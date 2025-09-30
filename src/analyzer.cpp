#include <limits.h>

#include <algorithm>
#include <math.h>
#include <numeric>
#include <omp.h>
#include <typeinfo>
#include "analyzer.h"
#include <chrono>
#include <iomanip>
#include "par.h"

//#define NUM_THREADS 4

int* Analyzer::getLevels() {
  return levels;
}

int Analyzer::getNumOfLevels() {
  return numOfLevels;
}

vector<vector<int>>& Analyzer::getLevelTable() {
  vector<vector<int>>& ref = levelTable;

  return ref;
}

DAG& Analyzer::getDAG() {
  DAG& ref = dag;

  return ref;
}

map<int, vector<int>>& Analyzer::getOriParents() {
  map<int, vector<int>>& ref = oriParents;

  return ref;
}

map<int, vector<double>>& Analyzer::getOriRowValues() {
  map<int, vector<double>>& ref = oriRowValues;

  return ref;
}

vector<vector<double>>& Analyzer::getValues() {
  vector<vector<double>>& ref = values;

  return ref;
}

vector<int>& Analyzer::getFlopsPerLevel() {
  vector<int>& ref = flopsPerLevel;

  return ref;
}

void Analyzer::setFlopsPerLevel(vector<int>& levelCost) {
  flopsPerLevel = levelCost;
}

map<int,int>& Analyzer::getFlopsBelowAvg() {
  map<int,int>& ref = flopsBelowAvg;

  return ref;
}

map<int,int>& Analyzer::getFlopsAboveAvg() {
  map<int,int>& ref = flopsAboveAvg;

  return ref;
}

bool Analyzer::getSingleLoopRows() {
  return singleLoopRows;
}

void Analyzer::setSingleLoopRows(bool looped) {
  singleLoopRows = looped;
}

float Analyzer::getALC() {
  return ALC;
}

float Analyzer::getAIR() {
  return AIR;
}

float Analyzer::getARL() {
  return ARL;
}

int Analyzer::getMMAD() {
  return MMAD;
}

int Analyzer::getMID() {
  return MID;
}

float Analyzer::getAMAD() {
  return AMAD;
}


#ifdef REWRITE_ENABLED
  vector<int>& Analyzer::getFlopsPerLevelRewrite() {
    vector<int>& ref = flopsPerLevelRewrite;

    return ref;
  }

  int Analyzer::findMaxLevelOfPredecessors(int row) {
    int maxLevel = 0;
    vector<int>& parents = dag[row].first;

    for(auto& parent : parents) {
      if(maxLevel < levels[parent]) 
        maxLevel = levels[parent];
    }
  
    return maxLevel;
  }

#endif

  // to avoid wrong calculations sort emptyLevels in descending order
  // before calling this function
  void Analyzer::updateWithEmptyLevels(vector<int>& emptyLevels) {
    numOfLevels -= emptyLevels.size();

/*auto erase_from_vector = [&](auto& vec) { 
  size_t j = 0;
  for(size_t i = 0; i < vec.size(); i++) {
    if(std::find(emptyLevels.begin(), emptyLevels.end(), i) == emptyLevels.end()) {
      vec[j++] = std::move(vec[i]);
    }
  }
  vec.resize(j);
};

auto start = std::chrono::steady_clock::now();
erase_from_vector(levelTable);
erase_from_vector(flopsPerLevel);
auto end = std::chrono::steady_clock::now();
cout << "chrono update:erase: " << chrono::duration<double>(end-start).count()*1000 << "\n";
*/
//auto start = std::chrono::steady_clock::now();
    for(auto& level : emptyLevels) {
//      cout << "erasing " << level <<  " at index " << level << "\n";
      levelTable.erase(levelTable.begin() + level);
      flopsPerLevel.erase(flopsPerLevel.begin() + level);
    }

//auto end = std::chrono::steady_clock::now();
//cout << "chrono update:erase: " << chrono::duration<double>(end-start).count()*1000 << "\n";

    // recalculate total and cost, ALC, ARL
    totalFLOPSPerLevel = accumulate(flopsPerLevel.begin(), flopsPerLevel.end(), 0.0);
    separateThinLevels();
    ARL = 0; ALC = 0;

    for(auto& level: flopsAboveAvg) {
      ARL += levelTable[level.first].size();
      ALC += level.second;
    }

    ARL /= flopsAboveAvg.size();
    ARL = ceil(ARL);

    ALC /= flopsAboveAvg.size();
    ALC = ceil(ALC);

 /*   ALC *= 2;
    ARL *= 2;*/
  }

void Analyzer::buildLevels() {
  if(matrixCSC != nullptr && matrixCSR != nullptr) {
    Part* L = matrixCSC->getL();
    vector<int>& rowPtrL = L->getRowPtr();
    vector<int>& colIdxL = L->getColIdx();
    vector<double>& valsL = L->getVals();
    vector<int>& rowPtrLCSR = matrixCSR->getL()->getRowPtr();
    int cols = L->getCols();
    int rows = L->getRows();

    int succLevels = posix_memalign((void **)&levels, 64, (cols-1)*sizeof(int));
    CHK_MEM_ERR(succLevels,"levels");

    std::fill(levels, levels + cols-1, 0);

    // create levels array:
    // index: row, value: level that the row is in
    int maxLevel = levels[0];
    for(int i = 0; i < cols-1; i++) {
      int currLevel = levels[rowPtrL[colIdxL[i]]];

      Connectivity& connectivity = dag[i];
      vector<int>& children = connectivity.second;

      for(int j = colIdxL[i]+1 ; j < colIdxL[i+1]; j++) {
        if(levels[rowPtrL[j]] <= currLevel) {
          levels[rowPtrL[j]] = currLevel+1;
          if(levels[rowPtrL[j]] > maxLevel)
              maxLevel = levels[rowPtrL[j]];
        }
     
        children.push_back(rowPtrL[j]);

        Connectivity& connectivityChild = dag[rowPtrL[j]];
        vector<int>& parents = connectivityChild.first;
        parents.push_back(i);

        vector<double>& rowValues = values[rowPtrL[j]];
        rowValues.push_back(valsL[j]);
      }
    }

    // fill in the value for the diagonal
    for(int i = 0; i < cols-1; i++) {
      vector<double>& rowValues = values[rowPtrL[colIdxL[i]]];
      rowValues.push_back(valsL[colIdxL[i]]);
    }

    // create levelTable
    levelTable = vector<vector<int>>(maxLevel + 1);
    for(int k = 0 ; k < cols-1 ; k++) {
      vector<int>& level = levelTable[levels[k]];
      level.push_back(k);
    }

    numOfLevels = maxLevel+1;
    levelTable.shrink_to_fit();
  } else {
    printf("%d: Wrong format for level building. It has to be CSC\n", __LINE__);
  }

  auto start = std::chrono::steady_clock::now();
  ARL = 0; AIR = 0;
  int size = levelTable.size();
  #ifdef PAR
   #pragma omp parallel for num_threads(NUM_THREADS) if(size > WL_SIZE) reduction(+:ARL,AIR)
  #endif
  //for(auto& level : levelTable) {
  for(int i = 0; i < levelTable.size(); i++) {
    vector<int>& level = levelTable[i];
    ARL += level.size();
    for(auto& row : level)
      AIR += dag[row].first.size();
  }
  auto end = std::chrono::steady_clock::now();
  #if defined(PAR) && (size > WL_SIZE)
    cout << "* chrono calculate buildLevels ARL;AIR: " << chrono::duration<double>(end-start).count()*1000 << "\n";
  #else
    cout << "chrono calculate buildLevels ARL;AIR: " << chrono::duration<double>(end-start).count()*1000 << "\n";
  #endif

  ARL /= levelTable.size();
  AIR /= matrixCSC->getNumOfRows();

  pro_ARL = ARL;
  pro_AIR = AIR;
}

void Analyzer::separateRows(vector<int>& loopedRows, vector<int>& unrolledRows) {
  // if row length (a.k.a number of parents < 5) put them in unrolledRows
  // rest is in loopedRows, will be executed in 1 loop as level-set
  vector<int>::iterator it = loopedRows.begin();
  while(it != loopedRows.end()) {
    if(dag[*it].first.size() < UNROLL_FACTOR) {
      unrolledRows.push_back(*it);
      loopedRows.erase(it);
    } else it++;
  }
}

void Analyzer::calculateFLOPS() {
  flopsPerLevel.resize(numOfLevels, 0);
  fill(flopsPerLevel.begin(), flopsPerLevel.end(),0);

  auto start = std::chrono::steady_clock::now();
      
//  int part = (int)ceil((double)numOfLevels/NUM_THREADS);
  // didnt use the fi-clause in omp pragma since still we need if-else structure, it doesntbring any benefits.
  #ifdef PAR
    if(numOfLevels > WL_SIZE) {
//      printf("calculateFLOPS: num_threads:%d\n", NUM_THREADS);
//      printf("numOfLevels: %d\n", numOfLevels);
      #pragma omp parallel num_threads(NUM_THREADS) reduction(+:totalFLOPSPerLevel)
      {
        int thread_id = omp_get_thread_num();
        int start = (thread_id * numOfLevels) / NUM_THREADS;
        if(thread_id == 0) start = 1;
        int end = ((thread_id + 1) * numOfLevels) / NUM_THREADS;
      
//        printf("thread id: %d start: %d, end: %d\n", thread_id, start, end);
        for(int i = start; i < end; i++) {
          int indegreeSum = 0;
          vector<int>& level = levelTable[i];
          for(auto& row : level)
            indegreeSum += dag[row].first.size()+1;
    
          flopsPerLevel[i] = (indegreeSum << 1) - level.size();
          totalFLOPSPerLevel += (indegreeSum << 1) - level.size();
        }
      }
    } else {
      for(int i = 1; i < numOfLevels; i++) {
        int indegreeSum = 0;
        vector<int>& level = levelTable[i];
        for(auto& row : level)
          indegreeSum += dag[row].first.size()+1;
    
        //flopsPerLevel.push_back((indegreeSum << 1) - level.size());
        flopsPerLevel[i] = (indegreeSum << 1) - level.size();
    //    indegreeSum = 0;
      }
    }
  #else
  for(int i = 1; i < numOfLevels; i++) {
    int indegreeSum = 0;
    vector<int>& level = levelTable[i];
    for(auto& row : level) {
      indegreeSum += dag[row].first.size()+1;
    }

    //flopsPerLevel.push_back((indegreeSum << 1) - level.size());
    flopsPerLevel[i] = (indegreeSum << 1) - level.size();
//    indegreeSum = 0;
  }
  #endif
  auto end = std::chrono::steady_clock::now();
    if(numOfLevels > WL_SIZE)
      cout << "* chrono calculateFLOPS : " << chrono::duration<double>(end-start).count()*1000 << "\n";
    else
      cout << "chrono calculateFLOPS : " << chrono::duration<double>(end-start).count()*1000 << "\n";


  //flopsPerLevel.shrink_to_fit();
  totalFLOPSPerLevel = accumulate(flopsPerLevel.begin(), flopsPerLevel.end(), 0.0);
  ALC = totalFLOPSPerLevel/numOfLevels;
  pro_ALC = ALC;

//  ALC *= 2;
}

int Analyzer::recalculateFLOPSFor(vector<int> loopedRows){
  int indegreeSum = 0;
  for(auto& row: loopedRows)
    indegreeSum += dag[row].first.size()+1;
    
  return ((indegreeSum << 1) - loopedRows.size());
}

void Analyzer::analyzeForCriteriaCostMap() {
  for(auto& level : flopsAboveAvg)
    ALC += level.second;

  ALC /= flopsAboveAvg.size();
  ALC = ceil(ALC);
}


//threeCriteria
void Analyzer::analyzeForCriteria3CRI() {
  AIR = ARL = MMAD = MID = 0;
  for(int i = 1; i < numOfLevels; i++) {
    vector<int>& level = levelTable[i];

    for(auto& row : level) {
      vector<int>& parents = dag[row].first;
      AIR += parents.size();

    }
  }

  AIR /= matrixCSC->getNumOfRows(); 
  AIR = ceil(AIR);
  
  for(auto& level : flopsBelowAvg) {
    vector<int>& currLevel = levelTable[level.first];
    ARL += currLevel.size();
  }

  ARL /= flopsBelowAvg.size();
  ARL = ceil(ARL);

  /*// ARL for all
 *   ARL = ceil(matrixCSC->getNumOfRows()/levelTable.size());*/

  cout << "\nALC: " << ALC << "\n";
  cout << "AIR: " << AIR << " ARL: " << ARL << "\n";
  }


// NEW3CRI
void Analyzer::analyzeForCriteria() {
  AIR = ARL = MMAD = MID = 0;
  AMAD = 0; // too costly to calculate

  auto start = chrono::high_resolution_clock::now();
  #ifdef PAR
  #pragma omp parallel for num_threads(NUM_THREADS) if(numOfLevels > WL_SIZE) reduction(+:AIR)
  #endif
  for(int i = 1; i < numOfLevels; i++) {
    vector<int>& level = levelTable[i];

    for(auto& row : level) {
      vector<int>& parents = dag[row].first;
      AIR += parents.size();


     /* REMOVED FROM ANALYSIS SINCE THEY DIDN'T PAY OUT
      *
      // find max distance between indegrees
      auto minMax = minmax_element(parents.begin(), parents.end());
      int indegreeDistance = *minMax.second - *minMax.first;

      // find max distance between memory accesses (indegrees & row itself)
      int memAccessDistance;
      if(row < *minMax.first)
        memAccessDistance = *minMax.second - row;
      else if(row > *minMax.second)
        memAccessDistance = row - *minMax.first;

      if(MID < indegreeDistance) {
        //cout << "MID: " << MID << " indegreeDistance: " << indegreeDistance << "\n";
        MID = indegreeDistance;
      }

      if(MMAD < memAccessDistance) {
//        cout << "MMAD: " << MMAD << " memAccessDistance: " << memAccessDistance << "\n";
        MMAD = memAccessDistance;
      }*/
    }
  }

  ALC = 0; ARL = 0;
  for(auto& level : flopsAboveAvg) {
    vector<int>& currLevel = levelTable[level.first];

    ARL += currLevel.size();
    ALC += level.second;
  }

  ALC /= flopsAboveAvg.size();
  ALC = ceil(ALC);

  AIR /= matrixCSC->getNumOfRows(); 
  AIR = ceil(AIR);
  
  ARL /= flopsAboveAvg.size();
  ARL = ceil(ARL);

  // ARL for all
  int ARL_ALL = ceil(matrixCSC->getNumOfRows()/levelTable.size());
  //printf("ARL_ALL:,%d\n", ARL_ALL);

//  cout << "\nALC: " << ALC << "\n";
//  cout << "AIR: " << AIR << " ARL: " << ARL << "\n";
//  cout << "MMAD: " << MMAD << " MID: " << MID << "\n\n";

//  printf("ARL BEFORE:, %.2f\n", ARL);
  char bash_cmd[50] = "getconf _NPROCESSORS_ONLN";
  FILE* pipe;
  int numThreads = 0;
  pipe = popen(bash_cmd, "r");
  if(pipe != NULL)
    fscanf(pipe, "%d", &numThreads);
  ARL = max(max((int)ARL, ARL_ALL), numThreads);
  
  float NEW_ALC = ((((int)AIR<<1)+1) * matrixCSC->getNumOfRows())/numOfLevels;
//  printf("ALC based on AIR =, (%.2f * %d) / %d:, %.2f\n", AIR,  matrixCSC->getNumOfRows(), numOfLevels, NEW_ALC);
  if(ALC < NEW_ALC)
    ALC = NEW_ALC;
//  printf("ALC:,%.2f\nARL:,%.2f\nAIR:,%.2f\n\n", ALC, ARL, AIR);

  float sum = 0.0;
  #ifdef PAR
  #pragma omp parallel for num_threads(NUM_THREADS) if(numOfLevels > WL_SIZE) reduction(+:sum)
  #endif
  for(int i = 1; i < numOfLevels; i++) {
    vector<int>& level = levelTable[i];

    for(auto& row : level) {
      vector<int>& parents = dag[row].first;
      sum += (parents.size() - AIR) * (parents.size() - AIR);
    }
  }

  sum /= matrixCSC->getNumOfRows(); 
  AIR_CV = sqrt(sum)/AIR;
  if(AIR_CV  < 0.5)
    AIR = 1.5 * AIR;


  auto end = chrono::high_resolution_clock::now();
  #if defined(PAR) || (numOfLevels > WL_SIZE)
    cout << fixed << setprecision(4) << "* chrono analyzeForCriteriaNEW3CRI:, " << chrono::duration<double>(end-start).count() * 1000 << "\n";
  #else
  cout << fixed << setprecision(4) << "chrono analyzeForCriteriaNEW3CRI:, " << chrono::duration<double>(end-start).count() * 1000 << "\n";
  #endif
}

// 3CRI_IMPROVED
void Analyzer::analyzeForCriteriaImproved() {
  AIR = ARL = MMAD = MID = 0;
  AMAD = 0; // too costly to calculate

  auto start = chrono::high_resolution_clock::now();

  ALC = 0; ARL = 0;
  for(auto& level : flopsAboveAvg) {
    vector<int>& currLevel = levelTable[level.first];

    ARL += currLevel.size();
    ALC += level.second;
  }

  ALC /= flopsAboveAvg.size();
  ALC = ceil(ALC);

  ARL /= flopsAboveAvg.size();
  ARL = ceil(ARL);

  AIR = ALC/ARL;

  // ARL for all
  int ARL_ALL = ceil(matrixCSC->getNumOfRows()/levelTable.size());

  char bash_cmd[50] = "getconf _NPROCESSORS_ONLN";
  FILE* pipe;
  int numThreads = 0;
  pipe = popen(bash_cmd, "r");
  if(pipe != NULL)
    fscanf(pipe, "%d", &numThreads);
  ARL = max(max((int)ARL, ARL_ALL), numThreads);
  
  ALC = AIR * ARL;

  auto end = chrono::high_resolution_clock::now();
  #if defined(PAR) || (numOfLevels > WL_SIZE)
    cout << fixed << setprecision(4) << "* chrono analyzeForCriteriaImproved:, " << chrono::duration<double>(end-start).count() * 1000 << "\n";
  #else
  cout << fixed << setprecision(4) << "chrono analyzeForCriteriaImproved:, " << chrono::duration<double>(end-start).count() * 1000 << "\n";
  #endif
}

// 2CRI
void Analyzer::analyzeForCriteria2() {
  AIR = ARL = MMAD = MID = 0;
  AMAD = 0; // too costly to calculate

  auto start = chrono::high_resolution_clock::now();
  #ifdef PAR
  #pragma omp parallel for num_threads(NUM_THREADS) if(numOfLevels > WL_SIZE) reduction(+:AIR)
  #endif
  for(int i = 1; i < numOfLevels; i++) {
    vector<int>& level = levelTable[i];

    for(auto& row : level) {
      vector<int>& parents = dag[row].first;
      AIR += parents.size();
    }
  }

  ALC = 0; ARL = 0;
  for(auto& level : flopsAboveAvg) {
    vector<int>& currLevel = levelTable[level.first];

    ARL += currLevel.size();
  }

  AIR /= matrixCSC->getNumOfRows(); 
  AIR = ceil(AIR);
  
  float rowCost = AIR * 8;
  printf("rowCost: %.2lf\n", rowCost);

  ARL /= flopsAboveAvg.size();
  ARL = ceil(ARL);

  // ARL for all
  int ARL_ALL = ceil(matrixCSC->getNumOfRows()/levelTable.size());

  char bash_cmd[50] = "getconf _NPROCESSORS_ONLN";
  FILE* pipe;
  int numThreads = 0;
  pipe = popen(bash_cmd, "r");
  if(pipe != NULL)
    fscanf(pipe, "%d", &numThreads);
  ARL = max(max((int)ARL, ARL_ALL), numThreads);
  
  printf("ARL calculated: %.2lf\n", ARL);
  ALC = rowCost * ARL;

  auto end = chrono::high_resolution_clock::now();
  #if defined(PAR) || (numOfLevels > WL_SIZE)
    cout << fixed << setprecision(4) << "* chrono analyzeForCriteria2:, " << chrono::duration<double>(end-start).count() * 1000 << "\n";
  #else
  cout << fixed << setprecision(4) << "chrono analyzeForCriteria2:, " << chrono::duration<double>(end-start).count() * 1000 << "\n";
  #endif
}

void Analyzer::analyzeForCriteriaCoeff() {
  /*float sum_ARL = 0.0, sum_ALC = 0.0;

  float sum = 0.0;
  auto start = chrono::high_resolution_clock::now();

   for(auto& level : flopsAboveAvg) {
    vector<int>& currLevel = levelTable[level.first];
    sum_ARL += (currLevel.size() - ARL) * (currLevel.size() - ARL);

    sum_ALC += (level.second - ALC) * (level.second - ALC);
  }

  sum_ARL /= flopsAboveAvg.size();
  ARL_CV = sqrt(sum_ARL)/ARL;

  sum_ALC /= flopsAboveAvg.size();
  ALC_CV = sqrt(sum_ALC)/ALC;*/

  float sum_ARL = 0.0, sum_ALC = 0.0, sum = 0.0;

  auto start = chrono::high_resolution_clock::now();
  pro_AIR = pro_ALC = pro_ARL = 0.0;
  #ifdef PAR
  #pragma omp parallel for num_threads(NUM_THREADS) if(numOfLevels > WL_SIZE) reduction(+:sum)
  #endif
  for(int i = 0; i < levelTable.size(); i++) {
    vector<int>& level = levelTable[i];
    pro_ARL += level.size();
    for(auto& row : level)
      pro_AIR += dag[row].first.size();
  }

  pro_ARL /= levelTable.size();
  pro_AIR /= matrixCSC->getNumOfRows();
  float total = accumulate(flopsPerLevel.begin(), flopsPerLevel.end(), 0.0);
  pro_ALC = total/numOfLevels;

  #ifdef PAR
  #pragma omp parallel for num_threads(NUM_THREADS) if(numOfLevels > WL_SIZE) reduction(+:sum)
  #endif
  for(int i = 1; i < numOfLevels; i++) {
    vector<int>& level = levelTable[i];

    sum_ARL += (level.size() - pro_ARL) * (level.size() - pro_ARL);
    sum_ALC += (flopsPerLevel[i] - pro_ALC) * (flopsPerLevel[i] - pro_ALC);

    for(auto& row : level) {
      vector<int>& parents = dag[row].first;
      sum += (parents.size() - pro_AIR) * (parents.size() - pro_AIR);
    }
  }

  sum_ARL /= flopsPerLevel.size();
  ARL_CV = sqrt(sum_ARL)/pro_ARL;

  sum_ALC /= flopsPerLevel.size();
  ALC_CV = sqrt(sum_ALC)/pro_ALC;

  sum /= matrixCSC->getNumOfRows(); 
  AIR_CV = sqrt(sum)/pro_AIR;

  printf("COEFF OF VAR.\nALC_CV:, %.2f\nARL_CV:, %.2f\nAIR_CV:, %.2f\n", ALC_CV, ARL_CV, AIR_CV);

  auto end = chrono::high_resolution_clock::now();
//  printf("num of levels size: %zu\n", flopsAboveAvg.size());
  #if defined(PAR) || (numOfLevels > WL_SIZE)
  cout << fixed << setprecision(4) << "* chrono analyzeForCriteriaCoeff:, " << chrono::duration<double>(end-start).count() * 1000 << "\n";
  #else
  cout << fixed << setprecision(4) << "chrono analyzeForCriteriaCoeff:, " << chrono::duration<double>(end-start).count() * 1000 << "\n";
  #endif
}

void Analyzer::separateThinLevels() {
  flopsAboveAvg.clear(); flopsBelowAvg.clear();

  #ifdef PAR
  #pragma omp parallel num_threads(2)
  {
    #pragma omp sections
    {
      #pragma omp section
      for(int i = 0; i < numOfLevels; i++) 
        if(flopsPerLevel[i] < ALC)
           flopsBelowAvg[i] = flopsPerLevel[i];

      #pragma omp section
      for(int i = 0; i < numOfLevels; i++) 
        if(flopsPerLevel[i] >= ALC)
           flopsAboveAvg[i] = flopsPerLevel[i];
    }
  }
  #else
  for(int i = 0; i < numOfLevels; i++) 
    if(flopsPerLevel[i] < ALC)
     flopsBelowAvg[i] = flopsPerLevel[i];
    else
     flopsAboveAvg[i] = flopsPerLevel[i];
  #endif
}

#ifdef REWRITE_ENABLED
  void Analyzer::saveOriginalValues() {
    for(int i = 0 ; i < dag.size(); i++) {
      vector<int> parents(dag[i].first);
      oriParents[i] = parents;

      vector<double> rowValues(values[i]);
      oriRowValues[i] = rowValues;
    }
  }
#endif

void Analyzer::printLevels() {
  cout << "printLevels:\n";
  if(matrixCSC != nullptr) {
    int rows = matrixCSC->getNumOfRows();

    cout << "num of rows: " << rows << "\n";
    for(int i = 0 ; i < rows ; i++)
      cout << i << ":, " << levels[i] << "\n";
    cout << "\n";
  } else {
    printf("%d: Wrong format while printing levels. It has to be CSC\n", __LINE__);
  }
}

void Analyzer::printLevelTable() {
  cout << "printLevelTable:\n";
  int cnt = 0;
  if(matrixCSC != nullptr) {
    for(auto& level : levelTable) {
      cout << "level " << cnt++ << "\n";
      for(auto& row : level)
        cout << row << ", ";
      cout << "\n\n";
    }
  } else {
    printf("%d: Wrong format while printing level table.\n", __LINE__);
  }
}

void Analyzer::printLevelSizes(){
  cout << "printLevelSizes:\n";
  cout << "num. of levels:, " << levelTable.size() << "\n";
  if(matrixCSC != nullptr) {
   for(int i = 0 ; i < levelTable.size(); i++)
      cout << i << ":, " << levelTable[i].size() << "\n";
  } else {
    printf("%d: Wrong format while printing level table.\n", __LINE__);
  }
}

void Analyzer::printDAG() {
  cout << "printDAG:\n";
  for(int i = 0 ; i < dag.size(); i++) {
    cout << "row " << i << ":\n";
    Connectivity& connectivity = dag[i];
    vector<int>& parents = connectivity.first;
    vector<int>& children = connectivity.second;

    cout << "parents:\n";
    for(auto& parent : parents)
      cout << parent << ", ";
    cout << "\n\n";

    cout << "children:\n";
    for(auto& child : children)
      cout << child << ", ";
    cout << "\n\n";
  }
}

void Analyzer::printValues() {
  cout << "printValues:\n";
  for(int i = 0 ; i < values.size(); i++) {
    cout << "row " << i << ":\n";
    for(auto& value : values[i]) {
      cout << value << ", ";
    cout << "\n";
    }

    cout << "\n\n";
  }
}

void Analyzer::printFLOPSPerLevel() {
  cout << "printFLOPSPerLevel:\n";
  cout << "total cost:, " << totalFLOPSPerLevel << "\n";
//  cout << "ALL LEVELS: avg. cost per level(ALC):, " << ALC << "\n";
//  cout << "ALL LEVELS: avg. # of rows per level(ARL):, " << ARL << "\n";
//  cout << "ALL LEVELS: avg. # of indegrees per row(AIR):, " << AIR << "\n";
  cout << "cost(FLOPS) per Level:\n";

  for(int i = 0; i < flopsPerLevel.size(); i++)
    cout << i << ":, " << flopsPerLevel[i] << "\n";
  cout << "\n";
}

void Analyzer::printFLOPSDivided() {
  cout << "printFLOPSDivided\n";
  cout << "avg. FLOPS per Level(ALC): " << ALC << "\n"; 
  cout << "flopsBelowAvg:\n";
  for(auto& level : flopsBelowAvg)
    cout << level.first << " :, " << level.second << "\n";

  cout << "flopsAboveAvg:\n";
  for(auto& level : flopsAboveAvg)
    cout << level.first << " :, " << level.second << "\n";
}

void Analyzer::printCriteria(){
  cout << "printCriteria\n";
  cout << "ALC:, " << ALC << "\n";
  cout << "ARL:, " << ARL << "\n";
  cout << "AIR:, " << AIR << "\n\n";
}

// reportType: BEFORE/AFTER
// reportStep: 0 (default), 1 (more), 2 (more more)
void Analyzer::report(string reportType, int reportStep) {
  cout << "\nREPORT:, " << reportType << "\n";
  if(reportType.compare("BEFORE") == 0){
    cout << "num. of levels:, " << getNumOfLevels() << "\n";
    printCriteria();
    analyzeForCriteriaCoeff();
    printFLOPSPerLevel();
//    printLevelTable();
    printLevelSizes();
//    printValues();
//    printDAG();
    validityCheck();
  }

  if(reportType.compare("AFTER") == 0){
    #ifdef REWRITE_ENABLED
      separateThinLevels();
   //   analyzeForCriteria();

      cout << "num. of levels:, " << getNumOfLevels() << "\n";
   //   printCriteria();
      analyzeForCriteriaCoeff();
//    printValues();
      //printFLOPSDivided();
//      printLevelTable();
   //   printDAG();
    #endif
    printFLOPSPerLevel();
    printLevelSizes();
    validityCheck();
  }

  if(reportStep == 1) {
    printLevelTable();
  } else if(reportStep == 2) {
    printLevels();
    printDAG();
    printValues();
  } else if(reportStep != 0)
    cout << "wrong report step. 0 (default), 1 (more), 2 (more more)\n";
  cout << "REPORT:, " << reportType << "\n\n";
}

void Analyzer::printDependencies() {
  cout << "printDependencies:\n";
  Part* L = matrixCSR->getL();
  int rows = L->getRows();
  cout << "rows: " << rows << "\n";

  for(int i = 0; i < rows-1 ; i++) {
    vector<int>& parents = dag[i].first;
    cout << "row " << i << " has " << parents.size() << " dependencies:\n";
    for(auto& parent : parents) {
      cout << parent << " ,";
    cout << "\n";
    }
  }
}

void Analyzer::validityCheck() {
    cout << "validity check:\n";
    cout << "numOfLevels: " << numOfLevels << "\n";
    cout << "FlopsPerLevel size: " << getFlopsPerLevel().size() << "\n";
    cout << "flopsBelowAvg size: " << getFlopsBelowAvg().size() << "\n";
    cout << "flopsAboveAvg size: " << getFlopsAboveAvg().size() << "\n";
//    cout <<  getFlopsPerLevel().size() << " =? " << getFlopsBelowAvg().size() + getFlopsAboveAvg().size() << "\n"; 
    cout << "level table size: " << getLevelTable().size() << "\n";

    int numLevels = getFlopsPerLevel().size();
    int numLevelsBelow = getFlopsBelowAvg().size();
    int numLevelsAbove = getFlopsAboveAvg().size();

    if((numLevels == numLevelsBelow + numLevelsAbove) && (numLevels == numOfLevels))
      cout << "[CHECK]: num. of levels equal for FlopsPerLevel, FlopsBelowAvg, FlopsAboveAvg\n";
      

    // rows in levelTable, levels, DAG have to be equal (a.k.a nnzs)
    int numRowsTable = 0;
    for(auto& level: levelTable)
      numRowsTable += level.size();

    int numRowsDAG = dag.size();

    cout << "numRowsTable: " << numRowsTable << "\n";
    cout << "numRowsDAG: " << numRowsDAG << "\n";
    cout << "matrix num of cols: " <<  matrixCSR->getL()->getRows()-1 << "\n";

    if((numRowsTable == numRowsDAG) && (numRowsDAG == matrixCSR->getL()->getRows()-1))
      cout << "[CHECK]: numRows equal for levelTable, DAG, matrix num. of cols. levels* cannot be checked (malloc'ed)\n";

    int numParents = 0, numChildren = 0; 
    for(int i = 0 ; i < dag.size(); i++) {
      Connectivity& connectivity = dag[i];
      numParents += connectivity.first.size();
      numChildren += connectivity.second.size();
    }

    cout << "numParents: " << numParents << "\n";
    cout << "numChildren: " << numChildren << "\n";
    cout << "matrix num of vals: " <<  matrixCSR->getL()->getNNZs() << "\n";

    // the equalities will fail since nnz count might change after transformation
    // still the number of parents & childrem must be the same
    if(matrixCSR->getL()->getNNZs() == numParents + numRowsDAG)
      cout << "[CHECK]: DAG numParents + DAG size equal to matrix num of vals\n";
    if(matrixCSR->getL()->getNNZs() == numChildren + numRowsDAG)
      cout << "[CHECK]: DAG numChildren + DAG size equal to matrix num of vals\n";

    for(int i = 0 ; i < dag.size(); i++) {
      Connectivity& connectivity = dag[i];
      for(auto& parent : connectivity.first){
        auto it = find(dag[parent].second.begin(), dag[parent].second.end(), i);
        if(it == dag[parent].second.end()){
          cout << i << " parent " << parent << " doesnt have " << i <<  " as a child\n";
          cout << i << "'s parents:\n";
          for(auto& item : connectivity.first) 
            cout << item << ",";
          cout << "\n";
          cout << "parent " << parent << "'s children:\n";
          for(auto& item : dag[parent].second) 
            cout << item << ",";
          cout << "\n";
        }
      }
    }

}
