#include <utility>
#include "strategy.h"
#include "par.h"

// <row,<oldLevel,newLevel>>
typedef map<int,pair<int,int>> ToBeRewritten;
typedef map<int,set<int>> RewritingMap;

// indegree & outdegree vectors of a row: dependencies and children
typedef pair< vector<int>,vector<int> > Connectivity;
typedef vector<Connectivity> DAG;

ToBeRewritten& RewritingStrategy::getToBeRewritten() {
  ToBeRewritten& ref = toBeRewritten;

  return ref;
}

RewritingMap& RewritingStrategy::getRewritingMap() {
  RewritingMap& ref = rewritingMap;

  return ref;
}

bool RewritingStrategy::isRewritten(int row) {
  return toBeRewritten.find(row) != toBeRewritten.end();
}

bool RewritingStrategy::isTopDown() {
  return (startPoint == StartPoint::TopDown);
}

bool RewritingStrategy::isScopeSelective() {
  return (scope == Scope::Selective);
}

set<int>& RewritingStrategy::getRewritingMapOf(int row) {
  set<int>& ref = rewritingMap[row];

  return ref;
}

// duplicate function with Analyzer class, don't want to do it static
int RewritingStrategy::findMaxLevelOfPredecessors(int row) {
  int maxLevel = 0;
  vector<int>& parents = dag[row].first;

//  #pragma omp parallel for reduction(max:maxLevel)
   for(auto& parent : parents) {
     if(maxLevel < levels[parent])
       maxLevel = levels[parent];
   }

   return maxLevel;
}

void RewritingStrategy::shiftRowsUp() {
  if(startPoint == StartPoint::TopDown) {
    for(auto& row : toBeRewritten) {
      createRewritingMapFor(row.first, toBeRewritten[row.first].second);
      updateLevelOf(row.first, toBeRewritten[row.first].second);
      #ifdef UPDATE_GRAPES
        updateGrapeRows(row.first);
      #endif
    }
  } else {
 //   cout << "bottom up\n";
    for(auto it = toBeRewritten.rbegin() ; it != toBeRewritten.rend() ; ++it) {
      createRewritingMapFor(it->first, toBeRewritten[it->first].second);
      updateLevelOf(it->first, toBeRewritten[it->first].second);
      #ifdef UPDATE_GRAPES
        updateGrapeRows(it->first);
      #endif
    }
  }

  /*cout << "leveltable level 0 size: " << levelTable[0].size() << "\n";
  for(auto& row : levelTable[0])
    cout << row << ": " << dag[row].first.size() << "\n";
  cout << "\n";*/
}

void RewritingStrategy::updateLevelOf(int row, int newLevel) {
  int oldLevel = levels[row];
  levels[row] = newLevel;

  vector<int>& oLevel = levelTable[oldLevel];
  vector<int>& nLevel = levelTable[newLevel];

  vector<int>::iterator it = lower_bound(oLevel.begin(), oLevel.end(), row);
  if(it != oLevel.end())
    oLevel.erase(it);
  else
    cout << __LINE__ << ": couldnt erase row " << row << " at level " << oldLevel << "\n";

  it = upper_bound(nLevel.begin(), nLevel.end(), row);
  // no need since iterator is equal to end()
/*      if(it == nLevel.end())
    nLevel.push_back(row);
  else*/
    nLevel.insert(it,row);
}


void RewritingStrategy::createRewritingMapFor(int row, int targetLevel) {
   vector<int> predsWithMaxLevel;
   int maxLevel = findPredecessorsWithMaxLevel(dag[row].first, predsWithMaxLevel);
   rewritingMap[row] = set<int>();
 
   while(maxLevel > targetLevel) {
     for(auto& pred : predsWithMaxLevel)
       expandPredsWith(pred, dag[row].first, row);
 
     rewritingMap[row].insert(predsWithMaxLevel.begin(), predsWithMaxLevel.end());
     predsWithMaxLevel.clear();
     maxLevel = findPredecessorsWithMaxLevel(dag[row].first, predsWithMaxLevel);
   }
  
   // maxLevel == targetLevel
   for(auto& pred : predsWithMaxLevel)
     expandPredsWith(pred, dag[row].first, row);
    
  // update children of final parents to contain row
  for(auto& parent: dag[row].first) {
    vector<int>& children = dag[parent].second;
    if(find(children.begin(), children.end(), row) == children.end()) {
      vector<int>::iterator it = lower_bound(children.begin(), children.end(), row);
      children.insert(it, row);
    }
  }

  rewritingMap[row].insert(predsWithMaxLevel.begin(), predsWithMaxLevel.end());
}

void RewritingStrategy::expandPredsWith(int row, vector<int>& preds, int child) {
  auto it = lower_bound(preds.begin(),preds.end(),row);
  // Replace row with first pred. of it instead of erasing and causing a shift.
  // Then push back the rest since order doesnt matter.
  if(it != preds.end()) {
    vector<int>& parents = dag[row].first;

    if(!parents.empty()) { 
      preds[it-preds.begin()] = parents[0];
      preds.insert(preds.end(),parents.begin()+1,parents.end());

      sort(preds.begin(), preds.end());
      auto it = unique(preds.begin(), preds.end());
      preds.resize(std::distance(preds.begin(),it));
    } else
      preds.erase(it);

    vector<int>& children = dag[row].second;
    it = find(children.begin(), children.end(), child);
    if(it != children.end())
      children.erase(it);
  } else
    cout << __LINE__ << ": row " << row << " doesnt exist in predecessors\n";
}

int RewritingStrategy::findPredecessorsWithMaxLevel(vector<int>& preds, vector<int>& predsWithMaxLevel) {
  int maxLevel = 0;
  for(int j = 0; j < preds.size(); j++) {
    if(maxLevel < levels[preds[j]])
      maxLevel = levels[preds[j]];
  }

  for(auto& pred : preds)
    if(levels[pred] == maxLevel)
      predsWithMaxLevel.push_back(pred);

  return maxLevel;
}


// after rewriting strategy is applied, some levels will be emptied.
// reflect the changes to related data structures (e.g. levelTable)
// by calling Analyzer::updateWithEmptyLevels
void RewritingStrategy::findEmptyLevels(vector<int>& emptyLevels) {
  int numOfLevels = levelTable.size();
  for(vector<vector<int>>::reverse_iterator rit = levelTable.rbegin() ; rit != levelTable.rend() ; ++rit) {
    if(rit->empty()) {
   //   cout << "empty level: " << numOfLevels-1-(rit-levelTable.rbegin()) << "\n";
      emptyLevels.push_back(numOfLevels-1-(rit-levelTable.rbegin()));
    }
  }
}


void RewriteLongerThan::policy() {
  if(startPoint == StartPoint::TopDown)
    for(int k = 0 ; k < rows ; k++) {
      if(dag[k].first.size() >= REWRITE_UP &&
          levels[k] > 1 && numOfRewrittenRows-- > 0)
         toBeRewritten[k] = make_pair(levels[k],-1);
    }
  else {
    for(int k = rows-1 ; k >= 0 ; k--) {
      if(dag[k].first.size() >= REWRITE_UP &&
         levels[k] > 1 && numOfRewrittenRows-- > 0)
         toBeRewritten[k] = make_pair(levels[k],-1);
    }
  }

  // we intend to rewrite rows * REWRITE_RATIO rows. The matrix might not have
  // that many rows with the criteria (>= REWRITE_UP) we set.
  numOfRewrittenRows = toBeRewritten.size();
  cout << "actual num. of rewritten rows: " << toBeRewritten.size() << "\n";
}


void RewriteShorterThan::policy() {
  if(startPoint == StartPoint::TopDown)
    for(int k = 0 ; k < rows ; k++) {
      if(dag[k].first.size() <= REWRITE_LOW &&
         levels[k] > 1 && numOfRewrittenRows-- > 0)
         toBeRewritten[k] = make_pair(levels[k],-1);
    }
  else
    for(int k = rows-1 ; k >= 0 ; k--) {
      if(dag[k].first.size() <= REWRITE_LOW &&
         levels[k] > 1 && numOfRewrittenRows-- > 0)
         toBeRewritten[k] = make_pair(levels[k],-1);
    }

  numOfRewrittenRows = toBeRewritten.size();
  cout << "actual num. of rewritten rows: " << toBeRewritten.size() << "\n";
}

void RewriteInBetween::policy() {
  if(startPoint == StartPoint::TopDown)
    for(int k = 0 ; k < rows ; k++) {
      if(dag[k].first.size() >= REWRITE_UP && dag[k].first.size() <= REWRITE_LOW &&
         levels[k] > 1 && numOfRewrittenRows-- > 0)
         toBeRewritten[k] = make_pair(levels[k],-1);
    }
  else {
    for(int k = rows-1 ; k >= 0 ; k--) {
      if(dag[k].first.size() >= REWRITE_UP && dag[k].first.size() <= REWRITE_LOW &&
         levels[k] > 1 && numOfRewrittenRows-- > 0)
         toBeRewritten[k] = make_pair(levels[k],-1);
    }
  }

  numOfRewrittenRows = toBeRewritten.size();
  cout << "actual num. of rewritten rows: " << toBeRewritten.size() << "\n";
}

// rewrite rows with the smallest # of in-degree at selected levels
// TODO: implement
void RewriteLowestIndegree::policy() {
}


/******************** RewriteByCostMap *********************/
void RewriteByCostMap::policy() {
  // start from level 1 and try rewriting to level 0 prioritizing the levels close to each other
  auto it_source = flopsBelowAvg.begin(); 
  auto it_target = it_source;
  it_source++;
  for(; it_source != flopsBelowAvg.end(); it_source++) {
    int level = it_source->first;

    if(it_source->first - it_target->first > REWRITING_DIST) {
      it_target = it_source;
      it_source++;
      level = it_source->first;
    }
      
    vector<int>& levelRows = levelTable[level];
    auto nextRow = levelRows.begin();

    createCostMapMax(it_source, it_target, nextRow);
  }
}

void RewriteByCostMap::policy(int start, int end, ToBeRewritten& rewritingList) {
  // start from level 1 and try rewriting to level 0 prioritizing the levels close to each other
  auto it_source = flopsBelowAvg.begin(); advance(it_source, start); 
  auto end_it = flopsBelowAvg.begin();    advance(end_it, end);
  auto it_target = it_source;
  it_source++;
  for(; it_source != end_it; it_source++) {
    int level = it_source->first;

    if(it_source->first - it_target->first > REWRITING_DIST) {
      it_target = it_source;
      it_source++;
      level = it_source->first;
    }
      
    vector<int>& levelRows = levelTable[level];
    auto nextRow = levelRows.begin();

    createCostMapMax(it_source, it_target, nextRow, rewritingList);
  }
}

void RewriteByCostMap::distPolicy() {
  ToBeRewritten rewritingList[NUM_THREADS];

  #pragma omp parallel num_threads(NUM_THREADS)
  {
    int thread_id = omp_get_thread_num();
    int start = (thread_id * flopsBelowAvg.size()) / NUM_THREADS;
    int end = ((thread_id + 1) * flopsBelowAvg.size()) / NUM_THREADS;

//    printf("thread id: %d start: %d, end: %d\n", thread_id, start, end);
    policy(start, end, rewritingList[thread_id]); 
  }

  // flatten out rewritingList array into toBeRewritten
  // needed for experimenting. If parallelizing pays out we'll make toBeRewritten an array permenantly.
  for(auto& list : rewritingList)
    toBeRewritten.merge(list);
}

// PAR: VERSION works on a level
void RewriteByCostMap::createCostMapMax(map<int,int>::iterator& sourceLevel, map<int,int>::iterator& targetLevel, vector<int>::iterator& nextRow, ToBeRewritten& rewritingList) {
  int level = sourceLevel->first;
  int levelEnd = targetLevel->first;
  vector<int>& levelRows = levelTable[level];

  for(; nextRow != levelRows.end() ; ++nextRow) {
     vector<int> predsWithMaxLevel;
     // do not alter the dag
     vector<int> parents = dag[*nextRow].first;
     int maxLevel = findPredecessorsWithMaxLevel(parents, predsWithMaxLevel);

     while(maxLevel > levelEnd) {
       for(auto& pred : predsWithMaxLevel)
         expandPredsWith(pred, parents, *nextRow);

       predsWithMaxLevel.clear();
       maxLevel = findPredecessorsWithMaxLevel(parents, predsWithMaxLevel);
     }
    
     // maxLevel == targetLevel
     for(auto& pred : predsWithMaxLevel)
       expandPredsWith(pred, parents, *nextRow);

     int numOfPreds = parents.size();
     int costRow = 0;
     if(numOfPreds > 0) {
       costRow = (numOfPreds << 1);
       if(numOfPreds >= 5 ) // looped
         costRow++;
     }

     if((levelCost[maxLevel] + costRow) <= avgCostPerLevel) {
       levelCost[maxLevel] += costRow;
       rewritingList[*nextRow] = make_pair(level,levelEnd);
     } else { // target level is full, continue with the next level as the target level
       // if rewriting stopped in the middle of a level remove already rewritten rows' cost from levelCost
       if(nextRow != levelRows.begin()) {
         float cost = 0;
         for(auto it = levelRows.begin(); it != nextRow; it++) {
           int numOfPreds = dag[*it].first.size();
           cost += (numOfPreds << 1) + 1; 
         }
         levelCost[level] -= cost;
       }

       targetLevel = sourceLevel;

       return;
     }
  }
}

// works on a level
void RewriteByCostMap::createCostMapMax(map<int,int>::iterator& sourceLevel, map<int,int>::iterator& targetLevel, vector<int>::iterator& nextRow) {
  int level = sourceLevel->first;
  int levelEnd = targetLevel->first;
  vector<int>& levelRows = levelTable[level];

  for(; nextRow != levelRows.end() ; ++nextRow) {
     vector<int> predsWithMaxLevel;
     // do not alter the dag
     vector<int> parents = dag[*nextRow].first;
     int maxLevel = findPredecessorsWithMaxLevel(parents, predsWithMaxLevel);

     while(maxLevel > levelEnd) {
       for(auto& pred : predsWithMaxLevel)
         expandPredsWith(pred, parents, *nextRow);

       predsWithMaxLevel.clear();
       maxLevel = findPredecessorsWithMaxLevel(parents, predsWithMaxLevel);
     }
    
     // maxLevel == targetLevel
     for(auto& pred : predsWithMaxLevel)
       expandPredsWith(pred, parents, *nextRow);

     int numOfPreds = parents.size();
     int costRow = 0;
     if(numOfPreds > 0) {
       costRow = (numOfPreds << 1);
       if(numOfPreds >= 5 ) // looped
         costRow++;
     }

     if((levelCost[maxLevel] + costRow) <= avgCostPerLevel) {
       levelCost[maxLevel] += costRow;
       toBeRewritten[*nextRow] = make_pair(level,levelEnd);
     } else { // target level is full, continue with the next level as the target level
       // if rewriting stopped in the middle of a level remove already rewritten rows' cost from levelCost
       if(nextRow != levelRows.begin()) {
         float cost = 0;
         for(auto it = levelRows.begin(); it != nextRow; it++) {
           int numOfPreds = dag[*it].first.size();
           cost += (numOfPreds << 1) + 1; 
         }
         levelCost[level] -= cost;
       }

       targetLevel = sourceLevel;

       return;
     }
  }
}

// version of expandPredsWith that does not alter the dag (children)
void RewriteByCostMap::expandPredsWith(int row, vector<int>& preds, int child) {
  auto it = lower_bound(preds.begin(),preds.end(),row);
  // Replace row with first pred. of it instead of erasing and causing a shift.
  // Then push back the rest since order doesnt matter.
  if(it != preds.end()) {
    vector<int>& parents = dag[row].first;

    if(!parents.empty()) { 
      preds[it-preds.begin()] = parents[0];
      preds.insert(preds.end(),parents.begin()+1,parents.end());

      sort(preds.begin(), preds.end());
      auto it = unique(preds.begin(), preds.end());
      preds.resize(std::distance(preds.begin(),it));
    } else
      preds.erase(it);
  } else
      cout << __LINE__ << ": row " << row << " doesnt exist in predecessors\n";
}

/******************** RewriteByThreeCriteria *********************/
void RewriteByThreeCriteria::policy() {
  auto it = flopsBelowAvg.begin();
  int targetLevel = it->first;
  it++;
  for(; it != flopsBelowAvg.end(); it++){
    vector<int>& levelRows = levelTable[it->first];
    int rewriteCount = 0;
    for(auto& row : levelRows){
      vector<int> predsWithMaxLevel;
      vector<int> parents = dag[row].first;
      int originalCost = (parents.size() << 1) + 1;
      int maxLevel = findPredecessorsWithMaxLevel(parents, predsWithMaxLevel);
       
      while(maxLevel > targetLevel) {
        for(auto& pred : predsWithMaxLevel)
          expandPredsWith(pred, parents, row);

        predsWithMaxLevel.clear();
        maxLevel = findPredecessorsWithMaxLevel(parents, predsWithMaxLevel);
      }

       // maxLevel == targetLevel
      for(auto& pred : predsWithMaxLevel)
        expandPredsWith(pred, parents, row);

       int numOfPreds = parents.size();
       int costRow = 0;
       if(numOfPreds > 0) {
         costRow = (numOfPreds << 1);
         if(numOfPreds >= 5 ) // looped
           costRow++;
       }

       if(levelCost[maxLevel]+costRow < avgCostPerLevel) {
         // relaxing the indegree constraint with rewriting distance == 1
//       if((parents.size() <= avgIndegreePerRow) || (it->first - maxLevel == 1)) {
         if((parents.size() <= avgIndegreePerRow) && (it->first - maxLevel < REWRITING_DIST)
             || (it->first - maxLevel == 1)) {
            rewriteCount++;

           levelCost[maxLevel] += costRow;
           levelSizeBelowAvg[maxLevel]++;
           levelCost[it->first] -= originalCost;
           levelSizeBelowAvg[it->first]--;

           toBeRewritten[row] = make_pair(levels[row],maxLevel);

           if(levelSizeBelowAvg[maxLevel] == avgLevelSize) {
             #ifdef FAILED_ANALYSIS
               failedRowsSize.push_back(maxLevel);
             #endif

             if(levelSizeBelowAvg[it->first] == 0 && next(it) != flopsBelowAvg.end())
               it++;
             while((levelSizeBelowAvg[it->first] >= avgLevelSize) && (it->first - maxLevel < REWRITING_DIST))
               it++;
             targetLevel = it->first;
             break;
           }
         }
         #ifdef FAILED_ANALYSIS
         else if(parents.size() > avgIndegreePerRow)
           failedRowsIndegree[row] = make_pair(it->first,parents.size());
         #endif
       } 
       #ifdef FAILED_ANALYSIS
         else 
           failedRowsCost[row] = make_pair(it->first,costRow);
       #endif
    } // for each row

    if(rewriteCount == 0) { // no rewrite happened for this level. change the target level
      while((levelSizeBelowAvg[it->first] >= avgLevelSize) && (it->first - targetLevel < REWRITING_DIST))
       it++;
      targetLevel = it->first;
    }
  }  // for each level
}

void RewriteByThreeCriteria::distPolicy() {
  ToBeRewritten rewritingList[NUM_THREADS];

  #pragma omp parallel num_threads(NUM_THREADS)
  {
    int thread_id = omp_get_thread_num();
    int start = (thread_id * flopsBelowAvg.size()) / NUM_THREADS;
    int end = ((thread_id + 1) * flopsBelowAvg.size()) / NUM_THREADS;

//    printf("flopsBelowAvg size:%d\n", flopsBelowAvg.size());
//    printf("thread id: %d start: %d, end: %d\n", thread_id, start, end);
    policy(start, end, rewritingList[thread_id]);
//    printf("thread %d finished\n", thread_id);
  }

  // flatten out rewritingList array into toBeRewritten
  // needed for experimenting. If parallelizing pays out we'll make toBeRewritten an array permenantly.
  for(auto& list : rewritingList)
    toBeRewritten.merge(list);
}

void RewriteByThreeCriteria::policy(int start, int end, ToBeRewritten& rewritingList) {
  auto it = flopsBelowAvg.begin(); advance(it, start);
  auto end_it = flopsBelowAvg.begin();    advance(end_it, end);
  auto it_target = it;

 // printf("start: %d, end:%d, INITIAL: %d, %d\n", start, end, it->first, end_it->first);
  int targetLevel = it->first;
  it++;

/*  if(end == flopsBelowAvg.size()) {
    end_it = flopsBelowAvg.end();
//    printf("last level: %d\n", flopsBelowAvg.rbegin()->first);
  }*/
  for(; it != end_it; it++) {

/*  auto it = flopsBelowAvg.begin();
  int targetLevel = it->first;
  it++;
  for(; it != flopsBelowAvg.end(); it++){*/
//    printf("level:%d last level:%d\n", it->first, end_it->first);

    vector<int>& levelRows = levelTable[it->first];

    int rewriteCount = 0;
    for(auto& row : levelRows){
      vector<int> predsWithMaxLevel;
      vector<int> parents = dag[row].first;
      int originalCost = (parents.size() << 1) + 1;
      int maxLevel = findPredecessorsWithMaxLevel(parents, predsWithMaxLevel);
       
      while(maxLevel > targetLevel) {
        for(auto& pred : predsWithMaxLevel)
          expandPredsWith(pred, parents, row);

        predsWithMaxLevel.clear();
        maxLevel = findPredecessorsWithMaxLevel(parents, predsWithMaxLevel);
      }

       // maxLevel == targetLevel
      for(auto& pred : predsWithMaxLevel)
        expandPredsWith(pred, parents, row);

       int numOfPreds = parents.size();
       int costRow = 0;
       if(numOfPreds > 0) {
         costRow = (numOfPreds << 1);
         if(numOfPreds >= 5 ) // looped
           costRow++;
       }

       if(levelCost[maxLevel]+costRow < avgCostPerLevel) {
         // relaxing the indegree constraint with rewriting distance == 1
//       if((parents.size() <= avgIndegreePerRow) || (it->first - maxLevel == 1)) {
         if((parents.size() <= avgIndegreePerRow) && (it->first - maxLevel < REWRITING_DIST)
             || (it->first - maxLevel == 1)) {
            rewriteCount++;

           levelCost[maxLevel] += costRow;
           levelSizeBelowAvg[maxLevel]++;
           levelCost[it->first] -= originalCost;
           levelSizeBelowAvg[it->first]--;

           //toBeRewritten[row] = make_pair(levels[row],maxLevel);
           rewritingList[row] = make_pair(levels[row],maxLevel);

           if(levelSizeBelowAvg[maxLevel] == avgLevelSize) {
             #ifdef FAILED_ANALYSIS
               failedRowsSize.push_back(maxLevel);
             #endif

  //             printf("levelSizeBelowAvg[maxLevel] == avgLevelSize\n");
             //if(levelSizeBelowAvg[it->first] == 0 && next(it) != flopsBelowAvg.end())
             if(levelSizeBelowAvg[it->first] == 0 && next(it) != end_it)
               it++;
             while((levelSizeBelowAvg[it->first] >= avgLevelSize) && 
                    (it->first - maxLevel < REWRITING_DIST) &&
                    next(it) != end_it)
               it++;
             targetLevel = it->first;
             break;
           }
         }
         #ifdef FAILED_ANALYSIS
         else if(parents.size() > avgIndegreePerRow)
           failedRowsIndegree[row] = make_pair(it->first,parents.size());
         #endif
       } 
       #ifdef FAILED_ANALYSIS
         else 
           failedRowsCost[row] = make_pair(it->first,costRow);
       #endif
    } // for each row

    if(rewriteCount == 0) { // no rewrite happened for this level. change the target level
  //    printf("rewriteCount: level:%d\n", it->first);
      while((levelSizeBelowAvg[it->first] >= avgLevelSize) && 
            (it->first - targetLevel < REWRITING_DIST) &&
            next(it) != end_it)
       it++;
      targetLevel = it->first;
    }
  }  // for each level
}

// version of expandPredsWith that does not alter the dag (children)
  void RewriteByThreeCriteria::expandPredsWith(int row, vector<int>& preds, int child) {
    auto it = lower_bound(preds.begin(),preds.end(),row);
    // Replace row with first pred. of it instead of erasing and causing a shift.
    // Then push back the rest since order doesnt matter.
    if(it != preds.end()) {
      vector<int>& parents = dag[row].first;

      if(!parents.empty()) {
        preds[it-preds.begin()] = parents[0];
        preds.insert(preds.end(),parents.begin()+1,parents.end());

        sort(preds.begin(), preds.end());
        auto it = unique(preds.begin(), preds.end());
        preds.resize(std::distance(preds.begin(),it));
      } else
        preds.erase(it);
    } else
        cout << __LINE__ << ": row " << row << " doesnt exist in predecessors\n";
  }

// rewrite rows on critical path (CP)
// TODO: implement
    void RewriteOnCP::policy() {
      int pathWithMaxLength = 0;
      int pathLength = 0;
      for(int i = 0 ; i < paths.size(); i++) {
        vector<int>& path = paths[i];
        if(pathLength < path.size()) {
          pathLength = path.size();
          pathWithMaxLength = i;
        }
      }
       
      vector<int>& criticalPath = paths[pathWithMaxLength];
      numOfRewrittenRows = criticalPath.size();
      for(auto& row : criticalPath)
        toBeRewritten[row] = make_pair(levels[row],-1); 
    }

    void RewriteOnCP::traversePath(vector<int>& path, int row) {
      int parentMaxLevel = row;
      int maxLevel = 0;
      if(levels[row] > 0)
        for(auto& parent : dag[row].first) {
          if(levels[parent] > maxLevel) {
            maxLevel = levels[parent];
            parentMaxLevel = parent; 
          }
        }

      path.push_back(parentMaxLevel);

      if(levels[parentMaxLevel] > 0)
        traversePath(path, parentMaxLevel);
    }

    void RewriteOnCP::buildCriticalPath(vector<int>& lastLevelRows) {
      #pragma omp parallel for schedule(static)
      for(int i = 0 ; i < lastLevelRows.size() ; i++) {
        int row = lastLevelRows[i];
        vector<int> path;
        path.push_back(row);
        traversePath(path, row);

        if(path.back() != 0)
          cout << __LINE__ << ": path doesnt end in 0\n";

        mutex.lock();
        paths.push_back(path);
        mutex.unlock();
      }
    }

    

// currently we're merging all levels with the min. # of rows for simplicity
void RewriteByLevel::policy() {
  auto it = levelLengthHistogram.begin();
  for(auto& level : it->second)
    levelsToBeRewritten.push_back(level);

  // if there's only 1 level with the min. # of rows, check for levels with the second min. # of rows
  if(levelsToBeRewritten.size() == 1) {
    it++;
    for(auto& level : it->second)
      levelsToBeRewritten.push_back(level);
  }

  sort(levelsToBeRewritten.begin(),levelsToBeRewritten.end());


  // auto itTargetLevel = levelsToBeRewritten.rbegin();
  auto itTargetLevel = levelsToBeRewritten.begin();
  // skip the first level since it's our targetLevel
  auto itrback = levelsToBeRewritten.rend();
  itrback--;
  for(auto it = levelsToBeRewritten.rbegin() ; it != itrback; it++) {
    // together with commented line before this for loop where itTargetLevel is set to levelsToBeRewritten.rbegin()
    // this if condition makes each level to be shifted to the previous level in levelsToBeRewritten
    /*if(++itTargetLevel == levelsToBeRewritten.rend())
      break; */

    for(auto& row : levelTable[*it])
      toBeRewritten[row] = make_pair(levels[row],*itTargetLevel); 
  }

    cout << __LINE__ << " num. of rewritten rows: " << toBeRewritten.size() << "\n";
    printRowsToBeRewritten();
}

void RewriteByLevel::prepare(map<int,int>& flopsBelowAvg) {
  // create a histogram of num. of levels with each num. of rows in flopsBelowAvg
  // merge all levels with the minimum num. of rows within themselves 
  for(auto it = flopsBelowAvg.begin() ; it != flopsBelowAvg.end(); ++it) {
    vector<int>& entry = levelLengthHistogram[levelTable[it->first].size()]; // key to histogram: level's size
    entry.push_back(it->first);
  }

/*  cout << "levelLengthHistogram:\n";
  for(auto& levelSize : levelLengthHistogram) {
    vector<int>& levels = levelSize.second;
    cout << "size: " << levelSize.first << "\n";
    for(auto& level : levels)
      cout << level << ", ";
    cout << "\n";
  }*/
}

  // this is for calculating for the manual approach (rewrite by level) when
  // calculation in rewrite module is unfunctional
  void RewriteByLevel::calculateFLOPSAfterRewrite(vector<int>& flopsPerLevelRewrite) {
cout << levelTable.size() << ", " << flopsPerLevelRewrite.size() << "\n";
for(int i = 0 ; i < levelTable.size(); i++) {
  vector<int>& level = levelTable[i];

  flopsPerLevelRewrite[i] = 0;
  for(auto& row : level) {
    if(isRewritten(row)) {
      int numOfPreds = dag[row].first.size();
      if(numOfPreds == 0)
        flopsPerLevelRewrite[i]++;
      else
        flopsPerLevelRewrite[i] += numOfPreds << 1;
    } else {
      flopsPerLevelRewrite[i] +=  (dag[row].first.size() << 1) + 1;
    }
  }
}

cout << "FLOPS after rewrite:\n";
for(auto& flops : flopsPerLevelRewrite)
  cout << flops << "\n";
  }

//          int costRow = numOfPreds <= 4 ? (numOfPreds << 1) + 1 : (numOfPreds << 1);
void RewriteByLevel::torso2_9levels_to10th(map<int,int>& flopsBelowAvg) {
  auto targetLevel = flopsBelowAvg.begin();

  auto it = flopsBelowAvg.begin();
  it++;
  for(; it != flopsBelowAvg.end() ; it++) {
    if(it->first % 10 == 0) {
      advance(targetLevel,10);
    } else if (it->first != targetLevel->first) {
      levelsToBeRewritten.push_back(it->first);

      for(auto& row : levelTable[it->first])
         toBeRewritten[row] = make_pair(it->first,targetLevel->first);
    }
  }

  cout << "done\n";
}

void RewriteByLevel::lung2_9levels_to10th() {
  // lung2 - 9 levels written to 10th
  levelsToBeRewritten.push_back(1);
  levelsToBeRewritten.push_back(2);
  int targetRow = 12;
  for(int i = 3; i < 452; i++) {
    if(i != targetRow)
      levelsToBeRewritten.push_back(i);
    else
      targetRow+=10;
  }
 
  levelsToBeRewritten.push_back(453);
  levelsToBeRewritten.push_back(454);
  levelsToBeRewritten.push_back(455);
  levelsToBeRewritten.push_back(477);
  levelsToBeRewritten.push_back(478);

  for(auto& row : levelTable[478])
     toBeRewritten[row] = make_pair(478,476);
  for(auto& row : levelTable[477])
     toBeRewritten[row] = make_pair(477,476);
  for(auto& row : levelTable[455])
     toBeRewritten[row] = make_pair(455,452);
  for(auto& row : levelTable[454])
     toBeRewritten[row] = make_pair(454,452);
  for(auto& row : levelTable[453])
     toBeRewritten[row] = make_pair(453,452);

  targetRow = 442;
  for(int i = 451; i > 11; i--) {
    if(targetRow != i)
      for(auto& row : levelTable[i])
        toBeRewritten[row] = make_pair(i,targetRow);
    else
      targetRow-=10;
  }

  for(int i = 11; i > 0; i--)
    for(auto& row : levelTable[i])
       toBeRewritten[row] = make_pair(i,0);
}

void RewriteByLevel::lung2_7levels_to8th() {
  // lung2 - 7 levels written to 8th
  levelsToBeRewritten.push_back(1);
  levelsToBeRewritten.push_back(2);
  levelsToBeRewritten.push_back(3);
  int targetRow = 10;
  for(int i = 5; i < 452; i++) {
    if(i != targetRow)
      levelsToBeRewritten.push_back(i);
    else
      targetRow+=8;
  }

  levelsToBeRewritten.push_back(453);
  levelsToBeRewritten.push_back(454);
  levelsToBeRewritten.push_back(455);
  levelsToBeRewritten.push_back(478);

  for(auto& row : levelTable[478])
     toBeRewritten[row] = make_pair(478,477);
  for(auto& row : levelTable[455])
     toBeRewritten[row] = make_pair(455,452);
  for(auto& row : levelTable[454])
     toBeRewritten[row] = make_pair(454,452);
  for(auto& row : levelTable[453])
     toBeRewritten[row] = make_pair(453,452);

  targetRow = 444;
  for(int i = 451; i > 3; i--) {
    if(targetRow != i)
      for(auto& row : levelTable[i])
        toBeRewritten[row] = make_pair(i,targetRow);
    else
      targetRow-=8;
  }
      
  for(auto& row : levelTable[3])
     toBeRewritten[row] = make_pair(3,0);
  for(auto& row : levelTable[2])
     toBeRewritten[row] = make_pair(2,0);
  for(auto& row : levelTable[1])
     toBeRewritten[row] = make_pair(1,0);
}

void RewriteByLevel::updateGrapeRows(int row) {
  vector<int>& children = dag[row].second;
 
  vector<int> toBeErased;
  for(int i = 0 ; i < children.size(); i++) {
    int childLevel = levels[children[i]];

    // child is already rewritten
    if(childLevel < levels[row]+1) {
      // parent (i.e. row) is not a parent of this child anymore (due to rewriting) remove it from children list
      if(find(dag[children[i]].first.begin(), dag[children[i]].first.end(), row) == dag[children[i]].first.end()) {
        toBeErased.push_back(i);
      }
    } else {
      int maxLevel = findMaxLevelOfPredecessors(children[i]);

      // this is a grape to be shifted up
      int targetLevel = maxLevel+1;
      if(targetLevel < childLevel) {
        bool update = true;
        bool isALevelToBeRewritten = (find(levelsToBeRewritten.begin(), levelsToBeRewritten.end(), targetLevel) != levelsToBeRewritten.end());
        if(isALevelToBeRewritten) {
              
          while(find(levelsToBeRewritten.begin(), levelsToBeRewritten.end(), targetLevel) != levelsToBeRewritten.end() &&
                targetLevel > childLevel)
            targetLevel++;
              
          // already scheduled for rewriting, update
          if(toBeRewritten.find(children[i]) != toBeRewritten.end()) {
            // if the new level is a level to be rewritten, we cannot shift the row there
            // instead just shift it to the nearest level that wont be rewritten
            if((targetLevel == childLevel) || (targetLevel == toBeRewritten[children[i]].second)) {
              toBeRewritten.erase(children[i]);    
              update = false;
            } else if(targetLevel < toBeRewritten[children[i]].second) {
              toBeRewritten.erase(children[i]);    

            } else if(targetLevel > toBeRewritten[children[i]].second) {
              toBeRewritten[children[i]].first = targetLevel;
            }
          } else {
            if(targetLevel == childLevel)
              update = false;
          }
        } else {
          if(toBeRewritten.find(children[i]) != toBeRewritten.end())
            if(targetLevel <= toBeRewritten[children[i]].second)
              toBeRewritten.erase(children[i]);    
        }

        if(update) {
          updateLevelOf(children[i], targetLevel);
          updateGrapeRows(children[i]);
        }
      }
    }
  }

  for(int i = 0 ; i < toBeErased.size() ; i++)
    toBeErased.erase(toBeErased.begin() + i);
}

/*bool RewriteByLevel::isALevelToBeRewritten(int level) {
  return (find(levelsToBeRewritten.begin(), levelsToBeRewritten.end(), level) != levelsToBeRewritten.end());
}

vector<int>& RewriteByLevel::getLevelsToBeRewritten() {
  vector<int>& ref = levelsToBeRewritten;

  return ref;
}*/

