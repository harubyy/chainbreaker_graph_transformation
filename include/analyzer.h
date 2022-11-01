#pragma once
#include <iostream>
#include <vector>
#include <map>
#include <utility>
#include <set>
#include "matrix.h"

#ifdef REWRITE_ENABLED
// TODO: make user defined
  // should be tuned according to DOP or some other metric(s)
  // REWRITE_UP >= REWRITE_DEPTH
//  #define REWRITE_UP      3       // upper bound to start rewriting
  #define REWRITE_DEPTH   0
  #define REWRITE_PERCENT 5
#endif


using namespace std;

// indegree & outdegree vectors of a row: dependencies and children
typedef pair< vector<int>,vector<int> > Connectivity;
typedef vector<Connectivity> DAG;

// TODO: alignment
// TODO: delete matrices CSR & CSC
class Analyzer {
  private:
    int* levels;
    vector< vector<int> > levelTable;
    DAG dag;
    vector< vector<double>> values;
    int numOfLevels;
    bool singleLoopRows;  // if at least 2 rows exist with nnzs > 10 in a level, loops will be merged, CSR loop will be used
    // Criteria for rewriting
    // 2*nnzs-rows FLOPS in total
    double totalFLOPSPerLevel;
    vector<int> flopsPerLevel;
    map<int,double> flopsBelowAvg;
    map<int,double> flopsAboveAvg;
    
    // ALC: avg. Level Cost
    // AIR: avg. indegrees(parents/dependencies) per row
    // ARL: avg. # of rows per level for flopsBelowAvg
    // MMAD: max. memory access distance of a row
    // MID: max. indegree distance (I'm curious)
    float ALC, AIR, ARL;
    int MMAD, MID;


    #ifdef REWRITE_ENABLED
      // levels to avoid altering levelTable by rewriting while traversing it.
      vector<int> flopsPerLevelRewrite;
    #endif

    Matrix* matrixCSR;
    Matrix* matrixCSC;


  public:
   Analyzer(Matrix* matrixCSR, Matrix* matrixCSC)
    : matrixCSR(matrixCSR), matrixCSC(matrixCSC) {}

   // TODO: any clearing to data structures needed?
    ~Analyzer() {
      matrixCSR = nullptr;
      matrixCSC = nullptr;

      if(levels != nullptr) {
        free(levels);
        levels = nullptr;
      }
    }

    int* getLevels();
    DAG& getDAG();
    int getNumOfLevels();
    bool getSingleLoopRows();
    vector< vector<int>>& getLevelTable();
    vector< vector<double>>& getValues();
    vector<int>& getFlopsPerLevel();
    void setFlopsPerLevel(vector<int>& levelCost);
    map<int,double>& getFlopsBelowAvg();
    map<int,double>& getFlopsAboveAvg();

    float getALC();
    float getAIR();
    float getARL();
    int getMMAD();
    int getMID();

    // compute AIR & ARL & MMAD & MID (ALC is calculated in calculateFLOPS)
    void analyzeForCriteria();

    void separateRows(int levelNum, int rowStartIndex, int rowEndIndex, vector<int>& loopedRows, vector<int>& unrolledRows);

    #ifdef REWRITE_ENABLED
      typedef map<int,pair<int,int>> ToBeRewritten;

      vector<int>& getFlopsPerLevelRewrite();
      int findMaxLevelOfPredecessors(int row);
      void correctAfterRewritingStrategy(vector<int>& emptyLevels);
    #endif

    void buildLevels();
    void calculateFLOPS();
    #ifdef REWRITE_ENABLED
      void calculateLevelsToBeRewritten();
    #endif

    void printLevels();
    void printLevelTable();
    void printLevelSizes();
    void printDAG();
    void printValues();
    void printFLOPSPerLevel();
    void report(string reportStep);
    #ifdef REWRITE_ENABLED
      void printDependencies();
    #endif
};
