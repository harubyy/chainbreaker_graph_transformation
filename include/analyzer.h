#pragma once
#include <iostream>
#include <vector>
#include <map>
#include <utility>
#include <set>
#include "matrix.h"

#define UNROLL_FACTOR 5

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
    map<int, vector<int>> oriParents; // original parents of all rows
    map<int, vector<double>> oriRowValues; // original row values of all rows
                                           
    vector< vector<double>> values;
    int numOfLevels;
    bool singleLoopRows;  // if at least 2 rows exist with nnzs > 10 in a level, loops will be merged, CSR loop will be used
    // Criteria for rewriting
    // 2*nnzs-rows FLOPS in total
    int totalFLOPSPerLevel;
    vector<int> flopsPerLevel;
    map<int,int> flopsBelowAvg;
    map<int,int> flopsAboveAvg;
    
    // ALC: avg. Level Cost
    // AIR: avg. indegrees(parents/dependencies) per row
    // ARL: avg. # of rows per level for flopsBelowAvg
    // MMAD: max. memory access distance of a row
    // AMAD: avg. memory access distance of a row
    // MID: max. indegree distance (I'm curious)
    float ALC, AIR, ARL;
    float ALC_CV, AIR_CV, ARL_CV;
    int MMAD, MID;
    float AMAD;


    #ifdef REWRITE_ENABLED
      // levels to avoid altering levelTable by rewriting while traversing it.
      vector<int> flopsPerLevelRewrite;
    #endif

    Matrix* matrixCSR;
    Matrix* matrixCSC;


  public:
   Analyzer(Matrix* matrixCSR, Matrix* matrixCSC, int cols)
    : matrixCSR(matrixCSR), matrixCSC(matrixCSC), dag(cols - 1), values(cols - 1) {}

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
    map<int,vector<int>>& getOriParents();
    map<int,vector<double>>& getOriRowValues();
    int getNumOfLevels();
    bool getSingleLoopRows();
    void setSingleLoopRows(bool looped);
    vector< vector<int>>& getLevelTable();
    vector< vector<double>>& getValues();
    vector<int>& getFlopsPerLevel();
    void setFlopsPerLevel(vector<int>& levelCost);
    map<int,int>& getFlopsBelowAvg();
    map<int,int>& getFlopsAboveAvg();

    float getALC();
    float getAIR();
    float getARL();
    int getMMAD();
    int getMID();
    float getAMAD();

    // compute AIR, ARL, MMAD, MID, AMAD (ALC is calculated in calculateFLOPS)
    void analyzeForCriteria();
    void analyzeForCriteriaCostMap();
    void analyzeForCriteriaCoeff();

    void separateRows(vector<int>& loopedRows, vector<int>& unrolledRows);
    void validityCheck();

    #ifdef REWRITE_ENABLED
      typedef map<int,pair<int,int>> ToBeRewritten;

      vector<int>& getFlopsPerLevelRewrite();
      int findMaxLevelOfPredecessors(int row);
    #endif

    void updateWithEmptyLevels(vector<int>& emptyLevels);
    void separateThinLevels();
    void buildLevels();
    void calculateFLOPS();
    int recalculateFLOPSFor(vector<int> loopedRows);
    #ifdef REWRITE_ENABLED
      void saveOriginalValues();
    #endif

    void printLevels();
    void printLevelTable();
    void printLevelSizes();
    void printDAG();
    void printValues();
    void printFLOPSPerLevel();
    void printFLOPSDivided();
    void printCriteria();
    void report(string reportType, int reportStep);
    void printDependencies();
};
