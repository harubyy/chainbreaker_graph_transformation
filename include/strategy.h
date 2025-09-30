#pragma once
#include <algorithm>
#include <vector>
#include <iterator>
#include <map>
#include <set>
#include <iostream>
#include <omp.h>
#include <shared_mutex>
#include <thread>
#include <chrono>
#include "analyzer.h"

#define REWRITE_UP    3
#define REWRITE_LOW   5
#define REWRITE_RATIO 0.05
//#define UPDATE_GRAPES  

// add random
enum class StartPoint {
  TopDown, BottomUp   };

enum class Scope   {
  Depth, Selective };

using namespace std;

// <row,<oldLevel,newLevel>>
typedef map<int,pair<int,int>> ToBeRewritten;
typedef map<int,set<int>> RewritingMap;

// indegree & outdegree vectors of a row: dependencies and children
typedef pair< vector<int>,vector<int> > Connectivity;
typedef vector<Connectivity> DAG;


// TODO: alignment
class RewritingStrategy {
  protected:
    int* levels;
    ToBeRewritten toBeRewritten;
    RewritingMap rewritingMap;
    int numOfRewrittenRows;
    vector< vector<int> >& levelTable;
    DAG& dag;

    int rows;
    StartPoint startPoint;
    Scope scope;


  public:
    RewritingStrategy(int rows, StartPoint startPoint, Analyzer* analyzer) : levelTable(analyzer->getLevelTable()), dag(analyzer->getDAG()) {
      this->levels = analyzer->getLevels();
      this->rows = rows;
      this->startPoint = startPoint;
      this->scope = Scope::Selective;
    }

    virtual void policy() = 0;
    ToBeRewritten& getToBeRewritten();
    RewritingMap& getRewritingMap();
    set<int>& getRewritingMapOf(int row);

    bool isRewritten(int row);
    bool isTopDown();
    bool isScopeSelective();

    void shiftRowsUp();
    void updateLevelOf(int row, int newLevel);
    void createRewritingMapFor(int row, int targetLevel);


    // duplicate function with Analyzer class, don't want to do it static
    int findMaxLevelOfPredecessors(int row);
    void expandPredsWith(int row, vector<int>& preds, int child);
    int findPredecessorsWithMaxLevel(vector<int>& preds, vector<int>& predsWithMaxLevel);

    // after rewriting strategy is applied, some levels will be emptied.
    // reflect the changes to related data structures (e.g. levelTable)
    // by calling Analyzer::updateWithEmptyLevels
    void findEmptyLevels(vector<int>& emptyLevels);

    // TODO: Fix this. It's wrong
    void updateGrapeRows(int row);

    void printRowsToBeRewritten() {
      cout << "num. of rewritten rows: " << toBeRewritten.size() << "\n";
      cout << "Rows to be rewritten:\n";
      for(auto& row : toBeRewritten)
        cout << "row : " << row.first << " at level " << row.second.first << " to " << row.second.second << "\n";
      cout << "\n";
    }

    void printRewritingMap() {
      cout << "Rewriting map for rows to be rewritten\n";
      for(auto& row : rewritingMap) {
        cout << "row: " << row.first << "\n";
        for(auto& item : row.second)
          cout << item << ",";
        cout << "\n";
      }

      cout << "\n\n";
    }

};

// RewriteLongerThan, RewriteShorterThan,
// RewriteInBetween, RewriteLowestIndegree
// should enable rows to be rewritten until (STOPPING CRITERIA)
//   -- they reach AIR, ARL, ALC or REWRITING_DIST (whichever is satisfied first)
// User defined parameters:
//   -- Length of the equation (a.k.a. # of parents or cost)
//   -- REWRITE_RATIO (can be selected from ALL, FLOPS_BELOW_AVG or FLOPS_ABOVE_AVG
// These methods should be invoked using appropriate indegree values (REWRITE_UP or REWRITE_LOW)
// TODO: add Random to StartPoint and select a level to scan or row
class RewriteLongerThan : public RewritingStrategy {
  public:
    RewriteLongerThan(int rows, StartPoint startPoint, Analyzer* analyzer) : 
                      RewritingStrategy(rows,startPoint, analyzer) {
      numOfRewrittenRows = rows * REWRITE_RATIO;
      cout << "intended num. of rewritten rows: " << numOfRewrittenRows << "\n";

      policy();
    }

    void policy();
};

class RewriteShorterThan : public RewritingStrategy {
  public:
    RewriteShorterThan(int rows, StartPoint startPoint, Analyzer* analyzer) : 
                      RewritingStrategy(rows,startPoint,analyzer) {
      numOfRewrittenRows = rows * REWRITE_RATIO;
      cout << "intended num. of rewritten rows: " << numOfRewrittenRows << "\n";

      policy();
    }

    void policy();
};

class RewriteInBetween : public RewritingStrategy {
  public:
    RewriteInBetween(int rows, StartPoint startPoint, Analyzer* analyzer) : 
                      RewritingStrategy(rows,startPoint,analyzer) {
      numOfRewrittenRows = rows * REWRITE_RATIO;
      cout << "intended num. of rewritten rows: " << numOfRewrittenRows << "\n";

      policy();
    }

    void policy();
};

// rewrite rows with the smallest # of in-degree at selectec levels
// TODO: implement
class RewriteLowestIndegree : public RewritingStrategy {
  public:
    RewriteLowestIndegree(int rows, StartPoint startPoint, Analyzer* analyzer) : 
                      RewritingStrategy(rows, startPoint, analyzer) {
      policy();
    }

    void policy();
};

class RewriteByCostMap : public RewritingStrategy {
  protected:
    float avgCostPerLevel;    // ALC
    int REWRITING_DIST = 20;

    vector<int> levelCost;
    map<int,int> flopsBelowAvg;


  public:
    RewriteByCostMap(int rows, StartPoint startPoint, Analyzer* analyzer) :
                      RewritingStrategy(rows, startPoint, analyzer) {
//     auto t1 = std::chrono::steady_clock::now();
       // copy the contents so that we can work on it
       levelCost = analyzer->getFlopsPerLevel();

       #ifdef REWRITE_ENABLED
         analyzer->separateThinLevels();
         analyzer->analyzeForCriteriaCostMap();
//         analyzer->analyzeForCriteriaCoeff();
//         analyzer->printCriteria();
         analyzer->saveOriginalValues();
       #endif

       this->avgCostPerLevel = analyzer->getALC();
       this->flopsBelowAvg = analyzer->getFlopsBelowAvg();

       auto start = chrono::high_resolution_clock::now();
       #ifdef PAR
         distPolicy();
         auto end = chrono::high_resolution_clock::now();
         cout << fixed <<  "* chrono policy: " << chrono::duration<double>(end-start).count() * 1000    << " :seconds\n";
       #else
         policy();
         auto end = chrono::high_resolution_clock::now();
         cout << fixed <<  "chrono policy: " << chrono::duration<double>(end-start).count() * 1000    << " :seconds\n";
       #endif


       analyzer->setFlopsPerLevel(levelCost);
    }
    
    // parallel exec.
    void distPolicy();
    void policy(int start, int end, ToBeRewritten& rewritingList);
    void createCostMapMax(map<int,int>::iterator& sourceLevel, map<int,int>::iterator& targetLevel, vector<int>::iterator& nextRow, ToBeRewritten& rewritingList);

    // serial exec.
    void policy();
    void createCostMapMax(map<int,int>::iterator& sourceLevel, map<int,int>::iterator& targetLevel, vector<int>::iterator& nextRow);

    // version of expandPredsWith that does not alter the dag (children)
    void expandPredsWith(int row, vector<int>& preds, int child);
};

class RewriteByThreeCriteria : public RewritingStrategy {
  protected:
    int REWRITING_DIST = 20;

    vector<int> levelCost;
    map<int,int> levelSizeBelowAvg;
    map<int,int> flopsBelowAvg;
    vector<int> numOfLevels;
//    #define FAILED_ANALYSIS ON
    #ifdef FAILED_ANALYSIS
      map<int,tuple<int,int>> failedRowsIndegree; // AIR: row level, row indegree
      map<int,tuple<int,int>> failedRowsCost;     // ALC: row level, row cost
      vector<int> failedRowsSize;                 // ARL: curr level, level size
    #endif

    float avgCostPerLevel, avgIndegreePerRow, avgLevelSize; // ALC, AIR, ARL

    public:
      RewriteByThreeCriteria(int rows, StartPoint startPoint, Analyzer* analyzer) :
                      RewritingStrategy(rows, startPoint, analyzer) {
          levelCost = analyzer->getFlopsPerLevel();
          #ifdef REWRITE_ENABLED
            analyzer->separateThinLevels();
            analyzer->analyzeForCriteria2();
            //analyzer->analyzeForCriteriaCoeff();
            analyzer->printCriteria();
            analyzer->saveOriginalValues();
          #endif
          this->flopsBelowAvg = analyzer->getFlopsBelowAvg();

          this->avgCostPerLevel = analyzer->getALC();
          this->avgIndegreePerRow = analyzer->getAIR();
          this->avgLevelSize = analyzer->getARL();

          // get # of rows for levels with flopsBelowAvg
          for(auto& level: flopsBelowAvg)
            levelSizeBelowAvg[level.first] = levelTable[level.first].size();
          auto start = chrono::high_resolution_clock::now();
          #ifdef PAR
            distPolicy();
            auto end = chrono::high_resolution_clock::now();
            cout << fixed <<  "* chrono policy: " << chrono::duration<double>(end-start).count() * 1000 << " :seconds\n";
          #else
            policy();
            auto end = chrono::high_resolution_clock::now();
            cout << fixed <<  "chrono policy: " << chrono::duration<double>(end-start).count() * 1000 << " :seconds\n";
          #endif
          

          analyzer->setFlopsPerLevel(levelCost);
          #ifdef FAILED_ANALYSIS
            cout << "failed AIR:," << analyzer->getAIR() << "\n";
            cout << "failedRowsIndegree: " << failedRowsIndegree.size() << "\n";
/*            cout << "indegrees:\n";
            for(auto& row : failedRowsIndegree)
              cout << get<1>(row.second) << "\n";*/
          
            cout << "failed ALC:," << analyzer->getALC() << "\n";
            cout << "failedRowsCost: " << failedRowsCost.size() << "\n";
/*            cout << "costs:\n";
            for(auto& row : failedRowsCost)
              cout << get<1>(row.second) << "\n";*/
          
            cout << "failedRowsSize: " << failedRowsSize.size() << "\n\n";
          #endif
      }

      // parallel exec.
      void distPolicy();
      void policy(int start, int end, ToBeRewritten& rewritingList);

      // serial exec.
      void policy();

      // version of expandPredsWith that does not alter the dag (children)
      void expandPredsWith(int row, vector<int>& preds, int child);
};

// rewrite rows on critical path (CP)
// TODO: implement
class RewriteOnCP : public RewritingStrategy {
  protected:
    vector< vector<int> > paths;
    shared_mutex mutex;

  public:
    RewriteOnCP(int rows, StartPoint startPoint, Analyzer* analyzer) : 
                      RewritingStrategy(rows, startPoint, analyzer) {
      policy();
    }

    void policy();
    void traversePath(vector<int>& path, int row);
    void buildCriticalPath(vector<int>& lastLevelRows);
};

// MANUAL(HANDCRAFTED) REWRITING METHOD
// methods using RewriteByLevel
// must traverse rowsToBeRewritten in REVERSE order
// hence RewriteByLevel is BottomUp by nature
class RewriteByLevel : public RewritingStrategy {
  protected:
    map<int,vector<int>> levelLengthHistogram;
    vector<int> levelsToBeRewritten;
    vector<int> flopsPerLevelRewrite;

  public:
    RewriteByLevel(int rows, StartPoint startPoint, Analyzer* analyzer) : 
                      RewritingStrategy(rows, startPoint, analyzer) {
//     prepare(flopsBelowAvg);    
//     policy();

//       lung2_9levels_to10th();
     torso2_9levels_to10th(analyzer->getFlopsBelowAvg());

     //  cout << "num. of rewritten rows: " << toBeRewritten.size() << "\n";
//       printRowsToBeRewritten();
    }
    
    void printLevelsToBeRewritten() {
      cout << "levelsToBeRewritten:\n";
      for(auto& level : levelsToBeRewritten)
        cout << level << ", ";
      cout << "\n";
    }

    // currently we're merging all levels with the min. # of rows for simplicity
    void policy();
    void prepare(map<int,int>& flopsBelowAvg);

    // this is for calculating for the manual approach (rewrite by level) when
    // calculation in rewrite module is unfunctional
    void calculateFLOPSAfterRewrite(vector<int>& flopsPerLevelRewrite);
    void torso2_9levels_to10th(map<int,int>& flopsBelowAvg);
    void lung2_9levels_to10th();
    void lung2_7levels_to8th();

    void updateGrapeRows(int row);

    //bool isALevelToBeRewritten(int level);
    //vector<int>& getLevelsToBeRewritten();
};
