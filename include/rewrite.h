#pragma once
#include <sys/time.h>
#include <cstring>
#include <vector>
#include <list>
#include <math.h>
#include <sstream>
#include <omp.h>
#include <sys/stat.h>
#include <chrono>
#include "matrix.h"
#include "analyzer.h"
#include "par.h"
#ifdef REWRITE_ENABLED
  #include "strategy.h"
#endif

using namespace std;

//#define PART_SIZE 30  // typically will be the size of a function executed by a thread

//#define NUM_THREADS 4
//#define LOWER_BOUND 40

// choose load balancing approach among threads for computing a level
#define BALANCE_ROWS  1
//#define BALANCE_ROWS_CV 1

// size of the table keeping track of rewriting in rewriteExecutor
#define TABLE_SIZE 200

// merge max of 10 consecutive single-threaded levels
// experimental, commented out in the code: rewrite.cpp
#define MERGE_DEPTH 10

// some stats are reported
//#define REPORT 1

// TODO: alignment for class variables
class Rewrite {
 private:
  Matrix* matrixCSR;
  Matrix* matrixCSC;

  Analyzer* analyzer;

  #ifdef REWRITE_ENABLED
    RewritingStrategy* rewritingStrategy;
  #endif

  vector<int> rowIndices;

  // keeps whether level holds only unrolled AND/OR rewritten rows or looped rows also (# of levels)
  vector<int> signatureLevel;
  // keeps whether part holds only unrolled AND/OR rewritten rows or looped rows also  (# of levels x # of parts per row)
  vector< vector<int> > signature; 
  // key: starting level, vector: total cost, signature, ending level for merged levels
  map<int, vector<int>> mergedLevels;
  // mergedLevels start and ending levels will be keys to mergedLevelsDist
  map<int, list<vector<int>>> mergedLevelsDist; // key: level, vector: looped rows, unrolled rows
  // when empty levels are erased merged level orders change. key: new lvl start, val: ori. lvl start
  map<int, int> levelLookUp;
  std::string const &fileName;
    
  chrono::time_point<chrono::steady_clock> ch_start,ch_end;
  chrono::time_point<chrono::steady_clock> ch_start2,ch_end2;
  chrono::duration<double> ch_ttime, ch_ttime2;

 public:
  #ifdef REWRITE_ENABLED
    Rewrite(Matrix* matrixCSR, Matrix* matrixCSC, std::string const &fileName, RewritingStrategy* rewritingStrategy, Analyzer* analyzer)
       : matrixCSR(matrixCSR), matrixCSC(matrixCSC), fileName(fileName), rewritingStrategy(rewritingStrategy), analyzer(analyzer) {
         if (mkdir(fileName.c_str(), 0777) == -1)
           cerr << "Error :  " << strerror(errno) << endl;
    }
  #else
    Rewrite(Matrix* matrixCSR, Matrix* matrixCSC, std::string const &fileName, Analyzer* analyzer)
        : matrixCSR(matrixCSR), matrixCSC(matrixCSC), analyzer(analyzer), fileName(fileName) {
          if (mkdir(fileName.c_str(), 0777) == -1)
            cerr << "Error :  " << strerror(errno) << endl;
        }
  #endif

  int dumpDataToMem(vector<double>& b);
  void writeMain(std::ostream &stream, int segment, int parentsSize);
  void writeFunc(std::ofstream &stream, vector< vector<int> >& tracker, int size, int maxNumOfThreads, int headerCounter);
  int writeMakefile(int execBlockCnt);
  int writeHeader(int partCounter, int partCounterStart, int partCounterEnd, vector<int>& threadCounts);
  void writeUtil();

  int writePart(int levelNum, int rowStartIndex, int rowEndIndex, int levelPart, vector<double>& b, std::ofstream& stream, int merged, int* sigType);

  int balanceLevel(int toBeBalanced, int& workloadPerThread);
  void mergeSingleThreadedLevels();
  void dumpLoop(int beginIndex, int endIndex, vector<int>& loopedRows, std::ofstream& stream);
  int dumpUnrolled(vector<int>& unrolledRows, vector<double> &b, std::ofstream& stream);
  #ifdef BALANCE_FLOPS
    int balanceRows(vector<int>& level);
  #elif BALANCE_ROWS_CV
    void collectFLOPS(vector<int>& level, vector<int>& flopsLevel);
    double calculateCV(int numOfSum, int start, vector<int>& level);
    void advanceTillBinLimit(vector<int>& level, int start, int& end, int binLimit);
    void binLastElement(vector<int>& level, int& start, int& end, vector<int>& rowDist, vector<int>& rowDistSum);
    double balanceRowsCV(vector<int>& level, vector<int>& rowDist, vector<int>& rowDistSum, double totalCV);
  #endif

  #ifdef REWRITE_ENABLED
    void rewriteInLoop(vector<int>& loopedRows, vector<double> &b);
    void calculateMultiplicants(int row, vector<double>& b, set<int>& rewritten, map<int,double>& multipliers, double rewritingMultiplicant, bool sign);
    void updateRewrittenRow(int row, vector<double>& b, set<int>& rewritten, map<int,double>& multipliers);
    void rewriteRow(int row, vector<double>& b);

  // for dumping transformed matrix into CSR format
  void buildTransformedMatrix(vector<int>& transformedRowPtr, vector<int>& transformedColIdx, vector<double>& transformedValues);
  #endif

  int rewriteExecutor(vector<double>& b);
  int rewrite();
};
