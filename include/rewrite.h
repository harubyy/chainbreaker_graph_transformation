#pragma once
#include <sys/time.h>
#include <cstring>
#include <vector>
#include <math.h>
#include <sstream>
#include <omp.h>
#include <sys/stat.h>
#include "matrix.h"
#include "analyzer.h"
#ifdef REWRITE_ENABLED
  #include "strategy.h"
#endif

using namespace std;

//#define PART_SIZE 30  // typically will be the size of a function executed by a thread
#define NUM_THREADS 12
#define LOWER_BOUND 10

// choose load balancing approach among threads for computing a level
#define BALANCE_ROWS  1
//#define BALANCE_ROWS_CV 1

// size of the table keeping track of rewriting in rewriteExecutor
#define TABLE_SIZE 200

// TODO: alignment for class variables
class Rewrite {
 private:
  Matrix* matrixCSR;
  Matrix* matrixCSC;

  // TODO: this shouldnt be here
  bool dumpCSR;

  #ifdef REWRITE_ENABLED
    map<int,vector<double>> rewrittenRowValues;
    map<int,vector<int>> rewrittenRowParents;
    map<int, double> rewrittenB;

    RewritingStrategy* rewritingStrategy;
  #endif

  Analyzer* analyzer;
  std::string const &fileName;

  vector<int> rowIndices;
  vector< vector<int> > startIndex;

  vector<int> signatureLevel;      // keeps whether level holds only unrolled AND/OR rewritten rows or looped rows also (# of levels)
  vector< vector<int> > signature; // keeps whether part holds only unrolled AND/OR rewritten rows or looped rows also  (# of levels x # of parts per row)

 public:
  #ifdef REWRITE_ENABLED
    Rewrite(Matrix* matrixCSR, Matrix* matrixCSC, std::string const &fileName, RewritingStrategy* rewritingStrategy, Analyzer* analyzer)
       : matrixCSR(matrixCSR), matrixCSC(matrixCSC), fileName(fileName), rewritingStrategy(rewritingStrategy), analyzer(analyzer) {
         if (mkdir(fileName.c_str(), 0777) == -1)
           cerr << "Error :  " << strerror(errno) << endl;

         //dumpCSR = true;
         dumpCSR = false;
    }
  #else
    Rewrite(Matrix* matrixCSR, Matrix* matrixCSC, std::string const &fileName, Analyzer* analyzer)
        : matrixCSR(matrixCSR), matrixCSC(matrixCSC), analyzer(analyzer), fileName(fileName) {
          if (mkdir(fileName.c_str(), 0777) == -1)
            cerr << "Error :  " << strerror(errno) << endl;
        }
  #endif

//  void allocateMemory(std::ostream &stream, int segment);
  int dumpDataToMem(vector<double>& b);
  void writeMain(std::ostream &stream, int segment, int parentsSize);
 // void writeFunc(int level, int runnablesStart, int numOfRunnables);
  void writeFunc(std::ofstream &stream, vector< vector<int> >& tracker, int size, int part, int maxNumOfThreads, int headerCounter);
  int writeMakefile(int execBlockCnt);
  int writeHeader(int partCounter, int partCounterStart, int partCounterEnd);
  void writeUtil();

  int writePart(int levelNum, int rowStartIndex, int rowEndIndex, int levelPart, vector<double>& b, std::ofstream& stream, int merged, int* sigType);

  int balanceLevel(int toBeBalanced, int& workloadPerThread);
  #ifdef BALANCE_FLOPS
    int balanceRows(vector<int>& level);
  #elif BALANCE_ROWS_CV
    void collectFLOPS(vector<int>& level, vector<int>& flopsLevel);
    double calculateCV(int numOfSum, int start, vector<int>& level);
    void advanceTillBinLimit(vector<int>& level, int start, int& end, int binLimit);
    void binLastElement(vector<int>& level, int& start, int& end, vector<int>& rowDist, vector<int>& rowDistSum);
    double balanceRowsCV(vector<int>& level, vector<int>& rowDist, vector<int>& rowDistSum, double totalCV);
    void reduceNumOfRows(vector<int>& level, vector<int>& rowDistReduced, vector<int>& rowDistSumReduced);
  #endif

  #ifdef REWRITE_ENABLED
    int rewriteRow(int row, vector<double>& b, string& currRow, int rewriteDepth);
    int rewriteRow(int row, vector<double>& b, string& currRow, set<int>& rewritten);
    void rewriteRow2(int row, vector<double>& b, set<int>& rewritten, map<int,double>& multipliers, double rewritingMultiplicant, bool sign);
    int rewriteLevel(int level, vector<double>& b, string& currRow, std::ostream &stream);
    int dumpString(int row, stringstream& currRow, set<int>& rewritten, map<int,double>& multipliers);

  // for dumping transformed matrix into CSR format
  void buildTransformedMatrix(vector<int>& transformedRowPtr, vector<int>& transformedColIdx, vector<double>& transformedValues);
  #endif

  int rewriteExecutor(vector<double>& b, vector<double> &x);
  int rewrite();


};
