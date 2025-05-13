#include <iostream>
#include <libufget.h>
#include <fstream>
#include <memory>
#include <chrono>
#include "matrix.h"
#include "rewrite.h"
#include "analyzer.h"
#ifdef REWRITE_ENABLED
  #include "strategy.h"
#endif
using namespace std;

int main(int argc, char *argv[]) {
  uf_collection_t *collection;
  int matID;

  collection = uf_collection_init();

  if (argc != 2) {
    cout << "Usage: programName matID\n";
    return 0;
  }

  matID = std::stoi(std::string(argv[1]));
  printf("UF Sparse Matrix Collection Database contains %d matrices.\n", uf_collection_num_matrices(collection));

  Matrix* matrixCSR = new Matrix(collection, matID);
  int successCSR = matrixCSR->convertToCSR();
  printf("success (0) for COO --> CSR:%d\n", successCSR);

  Matrix* matrixCSC = new Matrix(collection, matID);
  int successCSC = matrixCSC->convertToCSC();
  printf("success (0) for COO --> CSC:%d\n", successCSC);

  matrixCSR->extractLCSR();
//  cout << "################################################################ Printing after extractLCSR\n";
//  matrixCSR->print();

  matrixCSC->extractLCSC();

  int cols = matrixCSC->getL()->getCols();
  Analyzer* analyzer = new Analyzer(matrixCSR, matrixCSC, cols);
 auto t1 = std::chrono::steady_clock::now();
  analyzer->buildLevels();
  analyzer->calculateFLOPS();
 auto t2 = std::chrono::steady_clock::now();
 cout << "chrono build DAG: " << chrono::duration<double>(t2 - t1).count()*1000 << "\n";

  analyzer->report(string ("BEFORE"), 0);

  #ifdef REWRITE_ENABLED
 auto t3 = std::chrono::steady_clock::now();
//    RewriteByCostMap* rewritingStrategy = new RewriteByCostMap(matrixCSR->getL()->getRows() - 1, StartPoint::BottomUp, analyzer);
    RewriteByThreeCriteria* rewritingStrategy = new RewriteByThreeCriteria(matrixCSR->getL()->getRows() - 1, StartPoint::BottomUp, analyzer);

    rewritingStrategy->shiftRowsUp();
    vector<int> emptyLevels;
    rewritingStrategy->findEmptyLevels(emptyLevels);
    analyzer->updateWithEmptyLevels(emptyLevels);
 auto t4 = std::chrono::steady_clock::now();
 cout << "chrono transform DAG: " << chrono::duration<double>(t4 - t3).count()*1000 << "\n";

    rewritingStrategy->printRowsToBeRewritten();

 t1 = std::chrono::steady_clock::now();
    Rewrite rewriter(matrixCSR, matrixCSC, std::string(matrixCSR->getUF_matrix()->name), rewritingStrategy, analyzer);
  #else
 t1 = std::chrono::steady_clock::now();
    Rewrite rewriter(matrixCSR, matrixCSC, std::string(matrixCSR->getUF_matrix()->name), analyzer);
  #endif

  rewriter.rewrite();
 t2 = std::chrono::steady_clock::now();
 cout << "chrono rewrite module: " << chrono::duration<double>(t2 - t1).count()*1000 << "\n";
    
  analyzer->report(string ("AFTER"), 0);

  return 0;
}
