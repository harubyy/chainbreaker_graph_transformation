#include <ios>
#include <sys/time.h>
#include <cstring>
#include <vector>
#include <math.h>
#include <omp.h>
#include <sys/stat.h>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include "rewrite.h"
#include "util.h"

using namespace std;

// writes runX.c
void Rewrite::writeFunc(std::ofstream &stream, vector< vector<int> >& tracker, int size, int part, int maxNumOfThreads, int headerCounter) {
  if(!stream.is_open())
    std::cout << "Cannot open output file!\n";

  vector<int>& partStarts   = tracker[0];
  vector<int>& threadCounts = tracker[1];

/*  cout << "size: " << size << "\n";
  cout << "tracker table:\n";
  for(int i = 0; i < size ; i++)
    cout << i << ": " << partStarts[i] << ", " << threadCounts[i] << "\n";
  cout << "\n"; */

  stream << "#include \"calculators" + to_string(headerCounter) + ".h\"\n\n";

  int sign;
  if((part+1) * TABLE_SIZE > signatureLevel.size())
    sign = accumulate(signatureLevel.begin() + part * TABLE_SIZE, signatureLevel.end(), 0);
  else
    sign = accumulate(signatureLevel.begin() + part * TABLE_SIZE, signatureLevel.begin() + ((part+1) * TABLE_SIZE), 0);
  if(sign)
    stream << "void run" << part << "(double* x, double* b, int* parents, double* values, int* rowPtr, int* rowIndices) {\n";
  else
    stream << "void run" << part << "(double* x) {\n";

  int counter = partStarts[0];
  if (maxNumOfThreads > 1)
    stream << "\n#pragma omp parallel num_threads(" << maxNumOfThreads << ")\n{\n";

  int levelCounter = part * TABLE_SIZE;
  int threadCounter = 0;
  int mergeMark = 0, mergeBegin = 0, mergeCurrent = 0;
//  cout << "size: " << size << "\n";
  for(int j = 0 ; j < size ; j++) {
    if(threadCounts[j] > 1) {
      stream << "\n";
      stream << "  #pragma omp sections // " << j << ", " << to_string(threadCounts[j]) << "\n"
             << "  {\n";

 //     cout << "counter: " << counter << " threadCounts[j]: " << threadCounts[j] << " signature size: " << signature[levelCounter].size() << "\n";
      for(int i = counter ; i < counter+threadCounts[j] ; i++) {
        stream << "    #pragma omp section\n";
  //      cout << "levelCounter: " << levelCounter << " threadCounter: " << threadCounter << "\n";
        if(signature[levelCounter][threadCounter] == 0)
          stream << "    { calculate" << i << "(x); }\n";
        else
          stream << "    { calculate" << i << "(x, b, parents, values, rowPtr, rowIndices); }\n";

        threadCounter++;
      }
      stream << "  }\n";
    } else {
       if(!mergeBegin) {
         mergeBegin = 1;
         mergeMark = j;
         mergeCurrent = j;

         if(signature[levelCounter][0] == 0) {
            stream << "  #pragma omp single\n"
                   << "  calculate" << counter << "(x); \n";
          } else {
            stream << "  #pragma omp single\n"
                   << "  calculate" << counter << "(x, b, parents, values, rowPtr, rowIndices); \n";
          }

         stream << "\n";
       }

       if(mergeBegin == 1 && j == (mergeCurrent + 1))
         mergeCurrent = j; 

        // toggle the flag for the last merged
       if((mergeCurrent != j) || (mergeCurrent - mergeMark) == 9)
         mergeBegin = 0;
    }

    counter += threadCounts[j];
    levelCounter++;
    threadCounter = 0;
  }

  if (maxNumOfThreads > 1)
    stream << "  }\n";

  stream << "}\n";
    
  stream.close();
}

//void Rewrite::allocateMemory(std::ostream &stream, int numOfParts) {
void Rewrite::writeMain(std::ostream &stream, int numOfParts, int parentsSize) {
  Part* L = matrixCSR->getL();
  int rows = L->getRows();
  int vals = L->getNNZs();

  stream << "#include <stdlib.h>\n"
         << "#include <omp.h>\n"
         << "#include <math.h>\n"
         << "#include <stdio.h>\n"
         << "#include \"util.h\"\n\n";

  stream << "#define FILEPATH \"/tmp/" + fileName + ".bin\"\n"
         << "#define X_SIZE " +  to_string(rows-1) + "\n";

 #ifdef REWRITE_ENABLED
   if(analyzer->getSingleLoopRows()) {
     stream << "#define FILEPATH_B \"/tmp/" + fileName + "_b_TR.bin\"\n"
            << "#define FILEPATH_PARENTS \"/tmp/" + fileName + "_parents_TR.bin\"\n"
            << "#define FILEPATH_ROWPTR \"/tmp/" + fileName + "_rowPtr_TR.bin\"\n"
            << "#define FILEPATH_ROWINDICES \"/tmp/" + fileName + "_rowIndices_TR.bin\"\n"
            << "#define FILEPATH_VALUES \"/tmp/" + fileName + "_vals_TR.bin\"\n"
            << "#define FILESIZE_P (PARENTS_SIZE * sizeof(int))\n"
            << "#define FILESIZE_ROWPTR (X_SIZE * sizeof(int))\n"
            << "#define FILESIZE_ROWINDICES (" + to_string(rowIndices.size()) + " * sizeof(int))\n"
            << "#define FILESIZE_VAL (PARENTS_SIZE * sizeof(double))\n"
            << "#define PARENTS_SIZE " +  to_string(parentsSize) + "\n\n"; 
   } else {
     // just produce CSR format arrays for SYCL code
     if(dumpCSR) {
       stream << "#define FILEPATH_PARENTS \"/tmp/" + fileName + "_parents_TR.bin\"\n"
              << "#define FILEPATH_ROWPTR \"/tmp/" + fileName + "_rowPtr_TR.bin\"\n"
              << "#define FILEPATH_VALUES \"/tmp/" + fileName + "_vals_TR.bin\"\n"
              << "#define PARENTS_SIZE " +  to_string(parentsSize) + "\n\n"; 
     }
   }
 #else
   if(analyzer->getSingleLoopRows()) {
     stream << "#define FILEPATH_B \"/tmp/" + fileName + "_b.bin\"\n"
            << "#define FILEPATH_PARENTS \"/tmp/" + fileName + "_parents.bin\"\n"
            << "#define FILEPATH_ROWPTR \"/tmp/" + fileName + "_rowPtr.bin\"\n"
            << "#define FILEPATH_ROWINDICES \"/tmp/" + fileName + "_rowIndices.bin\"\n"
            << "#define FILEPATH_VALUES \"/tmp/" + fileName + "_vals.bin\"\n"
            << "#define FILESIZE_ROWINDICES (" + to_string(rowIndices.size()) + " * sizeof(int))\n"
            << "#define PARENTS_SIZE " +  to_string(vals) + "\n\n"; 
   } else {
     // just produce CSR format arrays for SYCL code
     if(dumpCSR) {
       stream << "#define FILEPATH_PARENTS \"/tmp/" + fileName + "_parents.bin\"\n"
              << "#define FILEPATH_ROWPTR \"/tmp/" + fileName + "_rowPtr.bin\"\n"
              << "#define FILEPATH_VALUES \"/tmp/" + fileName + "_vals.bin\"\n"
              << "#define PARENTS_SIZE " +  to_string(vals) + "\n\n"; 
     }
   }
 #endif

 stream << "#include <unistd.h>\n"
        << "#include <time.h>\n"
        << "\n"
        << "struct timespec t1, t2;\n"
        << "\n"
        << "#define TIME(...) fprintf(stdout, __VA_ARGS__)\n"
        << "#define START() clock_gettime(CLOCK_MONOTONIC, &t1);\n"
        << "#define STOP(event) \\\n"
        << " do {  \\\n"
        << "     double timeSpent = 0.0; \\\n"
        << "     struct timespec temp; \\\n"
        << "     clock_gettime(CLOCK_MONOTONIC, &t2); \\\n"
        << "     if (((&t2)->tv_nsec-(&t1)->tv_nsec)<0) { \\\n"
        << "         temp.tv_sec = (&t2)->tv_sec-(&t1)->tv_sec-1; \\\n"
        << "         temp.tv_nsec = 1000000000+(&t2)->tv_nsec-(&t1)->tv_nsec; \\\n"
        << "     } else { \\\n"
        << "         temp.tv_sec = (&t2)->tv_sec-(&t1)->tv_sec; \\\n"
        << "         temp.tv_nsec = (&t2)->tv_nsec-(&t1)->tv_nsec; \\\n"
        << "     } \\\n"
        << "     \\\n"
        << "     timeSpent = (temp.tv_nsec) / 1e6; \\\n"
        << "     timeSpent = timeSpent + temp.tv_sec*1000; \\\n"
        << "     TIME(\"%s: %4.2f ms\\n\",event,timeSpent); \\\n"
        << "    } while (0)\n\n";


  for(int i = 0; i < numOfParts; i++) {
    int sign;
    if((i+1) * TABLE_SIZE > signatureLevel.size())
      sign = accumulate(signatureLevel.begin() + i * TABLE_SIZE, signatureLevel.end(), 0);
    else
      sign = accumulate(signatureLevel.begin() + i * TABLE_SIZE, signatureLevel.begin() + ((i+1) * TABLE_SIZE), 0);

    if(sign == 0)
      stream << "void  run" << i << "(double* x);\n"; 
    else
      stream << "void run" << i << "(double* x, double* b, int* parents, double* values, int* rowPtr, int* rowIndices);\n";
  }

  stream << "\n";

  stream << "int main() {\n"
         << "  int fd;\n"                 // holds x (x_r & x_w)
         << "  int result;\n\n"
         << "  double* x;\n"
         << "  fd = open(FILEPATH, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);\n"
         << "  check_err_open_lseek(fd);\n"
         << "  result = lseek(fd, X_SIZE * sizeof(double) - 1, SEEK_SET);\n"
         << "  check_err_open_lseek(result);\n"
         << "  result = write(fd, \"\", 1);\n"
         << "  check_err_write(result, fd);\n\n"
         << "  x = mmap(0, X_SIZE * sizeof(double), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);\n"
         << "  check_fail_mmap(x, fd);\n\n";

 if(analyzer->getSingleLoopRows()) {
   stream << "  int fd_b;\n"             // holds b vector
          << "  int fd_parents;\n"       // holds parents vector
          << "  int fd_vals;\n"          // holds values vector
          << "  int fd_rowPtr;\n"        // holds rowPtr vector
          << "  int fd_rowIndices;\n\n"  // holds rowIndices vector
          << "  double* b;\n"
          << "  int* parents;\n"
          << "  int* rowPtr;\n"
          << "  int* rowIndices;\n"
          << "  double* values;\n\n"
          << "  fd_b = open(FILEPATH_B, O_RDWR);\n"
          << "  fd_parents = open(FILEPATH_PARENTS, O_RDWR);\n"
          << "  fd_rowPtr = open(FILEPATH_ROWPTR, O_RDWR);\n"
          << "  fd_rowIndices = open(FILEPATH_ROWINDICES, O_RDWR);\n"
          << "  fd_vals = open(FILEPATH_VALUES, O_RDWR);\n\n"
          << "  check_err_open_lseek(fd_b);\n"
          << "  check_err_open_lseek(fd_parents);\n"
          << "  check_err_open_lseek(fd_rowPtr);\n"
          << "  check_err_open_lseek(fd_rowIndices);\n"
          << "  check_err_open_lseek(fd_vals);\n\n"
          << "  parents = mmap(0, PARENTS_SIZE * sizeof(int), PROT_READ, MAP_SHARED, fd_parents, 0);\n"
          << "  values = mmap(0, PARENTS_SIZE * sizeof(double), PROT_READ, MAP_SHARED, fd_vals, 0);\n"
          << "  rowPtr = mmap(0, X_SIZE * sizeof(int), PROT_READ, MAP_SHARED, fd_rowPtr, 0);\n"
          << "  rowIndices = mmap(0, FILESIZE_ROWINDICES, PROT_READ, MAP_SHARED, fd_rowIndices, 0);\n"
          << "  b = mmap(0, X_SIZE * sizeof(double), PROT_READ, MAP_SHARED, fd_b, 0);\n\n"
          << "  check_fail_mmap(b, fd_b);\n"
          << "  check_fail_mmap((double*)rowPtr, fd_rowPtr);\n"
          << "  check_fail_mmap((double*)rowIndices, fd_rowIndices);\n"
          << "  check_fail_mmap((double*)parents, fd_parents);\n"
          << "  check_fail_mmap(values, fd_vals);\n\n";
 }

  stream << "START();\n";   
  for(int i = 0; i < numOfParts; i++) {
    int sign;
    if((i+1) * TABLE_SIZE > signatureLevel.size())
      sign = accumulate(signatureLevel.begin() + i * TABLE_SIZE, signatureLevel.end(), 0);
    else
      sign = accumulate(signatureLevel.begin() + i * TABLE_SIZE, signatureLevel.begin() + ((i+1) * TABLE_SIZE), 0);

    if(sign == 0)
      stream << "  run" << i << "(x);\n"; 
    else
      stream << "  run" << i << "(x, b, parents, values, rowPtr, rowIndices);\n"; 
  }

  stream << "STOP(\"execution\");\n\n";   

  stream << "  int errCnt = 0;\n"
         << "  for (int i = 0; i < " << (L->getRows()-1) << " ; i++) {\n"
         << "    if(((fabs(1.0000000000 - x[i])/1.0000000000) >= 1e-2) && fabs(x[i]) >= 0) {\n"
         << "      errCnt++;\n"
         << "      printf(\"x[%d]: %.5f\\n\",i,x[i]);\n"
         << "    }\n"
         << "  }\n\n";

  stream << "  if(!errCnt)\n"
         << "    printf(\"Chainbreaker passed!\\n\");\n"
         << "  else\n"
         << "    printf(\"Chainbreaker failed! errCnt:%d\\n\",errCnt);\n\n";

  if(analyzer->getSingleLoopRows())
    stream << "  if(munmap(b, X_SIZE * sizeof(double)) == -1 || \n \
                    munmap(parents, PARENTS_SIZE * sizeof(int)) == -1 || munmap(values, PARENTS_SIZE * sizeof(double)) == -1 ||  \n \
                    munmap(rowPtr, X_SIZE * sizeof(int)) == -1 || munmap(rowIndices, FILESIZE_ROWINDICES) == -1)\n"
           << "    perror(\"Error un-mmapping the file\");\n\n"
           << "  close(fd_b);\n"
           << "  close(fd_parents);\n"
           << "  close(fd_rowPtr);\n"
           << "  close(fd_rowIndices);\n"
           << "  close(fd_vals);\n";

  stream << "  if(munmap(x, X_SIZE * sizeof(double)) == -1)\n"
         << "    perror(\"Error un-mmapping the file\");\n\n"
         << "  close(fd);\n"
         << "  return 0;\n"
         << "}";
}

void Rewrite::writeUtil() {
  std::ofstream stream(fileName + "/util.h");
  if(!stream.is_open())
    std::cout << "Cannot open output file!\n";

  stream << "#include <sys/types.h>\n"
         << "#include <sys/stat.h>\n"
         << "#include <unistd.h>\n"
         << "#include <fcntl.h>\n"
         << "#include <sys/mman.h>\n\n"
         << "void check_err_open_lseek(int result) {\n"
         << "  if(result == -1) {\n"
         << "    perror(\"Error opening file/lseek for reading\");\n"
         << "    exit(EXIT_FAILURE);\n"
         << "  }\n"
         << "}\n\n"
         << "void check_fail_mmap(double* map, int fd) {\n"
         << "  if(map == MAP_FAILED) {\n"
         << "    close(fd);\n"
         << "    perror(\"Error mmapping the file\");\n"
         << "    exit(EXIT_FAILURE);\n"
         << "  }\n"
         << "}\n\n"
         << "void check_err_write(int result, int fd) {\n"
         << "  if(result != 1) {\n"
         << "    close(fd);\n"
         << "    perror(\"Error writing last byte of the file\");\n"
         << "    exit(EXIT_FAILURE);\n"
         << "  }\n"
         << "}\n\n";

  stream.close();
}

#ifdef REWRITE_ENABLED
int Rewrite::rewriteRow(int row, vector<double> &b, string& currRow, int rewriteDepth) {
  int flops = 0;

  DAG& dag = analyzer->getDAG();
  vector<vector<double>>& values = analyzer->getValues();
  vector<double>& rowValues = values[row];

  vector<int>& parents = rewritingStrategy->isRewritten(row) ? rewritingStrategy->getInitialParentsOf(row) : dag[row].first;

  if(parents.empty()) {
    if(rewriteDepth == REWRITE_DEPTH)
      currRow.append("x[" + to_string(row) + "]");
    else {
      //currRow.append("b[" + to_string(row) + "] / " + to_string(rowValues.back()));
      currRow.append(to_string(b[row]) + " / " + to_string(rowValues.back()));
      flops++;
    }
  } else {
    currRow += "((" + to_string(b[row]) + " - (";
    flops++;
    for(int j = 0 ; j < parents.size() ; j++) {
      flops++;
      if(rewriteDepth == REWRITE_DEPTH) {
        currRow += to_string(rowValues[j]) + "*x[" + to_string(parents[j]) + "]";
      } else {
        currRow += to_string(rowValues[j]) + "*(";
        flops += rewriteRow(parents[j], b, currRow, rewriteDepth+1);
        currRow += string(")");
      }

      if(j != parents.size() - 1) {
        currRow += string(" + ");
        flops++;
      }
    }

    currRow.append(string(")) / ") + to_string(rowValues.back()) + ")");
    flops++;
  }

  //cout << "row: " << row << " added " << flops << "\n";
  return flops;
}

int Rewrite::dumpString(int row, stringstream& currRow, set<int>& rewritten, map<int,double>& multipliers) {
  int flops = (multipliers.size() << 1);

//  cout << "dumpString multipliers[-1]: " << multipliers[-1] << "\nmultipliers size: " << multipliers.size() << "\n";
  currRow << multipliers[-1];
  multipliers.erase(-1);

  for(auto it = multipliers.begin() ; it != multipliers.end(); it++)
    if(rewritten.find(it->first) == rewritten.end()) {
      currRow << " - " << it->second <<  " * x[" << it->first <<  "]";
    }

  currRow << ";\n";
//  cout << "generated string: " << currRow << "\n";

  return flops;
}

// doesnt need to return flops
// put -1 as key for constant to multipliers
// map<int,double> multipliers;  // the constant will have the key -1
// rewritingMultiplicant is 1 (this is for rewritten)
// sign: + is true, - is false
void Rewrite::rewriteRow2(int row, vector<double> &b, set<int>& rewritten, map<int,double>& multipliers, double rewritingMultiplicant, bool sign) {

  DAG& dag = analyzer->getDAG();
  vector<vector<double>>& values = analyzer->getValues();
  vector<double>& rowValues = values[row];

  vector<int>& parents = rewritingStrategy->isRewritten(row) ? rewritingStrategy->getInitialParentsOf(row) : dag[row].first;

  /*cout.fixed;
  cout.precision(10);
  if(row == 5) cout << "(0) row: " << row << " rewritingMultiplicant: " << rewritingMultiplicant <<  " sign: " << (int)sign << "\n"; */
  if(!parents.empty()) {
    for(int j = 0 ; j < parents.size() ; j++) {
      auto it = rewritten.find(parents[j]);
      if(it != rewritten.end()) {
//        if(row == 5) cout << fixed << setprecision(10) << "(1) rewritten row: " << parents[j] << " rowValues[" << j << "]: " << rowValues[j] << " rowValues.back(): " << rowValues.back() << "\n";
        rewriteRow2(parents[j], b, rewritten, multipliers, rewritingMultiplicant * (rowValues[j] / rowValues.back()), !sign);
      } else {
  //if(row == 5)      cout << fixed << setprecision(10) << "(2) sign: " << sign << " rowValues[" << j  << "]: " << rowValues[j] << " rowValues.back(): " << rowValues.back() <<  " for x" << parents[j] << ": " << rewritingMultiplicant * (rowValues[j] /rowValues.back()) << "\n"; 
//        if(multipliers.find(parents[j]) != multipliers.end())
          if(sign == true) {
    //if(row == 5)       cout << " (3) adding to multipliers[" << parents[j] << "]: " << multipliers[parents[j]] << "\n";
           multipliers[parents[j]] +=  rewritingMultiplicant * (rowValues[j] /rowValues.back());
          } else {
      //if(row == 5)     cout << " (4) removing from multipliers[" << parents[j] << "]: " << multipliers[parents[j]] << "\n";
           multipliers[parents[j]] -=  rewritingMultiplicant * (rowValues[j] /rowValues.back());
          }
  //      else
    //      multipliers[parents[j]] = rewritingMultiplicant * (rowValues[j] /rowValues.back());
      }
    }
  
//if(row == 5)    cout << fixed << setprecision(10) << " (5) sign: " << sign << " b[" <<  row << "] (" << b[row] << ") / " << rowValues.back()  << "  rowValues.back() to the constant slot\n\n";
//if(row == 5)    cout << "row: " << row << " (6) rewritingMultiplicant: " << rewritingMultiplicant << "\n";
    if(sign == true)
      multipliers[-1] += rewritingMultiplicant * (b[row] / rowValues.back());
    else
      multipliers[-1] -= rewritingMultiplicant * (b[row] / rowValues.back());
  } else { // parents empty means I rewritten to level 0.
//if(row == 5)    cout << "(7) row: " << row << " rewritingMultiplicant: " << rewritingMultiplicant << " sign: " << sign << "\n";
//if(row == 5)    cout << fixed << setprecision(10) << "setting b[" << row << "]: " << b[row] << " / rowValues.back(): " << rowValues.back() << " for x" << row << "\n\n";
    if(sign == true)
      multipliers[-1] += rewritingMultiplicant * (b[row] / rowValues.back());
    else
      multipliers[-1] -= rewritingMultiplicant * (b[row] / rowValues.back());
  }

//if(row == 5)  cout << fixed << setprecision(10) << "(8) multipliers[-1]: " << multipliers[-1] << "\n";
}

int Rewrite::rewriteRow(int row, vector<double> &b, string& currRow, set<int>& rewritten) {
  int flops = 0;

  DAG& dag = analyzer->getDAG();
  vector<vector<double>>& values = analyzer->getValues();
  vector<double>& rowValues = values[row];

  vector<int>& parents = rewritingStrategy->isRewritten(row) ? rewritingStrategy->getInitialParentsOf(row) : dag[row].first;

  if(parents.empty()) {
    auto it = rewritten.find(row);
    if(it == rewritten.end())
      currRow.append("x[" + to_string(row) + "]");
    else {
      //currRow.append("b[" + to_string(row) + "] / " + to_string(rowValues.back()));
      //currRow.append(to_string(b[row]) + " / " + to_string(rowValues.back()));
      currRow.append(to_string(b[row] / rowValues.back()));
      flops++;
    }
  } else {
    auto it = rewritten.find(row);
    if(it == rewritten.end())
      currRow.append("x[" + to_string(row) + "]");
    else {
      //currRow += "((b[" + to_string(row) + "] - (";
      currRow += "((" + to_string(b[row]) + " - (";
      flops++;
      for(int j = 0 ; j < parents.size() ; j++) {
        flops++;
        auto it = rewritten.find(parents[j]);
        if(it != rewritten.end()) {
            currRow += to_string(rowValues[j]) + "*(";
     //       cout << "calling for parent " << parents[j] << "\n";
            flops += rewriteRow(parents[j], b, currRow, rewritten);
            currRow += string(")");
        } else
            currRow += to_string(rowValues[j]) + "*x[" + to_string(parents[j]) + "]";
  
        if(j != parents.size() - 1) {
          currRow += string(" + ");
          flops++;
        }
      }
  
      currRow.append(string(")) / ") + to_string(rowValues.back()) + ")");
      flops++;
    }
  }

  return flops;
}

// rewrite the whole level
int Rewrite::rewriteLevel(int level, vector<double> &b, string& currRow, std::ostream &stream) {
  int flops = 0;
  DAG& dag = analyzer->getDAG();
  vector<vector<double>>& values = analyzer->getValues();
  vector<vector<int>>& levelTable = analyzer->getLevelTable();

  for(auto& row : levelTable[level]) {
    set<int>& rewritten = rewritingStrategy->getRewritingMapOf(row);
    //int row = level[i];
    vector<int>& parents = dag[row].first;
    vector<double>& rowValues = values[row];

    if(parents.empty()) {
      stream << "  x[" << row << "] = " << b[row] << " / " << rowValues.back() << ";\n";
      flops++;
      continue;
    }

    ToBeRewritten& toBeRewritten = rewritingStrategy->getToBeRewritten();
    int* levels = analyzer->getLevels();

    stream << "  x[" << row << "] = (" << b[row] << " - (";
    flops++;

    string currRow;
    vector<int>& initialParents = rewritingStrategy->getInitialParentsOf(row);
    for(int j = 0 ; j < initialParents.size() ; j++) {
      currRow += to_string(rowValues[j]) + "*";
      flops++;

      if(rewritingStrategy->isScopeSelective()) {
        set<int>& rewritten = rewritingStrategy->getRewritingMapOf(row);
        flops += rewriteRow(initialParents[j], b, currRow, rewritten);
      } else
        flops += rewriteRow(initialParents[j], b, currRow, 0);

      if(j != initialParents.size() - 1) {
        currRow += string(" + ");
        flops++;
      }
    }

    stream << currRow;
    stream << ")) / " << rowValues.back() << ";\n";
    flops++;
   // flops += rewriteRow(row, b, currRow, targetLevel);
  }

  return flops;
}
#endif

/*int Rewrite::writeLevel(vector<int>& level) {
  vector<int>& rowPtrL =  matrixCSR->getL()->getRowPtr();
  // TODO: these will be string to stream except for level.size()
  #pragma omp parallel for schedule(static)
  for(int i = 0 ; i < level.size() ; i++) {
      // we'll denote mmaped rowValues and flatData as
      //                     values    and parents
      stream << "  double x" << row << " = 0;\n"
             << "  for(int j = " << rowPtrL[row] << "; j < " << (rowPtrL[row + 1] - 1) << " ; j++)\n"
             << "    x" << row << " += values[j] * x_r[parents[j]];\n\n"; 
      stream << "  x_w[" << row << "] = (b[" << row << "]-x" << row << ")/values[" << (rowPtrL[row + 1]-1) << "];\n\n";
  
      flops += ((rowPtrL[row + 1] - rowPtrL[row] - 1) << 1) + 1;
  }
}*/

int Rewrite::writePart(int levelNum, int rowStartIndex, int rowEndIndex, int levelPart, vector<double> &b, std::ofstream& stream, int merged, int* sigType) {
  int flops = 0;

  vector<vector<int>>& levelTable = analyzer->getLevelTable();
  vector<int>& level = levelTable[levelNum];

  DAG& dag = analyzer->getDAG();
  vector<vector<double>>& values = analyzer->getValues();

  vector<int>& rowPtrL =  matrixCSR->getL()->getRowPtr();

  stream.precision(5);

  vector<int> unrolledRows;
  vector<int> loopedRows(level.begin() + rowStartIndex, level.begin() + rowEndIndex);
  #ifdef REWRITE_ENABLED
    // TODO: this is not efficient, construct rewritten rows per level from the beginning
    vector<int> rewrittenRows;
  #endif

  analyzer->separateRows(levelNum, rowStartIndex, rowEndIndex, loopedRows, unrolledRows);
 
  if(!loopedRows.empty()) {
  //if(loopedRows.size() > 1) {
    /*cout << "level: " << levelNum << " startIndex: " << rowStartIndex << " endIndex: " << rowEndIndex << "\n";
    for(auto& row : loopedRows)
      cout << row << ", ";
    cout << "\n";*/

    if(rowStartIndex == 0) {
      signatureLevel.push_back(1);
      signature.push_back(vector<int>());
    } else
      signatureLevel[levelNum] = 1;

    signature[levelNum].push_back(1);  // if we OR signature[levelNum] == signatureLevel[levelNum]

    // remove rewritten rows from loopedRows, fill in rewrittenRows
    vector<int>::iterator it = loopedRows.begin();
    while(it != loopedRows.end()) {
      #ifdef REWRITE_ENABLED
        if(rewritingStrategy->isRewritten(*it)) {
          rewrittenRows.push_back(*it);
          it = loopedRows.erase(it);
        } else {
      #endif
        flops += ((dag[*it].first.size()) << 1) + 1;
        it++;
      #ifdef REWRITE_ENABLED
       }
      #endif
    }
  } else {
    if(rowStartIndex == 0) {
      signatureLevel.push_back(0);
      signature.push_back(vector<int>());
    }
      
    signature[levelNum].push_back(0);
  }

  #ifdef REWRITE_ENABLED
    if(!unrolledRows.empty()) {
      // remove rewritten rows from unrolledRows, fill in rewrittenRows
      vector<int>::iterator it = unrolledRows.begin();
      while(it != unrolledRows.end()) {
        if(rewritingStrategy->isRewritten(*it)) {
          rewrittenRows.push_back(*it);
          it = unrolledRows.erase(it);
        } else
          it++;
      }
    }
  #endif

  rowIndices.insert(rowIndices.end(), loopedRows.begin(), loopedRows.end());
  vector<int>& currLevelStartIndex = startIndex[levelNum];
  if(currLevelStartIndex.empty()) // this will work even if loopedRows is empty since start & end will be the same for the loop
    currLevelStartIndex.push_back(startIndex[levelNum-1].back() + loopedRows.size());
  else
    currLevelStartIndex.push_back(currLevelStartIndex.back() + loopedRows.size());

  /*cout << "rowIndices:\n";
  for(auto& index : rowIndices)
    cout << index << ", ";
  cout << "\n";*/


  if(!loopedRows.empty()) {
    //if(merged && *sigType == 1)
      *sigType = 2;

    // TODO: we actually can remove rowIndices from loopedRows with size == 1. That'll complicate things though.
    if(!merged)
      stream << "void calculate" << levelPart << "(double* x, double* b, int* parents, double* values, int* rowPtr, int* rowIndices) {\n";

    if(loopedRows.size() > 1) {
      if(currLevelStartIndex.size() == 1)
        stream << "  for(int i = " << startIndex[levelNum-1].back() << " ; i < " << currLevelStartIndex.back() << " ; i++) {\n";
      else
        stream << "  for(int i = " << currLevelStartIndex[currLevelStartIndex.size()-2] << " ; i < " << currLevelStartIndex.back() << " ; i++) {\n";

      stream << "    int row = rowIndices[i];\n";
    } else
        stream << "    int row = " << loopedRows[0] << ";\n";
          
    stream << "    double xi = 0;\n"
           << "    for (int j = rowPtr[row]; j < rowPtr[row+1]-1; j++)\n"
           << "      xi += values[j] * x[parents[j]];\n\n"
           << "    x[row] = (b[row]-xi)/values[rowPtr[row+1]-1];\n";

    if(loopedRows.size() > 1)
      stream << "  }\n\n";
  }

  if(!unrolledRows.empty()) {
 //   if(loopedRows.empty())
    if(loopedRows.empty() && !merged) {
      stream << "  void calculate" << levelPart << "(double* x) {\n";
      stream << "//void calculate" << levelPart << "(double* x, double* b, int* parents, double* values, int* rowPtr, int* rowIndices) {\n";
    }

    for(auto& row : unrolledRows) {

      vector<int>& parents = dag[row].first;
      vector<double>& rowValues = values[row];

      if(parents.empty()) {
  //      stream << setprecision(10) << "  x[" << row << "] = " << b[row] / rowValues.back() << ";\n";
        stream << std::scientific << "  x[" << row << "] = " << b[row] / rowValues.back() << ";\n";
        flops++;
      } else if(parents.size() == 1) {
        stream << std::scientific << "  x[" << row << "] = (" << b[row] << "-(" << rowValues[0] << ") * x[" << parents[0] << "])/" << rowValues.back() << ";\n";
        flops+= 3;
      } else if(parents.size() == 2) {
        stream << std::scientific << "  x[" << row << "] = (" << b[row] << "-((" \
               << rowValues[0] << ") * x[" << parents[0] << "]+(" \
               << rowValues[1] << ") * x[" << parents[1] << "]))/" << rowValues.back() << ";\n";
        flops+= 5;
      } else if(parents.size() == 3) {
        stream << std::scientific << "  x[" << row << "] = (" << b[row] << "-((" \
               << rowValues[0] << ") * x[" << parents[0] << "]+(" \
               << rowValues[1] << ") * x[" << parents[1] << "]+(" \
               << rowValues[2] << ") * x[" << parents[2] << "]))/" << rowValues.back() << ";\n";
        flops+= 7;
      } else if(parents.size() == 4) {
        stream << std::scientific << "  x[" << row << "] = (" << b[row] << "-((" \
               << rowValues[0] << ") * x[" << parents[0] << "]+(" \
               << rowValues[1] << ") * x[" << parents[1] << "]+(" \
               << rowValues[2] << ") * x[" << parents[2] << "]+(" \
               << rowValues[3] << ") * x[" << parents[3] << "]))/" << rowValues.back() << ";\n";
        flops+= 9;
      }
    }
  }

  #ifdef REWRITE_ENABLED
  if(!rewrittenRows.empty()) {
    if(loopedRows.empty() && unrolledRows.empty())
      if(!merged) {
        stream << "  void calculate" << levelPart << "(double* x) {\n";
        stream << "//void calculate" << levelPart << "(double* x, double* b, int* parents, double* values, int* rowPtr, int* rowIndices) {\n";
      }
    }

    ToBeRewritten& toBeRewritten = rewritingStrategy->getToBeRewritten();
    int* levels = analyzer->getLevels();

//    cout << "# of rewritten rows: " << rewrittenRows.size() << "\n";
    for(auto& row : rewrittenRows) {
      if(!rewritingStrategy->isRewritten(row))
        cout << "row: " << row << " doesnt appear in rewrittenStrategy's list\n"; 

      vector<int>& parents = dag[row].first;
      stream << "  x[" << row << "] = ";

      vector<int>& initialParents = rewritingStrategy->getInitialParentsOf(row);
      set<int>& rewritten = rewritingStrategy->getRewritingMapOf(row);
      map<int,double> multipliers;
      for(auto& row : rewritten)
        multipliers[row] = 0;
      multipliers[-1] = 0;
      double rewritingMultiplicant = 1;

      rewriteRow2(row, b, rewritten, multipliers, rewritingMultiplicant, true);
    
  auto t9 = std::chrono::steady_clock::now();
      // get updated matrix values for rewritten rows
      rewrittenB[row] = multipliers[-1];

      rewrittenRowParents[row] = vector<int>();
      rewrittenRowValues[row] = vector<double>();
      vector<double>& rewrittenValues = rewrittenRowValues[row];
      vector<int>& rewrittenParents = rewrittenRowParents[row];
  auto t10 = std::chrono::steady_clock::now();
  cout << " chrono updated values: " << std::chrono::duration_cast<std::chrono::duration<double>>(t10 - t9).count()*1000 << "\n";

      stringstream currRow;

      flops += dumpString(row, currRow, rewritten, multipliers);
      stream << currRow.str();
    
/*      ToBeRewritten& toBeRewritten = rewritingStrategy->getToBeRewritten();
      cout << "multipliers size: " << multipliers.size() << "\n";
      cout << "rewriting distance: " << toBeRewritten[row].second - toBeRewritten[row].first << "\n";*/
   t9 = std::chrono::steady_clock::now();
      for(auto it = multipliers.begin() ; it != multipliers.end(); it++) {
        if(rewritten.find(it->first) == rewritten.end()) {
          rewrittenValues.push_back(it->second);
          rewrittenParents.push_back(it->first);
        }
      }

  t10 = std::chrono::steady_clock::now();
  cout << " chrono updated values: " << std::chrono::duration_cast<std::chrono::duration<double>>(t10 - t9).count()*1000 << "\n";
    }
  #endif

  if(!merged)
    stream << "}\n\n";

  return flops;
}

int Rewrite::writeMakefile(int execBlockCnt) {
  std::ofstream stream(fileName + "/Makefile");
  if(!stream.is_open()) {
    std::cout << "Cannot open output file!\n";
    return -1;
  }

  stream << "CXX      := clang\n"
  //<< "CXXFLAGS := -O3 -ffast-math -funroll-loops -march=native -I/kuacc/apps/llvm-omp/include\n"
  << "CXXFLAGS := -O3 -flto=thin -ffast-math -march=native \n"
  << "LDFLAGS  := -fopenmp\n"
  << "SRC      := $(wildcard calculate*.c)\n"
  << "OBJECTS  := $(SRC:%.c=%.o)\n"
  << "SRC_RUN      := $(wildcard run*.c)\n"
  << "OBJECTS_RUN  := $(SRC_RUN:%.c=%.o)\n\n"
  << "%.o: %.c\n"
  << "\t$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ -c $<\n\n"
  << "all: $(OBJECTS) $(OBJECTS_RUN)\n"
  << "\t$(CXX) $(CXXFLAGS) -c main.c \n"
  << "\t$(CXX) $(CXXFLAGS) $(OBJECTS) $(OBJECTS_RUN) main.o $(LDFLAGS) -o chainbreaker\n\n";

  stream << "clean:\n"
         << "\t$(RM) $(OBJECTS) $(OBJECTS_RUN) main.o chainbreaker";

  stream << "\n\n";

  stream.close();
  return 0;
}

int Rewrite::writeHeader(int headerCounter, int partCounterStart, int partCounterEnd) {
  std::ofstream stream(fileName + "/calculators" + to_string(headerCounter) + ".h");
  if(!stream.is_open()) {
    std::cout << "Cannot open output file!\n";
    return -1;
  }

  int startingLevel = headerCounter * TABLE_SIZE;
  int counter = 0;
  stream << "#pragma once\n";
  for(int i = partCounterStart; i < partCounterEnd; i++) {
    if(signature[startingLevel][counter] == 0)
      stream << "void calculate" << to_string(i) << "(double* x);\n";
    else
      stream << "void calculate" << to_string(i) << "(double* x, double* b, int* parents, double* values, int* rowPtr, int* rowIndices);\n";

    counter++;
    if(counter == signature[startingLevel].size()) {
      startingLevel++;
      counter = 0;
    }
  }

  stream.close();
  return 0;
}

int Rewrite::balanceLevel(int toBeBalanced, int& workloadPerThread) {
    int numThreads = NUM_THREADS;
    workloadPerThread = toBeBalanced+NUM_THREADS-1;
    workloadPerThread /= NUM_THREADS;

 //     cout << "toBeBalanced size: " << toBeBalanced << "\n";
    if(workloadPerThread < LOWER_BOUND) {
      numThreads = toBeBalanced+LOWER_BOUND-1;
      numThreads /= LOWER_BOUND;
      workloadPerThread = toBeBalanced+numThreads-1;
      workloadPerThread /= numThreads;
 //       cout << "updated numThreads: " << numThreads << "\n";
 //       cout << "updated workloadPerThread: " << workloadPerThread << "\n";
    } else
      numThreads = (const int)ceil((double)(toBeBalanced/workloadPerThread));

    return numThreads;
}

#ifdef BALANCE_FLOPS
int Rewrite::balanceRows(vector<int>& level) {
//    int flopsPerThread = (const int)ceil((double)flopsPerLevel[j]/NUM_THREADS);
  
  // rowDistBelowAvg and rowdistAboveAvg keep row distribution to threads. rowDistBelowAvg
  // cooporate with flopSum which keeps track of sum of flops in rows assigned to a thread.
  // [GOAL]: distribute rows to threads where sum of flops assigned to each thread
  // is equally likely. A row can be assigned to 1 thread at most, multiple rows can
  // be assigned to a thread. We couldn't use a map instead of these 3 vectors since
  // sum of flops assigned to a thread does not need to satisfy uniqueness propoerty.
  // As needed new slots can be opened until we hit NUM_THREADS. Upon hitting, if there's
  // still a row unassigned, the one sitting at the end of rowDistBelowAvg is used : this
  // could be improved but overhead will become too costly if we keep rowDistBelowAvg and flopSum
  // sorted and try selecting the best candidate row to add on top.
  //

  // TODO: optimize for levels with all rows with same length: use BALANCE_ROWS for these
  // TODO: refactor this section

  cout << "\nnum of rows in level " << j << " : " << level.size() << "\n";
  cout << "flopsPerLevel: " << flopsPerLevel[j] << "\n";

  int avgFlopsPerLevel = flopsPerLevel[j]/level.size();
  cout << "avgFlopsPerLevel: " << avgFlopsPerLevel << "\n";

  vector<vector<int>> rowDistBelowAvg;
  vector<vector<int>> rowDistAboveAvg;
  vector<int> flopSum(NUM_THREADS,0);

  DAG& dag = analyzer->getDAG();

  for(int k = 0 ; k < level.size(); k++) {
    //int rowFlops = ((rowPtrL[level[k] + 1] - rowPtrL[level[k]]) << 1) - 1;
    vector<int>& parents = dag[level[k]].first;
    int rowFlops = ((parents.size()+1) << 1) - 1;
    cout << "row: " << level[k] << " with flops: " << rowFlops << "\n";
    if((rowDistAboveAvg.size() + rowDistBelowAvg.size()) > NUM_THREADS) {
      /*if(rowDistBelowAvg.empty()) {
        if(rowFlops <= avgFlopsPerLevel) {
        }
      } else*/
      // find the min
      int minIndex = 0;
      for(int m = 0 ; m < flopSum.size(); m++)
        if(flopSum[m] < flopSum[minIndex])
          minIndex = m;

      if(rowFlops <= avgFlopsPerLevel) {
        vector<int>& toBeMerged = rowDistBelowAvg[minIndex];
         toBeMerged.push_back(level[k]);
        if(rowFlops + flopSum[minIndex] > avgFlopsPerLevel) {
          cout << "rowDistBelowAvg size: " << rowDistBelowAvg.size() << "\n";
          rowDistAboveAvg.push_back(toBeMerged);
          if(rowDistBelowAvg.size() == 1) {
            rowDistBelowAvg.erase(rowDistBelowAvg.end() - 1);
            flopSum.erase(flopSum.end() - 1);
          }  else {
            rowDistBelowAvg.erase(rowDistBelowAvg.begin() + minIndex);
            flopSum.erase(flopSum.begin() + minIndex);
          }
        }
      } else {
        // find the second min
        int sndMinIndex = 0;
        for(int m = 0 ; m < flopSum.size(); m++)
          if(flopSum[m] < flopSum[sndMinIndex] && m != minIndex)
            sndMinIndex = m;

        vector<int>& toBeMerged = rowDistBelowAvg[minIndex];
        vector<int>& toBeMerged2 = rowDistBelowAvg[sndMinIndex];
        toBeMerged.insert(toBeMerged.end(), toBeMerged2.begin(), toBeMerged2.end());

        if(flopSum[minIndex] + flopSum[sndMinIndex] < avgFlopsPerLevel)
          flopSum[minIndex] += flopSum[sndMinIndex];
        else {   
          rowDistAboveAvg.push_back(toBeMerged);

          flopSum.erase(flopSum.begin() + minIndex);
          rowDistBelowAvg.erase(rowDistBelowAvg.begin() + minIndex);
          if(minIndex < sndMinIndex)
            sndMinIndex--;
        }
          
        flopSum.erase(flopSum.begin() + sndMinIndex);
        rowDistBelowAvg.erase(rowDistBelowAvg.begin() + sndMinIndex);
      }
    } else {
      if(rowFlops <= avgFlopsPerLevel) {
        if(rowDistBelowAvg.empty()) {
          rowDistBelowAvg.push_back(vector<int>(1,level[k])); 
          flopSum[0] = rowFlops;
        } else {
          int i = 0;
          while((i < rowDistBelowAvg.size()) && (flopSum[i] + rowFlops > avgFlopsPerLevel))
            i++;
       
          // we couldn't find a suitable candidate. Use the last slot in rowDistBelowAvg
          // and move it to rowDistAboveAvg. not the best approach but keeps from increasing
          // the overhead
          if(i != rowDistBelowAvg.size()) { 
            flopSum[i] += rowFlops;
            vector<int>& newSlot = rowDistBelowAvg[i];
            newSlot.push_back(level[k]);
          } else {
            if(rowDistBelowAvg.empty())
              rowDistAboveAvg.push_back(vector<int>()); 
            else
              rowDistAboveAvg.push_back(rowDistBelowAvg.back()); 

            vector<int>& newSlot = rowDistAboveAvg.back();
            newSlot.push_back(level[k]);
     //       cout << "rowDistAboveAvg has " << newSlot[0] << "\n";

            if(!rowDistBelowAvg.empty()) {
              rowDistBelowAvg.erase(rowDistBelowAvg.end()-1);
              flopSum.erase(flopSum.end()-1);
            }
          }
        }
      } else {
         rowDistAboveAvg.push_back(vector<int>()); 
         vector<int>& newSlot = rowDistAboveAvg.back();
         newSlot.push_back(level[k]);
      }
    }
  }

  cout << "numThreads after BALANCE_FLOPS: " << rowDistAboveAvg.size() << ", " << rowDistBelowAvg.size() << "\n";

  for(int m = 0 ; m < rowDistBelowAvg.size() ; m++) {
    cout << m << ":\n";
    for(auto& slot : rowDistBelowAvg[m])
      cout << slot << ", ";
    cout << "\n";
  }

  for(int m = 0 ; m < rowDistAboveAvg.size() ; m++) {
    cout << rowDistBelowAvg.size() + m << ":\n";
    for(auto& slot : rowDistAboveAvg[m])
      cout << slot << ", ";
    cout << "\n";
  }

  flopSum.shrink_to_fit();
  numThreads = rowDistAboveAvg.size() + rowDistBelowAvg.size();
  vector<vector<int>> rowDist(rowDistBelowAvg);
  rowDist.insert(rowDist.end(), rowDistAboveAvg.begin(),rowDistAboveAvg.end());
}
#endif

#ifdef BALANCE_ROWS_CV
// calculate coeff. of variation(CV) = stdev/mean
void Rewrite::collectFLOPS(vector<int>& level, vector<int>& flopsLevel) {
  DAG& dag = analyzer->getDAG();

 for(auto& row : level)
  flopsLevel.push_back(((dag[row].first.size()+1) << 1) - 1);

 /*for(auto& row : level)
   cout << row << ", ";
 cout << "\n";*/
}

double Rewrite::calculateCV(int numOfSum, int start, vector<int>& level) {
  /*for(auto& i : level)
    cout << i << ", ";
  cout << "\n";*/

  double mean = accumulate(level.begin(),level.end(),0.0)/(double)numOfSum;

  double sum = 0;
  for(int i = start; i < numOfSum; i++ )
    sum += (level[i] - mean) * (level[i] - mean);

  sum /= (numOfSum-1);

 double CV = sqrt(sum)/mean;
// cout << "start: " << start << " numOfSum: " << numOfSum << "\n";
// cout << "CV: " << CV << "\n";

 return CV;
}

void Rewrite::advanceTillBinLimit(vector<int>& level, int start, int& end, int binLimit) {
  int sum = accumulate(level.begin()+start, level.begin()+end, 0);

  while(end < level.size() && sum < binLimit) {
    sum += level[end];
    end = end+1;
  }
}

void Rewrite::binLastElement(vector<int>& level, int& start, int& end, vector<int>& rowDist, vector<int>& rowDistSum) {
  // if we reached the end consider whether to add the last element to the bin or not
  double CV_start_end_1 = 0;
  int currSum = 0;

  if(start == level.size()-1) {
    rowDist.push_back(1);
    rowDistSum.push_back(level.back());
    start++;

   /* cout << "rowDist:\n";
    for(auto& i : rowDist)
      cout << i << ", ";
    cout << "\n";

    cout << "rowDistSum:\n";
    for(auto& i : rowDistSum)
      cout << i << ", ";
    cout << "\n";*/

    return;
  }

  if(end == start + 1)
    currSum = accumulate(level.begin()+start, level.begin()+end, 0);
  else
    currSum = accumulate(level.begin()+start, level.begin()+(end-1), 0);

//  cout << "checking binLastElement:\n";
//  cout << "start: " << start << " end: " << end << " currSum: " << currSum << "\n";
  vector<int> currRowDistSum(rowDistSum);
  currRowDistSum.push_back(currSum);

/*  cout << "currRowDistSum:\n";
  for(auto& i : currRowDistSum)
    cout << i << ", ";
  cout << "\n";

  cout << "end: " << end << " level.size(): " << level.size() << "\n";*/

 // if(end == level.size())
 //   currRowDistSum.push_back(level[end-1]);

  if(currRowDistSum.size() > 1)
    CV_start_end_1 = calculateCV(currRowDistSum.size(), 0, currRowDistSum);

//  cout << "CV_start_end_1: " << CV_start_end_1 << "\n";
 
  if(end == start + 1)
    end++;
    
  currSum = accumulate(level.begin()+start, level.begin()+end, 0);

// currSum += *(level.end());
// cout << "currSum: " << currSum << " *level.end(): " << *(level.end()) << "\n";
  
  vector<int> currRowDistSumEnd(rowDistSum);
  currRowDistSumEnd.push_back(currSum);

  double CV_start_end = 0;
  if(currRowDistSum.size() > 1)
    CV_start_end = calculateCV(currRowDistSumEnd.size(), 0, currRowDistSumEnd);

/*  cout << "currSum: " << currSum << "\n";  

  cout << "currRowDistSumEnd:\n";
  for(auto& i : currRowDistSumEnd)
    cout << i << ", ";
  cout << "\n";

  cout << "CV_start_end: " << CV_start_end << "\n";*/
  
  if(CV_start_end_1 < CV_start_end) {
    rowDist.push_back(end-1-start);    
    rowDistSum.push_back(accumulate(level.begin()+start, level.begin()+(end-1), 0));    
    start = end-1;
  } else {
    rowDist.push_back(end-start);    
    rowDistSum.push_back(accumulate(level.begin()+start, level.begin()+end, 0));    
    start = end;
    end = end+1;
  }

/*  cout << "rowDist:\n";
  for(auto& i : rowDist)
    cout << i << ", ";
  cout << "\n";

  cout << "rowDistSum:\n";
  for(auto& i : rowDistSum)
    cout << i << ", ";
  cout << "\n";*/
}

double Rewrite::balanceRowsCV(vector<int>& level, vector<int>& rowDist, vector<int>& rowDistSum, double totalCV) {
//  double totalCV = calculateCV(level.size(), 0, level);
  //int binLimit = (*max_element(begin(level), end(level))) << 1;
  int binLimit = (*max_element(begin(level), end(level)));
  int start = 0; int end = 1; 

/*  cout << "balanceRowsCV: " << totalCV << ", ";
  cout << "CV before: " << totalCV << "\n";
  cout << "binLimit: " << binLimit << "\n";*/

  if(totalCV > 0.05) {
    // construct the first bin. This is not ideal but we need an initial bin and dont know the rest of the values yet
    // hence we fill the first bin no matter what
    advanceTillBinLimit(level, start, end, binLimit);
//    cout << "new end: " << end << "\n";

    if(end == start + 1){
      rowDist.push_back(1);    
      rowDistSum.push_back(level[start]);
      start = end;
      end++;
    } else {
      rowDist.push_back(end-1-start);    
      rowDistSum.push_back(accumulate(level.begin()+start, level.begin()+(end-1), 0));    
      start = end-1;
    }

/*    cout << "new start: " << start << " new end: " << end << "\n";

    cout << "rowDist:\n";
    for(auto& i : rowDist)
      cout << i << ", ";
    cout << "\n";

    cout << "rowDistSum:\n";
    for(auto& i : rowDistSum)
      cout << i << ", ";
    cout << "\n";*/

    //while(end <= level.size()) {
    while(start < level.size()) {
      advanceTillBinLimit(level, start, end, binLimit);
//      cout << "new end: " << end << "\n";
      binLastElement(level, start, end, rowDist, rowDistSum);
    }
  } else {
    // if variation is already < 0.05 do not bother running the algorithm
    rowDist.push_back(level.size());    
    rowDistSum.push_back(accumulate(level.begin(), level.end(), 0));    
  }

  double totalCVAfter = totalCV;
  if(rowDistSum.size() > 1)
    totalCVAfter = calculateCV(rowDistSum.size(), 0, rowDistSum);

  // if we increased the CV after running the algorithm, revert back
  // this can happen due to non-ideal condition for constructing the
  // first bin
  if((totalCVAfter > totalCV) && (level.size() <= NUM_THREADS)) {
    cout << totalCV << "\n";
    return totalCV;
  }

  cout << totalCVAfter << "\n";
  return totalCVAfter;
}

void Rewrite::reduceNumOfRows(vector<int>& level, vector<int>& rowDistReduced, vector<int>& rowDistSumReduced) {
  int partitionSize = (int)(ceil((double)level.size() / NUM_THREADS));
  int i = 0;

/*  cout << "reduceNumOfRows: levels \n";
  for(auto& i : level)
    cout << i << ", ";
  cout << "\n";

  cout << "partitionSize: " << partitionSize << "\n";*/
  while((i + partitionSize) < level.size()) {
    int sum = accumulate(level.begin()+i, level.begin()+i+partitionSize, 0);
    
    rowDistSumReduced.push_back(sum);
    rowDistReduced.push_back(partitionSize);
    i += partitionSize;
  }

  if(level.size()-i >= 1) {
    int sum = accumulate(level.begin()+i, level.end(), 0);
    rowDistSumReduced.push_back(sum);
    rowDistReduced.push_back(level.size()-i);
  }

/*  cout << "rowDistReduced:\n";
  for(auto& i : rowDistReduced)
    cout << i << ", ";
  cout << "\n";

  cout << "rowDistSumReduced:\n";
  for(auto& i : rowDistSumReduced)
    cout << i << ", ";
  cout << "\n";*/
}
#endif

int Rewrite::rewriteExecutor(vector<double> &b, vector<double> &x) {
  Part *L = matrixCSR->getL();
  writeUtil();

  vector<vector<int>>& levelTable = analyzer->getLevelTable();
  auto& flopsPerLevel = analyzer->getFlopsPerLevel();

  for(int i = 0 ; i < flopsPerLevel.size(); i++)
    startIndex.push_back(vector<int>());
  startIndex[0].push_back(0);

  int rowStartIndex = 0;
  int partCounter = 0;
  int headerCounter = 0;
  int partCounterStart = 0;

  int maxNumOfThreads = 1;

  // Tracker keeps start of each part & num of threads
  // part: workload of a thread
  vector< vector<int> > tracker;
  for(int i = 0 ; i < 2; i++)
    tracker.push_back(vector<int>(TABLE_SIZE,0));

//  #pragma omp parallel for private(rowStartIndex) schedule(static)

  // assign only 1 thread to level 0. Thread creating overhead is too much for such
  // easy & fast computation
  int flops = 0;
 // cout << "level: " << 0 << "\n";
    #ifdef BALANCE_ROWS
      std::ofstream streamZero(fileName + "/calculate" + to_string(partCounter) + ".c");
      if(!streamZero.is_open())
        std::cout << "Cannot open output file!\n";

      int flopsBefore = flops;
      int sig = 1;
      flops += writePart(0, 0, levelTable[0].size(), partCounter++, b, streamZero, 0, &sig);
      streamZero.close();
    #endif

  vector<int>& partStarts =   tracker[0];
  vector<int>& threadCounts = tracker[1];
  partStarts[0] = partCounter-1;
  threadCounts[0] = 1;

  int size = (0 == flopsPerLevel.size()-1) ? (flopsPerLevel.size() % TABLE_SIZE) : TABLE_SIZE; 

  if((0 == (TABLE_SIZE-1)) || (0 == flopsPerLevel.size()-1)) {
 //   cout << "j: " << 0 << " part: " << 0 << "\n"; 
    std::ofstream stream(fileName + "/run0.c");
    writeHeader(0, 0, 1);
    writeFunc(stream, tracker, size, 0, 1, headerCounter);
  }

  int prevSingleThreadCnt = 0, merged = 0, sigType = 1; // type1 : short one, type2: long one
  int prevPartCounter;
  for(int j = 1 ; j < flopsPerLevel.size(); j++) {
    int numThreads = NUM_THREADS;
    int workloadPerThread;
    vector<int>& level = levelTable[j];
    vector<int> flopsLevel;
    int flops = 0;

    #ifdef BALANCE_ROWS
      numThreads = balanceLevel(level.size(), workloadPerThread);
    #endif

    // new single threaded level found, if prev. level is also single threaded, merge
    if(numThreads == 1)
      prevSingleThreadCnt++;

    // merging finished before 10 files
    if((numThreads > 1) && (merged == 1) && (prevSingleThreadCnt < 10)) {
      closeParanthesis(fileName, prevPartCounter, sigType);

      prevSingleThreadCnt = merged = 0;
      sigType = 1;
    }

    if(prevSingleThreadCnt == 1)
      prevPartCounter = partCounter;

    if(prevSingleThreadCnt == 2) {
      cout << "prevPartCounter: " << prevPartCounter << " partCounter: " << partCounter << "\n";
      removeParanthesis(fileName, prevPartCounter);
    }

    if(numThreads > maxNumOfThreads)
      maxNumOfThreads = numThreads;

    int workloadCounter = 0;
//    cout << "level: " << j << "\n";
//    cout << "num threads: " << numThreads << "\n";
    #ifdef BALANCE_ROWS
      std::ofstream streamLevel;

      if(prevSingleThreadCnt >= 2 &&  prevSingleThreadCnt <= 10) { // write into prev file
        cout << "prevPartCounter: " << prevPartCounter << " partCounter: " << partCounter << "\n";
        streamLevel.open(fileName + "/calculate" + to_string(prevPartCounter) + ".c", std::ofstream::app);
        merged = 1;
      } else
        streamLevel.open(fileName + "/calculate" + to_string(partCounter) + ".c");

      if(!streamLevel.is_open())
        std::cout << "Cannot open output file!\n";

      for(int i = 0; i < numThreads-1; i++) {
        int flopsBefore = flops;
        flops += writePart(j, i*workloadPerThread, i*workloadPerThread+workloadPerThread, partCounter++, b, streamLevel, merged, &sigType);
    }
    
      int flopsBefore = flops;
      flops += writePart(j, (numThreads-1)*workloadPerThread, level.size(), partCounter++, b, streamLevel, merged, &sigType);

      streamLevel.close();

      if(prevSingleThreadCnt == 10) { // merge is done. Wait for a new opportunity
        closeParanthesis(fileName, prevPartCounter, sigType);

        prevSingleThreadCnt = merged = 0;
        sigType = 1;
      }
    #endif

//    cout << " part counter became: " << partCounter << "\n";;
    vector<int>& partStarts =   tracker[0];
    vector<int>& threadCounts = tracker[1];
    partStarts[(j % TABLE_SIZE)] = partCounter-numThreads;
    threadCounts[(j % TABLE_SIZE)] = numThreads;

    int size = (j == flopsPerLevel.size()-1) ? (flopsPerLevel.size() % TABLE_SIZE) : TABLE_SIZE; 

    if(((j % TABLE_SIZE) == (TABLE_SIZE-1)) || (j == flopsPerLevel.size()-1)) {
      int part = floor((((double)j)/TABLE_SIZE));
   //   cout << "j: " << j << " part: " << part << "\n"; 
      std::ofstream stream(fileName + "/run" + to_string(part) + ".c");
      writeHeader(headerCounter, partCounterStart, partCounter);
      writeFunc(stream, tracker, size, part, maxNumOfThreads, headerCounter);
      headerCounter++;
      partCounterStart = partCounter;
    }
  }

  return 0;
}

// TODO: update code according ALB updates (vector -> int*, validation check code, etc)
int Rewrite::rewrite() {
  int rows = matrixCSR->getNumOfRows();
  int cols = matrixCSR->getNumOfCols();
  int nnzs = matrixCSR->getNumOfVals();

  int *rowPtr = matrixCSR->getRowPtr();
  int *colIdx = matrixCSR->getColIdx();
  double *vals = matrixCSR->getVals();
  struct timeval t1{}, t2{};
  
  if(rows != cols) {
    printf("This is not a square matrix.n");
    return -1;
  }

  vector<double> xRef(cols, 1);
  vector<double> b(cols, 0);
  
  DAG& dag = analyzer->getDAG();
  vector<vector<double>>& values = analyzer->getValues();

  Part* LCSR = matrixCSR->getL();

  /////////////////////////////// PRINT SOME BASIC STATS ///////////////////////////////
  cout << "NUM_THREADS:" << NUM_THREADS << "\n";
  cout << "LOWER_BOUND: " << LOWER_BOUND << "\n";
  #ifdef BALANCE_ROWS
    cout << "BALANCE METHOD: BALANCE_ROWS\n";
  #else
    cout << "BALANCE METHOD: BALANCE_ROWS_CV\n";
  #endif
  cout << "TABLE_SIZE: " << TABLE_SIZE << "\n";
  /////////////////////////////// PRINT SOME BASIC STATS ///////////////////////////////

  auto t5 = std::chrono::steady_clock::now();
  // fill in the B values, these are needed for writePart
  // B values for rewritten rows will be filled in by writePart
  for(int i = 0; i < rows; i++) {
    #ifdef REWRITE_ENABLED
      vector<int>& parents = rewritingStrategy->isRewritten(i) ? rewritingStrategy->getInitialParentsOf(i) : dag[i].first;
    #else
      vector<int>& parents = dag[i].first;
    #endif

    vector<double>& rowValues = values[i];

    b[i] += rowValues.back() * xRef[i];
    for(int j = 0 ; j < parents.size() ; j++)
      b[i] += rowValues[j] * xRef[parents[j]];
  }
  auto t6 = std::chrono::steady_clock::now();
  cout << "chrono calculate B: " << std::chrono::duration_cast<std::chrono::duration<double>>(t6 - t5).count()*1000 << "\n";

  vector<double> x(cols, 0);
  rewriteExecutor(b, x);
  
  int parentsSize = dumpDataToMem(b);
  auto& flopsPerLevel = analyzer->getFlopsPerLevel();
  std::ofstream stream(fileName + "/main.c");
  if(!stream.is_open())
    std::cout << "Cannot open output file!\n";
  writeMain(stream, (int)(ceil((((double)flopsPerLevel.size())/TABLE_SIZE))), parentsSize);
  stream.close();
    
  writeMakefile(analyzer->getNumOfLevels());

  return 0;
}

int Rewrite::dumpDataToMem(vector<double>& b) {
  auto t7 = std::chrono::steady_clock::now();

  Part* LCSR = matrixCSR->getL();
  string extension(".bin");
  string extension_TR("_TR.bin");

  string pathVals ("/tmp/" + fileName + "_vals");
  string pathParents ("/tmp/" + fileName + "_parents");
  string pathRowPtr ("/tmp/" + fileName + "_rowPtr");

  int returnVal = 0;

  // dumpCSR is enabled so that we use external lib
  if(analyzer->getSingleLoopRows()) {
    string pathB ("/tmp/" + fileName + "_b");
    string pathRowIndices ("/tmp/" + fileName + "_rowIndices");

    #ifdef REWRITE_ENABLED
      vector<int> transformedRowPtr, transformedColIdx;
      vector<double> transformedValues;
      buildTransformedMatrix(transformedRowPtr, transformedColIdx, transformedValues);
      returnVal = transformedValues.size();

      // TODO: cant we update them in writePart? do we need to accumulate?
      for(auto& row : rewrittenB)
        b[row.first] = row.second;

      dumpToMem(b, pathB.append(extension_TR).c_str(), LCSR->getRows()-1);
      dumpToMem(transformedValues, pathVals.append(extension_TR).c_str(), transformedValues.size());
      dumpToMem(transformedRowPtr, pathRowPtr.append(extension_TR).c_str(), transformedRowPtr.size());
      dumpToMem(transformedColIdx, pathParents.append(extension_TR).c_str(), transformedColIdx.size());
      dumpToMem(rowIndices, pathRowIndices.append(extension_TR).c_str(), rowIndices.size());

      // do we want the original matrix in addition for SYCL code?
      if(dumpCSR) {
        dumpToMem(b, pathB.append(extension).c_str(), LCSR->getRows()-1);
        dumpToMem(LCSR->getVals(), pathVals.append(extension).c_str(), LCSR->getNNZs());
        dumpToMem(LCSR->getRowPtr(), pathRowPtr.append(extension).c_str(), LCSR->getRowPtr().size());
        dumpToMem(LCSR->getColIdx(), pathParents.append(extension).c_str(), LCSR->getColIdx().size());
        dumpToMem(rowIndices, pathRowIndices.append(extension).c_str(), rowIndices.size());
      }
    #else
      dumpToMem(b, pathB.append(extension).c_str(), LCSR->getRows()-1);
      dumpToMem(LCSR->getVals(), pathVals.append(extension).c_str(), LCSR->getNNZs());
      dumpToMem(LCSR->getRowPtr(), pathRowPtr.append(extension).c_str(), LCSR->getRowPtr().size());
      dumpToMem(LCSR->getColIdx(), pathParents.append(extension).c_str(), LCSR->getColIdx().size());
      dumpToMem(rowIndices, pathRowIndices.append(extension).c_str(), rowIndices.size());

      cout << "Matrix dimensions:\n";
      cout << "rowPtr: " << LCSR->getRowPtr().size() << "\n";
      cout << "colIdx: " << LCSR->getColIdx().size() << "\n";
      cout << "values: " << LCSR->getVals().size() << "\n";
    #endif
  } else if(dumpCSR) {
    // dump CSR format array for external lib code
    dumpToMem(LCSR->getVals(), pathVals.c_str(), LCSR->getNNZs());
    dumpToMem(LCSR->getRowPtr(), pathRowPtr.c_str(), LCSR->getRowPtr().size());
    dumpToMem(LCSR->getColIdx(), pathParents.c_str(), LCSR->getColIdx().size());

    cout << "Matrix dimensions:\n";
    cout << "rowPtr: " << LCSR->getRowPtr().size() << "\n";
    cout << "colIdx: " << LCSR->getColIdx().size() << "\n";
    cout << "values: " << LCSR->getVals().size() << "\n";

    #ifdef REWRITE_ENABLED
      vector<int> transformedRowPtr, transformedColIdx;
      vector<double> transformedValues;
      buildTransformedMatrix(transformedRowPtr, transformedColIdx, transformedValues);

      dumpToMem(transformedValues, pathVals.append(extension_TR).c_str(), transformedValues.size());
      dumpToMem(transformedRowPtr, pathRowPtr.append(extension_TR).c_str(), transformedRowPtr.size());
      dumpToMem(transformedColIdx, pathParents.append(extension_TR).c_str(), transformedColIdx.size());
    #endif
  }
  auto t8 = std::chrono::steady_clock::now();
  cout << "chrono dump data: " << std::chrono::duration_cast<std::chrono::duration<double>>(t8 - t7).count()*1000 << "\n";


  return returnVal;
}

#ifdef REWRITE_ENABLED
void Rewrite::buildTransformedMatrix(vector<int>& transformedRowPtr, vector<int>& transformedColIdx, vector<double>& transformedValues){
  DAG& dag = analyzer->getDAG();
  int rows = matrixCSR->getNumOfRows();
  vector<vector<double>>& values = analyzer->getValues();

  transformedRowPtr.push_back(0);

  for(int i = 0; i < rows; i++) {
    if(rewritingStrategy->isRewritten(i)) {
      vector<double>& rewrittenValues = rewrittenRowValues[i];
      transformedValues.insert(transformedValues.end(), rewrittenValues.begin(), rewrittenValues.end());
      transformedValues.push_back(1.0); // own value
    } else {
      vector<double>& rowValues = values[i];
      transformedValues.insert(transformedValues.end(), rowValues.begin(), rowValues.end());
    }

    vector<int>& parents = rewritingStrategy->isRewritten(i) ? rewrittenRowParents[i] : dag[i].first;
  //  cout << "row: " << i << " parents size: " << parents.size() << " last rowptr: " << transformedRowPtr.back() << "\n";
    transformedRowPtr.push_back(transformedRowPtr.back() + parents.size() + 1);

    transformedColIdx.insert(transformedColIdx.end(), parents.begin(), parents.end());
    transformedColIdx.push_back(i);
  }

  transformedColIdx.shrink_to_fit();
  transformedRowPtr.shrink_to_fit();
  transformedValues.shrink_to_fit();

  cout << "Transformed matrix dimensions:\n";
  cout << "rowPtr: " << transformedRowPtr.size() << "\n";
  cout << "colIdx: " << transformedColIdx.size() << "\n";
  cout << "values: " << transformedValues.size() << "\n";

}
#endif
