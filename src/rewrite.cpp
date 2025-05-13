#include <ios>
#include <iterator>
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
void Rewrite::writeFunc(std::ofstream &stream, vector< vector<int> >& tracker, int size, int maxNumOfThreads, int headerCounter) {
  if(!stream.is_open())
    std::cout << __LINE__ << " Cannot open output file!\n";

  vector<int>& partStarts   = tracker[0];
  vector<int>& threadCounts = tracker[1];

  stream << "#include \"calculators" + to_string(headerCounter) + ".h\"\n\n";

  int sign;
  if((headerCounter+1) * TABLE_SIZE > signatureLevel.size())
    sign = accumulate(signatureLevel.begin() + headerCounter * TABLE_SIZE, signatureLevel.end(), 0);
  else
    sign = accumulate(signatureLevel.begin() + headerCounter * TABLE_SIZE, signatureLevel.begin() + ((headerCounter+1) * TABLE_SIZE), 0);
  if(sign)
    stream << "void run" << headerCounter << "(double* x, double* b, int* parents, double* values, int* rowPtr, int* rowIndices) {\n";
  else
    stream << "void run" << headerCounter << "(double* x) {\n";

  int counter = partStarts[0];
  if (maxNumOfThreads > 1)
    stream << "\n#pragma omp parallel num_threads(" << maxNumOfThreads << ")\n{\n";

  int levelCounter = headerCounter * TABLE_SIZE, threadCounter = 0;
  for(int j = 0 ; j < size ; j++) {
    if(threadCounts[j] > 1) {
      stream << "\n";
      stream << "  #pragma omp sections // " << j << ", " << to_string(threadCounts[j]) << "\n"
             << "  {\n";

      for(int i = counter ; i < counter+threadCounts[j] ; i++) {
        stream << "    #pragma omp section\n";
        if(signature[levelCounter][threadCounter] == 0)
          stream << "    { calculate" << i << "(x); }\n";
        else
          stream << "    { calculate" << i << "(x, b, parents, values, rowPtr, rowIndices); }\n";

        threadCounter++;
      }
      stream << "  }\n";
    } else {
       if(threadCounts[j] == 1) {
        if(signature[levelCounter][threadCounter] == 0)
         stream << "  #pragma omp single // " << j << "\n"
                << "    { calculate" << counter << "(x); }\n";
       else
         stream << "  #pragma omp single // " << j << "\n"
                << "    { calculate" << counter << "(x, b, parents, values, rowPtr, rowIndices); }\n";
       } else {
         auto it = mergedLevels.find(levelCounter);
  
         if(it != mergedLevels.end()) {
           if((it->second)[1] == 1)
             stream << "  #pragma omp single // " << j << "\n"
                    << "  calculate" << counter << "(x); \n";
           else
             stream << "  #pragma omp single // " << j << "\n"
                    << "  calculate" << counter << "(x, b, parents, values, rowPtr, rowIndices); \n";
         }
       }
         
       stream << "\n";
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
       stream << "#define FILEPATH_PARENTS \"/tmp/" + fileName + "_parents_TR.bin\"\n"
              << "#define FILEPATH_ROWPTR \"/tmp/" + fileName + "_rowPtr_TR.bin\"\n"
              << "#define FILEPATH_VALUES \"/tmp/" + fileName + "_vals_TR.bin\"\n"
              << "#define PARENTS_SIZE " +  to_string(parentsSize) + "\n\n"; 
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
       stream << "#define FILEPATH_PARENTS \"/tmp/" + fileName + "_parents.bin\"\n"
              << "#define FILEPATH_ROWPTR \"/tmp/" + fileName + "_rowPtr.bin\"\n"
              << "#define FILEPATH_VALUES \"/tmp/" + fileName + "_vals.bin\"\n"
              << "#define PARENTS_SIZE " +  to_string(vals) + "\n\n"; 
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
         << "    if(((fabs(1.000000 - x[i])/1.000000) >= 1e-2) && fabs(x[i]) >= 0) {\n"
         << "      errCnt++;\n"
         << "      //printf(\"x[%d]: %.5f\\n\",i,x[i]);\n"
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
    std::cout << __LINE__ << " Cannot open output file!\n";

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
// when calling: put -1 as key for constant to multipliers
// rewritingMultiplicant is 1 (this is for rewritten)
// sign: + is true, - is false
void Rewrite::calculateMultiplicants(int row, vector<double> &b, set<int>& rewritten, map<int,double>& multipliers, double rewritingMultiplicant, bool sign) {
  DAG& dag = analyzer->getDAG();
  map<int, vector<int>>& oriParents = analyzer->getOriParents();
  vector<int>& parents = oriParents[row];

  map<int, vector<double>>& oriRowValues = analyzer->getOriRowValues();
  vector<double>& rowValues = oriRowValues[row];

  if(!parents.empty()) {
    for(int j = 0 ; j < parents.size() ; j++) {
      auto it = rewritten.find(parents[j]);
      if(it != rewritten.end()) {
  /*      if(parents[j] == 0 || parents[j] == 1 || parents[j] == 2 || parents[j] == 3) {
          cout << __LINE__ << " rewriting " << parents[j] << "\n";
          cout << "sign:" << !sign << "\nrewritingMultiplicant: " << rewritingMultiplicant << " * rowValues[" << j << "] (which is " << rowValues[j] << ") /" << rowValues.back() << " : " << rewritingMultiplicant * rowValues[j] /rowValues.back() << "\n";
        }*/
        calculateMultiplicants(parents[j], b, rewritten, multipliers, rewritingMultiplicant * (rowValues[j] / rowValues.back()), !sign);  
      } else {
   /*     if(parents[j] == 0 || parents[j] == 1 || parents[j] == 2 || parents[j] == 3) 
          cout << __LINE__ << " sign:" << sign << " multiplicant[" << parents[j] << "] get (+/-):" <<  "\nrewritingMultiplicant: " << rewritingMultiplicant << " * rowValues[" << j << "] (which is " << rowValues[j] << ") /" << rowValues.back() << " : " << rewritingMultiplicant * rowValues[j] /rowValues.back() << "\n";*/
        if(sign == true)
          multipliers[parents[j]] +=  rewritingMultiplicant * (rowValues[j] /rowValues.back());
        else 
          multipliers[parents[j]] -=  rewritingMultiplicant * (rowValues[j] /rowValues.back());
      }
    }

/*    if(row == 3 || row == 0 || row == 1 || row == 2 || row == 4)
      cout << __LINE__ << " sign:" << sign << " multiplicant[-1] get (+/-):" <<  "\nrewritingMultiplicant: " << rewritingMultiplicant << " * b[" << row << "] (which is " << b[row] << ") /" << rowValues.back() << " : " << rewritingMultiplicant * (b[row] / rowValues.back()) << "\n"; */

    if(sign == true)
      multipliers[-1] += rewritingMultiplicant * (b[row] / rowValues.back());
    else
      multipliers[-1] -= rewritingMultiplicant * (b[row] / rowValues.back());
  } else { // parents empty (rewriting to level 0)
/*        if(row == 3 || row == 0 || row == 1 || row == 2 || row == 4)
          cout << __LINE__ << " sign:" << sign << " multiplicant[-1] get (+/-):" <<  "\nrewritingMultiplicant: " << rewritingMultiplicant << " * b[" << row << "] (which is " << b[row] << ") /" << rowValues.back() << " : " << rewritingMultiplicant * (b[row] / rowValues.back()) << "\n";*/

    if(sign == true)
      multipliers[-1] += rewritingMultiplicant * (b[row] / rowValues.back());
    else
      multipliers[-1] -= rewritingMultiplicant * (b[row] / rowValues.back());
  }

/*    cout << "multipliers:\n";
    for(auto& multip : multipliers)
      cout << multip.first << ": " << multip.second << "\n";
    cout << "\n";
  }
*/
}

void Rewrite::updateRewrittenRow(int row, vector<double>& b, set<int>& rewritten, map<int,double>& multipliers) {
  DAG& dag = analyzer->getDAG();
  vector<vector<double>>& values = analyzer->getValues();
  vector<double>& rowValues = values[row];
  vector<int>& parents = dag[row].first;

  b[row] = multipliers[-1];
  multipliers.erase(-1);
  parents.clear();
  rowValues.clear();

  for(auto it = multipliers.begin() ; it != multipliers.end(); it++) {
    if(rewritten.find(it->first) == rewritten.end()) {
      parents.push_back(it->first);
      rowValues.push_back(it->second);
    }
  }

  // now itself will be 1.00 since rowValues come as divided by this value anyway (Lx = b)
  rowValues.push_back(1.00);

  if(parents.empty())
    b[row] = 1.00;
}

#endif

void Rewrite::dumpLoop(int beginIndex, int endIndex, vector<int>& loopedRows, std::ofstream& stream) {
  if(endIndex-beginIndex > 1) {
    stream << "  for(int i = " << beginIndex << " ; i < " << endIndex << "; i++) {\n"
           << "    int row = rowIndices[i];\n";
  } else // single looped row
    stream << "    { int row = " << loopedRows[0] << ";\n";

  stream << "      double xi = 0;\n"
         << "      for (int j = rowPtr[row]; j < rowPtr[row+1]-1; j++)\n"
         << "        xi += values[j] * x[parents[j]];\n\n"
         << "      x[row] = (b[row]-xi)/values[rowPtr[row+1]-1]; }\n\n";
}

int Rewrite::dumpUnrolled(vector<int>& unrolledRows, vector<double>& b, std::ofstream& stream) {
  int flops = 0;
  DAG& dag = analyzer->getDAG();
  vector<vector<double>>& values = analyzer->getValues();
 
  for(auto& row : unrolledRows) {
    vector<int>& parents = dag[row].first;
    vector<double>& rowValues = values[row];

//    cout << "row: " << row << " num. parents: " << parents.size() << "\n";
    if(!parents.empty()) {
      flops += (parents.size() << 1) + 1;
      #ifdef REWRITE_ENABLED
        if(rewritingStrategy->isRewritten(row)) flops--;
      #endif
    }

    if(parents.empty()) {
      #ifdef REWRITE_ENABLED
        if(rewritingStrategy->isRewritten(row))
          //stream << std::fixed << std::setprecision(6) << "  x[" << row << "] = " << b[row] << ";\n";
          stream << std::fixed << std::setprecision(6) << "  x[" << row << "] = " << rowValues.back() << ";\n";
        else
      #endif
        stream << std::fixed << std::setprecision(6) << "  x[" << row << "] = " << b[row] / rowValues.back() << ";\n";
    } else if(parents.size() == 1) {
      #ifdef REWRITE_ENABLED
        if(rewritingStrategy->isRewritten(row))
          stream << std::fixed << std::setprecision(6) << "  x[" << row << "] = " << b[row] << " - " << rowValues[0] << " * x[" << parents[0] << "];\n";
        else
      #endif
          stream << std::fixed << std::setprecision(6) << "  x[" << row << "] = (" << b[row] << " - " << rowValues[0] << " * x[" << parents[0] << "])/" << rowValues.back() << ";\n";
    } else if(parents.size() == 2) {
      #ifdef REWRITE_ENABLED
        if(rewritingStrategy->isRewritten(row)) {
          stream << std::fixed << std::setprecision(6) << "  x[" << row << "] = " << b[row] << "-(" \
                 << rowValues[0] << " * x[" << parents[0] << "] + " \
                 << rowValues[1] << " * x[" << parents[1] << "]);\n";
        } else
      #endif
          stream << std::fixed << std::setprecision(6) << "  x[" << row << "] = (" << b[row] << "-(" \
                 << rowValues[0] << " * x[" << parents[0] << "] + " \
                 << rowValues[1] << " * x[" << parents[1] << "]))/" << rowValues.back() << ";\n";
    } else if(parents.size() == 3) {
      #ifdef REWRITE_ENABLED
        if(rewritingStrategy->isRewritten(row)) {
          stream << std::fixed << std::setprecision(6) << "  x[" << row << "] = " << b[row] << "-(" \
                 << rowValues[0] << " * x[" << parents[0] << "] + " \
                 << rowValues[1] << " * x[" << parents[1] << "] +" \
                 << rowValues[2] << " * x[" << parents[2] << "]);\n";
        } else
      #endif
          stream << std::fixed << std::setprecision(6) << "  x[" << row << "] = (" << b[row] << "-(" \
                 << rowValues[0] << " * x[" << parents[0] << "] + " \
                 << rowValues[1] << " * x[" << parents[1] << "] + " \
                 << rowValues[2] << " * x[" << parents[2] << "]))/" << rowValues.back() << ";\n";
    } else if(parents.size() == 4) {
      #ifdef REWRITE_ENABLED
        if(rewritingStrategy->isRewritten(row)) {
          stream << std::fixed << std::setprecision(6) << "  x[" << row << "] = " << b[row] << "-(" \
                 << rowValues[0] << " * x[" << parents[0] << "] +" \
                 << rowValues[1] << " * x[" << parents[1] << "] +" \
                 << rowValues[2] << " * x[" << parents[2] << "] +" \
                 << rowValues[3] << " * x[" << parents[3] << "]);\n";
        } else
      #endif
          stream << std::fixed << std::setprecision(6) << "  x[" << row << "] = (" << b[row] << "-(" \
                 << rowValues[0] << " * x[" << parents[0] << "] +" \
                 << rowValues[1] << " * x[" << parents[1] << "] +" \
                 << rowValues[2] << " * x[" << parents[2] << "] +" \
                 << rowValues[3] << " * x[" << parents[3] << "]))/" << rowValues.back() << ";\n";
    }
  }

  stream << "\n";
  return flops;
}

#ifdef REWRITE_ENABLED
void Rewrite::rewriteInLoop(vector<int>& loopedRows, vector<double> &b) {
  int flops = 0; double total_time = 0.0;
  DAG& dag = analyzer->getDAG();

  vector<int>::iterator it = loopedRows.begin();
  while(it != loopedRows.end()) {
      if(rewritingStrategy->isRewritten(*it)) {
        rewriteRow(*it, b);
      }
    it++;
  }
}
#endif

int Rewrite::writePart(int levelNum, int rowStartIndex, int rowEndIndex, int levelPart, vector<double> &b, std::ofstream& stream, int signDump, int* sigType) {
  int flops = 0, beginLevel, endLevel, loopEnd, loopBegin = rowIndices.size();
  bool loopedRowsExist = false;

  vector<vector<int>>& levelTable = analyzer->getLevelTable();
  vector<int>& level = levelTable[levelNum];
  vector<int> unrolledRows, loopedRows;
      
  if(levelLookUp.find(levelNum) != levelLookUp.end()) {
    beginLevel = levelLookUp[levelNum];
    endLevel = mergedLevels[beginLevel][2];

    int prev = -1;
    list<int> removedLevels;
    for(int i = beginLevel; i <= endLevel; i++) {
      list<vector<int>>& loopedUnrolled = mergedLevelsDist[i];

      #ifdef REWRITE_ENABLED
        ch_start = std::chrono::steady_clock::now();
          rewriteInLoop(loopedUnrolled.front(), b);  // loopedRows
        ch_end = std::chrono::steady_clock::now();
        ch_ttime += (ch_end-ch_start);
      #endif
      analyzer->separateRows(loopedUnrolled.front(), loopedUnrolled.back());

      if(!loopedUnrolled.front().empty()) {
        loopedRowsExist = true;
        rowIndices.insert(rowIndices.end(), loopedUnrolled.front().begin(), loopedUnrolled.front().end());
      }

      // found prev: loops exist AND no unrolled
      if(prev == -1 && !loopedUnrolled.front().empty() && loopedUnrolled.back().empty()) {
        prev = i;
        // update the ending level since 
        mergedLevels[beginLevel][2] = prev;
      }

      // merge loops to prev level's loops if there's no prev unrolled
      if(prev > -1 && prev < i && !loopedUnrolled.front().empty()) {
        list<vector<int>>& prevLoopedUnrolled = mergedLevelsDist[prev];

        vector<int>& ref = loopedUnrolled.front();   // curr loops
        vector<int> loopsMerged(prevLoopedUnrolled.front());        // new vector out of prev loops since cannot resize

        loopsMerged.insert(loopsMerged.end(), ref.begin(), ref.end());

        prevLoopedUnrolled.begin()->clear();
        prevLoopedUnrolled.erase(prevLoopedUnrolled.begin());  // remove old prev level loops & insert extended one
        prevLoopedUnrolled.push_front(loopsMerged);
        ref.clear();

        if(!loopedUnrolled.back().empty())
          prev = -1;

        if(loopedUnrolled.front().empty() && loopedUnrolled.back().empty())
          removedLevels.push_back(i);
      }
    }

    for(auto& level: removedLevels)
      mergedLevelsDist.erase(level);

  } else {
    loopedRows.insert(loopedRows.begin(), level.begin() + rowStartIndex, level.begin() + rowEndIndex);

    #ifdef REWRITE_ENABLED
      ch_start = std::chrono::steady_clock::now();
        rewriteInLoop(loopedRows, b);
      ch_end = std::chrono::steady_clock::now();
      ch_ttime += (ch_end-ch_start);
    #endif
    analyzer->separateRows(loopedRows, unrolledRows);

    if(!loopedRows.empty()) {
      loopedRowsExist = true;
      rowIndices.insert(rowIndices.end(), loopedRows.begin(), loopedRows.end());
    }
  }

  if(loopedRowsExist) {
      analyzer->setSingleLoopRows(true);

    if(rowStartIndex == 0) {
      signatureLevel.push_back(1);
      signature.push_back(vector<int>());
    } else
      signatureLevel[levelNum] = 1;

    // update signature of starting merged level if its signature is type 1
    if(signDump != 0) {
      if(mergedLevels.find(beginLevel) != mergedLevels.end()) {
        if(mergedLevels[beginLevel][1] == 1)
          mergedLevels[beginLevel][1] = 2;
      }
    }

    signature[levelNum].push_back(1);  // if we OR signature[levelNum] == signatureLevel[levelNum]

    *sigType = 2;
    if(signDump == 2)  // beginning of single-threaded OR multi-threaded level
      stream << "  void calculate" << levelPart << "(double* x, double* b, int* parents, double* values, int* rowPtr, int* rowIndices) {\n";

  } else {
    if(rowStartIndex == 0) {
      signatureLevel.push_back(0);
      signature.push_back(vector<int>());
    }
      
    signature[levelNum].push_back(0);

    if(signDump == 2)  // beginning of single-threaded OR multi-threaded level
      stream << "  void calculate" << levelPart << "(double* x) {\n";
  }

  ch_start2 = std::chrono::steady_clock::now();

  int index = 0;
  if(levelLookUp.find(levelNum) != levelLookUp.end()) {
    endLevel = mergedLevels[beginLevel][2];
    map<int, list<vector<int>>>::iterator it = mergedLevelsDist.find(beginLevel);
    for(; it != mergedLevelsDist.end(); ++it) {
      int i = it->first; // get the level
      list<vector<int>>& loopedUnrolled = mergedLevelsDist[i];

      if(!loopedUnrolled.front().empty()) {
        loopEnd = loopBegin + loopedUnrolled.front().size();
        dumpLoop(loopBegin, loopEnd, loopedUnrolled.front(), stream);
        loopBegin = loopEnd;
     
        if(levelNum > 0)
          flops += analyzer->recalculateFLOPSFor(loopedUnrolled.front());
      }

      if(!loopedUnrolled.back().empty())
        flops += dumpUnrolled(loopedUnrolled.back(), b, stream);

      if(i == endLevel)
        break;
    }
  } else { 

    if(loopedRowsExist) {
      loopEnd = loopBegin + loopedRows.size();
      dumpLoop(loopBegin, loopEnd, loopedRows, stream);
      loopBegin = loopEnd;
        
      if(levelNum > 0)
        flops += analyzer->recalculateFLOPSFor(loopedRows);
    }

    if(!unrolledRows.empty())
      flops += dumpUnrolled(unrolledRows, b, stream);
  }
  
  stream << "}\n\n";

  ch_end2 = std::chrono::steady_clock::now();
  ch_ttime2 += (ch_end2 - ch_start2);

  return flops;
}

#ifdef REWRITE_ENABLED
void Rewrite::rewriteRow(int row, vector<double>& b) {
    if(!rewritingStrategy->isRewritten(row))
      cout << "row: " << row << " doesnt appear in rewrittenStrategy's list\n"; 

    set<int>& rewritten = rewritingStrategy->getRewritingMapOf(row);
    map<int,double> multipliers;
    for(auto& row : rewritten)
      multipliers[row] = 0;
    multipliers[-1] = 0;
    double rewritingMultiplicant = 1;

    calculateMultiplicants(row, b, rewritten, multipliers, rewritingMultiplicant, true);
    updateRewrittenRow(row, b, rewritten, multipliers);
}
#endif

int Rewrite::writeMakefile(int execBlockCnt) {
  std::ofstream stream(fileName + "/Makefile");
  if(!stream.is_open()) {
    std::cout << __LINE__ << " Cannot open output file!\n";
    return -1;
  }

  stream << "CXX      := clang\n"
  //<< "CXXFLAGS := -O3 -ffast-math -funroll-loops -march=native -I/kuacc/apps/llvm-omp/include\n"
 // << "CXXFLAGS := -O3 -flto=thin -ffast-math -march=native \n"
  << "CXXFLAGS := -O3 -flto=thin  -march=native -fopenmp\n"
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

int Rewrite::writeHeader(int headerCounter, int partCounterStart, int partCounterEnd, vector<int>& threadCounts) {
  std::ofstream stream(fileName + "/calculators" + to_string(headerCounter) + ".h");
  if(!stream.is_open()) {
    std::cout << __LINE__ << " Cannot open output file!\n";
    return -1;
  }

  int levelCounter = headerCounter * TABLE_SIZE;

  // if all signatures of a level are consumed, current level is finished
  // the latter counts in threadCounts, levelCounter doesn't start @0
  int counter = 0, cntLevelStart = 0;

  stream << "#pragma once\n";
  for(int i = partCounterStart; i < partCounterEnd; i++) {
     if(threadCounts[cntLevelStart] >= 1) {
        if(signature[levelCounter][counter] == 0)
          stream << "void calculate" << to_string(i) << "(double* x);\n";
        else
          stream << "void calculate" << to_string(i) << "(double* x, double* b, int* parents, double* values, int* rowPtr, int* rowIndices);\n";
     } else {
       auto it = mergedLevels.find(levelCounter);

       if(it != mergedLevels.end())
        if((it->second)[1] == 1)
          stream << "void calculate" << to_string(i) << "(double* x);\n";
        else
          stream << "void calculate" << to_string(i) << "(double* x, double* b, int* parents, double* values, int* rowPtr, int* rowIndices);\n";
     }

    counter++;

    if(counter == signature[levelCounter].size()) {
      levelCounter++;
      cntLevelStart++;
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

    if(workloadPerThread < LOWER_BOUND) {
      numThreads = toBeBalanced+LOWER_BOUND-1;
      numThreads /= LOWER_BOUND;
      workloadPerThread = toBeBalanced+numThreads-1;
      workloadPerThread /= numThreads;
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
#endif

void Rewrite::mergeSingleThreadedLevels() {
  vector<vector<int>>& levelTable = analyzer->getLevelTable();
  auto& flopsPerLevel = analyzer->getFlopsPerLevel();
  int* levels = analyzer->getLevels();
  int numThreads = 1, workloadPerThread = 1;
  vector<int> emptyLevels;

  for(int i = 0 ; i < flopsPerLevel.size()-1; i++) { // cannot merge the last level
    vector<int>& level = levelTable[i];
          
    int numThreads = (level.size() > 1) ? balanceLevel(level.size(), workloadPerThread) : 1;
    if(numThreads == 1) {
      if(mergedLevels.find(i) == mergedLevels.end()) {
        mergedLevels[i] = vector<int>(3, 0);
        mergedLevels[i][0] = flopsPerLevel[i];   // total cost
        mergedLevels[i][1] = 1;   // default to short type sig
        mergedLevels[i][2] = -1;  // ending level

        if(mergedLevels.size() == 1)
          levelLookUp[i] = i;   // new starting level # (after erasing empty levels it changes)

        mergedLevelsDist[i] = list<vector<int>>();
        vector<int>& level = levelTable[i];
        list<vector<int>>& loopedUnrolled = mergedLevelsDist[i];
        vector<int> rows(level.begin(), level.end());
        loopedUnrolled.push_back(rows);
        loopedUnrolled.push_back(vector<int>());

        int currLevel = i+1;
        while(currLevel < flopsPerLevel.size()) {
          vector<int>& level = levelTable[currLevel];
          numThreads = (level.size() > 1) ? balanceLevel(level.size(), workloadPerThread) : 1;
          vector<int>& currMerged = mergedLevels[i];

          // these are the regions that rewriting method couldnt be applied,
          // hence we fall back to the traditional sptrsv even without the levels
          // the classical looped calculation + unrolling
          if(numThreads == 1) { // merge consecutive single threaded levels until multi-threaded one is reached
            vector<int>& levelBegin = levelTable[i];
            vector<int>& levelEnd = levelTable[currLevel];
            for(auto& row : levelEnd)
              levels[row] = i;

            levelBegin.insert(levelBegin.end(),levelEnd.begin(),levelEnd.end());
            
            currMerged[0] += flopsPerLevel[currLevel];
            emptyLevels.push_back(currLevel);

            list<vector<int>>& loopedUnrolled = mergedLevelsDist[currLevel];
            vector<int> rows(level.begin(), level.end());
            loopedUnrolled.push_back(rows);
            loopedUnrolled.push_back(vector<int>());

            currLevel++;

            // NOT USED: MERGE_DEPTH (default 10) is reached or last level merged
            //if((currLevel == i + MERGE_DEPTH) || (currLevel == flopsPerLevel.size())) {
            if((currMerged[0] >= analyzer->getALC()) || (currLevel == flopsPerLevel.size())) {
              currMerged[2] = currLevel-1;
              flopsPerLevel[i] = currMerged[0];

              // update new starting level
              if(mergedLevels.size() > 1) {
                map<int, vector<int>>::reverse_iterator rit = mergedLevels.rbegin(); rit++;
                int start = rit->first;  // 0
                int end = (rit->second)[2]; // 455
                int newStartLevel = levelLookUp.rbegin()->first + (i - end);
                levelLookUp[newStartLevel] = i;
              }

              break;
            }
          } else { // multi-threaded level encountered
            if(currLevel != i+1) { 
              currMerged[2] = currLevel-1;
              currLevel++;
              flopsPerLevel[i] = currMerged[0];

              // update new starting level
              if(mergedLevels.size() > 1) {
                map<int, vector<int>>::reverse_iterator rit = mergedLevels.rbegin(); rit++;
                int start = rit->first;  // 0
                int end = (rit->second)[2]; // 455
                int newStartLevel = levelLookUp.rbegin()->first + (i - end);
                levelLookUp[newStartLevel] = i;
              }
            } 
            // we found only 1 single-threaded level:
            // starting & ending levels are the same: REMOVE IT
            else {
              mergedLevels[i].clear();
              mergedLevels.erase(i);
              mergedLevelsDist.erase(i);
              levelLookUp.erase(i);
            }

            break; 
          }
        } // while

        i = currLevel-1;
      }
    } // numThreads == 1
  }

  // needed since emptyLevels are erased using updateWithEmptyLevels()
  reverse(emptyLevels.begin(), emptyLevels.end());
  analyzer->updateWithEmptyLevels(emptyLevels);

  #ifdef REPORT
    cout << "empty levels:\n";
    for(auto& i : emptyLevels)
      cout << i << ", ";
    cout << "\n";

    cout << "\nmergedLevels:\n";
    for(map<int,vector<int>>::iterator iter = mergedLevels.begin(); iter != mergedLevels.end(); iter++) {
      vector<int>& ref = iter->second;
      cout << "starting level: " << iter->first << "\nending level: " << ref[2] << "\n\n";
    }

  cout << "\n\nmergedLevelsDist:\n";
  for(map<int,list<vector<int>>>::iterator iter = mergedLevelsDist.begin(); iter != mergedLevelsDist.end(); iter++) {
    list<vector<int>>& ref = iter->second;
    cout << "level: " << iter->first << "\n";

    cout << "looped rows: " << ref.front().size() << "\n";
    for(auto& row : ref.front())
      cout << row << ", ";
    cout << "\nunrolled rows: " << ref.back().size() << "\n";
    for(auto& row : ref.back())
      cout << row << ", ";
    cout << "\n";
  }
  #endif
}

int Rewrite::rewriteExecutor(vector<double> &b) {
  Part *L = matrixCSR->getL();
  analyzer->setSingleLoopRows(false);
  writeUtil();

  vector<vector<int>>& levelTable = analyzer->getLevelTable();
  auto& flopsPerLevel = analyzer->getFlopsPerLevel();
  int partCounter = 0, headerCounter = 0, partCounterStart = 0, maxNumOfThreads = 1;

  mergeSingleThreadedLevels();

  // Tracker keeps start of each part(tracker[0]) & num of thread(tracker[1]), part: workload of a thread
  vector< vector<int> > tracker(2, vector<int>(TABLE_SIZE,0));

  double chrono_par = 0.0, chrono_serial = 0.0;

  #ifdef PAR
  // parallel file generation (runX.c, calculatorsX.h)
  // define num of threads to be used so that each get a TABLE_SIZE, max num of threads can be NUM_THREADS
  int parFileGenThreadNum = flopsPerLevel.size()/TABLE_SIZE + ((flopsPerLevel.size() % TABLE_SIZE) > 0);

  //printf("flopsPerLevel.size(): %d, parFileGenThreadNum:%d\n", flopsPerLevel.size(), parFileGenThreadNum);

  int th_cnt = 0, numParts = parFileGenThreadNum;
  int* partCounterStarts = NULL, *sizes = NULL, *parts = NULL;
  map<int, vector<vector<int>> >* trackerAll = NULL;

  auto t_start = std::chrono::steady_clock::now();
  if(parFileGenThreadNum > 1) {
    partCounterStarts = (int*) malloc((parFileGenThreadNum + 1) * sizeof(int));
    sizes = (int*)malloc(parFileGenThreadNum * sizeof(int));
    partCounterStarts[0] = 0;
    trackerAll = new map<int, vector<vector<int>> >();
  }

  auto t_end = std::chrono::steady_clock::now();
  chrono_par += std::chrono::duration<double>(t_end - t_start).count();

  if(parFileGenThreadNum > NUM_THREADS) 
    parFileGenThreadNum = NUM_THREADS;
  #endif

  int sigType = 1; // type1 : short one, type2: long one
  for(int j = 0 ; j < flopsPerLevel.size(); j++) {
  int numThreads = 1, workloadPerThread = 1, flops = 0;
    vector<int>& level = levelTable[j];
    vector<int> flopsLevel;

    #ifdef BALANCE_ROWS
      numThreads = (levelLookUp.find(j) == levelLookUp.end()) ?  balanceLevel(level.size(), workloadPerThread): 1;
      #ifdef REPORT
        cout << "level " << j << " numThreads: " << numThreads << "\n";
      #endif
    #endif
      
    if(numThreads > maxNumOfThreads)
      maxNumOfThreads = numThreads;

    #ifdef BALANCE_ROWS
    std::ofstream streamLevel;
    streamLevel.open(fileName + "/calculate" + to_string(partCounter) + ".c");

    if(!streamLevel.is_open())
      std::cout << __LINE__ << " Cannot open output file!\n";

    // merged = 0 (not merged), =1 (merged), =2 (merged but beginning level)
    if(levelLookUp.find(j) != levelLookUp.end()) {
      //flops += writePart(0, 0, level.size(), partCounter++, b, streamLevel, 2, &sigType);
      flopsPerLevel[j] = writePart(j, 0, level.size(), partCounter++, b, streamLevel, 2, &sigType);
    } else {
      if(numThreads == 1) // single threaded 1 level (not merged)
        flopsPerLevel[j] = writePart(j, 0, level.size(), partCounter++, b, streamLevel, 2, &sigType);
      else { 
        flops += writePart(j, 0, workloadPerThread, partCounter++, b, streamLevel, 2, &sigType);
        for(int i = 1; i < numThreads-1; i++) {
          flops += writePart(j, i * workloadPerThread, (i+1) * workloadPerThread, partCounter++, b, streamLevel, 2, &sigType);
        }
    
        flops += writePart(j, (numThreads-1) * workloadPerThread, level.size(), partCounter++, b, streamLevel, 2, &sigType);
        flopsPerLevel[j] = flops;
      }
    }

    streamLevel.close();
    #endif
      
    vector<int>& partStarts =   tracker[0];
    vector<int>& threadCounts = tracker[1];
    partStarts[(j % TABLE_SIZE)] = partCounter-numThreads;
    threadCounts[(j % TABLE_SIZE)] = numThreads;

    #ifdef PAR
    if(parFileGenThreadNum > 1) {
//      printf("line: %d j: %d num of levels: %d\n", __LINE__, j, flopsPerLevel.size());
      auto t15 = std::chrono::steady_clock::now();
      if(((j % TABLE_SIZE) == (TABLE_SIZE-1)) || (j == flopsPerLevel.size()-1)) {
        partCounterStarts[th_cnt+1] = partCounter;
        sizes[th_cnt] = (j == flopsPerLevel.size()-1) ? (flopsPerLevel.size() % TABLE_SIZE) : TABLE_SIZE;
        vector<vector<int>> tckr(tracker);
        (*trackerAll)[th_cnt] = tckr;
        th_cnt++;
      } 
      auto t16 = std::chrono::steady_clock::now();
      chrono_par += std::chrono::duration<double>(t16 - t15).count();
    } else {
    #endif
      auto t13 = std::chrono::steady_clock::now();
      int size = (j == flopsPerLevel.size()-1) ? (flopsPerLevel.size() % TABLE_SIZE) : TABLE_SIZE; 
  
      if(((j % TABLE_SIZE) == (TABLE_SIZE-1)) || (j == flopsPerLevel.size()-1)) {
        std::ofstream stream(fileName + "/run" + to_string(headerCounter) + ".c");
        writeHeader(headerCounter, partCounterStart, partCounter, threadCounts);
 //       printf("SINGLE: size:%d, part:%d, headCounter:%d\n", size, part, headerCounter);
        writeFunc(stream, tracker, size, maxNumOfThreads, headerCounter);
        headerCounter++;
        partCounterStart = partCounter;
      }
    
      auto t14 = std::chrono::steady_clock::now();
      chrono_serial += std::chrono::duration<double>(t14 - t13).count();
    #ifdef PAR
    }
    #endif
  }

  #ifdef PAR
  auto t11 = std::chrono::steady_clock::now();
  if(parFileGenThreadNum > 1) {
    #pragma omp parallel for schedule(static,2) num_threads(parFileGenThreadNum)
    for(int i = 0; i < numParts; i++) {
      writeHeader(i, partCounterStarts[i], partCounterStarts[i+1], tracker[1]);
      std::ofstream stream(fileName + "/run" + to_string(i) + ".c");
      writeFunc(stream, (*trackerAll)[i], sizes[i], maxNumOfThreads, i);
    }

    free(sizes); free(partCounterStarts); delete trackerAll;
    sizes = NULL; partCounterStarts = NULL;
  }
  
  auto t12 = std::chrono::steady_clock::now();
  chrono_par += std::chrono::duration<double>(t12 - t11).count();

  if(parFileGenThreadNum > 1)
    cout << "* chrono dump runX.c calculatorsX.h: " << chrono_par * 1000 << "\n";
  else
 #endif
    cout << "chrono dump runX.c calculatorsX.h: " << chrono_serial * 1000 << "\n";

  #ifdef REPORT
    cout << "lookup table:\n";
    for(auto& item : levelLookUp)
      cout << item.first << ": " << item.second << "\n";
  
    cout << "leveltable, flopsPerLevel size: " << levelTable.size() << ", " << flopsPerLevel.size() << " NumOfLevels:" << analyzer->getNumOfLevels() << "\n";
  
    cout << "\nmergedLevels:\n";
    for(map<int,vector<int>>::iterator iter = mergedLevels.begin(); iter != mergedLevels.end(); iter++) {
      vector<int>& ref = iter->second;
      cout << "starting level: " << iter->first << "\nending level: " << ref[2] << "\n\n";
    }

  /*cout << "\n\nmergedLevelsDist:\n";
  for(map<int,list<vector<int>>>::iterator iter = mergedLevelsDist.begin(); iter != mergedLevelsDist.end(); iter++) {
    list<vector<int>>& ref = iter->second;
    cout << "level: " << iter->first << "\n";

    cout << "looped rows: " << ref.front().size() << "\n";
    for(auto& row : ref.front())
      cout << row << ", ";
    cout << "\nunrolled rows: " << ref.back().size() << "\n";
    for(auto& row : ref.back())
      cout << row << ", ";
    cout << "\n";
  }*/
  #endif
  
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
    printf("This is not a square matrix.\n");
    return -1;
  }

  vector<double> xRef(cols, 1.000000);
  vector<double> b(cols, 0.000000);
  
  DAG& dag = analyzer->getDAG();
  vector<vector<double>>& values = analyzer->getValues();

  Part* LCSR = matrixCSR->getL();

  /////////////////////////////// PRINT SOME BASIC STATS ///////////////////////////////
  #ifdef REPORT
    cout << "NUM_THREADS:" << NUM_THREADS << "\n";
    cout << "LOWER_BOUND: " << LOWER_BOUND << "\n";
    #ifdef BALANCE_ROWS
      cout << "BALANCE METHOD: BALANCE_ROWS\n";
    #else
      cout << "BALANCE METHOD: BALANCE_ROWS_CV\n";
    #endif
    cout << "TABLE_SIZE: " << TABLE_SIZE << "\n";
  #endif
  /////////////////////////////// PRINT SOME BASIC STATS ///////////////////////////////

  auto t5 = std::chrono::steady_clock::now();
  // fill in the B values, these are needed for writePart
  // B values for rewritten rows will be filled in by writePart
 
  #ifdef REWRITE_ENABLED
    map<int, vector<int>>& oriParents = analyzer->getOriParents();
    map<int, vector<double>>& oriRowValues = analyzer->getOriRowValues();
  #endif

  #ifdef PAR
  #pragma omp parallel for num_threads(NUM_THREADS)
  #endif
  for(int i = 0; i < rows; i++) {
  #ifdef REWRITE_ENABLED
    vector<int>& parents = oriParents[i];
    vector<double>& rowValues = oriRowValues[i];
  #else
    vector<int>& parents = dag[i].first;
    vector<double>& rowValues = values[i];
  #endif

    b[i] += rowValues.back() * xRef[i];
    for(int j = 0 ; j < parents.size() ; j++) {
      b[i] += rowValues[j] * xRef[parents[j]];
    }
  }


  auto t6 = std::chrono::steady_clock::now();
  #ifdef PAR
  cout << "* chrono calculate B: " << chrono::duration<double>(t6 - t5).count()*1000 << "\n";
  #else
  cout << "chrono calculate B: " << chrono::duration<double>(t6 - t5).count()*1000 << "\n";
  #endif

  // VERIFICATION OF XREF
//  cout << "VERIFICATION BEFORE\n";
//  verifyX(rows, dag, values, b, xRef);

  ch_ttime = std::chrono::duration<double>(0); 
  ch_ttime2 = std::chrono::duration<double>(0); 
  auto start = std::chrono::steady_clock::now();
  rewriteExecutor(b);
  auto end = std::chrono::steady_clock::now();
#ifdef REWRITE_ENABLED
  cout << "chrono rewrite: " << ch_ttime.count()*1000 << "\n";
#endif
  cout << "chrono calculate rewriteExecutor: " << chrono::duration<double>(end-start).count()*1000 << "\n";
  cout << "chrono dump calculateX.c: " << ch_ttime2.count()*1000 << "\n";

  int parentsSize = dumpDataToMem(b);
  auto& flopsPerLevel = analyzer->getFlopsPerLevel();

  start = std::chrono::steady_clock::now();
  std::ofstream stream(fileName + "/main.c");
  if(!stream.is_open())
    std::cout << __LINE__ << " Cannot open output file!\n";

  writeMain(stream, (int)(ceil((((double)flopsPerLevel.size())/TABLE_SIZE))), parentsSize);
  stream.close();
    
  writeMakefile(analyzer->getNumOfLevels());
  end = std::chrono::steady_clock::now();
  cout << "chrono dump main & makefile: " << chrono::duration<double>(end-start).count()*1000 << "\n";

  return 0;
}

int Rewrite::dumpDataToMem(vector<double>& b) {
  auto t7 = std::chrono::steady_clock::now();

  Part* LCSR = matrixCSR->getL();
  int returnVal = 0;

  #ifdef REWRITE_ENABLED
    string pathVals_TR ("/tmp/" + fileName + "_vals_TR.bin");
    string pathParents_TR ("/tmp/" + fileName + "_parents_TR.bin");
    string pathRowPtr_TR ("/tmp/" + fileName + "_rowPtr_TR.bin");

    vector<int> transformedRowPtr, transformedColIdx;
    vector<double> transformedValues;
    buildTransformedMatrix(transformedRowPtr, transformedColIdx, transformedValues);
    returnVal = transformedValues.size();
    
    // if loops exist
    if(analyzer->getSingleLoopRows()) {
      string pathB ("/tmp/" + fileName + "_b_TR.bin");
      string pathRowIndices ("/tmp/" + fileName + "_rowIndices_TR.bin");

      #ifdef PAR
      #pragma omp parallel sections num_threads(2)
      {
        #pragma omp section
        { dumpToMem(b, pathB.c_str(), LCSR->getRows()-1); }

        #pragma omp section
        { dumpToMem(rowIndices, pathRowIndices.c_str(), rowIndices.size()); }
      }
      #else
        dumpToMem(b, pathB.c_str(), LCSR->getRows()-1);
        dumpToMem(rowIndices, pathRowIndices.c_str(), rowIndices.size());
      #endif

    }

    #ifdef PAR
    #pragma omp parallel sections num_threads(3)
    {
      #pragma omp section
      { dumpToMem(transformedValues, pathVals_TR.c_str(), transformedValues.size());    }

      #pragma omp section
      { dumpToMem(transformedRowPtr, pathRowPtr_TR.c_str(), transformedRowPtr.size());  }

      #pragma omp section
      { dumpToMem(transformedColIdx, pathParents_TR.c_str(), transformedColIdx.size()); }
    }
    #else
      dumpToMem(transformedValues, pathVals_TR.c_str(), transformedValues.size());
      dumpToMem(transformedRowPtr, pathRowPtr_TR.c_str(), transformedRowPtr.size());
      dumpToMem(transformedColIdx, pathParents_TR.c_str(), transformedColIdx.size());
    #endif
  #endif

  #if !defined(REWRITE_ENABLED) || defined(DUMP_CSR) 
    string pathVals ("/tmp/" + fileName + "_vals.bin");
    string pathParents ("/tmp/" + fileName + "_parents.bin");
    string pathRowPtr ("/tmp/" + fileName + "_rowPtr.bin");

    // if loops exist
    if(analyzer->getSingleLoopRows()) {
      string pathB ("/tmp/" + fileName + "_b.bin");
      string pathRowIndices ("/tmp/" + fileName + "_rowIndices.bin");

      #ifdef PAR
      #pragma omp parallel sections num_threads(2)
      {
        #pragma omp section
        { dumpToMem(b, pathB.c_str(), LCSR->getRows()-1); }

        #pragma omp section
        { dumpToMem(rowIndices, pathRowIndices.c_str(), rowIndices.size()); }
      }
      #else
        dumpToMem(b, pathB.c_str(), LCSR->getRows()-1);
        dumpToMem(rowIndices, pathRowIndices.c_str(), rowIndices.size());
      #endif
    }

    #ifdef PAR
    #pragma omp parallel sections num_threads(3)
    {
      #pragma omp section
      { dumpToMem(LCSR->getVals(), pathVals.c_str(), LCSR->getNNZs()); }

      #pragma omp section
      { dumpToMem(LCSR->getRowPtr(), pathRowPtr.c_str(), LCSR->getRowPtr().size());  }

      #pragma omp section
      { dumpToMem(LCSR->getColIdx(), pathParents.c_str(), LCSR->getColIdx().size()); }
    }
    #else
      dumpToMem(LCSR->getVals(), pathVals.c_str(), LCSR->getNNZs());
      dumpToMem(LCSR->getRowPtr(), pathRowPtr.c_str(), LCSR->getRowPtr().size());
      dumpToMem(LCSR->getColIdx(), pathParents.c_str(), LCSR->getColIdx().size());
    #endif

    cout << "Matrix dimensions:\n";
    cout << "rowPtr: " << LCSR->getRowPtr().size() << "\n";
    cout << "colIdx: " << LCSR->getColIdx().size() << "\n";
    cout << "values: " << LCSR->getVals().size() << "\n";
  #endif

  auto t8 = std::chrono::steady_clock::now();
  #ifdef PAR
  cout << "* chrono dump data: " << chrono::duration<double>(t8 - t7).count()*1000 << "\n";
  #else
  cout << "chrono dump data: " << chrono::duration<double>(t8 - t7).count()*1000 << "\n";
  #endif


  return returnVal;
}

#ifdef REWRITE_ENABLED
void Rewrite::buildTransformedMatrix(vector<int>& transformedRowPtr, vector<int>& transformedColIdx, vector<double>& transformedValues){
  DAG& dag = analyzer->getDAG();
  int rows = matrixCSR->getNumOfRows();
  vector<vector<double>>& values = analyzer->getValues();

  transformedRowPtr.push_back(0);

  for(int i = 0; i < rows; i++) {
    vector<double>& rowValues = values[i];
    transformedValues.insert(transformedValues.end(), rowValues.begin(), rowValues.end());

    vector<int>& parents = dag[i].first;
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
