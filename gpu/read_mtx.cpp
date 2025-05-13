#include <iostream>
#include <string.h>
#include <stdint.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;

void check_err_open_lseek(int result) {
  if(result == -1) {
    perror("Error opening file/lseek for reading");
    exit(EXIT_FAILURE);
  }
}

void check_fail_mmap(double* map, int fd) {
  if(map == MAP_FAILED) {
    close(fd);
    perror("Error mmapping the file");
    exit(EXIT_FAILURE);
  }
}

void readCSRData(const char* path, const char* matrix_name, const char* transformed, int** rowPtrM, int** colIdxM, double** valsM, int* rows, int* cols, int* nnzs) {
  std::string path_prefix(path);
  path_prefix += matrix_name;
  std::string extension(".bin");

  std::string f_row_ptr(path_prefix);
  f_row_ptr += "_rowPtr";
  if(!std::string(transformed).compare("true"))
    f_row_ptr += "_TR";
  f_row_ptr += extension;

  std::string f_col_idx(path_prefix);
  f_col_idx += "_parents";
  if(!std::string(transformed).compare("true"))
    f_col_idx += "_TR";
  f_col_idx += extension;

  std::string f_values(path_prefix);
  f_values += "_vals";
  if(!std::string(transformed).compare("true"))
    f_values += "_TR";
  f_values += extension;

  printf("reading from: %s, %s, %s\n", f_row_ptr.c_str(), f_col_idx.c_str(), f_values.c_str());

  int fd_parents, fd_vals, fd_rowPtr;

  fd_parents = open(f_col_idx.c_str(), O_RDWR);
  fd_rowPtr = open(f_row_ptr.c_str(), O_RDWR);
  fd_vals = open(f_values.c_str(), O_RDWR);

  check_err_open_lseek(fd_parents);
  check_err_open_lseek(fd_rowPtr);
  check_err_open_lseek(fd_vals);

  struct stat buf;
  fstat(fd_parents, &buf);
  off_t parents_size = buf.st_size;

  std::cout << "parents_size: " << parents_size <<  " became(nnzs): " << parents_size/4 <<  "\n";

  fstat(fd_rowPtr, &buf);
  off_t x_size = buf.st_size;
  *rows = *cols = x_size/4-1;

  std::cout << "x_size: " << x_size <<  " became: " << x_size/4 <<  "\n";
  std::cout << "rows, cols: " << *rows << "\n";

  fstat(fd_vals, &buf);
  off_t x_size_ = buf.st_size;
  *nnzs =  x_size_/8;

  std:cout << "nnzs: " << *nnzs << "\n";

//  int* rowPtrM; int* colIdxM; double* valsM;
 * colIdxM = static_cast<int*>(mmap(0, *nnzs * sizeof(int), PROT_READ, MAP_SHARED, fd_parents, 0));
  *rowPtrM = static_cast<int*>(mmap(0, ((*rows) + 1) * sizeof(int), PROT_READ, MAP_SHARED, fd_rowPtr, 0));
  *valsM = static_cast<double*>(mmap(0, *nnzs * sizeof(double), PROT_READ, MAP_SHARED, fd_vals, 0));

  check_fail_mmap((double*)rowPtrM, fd_rowPtr);
  check_fail_mmap((double*)colIdxM, fd_parents);
  check_fail_mmap((double*)valsM, fd_vals);

/*  L->allocateForCSR(rows, nnzs);

  L->setNumOfCols(cols);
  L->setNumOfRows(rows);
  L->setNumOfVals(nnzs);
  L->setDataCSR(rowPtrM, colIdxM, valsM);*/
}

