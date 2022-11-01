#pragma once
#include <iostream>
#include <stdint.h>
#include <libufget.h>
#include <vector>

using namespace std;

#define CHK_MEM_ERR(succ,event) \
  do { \
   if(succ != 0) \
     printf("posix_memalign is not successfull!\n"); \
     \
   if(succ == ENOMEM) \
     printf("%s: insufficient memory\n", event); \
   else if(succ == EINVAL) \
     printf("%s: alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n", event); \
  } while (0)


typedef enum {
  COO = 0,
  CSR = 1,
  CSC = 2
} FORMAT;

class Part {
  private:
    vector<int> rowPtr;
    vector<int> colIdx;
    vector<double> vals;

    // TODO check for type of nnzs
    int rows, cols, nnzs;

  public:
    vector<int>& getRowPtr() {
      vector<int>& ref = rowPtr;

      return ref;
    }

    vector<int>& getColIdx() {
      vector<int>& ref = colIdx;

      return ref;
    }

    vector<double>& getVals() {
      vector<double>& ref = vals;

      return ref;
    }

    int getRows() {
      return rows;
    }

    int getCols() {
      return cols;
    }

    int getNNZs() {
      return nnzs;
    }

    void setDimensions(int rows, int cols, int nnzs) {
      this->rows = rows;
      this->cols = cols;
      this->nnzs = nnzs;
    }

    void print() {
      cout << "rows: " << rows << " cols: " << cols << " nnzs: " << nnzs << "\n";
      cout << "sizes of vectors: " << rowPtr.size() << " " << colIdx.size() << " " << vals.size() << "\n";
      printf("RowPtr:\n");
      for (int i = 0; i < rows; i++)
        cout << rowPtr[i] << ", ";
      cout << "\n";
  
      printf("Cols:\n");
      for (int i = 0; i < cols; i++)
        cout << colIdx[i] << ", ";
      cout << "\n";
  
      printf("vals:\n");
      for (int i = 0; i < nnzs; i++)
        cout << vals[i] << ", ";
      cout << "\n";
  
    }
};

// TODO: check for alignment
class Matrix {
  int rows, cols, nnzs;
  int* rowPtr;
  int* colIdx;
  double* vals;

  uf_matrix_t* uf_mat;
  uf_field_t field;
  int matID;

  Part* L;
  Part* U;
  FORMAT format;

 public:
  Matrix(uf_collection_t *collection, int matID) {
    uf_mat = uf_collection_get_by_id(collection, matID);
    printf("Analyzing Matrix[%4d] %s/%s\n", uf_mat->id, uf_mat->group_name, uf_mat->name);

    uf_matrix_coord_int32(&field, &rows, &cols, &nnzs, &rowPtr, &colIdx, (void **) &vals, uf_mat);
    cout << rows << ", " << cols << ", " << nnzs << endl;

    L = U = nullptr;
    format = COO;
  }

  ~Matrix() {
    free(rowPtr);
    free(colIdx);
    free(vals);
    delete uf_mat;

    delete L;
    delete U;
  }

  uf_matrix_t* getUF_matrix();
  int getNumOfRows();
  int getNumOfCols();
  int getNumOfVals();

  int* getRowPtr();
  int* getColIdx();
  double* getVals();

  int getID();
  FORMAT getFormat();
  Part* getL();

  int convertToCSR();
  int convertToOptimizedCSR();
  int convertToCSC();
  void extractLCSR();
  void extractLCSC();

  // TODO:
  void print() {
    printf("Matrix ID:%4d\nMatrix Name:%s\n", uf_mat->id, uf_mat->name);

    int rowCount, colCount, valCount;

    switch (format) {
      case COO:rowCount = colCount = valCount = nnzs;
        break;
      case CSR:rowCount = rows + 1;
        colCount = valCount = nnzs;
        break;
      case CSC:colCount = cols + 1;
        rowCount = valCount = nnzs;
        break;
      default:rowCount = colCount = valCount = 0;
        break;
    }

    if (L != NULL) {
      cout << "L exists for " << uf_mat->name << "\n";
      L->print();
    }

    if (U != NULL) {
      cout << "U exists for " << uf_mat->name << "\n";
      U->print();
    }
  }
};
