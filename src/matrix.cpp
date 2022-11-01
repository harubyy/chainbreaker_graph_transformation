#include "matrix.h"

uf_matrix_t* Matrix::getUF_matrix() {
  return uf_mat;
}

int Matrix::getNumOfRows() {
  return rows;
}

int Matrix::getNumOfCols() {
  return cols;
}

int Matrix::getNumOfVals() {
  return nnzs;
}

int* Matrix::getRowPtr() {
  return rowPtr;
}

int* Matrix::getColIdx() {
  return colIdx;
}

double* Matrix::getVals() {
  return vals;
}

int Matrix::getID() {
  return matID;
}

FORMAT Matrix::getFormat() {
  return format;
}

Part* Matrix::getL() {
  return L;
}

int Matrix::convertToCSR() {
  format = uf_matrix_coord_to_csr_int32(field, rows, cols, nnzs, &rowPtr, &colIdx, (void **) &vals) ? COO : CSR;
  return !format;
}

int Matrix::convertToCSC() {
  format = uf_matrix_coord_to_csc_int32(field, rows, cols, nnzs, &rowPtr, &colIdx, (void **) &vals) ? COO : CSC;
  return !format;
}

void Matrix::extractLCSC() {
  cout << "extract L for CSC format\n";

  if (L == nullptr) {
    L = new Part();

    // here actually we have rowIdx, colPtr and vals since this is CSC but we'll use the same class.
    // TODO: maybe we can generalize the pointer names so that it's not confusing for CSC
    vector<int>& rowPtrL = L->getRowPtr();
    vector<int>& colIdxL = L->getColIdx();
    vector<double>& valsL = L->getVals();

    rowPtrL.reserve(nnzs);
    colIdxL.reserve(cols + 1);
    valsL.reserve(nnzs);

    int nnzL = 0;
    int nnzPtr = 0;
    colIdxL.push_back(0);
    int count = 0;
    for (int i = 0; i < cols; i++) {
      count = 0;
      for (int j = colIdx[i]; j < colIdx[i + 1]; j++) {
  //        cout << "rowPtr[" << j << "]= " << rowPtr[j] << ", " << "i= " << i << ", nnzPtr= " << nnzPtr << "\n";
        if (rowPtr[j] >= i) {
          rowPtrL.push_back(rowPtr[j]);
          if((rowPtr[j] == i) && (rowPtr[j] == 0)) {
          //if(rowPtr[j] == i) {
            valsL.push_back(1.0); // vals[j];
          } else {
            valsL.push_back(vals[j]); // vals[j];
          }
          nnzPtr++;
          count++;
        }
      }
  
      colIdxL.insert(colIdxL.begin() + (i + 1), colIdxL[i] + count);
    }
  
    colIdxL.insert(colIdxL.begin() + cols, colIdxL[cols - 1] + count);
  
    printf("A's unit-lower triangular L: ( %i, %i ) nnz = %i\n", rows, cols, nnzPtr);
  
    L->setDimensions(nnzPtr, cols+1, nnzPtr);
  
    rowPtrL.shrink_to_fit();
    valsL.shrink_to_fit();
    // do this as well just in case for reserve allocating more than necessary at the beginning
    colIdxL.shrink_to_fit();
  
    format = CSC;
  } else {
    cout << "L already exists.\n";
    return;
  }
}

void Matrix::extractLCSR() {
  cout << "extract L for CSR format\n";

  if (L == nullptr) {
    L = new Part();

    // here actually we have rowIdx, colPtr and vals since this is CSC but we'll use the same class.
    // TODO: maybe we can generalize the pointer names so that it's not confusing for CSC
    vector<int>& rowPtrL = L->getRowPtr();
    vector<int>& colIdxL = L->getColIdx();
    vector<double>& valsL = L->getVals();

    rowPtrL.reserve(rows + 1);
    colIdxL.reserve(nnzs);
    valsL.reserve(nnzs);

    int nnzL = 0;
    int nnzPtr = 0;
    rowPtrL.push_back(0);
    int count = 0;
    for (int i = 0; i < rows; i++) {
      count = 0;
      for (int j = rowPtr[i]; j < rowPtr[i + 1]; j++) {
        if (colIdx[j] <= i) {
          colIdxL.push_back(colIdx[j]);
          if((colIdx[j] == i) && (colIdx[j] == 0)) {
          //if(colIdx[j] == i) {
            valsL.push_back(1.0); // vals[j];
          } else {
            valsL.push_back(vals[j]); // vals[j];
          }
          nnzPtr++;
          count++;
        }
      }
  
      rowPtrL.insert(rowPtrL.begin() + (i + 1), rowPtrL[i] + count);
    }
  
    rowPtrL.insert(rowPtrL.begin() + rows, rowPtrL[rows - 1] + count);
  
    printf("A's unit-lower triangular L: ( %i, %i ) nnz = %i\n", rows, cols, nnzPtr);
  
    L->setDimensions(rows + 1, nnzPtr, nnzPtr);
  
    rowPtrL.shrink_to_fit();
    valsL.shrink_to_fit();
    // do this as well just in case for reserve allocating more than necessary at the beginning
    colIdxL.shrink_to_fit();
  
    format = CSR;
  } else {
    cout << "L already exists.\n";
    return;
  }
}

