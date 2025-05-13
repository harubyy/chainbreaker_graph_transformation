// Compile with:
// nvcc -lcusparse -o gpu_solver gpu_solver.cu

#include <cuda_runtime.h>
#include <cusparse.h>
#include <iostream>
#include <vector>

#define CHECK_CUDA(err) if(err != cudaSuccess) { \
  std::cerr << "CUDA error: " << cudaGetErrorString(err) << "\n"; std::exit(EXIT_FAILURE); }
#define CHECK_CUSPARSE(err) if(err != CUSPARSE_STATUS_SUCCESS) { \
  std::cerr << "cuSPARSE error\n"; std::exit(EXIT_FAILURE); }

using namespace std;

void readCSRData(const char* path, const char* matrix_name, const char* transformed, int** rowPtrM, int** colIdxM, double** valsM, int* rows, int* cols, int* nnzs);

int main(int argc, char* argv[]) {
  if(argc != 2 && argc != 4) {
    cout << "Usage: programName matID\n";
    cout << "Usage: programName path_to_bin_files matrix_name true/false (transformed/not)\n";
    return 0;
  }

  cout << "injecting CSR data\n";
  int* rowPtr; int* colIdx; double* vals;
  int rows, cols, nnzs;
  readCSRData(argv[1], argv[2], argv[3], &rowPtr, &colIdx, &vals, &rows, &cols, &nnzs);


  // Lower triangular CSR matrix
/*  int h_csrRowPtr[] = {0, 1, 3, 4, 6, 7};
  int h_csrColInd[] = {0, 0, 1, 2, 0, 3, 4};
  float h_csrVal[]   = {3, 2, -2, 1, 4, 1, 3};
  int rows = 5, nnzs = 7;*/

  //float h_xRef[] = {1, 1, 1, 1, 1};
//  float h_xRef[rows];
  float* h_xRef = (float*)malloc(rows * sizeof(float));
  memset(h_xRef, 1.0, rows * sizeof(float));
  float h_b[rows];

  // Manually compute b = L * xRef
  for(int i = 0; i < rows; ++i) {
    float sum = 0;
    for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j)
      sum += vals[j] * h_xRef[colIdx[j]];
    h_b[i] = sum;
  }

  // Device memory
  int *d_csrRowPtr, *d_csrColInd;
  float *d_csrVal, *d_b, *d_x;
  CHECK_CUDA(cudaMalloc((void**)&d_csrRowPtr, (rows+1)*sizeof(int)));
  CHECK_CUDA(cudaMalloc((void**)&d_csrColInd, nnzs*sizeof(int)));
  CHECK_CUDA(cudaMalloc((void**)&d_csrVal,    nnzs*sizeof(float)));
  CHECK_CUDA(cudaMalloc((void**)&d_b,         rows*sizeof(float)));
  CHECK_CUDA(cudaMalloc((void**)&d_x,         rows*sizeof(float)));

  CHECK_CUDA(cudaMemcpy(d_csrRowPtr, rowPtr, (rows+1)*sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA(cudaMemcpy(d_csrColInd, colIdx, nnzs*sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA(cudaMemcpy(d_csrVal,    vals,    nnzs*sizeof(float), cudaMemcpyHostToDevice));
  CHECK_CUDA(cudaMemcpy(d_b,         h_b,         rows*sizeof(float), cudaMemcpyHostToDevice));

  // cuSPARSE handle & descriptor
  cusparseHandle_t handle;
  cusparseSpMatDescr_t matA;
  cusparseDnVecDescr_t vecX, vecB;
  void* dBuffer = nullptr;
  size_t bufferSize = 0;

  CHECK_CUSPARSE(cusparseCreate(&handle));

  CHECK_CUSPARSE(cusparseCreateCsr(&matA, rows, rows, nnzs,
    d_csrRowPtr, d_csrColInd, d_csrVal,
    CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
    CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F));

  CHECK_CUSPARSE(cusparseCreateDnVec(&vecX, rows, d_x, CUDA_R_32F));
  CHECK_CUSPARSE(cusparseCreateDnVec(&vecB, rows, d_b, CUDA_R_32F));

  // Create SpSV info
  cusparseSpSVDescr_t spsvDescr;
  CHECK_CUSPARSE(cusparseSpSV_createDescr(&spsvDescr));

  float alpha = 1.0f;
  CHECK_CUSPARSE(cusparseSpSV_bufferSize(
    handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
    &alpha, matA, vecB, vecX, CUDA_R_32F,
    CUSPARSE_SPSV_ALG_DEFAULT, spsvDescr, &bufferSize));

  CHECK_CUDA(cudaMalloc(&dBuffer, bufferSize));

  CHECK_CUSPARSE(cusparseSpSV_analysis(
    handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
    &alpha, matA, vecB, vecX, CUDA_R_32F,
    CUSPARSE_SPSV_ALG_DEFAULT, spsvDescr, dBuffer));

  CHECK_CUSPARSE(cusparseSpSV_solve(
    handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
    &alpha, matA, vecB, vecX, CUDA_R_32F,
    CUSPARSE_SPSV_ALG_DEFAULT, spsvDescr));

  float h_x[rows];
  CHECK_CUDA(cudaMemcpy(h_x, d_x, rows*sizeof(float), cudaMemcpyDeviceToHost));

  std::cout << "Computed solution x:\n";
  for(int i = 0; i < 40; ++i)
    std::cout << h_x[i] << " ";
  std::cout << "\nExpected (xRef):\n";
  for(int i = 0; i < 40; ++i)
    std::cout << h_xRef[i] << " ";
  std::cout << "\n";

  cudaFree(d_csrRowPtr); cudaFree(d_csrColInd); cudaFree(d_csrVal);
  cudaFree(d_b); cudaFree(d_x); cudaFree(dBuffer);
  cusparseDestroySpMat(matA);
  cusparseDestroyDnVec(vecX); cusparseDestroyDnVec(vecB);
  cusparseSpSV_destroyDescr(spsvDescr);
  cusparseDestroy(handle);
  return 0;
}

