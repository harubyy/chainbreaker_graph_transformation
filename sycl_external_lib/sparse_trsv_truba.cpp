/*******************************************************************************
* Copyright 2019-2021 Intel Corporation.
*
* This software and the related documents are Intel copyrighted  materials,  and
* your use of  them is  governed by the  express license  under which  they were
* provided to you (License).  Unless the License provides otherwise, you may not
* use, modify, copy, publish, distribute,  disclose or transmit this software or
* the related documents without Intel's prior written permission.
*
* This software and the related documents  are provided as  is,  with no express
* or implied  warranties,  other  than those  that are  expressly stated  in the
* License.
*******************************************************************************/

/*
*
*  Content:
*       This example demonstrates use of oneAPI Math Kernel Library (oneMKL)
*       DPCPP API oneapi::mkl::sparse::trsv to perform
*
*       sparse triangular solve on a SYCL device (Host, CPU, GPU).
*
*       op(A) * y = x
*
*       where op() is defined by one of
*oneapi::mkl::transpose::{nontrans,trans,conjtrans}
*
*
*       The supported floating point data types for gemm matrix data are:
*           float
*           double
*
*
*******************************************************************************/


/*
 * namespace oneapi::mkl:sparse {
    void set_csr_data (
        oneapi::mkl::sparse::matrix_handle_t handle,
        const intType num_rows,
        const intType num_cols,
        oneapi::mkl::index_base index,
        cl::sycl::buffer<intType, 1> &row_ptr,
        cl::sycl::buffer<intType, 1> &col_ind,
        cl::sycl::buffer<fp, 1> &val)
 }
 */

// stl includes
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <vector>
#include <string>

#include "mkl.h"
#include "oneapi/mkl.hpp"
#include <CL/sycl.hpp>

/*
 * COMPILE:
 * dpcpp -I${MKLROOT}/include -I${MKLROOT}/examples/dpcpp  -DMKL_ILP64 -DSYCL_DEVICES_host  -c  sparse_trsv.cpp -o sptrsv.o
 *
 * LINK:
 * dpcpp sptrsv.o -fsycl-device-code-split=per_kernel \
"${MKLROOT}/lib/intel64"/libmkl_sycl.a -Wl,-export-dynamic -Wl,--start-group \
"${MKLROOT}/lib/intel64"/libmkl_intel_ilp64.a \
"${MKLROOT}/lib/intel64"/libmkl_tbb_thread.a   \
"${MKLROOT}/lib/intel64"/libmkl_core.a -Wl,--end-group -L${TBBROOT}/lib/intel64/gcc4.8 -ltbb -lsycl -lOpenCL \
-lpthread -lm -ldl -o sptrsv
 */

// local includes
//#include "../common/common_for_examples.hpp"
#include "common_for_examples.hpp"
#include "common_for_sparse_examples.hpp"
//#include "common_for_sparse_examples.hpp"

//
// Main example for Sparse Triangular Solver consisting of
// initialization of A matrix, x and y vectors.
// Then the following system is solved
//
// op(A) * y = x
//
// and finally the results are post processed.
//

#include "util.hpp"
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>

int* parents;
int* rowPtr;
double* values;

void read_matrix(int x_size, int parents_size, const char* matrix_name, const char* transformed) {
  std::string path_prefix("/truba/home/buyilmaz/chainbreaker/Results/original/");
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

  int fd_parents;
  int fd_vals;
  int fd_rowPtr;


  //std::cout << f_row_ptr << " " << f_col_idx << " " << f_values << "\n";
  fd_parents = open(f_col_idx.c_str(), O_RDWR);
  fd_rowPtr = open(f_row_ptr.c_str(), O_RDWR);
  fd_vals = open(f_values.c_str(), O_RDWR);

  check_err_open_lseek(fd_parents);
  check_err_open_lseek(fd_rowPtr);
  check_err_open_lseek(fd_vals);

  parents = static_cast<int*>(mmap(0, parents_size * sizeof(int), PROT_READ, MAP_SHARED, fd_parents, 0));
  values = static_cast<double*>(mmap(0, parents_size * sizeof(double), PROT_READ, MAP_SHARED, fd_vals, 0));
  rowPtr = static_cast<int*>(mmap(0, x_size * sizeof(int), PROT_READ, MAP_SHARED, fd_rowPtr, 0));

  check_fail_mmap((double*)rowPtr, fd_rowPtr);
  check_fail_mmap((double*)parents, fd_parents);
  check_fail_mmap(values, fd_vals);

/*  std::cout << "parent_size: " << parents_size << " x_size: " << x_size << "\n";
  std::cout << "values:\n";
  for(int i = 0; i < parents_size; i++)
    std::cout << values[i] << ", ";
  std::cout << "\n";*/

//  std::cout << "rowPtr[x_size-1]: " << rowPtr[x_size-1] << "\n\n";
}

template <typename fp, typename intType>
//int run_sparse_triangular_solve_example(const cl::sycl::device &dev)
int run_sparse_triangular_solve_example(const cl::sycl::device &dev, int x_size, int parents_size)
{
    // Initialize data for Sparse Triangular Solve
    oneapi::mkl::transpose transA = oneapi::mkl::transpose::nontrans;

    // Matrix data size
    intType size  = 4;
    intType nrows = size * size * size;

    // Input matrix in CSR format
    std::vector<intType, mkl_allocator<intType, 64>> ia;
    std::vector<intType, mkl_allocator<intType, 64>> ja;
    std::vector<fp, mkl_allocator<fp, 64>> a;

  //std::cout << "parent_size: " << parents_size << " x_size: " << x_size << "\n";
    /*ia.resize(nrows + 1);   // rowPtr 64
    ja.resize(27 * nrows);  // colIdx 64 * 27
    a.resize(27 * nrows);   // values 64 * 27*/
    ia.resize(x_size+1);   // rowPtr
    ja.resize(parents_size);  // colIdx
    a.resize(parents_size);   // values

//    generate_sparse_matrix<fp, intType>(size, ia, ja, a);

    ia.insert(ia.begin(), rowPtr, rowPtr + (x_size+1));
    ja.insert(ja.begin(), parents, parents + parents_size);
    a.insert(a.begin(), values, values + parents_size);

    // Vectors x and y
    std::vector<fp, mkl_allocator<fp, 64>> x;
    std::vector<fp, mkl_allocator<fp, 64>> y;
    std::vector<fp, mkl_allocator<fp, 64>> z;
    x.resize(x_size);
    y.resize(x_size);
    z.resize(x_size);

    // Init vectors x and y
    //for (int i = 0; i < nrows; i++) {
    for (int i = 0; i < x_size; i++) {
        x[i] = set_fp_value(fp(1.0), fp(0.0));
        y[i] = set_fp_value(fp(0.0), fp(0.0));
        z[i] = set_fp_value(fp(0.0), fp(0.0));
    }

    // Catch asynchronous exceptions
    auto exception_handler = [](cl::sycl::exception_list exceptions) {
        for (std::exception_ptr const &e : exceptions) {
            try {
                std::rethrow_exception(e);
            }
            catch (cl::sycl::exception const &e) {
                std::cout << "Caught asynchronous SYCL "
                             "exception during sparse::trsv:\n"
                          << e.what() << std::endl;
            }
        }
    };

    //
    // Execute Triangulat Solve
    //

    // create execution queue and buffers of matrix data
    cl::sycl::queue main_queue(dev, exception_handler);

    /*cl::sycl::buffer<intType, 1> ia_buffer(ia.data(), (nrows + 1));
    cl::sycl::buffer<intType, 1> ja_buffer(ja.data(), (ia[nrows]));
    cl::sycl::buffer<fp, 1> a_buffer(a.data(), (ia[nrows]));
    cl::sycl::buffer<fp, 1> x_buffer(x.data(), x.size());
    cl::sycl::buffer<fp, 1> y_buffer(y.data(), y.size());*/
    cl::sycl::buffer<intType, 1> ia_buffer(ia.data(), x_size+1);
    cl::sycl::buffer<intType, 1> ja_buffer(ja.data(), (ia[x_size]));
    cl::sycl::buffer<fp, 1> a_buffer(a.data(), (ia[x_size]));
    cl::sycl::buffer<fp, 1> x_buffer(x.data(), x.size());
    cl::sycl::buffer<fp, 1> y_buffer(y.data(), y.size());

    // create and initialize handle for a Sparse Matrix in CSR format
    oneapi::mkl::sparse::matrix_handle_t handle;

    try {
        oneapi::mkl::sparse::init_matrix_handle(&handle);

        oneapi::mkl::sparse::set_csr_data(handle, x_size, x_size, oneapi::mkl::index_base::zero,
                                          ia_buffer, ja_buffer, a_buffer);
        //oneapi::mkl::sparse::set_csr_data(handle, nrows, nrows, oneapi::mkl::index_base::zero,
        //                                  ia_buffer, ja_buffer, a_buffer);

        // measure the operation
        auto start = std::chrono::steady_clock::now();


        oneapi::mkl::sparse::optimize_trsv(main_queue, oneapi::mkl::uplo::lower,
                                           oneapi::mkl::transpose::nontrans,
                                           oneapi::mkl::diag::unit, handle);

        // add oneapi::mkl::sparse::trsv to execution queue
        oneapi::mkl::sparse::trsv(main_queue, oneapi::mkl::uplo::lower,
                                  oneapi::mkl::transpose::nontrans, oneapi::mkl::diag::unit,
                                  handle, x_buffer, y_buffer);

        // measure the operation
        auto end = std::chrono::steady_clock::now();
        std::cout << "chainbreaker: " << std::chrono::duration_cast<std::chrono::duration<double>>(end-start).count()*1000 << " miliseconds\n";

        oneapi::mkl::sparse::release_matrix_handle(&handle);
    }
    catch (cl::sycl::exception const &e) {
        std::cout << "\t\tCaught synchronous SYCL exception:\n" << e.what() << std::endl;
        oneapi::mkl::sparse::release_matrix_handle(&handle);
        return 1;
    }
    catch (std::exception const &e) {
        std::cout << "\t\tCaught std exception:\n" << e.what() << std::endl;
        oneapi::mkl::sparse::release_matrix_handle(&handle);
        return 1;
    }

    //
    // Post Processing
    //

    std::cout << "\n\t\tsparse::trsv parameters:\n";
    std::cout << "\t\t\ttransA = "
              << (transA == oneapi::mkl::transpose::nontrans ?
                          "nontrans" :
                          (transA == oneapi::mkl::transpose::trans ? "trans" : "conjtrans"))
              << std::endl;
    //std::cout << "\t\t\tnrows = " << nrows << std::endl;
    std::cout << "\t\t\tnrows = " << x_size << std::endl;

    auto res = y_buffer.template get_access<cl::sycl::access::mode::read>();
    //for (intType row = 0; row < nrows; row++) {
    for (intType row = 0; row < x_size; row++) {
        fp tmp      = x[row];
        fp diag_val = set_fp_value(fp(0.0), fp(0.0));
        for (intType i = ia[row]; i < ia[row + 1]; i++) {
            if (ja[i] < row) {
//                std::cout << tmp << ", " << a[i] << ", " << z[ja[i]] << "\n";
                tmp -= a[i] * z[ja[i]];
            }
            else if (ja[i] == row) {
                //diag_val = a[i];
                diag_val = 1.0;
            }
        }
        z[row] = tmp / diag_val;
//        std::cout << "z:[row]: " << z[row] << "\n";
//        std::cout << res[row] << ", " << z[row] << ", " << row << "\n";
        if(!check_result(res[row], z[row], row)) return 1;
    }

    fp diff = set_fp_value(fp(0.0), fp(0.0));
    //for (intType i = 0; i < nrows; i++) {
    for (intType i = 0; i < x_size; i++) {
        diff += (z[i] - res[i]) * (z[i] - res[i]);
    }

    std::cout << "\n\t\t sparse::trsv residual:\n"
              << "\t\t\t" << diff << "\n\tFinished" << std::endl;

    return 0;
}

//
// Description of example setup, apis used and supported floating point type
// precisions
//
void print_example_banner()
{

    std::cout << "" << std::endl;
    std::cout << "###############################################################"
                 "#########"
              << std::endl;
    std::cout << "# Sparse Triangular Solve Example: " << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# y = op(A)^(-1) * x" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# where A is a sparse matrix in CSR format, x and y are "
                 "dense vectors"
              << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Using apis:" << std::endl;
    std::cout << "#   sparse::trsv" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Supported floating point type precisions:" << std::endl;
    std::cout << "#   float" << std::endl;
    std::cout << "#   double" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "###############################################################"
                 "#########"
              << std::endl;
    std::cout << std::endl;
}

//
// Main entry point for example.
//
// Dispatches to appropriate device types as set at build time with flag:
// -DSYCL_DEVICES_host -- only runs host implementation
// -DSYCL_DEVICES_cpu -- only runs SYCL CPU implementation
// -DSYCL_DEVICES_gpu -- only runs SYCL GPU implementation
// -DSYCL_DEVICES_all (default) -- runs on all: host, cpu and gpu devices
//
//  For each device selected and each supported data type, MatrixMultiplyExample
//  runs is with all supported data types
//

int main(int argc, char **argv)
{

  // x_size, parents_size, matrix_name, transformed (bool)
  read_matrix(atoi(argv[1]), atoi(argv[2]), argv[3], argv[4]);


    //print_example_banner();

    std::list<my_sycl_device_types> list_of_devices;
    set_list_of_devices(list_of_devices);

    int status = 0;
    for (auto it = list_of_devices.begin(); it != list_of_devices.end(); ++it) {

        cl::sycl::device my_dev;
        bool my_dev_is_found = false;
        get_sycl_device(my_dev, my_dev_is_found, *it);

        if (my_dev_is_found) {
            std::cout << "Running tests on " << sycl_device_names[*it] << ".\n";

/*            std::cout << "\tRunning with single precision real data type:" << std::endl;
              std::cout << argv[3] << ", ";
            //status = run_sparse_triangular_solve_example<float, std::int32_t>(my_dev);
            status = run_sparse_triangular_solve_example<float, std::int32_t>(my_dev, atoi(argv[1]), atoi(argv[2]));
            if(status != 0) return status;  */

            if (my_dev.get_info<cl::sycl::info::device::double_fp_config>().size() != 0) {
                std::cout << "\tRunning with double precision real data type:" << std::endl;
                std::cout << argv[3] << ", ";
                //status = run_sparse_triangular_solve_example<double, std::int32_t>(my_dev);
                status = run_sparse_triangular_solve_example<double, std::int32_t>(my_dev, atoi(argv[1]), atoi(argv[2]));
                if(status != 0) return status;
            }
        }
        else {
#ifdef FAIL_ON_MISSING_DEVICES
            std::cout << "No " << sycl_device_names[*it]
                      << " devices found; Fail on missing devices "
                         "is enabled.\n";
            return 1;
#else
            std::cout << "No " << sycl_device_names[*it] << " devices found; skipping "
                      << sycl_device_names[*it] << " tests.\n";
#endif
        }
    }

    return status;
}
