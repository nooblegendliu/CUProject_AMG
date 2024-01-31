#ifndef cuAMG_CORE_MATRIX_CUH_
#define cuAMG_CORE_MATRIX_CUH_

#include <cuda.h>
#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>
#include <thrust/device_vector.h>

namespace gpu {
	class csrMatrix {
	public:
		int n_rows;	// number of rows
		int n_cols;	// number of columns
		int nnz;	// number of non-zero elements

		// CUSPARSE library handle and statues;
		cusparseHandle_t cusparse_handle;
		cusparseStatus_t cusparse_status;

		// CUBLAS library handle and statues;
		cublasHandle_t cublas_handle;
		cublasStatus_t cublas_status;

		// cusparse matrix descriptor
		cusparseMatDescr_t descrA;
		cusparseOperation_t transA;

		// vector of csr matrix
		thrust::device_vector<int> row_ptr;
		thrust::device_vector<int> col_ind;
		thrust::device_vector<double> values;

	public:
		// constructor
		explicit csrMatrix() = delete;


	};
}









#endif