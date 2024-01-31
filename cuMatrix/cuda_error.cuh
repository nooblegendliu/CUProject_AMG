#ifndef cuAMG_UTILS_ERROR_CUH
#define cuAMG_UTILS_ERROR_CUH

#include <stdexcept>
#include <cuda_runtime.h>

inline void cudaCheckError(cudaError_t error) {
	if (error != cudaSuccess) {
		std::string msg{"CUDA error "};
		msg += cudaGetErrorName(error);
		msg += ": ";
		msg += cudaGetErrorString(error);
		throw std::runtime_error{msg};
	}
}

template <typename F>
void cudaChecked(F func) {
	func();
	cudaDeviceSynchronize();
	cudaCheckError(cudaGetLastError());
}

#endif