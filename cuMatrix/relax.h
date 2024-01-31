#ifndef AMG_UTILS_RELAX_HPP_
#define AMG_UTILS_RELAX_HPP_

#include <cfloat>

#include "vector.h"
#include "matrix.h"
#include "level.h"
#include <Eigen/Sparse>

namespace AMG {
	void jacobi(CSRMatrix* A, Vector& b, Vector& x, Vector& tmp,
		int num_sweeps = 2, double omega = 1.0, bool Print = false);
	void gaussSeidel(CSRMatrix* A, Vector& b, Vector& x,
		int num_sweeps = 2, double omega = 1.0, bool Print = false);


	void sor(CSRMatrix* A, Vector& b, Vector& x, Vector& tmp,
		int num_sweeps = 2, double omega = 1.0, bool Print = false);
	void ssor(CSRMatrix* A, Vector& b, Vector& x, Vector& tmp,
		int num_sweeps = 2, double omega = 1.0, bool Print = false);

	void Direct_solve(CSRMatrix* A, Vector& b, Vector& x);

}	// namespace AMG





#endif

