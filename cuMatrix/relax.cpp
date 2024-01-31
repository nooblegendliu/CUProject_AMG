#include "relax.h"
#include <iostream>
#include <iomanip>
namespace AMG {
	void jacobi(CSRMatrix* A, Vector& b, Vector& x, Vector& tmp,
		int num_sweeps, double omega, bool Print)
	{
		int row_start, row_end;
		double diag, row_sum;

		for (int iter = 0; iter < num_sweeps; iter++) {
			for (int i = 0; i < A->n_rows; i++) {
				tmp[i] = x[i];
			}

			for (int i = 0; i < A->n_rows; i++) {
				row_start = A->idx1[i];
				row_end = A->idx1[i + 1];
				if (row_start == row_end) continue;
				row_sum = 0.0;
				diag = 0.0;
				for (int j = row_start; j < row_end; j++) {
					int col = A->idx2[j];
					if (i == col)
						diag = A->vals[j];
					else
						row_sum += A->vals[j] * tmp[col];
				}
				if (fabs(diag) > zero_tol)
					x[i] = ((1.0 - omega) * tmp[i]) + (omega * (b[i] - row_sum) / diag);
			}
		}
		if (Print) {
			double res = 0.0;
			for (int i = 0; i < A->n_rows; i++) {
				row_start = A->idx1[i];
				row_end = A->idx1[i + 1];
				if (row_start == row_end) continue;
				row_sum = 0.0;
				for (int j = row_start; j < row_end; j++) {
					int col = A->idx2[j];
					row_sum += A->vals[j] * x[col];
				}
				res += (b[i] - row_sum) * (b[i] - row_sum);
			}
			std::cout << "Gauss-Seidel Residual: " << std::setprecision(12) << sqrt(res) << std::endl;
		}
	}	// jacobi
	void gaussSeidel(CSRMatrix* A, Vector& b, Vector& x,
		int num_sweeps, double omega, bool Print)
	{
		int row_start, row_end;
		double diag, row_sum;

		for (int iter = 0; iter < num_sweeps; iter++)
		{
			for (int i = 0; i < A->n_rows; i++) {
				row_start = A->idx1[i];
				row_end = A->idx1[i + 1];
				if (row_start == row_end) continue;
				row_sum = 0.0;
				diag = 0.0;
				for (int j = row_start; j < row_end; j++) {
					int col = A->idx2[j];
					if (i == col)
						diag = A->vals[j];
					else
						row_sum += A->vals[j] * x[col];
				}
				if (fabs(diag) > zero_tol)
					x[i] = ((1.0 - omega) * x[i]) + (omega * (b[i] - row_sum) / diag);
			}
		}
		if (Print) {
			double res = 0.0;
			for (int i = 0; i < A->n_rows; i++) {
				row_start = A->idx1[i];
				row_end = A->idx1[i + 1];
				if (row_start == row_end) continue;
				row_sum = 0.0;
				for (int j = row_start; j < row_end; j++) {
					int col = A->idx2[j];
					row_sum += A->vals[j] * x[col];
				}
				res += (b[i] - row_sum) * (b[i] - row_sum);
			}
			std::cout << "Gauss-Seidel Residual: " << std::setprecision(12) << sqrt(res) << std::endl;
		}
	}	// GaussSeidel

	void sor(CSRMatrix* A, Vector& b, Vector& x, Vector& tmp,
		int num_sweeps, double omega, bool Print)
	{
		int row_start, row_end, col;
		double diag_inv;
		double orig_x = 0.0;
		for (int iter = 0; iter < num_sweeps; iter++)
		{
			for (int i = 0; i < A->n_rows; i++)
			{
				orig_x = x[i];
				x[i] = b[i];
				row_start = A->idx1[i];
				row_end = A->idx1[i + 1];
				if (row_start == row_end) continue;
				diag_inv = 0.0;
				for (int j = row_start; j < row_end; j++) {
					col = A->idx2[j];
					if (i == col)
						diag_inv = omega / A->vals[j];
					else
						x[i] -= A->vals[j] * x[col];
				}
				x[i] = (1 - omega) * orig_x + diag_inv * x[i];
			}	// loop : rows
		}	// loop: iter
		if (Print) {
			double res = 0.0;
			double row_sum;
			for (int i = 0; i < A->n_rows; i++) {
				row_start = A->idx1[i];
				row_end = A->idx1[i + 1];
				if (row_start == row_end) continue;
				row_sum = 0.0;
				for (int j = row_start; j < row_end; j++) {
					int col = A->idx2[j];
					row_sum += A->vals[j] * x[col];
				}
				res += (b[i] - row_sum) * (b[i] - row_sum);
			}
			std::cout << "Gauss-Seidel Residual: " << std::setprecision(12) << sqrt(res) << std::endl;
		}
	}// sor
	void ssor(CSRMatrix* A, Vector& b, Vector& x, Vector& tmp, int num_sweeps,
		double omega, bool Print)
	{
		int row_start, row_end;
		double diag_inv;
		double orig_x = 0;

		for (int iter = 0; iter < num_sweeps; iter++)
		{
			for (int i = 0; i < A->n_rows; i++)
			{
				orig_x = x[i];
				x[i] = b[i];
				row_start = A->idx1[i];
				row_end = A->idx1[i + 1];
				if (row_start == row_end) continue;

				diag_inv = omega / A->vals[row_start];
				for (int j = row_start + 1; j < row_end; j++)
				{
					x[i] -= A->vals[j] * x[A->idx2[j]];
				}
				x[i] = diag_inv * x[i] + (1 - omega) * orig_x;
			}

			for (int i = A->n_rows - 1; i >= 0; i--)
			{
				orig_x = x[i];
				x[i] = b[i];
				row_start = A->idx1[i];
				row_end = A->idx1[i + 1];
				if (row_start == row_end) continue;

				diag_inv = omega / A->vals[row_start];
				for (int j = row_start + 1; j < row_end; j++)
				{
					x[i] -= A->vals[j] * x[A->idx2[j]];
				}
				x[i] = diag_inv * x[i] + (1 - omega) * orig_x;
			}
		}
		double row_sum;
		if (Print) {
			double res = 0.0;
			for (int i = 0; i < A->n_rows; i++) {
				row_start = A->idx1[i];
				row_end = A->idx1[i + 1];
				if (row_start == row_end) continue;
				row_sum = 0.0;
				for (int j = row_start; j < row_end; j++) {
					int col = A->idx2[j];
					row_sum += A->vals[j] * x[col];
				}
				res += (b[i] - row_sum) * (b[i] - row_sum);
			}
			std::cout << "Gauss-Seidel Residual: " << std::setprecision(12) << sqrt(res) << std::endl;
		}
	}	// ssor

	void Direct_solve(CSRMatrix* A, Vector& b, Vector& x)
	{
		Eigen::SparseMatrix<double> Eigen_mat(A->n_rows, A->n_cols);
		std::vector<Eigen::Triplet<double>> triplets;
		triplets.reserve(A->nnz);
		for (int i = 0; i < A->n_rows; ++i)
		{
			for (int i = 0; i < A->n_rows; ++i) {
				for (int j = A->idx1[i]; j < A->idx1[i + 1]; ++j) {
					int col = A->idx2[j];
					double val = A->vals[j];
					triplets.push_back(Eigen::Triplet<double>(i, col, val));
				}
			}
		}

		Eigen_mat.setFromTriplets(triplets.begin(), triplets.end());
		Eigen::VectorXd Eigen_b(b.size());
		for (int i = 0; i < b.size(); ++i)
		{
			Eigen_b(i) = b[i];
		}

		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
		solver.analyzePattern(Eigen_mat);
		solver.factorize(Eigen_mat);
		if (solver.info() != Eigen::Success) {
			throw std::runtime_error("Failed to perfom LU factorization");
		}

		Eigen::VectorXd Eigen_x = solver.solve(Eigen_b);
		if (solver.info() != Eigen::Success) {
			throw std::runtime_error("Failed to solve the linear system");
		}

		for (int i = 0; i < x.size(); ++i)
		{
			x[i] = Eigen_x(i);
		}
	}



}	// namespace AMG