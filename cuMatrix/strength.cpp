#include "matrix.h"

using namespace AMG;
/* declarations */
CSRMatrix* Classical_strength(CSRMatrix* A, double theta, int num_variables, int* variables);
CSRMatrix* Modclassical_strength(CSRMatrix* A, double theta, int num_variables, int* variables);


/* definitions */
CSRMatrix* Classical_strength(CSRMatrix* A, double theta, int num_variables, int* variables)
{
    int start, end, col;
    double val;
    double row_scale;
    double threshold;
    double diag;

    if (!A->sorted)
    {
        A->sort();
    }
    if (!A->diag_first)
    {
        A->move_diag();
    }

    CSRMatrix* S = new CSRMatrix(A->n_rows, A->n_cols);
    S->idx2.resize(A->nnz);
    S->vals.resize(A->nnz);

    S->nnz = 0;
    S->idx1[0] = 0;
    for (int i = 0; i < A->n_rows; i++)
    {
        // Always add the diagonal 
        start = A->idx1[i];
        end = A->idx1[i + 1];
        if (end - start)
        {
            if (A->idx2[start] == i)
            {
                diag = A->vals[start];
                S->idx2[S->nnz] = A->idx2[start];
                S->vals[S->nnz] = diag;
                S->nnz++;
                start++;
            }
            else
            {
                diag = 0.0;
            }

            if (num_variables == 1)
            {
                if (diag < 0.0) // find max off-diag value in row
                {
                    row_scale = -RAND_MAX;
                    for (int j = start; j < end; j++)
                    {
                        val = A->vals[j];
                        if (val > row_scale)
                        {
                            row_scale = val;
                        }
                    }
                }
                else // find min off-diag value in row
                {
                    row_scale = RAND_MAX;
                    for (int j = start; j < end; j++)
                    {
                        val = A->vals[j];
                        if (val < row_scale)
                        {
                            row_scale = val;
                        }
                    }
                }
            }
            else
            {
                if (diag < 0.0) // find max off-diag value in row
                {
                    row_scale = -RAND_MAX;
                    for (int j = start; j < end; j++)
                    {
                        col = A->idx2[j];
                        if (variables[i] == variables[col])
                        {
                            val = A->vals[j];
                            if (val > row_scale)
                            {
                                row_scale = val;
                            }
                        }
                    }
                }
                else // find min off-diag value in row
                {
                    row_scale = RAND_MAX;
                    for (int j = start; j < end; j++)
                    {
                        col = A->idx2[j];
                        if (variables[i] == variables[col])
                        {
                            val = A->vals[j];
                            if (val < row_scale)
                            {
                                row_scale = val;
                            }
                        }
                    }
                }

            }

            // Multiply row magnitude by theta
            threshold = row_scale * theta;

            // Add off-diagonals greater than threshold
            if (num_variables == 1)
            {
                if (diag < 0)
                {
                    for (int j = start; j < end; j++)
                    {
                        val = A->vals[j];
                        if (val > threshold)
                        {
                            S->idx2[S->nnz] = A->idx2[j];
                            S->vals[S->nnz] = val;
                            S->nnz++;
                        }
                    }
                }
                else
                {
                    for (int j = start; j < end; j++)
                    {
                        val = A->vals[j];
                        if (val < threshold)
                        {
                            S->idx2[S->nnz] = A->idx2[j];
                            S->vals[S->nnz] = val;
                            S->nnz++;
                        }
                    }
                }
            }
            else
            {
                if (diag < 0)
                {
                    for (int j = start; j < end; j++)
                    {
                        col = A->idx2[j];
                        if (variables[i] == variables[col])
                        {
                            val = A->vals[j];
                            if (val > threshold)
                            {
                                S->idx2[S->nnz] = col;
                                S->vals[S->nnz] = val;
                                S->nnz++;
                            }
                        }
                    }
                }
                else
                {
                    for (int j = start; j < end; j++)
                    {
                        col = A->idx2[j];
                        if (variables[i] == variables[col])
                        {
                            val = A->vals[j];
                            if (val < threshold)
                            {
                                S->idx2[S->nnz] = col;
                                S->vals[S->nnz] = val;
                                S->nnz++;
                            }
                        }
                    }
                }
            }
        }
        S->idx1[i + 1] = S->nnz;
    }
    S->idx2.resize(S->nnz);
    S->idx2.shrink_to_fit();
    S->vals.resize(S->nnz);
    S->vals.shrink_to_fit();

    return S;
}
CSRMatrix* Modclassical_strength(CSRMatrix* A, double theta, int num_variables, int* variables)
{
    int start, end, col;
    double val;
    double row_scale;
    double threshold;
    double diag;

    if (!A->sorted)
    {
        A->sort();
    }
    if (!A->diag_first)
    {
        A->move_diag();
    }

    CSRMatrix* S = new CSRMatrix(A->n_rows, A->n_cols);
    S->idx2.resize(A->nnz);
    S->vals.resize(A->nnz);

    S->nnz = 0;
    S->idx1[0] = 0;
    for (int i = 0; i < A->n_rows; i++) {
        start = A->idx1[i];
        end = A->idx1[i + 1];
        if (end - start)
        {
            // always add diagonal
            if (A->idx2[start] == i)
            {
                diag = A->vals[start];
                S->idx2[S->nnz] = A->idx2[start];
                S->vals[S->nnz] = diag;
                S->nnz++;
                start++;
            }
			else
			{
				diag = 0.0;
                S->idx2[S->nnz] = A->idx2[start];
                S->vals[S->nnz] = diag;
                S->nnz++;
                start++;
			}

            row_scale = -DBL_MAX;
            for (int j = start; j < end; j++) {
                val = A->vals[j];
                if (-val > row_scale)
                {
                    row_scale = -val;
                }
            }

            threshold = row_scale * theta;
            for (int j = start; j < end; j++) {
                val = A->vals[j];
                if (-val >= threshold)
                {
                    S->idx2[S->nnz] = A->idx2[j];
                    S->vals[S->nnz] = val;
                    S->nnz++;
                }
            }
        }
        S->idx1[i + 1] = S->nnz;
    }
    S->idx2.resize(S->nnz);
    S->idx2.shrink_to_fit();
    S->vals.resize(S->nnz);
    S->vals.shrink_to_fit();

    return S;
}

CSRMatrix* CSRMatrix::strength(strength_t strength_type,
    double theta, int num_variables, int* variables)
{
    switch (strength_type)
    {
    case Classical:
        return Classical_strength(this, theta, num_variables, variables);
    case ModClassical:
        return Modclassical_strength(this, theta, num_variables, variables);
    default:
        return Modclassical_strength(this, theta, num_variables, variables);
    }
    return NULL;
}













