#ifndef cuAMG_CORE_TYPES_H_
#define cuAMG_CORE_TYPES_H_

#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>

#include <vector>
#include <map>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <functional>
#include <set>

#include <cstdint>
#include <vector>
#include <stdexcept>

#include <random>

#define zero_tol 1e-16

namespace AMG {
	constexpr int TmpSelection = 4;
	constexpr int NewSelection = 3;
	constexpr int NewUnselection = 2;
	constexpr int Selected = 1;
	constexpr int Unselected = 0;
	constexpr int Unassigned = -1;
	constexpr int NoNeighbors = -2;

	using data_t = double;
	using index_t = int;
	enum strength_t { Classical, Modclassical };
	enum format_t { COO, CSR, CSC, BCOO, BSR, BSC };
	enum coarsen_t { RS, CLJP, Falgout, PMIS, HMIS };
	enum interp_t { Direct, ModClassical, Extended };
	enum agg_t { MIS };
	enum prolong_t { JacobiProlongation };
	enum relax_t { Jacobi, GaussSeidel, SOR, SSOR };

}	// namespace AMG
#endif
