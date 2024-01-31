#ifndef AMG_ML_LEVEL_HPP
#define AMG_ML_LEVEL_HPP

#include "types.h"
#include "matrix.h"
#include "vector.h"

namespace AMG {
	class Level
	{
	public:
		Level() : A{ nullptr }, P{ nullptr } {}
		~Level() {
			delete A;
			delete P;
		}

		CSRMatrix* A;
		CSRMatrix* P;
		Vector x;
		Vector b;
		Vector tmp;
	};	// class level




}	// namespace AMG
#endif


