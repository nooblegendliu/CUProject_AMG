#ifndef AMG_ML_MULTILEVEL_H_
#define AMG_ML_MULTILEVEL_H_

#include "types.h"
#include "matrix.h"
#include "vector.h"
#include "level.h"

#include "relax.h"


namespace AMG {

	class Multilevel
	{
	public:
		Multilevel(double _strong_threshold,
			strength_t _strength_type,
			relax_t _relax_type)
		{
			strong_threshold = _strong_threshold;
			strength_type = _strength_type;
			relax_type = _relax_type;
			// default settings
			num_smooth_sweeps = 5;
			relax_weight = 1.0;
			max_coarse = 500;
			max_levels = 25;

			store_residual = true;
		}

		virtual ~Multilevel()
		{
			for (std::vector<Level*>::iterator it = levels.begin();
				it != levels.end(); ++it)
			{
				delete* it;
			}
		}

		virtual void setup(CSRMatrix* Af) = 0;

		void setup_helper(CSRMatrix* Af)
		{
			printf("Strength %d\n", strength_type);
			int last_level = 0;

			levels.emplace_back(new Level());
			levels[0]->A = Af;
			levels[0]->A->sort();
			levels[0]->x.resize(Af->n_rows);
			levels[0]->b.resize(Af->n_rows);
			levels[0]->tmp.resize(Af->n_rows);
			levels[0]->P = nullptr;

			while (levels[last_level]->A->n_rows > max_coarse &&
				(max_levels == -1 || levels.size() < max_levels))
			{
				extend_hierarchy();
				last_level++;
				if (levels[last_level - 1]->A->n_rows == levels[last_level - 2]->A->n_rows) {
					break;
				}
			}

			num_levels = levels.size();
		}

		virtual void extend_hierarchy() = 0;

		









		void cycle(Vector& x, Vector& b, int level)
		{
			CSRMatrix* A = levels[level]->A;
			CSRMatrix* P = levels[level]->P;
			Vector& tmp = levels[level]->tmp;

			if (level == num_levels - 1)
			{
				Direct_solve(A, b, x);
			}
			else {
				levels[level + 1]->x.set_const_value(0.0);

				// relax
				switch (relax_type)
				{
				case Jacobi:
					jacobi(A, b, x, tmp, num_smooth_sweeps, relax_weight);
					break;
				case GaussSeidel:
					gaussSeidel(A, b, x, num_smooth_sweeps, relax_weight);
					break;
				case SOR:
					sor(A, x, b, tmp, num_smooth_sweeps, relax_weight);
					break;
				case SSOR:
					ssor(A, x, b, tmp, num_smooth_sweeps, relax_weight);
					break;
				default:
					jacobi(A, b, x, tmp, num_smooth_sweeps, relax_weight);
					break;
				}

				// residual
				A->residual(x, b, tmp);

				// restrict
				P->mult_T(tmp, levels[level + 1]->b);






			}
		}

	public:
		relax_t relax_type;
		strength_t strength_type;

		int num_smooth_sweeps;	//光滑迭代次数
		int max_coarse;			//最大最粗层级
		int max_levels;         //最大粗化层数

		double strong_threshold;	//强连接阈值
		double relax_weight;		//松弛因子

		bool store_residual;			//是否存储残差
		std::vector<double> residuals;	//残差

		std::vector<Level*> levels;	//层级
		int num_levels;				//层数
	};












}	// namespace AMG
#endif
