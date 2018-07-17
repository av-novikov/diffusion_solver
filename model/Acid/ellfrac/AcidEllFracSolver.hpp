#ifndef ACIDELLFRACSOLVER_HPP_
#define ACIDELLFRACSOLVER_HPP_

#include <fstream>

#include "model/AbstractSolver.hpp"
#include "model/Acid/ellfrac/AcidEllFracModel.hpp"
#include "method/ParalutionInterface.h"

namespace acidellfrac
{
	static const int poro_stencil = 7;
	static const int frac_stencil = 7;

	class AcidEllFracSolver : public AbstractSolver<AcidEllFrac>
	{
	protected:
		void control();
		void solveStep();
		void writeData();
		void doNextStep();

		double** jac;
		double* y;
		double* x;

		std::array<double, var_poro_size> averVal, averValPrev, dAverVal;
		double err_newton;
		int step_idx;

		std::ofstream S, P, qcells, pvd_frac, pvd_poro;
		ParSolver solver;

		int strNum;
		int* ind_i;
		int* ind_j;
		double* a;
		int* ind_rhs;
		double* rhs;
		int* cols;
		// Number of non-zero elements in sparse matrix
		int elemNum;
		int options[4];
		int repeat;

		void checkStability();
		void computeJac();
		void fill();
		void copySolution(const paralution::LocalVector<double>& sol);
		void fillIndices()
		{
			int counter = 0;
			for (int i = 0; i < strNum; i++)
				ind_rhs[i] = i;
		};
	public:
		AcidEllFracSolver(AcidEllFrac* _model);
		~AcidEllFracSolver();

		void start();
	};
}

#endif /* ACIDELLFRACSOLVER_HPP_ */
