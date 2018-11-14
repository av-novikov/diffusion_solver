#ifndef ACIDRECFRACSOLVER_HPP_
#define ACIDRECFRACSOLVER_HPP_

#include <fstream>

#include "model/AbstractSolver.hpp"
#include "model/Acid/recfrac/AcidRecFracModel.hpp"
#include "method/ParalutionInterface.h"

namespace acidrecfrac
{
	static const int poro_stencil = 7;
	static const int frac_stencil = 7;

	class AcidRecFracSolver : public AbstractSolver<AcidRecFrac>
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
		double err_newton_first, err_newton;
		int step_idx;

		std::ofstream S, P, qcells, pvd_frac, pvd_poro, trans;
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

        void analyzeNewtonConvergence();
		void checkVariables();
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
		AcidRecFracSolver(AcidRecFrac* _model);
		~AcidRecFracSolver();

		void start();
	};
}

#endif /* ACIDRECFRACSOLVER_HPP_ */
