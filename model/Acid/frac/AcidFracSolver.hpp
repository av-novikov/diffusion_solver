#ifndef ACIDFRACSOLVER_HPP_
#define ACIDFRACSOLVER_HPP_

#include <fstream>

#include "model/AbstractSolver.hpp"
#include "model/Acid/frac/AcidFracModel.hpp"
#include "method/ParalutionInterface.h"

namespace acidfrac
{
	static const int poro_stencil = 3;
	static const int frac_stencil = 7;

	class AcidFracSolver : public AbstractSolver<AcidFrac>
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

		std::ofstream S, P, qcells;
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
		AcidFracSolver(AcidFrac* _model);
		~AcidFracSolver();

		void start();
	};
}

#endif /* ACIDFRACSOLVER_HPP_ */
