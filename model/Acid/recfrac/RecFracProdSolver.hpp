#ifndef RECFRACPRODSOLVER_HPP_
#define RECFRACPRODSOLVER_HPP_

#include <fstream>

#include "model/AbstractSolver.hpp"
#include "model/Acid/recfrac/RecFracProd.hpp"
#include "method/ParalutionInterface.h"

namespace acidrecfrac_prod
{
	static const int stencil = 5;

	class RecFracProdSolver : public AbstractSolver<RecFracProd>
	{
	protected:
        
		void control();
		bool solveSmartStep();
        void writeData();
        bool doNextSmartStep();

        // garbage
        void solveStep() {};
        void doNextStep() {};

		double** jac;
		double* y;
		double* x;

		std::array<double, var_size> averVal, averValPrev, dAverVal;
		double err_newton_first, err_newton;
		int step_idx;

		std::ofstream S, P, qcells, pvd;
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
        RecFracProdSolver(RecFracProd* _model);
		~RecFracProdSolver();

		void start();
	};
}

#endif /* RECFRACPRODSOLVER_HPP_ */
