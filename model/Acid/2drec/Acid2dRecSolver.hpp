#ifndef ACID2DRECSOLVER_HPP_
#define ACID2DRECSOLVER_HPP_

#include <fstream>

#include "model/AbstractSolver.hpp"
#include "model/Acid/2drec/Acid2dRecModel.hpp"
#include "method/ParalutionInterface.h"
#include "method/HypreInterface.hpp"

namespace acid2drec
{
	static const int stencil = 5;

	template<class SolType>
	class Acid2dRecSolver : public AbstractSolver<Acid2dRecModel>
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
		double rel_tol, abs_tol;

		std::ofstream S, P, qcells, pvd;
		SolType solver;

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
		template<class VecType>
		void copySolution(const VecType& sol);
		void fillIndices()
		{
			int counter = 0;
			for (int i = 0; i < strNum; i++)
				ind_rhs[i] = i;
		};

	public:
		Acid2dRecSolver(Acid2dRecModel* _model);
		~Acid2dRecSolver();

		void start();
	};
}

#endif /* ACID2DRECSOLVER_HPP_ */
