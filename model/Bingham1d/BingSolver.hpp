#ifndef BINGSOLVER_HPP_
#define BINGSOLVER_HPP_

#include <iostream>
#include <cstdlib>

#include "model/AbstractSolver.hpp"
#include "method/sweep.h"
#include "model/Bingham1d/Bingham1d.hpp"

namespace bing1d
{
	class Bing1dSolver : public AbstractSolver<Bingham1d>, public Sweep
	{
	private:
		void MiddleAppr(int current, int MZ, int key);
		void LeftBoundAppr(int MZ, int key);
		void RightBoundAppr(int MZ, int key);
	protected:

		void construction_from_fz(int N, int n, int key);
		void control();
		void doNextStep();
		void solveStep();
		void writeData();

		std::ofstream plot_Pdyn;
		std::ofstream plot_qcells;

		double** jac;
	public:
		Bing1dSolver(Bingham1d* _model);
		~Bing1dSolver();
	};
};

#endif /* BINGSOLVER_HPP_ */
