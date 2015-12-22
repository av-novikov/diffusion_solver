#ifndef GAS1DSOLVER_H_
#define GAS1DSOLVER_H_

#include <iostream>
#include <cstdlib>
#include <map>

#include "model/AbstractSolver.hpp"
#include "method/sweep.h"
#include "model/Gas1D/Gas1D.h"

namespace gas1D
{
	class Gas1DSolver : public AbstractSolver<Gas1D>, public Sweep
	{
	private:
		void MiddleAppr(int current, int MZ, int key);
		void LeftBoundAppr(int MZ, int key);
		void RightBoundAppr (int MZ, int key);

	protected:		
		void construct_solution();
		void construction_from_fz(int N, int n, int key);
		void control();
		void doNextStep();
		void writeData();

		std::ofstream plot_P;
		std::ofstream plot_Q;

	public:
		Gas1DSolver(Gas1D* _model);
		Gas1DSolver(Gas1D* _model, int i);
		~Gas1DSolver();

		void fill();
		void start();
	};
};

#endif /* GAS1DSOLVER_H_ */