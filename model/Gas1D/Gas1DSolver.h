#ifndef GAS1DSOLVER_H_
#define GAS1DSOLVER_H_

#include <iostream>
#include <cstdlib>
#include <map>

#include "model/AbstractSolver.hpp"
#include "method/sweep.h"

namespace gas1D
{
	class Gas1D;
	class Gas1D_simple;

	template <class modelType>
	class Gas1DSolver : public AbstractSolver<modelType>, public Sweep
	{
	private:
		void MiddleAppr(int current, int MZ, int key);
		void LeftBoundAppr(int MZ, int key);
		void RightBoundAppr(int MZ, int key);

	protected:		
		void construct_solution();
		void construction_from_fz(int N, int n, int key);
		void control();
		void doNextStep();
		void writeData();

		std::ofstream plot_P;
		std::ofstream plot_Q;

	public:
		Gas1DSolver(modelType* _model);
		Gas1DSolver(modelType* _model, int i);
		~Gas1DSolver();

		void fill();
		void start();
	};

	typedef Gas1DSolver<Gas1D> Gas1DSol;
	typedef Gas1DSolver<Gas1D_simple> Gas1DSolSimp;
};

#endif /* GAS1DSOLVER_H_ */