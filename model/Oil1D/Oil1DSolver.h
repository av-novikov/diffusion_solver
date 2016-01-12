#ifndef OIL1DSOLVER_H_
#define OIL1DSOLVER_H_

#include "model/AbstractSolver.hpp"
#include "method/sweep.h"
#include "model/Oil1D/Oil1D.h"

namespace oil1D
{
	class Oil1DSolver : public AbstractSolver<Oil1D>, public Sweep
	{
	private:
		void MiddleAppr(int current, int MZ, int key);
		void LeftBoundAppr(int MZ, int key);
		void RightBoundAppr (int MZ, int key);

	protected:
		void construction_from_fz(int N, int n, int key);
		void control();
		void doNextStep();
		void writeData();

		std::ofstream plot_Pdyn;

	public:
		Oil1DSolver(Oil1D* _model);
		~Oil1DSolver();
	};
};

#endif /* OIL1DSOLVER_H_ */