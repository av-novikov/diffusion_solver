#ifndef OIL1DNITSOLVER_H_
#define OIL1DNITSOLVER_H_

#include "model/AbstractSolver.hpp"
#include "method/sweep.h"
#include "model/Oil1D_NIT/Oil1D_NIT.h"

namespace oil1D_NIT
{
	class Oil1DNITSolver : public AbstractSolver<Oil1D_NIT>, public Sweep
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

		double T_dim;

		std::ofstream plot_Pdyn;
		std::ofstream plot_Tdyn;
		std::ofstream plot_qcells;

	public:
		Oil1DNITSolver(Oil1D_NIT* _model);
		~Oil1DNITSolver();
	};
};

#endif /* OIL1DNITSOLVER_H_ */