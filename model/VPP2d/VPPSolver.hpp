#ifndef VPPSOLVER_HPP_
#define VPPSOLVER_HPP_

#include <iostream>
#include <cstdlib>
#include <map>

#include "model/AbstractSolver.hpp"
#include "method/sweep.h"
#include "model/VPP2d/VPP2d.hpp"

namespace vpp2d
{
	class VPPSolver : public AbstractSolver<VPP2d>, public Sweep
	{
	private:
		void MiddleAppr(int current, int MZ, int key);
		void LeftBoundAppr(int MZ, int key);
		void RightBoundAppr(int MZ, int key);

		void setLeftAppr(const int i, const int idx);
		void setMiddleAppr(const int i, const int idx);
		void setRightAppr(const int i, const int idx);

	protected:

		void construction_from_fz(int N, int n, int key);
		void control();
		void doNextStep();
		void solveStep();
		void writeData();

		std::ofstream plot_P;
		std::ofstream plot_S;
		std::ofstream plot_C;
		std::ofstream plot_qcells;

		int n;
		bool isChange;
		TVector<double> dq;
		TVector<double> q;
		TMatrix<double> dpdq;
		TMatrix<double> mat;
		TVector<double> b;

		void fillq();
		void fillDq();
		void filldPdQ(double mult);
		void solveSystem();
		void solveDq(double mult);

		inline void printWellRates()
		{
			double DQ = model->Q_sum;
			std::map<int, double>::iterator it;
			int k = 0;
			for (it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
			{
				std::cout << "Rate in " << it->first << " = " << it->second * model->Q_dim * 86400.0 << "\t";
				std::cout << "Press in " << it->first << " = " << model->cells[it->first].u_next.p << std::endl;
				DQ -= it->second;
				k++;
			}
			std::cout << "Summary rate deviation = " << DQ * model->Q_dim * 86400.0 << std::endl;
			std::cout << std::endl;
		};

		double** jac;

	public:
		VPPSolver(VPP2d* _model);
		~VPPSolver();
	};
};

#endif /* VPPSOLVER_HPP_ */
