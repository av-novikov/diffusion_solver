#ifndef OILRZSOLVER_H_
#define OILRZSOLVER_H_

#include <iostream>
#include <cstdlib>
#include <map>

#include "model/AbstractSolver.hpp"
#include "method/sweep.h"
#include "model/Oil_RZ/Oil_RZ.h"

namespace oil_rz
{
	class OilRZSolver : public AbstractSolver<Oil_RZ>, public Sweep
	{
	private:
		void MiddleAppr(int current, int MZ, int key);
		void LeftBoundAppr(int MZ, int key);
		void RightBoundAppr (int MZ, int key);

	protected:
		void TopAppr(int i, int key);
		void BottomAppr(int i, int key);
		
		void construction_from_fz(int N, int n, int key);
		void control();
		void doNextStep();
		void solveStep();
		void writeData();

		std::ofstream plot_Pdyn;
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
			std::map<int,double>::iterator it;
			int k = 0;
			for(it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
			{
				std::cout << "Rate in " << it->first << " = " << it->second * model->Q_dim * 86400.0 << "\t";
				std::cout << "Press in " << it->first << " = " << model->cells[ it->first ].u_next.p << std::endl;
				DQ -= it->second;
				k++;
			}
			std::cout << "Summary rate deviation = " << DQ * model->Q_dim * 86400.0 << std::endl;
			std::cout << std::endl;
		};

	public:
		OilRZSolver(Oil_RZ* _model);
		~OilRZSolver();
	};
};

#endif /* OILRZSOLVER_H_ */