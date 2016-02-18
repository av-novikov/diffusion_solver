#ifndef PAR3DSOLVER_H_
#define PAR3DSOLVER_H_

#include <iostream>
#include <cstdlib>
#include <map>

#include "model/AbstractSolver.hpp"
#include "method/ParalutionInterface.h"
#include "method/sweep.h"

#include "model/3D/GasOil_3D/GasOil_3D.h"

namespace gasOil_3d
{
	class Par3DSolver : public AbstractSolver<GasOil_3D>, public ParSolver
	{
	protected:
		std::ofstream plot_Sdyn;
		std::ofstream plot_Pdyn;
		std::ofstream plot_qcells;

		void control();
		void doNextStep();
		void solveStep();
		void writeData();

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

		int counter;
		void Solve();
		void fillLeft();
		void fillMiddle();
		void fillRight();
		void fillTop(int i);
		void fillBottom(int i);

	public:
		Par3DSolver(GasOil_3D* _model);
		~Par3DSolver();
	};
}

#endif /* PAR3DSOLVER_H_ */
