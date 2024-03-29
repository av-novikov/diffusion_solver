#ifndef GASOIL3DSOLVER_H_
#define GASOIL3DSOLVER_H_

#include <iostream>
#include <cstdlib>
#include <map>

#include "model/AbstractSolver.hpp"
#include "method/sweep.h"
#include "model/3D/GasOil_3D/GasOil_3D.h"

namespace gasOil_3d
{
	class GasOil3DSolver : public AbstractSolver<GasOil_3D>, public Sweep
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

		std::ofstream plot_Sdyn;
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
		//void fillGrad(double mult);
		//void fillGess();

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
		inline int getCalcIdx(int idx)
		{
			if(idx < 0)
				return 2 * model->cellsNum_phi * (model->cellsNum_z+2) + idx;
			else if(idx > 2 * model->cellsNum_phi * (model->cellsNum_z+2) )
				return idx - 2 * model->cellsNum_phi * (model->cellsNum_z+2);
			else
				return idx;
		};

	public:
		GasOil3DSolver(GasOil_3D* _model);
		~GasOil3DSolver();
	};
};

#endif /* GASOIL3DSOLVER_H_ */