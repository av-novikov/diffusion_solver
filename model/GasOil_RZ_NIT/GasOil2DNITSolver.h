#ifndef GASOIL2DNITSOLVER_H_
#define GASOIL2DNITSOLVER_H_

#include <iostream>

#include "model/AbstractSolver.hpp"
#include "method/sweep.h"
#include "model/GasOil_RZ_NIT/GasOil_RZ_NIT.h"

namespace gasOil_rz_NIT
{
	class GasOil2DNITSolver : public AbstractSolver<GasOil_RZ_NIT>, public Sweep
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
		void writeData();

		std::ofstream plot_Sdyn;
		std::ofstream plot_Pdyn;
		std::ofstream plot_Tdyn;
		
		double T_dim;

	public:
		GasOil2DNITSolver(GasOil_RZ_NIT* _model);
		~GasOil2DNITSolver();
	};

};

#endif /* GASOIL2DNITSOLVER_H_ */