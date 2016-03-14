#ifndef PARPERFNITSOLVER_H_
#define PARPERFNITSOLVER_H_

#include <iostream>
#include <cstdlib>
#include <map>

#include "model/cells/stencils/Stencil.h"
#include "model/AbstractSolver.hpp"
#include "method/ParalutionInterface.h"
#include "method/sweep.h"

#include "model/3D/Perforation/GasOil_Perf_NIT.h"

namespace gasOil_perf_nit
{
	class ParPerfNITSolver : public AbstractSolver<GasOil_Perf_NIT>
	{
	protected:
		std::ofstream plot_Tdyn;
		std::ofstream plot_Sdyn;
		std::ofstream plot_Pdyn;
		std::ofstream plot_qcells;

		double T_dim;
		
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

		void fill(int key);
		void fillIndices(int key);
		void copySolution(const paralution::LocalVector<double>& sol, int key);

		// Numerical stencils for matrix filling
		UsedStencils<GasOil_Perf_NIT>* stencils;

		// Sparse matrix solvers
		ParSolver pres_solver;
		ParSolver temp_solver;

		// Coordinate form of sparse matrix & dense vector
		int* ind_i;
		int* ind_j;
		double* a;
		int* ind_rhs;
		double* rhs;
		//
		int* tind_i;
		int* tind_j;
		double* ta;
		int* tind_rhs;
		double* trhs;

		// Number of non-zero elements in sparse matrix
		int presElemNum;
		int tempElemNum;

	public:
		ParPerfNITSolver(GasOil_Perf_NIT* _model);
		~ParPerfNITSolver();

		void start();
	};
}

#endif /* PARPERFNITSOLVER_H_ */
