#ifndef GASOILELLIPTICSOLVER_HPP_
#define GASOILELLIPTICSOLVER_HPP_

#include <iostream>
#include <cstdlib>
#include <map>

#include "method/sweep.h"
#include "model/AbstractSolver.hpp"
#include "method/ParalutionInterface.h"
#include "model/GasOil_Elliptic/GasOil_Elliptic.hpp"

namespace gasOil_elliptic
{
	class GasOilEllipticSolver : public AbstractSolver<GasOil_Elliptic>
	{
	protected:
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

		void control();
		void doNextStep();
		void solveStep();
		void writeData();

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

		// Sparse matrix solver
		ParSolver solver;

		void fill();
		void fillIndices();
		void copySolution(const paralution::LocalVector<double>& sol);

		// Coordinate form of sparse matrix & dense vector
		int* ind_i;
		int* ind_j;
		double* a;
		int* ind_rhs;
		double* rhs;
		// Number of non-zero elements in sparse matrix
		int elemNum;

	public:
		GasOilEllipticSolver(GasOil_Elliptic* _model);
		~GasOilEllipticSolver();

		void start();
	};
};

#endif /* GASOILELLIPTICSOLVER_HPP_ */
