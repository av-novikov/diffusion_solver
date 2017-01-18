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
		inline const std::initializer_list<int> getMatrixStencil(const Cell& cell)
		{
			int stencil_idx[7];

			if (cell.type == MIDDLE)
			{
				stencil_idx[0] = cell.num;
				// Special neighbor search for center cells
				if (cell.num % ((model->cellsNum_mu + 2) * (model->cellsNum_z + 2)) > model->cellsNum_z + 1)
					stencil_idx[1] = model->getCellIdx(cell.num, cell.num - model->cellsNum_z - 2);
				else
				{
					int nu_idx = cell.num / ((model->cellsNum_z + 2) * (model->cellsNum_mu + 2));
					stencil_idx[1] = model->getCellIdx(cell.num, model->cellsNum + cell.num - (2 * nu_idx + 1) * (model->cellsNum_z + 2) * (model->cellsNum_mu + 2));
				}
				stencil_idx[2] = model->getCellIdx(cell.num, cell.num + model->cellsNum_z + 2);
				stencil_idx[3] = model->getCellIdx(cell.num, cell.num - 1);
				stencil_idx[4] = model->getCellIdx(cell.num, cell.num + 1);

				if (cell.num < (model->cellsNum_mu + 2) * (model->cellsNum_z + 2))
					stencil_idx[5] = model->getCellIdx(cell.num, cell.num +
					(model->cellsNum_mu + 2) * (model->cellsNum_z + 2) * (model->cellsNum_nu - 1));
				else
					stencil_idx[5] = model->getCellIdx(cell.num, cell.num -
					(model->cellsNum_mu + 2) * (model->cellsNum_z + 2));
				if (cell.num < (model->cellsNum_mu + 2) * (model->cellsNum_z + 2) * (model->cellsNum_nu - 1))
					stencil_idx[6] = model->getCellIdx(cell.num, cell.num +
					(model->cellsNum_mu + 2) * (model->cellsNum_z + 2));
				else
					stencil_idx[6] = model->getCellIdx(cell.num, cell.num -
					(model->cellsNum_mu + 2) * (model->cellsNum_z + 2) * (model->cellsNum_nu - 1));
				return{ stencil_idx[0], stencil_idx[1], stencil_idx[2], stencil_idx[3], 
						stencil_idx[4], stencil_idx[5], stencil_idx[6] };
			}
			else if (cell.type == MIDDLE_SIDE)
			{
				stencil_idx[0] = cell.num;
				stencil_idx[1] = model->getCellIdx(cell.num, cell.num + model->cellsNum_z + 2);
				stencil_idx[2] = model->getCellIdx(cell.num, cell.num - 1);
				stencil_idx[3] = model->getCellIdx(cell.num, cell.num + 1);

				if (cell.num < (model->cellsNum_mu + 2) * (model->cellsNum_z + 2))
					stencil_idx[4] = model->getCellIdx(cell.num, cell.num +
						(model->cellsNum_mu + 2) * (model->cellsNum_z + 2) * (model->cellsNum_nu - 1));
				else
					stencil_idx[4] = model->getCellIdx(cell.num, cell.num -
						(model->cellsNum_mu + 2) * (model->cellsNum_z + 2));
				if (cell.num < (model->cellsNum_mu + 2) * (model->cellsNum_z + 2) * (model->cellsNum_nu - 1))
					stencil_idx[5] = model->getCellIdx(cell.num, cell.num +
						(model->cellsNum_mu + 2) * (model->cellsNum_z + 2));
				else
					stencil_idx[5] = model->getCellIdx(cell.num, cell.num -
						(model->cellsNum_mu + 2) * (model->cellsNum_z + 2) * (model->cellsNum_nu - 1));
				return{ stencil_idx[0], stencil_idx[1], stencil_idx[2], stencil_idx[3],
					stencil_idx[4], stencil_idx[5] };
			}
			else if (cell.type == RIGHT)
			{
				stencil_idx[0] = cell.num;
				stencil_idx[1] = cell.num - model->cellsNum_z - 2;
				return{ stencil_idx[0], stencil_idx[1] };
			}
			else if (cell.type == TOP)
			{
				stencil_idx[0] = cell.num;
				stencil_idx[1] = cell.num + 1;
				return{ stencil_idx[0], stencil_idx[1] };
			}
			else if (cell.type == BOTTOM)
			{
				stencil_idx[0] = cell.num;
				stencil_idx[1] = cell.num - 1;
				return{ stencil_idx[0], stencil_idx[1] };
			}
			else
			{
				stencil_idx[0] = model->cellsNum + cell.num;
				stencil_idx[1] = model->nebrMap[cell.num].first;
				stencil_idx[2] = model->nebrMap[cell.num].second;
				return{ stencil_idx[0], stencil_idx[1], stencil_idx[2] };
			}
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
