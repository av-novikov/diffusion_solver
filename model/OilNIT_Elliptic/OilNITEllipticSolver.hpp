#ifndef OILNITELLIPTICSOLVER_HPP_
#define OILNITELLIPTICSOLVER_HPP_

#include <iostream>
#include <cstdlib>
#include <map>

#include "method/sweep.h"
#include "model/AbstractSolver.hpp"
#include "method/ParalutionInterface.h"
#include "model/OilNIT_Elliptic/OilNIT_Elliptic.hpp"

namespace oilnit_elliptic
{
	class OilNITEllipticSolver : public AbstractSolver<OilNIT_Elliptic>
	{
	protected:
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
			double DQ = model->Q_sum * model->Q_dim * 86400.0;
			std::map<int, double>::iterator it;
			int k = 0;
			for (it = model->Qcell_ellipse.begin(); it != model->Qcell_ellipse.end(); ++it)
			{
				std::cout << "Rate in " << it->first << " = " << it->second * model->Q_dim * 86400.0 << "\t";
				std::cout << "Press in " << it->first << " = " << model->cells[it->first].u_next.p * model->P_dim / BAR_TO_PA << std::endl;
				DQ -= it->second;
				k++;
			}
			std::cout << "Summary rate deviation = " << DQ * model->Q_dim * 86400.0 << std::endl;
			std::cout << std::endl;
		};
		inline const std::vector<int> getMatrixStencil(const Cell& cell)
		{
			std::vector<int> stencil_idx;

			if (cell.type == MIDDLE)
			{
				stencil_idx.resize(7);
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
				return stencil_idx;
			}
			else if (cell.type == MIDDLE_SIDE)
			{
				stencil_idx.resize(6);
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
				return stencil_idx;
			}
			else if (cell.type == RIGHT)
			{
				stencil_idx.resize(2);
				stencil_idx[0] = cell.num;
				stencil_idx[1] = cell.num - model->cellsNum_z - 2;
				return stencil_idx;
			}
			else if (cell.type == TOP)
			{
				stencil_idx.resize(2);
				stencil_idx[0] = cell.num;
				stencil_idx[1] = cell.num + 1;
				return stencil_idx;
			}
			else if (cell.type == BOTTOM)
			{
				stencil_idx.resize(2);
				stencil_idx[0] = cell.num;
				stencil_idx[1] = cell.num - 1;
				return stencil_idx;
			}
			else
			{
				stencil_idx.resize(2);
				stencil_idx[0] = model->cellsNum + cell.num;
				stencil_idx[1] = model->nebrMap[cell.num].first;
				return stencil_idx;
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
		OilNITEllipticSolver(OilNIT_Elliptic* _model);
		~OilNITEllipticSolver();

		void start();
	};
};

#endif /* OILNITELLIPTICSOLVER_HPP_ */
