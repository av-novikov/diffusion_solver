#ifndef BLACKOILNITELLIPTICSOLVER_HPP_
#define BLACKOILNITELLIPTICSOLVER_HPP_

#include <iostream>
#include <cstdlib>
#include <map>

#include "method/sweep.h"
#include "model/AbstractSolver.hpp"
#include "model/BlackOilNIT_Elliptic/BlackOilNIT_Elliptic.hpp"

namespace blackoilnit_elliptic
{
	template <typename solType>
	class BlackOilNITEllipticSolver : public AbstractSolver<BlackOilNIT_Elliptic>
	{
	public:
		typedef typename solType::Vector Vector;
	protected:
		std::ofstream plot_Pdyn;
		std::ofstream plot_Tdyn;
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
		void solveTempStep();
		void writeData();

		inline void printWellRates()
		{
			double DQ = model->Q_sum_quater * model->Q_dim * 86400.0;
			std::map<int, double>::iterator it;
			int k = 0;
			for (it = model->Qcell_ellipse.begin(); it != model->Qcell_ellipse.end(); ++it)
			{
				std::cout << "Rate in " << it->first << " = " << it->second * model->Q_dim * 86400.0 << "\t";
				std::cout << "Press in " << it->first << " = " << model->wellCells[it->first].u_next.p * model->P_dim / BAR_TO_PA << std::endl;
				DQ -= it->second * model->Q_dim * 86400.0;
				k++;
			}
			std::cout << "Summary rate deviation = " << DQ * model->Q_dim * 86400.0 << std::endl;
			std::cout << std::endl;
		};
		inline std::vector<int> getMatrixStencil(const Cell& cell, const int val)
		{
			std::vector<int> stencil_idx;

			if (cell.type == Type::MIDDLE)
			{
				stencil_idx.resize(stencil);
				stencil_idx[0] = cell.num;
				// Special neighbor search for center cells
				if (cell.num % ((model->cellsNum_mu + 2) * (model->cellsNum_z + 2)) > model->cellsNum_z + 1)
					stencil_idx[1] = model->getCellIdx(cell.num, cell.num - model->cellsNum_z - 2);
				else
				{
					int nu_idx = cell.num / ((model->cellsNum_z + 2) * (model->cellsNum_mu + 2));
					stencil_idx[1] = model->getCellIdx(cell.num, model->cellsNum + cell.num - 2 * nu_idx * (model->cellsNum_z + 2) * (model->cellsNum_mu + 2));
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
			else if (cell.type == Type::MIDDLE_SIDE)
			{
				stencil_idx.resize(stencil - 1);
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
			else if (cell.type == Type::RIGHT)
			{
				stencil_idx.resize(Rstencil);
				stencil_idx[0] = cell.num;
				stencil_idx[1] = cell.num - model->cellsNum_z - 2;
				return stencil_idx;
			}
			else if (cell.type == Type::TOP)
			{
				stencil_idx.resize(Vstencil);
				stencil_idx[0] = cell.num;
				stencil_idx[1] = cell.num + 1;
				return stencil_idx;
			}
			else if (cell.type == Type::BOTTOM)
			{
				stencil_idx.resize(Vstencil);
				stencil_idx[0] = cell.num;
				stencil_idx[1] = cell.num - 1;
				return stencil_idx;
			}
			else
			{
				if (val == TEMP)
				{
					stencil_idx.resize(TLstencil);
					stencil_idx[0] = model->cellsNum + cell.num;
					stencil_idx[1] = model->nebrMap[cell.num].first;
					stencil_idx[2] = model->nebrMap[cell.num].second;
					return stencil_idx;
				}
				else if (val == PRES)
				{
					stencil_idx.resize(Lstencil);
					stencil_idx[0] = model->cellsNum + cell.num;
					stencil_idx[1] = model->nebrMap[cell.num].first;
					return stencil_idx;
				}
			}
		};

		// Sparse matrix solver
		solType pres_solver, temp_solver;

		void fill(const int val);
		void fillIndices();
		void copySolution(const Vector& sol, const int val);

		// Coordinate form of sparse matrix & dense vector
		int *ind_i, *ind_j, *tind_i, *tind_j, *ind_rhs, *cols, *t_cols;
		double *a, *rhs;
		// Number of non-zero elements in sparse matrix
		int elemNum, telemNum;

		std::vector<std::vector<double>> rateRatios;
		std::vector<double> rateRatiosAmongIntervals;

	public:
		BlackOilNITEllipticSolver(BlackOilNIT_Elliptic* _model);
		~BlackOilNITEllipticSolver();

		void start();
	};
};

#endif /* BLACKOILNITELLIPTICSOLVER_HPP_ */
