#ifndef ACID2DSOLVER_HPP_
#define ACID2DSOLVER_HPP_

#include "model/Basic2d/Basic2dSolver.hpp"
#include "model/Acid/2d/Acid2d.hpp"
#include "method/ParalutionInterface.h"

#include <array>

namespace acid2d
{
	class Acid2dSolver : public basic2d::Basic2dSolver<acid2d::Acid2d>
	{
	public:
		typedef ParSolver::Vector Vector;
	private:
		void MiddleAppr(int current, int MZ, int key);
		void LeftBoundAppr(int MZ, int key);
		void RightBoundAppr(int MZ, int key);
	protected:
		void solveStep();
		void writeData();
		void checkStability();

		void construction_from_fz(int N, int n, int key);
		static const int var_size = Variable::size - 1;
		static const int size = Variable::size;

		std::array<double, var_size> averVal, averValPrev, dAverVal;
		double err_newton;

		std::ofstream mat_a, mat_b, mat_c, rhs_os;
		void writeMatrixes();

		std::ofstream S, P, qcells;

		// Coordinate form of sparse matrix & dense vector
		int *ind_i, *ind_j, *ind_rhs, *cols;
		double *a, *rhs;
		// Number of non-zero elements in sparse matrix
		int elemNum;
		ParSolver solver;

		void fill();
		void fillIndices();
		void copySolution(const Vector& sol);
		inline std::vector<int> getMatrixStencil(const Cell& cell)
		{
			std::vector<int> stencil_idx;

			if (cell.type == Type::MIDDLE)
			{
				stencil_idx.resize(stencil);

				int neighbor[stencil-1];
				model->getNeighborIdx(cell.num, neighbor);

				stencil_idx[0] = cell.num;
				for (int i = 0; i < stencil-1; i++)
					stencil_idx[i+1] = neighbor[i];
			}
			else if (cell.type == Type::WELL_LAT)
			{
				stencil_idx.resize(Lstencil);

				stencil_idx[0] = cell.num;
				stencil_idx[1] = cell.num + model->cellsNum_z + 2;
				stencil_idx[2] = cell.num + 2 * model->cellsNum_z + 4;
			}
			else if (cell.type == Type::RIGHT)
			{
				stencil_idx.resize(Rstencil);

				stencil_idx[0] = cell.num;
				stencil_idx[1] = cell.num - model->cellsNum_z - 2;
				stencil_idx[2] = cell.num - 2 * model->cellsNum_z - 4;
			}
			else if (cell.type == Type::TOP)
			{
				stencil_idx.resize(Vstencil);

				stencil_idx[0] = cell.num;
				stencil_idx[1] = cell.num + 1;
			}
			else if (cell.type == Type::BOTTOM)
			{
				stencil_idx.resize(Vstencil);

				stencil_idx[0] = cell.num;
				stencil_idx[1] = cell.num - 1;
			}

			return stencil_idx;
		};
	public:
		Acid2dSolver(acid2d::Acid2d* _model);
		~Acid2dSolver();

		void start();
	};
}

#endif /* ACID2DSOLVER_HPP_ */