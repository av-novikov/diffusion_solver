#ifndef WAXNIT1DSOLVER_HPP_
#define WAXNIT1DSOLVER_HPP_

#include "model/Basic1d/Basic1dSolver.hpp"
#include "model/WaxNIT/1d/WaxNIT1d.hpp"
#include "method/ParalutionInterface.h" 

namespace wax_nit1d
{
	class WaxNIT1dSolver : public basic1d::Basic1dSolver<WaxNIT1d>
	{
	public:
		typedef ParSolver::Vector Vector;
	private:
		//void MiddleAppr(int current, int MZ, int key);
		//void LeftBoundAppr(int MZ, int key);
		//void RightBoundAppr(int MZ, int key);
	protected:
		void solveStep();
		void writeData();
		//void construction_from_fz(int N, int n, int key);
		void checkStability();
		
		static const int var_size = Variable::size;
		static const int size = Variable::size;

		std::array<double, var_size> averVal, averValPrev, dAverVal;
		double err_newton;
		int step_idx;

		std::ofstream S;
		std::ofstream P;
		std::ofstream poro;
		std::ofstream qcells;
		std::ofstream pvd;
		
		//std::ofstream mat_a, mat_b, mat_c, rhs_os;
		//void writeMatrixes();

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

				int neighbor[stencil - 1];
				model->getNeighborIdx(cell.num, neighbor);

				stencil_idx[0] = cell.num;
				for (int i = 0; i < stencil - 1; i++)
					stencil_idx[i + 1] = neighbor[i];
			}
			else if (cell.type == Type::WELL_LAT)
			{
				stencil_idx.resize(Lstencil);

				stencil_idx[0] = cell.num;
				stencil_idx[1] = cell.num + 1;
				stencil_idx[2] = cell.num + 2;
			}
			else if (cell.type == Type::RIGHT)
			{
				stencil_idx.resize(Rstencil);

				stencil_idx[0] = cell.num;
				stencil_idx[1] = cell.num - 1;
				stencil_idx[2] = cell.num - 2;
			}

			return stencil_idx;
		};
	public:
		WaxNIT1dSolver(WaxNIT1d* _model);
		~WaxNIT1dSolver();

		void start();
	};
};

#endif /* WAXNIT1DSOLVER_HPP_ */