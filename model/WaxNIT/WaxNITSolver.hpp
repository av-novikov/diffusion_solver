#ifndef WAXNITSOLVER_HPP_
#define WAXNITSOLVER_HPP_

#include "model/Basic2d/Basic2dSolver.hpp"
#include "model/WaxNIT/WaxNIT.hpp"
#include "method/ParalutionInterface.h" 

namespace wax_nit
{
	class WaxNITSolver : public basic2d::Basic2dSolver<WaxNIT>
	{
	public:
		typedef ParSolver::Vector Vector;
	private:
		//void MiddleAppr(int current, int MZ, int key);
		//void LeftBoundAppr(int MZ, int key);
		//void RightBoundAppr(int MZ, int key);
		void setWaxTemperature();
	protected:
		void solveStep();
		void writeData();
		//void construction_from_fz(int N, int n, int key);
		void checkStability();
		
		static const int var_size = Variable::size - 2;
		static const int size = Variable::size;

		std::array<double, var_size> averVal, averValPrev, dAverVal;
		double err_newton;
		int step_idx;

		std::ofstream S;
		std::ofstream P;
		std::ofstream poro;
		std::ofstream T;
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
		WaxNITSolver(WaxNIT* _model);
		~WaxNITSolver();

		void start();
	};
};

#endif /* WAXNITSOLVER_HPP_ */