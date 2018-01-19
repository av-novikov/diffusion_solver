#ifndef WAXNITSOLVER_HPP_
#define WAXNITSOLVER_HPP_

#include "model/Basic2d/Basic2dSolver.hpp"
#include "model/WaxNIT/WaxNIT.hpp"

namespace wax_nit
{
	class WaxNITSolver : public basic2d::Basic2dSolver<WaxNIT>
	{
	private:
		void MiddleAppr(int current, int MZ, int key);
		void LeftBoundAppr(int MZ, int key);
		void RightBoundAppr(int MZ, int key);
	protected:
		void solveStep();
		void writeData();
		void construction_from_fz(int N, int n, int key);
		void checkStability();
		
		static const int var_size = Variable::size - 1;
		static const int size = Variable::size;

		std::array<double, var_size> averVal, averValPrev, dAverVal;
		double err_newton;

		std::ofstream S;
		std::ofstream P;
		std::ofstream poro;
		std::ofstream T;
		std::ofstream qcells;
		
		std::ofstream mat_a, mat_b, mat_c, rhs_os;
		void writeMatrixes();
	public:
		WaxNITSolver(WaxNIT* _model);
		~WaxNITSolver();
	};
};

#endif /* WAXNITSOLVER_HPP_ */