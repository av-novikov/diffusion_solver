#ifndef ACID2DSOLVER_HPP_
#define ACID2DSOLVER_HPP_

#include "model/Basic2d/Basic2dSolver.hpp"
#include "model/Acid/2d/Acid2d.hpp"

#include <array>

namespace acid2d
{
	class Acid2dSolver : public basic2d::Basic2dSolver<acid2d::Acid2d>
	{
	private:
		void MiddleAppr(int current, int MZ, int key);
		void LeftBoundAppr(int MZ, int key);
		void RightBoundAppr(int MZ, int key);
	protected:
		void solveStep();
		void writeData();
		void checkStability() {};

		void construction_from_fz(int N, int n, int key);
		static const int var_size = Variable::size;
		static const int size = Variable::size;

		std::array<double, var_size> averVal, averValPrev, dAverVal;
		double err_newton;

		std::ofstream mat_a, mat_b, mat_c, rhs;
		void writeMatrixes();

		std::ofstream S;
		std::ofstream P;
		std::ofstream qcells;
	public:
		Acid2dSolver(acid2d::Acid2d* _model);
		~Acid2dSolver();
	};
}

#endif /* ACID2DSOLVER_HPP_ */