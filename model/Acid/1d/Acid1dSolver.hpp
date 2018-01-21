#ifndef ACID1DSOLVER_HPP_
#define ACID1DSOLVER_HPP_

#include "model/Basic1d/Basic1dSolver.hpp"
#include "model/Acid/1d/Acid1d.hpp"

#include <array>

namespace acid1d
{
	class Acid1dSolver : public basic1d::Basic1dSolver<acid1d::Acid1d>
	{
	private:
		void MiddleAppr(int current, int MZ, int key);
		void LeftBoundAppr(int MZ, int key);
		void RightBoundAppr(int MZ, int key);
	protected:
		void solveStep();
		void writeData();
		void checkStability();

		void construction_from_fz(int N, int n, int key);
		static const int var_size = Variable::size;
		static const int size = Variable::size;

		std::array<double, var_size> averVal, averValPrev, dAverVal;
		double err_newton;

		std::ofstream mat_a, mat_b, mat_c, rhs_os;
		void writeMatrixes();

		std::ofstream S, P, qcells;
	public:
		Acid1dSolver(acid1d::Acid1d* _model);
		~Acid1dSolver();

		void start();
	};
};

#endif /* ACID1DSOLVER_HPP_ */