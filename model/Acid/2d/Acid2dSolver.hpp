#ifndef ACID2DSOLVER_HPP_
#define ACID2DSOLVER_HPP_

#include "model/Basic2d/Basic2dSolver.hpp"
#include "model/Acid/2d/Acid2d.hpp"

namespace acid2d
{
	class Acid2dSolver : public basic2d::Basic2dSolver<acid2d::Acid2d>
	{
	private:
		void MiddleAppr(int current, int MZ, int key) {};
		void LeftBoundAppr(int MZ, int key) {};
		void RightBoundAppr(int MZ, int key) {};
	protected:
		void solveStep() {};
		void writeData() {};
		void construction_from_fz(int N, int n, int key) {};
		void checkStability() {};

		static const int var_size = Variable::size - 1;
		static const int size = Variable::size;

		std::ofstream S;
		std::ofstream P;
		std::ofstream qcells;
	public:
		Acid2dSolver(acid2d::Acid2d* _model);
		~Acid2dSolver();
	};
}

#endif /* ACID2DSOLVER_HPP_ */