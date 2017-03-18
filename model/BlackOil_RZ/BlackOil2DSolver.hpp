#ifndef BLACKOIL2DSOLVER_HPP_
#define BLACKOIL2DSOLVER_HPP_

#include "model/Basic2d/Basic2dSolver.hpp"
#include "model/BlackOil_RZ/BlackOil_RZ.hpp"

namespace blackoil_rz
{
	class BlackOil2dSolver : public basic2d::Basic2dSolver<BlackOil_RZ>
	{
	private:
		void MiddleAppr(int current, int MZ, int key);
		void LeftBoundAppr(int MZ, int key);
		void RightBoundAppr(int MZ, int key);
	protected:
		void solveStep();
		void writeData();
		void construction_from_fz(int N, int n, int key);

		static const int var_size = Variable::size - 1;
		static const int size = Variable::size;

		std::ofstream S;
		std::ofstream P;
		std::ofstream qcells;
	public:
		BlackOil2dSolver(BlackOil_RZ* _model);
		~BlackOil2dSolver();
	};
};

#endif BLACKOIL2DSOLVER_HPP_