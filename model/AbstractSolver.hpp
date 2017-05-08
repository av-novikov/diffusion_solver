#ifndef ABSTRACTSOLVER_HPP_
#define ABSTRACTSOLVER_HPP_

#include <algorithm>
#include <cmath>
#include <iostream>
#include <array>

template <class modelType>
class AbstractSolver {
	protected:

		modelType* model;
		int size;
			
		int curTimePeriod;
		double Tt;
		double cur_t, cur_t_log;

		int idx1, idx2;

		double t_dim;

		int iterations;
		
		void copyIterLayer();
		void revertIterLayer();
		void copyTimeLayer();
		
		double convergance(int& ind, int& varInd);
		double averValue(int varInd);
		void averValue(std::array<double, modelType::var_size>& aver);
		
		virtual void writeData() = 0;
		virtual void control() = 0;
		virtual void doNextStep() = 0;
		virtual void solveStep() = 0;

		double NEWTON_STEP;
		double CHOP_MULT;
		double MAX_SAT_CHANGE;
		double CONV_W2, CONV_VAR;
		int MAX_ITER;

		
		virtual void checkStability();

	public:
		AbstractSolver(modelType* _model);
		virtual ~AbstractSolver();
		
		virtual void fill();
		virtual void start();
	
};

#endif /* ABSTRACTSOLVER_HPP_ */
