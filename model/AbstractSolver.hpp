#ifndef ABSTRACTSOLVER_HPP_
#define ABSTRACTSOLVER_HPP_

#include <algorithm>
#include <cmath>
#include <iostream>

#define TEMP 0
#define PRES 1
#define SAT 2

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
		void copyTimeLayer();
		
		double convergance(int& ind, int& varInd);
		double averValue(int varInd);
		
		virtual void construct_solution();
		virtual void writeData() = 0;
		virtual void control() = 0;
		virtual void doNextStep() = 0;

	public:
		AbstractSolver(modelType* _model);
		virtual ~AbstractSolver();
		
		virtual void fill();
		virtual void start();
	
};

#endif /* ABSTRACTSOLVER_HPP_ */
