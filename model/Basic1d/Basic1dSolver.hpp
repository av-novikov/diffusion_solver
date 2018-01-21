#ifndef BASIC1DSOLVER_HPP_
#define BASIC1DSOLVER_HPP_

#include <iostream>
#include <vector>
#include <map>

#include "model/AbstractSolver.hpp"
#include "method/sweep.h"
#include "model/Basic1d/Basic1d.hpp"

namespace basic1d
{
	template <class modelType>
	class Basic1dSolver : public AbstractSolver<modelType>, public Sweep
	{
	protected:
		void control()
		{
			writeData();

			if (cur_t >= model->period[curTimePeriod])
			{
				curTimePeriod++;
				model->ht = model->ht_min;
				model->setPeriod(curTimePeriod);
			}

			if (model->ht <= model->ht_max && iterations < 6)
				model->ht = model->ht * 1.5;
			else if (iterations > 6 && model->ht > model->ht_min)
				model->ht = model->ht / 1.5;

			if (cur_t + model->ht > model->period[curTimePeriod])
				model->ht = model->period[curTimePeriod] - cur_t;

			cur_t += model->ht;
		}
		void doNextStep()
		{
			solveStep();
		}
		bool isChange;
	public:
		Basic1dSolver(modelType* _model) : AbstractSolver<modelType>(_model) 
		{
			Initialize(model->cellsNum_x + 2, model->var_size);
			Tt = model->period[model->period.size() - 1];
		};
		~Basic1dSolver() {};
	};

};

#endif /* BASIC1DSOLVER_HPP_ */
