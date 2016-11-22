#include "model/AbstractSolver.hpp"
#include "util/utils.h"

#include "model/GasOil_RZ/GasOil_RZ.h"

#include <iomanip>

using namespace std;

template <class modelType>
AbstractSolver<modelType>::AbstractSolver(modelType* _model) : model(_model), size(_model->getCellsNum()), Tt(model->period[model->period.size()-1])
{
	newton_step = 1.0;
	isWellboreAffect = false;
	cur_t = cur_t_log = 0.0;
	curTimePeriod = 0;

	idx1 = int( (model->perfIntervals[0].first + model->perfIntervals[0].second) / 2 );
	//idx2 = idx1 + model->cellsNum_z + 1;

	t_dim = model->t_dim;
}


template <class modelType>
AbstractSolver<modelType>::~AbstractSolver()
{
}

template <class modelType>
void AbstractSolver<modelType>::start()
{
	int counter = 0;
	iterations = 8;

	model->setPeriod(curTimePeriod);
	while(cur_t < Tt)
	{
		control();
		if( model->isWriteSnaps )
			model->snapshot_all(counter++);
		doNextStep();
		copyTimeLayer();
		cout << "---------------------NEW TIME STEP---------------------" << endl;
		cout << setprecision(6);
		cout << "time = " << cur_t << endl;
	}
	if( model->isWriteSnaps )
		model->snapshot_all(counter++);
	writeData();
}

template <class modelType>
void AbstractSolver<modelType>::fill()
{
}

template <class modelType>
void AbstractSolver<modelType>::copyIterLayer()
{
	for (int i = 0; i < model->cells.size(); i++)
		model->cells[i].u_iter = model->cells[i].u_next;
}

template <class modelType>
void AbstractSolver<modelType>::revertIterLayer()
{
	for (int i = 0; i < model->cells.size(); i++)
		model->cells[i].u_next = model->cells[i].u_iter;
}

template <class modelType>
void AbstractSolver<modelType>::copyTimeLayer()
{
	for(int i = 0; i < model->cells.size(); i++)
		model->cells[i].u_prev = model->cells[i].u_iter = model->cells[i].u_next;
}

template <class modelType>
double AbstractSolver<modelType>::convergance(int& ind, int& varInd)
{
	double relErr = 0.0;
	double cur_relErr = 0.0;

	double var_next, var_iter;
		
	for(int i = 0; i < model->cells[0].varNum; i++)
	{
		for(int j = 0; j < model->cells.size(); j++)
		{
			var_next = model->cells[j].u_next.values[i];	var_iter = model->cells[j].u_iter.values[i];
			if(fabs(var_next) > EQUALITY_TOLERANCE)	
			{
				cur_relErr = fabs( (var_next - var_iter) / var_next );
				if(cur_relErr > relErr)
				{
					relErr = cur_relErr;
					ind  = j;
					varInd = i;
				}
			}
		}
	}
	
	return relErr;
}

template <>
double AbstractSolver<gasOil_rz::GasOil_RZ>::convergance(int& ind, int& varInd)
{
	double relErr = 0.0;
	double cur_relErr = 0.0;

	double var_next, var_iter;

	for (int j = 0; j < model->cells.size(); j++)
	{
		gasOil_rz::Cell& cell = model->cells[j];
		var_next = cell.u_next.p;	var_iter = cell.u_iter.p;
		if (fabs(var_next) > EQUALITY_TOLERANCE)
		{
			cur_relErr = fabs((var_next - var_iter) / var_next);
			if (cur_relErr > relErr)
			{
				relErr = cur_relErr;
				ind = j;
				varInd = 0;
			}
		}
		
		if (cell.u_next.SATUR)
		{
			var_next = cell.u_next.s;	var_iter = cell.u_iter.s;
		}
		else {
			var_next = cell.u_next.p_bub;	var_iter = cell.u_iter.p_bub;
		}
		if (fabs(var_next) > EQUALITY_TOLERANCE)
		{
			cur_relErr = fabs((var_next - var_iter) / var_next);
			if (cur_relErr > relErr)
			{
				relErr = cur_relErr;
				ind = j;
				varInd = 1;
			}
		}
	}

	return relErr;
}

template <class modelType>
double AbstractSolver<modelType>::averValue(const int varInd)
{
	double tmp = 0.0;

	for(const auto& cell : model->cells)
	{
		tmp += cell.u_next.values[varInd] * cell.V;
	}
	
	return tmp / model->Volume;
}

template <>
double AbstractSolver<gasOil_rz::GasOil_RZ>::averValue(const int varInd)
{
	double tmp = 0.0;

	if (varInd == 0)
		for (const auto& cell : model->cells)
		{
			tmp += cell.u_next.values[varInd] * cell.V;
		}
	else
		for(const auto& cell : model->cells)
		{
			if(cell.u_next.SATUR)
				tmp += cell.u_next.values[varInd] * cell.V;
			else
				tmp += cell.u_next.values[varInd+1] * cell.V;
		}

	return tmp / model->Volume;
}

template <class modelType>
void AbstractSolver<modelType>::construct_solution()
{
}

template class AbstractSolver<gasOil_rz::GasOil_RZ>;