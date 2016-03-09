#include "model/AbstractSolver.hpp"
#include "util/utils.h"

#include "model/Oil1D/Oil1D.h"
#include "model/Gas1D/Gas1D.h"
#include "model/Gas1D/Gas1D_simple.h"
#include "model/Oil1D_NIT/Oil1D_NIT.h"
#include "model/Oil_RZ/Oil_RZ.h"
#include "model/GasOil_RZ/GasOil_RZ.h"
#include "model/GasOil_RZ_NIT/GasOil_RZ_NIT.h"

#include "model/3D/GasOil_3D/GasOil_3D.h"
#include "model/3D/GasOil_3D_NIT/GasOil_3D_NIT.h"

#include "model/3D/Perforation/GasOil_Perf.h"
#include "model/3D/Perforation/GasOil_Perf_NIT.h"

using namespace std;

template <class modelType>
AbstractSolver<modelType>::AbstractSolver(modelType* _model) : model(_model), size(_model->getCellsNum()), Tt(model->period[model->period.size()-1])
{
	cur_t = cur_t_log = 0.0;
	curTimePeriod = 0;

	idx1 = int( (model->perfIntervals[0].first + model->perfIntervals[0].second) / 2 );
	//idx2 = idx1 + model->cellsNum_z + 1;

	t_dim = model->t_dim;
}

template <>
AbstractSolver<gasOil_perf::GasOil_Perf>::AbstractSolver(gasOil_perf::GasOil_Perf* _model) : model(_model), size(_model->getCellsNum()), Tt(model->period[model->period.size() - 1])
{
	cur_t = cur_t_log = 0.0;
	curTimePeriod = 0;

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

template <>
void AbstractSolver<gasOil_perf::GasOil_Perf>::copyIterLayer()
{
	for (int i = 0; i < model->cells.size(); i++)
		model->cells[i].u_iter = model->cells[i].u_next;

	for (int i = 0; i < model->tunnelCells.size(); i++)
		model->tunnelCells[i].u_iter = model->tunnelCells[i].u_next;
}

template <class modelType>
void AbstractSolver<modelType>::copyTimeLayer()
{
	for(int i = 0; i < model->cells.size(); i++)
		model->cells[i].u_prev = model->cells[i].u_iter = model->cells[i].u_next;
}

template <>
void AbstractSolver<gasOil_perf::GasOil_Perf>::copyTimeLayer()
{
	for (int i = 0; i < model->cells.size(); i++)
		model->cells[i].u_prev = model->cells[i].u_iter = model->cells[i].u_next;

	for (int i = 0; i < model->tunnelCells.size(); i++)
		model->tunnelCells[i].u_prev = model->tunnelCells[i].u_iter = model->tunnelCells[i].u_next;
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
double AbstractSolver<gasOil_perf::GasOil_Perf>::convergance(int& ind, int& varInd)
{
	double relErr = 0.0;
	double cur_relErr = 0.0;

	double var_next, var_iter;

	for (int i = 0; i < model->cells[0].varNum; i++)
	{
		for (int j = 0; j < model->cells.size(); j++)
		{
			var_next = model->cells[j].u_next.values[i];	var_iter = model->cells[j].u_iter.values[i];
			if (fabs(var_next) > EQUALITY_TOLERANCE)
			{
				cur_relErr = fabs((var_next - var_iter) / var_next);
				if (cur_relErr > relErr)
				{
					relErr = cur_relErr;
					ind = j;
					varInd = i;
				}
			}
		}

		for (int j = 0; j < model->tunnelCells.size(); j++)
		{
			var_next = model->tunnelCells[j].u_next.values[i];	var_iter = model->tunnelCells[j].u_iter.values[i];
			if (fabs(var_next) > EQUALITY_TOLERANCE)
			{
				cur_relErr = fabs((var_next - var_iter) / var_next);
				if (cur_relErr > relErr)
				{
					relErr = cur_relErr;
					ind = j;
					varInd = i;
				}
			}
		}
	}

	return relErr;
}

template <>
double AbstractSolver<gasOil_rz_NIT::GasOil_RZ_NIT>::convergance(int& ind, int& varInd)
{
	double relErr = 0.0;
	double cur_relErr = 0.0;

	double var_next, var_iter;
		
	for(int i = 1; i < model->cells[0].varNum; i++)
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
double AbstractSolver<gasOil_3d_NIT::GasOil_3D_NIT>::convergance(int& ind, int& varInd)
{
	double relErr = 0.0;
	double cur_relErr = 0.0;

	double var_next, var_iter;
		
	for(int i = 1; i < model->cells[0].varNum; i++)
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

template <class modelType>
double AbstractSolver<modelType>::averValue(int varInd)
{
	double tmp = 0.0;

	for(int i = 0; i < model->cells.size(); i++)
	{
		tmp += model->cells[i].u_next.values[varInd] * model->cells[i].V;
	}
	
	return tmp / model->Volume;
}

template <>
double AbstractSolver<gasOil_perf::GasOil_Perf>::averValue(int varInd)
{
	double tmp = 0.0;

	for (int i = 0; i < model->cells.size(); i++)
	{
		tmp += model->cells[i].u_next.values[varInd] * model->cells[i].V;
	}

	for (int i = 0; i < model->tunnelCells.size(); i++)
	{
		tmp += model->tunnelCells[i].u_next.values[varInd] * model->tunnelCells[i].V;
	}

	return tmp / model->Volume;
}

template <class modelType>
void AbstractSolver<modelType>::construct_solution()
{
}

template class AbstractSolver<oil1D::Oil1D>;
template class AbstractSolver<gas1D::Gas1D>;
template class AbstractSolver<gas1D::Gas1D_simple>;
template class AbstractSolver<oil1D_NIT::Oil1D_NIT>;
template class AbstractSolver<oil_rz::Oil_RZ>;
template class AbstractSolver<gasOil_rz_NIT::GasOil_RZ_NIT>;
template class AbstractSolver<gasOil_rz::GasOil_RZ>;

template class AbstractSolver<gasOil_3d::GasOil_3D>;
template class AbstractSolver<gasOil_3d_NIT::GasOil_3D_NIT>;

template class AbstractSolver<gasOil_perf::GasOil_Perf>;
template class AbstractSolver<gasOil_perf_nit::GasOil_Perf_NIT>;