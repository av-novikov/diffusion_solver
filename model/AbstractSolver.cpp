#include "model/AbstractSolver.hpp"
#include "util/utils.h"

#include "model/GasOil_RZ/GasOil_RZ.h"
#include "model/VPP2d/VPP2d.hpp"
#include "model/Bingham1d/Bingham1d.hpp"
#include "model/GasOil_Elliptic/GasOil_Elliptic.hpp"
#include "model/OilNIT_Elliptic/OilNIT_Elliptic.hpp"
#include "model/BlackOilNIT_Elliptic/BlackOilNIT_Elliptic.hpp"
#include "model/BlackOil_RZ/BlackOil_RZ.hpp"
#include "model/Acid/1d/Acid1d.hpp"
#include "model/Acid/2d/Acid2d.hpp"
#include "model/Acid/2dnit/Acid2dNIT.hpp"
#include "model/Acid/frac/AcidFracModel.hpp"
#include "model/WaxNIT/WaxNIT.hpp"

#include <iomanip>

using namespace std;

template <class modelType>
AbstractSolver<modelType>::AbstractSolver(modelType* _model) : model(_model), size(_model->getCellsNum()), Tt(model->period[model->period.size() - 1])
{
	NEWTON_STEP = 1.0;
	cur_t = cur_t_log = 0.0;
	curTimePeriod = 0;

	//idx1 = int((model->perfIntervals[0].first + model->perfIntervals[0].second) / 2);
	//idx2 = idx1 + model->cellsNum_z + 1;

	t_dim = model->t_dim;
}
AbstractSolver<gasOil_elliptic::GasOil_Elliptic>::AbstractSolver(gasOil_elliptic::GasOil_Elliptic* _model) : model(_model), size(_model->getCellsNum()), Tt(model->period[model->period.size() - 1])
{
	NEWTON_STEP = 1.0;
	cur_t = cur_t_log = 0.0;
	curTimePeriod = 0;

	t_dim = model->t_dim;
}
AbstractSolver<oilnit_elliptic::OilNIT_Elliptic>::AbstractSolver(oilnit_elliptic::OilNIT_Elliptic* _model) : model(_model), size(_model->getCellsNum()), Tt(model->period[model->period.size() - 1])
{
	NEWTON_STEP = 1.0;
	cur_t = cur_t_log = 0.0;
	curTimePeriod = 0;

	t_dim = model->t_dim;
}
AbstractSolver<blackoilnit_elliptic::BlackOilNIT_Elliptic>::AbstractSolver(blackoilnit_elliptic::BlackOilNIT_Elliptic* _model) : model(_model), size(_model->getCellsNum()), Tt(model->period[model->period.size() - 1])
{
	NEWTON_STEP = 1.0;
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
	for (auto& cell : model->cells)
		cell.u_iter = cell.u_next;
}
void AbstractSolver<gasOil_elliptic::GasOil_Elliptic>::copyIterLayer()
{
	for (auto& cell : model->cells)
		cell.u_iter = cell.u_next;

	for (auto& cell : model->wellCells)
		cell.u_iter = cell.u_next;
}
void AbstractSolver<oilnit_elliptic::OilNIT_Elliptic>::copyIterLayer()
{
	for (auto& cell : model->cells)
		cell.u_iter = cell.u_next;

	for (auto& cell : model->wellCells)
		cell.u_iter = cell.u_next;
}
void AbstractSolver<blackoilnit_elliptic::BlackOilNIT_Elliptic>::copyIterLayer()
{
	for (auto& cell : model->cells)
		cell.u_iter = cell.u_next;

	for (auto& cell : model->wellCells)
		cell.u_iter = cell.u_next;
}
void AbstractSolver<acidfrac::AcidFrac>::copyIterLayer()
{
	for (auto& cell : model->cells_frac)
		cell.u_iter = cell.u_next;

	for(auto& grid : model->poro_grids)
		for (auto& cell : grid.cells)
			cell.u_iter = cell.u_next;
}

template <class modelType>
void AbstractSolver<modelType>::revertIterLayer()
{
	for (int i = 0; i < model->cells.size(); i++)
		model->cells[i].u_next = model->cells[i].u_iter;
}
void AbstractSolver<acidfrac::AcidFrac>::revertIterLayer()
{
	for (auto& cell : model->cells_frac)
		cell.u_next = cell.u_iter;

	for (auto& grid : model->poro_grids)
		for (auto& cell : grid.cells)
			cell.u_next = cell.u_iter;
}
template <class modelType>
void AbstractSolver<modelType>::copyTimeLayer()
{
	for (auto& cell : model->cells)
		cell.u_prev = cell.u_iter = cell.u_next;
}
void AbstractSolver<gasOil_elliptic::GasOil_Elliptic>::copyTimeLayer()
{
	for (auto& cell : model->cells)
		cell.u_prev = cell.u_iter = cell.u_next;

	for (auto& cell : model->wellCells)
		cell.u_prev = cell.u_iter = cell.u_next;
}
void AbstractSolver<oilnit_elliptic::OilNIT_Elliptic>::copyTimeLayer()
{
	for (auto& cell : model->cells)
		cell.u_prev = cell.u_iter = cell.u_next;

	for (auto& cell : model->wellCells)
		cell.u_prev = cell.u_iter = cell.u_next;
}
void AbstractSolver<blackoilnit_elliptic::BlackOilNIT_Elliptic>::copyTimeLayer()
{
	for (auto& cell : model->cells)
		cell.u_prev = cell.u_iter = cell.u_next;

	for (auto& cell : model->wellCells)
		cell.u_prev = cell.u_iter = cell.u_next;
}
void AbstractSolver<acidfrac::AcidFrac>::copyTimeLayer()
{
	for (auto& cell : model->cells_frac)
		cell.u_prev = cell.u_iter = cell.u_next;

	for (auto& grid : model->poro_grids)
		for (auto& cell : grid.cells)
			cell.u_prev = cell.u_iter = cell.u_next;
}

template <class modelType>
double AbstractSolver<modelType>::convergance(int& ind, int& varInd)
{
	double relErr = 0.0;
	double cur_relErr = 0.0;

	double var_next, var_iter;
		
	for(int i = 0; i < modelType::var_size; i++)
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
double AbstractSolver<blackoil_rz::BlackOil_RZ>::convergance(int& ind, int& varInd)
{
	double relErr = 0.0;
	double cur_relErr = 0.0;
	double var_next, var_iter;

	for (int j = 0; j < model->cells.size(); j++)
	{
		blackoil_rz::Cell& cell = model->cells[j];

		for (int i = 0; i < model->var_size; i++)
		{
			if (i == 2 && !cell.u_next.SATUR) { var_next = cell.u_next.values[i+1];	var_iter = cell.u_iter.values[i+1]; }
			else { var_next = cell.u_next.values[i];	var_iter = cell.u_iter.values[i]; }

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
			else
				exit(-1);
		}
	}

	return relErr;
}
double AbstractSolver<wax_nit::WaxNIT>::convergance(int& ind, int& varInd)
{
	double relErr = 0.0;
	double cur_relErr = 0.0;
	double var_next, var_iter;

	for (int j = 0; j < model->cells.size(); j++)
	{
		wax_nit::Cell& cell = model->cells[j];

		for (int i = 0; i < model->var_size; i++)
		{
			if ((i == model->var_size - 2) && !cell.u_next.satur_gas) { var_next = cell.u_next.values[i + 2];	var_iter = cell.u_iter.values[i + 2]; }
			else if ((i == model->var_size - 1) && !cell.u_next.satur_wax) { var_next = cell.u_next.values[i + 2];	var_iter = cell.u_iter.values[i + 2]; }
			else { var_next = cell.u_next.values[i];	var_iter = cell.u_iter.values[i]; }

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
		}
	}

	return relErr;
}
double AbstractSolver<acid2d::Acid2d>::convergance(int& ind, int& varInd)
{
	double relErr = 0.0;
	double cur_relErr = 0.0;
	double var_next, var_iter;

	for (int j = 0; j < model->cells.size(); j++)
	{
		acid2d::Cell& cell = model->cells[j];

		for (int i = 0; i < model->var_size; i++)
		{
			if ((i == acid2d::Acid2d::var_size - 1) && !cell.u_next.SATUR) 
			{ var_next = cell.u_next.values[i + 1];	var_iter = cell.u_iter.values[i + 1]; }
			else 
			{ var_next = cell.u_next.values[i];	var_iter = cell.u_iter.values[i]; }

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
		}
	}

	return relErr;
}
double AbstractSolver<gasOil_elliptic::GasOil_Elliptic>::convergance(int& ind, int& varInd)
{
	double relErr = 0.0;
	double cur_relErr = 0.0;

	double var_next, var_iter;

	for (int j = 0; j < model->cells.size(); j++)
	{
		const auto& cell = model->cells[j];
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

	for (int j = 0; j < model->wellCells.size(); j++)
	{
		const auto& cell = model->wellCells[j];
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
double AbstractSolver<oilnit_elliptic::OilNIT_Elliptic>::convergance(int& ind, int& varInd)
{
	double relErr = 0.0;
	double cur_relErr = 0.0;

	double var_next, var_iter;

	for (int j = 0; j < model->cells.size(); j++)
	{
		const auto& cell = model->cells[j];
		var_next = cell.u_next.p;	var_iter = cell.u_iter.p;
		if (fabs(var_next) > EQUALITY_TOLERANCE)
		{
			cur_relErr = fabs((var_next - var_iter) / var_next);
			if (cur_relErr > relErr)
			{
				relErr = cur_relErr;
				ind = j;
				varInd = PRES;
			}
		}
	}

	for (int j = 0; j < model->wellCells.size(); j++)
	{
		const auto& cell = model->wellCells[j];
		var_next = cell.u_next.p;	var_iter = cell.u_iter.p;
		if (fabs(var_next) > EQUALITY_TOLERANCE)
		{
			cur_relErr = fabs((var_next - var_iter) / var_next);
			if (cur_relErr > relErr)
			{
				relErr = cur_relErr;
				ind = j;
				varInd = PRES;
			}
		}
	}

	return relErr;
}
double AbstractSolver<blackoilnit_elliptic::BlackOilNIT_Elliptic>::convergance(int& ind, int& varInd)
{
	double relErr = 0.0;
	double cur_relErr = 0.0;

	double var_next, var_iter;

	for (int j = 0; j < model->cells.size(); j++)
	{
		const auto& cell = model->cells[j];
		var_next = cell.u_next.p;	var_iter = cell.u_iter.p;
		if (fabs(var_next) > EQUALITY_TOLERANCE)
		{
			cur_relErr = fabs((var_next - var_iter) / var_next);
			if (cur_relErr > relErr)
			{
				relErr = cur_relErr;
				ind = j;
				varInd = PRES;
			}
		}
	}

	for (int j = 0; j < model->wellCells.size(); j++)
	{
		const auto& cell = model->wellCells[j];
		var_next = cell.u_next.p;	var_iter = cell.u_iter.p;
		if (fabs(var_next) > EQUALITY_TOLERANCE)
		{
			cur_relErr = fabs((var_next - var_iter) / var_next);
			if (cur_relErr > relErr)
			{
				relErr = cur_relErr;
				ind = j;
				varInd = PRES;
			}
		}
	}

	return relErr;
}
double AbstractSolver<acidfrac::AcidFrac>::convergance(int& ind, int& varInd)
{
	double relErr = 0.0;
	double cur_relErr = 0.0;

	double var_next, var_iter;

	for (int i = 0; i < acidfrac::var_frac_size; i++)
	{
		for (const auto& cell : model->cells_frac)
		{
			var_next = cell.u_next.values[i];	var_iter = cell.u_iter.values[i];
			if (fabs(var_next) > EQUALITY_TOLERANCE)
			{
				cur_relErr = fabs((var_next - var_iter) / var_next);
				if (cur_relErr > relErr)
				{
					relErr = cur_relErr;
					ind = cell.num;
					varInd = i;
				}
			}
		}
	}

	for (int i = 0; i < acidfrac::var_poro_size; i++)
	{
		for (const auto& grid : model->poro_grids)
		{
			for (const auto& cell : grid.cells)
			{
				var_next = cell.u_next.values[i];	var_iter = cell.u_iter.values[i];
				if (fabs(var_next) > EQUALITY_TOLERANCE)
				{
					cur_relErr = fabs((var_next - var_iter) / var_next);
					if (cur_relErr > relErr)
					{
						relErr = cur_relErr;
						ind = model->cellsNum + grid.start_idx + cell.num;
						varInd = i;
					}
				}
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
double AbstractSolver<acidfrac::AcidFrac>::averValue(const int varInd)
{
	return 0.0;
}
template <class modelType>
void AbstractSolver<modelType>::averValue(std::array<double, modelType::var_size>& aver)
{
	std::fill(aver.begin(), aver.end(), 0.0);

	for (const auto& cell : model->cells)
		for (int i = 0; i < modelType::var_size; i++)
			aver[i] += cell.u_next.values[i] * cell.V;

	for(auto& val : aver)
		val /= model->Volume;
}
void AbstractSolver<wax_nit::WaxNIT>::averValue(std::array<double, wax_nit::WaxNIT::var_size>& aver)
{
	std::fill(aver.begin(), aver.end(), 0.0);

	for (const auto& cell : model->cells)
		for (int i = 0; i < wax_nit::WaxNIT::var_size; i++)
		{
			if (i == model->var_size - 2)
			{
				if (!cell.u_next.satur_gas)
					aver[i] += cell.u_next.values[i + 2] * cell.V;
			}
			else if (i == model->var_size - 1)
			{
				if (!cell.u_next.satur_wax)
					aver[i] += cell.u_next.values[i + 2] * cell.V;
			}
			else
				aver[i] += cell.u_next.values[i] * cell.V;
		}
	for (auto& val : aver)
		val /= model->Volume;
}
void AbstractSolver<acid2d::Acid2d>::averValue(std::array<double, acid2d::Acid2d::var_size>& aver)
{
	std::fill(aver.begin(), aver.end(), 0.0);

	for (const auto& cell : model->cells)
	{
		for (int i = 0; i < acid2d::Acid2d::var_size - 1; i++)
		aver[i] += cell.u_next.values[i] * cell.V;

	if (cell.u_next.SATUR)
		aver[acid2d::Acid2d::var_size - 1] += cell.u_next.values[acid2d::Acid2d::var_size - 1] * cell.V;
	else
		aver[acid2d::Acid2d::var_size - 1] += cell.u_next.values[acid2d::Acid2d::var_size] * cell.V;
	}

	for (auto& val : aver)
		val /= model->Volume;
}
void AbstractSolver<acidfrac::AcidFrac>::averValue(std::array<double, acidfrac::AcidFrac::var_size>& aver)
{
	std::fill(aver.begin(), aver.end(), 0.0);

	for (const auto& cell : model->cells_frac)
	{
		aver[1] += cell.u_next.values[0] * cell.V;
		aver[4] += cell.u_next.values[1] * cell.V;
	}
	aver[1] /= model->Volume;	aver[4] /= model->Volume;

	for (const auto& grid : model->poro_grids)
		for (const auto& cell : grid.cells)
			for (int i = 0; i < acidfrac::var_poro_size; i++)
				aver[i] += cell.u_next.values[i] * cell.V / grid.Volume;
}
void AbstractSolver<blackoilnit_elliptic::BlackOilNIT_Elliptic>::averValue(std::array<double, blackoilnit_elliptic::BlackOilNIT_Elliptic::var_size>& aver)
{
	std::fill(aver.begin(), aver.end(), 0.0);

	for (const auto& cell : model->cells)
	{
		for (int i = 0; i < blackoilnit_elliptic::BlackOilNIT_Elliptic::var_size - 1; i++)
			aver[i] += cell.u_next.values[i] * cell.V;

		if (cell.u_next.SATUR)
			aver[blackoilnit_elliptic::BlackOilNIT_Elliptic::var_size - 1] += cell.u_next.values[blackoilnit_elliptic::BlackOilNIT_Elliptic::var_size - 1] * cell.V;
		else
			aver[blackoilnit_elliptic::BlackOilNIT_Elliptic::var_size - 1] += cell.u_next.values[blackoilnit_elliptic::BlackOilNIT_Elliptic::var_size] * cell.V;
	}

	for (auto& val : aver)
		val /= model->Volume;
}
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
double AbstractSolver<blackoil_rz::BlackOil_RZ>::averValue(const int varInd)
{
	double tmp = 0.0;

	if (varInd == 2)
		for (const auto& cell : model->cells)
		{
			if (cell.u_next.SATUR)
				tmp += cell.u_next.values[varInd] * cell.V;
			else
				tmp += cell.u_next.values[varInd + 1] * cell.V;
		}
	else
		for (const auto& cell : model->cells)
			tmp += cell.u_next.values[varInd] * cell.V;

	return tmp / model->Volume;
}
double AbstractSolver<gasOil_elliptic::GasOil_Elliptic>::averValue(const int varInd)
{
	double tmp = 0.0;

	if (varInd == 0)
	{
		for (const auto& cell : model->cells)
			tmp += cell.u_next.values[varInd] * cell.V;

		for (const auto& cell : model->wellCells)
			tmp += cell.u_next.values[varInd] * cell.V;
	}
	else
	{
		for (const auto& cell : model->cells)
		{
			if (cell.u_next.SATUR)
				tmp += cell.u_next.values[varInd] * cell.V;
			else
				tmp += cell.u_next.values[varInd + 1] * cell.V;
		}

		for (const auto& cell : model->wellCells)
		{
			if (cell.u_next.SATUR)
				tmp += cell.u_next.values[varInd] * cell.V;
			else
				tmp += cell.u_next.values[varInd + 1] * cell.V;
		}
	}

	return tmp / model->Volume;
}
double AbstractSolver<oilnit_elliptic::OilNIT_Elliptic>::averValue(const int varInd)
{
	double tmp = 0.0;

	for (const auto& cell : model->cells)
		tmp += cell.u_next.values[varInd] * cell.V;

	for (const auto& cell : model->wellCells)
		tmp += cell.u_next.values[varInd] * cell.V;

	return tmp / model->Volume;
}

template <class modelType>
void AbstractSolver<modelType>::checkStability()
{
}

template class AbstractSolver<gasOil_rz::GasOil_RZ>;
template class AbstractSolver<vpp2d::VPP2d>;
template class AbstractSolver<bing1d::Bingham1d>;
template class AbstractSolver<gasOil_elliptic::GasOil_Elliptic>;
template class AbstractSolver<oilnit_elliptic::OilNIT_Elliptic>;
template class AbstractSolver<blackoilnit_elliptic::BlackOilNIT_Elliptic>;
template class AbstractSolver<blackoil_rz::BlackOil_RZ>;
template class AbstractSolver<acid2d::Acid2d>;
template class AbstractSolver<acid2dnit::Acid2dNIT>;
template class AbstractSolver<acid1d::Acid1d>;
template class AbstractSolver<acidfrac::AcidFrac>;
template class AbstractSolver<wax_nit::WaxNIT>;