#include "model/cells/stencils/Stencil.h"

#include "model/3D/GasOil_3D/GasOil_3D.h"
#include "model/3D/Perforation/GasOil_Perf.h"

using namespace std;

/*--------------------MidStencil--------------------*/

template <class modelType>
MidStencil<modelType>::MidStencil(modelType* _model, FillFoo _foo) : model(_model), foo(_foo)
{
}

template <class modelType>
MidStencil<modelType>::~MidStencil()
{
}

template <class modelType>
void MidStencil<modelType>::setStorage(double* _matrix, int* _ind_i, int* _ind_j, double* _rhs)
{
	matrix = _matrix;
	ind_i = _ind_i;
	ind_j = _ind_j;
	rhs = _rhs;
}

template <class modelType>
void MidStencil<modelType>::fill(int cellIdx, int *counter)
{
}

template <class modelType>
void MidStencil<modelType>::fillIndex(int cellIdx, int *counter)
{
}

template <>
void MidStencil<gasOil_3d::GasOil_3D>::fill(int cellIdx, int *counter)
{
	int nebr [7];
	model->getStencilIdx(cellIdx, nebr);

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 7; j++)
		{
			// dP
			matrix[(*counter)++] = foo.mat[i][2 * j](cellIdx, nebr[j]);

			// dS
			matrix[(*counter)++] = foo.mat[i][2 * j + 1](cellIdx, nebr[j]);
		}

		rhs[2 * cellIdx + i] = -foo.rhs[i](cellIdx);
	}		
}

template <>
void MidStencil<gasOil_3d::GasOil_3D>::fillIndex(int cellIdx, int *counter)
{
	int nebr[7];
	model->getStencilIdx(cellIdx, nebr);

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 7; j++)
		{
			ind_i[*counter] = 2 * cellIdx + i;
			ind_j[(*counter)++] = 2 * nebr[j];
			
			ind_i[*counter] = 2 * cellIdx + i;
			ind_j[(*counter)++] = 2 * nebr[j] + 1;
		}
	}
}

template <>
void MidStencil<gasOil_perf::GasOil_Perf>::fill(int cellIdx, int *counter)
{
	gasOil_perf::Cell* nebr[7];
	model->getStencilIdx(cellIdx, nebr);

	if (nebr[0]->isUsed)
	{
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 7; j++)
			{
				matrix[(*counter)++] = foo.mat[i][2 * j](cellIdx, nebr[j]->num);
				matrix[(*counter)++] = foo.mat[i][2 * j + 1](cellIdx, nebr[j]->num);
			}

			rhs[2 * cellIdx + i] = -foo.rhs[i](cellIdx);
		}
	}
	else
	{
		matrix[(*counter)++] = 1.0;
		matrix[(*counter)++] = 1.0;

		rhs[2 * cellIdx] = 0.0;
		rhs[2 * cellIdx + 1] = 0.0;
	}
}

template <>
void MidStencil<gasOil_perf::GasOil_Perf>::fillIndex(int cellIdx, int *counter)
{
	gasOil_perf::Cell* nebr[7];
	model->getStencilIdx(cellIdx, nebr);

	if (nebr[0]->isUsed)
	{
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 7; j++)
			{
				ind_i[*counter] = 2 * cellIdx + i;
				if (nebr[j]->isUsed)
					ind_j[(*counter)++] = 2 * nebr[j]->num;
				else
					ind_j[(*counter)++] = 2 * (model->cellsNum + model->getCell(nebr[0]->num, nebr[j]->num).num);

				ind_i[*counter] = 2 * cellIdx + i;
				if (nebr[j]->isUsed)
					ind_j[(*counter)++] = 2 * nebr[j]->num + 1;
				else
					ind_j[(*counter)++] = 2 * (model->cellsNum + model->getCell(nebr[0]->num, nebr[j]->num).num) + 1;
			}

		}
	}
	else {
		ind_i[*counter] = 2 * cellIdx;
		ind_j[(*counter)++] = 2 * cellIdx;
		ind_i[*counter] = 2 * cellIdx + 1;
		ind_j[(*counter)++] = 2 * cellIdx + 1;
	}
}

/*--------------------LeftStencil--------------------*/

template <class modelType>
LeftStencil<modelType>::LeftStencil(modelType* _model, FillFoo _foo) : model(_model), foo(_foo)
{
}

template <class modelType>
LeftStencil<modelType>::~LeftStencil()
{
}

template <class modelType>
void LeftStencil<modelType>::setStorage(double* _matrix, int* _ind_i, int* _ind_j, double* _rhs)
{
	matrix = _matrix;
	ind_i = _ind_i;
	ind_j = _ind_j;
	rhs = _rhs;
}

template <class modelType>
void LeftStencil<modelType>::fill(int cellIdx, int *counter)
{
}

template <class modelType>
void LeftStencil<modelType>::fillIndex(int cellIdx, int *counter)
{
}

template <>
void LeftStencil<gasOil_3d::GasOil_3D>::fill(int cellIdx, int *counter)
{
	int nebr1Idx = cellIdx + model->cellsNum_z + 2;
	int nebr2Idx = cellIdx + 2 * model->cellsNum_z + 4;

	for (int i = 0; i < 2; i++)
	{
		matrix[(*counter)++] = foo.mat[i][0](cellIdx, cellIdx);
		matrix[(*counter)++] = foo.mat[i][1](cellIdx, cellIdx);

		matrix[(*counter)++] = foo.mat[i][2](cellIdx, nebr1Idx);
		matrix[(*counter)++] = foo.mat[i][3](cellIdx, nebr1Idx);

		if (i == 1)
		{
			matrix[(*counter)++] = foo.mat[i][4](cellIdx, nebr2Idx);
			matrix[(*counter)++] = foo.mat[i][5](cellIdx, nebr2Idx);
		}

		// Right side
		rhs[2 * cellIdx + i] = -foo.rhs[i](cellIdx);
	}
}

template <>
void LeftStencil<gasOil_3d::GasOil_3D>::fillIndex(int cellIdx, int *counter)
{
	int nebr1Idx = cellIdx + model->cellsNum_z + 2;
	int nebr2Idx = cellIdx + 2 * model->cellsNum_z + 4;

	for (int i = 0; i < 2; i++)
	{
		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * cellIdx;

		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * cellIdx + 1;

		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * nebr1Idx;

		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * nebr1Idx + 1;

		if (i == 1)
		{
			ind_i[*counter] = 2 * cellIdx + 1;
			ind_j[(*counter)++] = 2 * nebr2Idx;

			ind_i[*counter] = 2 * cellIdx + 1;
			ind_j[(*counter)++] = 2 * nebr2Idx + 1;
		}
	}
}

template <>
void LeftStencil<gasOil_perf::GasOil_Perf>::fill(int cellIdx, int *counter)
{
	int nebr1Idx = model->nebrMap[cellIdx].first;
	int nebr2Idx = model->nebrMap[cellIdx].second;

	for (int i = 0; i < 2; i++)
	{
		matrix[(*counter)++] = foo.mat[i][0](cellIdx, cellIdx);
		matrix[(*counter)++] = foo.mat[i][1](cellIdx, cellIdx);

		matrix[(*counter)++] = foo.mat[i][2](cellIdx, nebr1Idx);
		matrix[(*counter)++] = foo.mat[i][3](cellIdx, nebr1Idx);

		if (i == 1)
		{
			matrix[(*counter)++] = foo.mat[i][4](cellIdx, nebr2Idx);
			matrix[(*counter)++] = foo.mat[i][5](cellIdx, nebr2Idx);
		}

		// Right side
		rhs[2 * (cellIdx + model->cellsNum) + i] = -foo.rhs[i](cellIdx);
	}
}

template <>
void LeftStencil<gasOil_perf::GasOil_Perf>::fillIndex(int cellIdx, int *counter)
{
	const int cellMatIdx = cellIdx + model->cellsNum;
	const int nebr1Idx = model->nebrMap[cellIdx].first;
	const int nebr2Idx = model->nebrMap[cellIdx].second;
	
	for (int i = 0; i < 2; i++)
	{
		ind_i[*counter] = 2 * cellMatIdx + i;
		ind_j[(*counter)++] = 2 * cellMatIdx;

		ind_i[*counter] = 2 * cellMatIdx + i;
		ind_j[(*counter)++] = 2 * cellMatIdx + 1;

		ind_i[*counter] = 2 * cellMatIdx + i;
		ind_j[(*counter)++] = 2 * nebr1Idx;

		ind_i[*counter] = 2 * cellMatIdx + i;
		ind_j[(*counter)++] = 2 * nebr1Idx + 1;

		if (i == 1)
		{
			ind_i[*counter] = 2 * cellMatIdx + 1;
			ind_j[(*counter)++] = 2 * nebr2Idx;

			ind_i[*counter] = 2 * cellMatIdx + 1;
			ind_j[(*counter)++] = 2 * nebr2Idx + 1;
		}
	}
}

/*--------------------RightStencil--------------------*/

template <class modelType>
RightStencil<modelType>::RightStencil(modelType* _model, FillFoo _foo) : model(_model), foo(_foo)
{
}

template <class modelType>
RightStencil<modelType>::~RightStencil()
{
}

template <class modelType>
void RightStencil<modelType>::setStorage(double* _matrix, int* _ind_i, int* _ind_j, double* _rhs)
{
	matrix = _matrix;
	ind_i = _ind_i;
	ind_j = _ind_j;
	rhs = _rhs;
}

template <class modelType>
void RightStencil<modelType>::fill(int cellIdx, int *counter)
{
}

template <class modelType>
void RightStencil<modelType>::fillIndex(int cellIdx, int *counter)
{
}

template <>
void RightStencil<gasOil_3d::GasOil_3D>::fill(int cellIdx, int *counter)
{
	int nebrIdx = cellIdx - model->cellsNum_z - 2;
	for (int i = 0; i < 2; i++)
	{
		matrix[(*counter)++] = foo.mat[i][0](cellIdx, cellIdx);
		matrix[(*counter)++] = foo.mat[i][1](cellIdx, cellIdx);

		matrix[(*counter)++] = foo.mat[i][2](cellIdx, nebrIdx);
		matrix[(*counter)++] = foo.mat[i][3](cellIdx, nebrIdx);

		// Right side
		rhs[2 * cellIdx + i] = -foo.rhs[i](cellIdx);
	}
}

template <>
void RightStencil<gasOil_3d::GasOil_3D>::fillIndex(int cellIdx, int *counter)
{
	int nebrIdx = cellIdx - model->cellsNum_z - 2;
	
	for (int i = 0; i < 2; i++)
	{
		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * cellIdx;

		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * cellIdx + 1;

		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * nebrIdx;

		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * nebrIdx + 1;
	}
}

void RightStencil<gasOil_perf::GasOil_Perf>::fill(int cellIdx, int *counter)
{
	int nebrIdx = cellIdx - model->cellsNum_z - 2;
	for (int i = 0; i < 2; i++)
	{
		matrix[(*counter)++] = foo.mat[i][0](cellIdx, cellIdx);
		matrix[(*counter)++] = foo.mat[i][1](cellIdx, cellIdx);

		matrix[(*counter)++] = foo.mat[i][2](cellIdx, nebrIdx);
		matrix[(*counter)++] = foo.mat[i][3](cellIdx, nebrIdx);

		// Right side
		rhs[2 * cellIdx + i] = -foo.rhs[i](cellIdx);
	}
}

template <>
void RightStencil<gasOil_perf::GasOil_Perf>::fillIndex(int cellIdx, int *counter)
{
	int nebrIdx = cellIdx - model->cellsNum_z - 2;

	for (int i = 0; i < 2; i++)
	{
		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * cellIdx;

		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * cellIdx + 1;

		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * nebrIdx;

		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * nebrIdx + 1;
	}
}

/*--------------------TopStencil--------------------*/

template <class modelType>
TopStencil<modelType>::TopStencil(modelType* _model, FillFoo _foo) : model(_model), foo(_foo)
{
}

template <class modelType>
TopStencil<modelType>::~TopStencil()
{
}

template <class modelType>
void TopStencil<modelType>::setStorage(double* _matrix, int* _ind_i, int* _ind_j, double* _rhs)
{
	matrix = _matrix;
	ind_i = _ind_i;
	ind_j = _ind_j;
	rhs = _rhs;
}

template <class modelType>
void TopStencil<modelType>::fill(int cellIdx, int *counter)
{
}

template <class modelType>
void TopStencil<modelType>::fillIndex(int cellIdx, int *counter)
{
}

template <>
void TopStencil<gasOil_3d::GasOil_3D>::fill(int cellIdx, int *counter)
{
	int nebrIdx = cellIdx + 1;
	for (int i = 0; i < 2; i++)
	{
		matrix[(*counter)++] = foo.mat[i][0](cellIdx, cellIdx);
		matrix[(*counter)++] = foo.mat[i][1](cellIdx, cellIdx);

		matrix[(*counter)++] = foo.mat[i][2](cellIdx, nebrIdx);
		matrix[(*counter)++] = foo.mat[i][3](cellIdx, nebrIdx);

		// Right side
		rhs[2 * cellIdx + i] = -foo.rhs[i](cellIdx);
	}
}

template <>
void TopStencil<gasOil_3d::GasOil_3D>::fillIndex(int cellIdx, int *counter)
{
	int nebrIdx = cellIdx + 1;
	for (int i = 0; i < 2; i++)
	{
		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * cellIdx;

		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * cellIdx + 1;

		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * nebrIdx;

		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * nebrIdx + 1;
	}
}

template <>
void TopStencil<gasOil_perf::GasOil_Perf>::fill(int cellIdx, int *counter)
{
	int nebrIdx = cellIdx + 1;
	for (int i = 0; i < 2; i++)
	{
		matrix[(*counter)++] = foo.mat[i][0](cellIdx, cellIdx);
		matrix[(*counter)++] = foo.mat[i][1](cellIdx, cellIdx);

		matrix[(*counter)++] = foo.mat[i][2](cellIdx, nebrIdx);
		matrix[(*counter)++] = foo.mat[i][3](cellIdx, nebrIdx);

		// Right side
		rhs[2 * cellIdx + i] = -foo.rhs[i](cellIdx);
	}
}

template <>
void TopStencil<gasOil_perf::GasOil_Perf>::fillIndex(int cellIdx, int *counter)
{
	int nebrIdx = cellIdx + 1;
	for (int i = 0; i < 2; i++)
	{
		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * cellIdx;

		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * cellIdx + 1;

		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * nebrIdx;

		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * nebrIdx + 1;
	}
}

/*--------------------BotStencil--------------------*/

template <class modelType>
BotStencil<modelType>::BotStencil(modelType* _model, FillFoo _foo) : model(_model), foo(_foo)
{
}

template <class modelType>
BotStencil<modelType>::~BotStencil()
{
}

template <class modelType>
void BotStencil<modelType>::setStorage(double* _matrix, int* _ind_i, int* _ind_j, double* _rhs)
{
	matrix = _matrix;
	ind_i = _ind_i;
	ind_j = _ind_j;
	rhs = _rhs;
}

template <class modelType>
void BotStencil<modelType>::fill(int cellIdx, int *counter)
{
}

template <class modelType>
void BotStencil<modelType>::fillIndex(int cellIdx, int *counter)
{
}

template <>
void BotStencil<gasOil_3d::GasOil_3D>::fill(int cellIdx, int *counter)
{
	int nebrIdx = cellIdx - 1;
	for (int i = 0; i < 2; i++)
	{
		matrix[(*counter)++] = foo.mat[i][0](cellIdx, cellIdx);
		matrix[(*counter)++] = foo.mat[i][1](cellIdx, cellIdx);

		matrix[(*counter)++] = foo.mat[i][2](cellIdx, nebrIdx);
		matrix[(*counter)++] = foo.mat[i][3](cellIdx, nebrIdx);

		// Right side
		rhs[2 * cellIdx + i] = -foo.rhs[i](cellIdx);
	}
}

template <>
void BotStencil<gasOil_3d::GasOil_3D>::fillIndex(int cellIdx, int *counter)
{
	int nebrIdx = cellIdx - 1;
	for (int i = 0; i < 2; i++)
	{
		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * cellIdx;

		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * cellIdx + 1;

		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * nebrIdx;

		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * nebrIdx + 1;
	}
}

template <>
void BotStencil<gasOil_perf::GasOil_Perf>::fill(int cellIdx, int *counter)
{
	int nebrIdx = cellIdx - 1;
	for (int i = 0; i < 2; i++)
	{
		matrix[(*counter)++] = foo.mat[i][0](cellIdx, cellIdx);
		matrix[(*counter)++] = foo.mat[i][1](cellIdx, cellIdx);

		matrix[(*counter)++] = foo.mat[i][2](cellIdx, nebrIdx);
		matrix[(*counter)++] = foo.mat[i][3](cellIdx, nebrIdx);

		// Right side
		rhs[2 * cellIdx + i] = -foo.rhs[i](cellIdx);
	}
}

template <>
void BotStencil<gasOil_perf::GasOil_Perf>::fillIndex(int cellIdx, int *counter)
{
	int nebrIdx = cellIdx - 1;
	for (int i = 0; i < 2; i++)
	{
		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * cellIdx;

		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * cellIdx + 1;

		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * nebrIdx;

		ind_i[*counter] = 2 * cellIdx + i;
		ind_j[(*counter)++] = 2 * nebrIdx + 1;
	}
}
/*--------------------UsedStencils--------------------*/

template <class modelType>
UsedStencils<modelType>::UsedStencils(modelType* _model)
{
	middle = new MidStencil<modelType>(_model, _model->middleFoo);
	left = new LeftStencil<modelType>(_model, _model->leftFoo);
	right = new RightStencil<modelType>(_model, _model->rightFoo);
	top = new TopStencil<modelType>(_model, _model->topFoo);
	bot = new BotStencil<modelType>(_model, _model->botFoo);
}

template <class modelType>
UsedStencils<modelType>::~UsedStencils()
{
	delete middle;
	delete left;
	delete right;
	delete top;
	delete bot;
}

template <class modelType>
void UsedStencils<modelType>::setStorages(double* _matrix, int* _ind_i, int* _ind_j, double* _rhs)
{
	middle->setStorage(_matrix, _ind_i, _ind_j, _rhs);
	left->setStorage(_matrix, _ind_i, _ind_j, _rhs);
	right->setStorage(_matrix, _ind_i, _ind_j, _rhs);
	top->setStorage(_matrix, _ind_i, _ind_j, _rhs);
	bot->setStorage(_matrix, _ind_i, _ind_j, _rhs);
}

template class MidStencil<gasOil_3d::GasOil_3D>;
template class LeftStencil<gasOil_3d::GasOil_3D>;
template class RightStencil<gasOil_3d::GasOil_3D>;
template class TopStencil<gasOil_3d::GasOil_3D>;
template class BotStencil<gasOil_3d::GasOil_3D>;
template class UsedStencils<gasOil_3d::GasOil_3D>;

template class MidStencil<gasOil_perf::GasOil_Perf>;
template class LeftStencil<gasOil_perf::GasOil_Perf>;
template class RightStencil<gasOil_perf::GasOil_Perf>;
template class TopStencil<gasOil_perf::GasOil_Perf>;
template class BotStencil<gasOil_perf::GasOil_Perf>;
template class UsedStencils<gasOil_perf::GasOil_Perf>;
