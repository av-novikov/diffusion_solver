#include "model/AcidFrac/FractureModel.hpp"

using namespace frac;

FractureModel::FractureModel()
{
	x = new double[stencil * Variable::size];
	y = new double[var_size];

	jac = new double*[var_size];
	for (int i = 0; i < var_size; i++)
		jac[i] = new double[stencil * Variable::size];
}
FractureModel::~FractureModel()
{
	delete[] x, y, h;
	delete[] var;

	for (int i = 0; i < var_size; i++)
		delete[] jac[i];
	delete[] jac;
}
void FractureModel::setProps(const Properties& props)
{
	cellsNum_x = props.cellsNum_x;
	cellsNum_y = props.cellsNum_y;
	cellsNum_z = props.cellsNum_z;
	w2 = props.w2;	
	l2 = props.l2;	
	height = props.height;
}
void FractureModel::makeDimLess()
{
	R_dim = l2;
	P_dim = p_init;

	l2 /= R_dim;	w2 /= R_dim;	height /= R_dim;
}
void FractureModel::setInitialState()
{
	for(auto& cell : cells)
	{
		cell.u_prev.p = cell.u_iter.p = cell.u_next.p = p_init;
		cell.u_prev.c = cell.u_iter.c = cell.u_next.c = c_init;
	}
	var = new TapeVariable[stencil];
	h = new adouble[var_size];
}
void FractureModel::buildGrid()
{
	int counter = 0;
	const double hx = l2 / cellsNum_x;
	const double hy = w2 / cellsNum_y;
	const double hz = height / cellsNum_z;
	double x = 0.0, y = 0.0, z = 0.0;
	Type cur_type;
	for (int i = 0; i < cellsNum_x + 2; i++)
	{
		for (int j = 0; j < cellsNum_y + 1; j++)
		{
			for (int k = 0; k < cellsNum_z + 2; k++)
			{

				cells.push_back(Cell(counter++, x, y, z, hx, hy, hz, cur_type));
				if (k == 0 || k == cellsNum_z)
					z += hz / 2.0;
				else
					z += hz;
			}
			if (j == 0)
				y += hy / 2.0;
			else
				y += hy;
		}
		if (i == 0 || i == cellsNum_x)
			x += hx / 2.0;
		else
			x += hx;
	}
}
void FractureModel::setPeriod(int period)
{
}