#include "snapshotter/Snapshotter.h"
#include "util/utils.h"

#include "model/Oil1D/Oil1D.h"
#include "model/Oil1D_NIT/Oil1D_NIT.h"
#include "model/Oil_RZ/Oil_RZ.h"
#include "model/GasOil_RZ/GasOil_RZ.h"
#include "model/GasOil_RZ_NIT/GasOil_RZ_NIT.h"

using namespace std;

template <class modelType>
const string Snapshotter<modelType>::prefix = "snaps/";

template <class modelType>
Snapshotter<modelType>::Snapshotter()
{
}

template <class modelType>
Snapshotter<modelType>::~Snapshotter()
{
}

template <class modelType>
void Snapshotter<modelType>::setModel(modelType* _model)
{
	model = _model;

	t_dim = model->t_dim;
	r_dim = model->R_dim;
}

template <>
void Snapshotter<oil1D::Oil1D>::setModel(oil1D::Oil1D* _model)
{
	model = _model;

	t_dim = model->t_dim;
	r_dim = model->R_dim;

	nx = model->cellsNum_r + 2;
}

template <>
void Snapshotter<oil_rz::Oil_RZ>::setModel(oil_rz::Oil_RZ* _model)
{
	model = _model;

	t_dim = model->t_dim;
	r_dim = model->R_dim;
	T_dim = model->T_dim;

	nx = model->cellsNum_r + 2;
	ny = model->cellsNum_z + 2;
}

template <>
void Snapshotter<gasOil_rz::GasOil_RZ>::setModel(gasOil_rz::GasOil_RZ* _model)
{
	model = _model;

	t_dim = model->t_dim;
	r_dim = model->R_dim;
	T_dim = model->T_dim;

	nx = model->cellsNum_r + 2;
	ny = model->cellsNum_z + 2;
}

template <>
void Snapshotter<gasOil_rz_NIT::GasOil_RZ_NIT>::setModel(gasOil_rz_NIT::GasOil_RZ_NIT* _model)
{
	model = _model;

	t_dim = model->t_dim;
	r_dim = model->R_dim;
	T_dim = model->T_dim;

	nx = model->cellsNum_r + 2;
	ny = model->cellsNum_z + 2;
}

template <class modelType>
string Snapshotter<modelType>::replace(string filename, string from, string to)
{
	size_t start_pos = 0;
    while((start_pos = filename.find(from, start_pos)) != string::npos) 
	{
		filename.replace(start_pos, from.length(), to);
        start_pos += to.length();
    }
	return filename;
}

template <class modelType>
string Snapshotter<modelType>::getFileName(int i)
{
	string filename = pattern;
	return replace(filename, "%{STEP}", to_string(i));
}

template class Snapshotter<oil1D::Oil1D>;
template class Snapshotter<oil1D_NIT::Oil1D_NIT>;
template class Snapshotter<oil_rz::Oil_RZ>;
template class Snapshotter<gasOil_rz::GasOil_RZ>;
template class Snapshotter<gasOil_rz_NIT::GasOil_RZ_NIT>;