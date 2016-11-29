#include "snapshotter/Snapshotter.h"
#include "util/utils.h"

#include "model/GasOil_RZ/GasOil_RZ.h"

#include "model/Acid/2d/Acid2d.hpp"
#include "model/VPP2d/VPP2d.hpp"

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

	P_dim = model->P_dim;
	t_dim = model->t_dim;
	r_dim = model->R_dim;
}

template <>
void Snapshotter<gasOil_rz::GasOil_RZ>::setModel(gasOil_rz::GasOil_RZ* _model)
{
	model = _model;

	t_dim = model->t_dim;
	r_dim = model->R_dim;
	T_dim = model->T_dim;
	P_dim = model->P_dim;

	nx = model->cellsNum_r + 2;
	ny = model->cellsNum_z + 2;
}

template <>
void Snapshotter<acid2d::Acid2d>::setModel(acid2d::Acid2d* _model)
{
	model = _model;

	t_dim = model->t_dim;
	r_dim = model->R_dim;
	T_dim = model->T_dim;
	P_dim = model->P_dim;

	nx = model->cellsNum_r + 2;
	ny = model->cellsNum_z + 2;
}

template <>
void Snapshotter<vpp2d::VPP2d>::setModel(vpp2d::VPP2d* _model)
{
	model = _model;

	t_dim = model->t_dim;
	r_dim = model->R_dim;
	T_dim = model->T_dim;
	P_dim = model->P_dim;

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

template class Snapshotter<gasOil_rz::GasOil_RZ>;

template class Snapshotter<acid2d::Acid2d>;
template class Snapshotter<vpp2d::VPP2d>;