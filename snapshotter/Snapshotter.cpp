#include "snapshotter/Snapshotter.h"
#include "util/utils.h"

#include "model/GasOil_RZ/GasOil_RZ.h"

#include "model/Acid/2d/Acid2d.hpp"
#include "model/VPP2d/VPP2d.hpp"
#include "model/Bingham1d/Bingham1d.hpp"
#include "model/GasOil_Elliptic/GasOil_Elliptic.hpp"
#include "model/OilNIT_Elliptic/OilNIT_Elliptic.hpp"
#include "model/GasOilNIT_Elliptic/GasOilNIT_Elliptic.hpp"

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
void Snapshotter<gasOil_elliptic::GasOil_Elliptic>::setModel(gasOil_elliptic::GasOil_Elliptic* _model)
{
	model = _model;

	t_dim = model->t_dim;
	r_dim = model->R_dim;
	T_dim = model->T_dim;
	P_dim = model->P_dim;

	nx = model->cellsNum_mu + 2;
	ny = model->cellsNum_nu;
	nz = model->cellsNum_z + 2;
}
template <>
void Snapshotter<gasOilnit_elliptic::GasOilNIT_Elliptic>::setModel(gasOilnit_elliptic::GasOilNIT_Elliptic* _model)
{
	model = _model;

	T_dim = model->T_dim;
	t_dim = model->t_dim;
	r_dim = model->R_dim;
	T_dim = model->T_dim;
	P_dim = model->P_dim;

	nx = model->cellsNum_mu + 2;
	ny = model->cellsNum_nu;
	nz = model->cellsNum_z + 2;
}
template <>
void Snapshotter<oilnit_elliptic::OilNIT_Elliptic>::setModel(oilnit_elliptic::OilNIT_Elliptic* _model)
{
	model = _model;

	T_dim = model->T_dim;
	t_dim = model->t_dim;
	r_dim = model->R_dim;
	T_dim = model->T_dim;
	P_dim = model->P_dim;

	nx = model->cellsNum_mu + 2;
	ny = model->cellsNum_nu;
	nz = model->cellsNum_z + 2;
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

template <>
void Snapshotter<bing1d::Bingham1d>::setModel(bing1d::Bingham1d* _model)
{
	model = _model;

	t_dim = model->t_dim;
	r_dim = model->R_dim;
	P_dim = model->P_dim;

	nx = model->cellsNum_r + 2;
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
template class Snapshotter<bing1d::Bingham1d>;
template class Snapshotter<gasOil_elliptic::GasOil_Elliptic>;
template class Snapshotter<oilnit_elliptic::OilNIT_Elliptic>;
template class Snapshotter<gasOilnit_elliptic::GasOilNIT_Elliptic>;