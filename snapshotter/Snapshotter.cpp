#include "snapshotter/Snapshotter.h"
#include "util/utils.h"

#include "model/GasOil_RZ/GasOil_RZ.h"

#include "model/Acid/2dnit/Acid2dNIT.hpp"
#include "model/Acid/2d/Acid2d.hpp"
#include "model/Acid/1d/Acid1d.hpp"
#include "model/Acid/frac/AcidFracModel.hpp"
#include "model/Acid/ellfrac/AcidEllFracModel.hpp"
#include "model/Acid/recfrac/AcidRecFracModel.hpp"
#include "model/Acid/recfracmov/AcidRecFracMovModel.hpp"
#include "model/VPP2d/VPP2d.hpp"
#include "model/Bingham1d/Bingham1d.hpp"
#include "model/GasOil_Elliptic/GasOil_Elliptic.hpp"
#include "model/OilNIT_Elliptic/OilNIT_Elliptic.hpp"
#include "model/BlackOilNIT_Elliptic/BlackOilNIT_Elliptic.hpp"
#include "model/BlackOil_RZ/BlackOil_RZ.hpp"
#include "model/WaxNIT/2d/WaxNIT.hpp"
#include "model/WaxNIT/1d/WaxNIT1d.hpp"

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
template<>
void Snapshotter<acidfrac::AcidFrac>::setModel(acidfrac::AcidFrac* _model)
{
	model = _model;

	T_dim = model->T_dim;
	t_dim = model->t_dim;
	r_dim = model->R_dim;
	T_dim = model->T_dim;
	P_dim = model->P_dim;

	nx = model->cellsNum_x + 2;
	ny = model->cellsNum_y + 1;
	nz = model->cellsNum_z + 2;
}
template<>
void Snapshotter<acidellfrac::AcidEllFrac>::setModel(acidellfrac::AcidEllFrac* _model)
{
	model = _model;

	T_dim = model->T_dim;
	t_dim = model->t_dim;
	r_dim = model->R_dim;
	T_dim = model->T_dim;
	P_dim = model->P_dim;

	nx = model->cellsNum_nu + 2;
	ny = model->cellsNum_mu_frac + 1;
	nz = model->cellsNum_z + 2;
}
template<>
void Snapshotter<acidrecfracmov::AcidRecFracMov>::setModel(acidrecfracmov::AcidRecFracMov* _model)
{
	model = _model;

	T_dim = model->T_dim;
	t_dim = model->t_dim;
	r_dim = model->R_dim;
	T_dim = model->T_dim;
	P_dim = model->P_dim;

	nx = model->cellsNum_x + 2;
	ny = model->cellsNum_y_frac + 1;
	nz = model->cellsNum_z + 2;
}
template <class modelType>
void Snapshotter<modelType>::setModel(modelType* _model)
{
	model = _model;

	P_dim = model->P_dim;
	t_dim = model->t_dim;
	r_dim = model->R_dim;
}
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
void Snapshotter<blackoilnit_elliptic::BlackOilNIT_Elliptic>::setModel(blackoilnit_elliptic::BlackOilNIT_Elliptic* _model)
{
	model = _model;

	T_dim = model->T_dim;
	t_dim = model->t_dim;
	r_dim = model->R_dim;
	T_dim = model->T_dim;
	P_dim = model->P_dim;

	nx = model->cellsNum_nu + 2; 
	ny = model->cellsNum_mu + 1;
	nz = model->cellsNum_z + 2;
}
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
void Snapshotter<acid1d::Acid1d>::setModel(acid1d::Acid1d* _model)
{
	model = _model;

	t_dim = model->t_dim;
	r_dim = model->R_dim;
	T_dim = model->T_dim;
	P_dim = model->P_dim;

	nx = model->cellsNum_x + 2;
	ny = 3;
}
void Snapshotter<acid2dnit::Acid2dNIT>::setModel(acid2dnit::Acid2dNIT* _model)
{
	model = _model;

	t_dim = model->t_dim;
	r_dim = model->R_dim;
	T_dim = model->T_dim;
	P_dim = model->P_dim;

	nx = model->cellsNum_r + 2;
	ny = model->cellsNum_z + 2;
}
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
void Snapshotter<bing1d::Bingham1d>::setModel(bing1d::Bingham1d* _model)
{
	model = _model;

	t_dim = model->t_dim;
	r_dim = model->R_dim;
	P_dim = model->P_dim;

	nx = model->cellsNum_r + 2;
}
void Snapshotter<blackoil_rz::BlackOil_RZ>::setModel(blackoil_rz::BlackOil_RZ* _model)
{
	model = _model;

	t_dim = model->t_dim;
	r_dim = model->R_dim;
	T_dim = model->T_dim;
	P_dim = model->P_dim;

	nx = model->cellsNum_r + 2;
	ny = model->cellsNum_z + 2;
}
void Snapshotter<wax_nit::WaxNIT>::setModel(wax_nit::WaxNIT* _model)
{
	model = _model;

	T_dim = model->T_dim;
	t_dim = model->t_dim;
	r_dim = model->R_dim;
	T_dim = model->T_dim;
	P_dim = model->P_dim;

	nx = model->cellsNum_r + 2;
	ny = model->cellsNum_z + 2;
}
void Snapshotter<wax_nit1d::WaxNIT1d>::setModel(wax_nit1d::WaxNIT1d* _model)
{
	model = _model;

	t_dim = model->t_dim;
	r_dim = model->R_dim;
	T_dim = model->T_dim;
	P_dim = model->P_dim;

	nx = model->cellsNum_x + 2;
	ny = 3;
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
template <class modelType>
string Snapshotter<modelType>::getFileName(int i, const string name)
{
	string filename = pattern;
	return replace(replace(filename, "%{STEP}", to_string(i)), "%{NAME}", name);
}

template class Snapshotter<gasOil_rz::GasOil_RZ>;
template class Snapshotter<acid1d::Acid1d>;
template class Snapshotter<acid2d::Acid2d>;
template class Snapshotter<acid2dnit::Acid2dNIT>;
template class Snapshotter<acidfrac::AcidFrac>;
template class Snapshotter<acidellfrac::AcidEllFrac>;
template class Snapshotter<acidrecfrac::AcidRecFrac>;
template class Snapshotter<acidrecfracmov::AcidRecFracMov>;
template class Snapshotter<vpp2d::VPP2d>;
template class Snapshotter<bing1d::Bingham1d>;
template class Snapshotter<gasOil_elliptic::GasOil_Elliptic>;
template class Snapshotter<oilnit_elliptic::OilNIT_Elliptic>;
template class Snapshotter<blackoilnit_elliptic::BlackOilNIT_Elliptic>;
template class Snapshotter<blackoil_rz::BlackOil_RZ>;
template class Snapshotter<wax_nit::WaxNIT>;
template class Snapshotter<wax_nit1d::WaxNIT1d>;