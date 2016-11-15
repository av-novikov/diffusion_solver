#include "model/AbstractModel.hpp"
#include "model/cells/Variables.hpp"

#include "model/Oil1D/Oil1D.h"
#include "model/Gas1D/Gas1D.h"
#include "model/Gas1D/Gas1D_simple.h"
#include "model/Oil1D_NIT/Oil1D_NIT.h"
#include "model/Oil_RZ/Oil_RZ.h"
#include "model/Oil_RZ_NIT/Oil_RZ_NIT.h"
#include "model/GasOil_RZ/GasOil_RZ.h"
#include "model/GasOil_RZ_NIT/GasOil_RZ_NIT.h"

#include "model/3D/GasOil_3D/GasOil_3D.h"
#include "model/3D/GasOil_3D_NIT/GasOil_3D_NIT.h"

#include "model/3D/Perforation/GasOil_Perf.h"
#include "model/3D/Perforation/Oil_Perf_NIT.h"
#include "model/3D/Perforation/GasOil_Perf_NIT.h"

#include "model/Acid/2d/Acid2d.hpp"
#include "model/VPP2d/VPP2d.hpp"

using namespace std;

template <typename varType, typename propsType,
template <typename varType> class cellType, class modelType>
AbstractModel<varType, propsType, cellType, modelType>::AbstractModel()
{
	isWriteSnaps = false;
}

template <typename varType, typename propsType,
template <typename varType> class cellType, class modelType>
AbstractModel<varType, propsType, cellType, modelType>::~AbstractModel()
{
}

template <typename varType, typename propsType,
template <typename varType> class cellType, class modelType>
void AbstractModel<varType, propsType, cellType, modelType>::setInitialState()
{
	vector<cellType<varType> >::iterator it;
	for(it = cells.begin(); it != cells.end(); ++it)
		it->u_prev = it->u_iter = it->u_next = varInit;
}

template <typename varType, typename propsType,
template <typename varType> class cellType, class modelType>
void AbstractModel<varType, propsType, cellType, modelType>::load(propsType& props)
{
	setProps(props);

	buildGridLog();
	setPerforated();
	setInitialState();
}

template <typename varType, typename propsType,
template <typename varType> class cellType, class modelType>
int AbstractModel<varType, propsType, cellType, modelType>::getCellsNum()
{
	return cellsNum;
}

template <typename varType, typename propsType,
template <typename varType> class cellType, class modelType>
void AbstractModel<varType, propsType, cellType, modelType>::setSnapshotter(string type, modelType* model)
{
	if(type == "VTK") {
		snapshotter = new VTKSnapshotter<modelType>();
		snapshotter->setModel(model);
		isWriteSnaps = true;
	} else if(type == "GRDECL") {
		snapshotter = new GRDECLSnapshotter<modelType>();
		snapshotter->setModel(model);
		isWriteSnaps = true;
	} else {
		isWriteSnaps = false;
	}
}

template <typename varType, typename propsType,
template <typename varType> class cellType, class modelType>
void AbstractModel<varType, propsType, cellType, modelType>::setWellborePeriod(int period, double cur_t)
{
}

template <typename varType, typename propsType,
template <typename varType> class cellType, class modelType>
void AbstractModel<varType, propsType, cellType, modelType>::snapshot(int i)
{
	snapshotter->dump(i);
}

template <typename varType, typename propsType,
template <typename varType> class cellType, class modelType>
void AbstractModel<varType, propsType, cellType, modelType>::snapshot_all(int i)
{
	snapshotter->dump_all(i);
}

template class AbstractModel<Var1phase, oil1D::Properties, RadialCell, oil1D::Oil1D>;
template class AbstractModel<Var1phase, gas1D::Properties, RadialCell, gas1D::Gas1D>;
template class AbstractModel<Var1phase, gas1D::Properties, RadialCell, gas1D::Gas1D_simple>;
template class AbstractModel<Var1phaseNIT, oil1D_NIT::Properties, RadialCell, oil1D_NIT::Oil1D_NIT>;
template class AbstractModel<Var1phase, oil_rz::Properties, CylCell2D, oil_rz::Oil_RZ>;
template class AbstractModel<Var1phaseNIT, oil_rz_nit::Properties, CylCell2D, oil_rz_nit::Oil_RZ_NIT>;
template class AbstractModel<Var2phase, gasOil_rz::Properties, CylCell2D, gasOil_rz::GasOil_RZ>;
template class AbstractModel<Var2phaseNIT, gasOil_rz_NIT::Properties, CylCell2D, gasOil_rz_NIT::GasOil_RZ_NIT>;

template class AbstractModel<Var2phase, gasOil_3d::Properties, CylCell3D, gasOil_3d::GasOil_3D>;
template class AbstractModel<Var2phaseNIT, gasOil_3d_NIT::Properties, CylCell3D, gasOil_3d_NIT::GasOil_3D_NIT>;

template class AbstractModel<Var2phase, gasOil_perf::Properties, CylCellPerf, gasOil_perf::GasOil_Perf>;
template class AbstractModel<Var1phaseNIT, oil_perf_nit::Properties, CylCellPerf, oil_perf_nit::Oil_Perf_NIT>;
template class AbstractModel<Var2phaseNIT, gasOil_perf_nit::Properties, CylCellPerf, gasOil_perf_nit::GasOil_Perf_NIT>;

template class AbstractModel<VarSimpleAcid, acid2d::Properties, acid2d::TCell, acid2d::Acid2d>;
template class AbstractModel<VarSimpleVPP, vpp2d::Properties, vpp2d::TCell, vpp2d::VPP2d>;