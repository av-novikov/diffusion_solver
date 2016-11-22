#include "model/AbstractModel.hpp"
#include "model/cells/Variables.hpp"

#include "model/GasOil_RZ/GasOil_RZ.h"

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

template class AbstractModel<Var2phase, gasOil_rz::Properties, gasOil_rz::TCell, gasOil_rz::GasOil_RZ>;

template class AbstractModel<VarSimpleAcid, acid2d::Properties, acid2d::TCell, acid2d::Acid2d>;
template class AbstractModel<VarSimpleVPP, vpp2d::Properties, vpp2d::TCell, vpp2d::VPP2d>;