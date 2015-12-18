#include "model/AbstractModel.hpp"
#include "model/cells/Variables.hpp"

#include "model/Oil1D/Oil1D.h"
#include "model/Gas1D/Gas1D.h"
#include "model/Oil1D_NIT/Oil1D_NIT.h"
#include "model/Oil_RZ/Oil_RZ.h"
#include "model/GasOil_RZ/GasOil_RZ.h"
#include "model/GasOil_RZ_NIT/GasOil_RZ_NIT.h"

using namespace std;

template <typename varType, typename propsType,
template <typename varType> class cellType>
AbstractModel<varType, propsType, cellType>::AbstractModel()
{
}

template <typename varType, typename propsType,
template <typename varType> class cellType>
AbstractModel<varType, propsType, cellType>::~AbstractModel()
{
}

template <typename varType, typename propsType,
template <typename varType> class cellType>
void AbstractModel<varType, propsType, cellType>::setInitialState()
{
	vector<cellType<varType> >::iterator it;
	for(it = cells.begin(); it != cells.end(); ++it)
		it->u_prev = it->u_iter = it->u_next = varInit;

	//for(auto cell : cells)
	//	cell.u_prev = cell.u_iter = cell.u_next = varInit;
}

template <typename varType, typename propsType,
template <typename varType> class cellType>
void AbstractModel<varType, propsType, cellType>::load(propsType& props)
{
	setProps(props);

	buildGridLog();
	setPerforated();
	setInitialState();
}

template <typename varType, typename propsType,
template <typename varType> class cellType>
int AbstractModel<varType, propsType, cellType>::getCellsNum()
{
	return cellsNum;
}

template class AbstractModel<Var1phase, oil1D::Properties, RadialCell>;
template class AbstractModel<Var1phase, gas1D::Properties, RadialCell>;
template class AbstractModel<Var1phaseNIT, oil1D_NIT::Properties, RadialCell>;
template class AbstractModel<Var1phase, oil_rz::Properties, CylCell>;
template class AbstractModel<Var2phaseNIT, gasOil_rz_NIT::Properties, CylCell>;
template class AbstractModel<Var2phase, gasOil_rz::Properties, CylCell>;