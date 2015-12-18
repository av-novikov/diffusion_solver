#include "model/cells/RadialCell.hpp"
#include "model/cells/Variables.hpp"

template <typename varType>
RadialCell<varType>::RadialCell()
{
}

template <typename varType>
RadialCell<varType>::RadialCell(int _num, double _r, double _hr, double _hz) : AbstractCell<varType>(_num), r(_r), hr(_hr) 
{
	V = 2.0 * M_PI * r * hr * _hz;
}

template <typename varType>
RadialCell<varType>::~RadialCell()
{
}

template class RadialCell<Var1phase>;
template class RadialCell<Var1phaseNIT>;
template class RadialCell<Var2phase>;
template class RadialCell<Var2phaseNIT>;