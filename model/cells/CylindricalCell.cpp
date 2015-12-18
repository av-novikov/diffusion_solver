#include "model/cells/CylindricalCell.h"
#include "model/cells/Variables.hpp"

template <typename varType>
CylCell<varType>::CylCell()
{
}

template <typename varType>
CylCell<varType>::CylCell(int _num, double _r, double _z, double _hr, double _hz) : AbstractCell<varType>(_num), r(_r), z(_z), hr(_hr), hz(_hz)
{
	V = 2 * M_PI * r * hr * hz;
}

template <typename varType>
CylCell<varType>::~CylCell()
{
}

template class CylCell<Var1phase>;
template class CylCell<Var1phaseNIT>;
template class CylCell<Var2phase>;
template class CylCell<Var2phaseNIT>;