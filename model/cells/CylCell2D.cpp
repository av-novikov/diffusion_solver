#include "model/cells/CylCell2D.h"
#include "model/cells/Variables.hpp"

template <typename varType>
CylCell2D<varType>::CylCell2D()
{
}

template <typename varType>
CylCell2D<varType>::CylCell2D(int _num, double _r, double _z, double _hr, double _hz) : AbstractCell<varType>(_num), r(_r), z(_z), hr(_hr), hz(_hz)
{
	V = 2 * M_PI * r * hr * hz;
}

template <typename varType>
CylCell2D<varType>::~CylCell2D()
{
}

template class CylCell2D<Var1phase>;
template class CylCell2D<Var1phaseNIT>;
template class CylCell2D<Var2phase>;
template class CylCell2D<Var2phaseNIT>;