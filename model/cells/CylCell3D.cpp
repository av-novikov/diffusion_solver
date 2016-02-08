#include "model/cells/CylCell3D.h"
#include "model/cells/Variables.hpp"

template <typename varType>
CylCell3D<varType>::CylCell3D()
{
}

template <typename varType>
CylCell3D<varType>::CylCell3D(int _num, double _r, double _phi, double _z, double _hr, double _hphi, double _hz) : AbstractCell<varType>(_num), r(_r), phi(_phi), z(_z), hr(_hr), hphi(_hphi), hz(_hz)
{
	V = hphi * r * hr * hz;
}

template <typename varType>
CylCell3D<varType>::~CylCell3D()
{
}

template class CylCell3D<Var1phase>;
template class CylCell3D<Var1phaseNIT>;
template class CylCell3D<Var2phase>;
template class CylCell3D<Var2phaseNIT>;