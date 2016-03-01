#include "model/cells/CylCellPerf.h"
#include "model/cells/Variables.hpp"

template <typename varType>
CylCellPerf<varType>::CylCellPerf()
{
}

template <typename varType>
CylCellPerf<varType>::CylCellPerf(int _num, double _r, double _phi, double _z, double _hr, double _hphi, double _hz) : AbstractCell<varType>(_num), r(_r), phi(_phi), z(_z), hr(_hr), hphi(_hphi), hz(_hz)
{
	V = hphi * r * hr * hz;
	isUsed = true;
	isTunnel = false;
}

template <typename varType>
CylCellPerf<varType>::CylCellPerf(int _num, double _r, double _phi, double _z, double _hr, double _hphi, double _hz, int _tunNum) : AbstractCell<varType>(_num), r(_r), phi(_phi), z(_z), hr(_hr), hphi(_hphi), hz(_hz), tunNum(_tunNum)
{
	V = hphi * r * hr * hz;
	isUsed = true;
	isTunnel = true;
}

template <typename varType>
CylCellPerf<varType>::~CylCellPerf()
{
}

template class CylCellPerf<Var1phase>;
template class CylCellPerf<Var1phaseNIT>;
template class CylCellPerf<Var2phase>;
template class CylCellPerf<Var2phaseNIT>;