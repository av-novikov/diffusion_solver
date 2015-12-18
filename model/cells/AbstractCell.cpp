#include "model/cells/AbstractCell.hpp"
#include "model/cells/Variables.hpp"

template<>
const int AbstractCell<Var1phase>::varNum = 1;

template<>
const int AbstractCell<Var1phaseNIT>::varNum = 2;

template<>
const int AbstractCell<Var2phase>::varNum = 3;

template<>
const int AbstractCell<Var2phaseNIT>::varNum = 4;


template <typename varType>
AbstractCell<varType>::AbstractCell()
{
	V = 0.0;
	num = -1;
}

template <typename varType>
AbstractCell<varType>::AbstractCell(int _num) : num(_num)
{
}

template <typename varType>
AbstractCell<varType>::~AbstractCell()
{
}

template class AbstractCell<Var1phase>;
template class AbstractCell<Var1phaseNIT>;
template class AbstractCell<Var2phase>;
template class AbstractCell<Var2phaseNIT>;