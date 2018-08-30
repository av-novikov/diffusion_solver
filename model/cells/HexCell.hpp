#ifndef HEXCELL_HPP_
#define HEXCELL_HPP_

#include "model/cells/AbstractCell.hpp"

template <typename varType, typename PropsType = EmptyStruct>
class HexCell : public AbstractCell<varType>
{
public:
	double x, y, z;
	double hx, hy, hz;
	PropsType* props;

	HexCell() {};
	HexCell(int _num, double _x, double _y, double _z,
					double _hx, double _hy, double _hz) :
		AbstractCell<varType>(_num, -1), x(_x), y(_y), z(_z), hx(_hx), hy(_hy), hz(_hz)
	{
		V = hx * hy * hz;
	};
	HexCell(int _num, double _x, double _y, double _z,
					double _hx, double _hy, double _hz, const Type _type) :
		AbstractCell<varType>(_num, _type), x(_x), y(_y), z(_z), hx(_hx), hy(_hy), hz(_hz)
	{
		V = hx * hy * hz;
	};
	~HexCell() {};
};

#endif /* HEXCELL_HPP_ */
