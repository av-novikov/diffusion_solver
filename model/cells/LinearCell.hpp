#ifndef LINEAR_CELL_HPP_
#define LINEAR_CELL_HPP_

#include "model/cells/AbstractCell.hpp"

template <typename varType, typename PropsType = EmptyStruct>
class LinearCell : public AbstractCell<varType>
{
public:
	double x;
	double hx;

	PropsType* props;

	LinearCell() {};
	LinearCell(int _num, double _x, double _hx, double hz) :
		AbstractCell<varType>(_num), x(_x), hx(_hx)
	{
		V = hx * hz;
	};
	LinearCell(int _num, double _x,  double _hx, double hz, const Type _type) :
		AbstractCell<varType>(_num, _type), x(_x), hx(_hx)
	{
		V = hx * hz;
	};
	~LinearCell() {};
};

#endif /* LINEAR_CELL_HPP_ */