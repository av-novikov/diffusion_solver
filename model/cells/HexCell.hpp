#ifndef HEXCELL_HPP_
#define HEXCELL_HPP_

#include "model/cells/AbstractCell.hpp"
#include "model/cells/Variables.hpp"

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
    template <typename TVar, typename TProp,
        typename = typename std::enable_if<std::is_same<VarOilWater, varType>::value
                                        && std::is_same<JustAcid, TVar>::value 
                                        && std::is_same<PropsType, TProp>::value>::type>
    HexCell(const HexCell<TVar,TProp>& rhs, const int _num, const Type _type) : AbstractCell<varType>(_num, _type)
    {
        x = rhs.x;      y = rhs.y;      z = rhs.z;
        hx = rhs.hx;    hy = rhs.hy;    hz = rhs.hz;
        props = rhs.props;
        V = rhs.V;
    };
    ~HexCell() {};
};

#endif /* HEXCELL_HPP_ */
