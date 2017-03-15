#ifndef CYLCELL2D_H_
#define CYLCELL2D_H_

#define PREV 0
#define ITER 1
#define NEXT 2

#include "model/cells/AbstractCell.hpp"

template <typename varType>
class CylCell2D : public AbstractCell<varType>
{
public:
	double r;
	double z;

	double hr;
	double hz;

	CylCell2D() {};
	CylCell2D(int _num, double _r, double _z, double _hr, double _hz) :
				AbstractCell<varType>(_num, -1), r(_r), z(_z), hr(_hr), hz(_hz)
				{ V = 2 * M_PI * r * hr * hz; };
	CylCell2D(int _num, double _r, double _z, double _hr, double _hz, const Type _type) :
		AbstractCell<varType>(_num, _type), r(_r), z(_z), hr(_hr), hz(_hz)
				{ V = 2 * M_PI * r * hr * hz;	};
	~CylCell2D() {};
};

template <typename varType, typename PropsType = EmptyStruct>
class NewCylCell2D : public AbstractCell<varType>
{
public:
	double r;
	double z;

	double hr;
	double hz;

	PropsType* props;

	NewCylCell2D() {};
	NewCylCell2D(int _num, double _r, double _z, double _hr, double _hz) :
				AbstractCell<varType>(_num), r(_r), z(_z), hr(_hr), hz(_hz)
				{ V = 2 * M_PI * r * hr * hz; };
	NewCylCell2D(int _num, double _r, double _z, double _hr, double _hz, const Type _type) :
		AbstractCell<varType>(_num, _type), r(_r), z(_z), hr(_hr), hz(_hz)
	{
		V = 2 * M_PI * r * hr * hz;
	};
	~NewCylCell2D() {};
};

namespace std {
	template <typename varType>
	inline std::ostream& operator<<(std::ostream& os, const CylCell2D<varType>& cell)
	{
		os << "\nCell number:\t" << cell.num << std::endl;
		os << "Mass center:\t" << cell.r << "\t" << cell.z << std::endl;
		os << "Values:\t";
		for(int i = 0; i < varType::varNum; i++)
			os << cell.u_next.values[i] << "\t";
		os << std::endl;

		return os;
	}
}

#endif /* CYLCELL2D_H_ */