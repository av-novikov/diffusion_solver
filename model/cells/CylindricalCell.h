#ifndef CYLINDRICALCELL_H_
#define CYLINDRICALCELL_H_

#define VAR_NUM 4
#define PREV 0
#define ITER 1
#define NEXT 2

#define R_AXIS 0
#define Z_AXIS 1

#include "model/cells/AbstractCell.hpp"

template <typename varType>
class CylCell : public AbstractCell<varType>
{
public:
	double r;
	double z;

	double hr;
	double hz;

	CylCell();
	CylCell(int _num, double _r, double _z, double _hr, double _hz);
	~CylCell();
};

namespace std {
	template <typename varType>
	inline std::ostream& operator<<(std::ostream& os, const CylCell<varType>& cell)
	{
		os << "\nCell number:\t" << cell.num << std::endl;
		os << "Mass center:\t" << cell.r << "\t" << cell.z << std::endl;
		os << "Values:\t";
		for(int i = 0; i < VAR_NUM; i++)
			os << cell.u_next.values[i] << "\t";
		os << std::endl;

		return os;
	}
}

#endif /* CYLINDRICALCELL_H_ */