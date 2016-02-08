#ifndef CYLCELL3D_H_
#define CYLCELL3D_H_

#define PREV 0
#define ITER 1
#define NEXT 2

#define R_AXIS 0
#define Z_AXIS 1
#define PHI_AXIS 2

#include "model/cells/AbstractCell.hpp"

template <typename varType>
class CylCell3D : public AbstractCell<varType>
{
public:
	double r;
	double phi;
	double z;

	double hr;
	double hphi;
	double hz;

	CylCell3D();
	CylCell3D(int _num, double _r, double _phi, double _z, double _hr, double _hphi, double _hz);
	~CylCell3D();
};

namespace std {
	template <typename varType>
	inline std::ostream& operator<<(std::ostream& os, const CylCell3D<varType>& cell)
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

#endif /* CYLCELL3D_H_ */