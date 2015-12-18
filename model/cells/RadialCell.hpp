#ifndef RADIALCELL_HPP_
#define RADIALCELL_HPP_

#include "model/cells/AbstractCell.hpp"

template <typename varType>
class RadialCell : public AbstractCell<varType>
{
	public:
		double r;
		double hr;
		
		RadialCell();
		RadialCell(int _num, double _r, double _hr, double _hz);
		~RadialCell();			
};

#endif /* RADIALCELL_HPP_ */
