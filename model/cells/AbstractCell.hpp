#ifndef ABSTRACTCELL_HPP_
#define ABSTRACTCELL_HPP_

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>

template <typename varType>
class AbstractCell
{
	public:		
		int num;
		static const int varNum = varType::size;
		
		double V;
		
		varType u_next;
		varType u_iter;
		varType u_prev;
		
		AbstractCell() : V(0.0), num(-1) {};
		AbstractCell(int _num) : num(_num) {};
		~AbstractCell() {};
};

#endif /* ABSTRACTCELL_HPP_ */
