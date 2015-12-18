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
		static const int varNum;
		
		double V;
		
		varType u_next;
		varType u_iter;
		varType u_prev;
		
		AbstractCell();
		AbstractCell(int _num);
		~AbstractCell();
};

#endif /* ABSTRACTCELL_HPP_ */
