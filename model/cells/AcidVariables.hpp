#ifndef ACIDVARIABLES_HPP_
#define ACIDVARIABLES_HPP_

#include "adolc/adouble.h"

struct FirstAcid
{
	static const int size = 5;

	union {
		double values[size];
		struct
		{
			double m;
			double p;
			double sw;
			double xa;
			double xw;
		};
	};
};
struct TapeFirstAcid
{
	static const int size = 5;

	adouble m;
	adouble p;
	adouble sw;
	adouble xa;
	adouble xw;
};

#endif /* ACIDVARIABLES_HPP_ */
