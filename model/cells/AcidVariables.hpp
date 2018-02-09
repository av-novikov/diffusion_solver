#ifndef ACIDVARIABLES_HPP_
#define ACIDVARIABLES_HPP_

#include "adolc/adouble.h"

struct FirstAcid
{
	static const int size = 7;

	union {
		double values[size];
		struct
		{
			double m;
			double p;
			double sw;
			double xa;
			double xw;
			double so;
			double p_bub;
		};
	};

	bool SATUR;
};
struct TapeFirstAcid
{
	static const int size = 7;

	adouble m;
	adouble p;
	adouble sw;
	adouble xa;
	adouble xw;
	adouble so;
	adouble p_bub;
	bool SATUR;
};
struct JustAcid
{
	static const int size = 6;

	union {
		double values[size];
		struct
		{
			double m;
			double p;
			double sw;
			double xw;
			double xa;
			double xs;
		};
	};
};
struct TapeJustAcid
{
	static const int size = 6;

	adouble m;
	adouble p;
	adouble sw;
	adouble xw;
	adouble xa;
	adouble xs;
};
struct JustAcidNIT
{
	static const int size = 7;

	union {
		double values[size];
		struct
		{
			double m;
			double p;
			double sw;
			double xw;
			double xa;
			double xs;
			double t;
		};
	};
};
struct TapeJustAcidNIT
{
	static const int size = 7;

	adouble m;
	adouble p;
	adouble sw;
	adouble xw;
	adouble xa;
	adouble xs;
	adouble t;
};
struct VarFrac
{
	static const int size = 2;

	union {
		double values[size];
		struct
		{
			double p;
			double c;
		};
	};
};
struct TapeVarFrac
{
	static const int size = 2;
	adouble p;
	adouble c;
};

#endif /* ACIDVARIABLES_HPP_ */
