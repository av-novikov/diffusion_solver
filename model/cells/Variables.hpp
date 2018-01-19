#ifndef VARIABLES_HPP_
#define VARIABLES_HPP_

#include "adolc/adouble.h"

#define TEMP 0
#define PRES 1
#define SAT 2

struct Var1phase 
{
	static const int size = 1;

	union {
		double values [1];
		struct 
		{
			double p;
		};
	};
};
struct TapeVar1Phase
{
	static const int size = 1;
	adouble p;
};
struct Var1phaseNIT 
{
	static const int size = 2;

	union {
		double values [2];
		struct 
		{
			double t;
			double p;
		};
	};
};
struct TapeVar1PhaseNIT
{
	static const int size = 1;
	adouble t;
};
struct Var2phase 
{
	static const int size = 3;

	union {
		double values [3];
		struct 
		{
			double p;
			double s;
			double p_bub;
		};
	};
	bool SATUR;
};
struct TapeVarGasOil
{
	static const int size = 3;

	adouble p;
	adouble s;
	adouble p_bub;

	bool SATUR;
};
struct Var2phaseNIT
{
	static const int size = 4;

	union {
		double values [4];
		struct 
		{
			double t;
			double p;
			double s;
			double p_bub;
		};
	};
	bool SATUR;
};

struct VarBlackOil
{
	static const int size = 4;

	union {
		double values[4];
		struct
		{
			double p;
			double s_w;
			double s_o;
			double p_bub;
		};
	};
	bool SATUR;
};
struct TapeVarBlackOil
{
	static const int size = 4;

	adouble p;
	adouble s_w;
	adouble s_o;
	adouble p_bub;

	bool SATUR;
};
struct VarBlackOilNIT
{
	static const int size = 5;

	union {
		double values[5];
		struct
		{
			double t;
			double p;
			double s_w;
			double s_o;
			double p_bub;
		};
	};
	bool SATUR;
};

struct VarWaxNIT
{
	static const int size = 6;

	union {
		double values[6];
		struct
		{
			double m;
			double t;
			double p;
			double s_w;
			double s_o;
			double p_bub;
		};
	};
	bool SATUR;
};
struct TapeVarWaxNIT
{
	static const int size = 6;
	
	adouble m;
	adouble t;
	adouble p;
	adouble s_w;
	adouble s_o;
	adouble p_bub;
	bool SATUR;
};

template <typename DataType>
struct VarSimpleVPP
{
	static const int size = 3;

	union {
		DataType values[size];
		struct
		{
			DataType p;
			DataType s;
			DataType c;
		};
	};

	template <class OtherType>
	VarSimpleVPP<DataType>& operator=(const VarSimpleVPP<OtherType>& other)
	{
		for (int i = 0; i < size; i++)
			values[i] = other.values[i];

		return *this;
	};

	VarSimpleVPP() {};
	~VarSimpleVPP() {};
};
struct TapeVarSimpleVPP
{
	static const int size = 3;

	adouble p;
	adouble s;
	adouble c;
};

#endif /* VARIABLES_HPP_ */
