#ifndef VARIABLES_HPP_
#define VARIABLES_HPP_

struct Var1phase 
{
	union {
		double values [1];
		struct 
		{
			double p;
		};
	};
};

struct Var1phaseNIT 
{
	union {
		double values [2];
		struct 
		{
			double t;
			double p;
		};
	};
};

struct Var2phase 
{
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

struct Var2phaseNIT
{
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

struct VarSimpleAcid
{
	union {
		double values[5];
		struct
		{
			double m;
			double p;
			double s;
			double Ys;
			double Ya;
		};
	};
};

#endif /* VARIABLES_HPP_ */
