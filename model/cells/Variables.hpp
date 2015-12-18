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
	
	Var1phase& operator=(const Var1phase& a)
	{
		p = a.p;
		return *this;
	}
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
	
	Var1phaseNIT& operator=(const Var1phaseNIT& a)
	{
		p = a.p;
		t = a.t;
		return *this;
	}
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
	
	Var2phase& operator=(const Var2phase& a)
	{
		p = a.p;
		s = a.s;
		p_bub = a.p_bub;
		SATUR = a.SATUR;
		
		return *this;
	}
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
	
	Var2phaseNIT& operator=(const Var2phaseNIT& a)
	{
		p = a.p;
		s = a.s;
		p_bub = a.p_bub;
		t = a.t;
		SATUR = a.SATUR;
		
		return *this;
	}
};

#endif /* VARIABLES_HPP_ */
