#ifndef GAS1D_TEST_H_
#define GAS1D_TEST_H_

#include "tests/base-test.h"
#include "Scene.h"
#include "model/Gas1D/Gas1D.h"
#include "model/Gas1D/Gas1DSolver.h"

class Gas1D_Wrapped : public gas1D::Gas1D
{
	friend class Gas1D_Test;
protected:
	inline double getVisc(double p) const
	{
		if (leftBoundIsRate)
			return 0.00001 / P_dim / t_dim * p / cells[0].u_next.p;
		else
			return 0.00001 / P_dim / t_dim * p / Pwf;
	};
	inline double getVisc_dp(double p) const
	{
		if (leftBoundIsRate)
			return 0.00001 / P_dim / t_dim / cells[0].u_next.p;
		else
			return 0.00001 / P_dim / t_dim / Pwf;
	};
	inline double getZ(double p) const
	{
		if (leftBoundIsRate)
			return sqrt( p / cells[0].u_next.p );
		else
			return sqrt( p / Pwf );
	};
	inline double getZ_dp(double p) const
	{
		if (leftBoundIsRate)
			return 1.0 / 2.0 / sqrt( p * cells[0].u_next.p );
		else
			return 1.0 / 2.0 / sqrt( p * Pwf );
	};
};

class Gas1D_Test : public BaseTest<gas1D::Properties, Scene<Gas1D_Wrapped, gas1D::Gas1DSol, gas1D::Properties> >
{
protected:
	gas1D::Properties* getProps();

	double getStatRate() const;

public:
	void test();
};

#endif /* GAS1D_TEST_H_ */