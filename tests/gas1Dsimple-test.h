#ifndef GAS1DSIMPLE_TEST_H_
#define GAS1DSIMPLE_TEST_H_

#include "tests/base-test.h"
#include "Scene.h"
#include "model/Gas1D/Gas1D.h"
#include "model/Gas1D/Gas1DSolver.h"

class Gas1D_Simple_Test : public BaseTest<gas1D::Properties, Scene<gas1D::Gas1D_simple, gas1D::Gas1DSolSimp, gas1D::Properties> >
{
protected:
	gas1D::Properties* getProps();

public:
	void test();
};

#endif /* GAS1DSIMPLE_TEST_H_ */