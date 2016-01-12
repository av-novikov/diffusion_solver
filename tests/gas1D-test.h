#ifndef GAS1D_TEST_H_
#define GAS1D_TEST_H_

#include "tests/base-test.h"
#include "Scene.h"
#include "model/Gas1D/Gas1D.h"
#include "model/Gas1D/Gas1DSolver.h"

class Gas1D_Test : public BaseTest<gas1D::Properties, Scene<gas1D::Gas1D, gas1D::Gas1DSol, gas1D::Properties> >
{
protected:
	gas1D::Properties* getProps();

public:
	void test();
};

#endif /* GAS1D_TEST_H_ */