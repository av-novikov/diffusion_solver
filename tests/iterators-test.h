#ifndef ITERATORS_TEST_H_
#define ITERATORS_TEST_H_

#include "tests/base-test.h"
#include "Scene.h"
#include "model/3D/GasOil_3D/GasOil_3D.h"
#include "model/3D/GasOil_3D/Par3DSolver.h"

class Iterators_Test : public BaseTest<gasOil_3d::Properties, Scene<gasOil_3d::GasOil_3D, gasOil_3d::Par3DSolver, gasOil_3d::Properties> >
{
protected:
	gasOil_3d::Properties* getProps();

public:
	void test();
	void run();
};

#endif /* ITERATORS_TEST_H_ */
