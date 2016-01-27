#ifndef OIL1D_TEST_H_
#define OIL1D_TEST_H_

#include "tests/base-test.h"
#include "Scene.h"
#include "model/Oil1D_NIT/Oil1D_NIT.h"
#include "model/Oil1D_NIT/Oil1DNITSolver.h"

class Oil1D_Test : public BaseTest<oil1D_NIT::Properties, Scene<oil1D_NIT::Oil1D_NIT, oil1D_NIT::Oil1DNITSolver, oil1D_NIT::Properties> >
{
protected:
	oil1D_NIT::Properties* getProps();

public:
	void test();
	void stat_test();
	void non_stat_test();
	void jt_temp_test();
};

#endif /* OIL1D_TEST_H_ */