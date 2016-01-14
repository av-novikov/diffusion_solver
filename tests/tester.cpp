#include "gtest/gtest.h"

#include "tests/gas1D-test.h"
#include "tests/gas1Dsimple-test.h"

TEST(Gas1DTest, StationaryRate)
{
	Gas1D_Test test;
	test.run();
	test.test();
}

TEST(Gas1DSimpleTest, StationaryRate)
{
	Gas1D_Simple_Test test;
	test.run();
	test.test();
}