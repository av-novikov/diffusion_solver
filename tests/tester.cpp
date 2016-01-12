#include "gtest/gtest.h"

#include "tests/gas1D-test.h"
#include "tests/gas1Dsimple-test.h"

TEST(Gas1DTest, Gas1D)
{
	Gas1D_Test test;
	test.run();
	test.test();
}

TEST(Gas1DSimpleTest, Gas1DSimple)
{
	Gas1D_Simple_Test test;
	test.run();
	test.test();
}