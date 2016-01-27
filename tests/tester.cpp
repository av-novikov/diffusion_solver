#include "gtest/gtest.h"

#include "tests/oil1D-test.h"
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

TEST(Oil1D, StationaryRate)
{
	Oil1D_Test test;
	test.run();
	test.stat_test();
}

TEST(Oil1D, NonStationaryRate)
{
	Oil1D_Test test;
	test.run();
	test.non_stat_test();
}

TEST(Oil1D_NIT, JT_Temperature)
{
	Oil1D_Test test;
	test.run();
	test.jt_temp_test();
}