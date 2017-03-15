#ifndef BLACKOILRZ_PROPERTIES_HPP_
#define BLACKOILRZ_PROPERTIES_HPP_

#include <array>
#include "model/Basic2d/Properties.hpp"

namespace blackoil_rz
{
	const int mid = basic2d::mid;
	const int left = basic2d::left;
	const int right = basic2d::right;
	const int vertical = basic2d::vertical;

	struct Skeleton_Props : public basic2d::Skeleton_Props
	{
		// Initial values
		double p_init;
		double so_init;
		double sw_init;
	};
	struct Liquid_Props : public basic2d::Liquid_Props
	{
	};
	struct Gas_Props : public basic2d::Gas_Props
	{
	};

	struct Properties : public basic2d::Properties
	{
		std::vector<Skeleton_Props> props_sk;
		Liquid_Props props_oil;
		Liquid_Props props_wat;
		Gas_Props props_gas;
	};
};

#endif /* BLACKOILRZ_PROPERTIES_HPP_ */