#ifndef ACID2D_PROPERTIES_HPP_
#define ACID2D_PROPERTIES_HPP_

#include "model/Basic2d/Properties.hpp"

#include "adolc/adouble.h"
#include "adolc/taping.h"

namespace acid2d
{
	// ADOLC stencil ids
	const int mid = basic2d::mid;
	const int left = basic2d::left;
	const int right = basic2d::right;
	const int vertical = basic2d::vertical;

	struct Mineral
	{

	};
	struct Skeleton_Props : public basic2d::Skeleton_Props
	{
		std::vector<Mineral> minerals;

		// Initial values
		double m_init;
		double p_init;
		double s_init;
		double Ya_init;
		double Ys_init;
	};
	struct Component
	{
	};
	struct Liquid_Props : public basic2d::Liquid_Props
	{
		std::vector<Component> components;
	};
	struct Gas_Props : public basic2d::Gas_Props
	{
		std::vector<Component> components;
	};
	struct Properties : public basic2d::Properties
	{
		std::vector<Skeleton_Props> props_sk;
		Liquid_Props props_l;
		Gas_Props props_g;

		// Data set (saturation, relative oil permeability)
		std::vector< std::pair<double, double> > kr_l;
		// Data set (saturation, relative gas permeability)
		std::vector< std::pair<double, double> > kr_g;
	};
};

#endif /* ACID2D_PROPERTIES_HPP_ */