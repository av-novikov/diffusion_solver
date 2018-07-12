#ifndef POINT_HPP_
#define POINT_HPP_

#include <cmath>
#include <iostream>
#include <vector>
#include <new>

#include "adolc/drivers/drivers.h"
#include "adolc/adolc.h"

#include "util/utils.h"

namespace point
{
	enum PointType { NOTYPE, INNER, BORDER };
	struct CartesianPoint3d
	{
		PointType type;
		const int id;
		union
		{
			double coords[3];
			struct
			{
				double x;	double y;	double z;
			};
		};

		CartesianPoint3d() : id(-1), type(NOTYPE) {};
		CartesianPoint3d(const double _x, const double _y, const double _z) : id(-1), x(_x), y(_y), z(_z), type(NOTYPE) { };
		CartesianPoint3d(const int _id, const double _x, const double _y, const double _z) : id(_id), x(_x), y(_y), z(_z), type(NOTYPE) { };
		~CartesianPoint3d() {};
		CartesianPoint3d(const CartesianPoint3d& a) : id(a.id)
		{
			(*this) = a;
		};
		CartesianPoint3d& operator=(const CartesianPoint3d& rhs)
		{
			x = rhs.x, y = rhs.y, z = rhs.z;
			type = rhs.type;
			return *this;
		};
		CartesianPoint3d& operator/=(const double k)
		{
			x /= k;	y /= k;	z /= k;
			return *this;
		};
		CartesianPoint3d& operator+=(const CartesianPoint3d& rhs)
		{
			x += rhs.x;	y += rhs.y;	z += rhs.z;
			return *this;
		};
		CartesianPoint3d& operator-=(const CartesianPoint3d& rhs)
		{
			x -= rhs.x;	y -= rhs.y;	z -= rhs.z;
		};

		inline std::array<double, 3> getMetric()
		{
			return {1.0, 1.0, 1.0};
		};
		inline double norm() const { return sqrt(x * x + y * y + z * z); };
	};
	struct EllipticPoint3d
	{
		PointType type;
		const int id;
		union
		{
			double coords[3];
			struct
			{
				double mu;	double nu;	double z;
			};
		};
		static double a;

		EllipticPoint3d() : id(-1), type(NOTYPE) {};
		EllipticPoint3d(const double _mu, const double _nu, const double _z) : id(-1), mu(_mu), nu(_nu), z(_z), type(NOTYPE) { };
		EllipticPoint3d(const int _id, const double _mu, const double _nu, const double _z) : id(_id), mu(_mu), nu(_nu), z(_z), type(NOTYPE) { };
		~EllipticPoint3d() {};
		EllipticPoint3d(const CartesianPoint3d& a) : id(a.id)
		{
			(*this) = a;
		};
		EllipticPoint3d& operator=(const EllipticPoint3d& rhs)
		{
			mu = rhs.mu, nu = rhs.nu, z = rhs.z;
			type = rhs.type;
			return *this;
		};
		EllipticPoint3d& operator/=(const double k)
		{
			mu /= k;	nu /= k;	z /= k;
			return *this;
		};
		EllipticPoint3d& operator+=(const EllipticPoint3d& rhs)
		{
			mu += rhs.mu;	nu += rhs.nu;	z += rhs.z;
			return *this;
		};
		EllipticPoint3d& operator-=(const EllipticPoint3d& rhs)
		{
			mu -= rhs.mu;	nu -= rhs.nu;	z -= rhs.z;
		};

		inline double getH() const
		{
			return a * sqrt(sinh(mu) * sinh(mu) + sin(nu) * sin(nu));
		};
		inline std::array<double, 3> getMetric() const
		{
			const double h = getH();
			return{h * h, h * h, 1.0};
		};
		inline double dist() const
		{
			return sqrt(mu * mu + nu * nu + z * z);
		}
		inline double norm() const
		{ 
			const auto g = getMetric();
			return sqrt( g[0] * mu * mu + g[1] * nu * nu + g[2] * z * z ); 
		};
		inline CartesianPoint3d getCartesian() const
		{
			return CartesianPoint3d(id, a * cosh(mu) * cos(nu), a * sinh(mu) * sin(nu), z);
		};
		inline CartesianPoint3d getVectorCartesian(double u_mu, double u_nu, double u_z) const
		{
			const auto g = getMetric();
			const double u_x = a * (sinh(mu) * cos(nu) * u_mu - cosh(mu) * sin(nu) * u_nu) / sqrt(g[0]);
			const double u_y = a * (cosh(mu) * sin(nu) * u_mu - sinh(mu) * cos(nu) * u_nu) / sqrt(g[1]);
			return{ u_x, u_y, u_z };
		};
	};

	inline CartesianPoint3d getCartesianFromElliptic(const double mu, const double nu, const double z)
	{
		return{ EllipticPoint3d::a * cosh(mu) * cos(nu), 
				EllipticPoint3d::a * sinh(mu) * sin(nu), 
				z };
	};

	template <typename Point>
	inline std::ostream& operator<<(std::ostream& os, const Point& a)
	{
		os << a.coords[0] << " " << a.coords[1] << " " << a.coords[2] << std::endl;
		return os;
	}
	/*template <typename Point>
	inline bool operator==(const Point& a1, const Point& a2)
	{
		if ((fabs(a2.coords[0] - a1.coords[0]) > EQUALITY_TOLERANCE) ||
			(fabs(a2.coords[1] - a1.coords[1]) > EQUALITY_TOLERANCE) ||
			(fabs(a2.coords[2] - a1.coords[2]) > EQUALITY_TOLERANCE))
			return false;
		else
			return true;
	};*/
	template <typename Point>
	inline Point operator-(const Point& rhs)
	{
		return Point(-rhs.coords[0], -rhs.coords[1], -rhs.coords[2]);
	};
	template <typename Point>
	inline Point operator-(const Point& a1, const Point& a2)
	{
		return Point(a1.coords[0] - a2.coords[0], a1.coords[1] - a2.coords[1], a1.coords[2] - a2.coords[2]);
	};
	template <typename Point>
	inline Point operator+(const Point& rhs)
	{
		return Point(rhs.coords[0], rhs.coords[1], rhs.coords[2]);
	};
	template <typename Point>
	inline Point operator+(const Point& a1, const Point& a2)
	{
		return Point(a1.coords[0] + a2.coords[0], a1.coords[1] + a2.coords[1], a1.coords[2] + a2.coords[2]);
	};
	template <typename Point>
	inline Point operator*(const Point& a1, double k)
	{
		return Point(a1.coords[0] * k, a1.coords[1] * k, a1.coords[2] * k);
	};
	template <typename Point>
	inline Point operator*(double k, const Point& a1)
	{
		return a1 * k;
	};
	template <typename Point>
	inline Point operator/(const Point& a1, double k)
	{
		return Point(a1.coords[0] / k, a1.coords[1] / k, a1.coords[2] / k);
	};
	template <typename Point>
	inline Point operator/(const Point& a1, const Point& a2)
	{
		return Point(a1.coords[0] / a2.coords[0], a1.coords[1] / a2.coords[1], a1.coords[2] / a2.coords[2]);
	};
	template <typename Point>
	inline Point operator*(const Point& a1, const Point& a2)
	{
		return Point(a1.coords[0] * a2.coords[0], a1.coords[1] * a2.coords[1], a1.coords[2] * a2.coords[2]);
	};
	template <typename Point>
	inline double dot_product(const Point& a1, const Point& a2)
	{
		return a1.coords[0] * a2.coords[0] + a1.coords[1] * a2.coords[1] + a1.coords[2] * a2.coords[2];
	};
	template <typename Point>
	inline double dot_product(const Point& a1, const Point& a2, const Point& in)
	{
		const auto g = in.getMetric();
		return g[0] * a1.coords[0] * a2.coords[0] + g[1] * a1.coords[1] * a2.coords[1] + g[2] * a1.coords[2] * a2.coords[2];
	};

	/*template <typename Point>
	inline Point vector_product(const Point& a1, const Point& a2)
	{
		return{ a1.y * a2.z - a1.z * a2.y,
			a1.z * a2.x - a1.x * a2.z,
			a1.x * a2.y - a1.y * a2.x };
	};
	template <typename Point>
	inline double distance(const Point& a1, const Point& a2)
	{
		return sqrt(dot_product(a2 - a1, a2 - a1));
	};
	template <typename Point>
	inline double square(const Point& a1, const Point& a2, const Point& a3)
	{
		const double a = distance(a1, a2);
		const double b = distance(a2, a3);
		const double c = distance(a3, a1);
		const double p = (a + b + c) / 2.0;
		return sqrt(p * (p - a) * (p - b) * (p - c));
	};
	template <typename Point>
	inline double square(const Point& a1, const Point& a2, const Point& a3, const Point& a4)
	{
		return square(a1, a2, a3) + square(a3, a4, a1);
	};*/
};

#endif /* POINT_HPP_ */
