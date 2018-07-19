#ifndef ELLIPTICCELL_HPP_
#define ELLIPTICCELL_HPP_

#include "model/cells/AbstractCell.hpp"

template <typename varType, typename PropsType = EmptyStruct>
class EllipticCell : public AbstractCell<varType>
{
public:
	typedef TPoint<3> Point;

	static double a;

	double mu;
	double nu;
	double z;

	double hmu;
	double hnu;
	double hz;

	PropsType* props;
	bool isUsed;

	inline std::array<double, 3> getCartesian() const
	{
		return { a * cosh(mu) * cos(nu), a * sinh(mu) * sin(nu), z };
	};
	inline std::array<double, 3> getVectorCartesian(double u_mu, double u_nu, double u_z) const
	{
		const double u_x = a * (sinh(mu) * cos(nu) * u_mu - cosh(mu) * sin(nu) * u_nu) / getH(mu, nu);
		const double u_y = a * (cosh(mu) * sin(nu) * u_mu + sinh(mu) * cos(nu) * u_nu) / getH(mu, nu);
		return{ u_x, u_y, u_z };
	};

	EllipticCell() : isUsed(true) {};
	EllipticCell(int _num, double _mu, double _nu, double _z, double _hmu, double _hnu, double _hz, const Type _type) : AbstractCell<varType>(_num, _type), mu(_mu), nu(_nu), z(_z), hmu(_hmu), hnu(_hnu), hz(_hz)
	{
		V = a * a * (sinh(mu) * sinh(mu) + sin(nu) * sin(nu)) * hmu * hnu * hz;

		isUsed = true;
	};
	~EllipticCell() {};

	static double getH(const double mu, const double nu)
	{
		return a * sqrt(sinh(mu) * sinh(mu) + sin(nu) * sin(nu));
	};
	bool operator==(const EllipticCell& rhs) const
	{
		return (mu == rhs.mu && nu == rhs.nu && z == rhs.z);
	}

	std::array<size_t, 6> nebrs;
	std::array<double, 6> faces_dist;
};

template <typename varType, typename PropsType>
double EllipticCell<varType, PropsType>::a;

template <class ellipCell>
inline typename ellipCell::Point getCartesian(double mu, double nu, double z)
{
	return{ ellipCell::a * cosh(mu) * cos(nu),
			ellipCell::a * sinh(mu) * sin(nu),
			z };
};

#endif /* ELLIPTICCELL_HPP_ */
