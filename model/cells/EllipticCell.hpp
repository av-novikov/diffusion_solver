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

	bool isGhost;

	inline std::array<double, 3> getCartesian() const
	{
		return { a * cosh(mu) * cos(nu), a * sinh(mu) * sin(nu), z };
	};

	EllipticCell() : isGhost(false) {};
	EllipticCell(int _num, double _mu, double _nu, double _z, double _hmu, double _hnu, double _hz) : AbstractCell<varType>(_num), mu(_mu), nu(_nu), z(_z), hmu(_hmu), hnu(_hnu), hz(_hz)
	{
		V = a * a * (sinh(mu) * sinh(mu) + sin(nu) * sin(nu)) * hz;
		if (V == 0.0)
			isGhost = true;
		else
			isGhost = false;
	};
	~EllipticCell() {};
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
