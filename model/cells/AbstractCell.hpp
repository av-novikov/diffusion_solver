#ifndef ABSTRACTCELL_HPP_
#define ABSTRACTCELL_HPP_

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <array>

struct EmptyStruct
{};

template <int N> using TPoint = std::array<double, N>;

template <typename varType>
class AbstractCell
{
public: 
	enum class Type : int
	{MIDDLE, MIDDLE_SIDE, RIGHT, TOP, BOTTOM, WELL_LAT, WELL_TOP, WELL_BOT, WELL_SIDE, FRAC_IN, FRAC_BORDER, NOTYPE};
	public:		
		const int num;
		static const int varNum = varType::size;
		
		double V;

		const Type type;
		
		varType u_next;
		varType u_iter;
		varType u_prev;
		
		AbstractCell() : V(0.0), num(-1), type(Type::NOTYPE) {};
		AbstractCell(const int _num) : V(0.0), num(_num), type(Type::NOTYPE) {};
		AbstractCell(const int _num, const Type _type) : num(_num), type(_type) {};
		~AbstractCell() {};

		AbstractCell(const AbstractCell& cell) : num(cell.num), type(cell.type)
		{
			V = cell.V;
			u_next = cell.u_next;	u_iter = cell.u_iter;	u_prev = cell.u_prev;
		};
		AbstractCell& operator=(const AbstractCell& cell)
		{
			V = cell.V;
			u_next = cell.u_next;	u_iter = cell.u_iter;	u_prev = cell.u_prev;
			return *this;
		};
};

#endif /* ABSTRACTCELL_HPP_ */
