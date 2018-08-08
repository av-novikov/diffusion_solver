#ifndef NEWELLIPTICCELL_HPP_
#define NEWELLIPTICCELL_HPP_

#include <utility>
#include <unordered_map>
#include <functional>
#include <algorithm>

#include "model/cells/Point.hpp"
#include "model/cells/AbstractCell.hpp"

namespace new_cell
{
	template <class TPoint>
	struct Face
	{
		//TPoint c;
		double S;
	};
	struct AdjancedCellIdx
	{
		int first;
		int second;
		bool operator==(const AdjancedCellIdx& rhs) const
		{
			return (first == rhs.first && second == rhs.second) || (second == rhs.first && first == rhs.second);
		};
	};
	struct IdxHasher
	{
		static const int SIZE = 65536;
		int operator() (const AdjancedCellIdx& idx) const
		{
			//return hash<size_t>()(idx.first) + hash<size_t>()(idx.second);
			return SIZE * std::min(idx.first, idx.second) + std::max(idx.first, idx.second);
		};
	};

	const int NEBRS_NUM = 6;

	template <typename varType, typename PropsType = EmptyStruct>
	class EllipticCell : public AbstractCell<varType>
	{
	public:
		typedef point::EllipticPoint3d Point;
		typedef Face<Point> Face;

		//static double a;
		Point c, h;
		std::array<int, NEBRS_NUM> nebrs;
		std::array<char, NEBRS_NUM> nebrs_idx;
		std::array<double, NEBRS_NUM> faces_dist;
		bool isUsed;

		EllipticCell() : isUsed(true) { nebrs.fill(-1);	nebrs_idx.fill(-1);	};
		EllipticCell(int _num, double _mu, double _nu, double _z, double _hmu, double _hnu, double _hz, const Type _type) : AbstractCell<varType>(_num, _type), 
			c(_mu, _nu, _z), h(_hmu, _hnu, _hz)
		{
			const auto g = c.getMetric();
			V = sqrt(g[0] * g[1] * g[2]) * h.mu * h.nu * h.z;
			isUsed = true;
			nebrs.fill(-1);	nebrs_idx.fill(-1);
		};
		EllipticCell(int _num, const Point& _c, const Point& _h, const Type _type) : AbstractCell<varType>(_num, _type), c(_c), h(_h)
		{
			const auto g = c.getMetric();
			V = sqrt(g[0] * g[1] * g[2]) * h.mu * h.nu * h.z;
			isUsed = true;
			nebrs.fill(-1);	nebrs_idx.fill(-1);
		};
		~EllipticCell() {};

		bool operator==(const EllipticCell& rhs) const
		{
			return (mu == rhs.mu && nu == rhs.nu && z == rhs.z);
		}
	};
};

#endif /* NEWELLIPTICCELL_HPP_ */
