#ifndef ACIDELLFRACMODEL_HPP_
#define ACIDELLFRACMODEL_HPP_

#include <map>
#include "model/Acid/ellfrac/Properties.hpp"
#include "model/cells/AcidVariables.hpp"
#include "model/cells/NewEllipticCell.hpp"
#include "model/cells/HexCell.hpp"
#include "model/Basic1d/Basic1d.hpp"
#include "util/Interpolate.h"

namespace acidellfrac
{
	typedef JustAcid PoroVariable;
	typedef TapeJustAcid PoroTapeVariable;
	typedef JustAcid PoroVariable;
	typedef TapeJustAcid PoroTapeVariable;
	typedef VarFrac FracVariable;
	typedef TapeVarFrac FracTapeVariable;

	template <typename TVariable> using TCell = new_cell::EllipticCell<TVariable, Skeleton_Props>;
	typedef new_cell::EllipticCell<PoroVariable> PoroCell;
	typedef new_cell::EllipticCell<FracVariable> FracCell;
	typedef PoroCell::Type PoroType;
	typedef FracCell::Type FracType;
	typedef point::EllipticPoint3d Point;

	typedef CalciteReaction CurrentReaction;
	typedef CurrentReaction::REACTS REACTS;

	static const int var_poro_size = PoroVariable::size;
	static const int var_frac_size = FracVariable::size;
	typedef FracCell::Face Face;
	//typedef new_cell::AdjancedCellIdx CellIdx;
	typedef std::unordered_map<new_cell::AdjancedCellIdx,Face,new_cell::IdxHasher> FaceMap;
	const int NEBRS_NUM = new_cell::NEBRS_NUM;

	class AcidEllFrac
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class AbstractSolver;
		friend class AcidEllFracSolver; 
	public:
		static const int var_size = var_poro_size;
	protected:
		// Porous medium
		int skeletonsNum;
		std::vector<Skeleton_Props> props_sk;
		// Fracture
		FracProperties props_frac;
		// Fluids
		CurrentReaction reac;
		Water_Props props_w;
		Oil_Props props_o;
		Gas_Props props_g;
		
		// Grids
		int cellsNum_frac, cellsNum_poro, cellsNum_mu_frac, cellsNum_mu_poro, cellsNum_nu, cellsNum_z;
		double Volume_frac, Volume_poro;
		//std::vector<PoroGrid> poro_grids;
		std::vector<PoroCell> cells_poro;
		std::vector<FracCell> cells_frac;
		double re;
		FaceMap fmap_frac, fmap_poro, fmap_inter;
		
		// Temporary properties
		double ht, ht_min, ht_max;
		int periodsNum;
		double Q_sum, Pwf, c, height_perf;
		std::map<int, double> Qcell;
		std::vector<double> period, rate, pwf, cs;
		std::vector<bool> LeftBoundIsRate;
		bool leftBoundIsRate;
		bool rightBoundIsPres;
		// Snapshotter
		bool isWriteSnaps;
		Snapshotter<AcidEllFrac>* snapshotter;

		void buildFracGrid();
		void buildPoroGrid();
		void processGeometry();
		void setProps(Properties& props);
		void makeDimLess();
		void setInitialState();
		void setPerforated();
		void calculateTrans();
		// Service functions
		// Schemes
		PoroTapeVariable solvePoro(const PoroCell& cell);
		PoroTapeVariable solvePoroMid(const PoroCell& cell);
		PoroTapeVariable solvePoroLeft(const PoroCell& cell);
		PoroTapeVariable solvePoroRight(const PoroCell& cell);
		PoroTapeVariable solvePoroBorder(const PoroCell& cell);
		FracTapeVariable solveFrac(const FracCell& cell);
		FracTapeVariable solveFracIn(const FracCell& cell);
		FracTapeVariable solveFracMid(const FracCell& cell);
		FracTapeVariable solveFracBorder(const FracCell& cell);
		// Service calculations
		template <class TCell>
		inline int getUpwindIdx(const TCell& cell, const TCell& beta) const
		{
			if (cell.u_next.p < beta.u_next.p)
				return beta.num;
			else
				return cell.num;
		};
		template <>
		inline int getUpwindIdx(const FracCell& cell, const FracCell& beta) const
		{
			/*if (fabs(cell.c.mu - beta.c.mu) > cell.h.mu * EQUALITY_TOLERANCE)
				return (cell.c.mu > beta.c.mu) ? beta.num : cell.num;
			if (cell.u_next.p < beta.u_next.p)
				return beta.num;
			else
				return cell.num;*/
			return -1;
		};
		inline adouble getAverage(adouble p1, const double l1, adouble p2, const double l2) const
		{
			return (p1 * (adouble)l2 + p2 * (adouble)l1) / (adouble)(l1 + l2);
		};
		inline adouble getPoroTrans(const PoroCell& cell, const PoroTapeVariable& next, const char idx,
									const PoroCell& beta, const PoroTapeVariable& nebr) const
		{
			adouble k1, k2, S;
			const auto& props = props_sk.back();
			k1 = props.getPermCoseni(next.m, next.p);
			k2 = props.getPermCoseni(nebr.m, nebr.p);
			if (k1 == 0.0 && k2 == 0.0)
				return 0.0;
			S = fmap_poro.at({cell.num, beta.num}).S;
			return k1 * k2 * S / (k1 * cell.faces_dist[idx] + k2 * beta.faces_dist[cell.nebrs_idx[idx]]);
		};
		/*inline adouble getPoroAverage(adouble p1, const PoroCell& cell1, adouble p2, const PoroCell& cell2) const
		{
			return (p1 * (adouble)cell2.hx + p2 * (adouble)cell1.hx) / (adouble)(cell1.hx + cell2.hx);
		};
		inline adouble getFracAverage(const adouble p1, const FracCell& cell1, const adouble p2, const FracCell& cell2) const
		{
			if (cell1.x - cell2.x > EQUALITY_TOLERANCE)
			{
				//return (p1 * (adouble)cell2.hx + p2 * (adouble)cell1.hx) / (adouble)(cell1.hx + cell2.hx);
				if (cell1.x < cell2.x)
					return p1;
				else
					return p2;
			}
			else if (cell1.y - cell2.y > EQUALITY_TOLERANCE)
			{
				return (p1 * (adouble)cell2.hy + p2 * (adouble)cell1.hy) / (adouble)(cell1.hy + cell2.hy);
			}
			else
				return (p1 * (adouble)cell2.hz + p2 * (adouble)cell1.hz) / (adouble)(cell1.hz + cell2.hz);
		};*/
		inline adouble getReactionRate(const PoroTapeVariable& var, const Skeleton_Props& props) const
		{
			adouble m = props.getPoro(var.m, var.p);
			return var.sw * props_w.getDensity(var.p, var.xa, var.xw, var.xs) *
				(var.xa - props.xa_eqbm) *
				reac.getReactionRate(props.m_init, m) / reac.comps[REACTS::ACID].mol_weight;
		};
		inline const int getFracOut(const int idx) const
		{
			return int(idx / (cellsNum_mu_frac + 1)) * (cellsNum_mu_frac + 1) + cellsNum_mu_frac;
		};
		inline const int getFirstMuPoro(const int idx) const
		{
			const auto& out_cell = cells_frac[getFracOut(idx)];
			//assert(out_cell.type == FracType::FRAC_OUT);
			const int nebr_idx = int(idx / (cellsNum_mu_frac + 1)) * (cellsNum_mu_poro + 2);
			assert(out_cell.nebrs[1] < 0 || out_cell.nebrs[1] == nebr_idx);
			return nebr_idx;
		};
		inline adouble getFlowLeak(const FracCell& cell) const 
		{
			const int idx = getFirstMuPoro(cell.num);
			const auto& poro_cell = cells_poro[idx];
			assert(poro_cell.type == PoroType::WELL_LAT);
			const auto& poro_beta = cells_poro[idx + 1];
			const auto& next = x_poro[idx];
			const auto& nebr = x_poro[idx + 1];
			const auto& props = props_sk.back();
			return -props.getPermCoseni(next.m, next.p) * props_w.getKr(next.sw, next.m, &props) /
				props_w.getViscosity(next.p, next.xa, next.xw, next.xs) / next.m / next.sw *
				(nebr.p - next.p) / (poro_cell.faces_dist[1] + poro_beta.faces_dist[0]);
		};
		inline double getFracDistance(const int idx1, const int idx2) const
		{
			const auto& cell1 = cells_frac[idx1];
			const auto& cell2 = cells_frac[idx2];
			return point::distance(cell1.c, cell2.c, (cell1.c + cell2.c) / 2.0);
		};
		inline adouble getVelocity(const FracCell& cell, const PoroCell& beta) const
		{
			std::array<adouble, 2> res;
			const auto pt = (cell.c + beta.c) / 2.0;
			const double h = pt.getH();

			const auto& out_cell = beta;
			const double rat_cell = sinh(cell.c.mu) / sinh(out_cell.c.mu);
			const double rat_cur = sinh(pt.mu) / sinh(out_cell.c.mu);
			double alpha = -props_frac.w2_avg * props_frac.w2_avg / 2.0 / props_w.visc * (1.0 - rat_cell * rat_cell);
			adouble vel_y = getFlowLeak(cell) * (1.5 * rat_cur - 0.5 * rat_cur * rat_cur * rat_cur);

			res[0] = -(x_frac[cell.num].p - x_poro[beta.num].p) / point::distance(cell.c, beta.c, cell.c);
			const auto& pt_plus1 = cells_frac[cell.nebrs[3]].c;
			const auto& pt_plus2 = cells_poro[beta.nebrs[3]].c;
			const auto pt_plus = (pt_plus1 + pt_plus2) / 2.0;
			adouble p_plus = getAverage(x_frac[cell.nebrs[3]].p, point::distance(pt_plus1, pt_plus, pt_plus),
				x_poro[beta.nebrs[3]].p, point::distance(pt_plus2, pt_plus, pt_plus));
			const auto& pt_minus1 = cells_frac[cell.nebrs[2]].c;
			const auto& pt_minus2 = cells_poro[beta.nebrs[2]].c;
			const auto pt_minus = (pt_minus1 + pt_minus2) / 2.0;
			adouble p_minus = getAverage(x_frac[cell.nebrs[2]].p, point::distance(pt_minus1, pt_minus, pt_minus),
				x_poro[beta.nebrs[2]].p, point::distance(pt_minus2, pt_minus, pt_minus));
			res[1] = (p_plus - p_minus) / point::distance(pt_plus, pt_minus, (pt_plus + pt_minus) / 2.0);
			adouble vel_x = alpha * Point::a * (sinh(pt.mu) * cos(pt.nu) * res[0] - cosh(pt.mu) * sin(pt.nu) * res[1]) / h / h;
			return h * pt.getEllipticalVector(vel_x, vel_y)[0];
		};
		inline adouble getVelocity(const FracCell& cell, const FracCell& beta) const
		{
			std::array<adouble, 2> res;
			const auto pt = (cell.c + beta.c) / 2.0;
			const double h = pt.getH();

			const auto& out_cell = cells_poro[getFirstMuPoro(cell.num)];
			const double rat_cell = sinh(cell.c.mu) / sinh(out_cell.c.mu);
			const double rat_cur = sinh(pt.mu) / sinh(out_cell.c.mu);
			double alpha = -props_frac.w2_avg * props_frac.w2_avg / 2.0 / props_w.visc * (1.0 - rat_cell * rat_cell);
			adouble vel_y = getFlowLeak(cell) * (1.5 * rat_cur - 0.5 * rat_cur * rat_cur * rat_cur);

			if (fabs(cell.c.mu - beta.c.mu) > EQUALITY_TOLERANCE)
			{
				res[0] = (x_frac[cell.num].p - x_frac[beta.num].p) / getFracDistance(cell.num, beta.num);
				res[0] = cell.num > beta.num ? res[0] : -res[0];

				const Point& pt_plus1 = cells_frac[cell.nebrs[3]].c;
				const Point& pt_plus2 = cells_frac[cell.nebrs[3] + (beta.num - cell.num)].c;
				const Point pt_plus = (pt_plus1 + pt_plus2) / 2.0;
				adouble p_plus = getAverage(x_frac[cell.nebrs[3]].p,							point::distance(pt_plus1, pt_plus, pt_plus),
											x_frac[cell.nebrs[3] + (beta.num - cell.num)].p, point::distance(pt_plus2, pt_plus, pt_plus));
				const Point& pt_minus1 = cells_frac[cell.nebrs[2]].c;
				const Point& pt_minus2 = cells_frac[cell.nebrs[2] + (beta.num - cell.num)].c;
				const Point pt_minus = (pt_minus1 + pt_minus2) / 2.0;
				adouble p_minus = getAverage(x_frac[cell.nebrs[2]].p, point::distance(pt_minus1, pt_minus, pt_minus),
											x_frac[cell.nebrs[2] + (beta.num - cell.num)].p, point::distance(pt_minus2, pt_minus, pt_minus));
				res[1] = (p_plus - p_minus) / point::distance(pt_plus, pt_minus, (pt_plus + pt_minus) / 2.0);
				adouble vel_x = alpha * Point::a * (sinh(pt.mu) * cos(pt.nu) * res[0] - cosh(pt.mu) * sin(pt.nu) * res[1]) / h / h;
				return h * pt.getEllipticalVector(vel_x, vel_y)[0];
			}
			if (fabs(cell.c.nu - beta.c.nu) > EQUALITY_TOLERANCE)
			{
				res[1] = (x_frac[cell.num].p - x_frac[beta.num].p) / getFracDistance(cell.num, beta.num);
				res[1] = cell.num > beta.num ? -res[1] : res[1];

				const int idx = cell.nebrs[1] + (beta.num - cell.num) / (cellsNum_mu_frac + 1) * (cellsNum_mu_poro + 2);
				const Point& pt_plus1 = (cell.type != FracType::FRAC_OUT ? cells_frac[cell.nebrs[1]].c : cells_poro[cell.nebrs[1]].c);
				const Point& pt_plus2 = (cell.type != FracType::FRAC_OUT ? cells_frac[cell.nebrs[1] + (beta.num - cell.num)].c :
										cells_poro[idx].c);
				const Point pt_plus = (pt_plus1 + pt_plus2) / 2.0;
				const auto& x1 = (cell.type != FracType::FRAC_OUT ? x_frac[cell.nebrs[1]].p : x_poro[cell.nebrs[1]].p);
				const auto& x2 = (cell.type != FracType::FRAC_OUT ? x_frac[cell.nebrs[1] + (beta.num - cell.num)].p :
										x_poro[idx].p);
				adouble p_plus = getAverage(x1, point::distance(pt_plus1, pt_plus, pt_plus), x2, point::distance(pt_plus2, pt_plus, pt_plus));
				const Point& pt_minus1 = cells_frac[cell.nebrs[0]].c;
				const Point& pt_minus2 = cells_frac[cell.nebrs[0] + (beta.num - cell.num)].c;
				const Point pt_minus = (pt_minus1 + pt_minus2) / 2.0;
				adouble p_minus = getAverage(x_frac[cell.nebrs[0]].p, point::distance(pt_minus1, pt_minus, pt_minus),
					x_frac[cell.nebrs[0] + (beta.num - cell.num)].p, point::distance(pt_minus2, pt_minus, pt_minus));
				res[0] = (p_plus - p_minus) / point::distance(pt_plus, pt_minus, (pt_plus + pt_minus) / 2.0);
				adouble vel_x = alpha * Point::a * (sinh(pt.mu) * cos(pt.nu) * res[0] - cosh(pt.mu) * sin(pt.nu) * res[1]) / h / h;
				return h * pt.getEllipticalVector(vel_x, vel_y)[1];
			}
		};
		/*inline adouble getFlowLeakNew(const FracCell& cell)
		{
			const auto& grid = poro_grids[frac2poro[cells_frac[getRowOuter(cell.num)].num]];
			const auto& beta = cells_frac[grid.frac_nebr->num];
			assert(beta.type == FracType::FRAC_OUT);
			
			const auto& next = x_poro[grid.start_idx];
			const auto& nebr = x_frac[beta.num];
			
			return -props_frac.w2 * props_frac.w2 / 3.0 / props_w.visc * (next.p - nebr.p) / (grid.cells[0].x - beta.y);
		}
		inline adouble getQuadAppr(const std::array<adouble, 3> p, const std::array<double, 3> x, const double cur_x) const
		{
			const double den = (x[0] - x[1]) * (x[0] - x[2]) * (x[1] - x[2]);
			adouble a = (	p[2] * (x[0] - x[1]) +
							p[0] * (x[1] - x[2]) + 
							p[1] * (x[2] - x[0])				) / den;
			adouble b = (	p[2] * (x[1] * x[1] - x[0] * x[0]) + 
							p[0] * (x[2] * x[2] - x[1] * x[1]) + 
							p[1] * (x[0] * x[0] - x[2] * x[2])	) / den;
			adouble c = (	p[2] * x[0] * x[1] * (x[0] - x[1]) +
							x[2] * (p[0] * x[1] * (x[1] - x[2]) + p[1] * x[0] * (x[2] - x[0]))) / den;
			return a * cur_x * cur_x + b * cur_x + c;
		};
		inline void getPoroNeighborIdx(const int cur, int* const neighbor)
		{
			neighbor[0] = cur - 1;
			neighbor[1] = cur + 1;
		};
		inline double upwindIsCur(const PoroCell& cell, const PoroCell& beta)
		{
			if (cell.u_next.p < beta.u_next.p)
				return 0.0;
			else
				return 1.0;
		};*/
	public:
		// Dimensions
		double t_dim, R_dim, P_dim, T_dim, Q_dim, grav;

		AcidEllFrac();
		~AcidEllFrac();

		void load(Properties& props)
		{
			setProps(props);
			setSnapshotter("", this);
			buildFracGrid();
			buildPoroGrid();
			processGeometry();
			setPerforated();
			setInitialState();

			snapshotter->dump_all(0);
		};
		void snapshot_all(int i) { snapshotter->dump_all(i); }
		void setSnapshotter(std::string type, AcidEllFrac* model)
		{
			snapshotter = new VTKSnapshotter<AcidEllFrac>();
			snapshotter->setModel(model);
			isWriteSnaps = true;
		};

		void setPeriod(int period);

		adouble* h;
		PoroTapeVariable* x_poro;
		FracTapeVariable* x_frac;

		int getCellsNum() { return cellsNum_frac + cellsNum_poro; };
		double getRate(int cur) const;
	};
};

#endif /* ACIDELLFRACMODEL_HPP_ */