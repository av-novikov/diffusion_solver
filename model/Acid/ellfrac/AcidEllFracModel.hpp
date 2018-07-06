#ifndef ACIDELLFRACMODEL_HPP_
#define ACIDELLFRACMODEL_HPP_

#include <map>
#include "model/Acid/ellfrac/Properties.hpp"
#include "model/cells/AcidVariables.hpp"
#include "model/cells/EllipticCell.hpp"
#include "model/cells/HexCell.hpp"
#include "model/Basic1d/Basic1d.hpp"
#include "util/Interpolate.h"

namespace acidellfrac
{
	class AcidEllFracSolver
	{
	public:
		void start() {};
	};

	typedef JustAcid PoroVariable;
	typedef TapeJustAcid PoroTapeVariable;
	typedef JustAcid PoroVariable;
	typedef TapeJustAcid PoroTapeVariable;
	typedef VarFrac FracVariable;
	typedef TapeVarFrac FracTapeVariable;

	//struct PoroGrid;
	/*struct PoroGrid
	{
	int id, start_idx;
	std::vector<PoroCell> cells;
	double Volume;
	const Skeleton_Props* props_sk;
	double hx, hz;
	int cellsNum;
	const FracCell* frac_nebr;
	double width, trans;
	};*/
	template <typename TVariable> using TCell = EllipticCell<TVariable, Skeleton_Props>;
	typedef EllipticCell<PoroVariable> PoroCell;
	typedef EllipticCell<FracVariable> FracCell;
	typedef PoroCell::Type PoroType;
	typedef FracCell::Type FracType;
	typedef FracCell::Point Point;

	typedef CalciteReaction CurrentReaction;
	typedef CurrentReaction::REACTS REACTS;

	static const int var_poro_size = PoroVariable::size;
	static const int var_frac_size = FracVariable::size;

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
		size_t cellsNum_frac, cellsNum_poro, cellsNum_mu_frac, cellsNum_mu_poro, cellsNum_nu, cellsNum_z;
		double Volume_frac, Volume_poro;
		//std::vector<PoroGrid> poro_grids;
		std::vector<PoroCell> cells_poro;
		std::vector<FracCell> cells_frac;
		double re;
		std::map<int, int> frac2poro;
		std::map<int, int> border_nebrs;
		
		// Temporary properties
		double ht, ht_min, ht_max;
		int periodsNum;
		double Q_sum, Pwf, c, height_perf;
		std::map<int, double> Qcell;
		std::vector<double> period, rate, pwf, cs;
		std::vector<bool> LeftBoundIsRate;
		bool leftBoundIsRate;
		bool rightBoundIsPres;
		// Scheme
		/*TapeVariable* var;
		adouble* h;
		double* x;
		double* y;
		double** jac;*/
		// Snapshotter
		bool isWriteSnaps;
		Snapshotter<AcidEllFrac>* snapshotter;

		size_t getInitId2OutCell(const FracCell& cell);
		void buildGrid();
		void buildGridUniform();
		void buildPoroGrid();
		void setProps(Properties& props);
		void makeDimLess();
		void setInitialState();
		void setPerforated();
		void calculateTrans();
		// Schemes
		/*PoroTapeVariable solvePoro(const PoroCell& cell);
		PoroTapeVariable solvePoroMid(const PoroCell& cell);
		PoroTapeVariable solvePoroLeft(const PoroCell& cell);
		PoroTapeVariable solvePoroRight(const PoroCell& cell);
		FracTapeVariable solveFrac(const FracCell& cell);
		FracTapeVariable solveFracIn(const FracCell& cell);
		FracTapeVariable solveFracMid(const FracCell& cell);
		FracTapeVariable solveFracBorder(const FracCell& cell);
		FracTapeVariable solveFracOut(const FracCell& cell);
		// Service calculations
		inline adouble getPoroTrans(const PoroCell& cell, const PoroTapeVariable& next, const PoroCell& beta, const PoroTapeVariable& nebr) const
		{
			/*adouble k1, k2, S;
			const auto& props = *cell.props->props_sk;
			k1 = props.getPermCoseni(next.m, next.p);
			k2 = props.getPermCoseni(nebr.m, nebr.p);
			if (k1 == 0.0 && k2 == 0.0)
				return 0.0;
			S = cell.props->hz * cell.props->hx;
			return 2.0 * k1 * k2 * S / (k1 * beta.hx + k2 * cell.hx);
		};*/
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
		/*inline adouble getFlowLeak(const FracCell& cell)
		{
			const auto& grid = poro_grids[frac2poro[cells_frac[getRowOuter(cell.num)].num]];
			const auto& next = x_poro[grid.start_idx];
			const auto& nebr = x_poro[grid.start_idx + 1];
			return -grid.props_sk->getPermCoseni(next.m, next.p) * props_w.getKr(next.sw, next.m, grid.props_sk) /
				props_w.getViscosity(next.p, next.xa, next.xw, next.xs) / next.m / next.sw *
				(nebr.p - next.p) / (grid.cells[1].x - grid.cells[0].x);
		}
		inline adouble getFlowLeakNew(const FracCell& cell)
		{
			const auto& grid = poro_grids[frac2poro[cells_frac[getRowOuter(cell.num)].num]];
			const auto& beta = cells_frac[grid.frac_nebr->num];
			assert(beta.type == FracType::FRAC_OUT);
			
			const auto& next = x_poro[grid.start_idx];
			const auto& nebr = x_frac[beta.num];
			
			return -props_frac.w2 * props_frac.w2 / 3.0 / props_w.visc * (next.p - nebr.p) / (grid.cells[0].x - beta.y);
		}*/
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
		// Service functions
		/*inline void getPoroNeighborIdx(const int cur, int* const neighbor)
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
		};
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
			if (fabs(cell.y - beta.y) > cell.hy * EQUALITY_TOLERANCE)
				return (cell.y > beta.y) ? beta.num : cell.num;
			if (cell.u_next.p < beta.u_next.p)
				return beta.num;
			else
				return cell.num;
		};
		inline const int getRowOuter(const int idx) const
		{
			const int outer_idx = int(idx / (cellsNum_y + 1)) * (cellsNum_y + 1) + cellsNum_y;
			assert(cells_frac[outer_idx].type == FracType::FRAC_OUT);
			return outer_idx;
		};
		inline void getNeighborIdx(int cur, int* const neighbor)
		{
			neighbor[0] = cur - (cellsNum_y + 1) * (cellsNum_z + 2);
			neighbor[1] = cur + (cellsNum_y + 1) * (cellsNum_z + 2);
			neighbor[2] = cur - 1;
			neighbor[3] = cur + 1;
			neighbor[4] = cur - cellsNum_y - 1;
			neighbor[5] = cur + cellsNum_y + 1;
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
			buildGridUniform();
			//buildGrid();
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