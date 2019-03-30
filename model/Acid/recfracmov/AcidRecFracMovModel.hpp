#ifndef ACIDRECFRACMOVMODEL_HPP_
#define ACIDRECFRACMOVMODEL_HPP_

#include <map>
#include "model/Acid/recfracmov/Properties.hpp"
#include "model/cells/AcidVariables.hpp"
#include "model/cells/HexCell.hpp"
#include "model/Basic1d/Basic1d.hpp"
#include "util/Interpolate.h"

namespace acidrecfracmov
{
	enum Regime {INJECTION, STOP};

	typedef JustAcid PoroVariable;
	typedef TapeJustAcid PoroTapeVariable;
	typedef JustAcid PoroVariable;
	typedef TapeJustAcid PoroTapeVariable;
	typedef VarFrac FracVariable;
	typedef TapeVarFrac FracTapeVariable;

	typedef HexCell<PoroVariable, Skeleton_Props> PoroCell;
	typedef HexCell<FracVariable> FracCell;
	typedef PoroCell::Type PoroType;
	typedef FracCell::Type FracType;

	typedef DolomiteReaction CurrentReaction;
	typedef CurrentReaction::REACTS REACTS;

	static const int var_poro_size = PoroVariable::size;
	static const int var_frac_size = FracVariable::size;
	//typedef FracCell::Face Face;
	//typedef std::unordered_map<new_cell::AdjancedCellIdx,Face,new_cell::IdxHasher> FaceMap;
	const int NEBRS_NUM = 6;

	class AcidRecFracMov
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class AbstractSolver;
		friend class AcidRecFracMovSolver; 
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
		int cellsNum_frac, cellsNum_poro, cellsNum_y_frac, cellsNum_y_poro, cellsNum_x, cellsNum_z;
		double Volume_frac, Volume_poro;
		//std::vector<PoroGrid> poro_grids;
		std::vector<PoroCell> cells_poro;
		std::vector<FracCell> cells_frac;
		double re;
		double dist_x;
		double max_vel_x, max_vel_y, max_vel_z;
		
		// Temporary properties
		double ht, ht_min, ht_max;
		int periodsNum;
		double Q_sum, Pwf, c, height_perf;
		std::map<int, double> Qcell;
		std::vector<double> period, rate, pwf, cs;
		std::vector<bool> LeftBoundIsRate;
		double injected_sol_volume, injected_acid_volume, max_sol_volume;
		bool leftBoundIsRate;
		bool rightBoundIsPres;
		// Snapshotter
		bool isWriteSnaps;
		Snapshotter<AcidRecFracMov>* snapshotter;
		std::vector<double> trans;
		std::vector<double> widths;
		double width;

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
		PoroTapeVariable solvePoro(const PoroCell& cell, const Regime reg);
		PoroTapeVariable solvePoroMid(const PoroCell& cell);
		PoroTapeVariable solvePoroLeft(const PoroCell& cell, const Regime reg);
		PoroTapeVariable solvePoroRight(const PoroCell& cell);
		PoroTapeVariable solvePoroBorder(const PoroCell& cell);
		FracTapeVariable solveFrac(const FracCell& cell, const Regime reg);
		FracTapeVariable solveFracIn(const FracCell& cell);
		FracTapeVariable solveFracMid(const FracCell& cell, const Regime reg);
		FracTapeVariable solveFracBorder(const FracCell& cell);
		// Service calculations
		const int getSkeletonId(const PoroCell& cell) const
		{
			// Random
			/*const int i_idx = int(cell.num / ((cellsNum_y_poro + 2) * (cellsNum_z + 2)));
			const int k_idx = int((cell.num - i_idx * (cellsNum_y_poro + 2) * (cellsNum_z + 2)) / (cellsNum_y_poro + 2));
			int i_idx_modi = i_idx - 1;
			if (i_idx_modi < 0)
				i_idx_modi = 0;
			if (i_idx_modi == cellsNum_x)
				i_idx_modi = cellsNum_x - 1;

			int k_idx_modi = k_idx - 1;
			if (k_idx_modi < 0)
				k_idx_modi = 0;
			if (k_idx_modi == cellsNum_z)
				k_idx_modi = cellsNum_z - 1;

			return i_idx_modi * cellsNum_z + k_idx_modi;*/

			// Example
			const int i_idx = int(cell.num / ((cellsNum_y_poro + 2) * (cellsNum_z + 2)));
			const int k_idx = int((cell.num - i_idx * (cellsNum_y_poro + 2) * (cellsNum_z + 2)) / (cellsNum_y_poro + 2));
			
			int k_idx_modi = k_idx - 1;
			if (k_idx_modi < 0)
				k_idx_modi = 0;
			if (k_idx_modi == cellsNum_z)
				k_idx_modi = cellsNum_z - 1;
			return k_idx_modi;
		}
		template <class TCell>
		inline int getUpwindIdx(const TCell& cell, const TCell& beta) const
		{
			if (cell.u_next.p < beta.u_next.p)
				return beta.num;
			else
				return cell.num;
		};
		inline const size_t getFracNebr(const size_t poro_idx) const
		{
			return poro_idx / (cellsNum_y_poro + 2) * (cellsNum_y_frac + 1) + cellsNum_y_frac;
		};
		inline const size_t getPoroNebr(const size_t frac_idx)
		{
			return (frac_idx - cellsNum_y_frac) / (cellsNum_y_frac + 1) * (cellsNum_y_poro + 2);
		};
		inline adouble getAverage(adouble p1, const double l1, adouble p2, const double l2) const
		{
			return (p1 * (adouble)l2 + p2 * (adouble)l1) / (adouble)(l1 + l2);
		};
		inline adouble getPoroTrans(const PoroCell& cell, const PoroTapeVariable& next,	const PoroCell& beta, const PoroTapeVariable& nebr) const
		{
			adouble k1, k2;
			k1 = cell.props->getPermCoseni(next.m, next.p);
			k2 = beta.props->getPermCoseni(nebr.m, nebr.p);
			if (k1 == 0.0 && k2 == 0.0)
				return 0.0;

			if (fabs(cell.x - beta.x) > EQUALITY_TOLERANCE)
				return 2.0 * k1 * k2 * cell.hy * cell.hz / (k1 * beta.hx + k2 * cell.hx);
			if (fabs(cell.y - beta.y) > EQUALITY_TOLERANCE)
				return 2.0 * k1 * k2 * cell.hx * cell.hz / (k1 * beta.hy + k2 * cell.hy);
			if (fabs(cell.z - beta.z) > EQUALITY_TOLERANCE)
				return 2.0 * k1 * k2 * cell.hx * cell.hy / (k1 * beta.hz + k2 * cell.hz);
		};
		inline adouble getReactionRate(const PoroTapeVariable& var, const Skeleton_Props& props) const
		{
			adouble tmp = (var.xa - props.xa_eqbm) * props_w.getDensity(var.p, var.xa, var.xw, var.xs) / reac.comps[REACTS::ACID].mol_weight;
			adouble isAboveEQ = (tmp.value() > 0.0) ? (adouble)true : (adouble)false;
			adouble tmp1;
			condassign(tmp1, isAboveEQ, pow(tmp, reac.alpha), (adouble)0.0);
			return /*var.sw */	tmp * reac.getReactionRate(props.m_init, props.m_max, var.m);
		};
        inline adouble getReactionRateOutput(const PoroVariable& var, const Skeleton_Props& props) const
        {
            adouble tmp = (var.xa - props.xa_eqbm) * props_w.getDensity(var.p, var.xa, var.xw, var.xs) / reac.comps[REACTS::ACID].mol_weight;
            adouble isAboveEQ = (tmp.value() > 0.0) ? (adouble)true : (adouble)false;
            adouble tmp1;
            condassign(tmp1, isAboveEQ, pow(tmp, reac.alpha), (adouble)0.0);
            return /*var.sw */	tmp * reac.getReactionRate(props.m_init, props.m_max, var.m);
        };
		inline double getDarmkoller(const PoroCell& cell, const PoroVariable& var, const Skeleton_Props& props) const
		{
			std::array<double, 3> vel;
			vel = getPoroWaterVelocity(cell);

			//adouble tmp = (var.xa - props.xa_eqbm) * props_w.getDensity(var.p, var.xa, var.xw, var.xs) / reac.comps[REACTS::ACID].mol_weight;
			//adouble isAboveEQ = (tmp.value() > 0.0) ? (adouble)true : (adouble)false;
			//adouble tmp1;
			//condassign(tmp1, isAboveEQ, pow(tmp, reac.alpha - 1.0), (adouble)0.0);
            double vel_mag = sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
            double tmp22 = /*var.sw */ reac.getSpecificReactionRate().value();
            if (vel_mag > EQUALITY_TOLERANCE && fabs(var.xa - props.xa_eqbm) > 1.E-3)
                return tmp22 / vel_mag;
            else
                return 0.0;
		};
		inline const int getFracOut(const int idx) const
		{
			return int(idx / (cellsNum_y_frac + 1)) * (cellsNum_y_frac + 1) + cellsNum_y_frac;
		};
		inline const int getFirstPoro(const int idx) const
		{;
			return int(idx / (cellsNum_y_frac + 1)) * (cellsNum_y_poro + 2);
		};
		inline adouble getFlowLeak(const FracCell& cell) const 
		{
			const int idx = getFirstPoro(cell.num);
			const auto& poro_cell = cells_poro[idx];
			assert(poro_cell.type == PoroType::WELL_LAT);
			const auto& poro_beta = cells_poro[idx + 1];
			const auto& next = x_poro[idx];
			const auto& nebr = x_poro[idx + 1];
			const auto& props = *poro_cell.props;
			double tmp1 = nebr.p.value();
			double tmp2 = next.p.value();
            double tmp3 = props.getPermCoseni(next.m, next.p).value();
            double tmp4 = props_w.getKr(next.sw, next.m, &props).value();
            double tmp5 = 1.0 / (poro_cell.hy + poro_beta.hy);
            double tmp6 = next.sw.value();
            double tmp7 = next.m.value();
            if(next.sw.value() != 0.0)
			    return -props.getPermCoseni(next.m, next.p) * props_w.getKr(next.sw, next.m, &props) /
    				props_w.getViscosity(next.p, next.xa, next.xw, next.xs) / next.m / next.sw *
	    			(nebr.p - next.p) * 2.0 / (poro_cell.hy + poro_beta.hy);
            else
                return -props.getPermCoseni(next.m, next.p) / props_w.getViscosity(next.p, next.xa, next.xw, next.xs) / 
                next.m * (nebr.p - next.p) * 2.0 / (poro_cell.hy + poro_beta.hy);
		};
        inline std::array<double, 3> getPoroWaterVelocity(const PoroCell& cell) const
        {
            int neighbor[NEBRS_NUM];
            getPoroNeighborIdx(cell.num, neighbor);
            const auto& next = cell.u_next;

            double transmissivity = -cell.props->getPermCoseni(next.m, next.p).value() * 
                props_w.getKr(next.sw, next.m, cell.props).value() / props_w.getViscosity(next.p, next.xa, next.xw, next.xs).value();
            double vel_x, vel_y, vel_z;
            vel_x = transmissivity * (cells_poro[neighbor[3]].u_next.p - cells_poro[neighbor[2]].u_next.p) / (cells_poro[neighbor[3]].x - cells_poro[neighbor[2]].x);
            if (cell.type != PoroType::WELL_LAT)
                vel_y = transmissivity * (cells_poro[neighbor[1]].u_next.p - cells_poro[neighbor[0]].u_next.p) / (cells_poro[neighbor[1]].y - cells_poro[neighbor[0]].y);
            else
                vel_y = getFlowLeak(cells_frac[getFracNebr(cell.num)]).value();
            vel_z = transmissivity * (cells_poro[neighbor[5]].u_next.p - cells_poro[neighbor[4]].u_next.p) / (cells_poro[neighbor[5]].z - cells_poro[neighbor[4]].z);
            return{ vel_x, vel_y, vel_z };
        };
		/*inline double getFracDistance(const int idx1, const int idx2) const
		{
			const auto& cell1 = cells_frac[idx1];
			const auto& cell2 = cells_frac[idx2];
			return point::distance(cell1.c, cell2.c, (cell1.c + cell2.c) / 2.0);
		};*/
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
		};*/
		inline void getPoroNeighborIdx(const int cur, int* const neighbor) const
		{
			neighbor[0] = cur - 1;
			neighbor[1] = cur + 1;
			neighbor[2] = cur - (cellsNum_y_poro + 2) * (cellsNum_z + 2);
			neighbor[3] = cur + (cellsNum_y_poro + 2) * (cellsNum_z + 2);
			neighbor[4] = cur - (cellsNum_y_poro + 2);
			neighbor[5] = cur + (cellsNum_y_poro + 2);
		};
		inline void getFracNeighborIdx(const int cur, int* const neighbor)
		{
			neighbor[0] = cur - 1;
			if(cells_frac[cur].type == FracType::FRAC_OUT)
				neighbor[1] = getPoroNebr(cur);
			else
				neighbor[1] = cur + 1;
			neighbor[2] = cur - (cellsNum_y_frac + 1) * (cellsNum_z + 2);
			neighbor[3] = cur + (cellsNum_y_frac + 1) * (cellsNum_z + 2);
			neighbor[4] = cur - (cellsNum_y_frac + 1);
			neighbor[5] = cur + (cellsNum_y_frac + 1);
		};
		inline double upwindIsCur(const PoroCell& cell, const PoroCell& beta)
		{
			if (cell.u_next.p < beta.u_next.p)
				return 0.0;
			else
				return 1.0;
		};
	public:
		// Dimensions
		double t_dim, R_dim, P_dim, T_dim, Q_dim, grav;

		AcidRecFracMov();
		~AcidRecFracMov();

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
		void setSnapshotter(std::string type, AcidRecFracMov* model)
		{
			snapshotter = new VTKSnapshotter<AcidRecFracMov>();
			snapshotter->setModel(model);
			isWriteSnaps = true;
		};

		void setPeriod(int period);

		adouble* h;
		PoroTapeVariable* x_poro;
		FracTapeVariable* x_frac;

		int getCellsNum() { return cellsNum_frac + cellsNum_poro; };
		double getRate(const int idx) const;
	};
};

#endif /* ACIDRECFRACMOVMODEL_HPP_ */