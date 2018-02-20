#ifndef ACIDFRACMODEL_HPP_
#define ACIDFRACMODEL_HPP_

#include <map>
#include "model/Acid/frac/Properties.hpp"
#include "model/cells/AcidVariables.hpp"
#include "model/cells/LinearCell.hpp"
#include "model/cells/HexCell.hpp"
#include "model/Basic1d/Basic1d.hpp"
#include "util/Interpolate.h"

namespace acidfrac
{
	typedef JustAcid PoroVariable;
	typedef TapeJustAcid PoroTapeVariable;
	typedef JustAcid PoroVariable;
	typedef TapeJustAcid PoroTapeVariable;
	typedef VarFrac FracVariable;
	typedef TapeVarFrac FracTapeVariable;

	struct PoroGrid;
	template <typename TVariable> using TCell = LinearCell<TVariable, Skeleton_Props>;
	typedef LinearCell<PoroVariable, PoroGrid> PoroCell;
	typedef HexCell<FracVariable> FracCell;
	typedef PoroCell::Type PoroType;
	typedef FracCell::Type FracType;
	struct PoroGrid
	{
		int id, start_idx;
		std::vector<PoroCell> cells;
		double Volume;
		const Skeleton_Props* props_sk;
		double hx, hz;
		int cellsNum;
		const FracCell* frac_nebr;
	};

	typedef CalciteReaction CurrentReaction;
	typedef CurrentReaction::REACTS REACTS;

	static const int var_poro_size = PoroVariable::size;
	static const int var_frac_size = FracVariable::size;

	class AcidFrac
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class AbstractSolver;
		friend class AcidFracSolver; 
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
		std::vector<PoroGrid> poro_grids;
		std::vector<FracCell> cells_frac;
		int cellsNum, cellsNum_x, cellsNum_y, cellsNum_z, cellsPoroNum;
		double Volume;
		std::vector<int> cellsNum_y_1d;
		std::vector<double> xe;
		std::map<FracCell*, PoroGrid*> frac2poro;
		// Temporary properties
		double ht, ht_min, ht_max;
		int periodsNum;
		double Q_sum, Pwf, c, height_perf;
		std::map<int, double> Qcell;
		std::vector<double> period, rate, pwf, cs;
		bool leftBoundIsRate, rightBoundIsPres;
		// Scheme
		/*TapeVariable* var;
		adouble* h;
		double* x;
		double* y;
		double** jac;*/
		// Snapshotter
		bool isWriteSnaps;
		Snapshotter<AcidFrac>* snapshotter;

		size_t getInitId2OutCell(const FracCell& cell);
		void buildGrid();
		void build1dGrid(const FracCell& cell, const int grid_id);
		void setProps(Properties& props);
		void makeDimLess();
		void setInitialState();
		void setPerforated();
		// Schemes
		PoroTapeVariable solvePoro(const PoroCell& cell);
		PoroTapeVariable solvePoroMid(const PoroCell& cell);
		PoroTapeVariable solvePoroLeft(const PoroCell& cell);
		PoroTapeVariable solvePoroRight(const PoroCell& cell);
		FracTapeVariable solveFrac(const FracCell& cell);
		FracTapeVariable solveFracIn(const FracCell& cell);
		FracTapeVariable solveFracMid(const FracCell& cell);
		FracTapeVariable solveFracBorder(const FracCell& cell);
		FracTapeVariable solveFracOut(const FracCell& cell);
		// Service calculations
		inline adouble getPoroTrans(const PoroCell& cell, adouble m_cell, const PoroCell& beta, adouble m_beta) const
		{
			adouble k1, k2, S;
			const auto& props = *cell.props->props_sk;
			k1 = props.getPermCoseni(m_cell);
			k2 = props.getPermCoseni(m_beta);
			if (k1 == 0.0 && k2 == 0.0)
				return 0.0;
			S = cell.props->hz * cell.props->hx;
			return 2.0 * k1 * k2 * S / (k1 * beta.hx + k2 * cell.hx);
		};
		inline adouble getAverage(adouble p1, const PoroCell& cell1, adouble p2, const PoroCell& cell2) const
		{
			return (p1 * (adouble)cell2.hx + p2 * (adouble)cell1.hx) / (adouble)(cell1.hx + cell2.hx);
		};
		inline adouble getReactionRate(const PoroTapeVariable& var, const Skeleton_Props& props) const
		{
			return var.sw * props_w.getDensity(var.p, var.xa, var.xw, var.xs) *
				(var.xa - props.xa_eqbm) *
				reac.getReactionRate(props.m_init, var.m) / reac.comps[REACTS::ACID].mol_weight;
		};
		// Service functions
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
		};
		inline int getUpwindIdx(const PoroCell& cell, const PoroCell& beta) const
		{
			if (cell.u_next.p < beta.u_next.p)
				return beta.num;
			else
				return cell.num;
		};
	public:
		// Dimensions
		double t_dim, R_dim, P_dim, T_dim, Q_dim, grav;

		AcidFrac();
		~AcidFrac();

		void load(Properties& props)
		{
			setProps(props);
			setSnapshotter("", this);
			buildGrid();
			setPerforated();
			setInitialState();

			snapshotter->dump_all(0);
		};
		void snapshot_all(int i) { snapshotter->dump_all(i); }
		void setSnapshotter(std::string type, AcidFrac* model)
		{
			snapshotter = new VTKSnapshotter<AcidFrac>();
			snapshotter->setModel(model);
			isWriteSnaps = true;
		};

		void setPeriod(int period);

		adouble* h;
		PoroTapeVariable* x_poro;
		FracTapeVariable* x_frac;

		int getCellsNum() { return cellsNum + cellsPoroNum; };
	};
};

#endif /* ACIDFRACMODEL_HPP_ */