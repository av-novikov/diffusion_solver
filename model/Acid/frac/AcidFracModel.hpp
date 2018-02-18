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
	
	template <typename TVariable> using TCell = LinearCell<TVariable, Skeleton_Props>;
	typedef LinearCell<PoroVariable, Skeleton_Props> PoroCell;
	typedef HexCell<FracVariable> FracCell;
	typedef PoroCell::Type PoroType;
	typedef FracCell::Type FracType;

	struct PoroGrid
	{
		std::vector<PoroCell> cells;
		double Volume;
		const Skeleton_Props* props_sk;
		double hx, hz;
		int cellsNum;
		const FracCell* frac_nebr;
	};

	typedef CalciteReaction CurrentReaction;
	typedef CurrentReaction::REACTS REACTS;

	class AcidFrac
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class VTKSnapshotter;
		friend class AcidFracSolver; 
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
		int cellsNum, cellsNum_x, cellsNum_y, cellsNum_z;
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
		void build1dGrid(const FracCell& cell);
		void setProps(Properties& props);
		void makeDimLess();
		void setInitialState();
		void setPerforated();

		// Service functionss

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
	};
};

#endif /* ACIDFRACMODEL_HPP_ */