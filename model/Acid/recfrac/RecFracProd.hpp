#ifndef RECFRACPROD_HPP_
#define RECFRACPROD_HPP_

#include <map>
#include "model/Acid/recfrac/Properties.hpp"
#include "model/cells/AcidVariables.hpp"
#include "model/cells/HexCell.hpp"
#include "model/Basic1d/Basic1d.hpp"
#include "util/Interpolate.h"
#include "model/Acid/recfrac/AcidRecFracModel.hpp"

namespace acidrecfrac_prod
{
    typedef acidrecfrac::Properties Properties;
    typedef acidrecfrac::Water_Props Water_Props;
    typedef acidrecfrac::Oil_Props Oil_Props;
    typedef acidrecfrac::Skeleton_Props Skeleton_Props;
    typedef VarOilWater Variable;
    typedef TapeOilWater TapeVariable;
    typedef acidrecfrac::PoroCell PoroCell;
    typedef HexCell<Variable, Skeleton_Props> Cell;
    typedef Cell::Type Type;
    const int NEBRS_NUM = acidrecfrac::NEBRS_NUM;

    static const int var_size = Variable::size;
	class RecFracProd
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class AbstractSolver;
		friend class RecFracProdSolver; 
	public:
        static const int var_size = var_size;
	protected:
		std::vector<Skeleton_Props> props_sk;
		// Fluids
		Water_Props props_w;
		Oil_Props props_o;
		
		// Grids
        double Volume;
		std::vector<Cell> cells;
        double x_size, y_size, z_size;
        double prev_x_size, prev_y_size, prev_z_size;
        int cellsNum, cellsNum_x, cellsNum_y, cellsNum_z;
        int prev_cellsNum_x, prev_cellsNum_y, prev_cellsNum_z;

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
		Snapshotter<RecFracProd>* snapshotter;

		void buildGrid(std::vector<PoroCell>& cells_poro);
		void processGeometry();
		void setProps(Properties& props);
		void makeDimLess();
		void setInitialState(const std::vector<PoroCell>& cells_poro);
		void setPerforated();
		// Service calculations
		template <class TCell>
		inline int getUpwindIdx(const TCell& cell, const TCell& beta) const
		{
			if (cell.u_next.p < beta.u_next.p)
				return beta.num;
			else
				return cell.num;
		};
		inline adouble getAverage(adouble p1, const double l1, adouble p2, const double l2) const
		{
			return (p1 * (adouble)l2 + p2 * (adouble)l1) / (adouble)(l1 + l2);
		};
        inline adouble getPoroTrans(const Cell& cell, const Cell& beta) const
		{
			adouble k1, k2;
			k1 = cell.props->getPermCoseni(cell.u_next.m, cell.u_next.p);
			k2 = beta.props->getPermCoseni(cell.u_next.m, cell.u_next.p);
			if (k1 == 0.0 && k2 == 0.0)
				return 0.0;

			if (fabs(cell.x - beta.x) > EQUALITY_TOLERANCE)
				return 2.0 * k1 * k2 * cell.hy * cell.hz / (k1 * beta.hx + k2 * cell.hx);
			if (fabs(cell.y - beta.y) > EQUALITY_TOLERANCE)
				return 2.0 * k1 * k2 * cell.hx * cell.hz / (k1 * beta.hy + k2 * cell.hy);
			if (fabs(cell.z - beta.z) > EQUALITY_TOLERANCE)
				return 2.0 * k1 * k2 * cell.hx * cell.hy / (k1 * beta.hz + k2 * cell.hz);
            return 0.0;
		};
		inline double upwindIsCur(const Cell& cell, const Cell& beta)
		{
			if (cell.u_next.p < beta.u_next.p)
				return 0.0;
			else
				return 1.0;
		};
        const int getSkeletonId(const Cell& cell) const
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
            const int i_idx = int(cell.num / ((cellsNum_y + 2) * (cellsNum_z + 2)));
            const int k_idx = int((cell.num - i_idx * (cellsNum_y + 2) * (cellsNum_z + 2)) / (cellsNum_y + 2));

            int k_idx_modi = k_idx - 1;
            if (k_idx_modi < 0)
                k_idx_modi = 0;
            if (k_idx_modi == cellsNum_z)
                k_idx_modi = cellsNum_z - 1;
            return k_idx_modi;
        }
        inline void getPoroNeighborIdx(const int cur, int* const neighbor) const
        {
            neighbor[0] = cur - 1;
            neighbor[1] = cur + 1;
            neighbor[2] = cur - (cellsNum_y + 2) * (cellsNum_z + 2);
            neighbor[3] = cur + (cellsNum_y + 2) * (cellsNum_z + 2);
            neighbor[4] = cur - (cellsNum_y + 2);
            neighbor[5] = cur + (cellsNum_y + 2);
        };

        TapeVariable solvePoroMid(const Cell& cell);
        TapeVariable solvePoroLeft(const Cell& cell);
        TapeVariable solvePoroRight(const Cell& cell);
        TapeVariable solvePoroBorder(const Cell& cell);
        TapeVariable solveFrac(const Cell& cell);
	public:
		// Dimensions
		double t_dim, R_dim, P_dim, T_dim, Q_dim, grav;

        RecFracProd();
		~RecFracProd();

		void load(Properties& props, std::vector<PoroCell>& cells_poro)
		{
			setProps(props);
			setSnapshotter("", this);
			buildGrid(cells_poro);
			processGeometry();
			setPerforated();
			setInitialState(cells_poro);
		};
		void snapshot_all(int i) { snapshotter->dump_all(i); }
		void setSnapshotter(std::string type, RecFracProd* model)
		{
			snapshotter = new VTKSnapshotter<RecFracProd>();
			snapshotter->setModel(model);
			isWriteSnaps = true;
		};
		void setPeriod(int period);

		adouble* h;
		TapeVariable* x;

		int getCellsNum() { return cellsNum; };
		double getRate(const Cell& cell) const;
	};
};

#endif /* RECFRACPROD_HPP_ */