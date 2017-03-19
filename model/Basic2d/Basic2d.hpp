#ifndef BASIC2D_HPP_
#define BASIC2D_HPP_

#include "model/cells/Variables.hpp"
#include "model/GasOil_RZ/Properties.hpp"
#include "model/AbstractModel.hpp"

#include <vector>
#include <map>
#include <cassert>

namespace basic2d
{
	static const int R_AXIS = 0;
	static const int Z_AXIS = 1;

	static const int stencil = 5;
	static const int Lstencil = 3;
	static const int Rstencil = 3;
	static const int Vstencil = 2;

	template <typename varType, typename propsType, typename sk_propsType,
	template <typename varType> class cellType, class modelType>
	class Basic2d : public AbstractModel<varType, propsType, cellType, modelType>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class AbstractMethod;
		template<typename> friend class Basic2dSolver;
	public:
		typedef cellType<varType> Cell;
		typename typedef Cell::Type Type;

		static const int var_size;
	protected:

		std::vector<sk_propsType> props_sk;

		// Continuum properties
		int skeletonsNum;
		// Number of cells in radial direction
		int cellsNum_r;
		// Number of cells in vertical direction
		int cellsNum_z;

		void setBasicProps(propsType& props)
		{
			leftBoundIsRate = props.leftBoundIsRate;
			rightBoundIsPres = props.rightBoundIsPres;

			// Setting grid properties
			r_w = props.r_w;
			r_e = props.r_e;
			cellsNum_r = props.cellsNum_r;
			cellsNum_z = props.cellsNum_z;
			cellsNum = (cellsNum_r + 2) * (cellsNum_z + 2);

			// Setting skeleton properties
			perfIntervals = props.perfIntervals;

			skeletonsNum = props.props_sk.size();
			checkSkeletons(props.props_sk);
			for (int j = 0; j < skeletonsNum; j++)
			{
				props_sk[j].perm_r = MilliDarcyToM2(props_sk[j].perm_r);
				props_sk[j].perm_z = MilliDarcyToM2(props_sk[j].perm_z);
			}

			periodsNum = props.timePeriods.size();
			for (int i = 0; i < periodsNum; i++)
			{
				period.push_back(props.timePeriods[i]);
				if (leftBoundIsRate)
					rate.push_back(props.rates[i] / 86400.0);
				else
					pwf.push_back(props.pwf[i]);
				for (int j = 0; j < skeletonsNum; j++)
				{
					if (props_sk[j].radiuses_eff[i] > props.r_w)
						props_sk[j].perms_eff.push_back(MilliDarcyToM2(props.props_sk[j].perm_r * log(props.props_sk[j].radiuses_eff[i] / props.r_w) / (log(props.props_sk[j].radiuses_eff[i] / props.r_w) + props.props_sk[j].skins[i])));
					else
						props_sk[j].perms_eff.push_back(MilliDarcyToM2(props.props_sk[j].perm_r));
				}
			}

			// Temporal properties
			ht = props.ht;
			ht_min = props.ht_min;
			ht_max = props.ht_max;

			alpha = props.alpha;
			depth_point = props.depth_point;
		};
		void makeBasicDimLess()
		{
			// Main units
			R_dim = r_w;
			t_dim = 3600.0;
			P_dim = props_sk[0].p_init;

			// Temporal properties
			ht /= t_dim;
			ht_min /= t_dim;
			ht_max /= t_dim;

			// Grid properties
			r_w /= R_dim;
			r_e /= R_dim;

			// Skeleton properties
			for (int i = 0; i < skeletonsNum; i++)
			{
				props_sk[i].perm_r /= (R_dim * R_dim);
				props_sk[i].perm_z /= (R_dim * R_dim);
				props_sk[i].d_pore_r /= R_dim;
				props_sk[i].d_pore_z /= R_dim;

				props_sk[i].beta /= (1.0 / P_dim);
				props_sk[i].dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
				props_sk[i].h1 = (props_sk[i].h1 - depth_point) / R_dim;
				props_sk[i].h2 = (props_sk[i].h2 - depth_point) / R_dim;
				props_sk[i].height /= R_dim;
				props_sk[i].p_init /= P_dim;
				props_sk[i].p_out /= P_dim;
				props_sk[i].p_ref /= P_dim;

				for (int j = 0; j < periodsNum; j++)
				{
					props_sk[i].perms_eff[j] /= (R_dim * R_dim);
					props_sk[i].radiuses_eff[j] /= R_dim;
				}
			}

			Q_dim = R_dim * R_dim * R_dim / t_dim;
			for (int i = 0; i < periodsNum; i++)
			{
				period[i] /= t_dim;
				if (leftBoundIsRate)
					rate[i] /= Q_dim;
				else
					pwf[i] /= P_dim;
			}

			// Rest properties
			alpha /= t_dim;
			//depth_point = 0.0;
		};
		void checkSkeletons(const vector<sk_propsType>& props)
		{
			vector<sk_propsType>::const_iterator it = props.begin();
			double tmp;
			int indxs = 0;

			assert(it->h2 - it->h1 == it->height);
			indxs += it->cellsNum_z;
			tmp = it->h2;
			++it;

			while (it != props.end())
			{
				assert(it->h1 == tmp);
				assert(it->h2 - it->h1 == it->height);
				indxs += it->cellsNum_z;
				tmp = it->h2;
				++it;
			}
			assert(indxs == cellsNum_z);
		}
		void buildGridLog()
		{
			cells.reserve(cellsNum);

			Volume = 0.0;
			int counter = 0;
			int skel_idx = 0, cells_z = 0;

			double r_prev = r_w;
			double logMax = log(r_e / r_w);
			double logStep = logMax / (double)cellsNum_r;

			double hz = 0.0;//(props_sk[skel_idx].h2 - props_sk[skel_idx].h1) / (double)(props_sk[skel_idx].cellsNum_z);
			double cm_z = props_sk[skel_idx].h1;
			double hr = r_prev * (exp(logStep) - 1.0);
			double cm_r = r_w;

			// Left border
			cells.push_back(Cell(counter++, cm_r, cm_z, 0.0, 0.0, Type::WELL_LAT));
			for (int i = 0; i < cellsNum_z; i++)
			{
				hz = (props_sk[skel_idx].h2 - props_sk[skel_idx].h1) / (double)(props_sk[skel_idx].cellsNum_z);
				cm_z += (cells[cells.size() - 1].hz + hz) / 2.0;

				cells.push_back(Cell(counter++, cm_r, cm_z, 0.0, hz, Type::WELL_LAT));
				cells_z++;

				if (cells_z >= props_sk[skel_idx].cellsNum_z)
				{
					cells_z = 0;
					skel_idx++;
				}
			}
			cells.push_back(Cell(counter++, cm_r, cm_z + hz / 2.0, 0.0, 0.0, Type::WELL_LAT));

			// Middle cells
			for (int j = 0; j < cellsNum_r; j++)
			{
				skel_idx = 0;	cells_z = 0;
				cm_z = props_sk[0].h1;
				cm_r = r_prev * (exp(logStep) + 1.0) / 2.0;
				hr = r_prev * (exp(logStep) - 1.0);

				cells.push_back(Cell(counter++, cm_r, cm_z, hr, 0.0, Type::TOP));
				for (int i = 0; i < cellsNum_z; i++)
				{
					hz = (props_sk[skel_idx].h2 - props_sk[skel_idx].h1) / (double)(props_sk[skel_idx].cellsNum_z);
					cm_z += (cells[cells.size() - 1].hz + hz) / 2.0;

					cells.push_back(Cell(counter++, cm_r, cm_z, hr, hz, Type::MIDDLE));
					Volume += cells[cells.size() - 1].V;
					cells_z++;

					if (cells_z >= props_sk[skel_idx].cellsNum_z)
					{
						cells_z = 0;
						skel_idx++;
					}
				}
				cells.push_back(Cell(counter++, cm_r, cm_z + hz / 2.0, hr, 0.0, Type::BOTTOM));

				r_prev = r_prev * exp(logStep);
			}

			// Right border
			cm_z = props_sk[0].h1;
			cm_r = r_e;

			cells.push_back(Cell(counter++, cm_r, cm_z, 0.0, 0.0, Type::RIGHT));
			skel_idx = 0;	cells_z = 0;
			for (int i = 0; i < cellsNum_z; i++)
			{
				hz = (props_sk[skel_idx].h2 - props_sk[skel_idx].h1) / (double)(props_sk[skel_idx].cellsNum_z);
				cm_z += (cells[cells.size() - 1].hz + hz) / 2.0;

				cells.push_back(Cell(counter++, cm_r, cm_z, 0.0, hz, Type::RIGHT));
				cells_z++;

				if (cells_z >= props_sk[skel_idx].cellsNum_z)
				{
					cells_z = 0;
					skel_idx++;
				}
			}
			cells.push_back(Cell(counter++, cm_r, cm_z + hz / 2.0, 0.0, 0.0, Type::RIGHT));
		}
		void setPerforated()
		{
			height_perf = 0.0;
			vector<pair<int, int> >::iterator it;
			for (it = perfIntervals.begin(); it != perfIntervals.end(); ++it)
			{
				for (int i = it->first; i <= it->second; i++)
				{
					Qcell[i] = 0.0;
					height_perf += cells[i].hz;
				}
			}
		};
		void setRateDeviation(int num, double ratio)
		{
			Qcell[num] += Q_sum * ratio;
		}
		double solveH()
		{
			double H = 0.0;
			double p1, p0;

			std::map<int, double>::iterator it = Qcell.begin();
			for (int i = 0; i < Qcell.size() - 1; i++)
			{
				p0 = cells[it->first].u_next.p;
				p1 = cells[(++it)->first].u_next.p;

				H += (p1 - p0) * (p1 - p0) / 2.0;
			}

			return H;
		}
		inline double getNablaP(Cell& cell, int varNum, int axis)
		{
			Cell* nebr1;
			Cell* nebr2;
			double h, r_eff;

			switch (axis)
			{
			case R_AXIS:
				nebr1 = &cells[cell.num - cellsNum_z - 2];
				nebr2 = &cells[cell.num + cellsNum_z + 2];

				r_eff = cell.props->radius_eff;
				if ((nebr1->r < r_eff) && (nebr2->r > r_eff))
				{
					if (cell.r > r_eff)
						nebr1 = &cell;
					else
						nebr2 = &cell;
				}
				h = nebr2->r - nebr1->r;

				break;
			case Z_AXIS:
				nebr1 = &cells[cell.num - 1];
				nebr2 = &cells[cell.num + 1];
				h = nebr2->z - nebr1->z;
				break;
			}

			switch (varNum)
			{
			case PREV:
				return (nebr2->u_prev.p - nebr1->u_prev.p) / h;
			case ITER:
				return (nebr2->u_iter.p - nebr1->u_iter.p) / h;
			case NEXT:
				return (nebr2->u_next.p - nebr1->u_next.p) / h;
			}
		};

		// Service functions
		inline double upwindIsCur(int cur, int beta)
		{
			if (cells[cur].u_next.p < cells[beta].u_next.p)
				return 0.0;
			else
				return 1.0;
		};
		inline int getUpwindIdx(int cur, int beta) const
		{
			if (cells[cur].u_next.p < cells[beta].u_next.p)
				return beta;
			else
				return cur;
		};
		inline void getNeighborIdx(int cur, int* const neighbor)
		{
			neighbor[0] = cur - cellsNum_z - 2;
			neighbor[1] = cur + cellsNum_z + 2;
			neighbor[2] = cur - 1;
			neighbor[3] = cur + 1;
		};
		inline void getStencil(const Cell& cell, Cell** const neighbor)
		{
			neighbor[0] = &cell;
			switch (cell.type)
			{
			case Type::MIDDLE:
				neighbor[1] = &cells[cell.num - cellsNum_z - 2];
				neighbor[2] = &cells[cell.num + cellsNum_z + 2];
				neighbor[3] = &cells[cell.num - 1];
				neighbor[4] = &cells[cell.num + 1];
				break;
			case Type::WELL_LAT:
				neighbor[1] = &cells[cell.num + cellsNum_z + 2];
				neighbor[2] = &cells[cell.num + 2 * cellsNum_z + 4];
				break;
			case Type::RIGHT:
				neighbor[1] = &cells[cell.num - cellsNum_z - 2];
				neighbor[2] = &cells[cell.num - 2 * cellsNum_z - 4];
				break;
			case Type::TOP:
				neighbor[1] = &cells[cell.num + 1];
			case Type::BOTTOM:
				neighbor[1] = &cells[cell.num - 1];
			}
		};
		inline int getSkeletonIdx(const Cell& cell) const
		{
			int idx = 0;
			while (idx < props_sk.size())
			{
				if (cell.z <= props_sk[idx].h2 + EQUALITY_TOLERANCE)
					return idx;
				idx++;
			}
			exit(-1);
		};

	public:
		Basic2d() {};
		~Basic2d() {};

		void setPeriod(int period)
		{
			if (leftBoundIsRate)
			{
				Q_sum = rate[period];

				if (period == 0 || rate[period - 1] < EQUALITY_TOLERANCE) {
					std::map<int, double>::iterator it;
					for (it = Qcell.begin(); it != Qcell.end(); ++it)
						it->second = Q_sum * cells[it->first].hz / height_perf;
				}
				else {
					std::map<int, double>::iterator it;
					for (it = Qcell.begin(); it != Qcell.end(); ++it)
						it->second = it->second * Q_sum / rate[period - 1];
				}
			}
			else
			{
				Pwf = pwf[period];
				Q_sum = 0.0;
			}

			for (int i = 0; i < skeletonsNum; i++)
			{
				props_sk[i].radius_eff = props_sk[i].radiuses_eff[period];
				props_sk[i].perm_eff = props_sk[i].perms_eff[period];
				props_sk[i].skin = props_sk[i].skins[period];
			}
		}
	};
}

#endif /* BASIC2D_HPP_ */
