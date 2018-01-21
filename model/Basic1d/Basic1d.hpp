#ifndef BASIC1D_HPP_
#define BASIC1D_HPP_

#include "model/cells/Variables.hpp"
#include "model/Basic1d/Properties.hpp"
#include "model/AbstractModel.hpp"

#include <vector>
#include <map>
#include <cassert>

namespace basic1d
{
	static const int R_AXIS = 0;

	static const int stencil = 3;
	static const int Lstencil = 3;
	static const int Rstencil = 3;

	template <typename varType, typename propsType, typename sk_propsType,
	template <typename varType> class cellType, class modelType>
	class Basic1d : public AbstractModel<varType, propsType, cellType, modelType>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class Basic1dSolver;
	public:
		typedef cellType<varType> Cell;
		typename typedef Cell::Type Type;
	protected:

		sk_propsType props_sk;

		// Continuum properties
		int skeletonsNum;
		// Number of cells in radial direction
		int cellsNum_x;

		void setPerforated()
		{
			Qcell[0] = 0.0;
		};
		void setBasicProps(propsType& props)
		{
			leftBoundIsRate = props.leftBoundIsRate;
			rightBoundIsPres = props.rightBoundIsPres;

			// Setting grid properties
			r_w = props.r_w;
			r_e = props.r_e;
			cellsNum_x = props.cellsNum_x;
			cellsNum = cellsNum_x + 2;

			props_sk.perm = MilliDarcyToM2(props_sk.perm);

			periodsNum = props.timePeriods.size();
			for (int i = 0; i < periodsNum; i++)
			{
				period.push_back(props.timePeriods[i]);
				if (leftBoundIsRate)
					rate.push_back(props.rates[i] / 86400.0);
				else
					pwf.push_back(props.pwf[i]);
				if (props_sk.radiuses_eff[i] > props.r_w)
					props_sk.perms_eff.push_back(MilliDarcyToM2(props.props_sk.perm * log(props.props_sk.radiuses_eff[i] / props.r_w) / (log(props.props_sk.radiuses_eff[i] / props.r_w) + props.props_sk.skins[i])));
				else
					props_sk.perms_eff.push_back(MilliDarcyToM2(props.props_sk.perm));
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
			R_dim = r_e / 100.0;
			t_dim = 3600.0;
			P_dim = props_sk.p_init;

			// Temporal properties
			ht /= t_dim;
			ht_min /= t_dim;
			ht_max /= t_dim;

			// Grid properties
			r_w /= R_dim;
			r_e /= R_dim;

			// Skeleton properties
			props_sk.perm /= (R_dim * R_dim);

			props_sk.beta /= (1.0 / P_dim);
			props_sk.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
			props_sk.h1 = (props_sk.h1 - depth_point) / R_dim;
			props_sk.h2 = (props_sk.h2 - depth_point) / R_dim;
			props_sk.height /= R_dim;
			props_sk.p_init /= P_dim;
			props_sk.p_out /= P_dim;
			props_sk.p_ref /= P_dim;

			for (int j = 0; j < periodsNum; j++)
			{
				props_sk.perms_eff[j] /= (R_dim * R_dim);
				props_sk.radiuses_eff[j] /= R_dim;
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

			grav /= (R_dim / t_dim / t_dim);
			// Rest properties
			alpha /= t_dim;
			//depth_point = 0.0;
		};
		void buildGridLog()
		{
			cells.reserve(cellsNum);

			Volume = 0.0;
			int counter = 0;

			double r_prev = r_w;
			double logMax = log(r_e / r_w);
			double logStep = logMax / (double)cellsNum_x;

			double hz = props_sk.h2 - props_sk.h1;
			//double cm_z = props_sk[skel_idx].h1;
			double hr = r_prev * (exp(logStep) - 1.0);
			double cm_r = r_w;

			// Left border
			cells.push_back(Cell(counter++, cm_r, 0.0, hz, Type::WELL_LAT));
			// Middle cells
			for (int j = 0; j < cellsNum_x; j++)
			{
				cm_r = r_prev * (exp(logStep) + 1.0) / 2.0;
				hr = r_prev * (exp(logStep) - 1.0);

				cells.push_back(Cell(counter++, cm_r, hr, hz, Type::MIDDLE));
				Volume += cells[cells.size() - 1].V;
				r_prev = r_prev * exp(logStep);
			}

			// Right border
			cm_r = r_e;
			cells.push_back(Cell(counter++, cm_r, 0.0, hz, Type::RIGHT));
		}
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
		inline double getNablaP(Cell& cell, int varNum)
		{
			Cell* nebr1;
			Cell* nebr2;
			double h, r_eff;

			nebr1 = &cells[cell.num - 1];
			nebr2 = &cells[cell.num + 1];

			r_eff = cell.props->radius_eff;
			if ((nebr1->x < r_eff) && (nebr2->x > r_eff))
			{
				if (cell.x > r_eff)
					nebr1 = &cell;
				else
					nebr2 = &cell;
			}
			h = nebr2->x - nebr1->x;
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
			neighbor[0] = cur - 1;
			neighbor[1] = cur + 1;
		};
		inline void getStencil(const Cell& cell, Cell** const neighbor)
		{
			neighbor[0] = &cell;
			switch (cell.type)
			{
			case Type::MIDDLE:
				neighbor[1] = &cells[cell.num - 1];
				neighbor[2] = &cells[cell.num + 1];
				break;
			case Type::WELL_LAT:
				neighbor[1] = &cells[cell.num + 1];
				neighbor[2] = &cells[cell.num + 2];
				break;
			case Type::RIGHT:
				neighbor[1] = &cells[cell.num - 1];
				neighbor[2] = &cells[cell.num - 2];
				break;
			}
		};

	public:
		Basic1d() {};
		~Basic1d() {};

		virtual void setPeriod(int period)
		{
			if (leftBoundIsRate)
			{
				Q_sum = rate[period];

				if (period == 0 || rate[period - 1] < EQUALITY_TOLERANCE) {
					std::map<int, double>::iterator it;
					for (it = Qcell.begin(); it != Qcell.end(); ++it)
						it->second = Q_sum;
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
			props_sk.radius_eff = props_sk.radiuses_eff[period];
			props_sk.perm_eff = props_sk.perms_eff[period];
			props_sk.skin = props_sk.skins[period];
		}
	};
}

#endif /* BASIC1D_HPP_ */
