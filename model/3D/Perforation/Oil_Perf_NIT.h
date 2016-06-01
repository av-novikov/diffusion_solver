#ifndef OIL_PERF_NIT_H_
#define OIL_PERF_NIT_H_

#include <vector>
#include <map>
#include <string>

#include "model/cells/Iterators.h"
#include "model/cells/Variables.hpp"
#include "model/cells/CylCellPerf.h"
#include "model/AbstractModel.hpp"
#include "util/Interpolate.h"
#include "util/utils.h"

namespace oil_perf_nit
{
	typedef CylCellPerf<Var1phaseNIT> Cell;
	typedef Iterator<CylCellPerf<Var1phaseNIT> > Iterator;

	struct Skeleton_Props
	{
		// Porosity in STC
		double m; 
		// Density of skeleton matter in STC [kg/m3]
		double dens_stc;
		// Compessibility [1/Pa]
		double beta;
		// Permeability along radial direction [mD]
		double perm_r;
		// Permeability along vertical direction [mD]
		double perm_z;

		// Permeability of colmatage zone [mD]
		std::vector<double> perms_eff;
		// Radius of colmatage zone [m]
		std::vector<double> radiuses_eff;
		// Vector of skins
		std::vector<double> skins;
		double perm_eff;
		double radius_eff;
		double skin;

		// Top and bottom depth of perforation
		double h1, h2;
		// Height of formation [m]
		double height;

		int cellsNum_z;

		double p_out;
		double p_init;
		double t_init;

		// Thermal properties

		// Thermal diffusivity coefficient [m2/sec]
		double kappa_eff;
		// Mass heat capacity [J/kg/K]
		double c;
		// Thermal conductivity coefficient [W/m/K]
		double lambda_r;
		// Thermal conductivity coefficient [W/m/K]
		double lambda_z;
	};

	struct Fluid_Props
	{
		// Viscosity [cP]
		double visc;
		// Density of fluid in STC [kg/m3]
		double dens_stc;
		// Volume factor for well bore
		double b_bore;
		// Compessibility [1/Pa]
		double beta;
		// Relative fluid permeability
		Interpolate* kr;
		// Fluid volume factor
		Interpolate* b;

		// Thermal properties

		// Mass heat capacity [J/kg/K]
		double c;
		// Thermal conductivity coefficient [W/m/K]
		double lambda;
		// Joule-thompson coefficient [K/Pa]
		double jt;
		// Adiabatic coefficient [K/Pa]
		double ad;
	};

	struct Properties
	{
		// Vector of start times of periods [sec]
		std::vector<double> timePeriods;
		// Vector of rates [m3/day]
		std::vector<double> rates;
		// Vector of BHPs [Pa]
		std::vector<double> pwf;

		// If left boundary condition would be 2nd type
		bool leftBoundIsRate;
		// If right boundary condition would be 1st type
		bool rightBoundIsPres;
	
		// Perforated tunnels
		/* Cell number && number of cells in depth */
		std::vector<std::pair<int,int> > perfTunnels;
		// Time step limits
		// Initial time step [sec]
		double ht;
		// Minimal time step [sec]
		double ht_min;
		// Maximum time step [sec]
		double ht_max;
		// During the time flow rate decreases 'e' times in well test [sec] 
		double alpha;

		// Inner radius of well [m]
		double r_w;
		// Radius of formation [m]
		double r_e;
	
		// Number of cells in radial direction
		int cellsNum_r;
		// Number of cells in axial direction
		int cellsNum_phi;
		// Number of cells in vertical direction
		int cellsNum_z;

		std::vector<Skeleton_Props> props_sk;
		Fluid_Props props_oil;

		double depth_point;
	};

	class Oil_Perf_NIT : public AbstractModel<Var1phaseNIT, Properties, CylCellPerf, Oil_Perf_NIT>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class GRDECLSnapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class AbstractSolver;
		friend class OilPerfNITSolver;
		template<typename> friend class MidStencil;
		template<typename> friend class LeftStencil;
		template<typename> friend class RightStencil;
		template<typename> friend class TopStencil;
		template<typename> friend class BotStencil;
		template<typename> friend class UsedStencils;

	protected:

		// Middle iterator
		Iterator* midIter;
		Iterator* midBegin;
		Iterator* midEnd;

		// Left iterator
		Iterator* leftIter;
		Iterator* leftBegin;
		Iterator* leftEnd;

		// Right iterator
		Iterator* rightIter;
		Iterator* rightBegin;
		Iterator* rightEnd;

		// Continuum properties
		int skeletonsNum;
		std::vector<Skeleton_Props> props_sk;
		Fluid_Props props_oil;

		// Number of cells in radial direction
		int cellsNum_r;
		// Number of cells in axial direction
		int cellsNum_phi;
		// Number of cells in vertical direction
		int cellsNum_z;

		// Structure to construct tunnelss
		std::vector<std::pair<int, int> > perfTunnels;
		// Border cells in tunnels
		std::vector<Cell> tunnelCells;
		void buildTunnels();
		void setUnused();
		std::map<int,int> tunnelNebrMap;
		std::map<int,std::pair<int, int> > nebrMap;

		inline Cell& getCell(int num)
		{
			Cell& cell = cells[num];
			//assert(cell.isUsed);

			return cell;
		};

		inline Cell& getCell(int num, int beta)
		{
			Cell& cell = cells[num];
			//assert(cell.isUsed);

			Cell& nebr = cells[beta];
			if (nebr.isUsed)
				return nebr;
			else
				return tunnelCells[tunnelNebrMap[num]];
		};

		inline const Cell& getCell(const int num) const
		{
			const Cell& cell = cells[num];
			//assert(cell.isUsed);

			return cell;
		};

		inline const Cell& getCell(const int num, const int beta) const
		{
			const Cell& cell = cells[num];
			//assert(cell.isUsed);

			const Cell& nebr = cells[beta];
			if (nebr.isUsed)
				return nebr;
			else
				return tunnelCells[ tunnelNebrMap.at(num) ];
		};

		// BHP will be converted to the depth
		double depth_point;
		// During the time flow rate decreases 'e' times in well test [sec] 
		double alpha;

		void setInitialState();
		// Set all properties
		void setProps(Properties& props);
		// Make all properties dimensionless
		void makeDimLess();
		// Build grid
		void buildGridLog();
		// Set perforated cells
		void setPerforated();
		// Set some deviation to rate distribution
		void setRateDeviation(int num, double ratio);
		// Check formations properties
		void checkSkeletons(const std::vector<Skeleton_Props>& props);

		// Service functions
		inline double upwindIsCur(const Cell* cell, const Cell* nebr) const
		{
			if(cell->u_next.p < nebr->u_next.p)
				return 0.0;
			else
				return 1.0;
		};
		inline const Cell* getUpwindIdx(const Cell* cell, const Cell* nebr) const
		{
			assert(cell->isUsed);
			assert(nebr->isUsed);
			if(cell->u_next.p < nebr->u_next.p)
				return nebr;
			else
				return cell;
		};
		inline void getNeighborIdx(Cell& cur, Cell** const neighbor)
		{
			assert(cur.isUsed);
			neighbor[0] = &getCell(cur.num, cur.num - cellsNum_z - 2);
			neighbor[1] = &getCell(cur.num, cur.num + cellsNum_z + 2);
			neighbor[2] = &getCell(cur.num, cur.num - 1);
			neighbor[3] = &getCell(cur.num, cur.num + 1);
			if (cur.num < (cellsNum_r + 2) * (cellsNum_z + 2))
				neighbor[4] = &getCell(cur.num, cur.num + (cellsNum_r + 2) * (cellsNum_z + 2) * (cellsNum_phi - 1));
			else
				neighbor[4] = &getCell(cur.num, cur.num - (cellsNum_r + 2) * (cellsNum_z + 2));
			if (cur.num < (cellsNum_r + 2) * (cellsNum_z + 2) * (cellsNum_phi - 1))
				neighbor[5] = &getCell(cur.num, cur.num + (cellsNum_r + 2) * (cellsNum_z + 2));
			else
				neighbor[5] = &getCell(cur.num, cur.num - (cellsNum_r + 2) * (cellsNum_z + 2) * (cellsNum_phi - 1));
		};
		inline void getStencilIdx(int cur, Cell** const neighbor)
		{
			neighbor[0] = &getCell(cur);
			neighbor[1] = &getCell(cur - cellsNum_z - 2);
			neighbor[2] = &getCell(cur + cellsNum_z + 2);
			neighbor[3] = &getCell(cur - 1);
			neighbor[4] = &getCell(cur + 1);
			if (cur < (cellsNum_r + 2) * (cellsNum_z + 2))
				neighbor[5] = &getCell(cur + (cellsNum_r + 2) * (cellsNum_z + 2) * (cellsNum_phi - 1));
			else
				neighbor[5] = &getCell(cur - (cellsNum_r + 2) * (cellsNum_z + 2));
			if (cur < (cellsNum_r + 2) * (cellsNum_z + 2) * (cellsNum_phi - 1))
				neighbor[6] = &getCell(cur + (cellsNum_r + 2) * (cellsNum_z + 2));
			else
				neighbor[6] = &getCell(cur - (cellsNum_r + 2) * (cellsNum_z + 2) * (cellsNum_phi - 1));
		};
		inline int getIdx(int i)
		{
			if( i < 0 )
				return cellsNum + i;
			else if( i > cellsNum)
				return i - cellsNum;
			else
				return i;
		}
		inline int getSkeletonIdx(const Cell& cell) const
		{
			if (!cell.isTunnel)
			{
				int idx = 0;
				while (idx < props_sk.size())
				{
					if (cell.z <= props_sk[idx].h2 + EQUALITY_TOLERANCE)
						return idx;
					idx++;
				}
				exit(-1);
			}
			else
				return getSkeletonIdx( getCell(perfTunnels[cell.tunNum].first) );
		};

		// Solving coefficients
		inline double getPoro(double p, Cell& cell) const
		{
			const int idx = getSkeletonIdx(cell);
			return props_sk[idx].m * (1.0 + props_sk[idx].beta * p );
		};
		inline double getPoro_dp(Cell& cell) const
		{
			const int idx = getSkeletonIdx(cell);
			return props_sk[idx].m * props_sk[idx].beta;
		};
		inline double getTrans(Cell& cell, Cell& beta)
		{
			double k1, k2, S;
			const int idx1 = getSkeletonIdx(cell);
			const int idx2 = getSkeletonIdx(beta);

			if( fabs(cell.z - beta.z) > EQUALITY_TOLERANCE ) {
				k1 = props_sk[idx1].perm_z;
				k2 = props_sk[idx2].perm_z;
				if(k1 == 0.0 && k2 == 0.0)
					return 0.0;
				S = cell.hphi * cell.r * cell.hr;
				return 2.0 * k1 * k2 * S / (k1 * beta.hz + k2 * cell.hz);
			} else if( fabs(cell.r - beta.r) > EQUALITY_TOLERANCE) {
				k1 = (cell.r > props_sk[idx1].radius_eff ? props_sk[idx1].perm_r : props_sk[idx1].perm_eff);
				k2 = (beta.r > props_sk[idx2].radius_eff ? props_sk[idx2].perm_r : props_sk[idx2].perm_eff);
				if(k1 == 0.0 && k2 == 0.0)
					return 0.0;
				S = cell.hphi * cell.hz * (cell.r + sign(beta.num - cell.num) * cell.hr / 2.0);
				return 2.0 * k1 * k2 * S / (k1 * beta.hr + k2 * cell.hr);
			} else if(fabs(cell.phi - beta.phi) > EQUALITY_TOLERANCE) {
				k1 = (cell.r > props_sk[idx1].radius_eff ? props_sk[idx1].perm_r : props_sk[idx1].perm_eff);
				if(k1 == 0.0)
					return 0.0;
				S = cell.hr * cell.hz;
				return 2.0 * k1 * S / (beta.hr + cell.hr);
			}
		};
		inline double getPerm_r(const Cell& cell) const
		{
			const int idx = getSkeletonIdx(cell);
			return (cell.r > props_sk[idx].radius_eff ? props_sk[idx].perm_r : props_sk[idx].perm_eff);
		};
		inline double getBoreB_oil(double p) const
		{
			return props_oil.b_bore;
			//return getB_oil(p, p_bub, SATUR);
		};

		inline double getRho(double p) const
		{
			return props_oil.dens_stc * (1.0 + props_oil.beta * p);
		};
		inline double getRho_dp() const
		{
			return props_oil.dens_stc * props_oil.beta;
		};
		inline double getRho(Cell& cell1, Cell& cell2) const
		{
			if (fabs(cell1.z - cell2.z) > EQUALITY_TOLERANCE)
				return (cell1.hz * getRho(cell2.u_next.p) + cell2.hz * getRho(cell1.u_next.p)) / (cell1.hz + cell2.hz);
			else if (fabs(cell1.r - cell2.r) > EQUALITY_TOLERANCE)
				return (cell1.hr * getRho(cell2.u_next.p) + cell2.hr * getRho(cell1.u_next.p)) / (cell1.hr + cell2.hr);
			else if (fabs(cell1.phi - cell2.phi) > EQUALITY_TOLERANCE)
				return (cell1.hphi * getRho(cell2.u_next.p) + cell2.hphi * getRho(cell1.u_next.p)) / (cell1.hphi + cell2.hphi);
		};
		inline double getRho_dp(Cell& cell1, Cell& cell2) const
		{
			if (fabs(cell1.z - cell2.z) > EQUALITY_TOLERANCE)
				return cell2.hz * getRho_dp() / (cell1.hz + cell2.hz);
			else if (fabs(cell1.r - cell2.r) > EQUALITY_TOLERANCE)
				return cell2.hr * getRho_dp() / (cell1.hr + cell2.hr);
			else if (fabs(cell1.phi - cell2.phi) > EQUALITY_TOLERANCE)
				return cell2.hphi * getRho_dp() / (cell1.hphi + cell2.hphi);
		};
		inline double getRho_dp_beta(Cell& cell1, Cell& cell2) const
		{
			if (fabs(cell1.z - cell2.z) > EQUALITY_TOLERANCE)
				return cell1.hz * getRho_dp() / (cell1.hz + cell2.hz);
			else if (fabs(cell1.r - cell2.r) > EQUALITY_TOLERANCE)
				return cell1.hr * getRho_dp() / (cell1.hr + cell2.hr);
			else if (fabs(cell1.phi - cell2.phi) > EQUALITY_TOLERANCE)
				return cell1.hphi * getRho_dp() / (cell1.hphi + cell2.hphi);
		};

		// Thermal functions
		inline double getNablaP(Cell& cell, int varNum, int axis)
		{
			Cell* nebr1;
			Cell* nebr2;
			double h;
			double r_eff;
			int idx = getSkeletonIdx(cell);

			switch(axis)
			{
			case R_AXIS:
				nebr1 = &getCell(cell.num, cell.num - cellsNum_z - 2);
				nebr2 = &getCell(cell.num, cell.num + cellsNum_z + 2);
				
				r_eff = props_sk[idx].radius_eff;
				if ( (nebr1->r < r_eff) && (nebr2->r > r_eff) )
				{
					if (cell.r > r_eff)
						nebr1 = &cell;
					else
						nebr2 = &cell;
				}
				h = nebr2->r - nebr1->r;
				break;
			case PHI_AXIS:
				nebr1 = &getCell(cell.num, getIdx(cell.num - (cellsNum_r + 2) * (cellsNum_z + 2)) );
				nebr2 = &getCell(cell.num, getIdx(cell.num + (cellsNum_r + 2) * (cellsNum_z + 2)) );
				h = nebr1->r * (nebr2->phi - nebr1->phi);
				if(abs(nebr2->phi - nebr1->phi) > 2.0 * cell.hphi + EQUALITY_TOLERANCE)
					h = nebr1->r * (nebr2->phi - nebr1->phi + 2.0 * M_PI);
				break;
			case Z_AXIS:
				nebr1 = &getCell(cell.num, cell.num - 1);
				nebr2 = &getCell(cell.num, cell.num + 1);
				h = nebr2->z - nebr1->z;
				break;
			}

			switch(varNum)
			{
			case PREV:
				return (nebr2->u_prev.p - nebr1->u_prev.p ) / h;
			case ITER:
				return (nebr2->u_iter.p - nebr1->u_iter.p ) / h;
			case NEXT:
				return (nebr2->u_next.p - nebr1->u_next.p ) / h;
			}
		};
		inline double getOilVelocity(Cell& cell, int varNum, int axis)
		{
			const int idx = getSkeletonIdx( cell );

			Var1phaseNIT* var;
			switch(varNum)
			{
			case PREV:
				var = &cell.u_prev;
				break;
			case ITER:
				var = &cell.u_iter;
				break;
			case NEXT:
				var = &cell.u_next;
				break;
			}

			switch(axis)
			{
			case R_AXIS:
				return -getPerm_r(cell) / props_oil.visc * getNablaP(cell, varNum, axis);
			case PHI_AXIS:
				return -getPerm_r(cell) / props_oil.visc * getNablaP(cell, varNum, axis);
			case Z_AXIS:
				return -props_sk[idx].perm_z / props_oil.visc * getNablaP(cell, varNum, axis);
			}
		};

		inline double getCn(Cell& cell) const
		{
			const int idx = getSkeletonIdx(cell);
			return getPoro(cell.u_next.p, cell) * getRho(cell.u_next.p) * props_oil.c +
				(1.0 - getPoro(cell.u_next.p, cell)) * props_sk[idx].dens_stc * props_sk[idx].c;
		};
		inline double getAd(Cell& cell) const
		{
			return getPoro(cell.u_next.p, cell) * getRho(cell.u_next.p) * props_oil.c * props_oil.ad;
		};
		inline double getLambda(Cell& cell, int axis)
		{
			const int idx = getSkeletonIdx(cell);
			switch (axis)
			{
			case R_AXIS:
				return getPoro(cell.u_next.p, cell) * props_oil.lambda +
					(1.0 - getPoro(cell.u_next.p, cell)) * props_sk[idx].lambda_r;
			case PHI_AXIS:
				return getPoro(cell.u_next.p, cell) * props_oil.lambda +
					(1.0 - getPoro(cell.u_next.p, cell)) * props_sk[idx].lambda_r;
			case Z_AXIS:
				return getPoro(cell.u_next.p, cell) * props_oil.lambda +
					(1.0 - getPoro(cell.u_next.p, cell)) * props_sk[idx].lambda_z;
			}
		};
		inline double getLambda(Cell& cell1, Cell& cell2)
		{
			if (fabs(cell1.z - cell2.z) > EQUALITY_TOLERANCE)
				return (cell1.hz * getLambda(cell2, Z_AXIS) + cell2.hz * getLambda(cell1, Z_AXIS)) / (cell1.hz + cell2.hz);
			else if (fabs(cell1.r - cell2.r) > EQUALITY_TOLERANCE)
				return (cell1.hr * getLambda(cell2, R_AXIS) + cell2.hr * getLambda(cell1, R_AXIS)) / (cell1.hr + cell2.hr);
			else if (fabs(cell1.phi - cell2.phi) > EQUALITY_TOLERANCE)
				return (cell1.hphi * getLambda(cell2, PHI_AXIS) + cell2.hphi * getLambda(cell1, PHI_AXIS)) / (cell1.hphi + cell2.hphi);
		};
		inline double getJT(Cell& cell, int varNum, int axis)
		{
			Var1phaseNIT* var;
			switch (varNum)
			{
			case PREV:
				var = &cell.u_prev;
				break;
			case ITER:
				var = &cell.u_iter;
				break;
			case NEXT:
				var = &cell.u_next;
				break;
			}

			return getRho(var->p) * props_oil.c * props_oil.jt * getOilVelocity(cell, varNum, axis);
		};
		inline double getA(Cell& cell, int varNum, int axis)
		{
			Var1phaseNIT* var;
			switch (varNum)
			{
			case PREV:
				var = &cell.u_prev;
				break;
			case ITER:
				var = &cell.u_iter;
				break;
			case NEXT:
				var = &cell.u_next;
				break;
			}

			return getRho(var->p) * props_oil.c * getOilVelocity(cell, varNum, axis);
		};

		// First eqn
		double solve_eq(int cur);
		double solve_eq_dp(int cur, int beta);
		double solve_eq_dp_beta(int cur, int beta);

		/*-------------- Left cells ------------------*/

		inline double solve_eqLeft(int cur)
		{
			Cell& cell = tunnelCells[cur];
			Cell& nebr = getCell(nebrMap[cur].first);
			const Var1phaseNIT& next = cell.u_next;
			const Var1phaseNIT& upwd = getUpwindIdx(&cell, &nebr)->u_next;

			if (leftBoundIsRate)
				return getTrans(cell, nebr) / props_oil.visc / getBoreB_oil(next.p) * (nebr.u_next.p - next.p) - Qcell[cur];
			else
				return next.p - Pwf;
		}

		inline double solve_eqLeft_dp(int cur, int beta)
		{
			Cell& cell = tunnelCells[cur];
			Cell& nebr = getCell(nebrMap[cur].first);
			const Var1phaseNIT& next = cell.u_next;
			const Var1phaseNIT& upwd = getUpwindIdx(&cell, &nebr)->u_next;

			if (leftBoundIsRate)
				return -getTrans(cell, nebr) / getBoreB_oil(next.p) / props_oil.visc;
			else
				return 1.0;
		}

		inline double solve_eqLeft_dp_beta(int cur, int beta)
		{
			Cell& cell = tunnelCells[cur];
			Cell& nebr = getCell(nebrMap[cur].first);
			const Var1phaseNIT& next = cell.u_next;
			const Var1phaseNIT& upwd = getUpwindIdx(&cell, &nebr)->u_next;

			if (leftBoundIsRate)
				return getTrans(cell, nebr) / getBoreB_oil(next.p) / props_oil.visc;
			else
				return 0.0;
		}

		/*-------------- Right cells ------------------*/

		inline double solve_eqRight(int cur)
		{
			const Cell& cell = getCell(cur);

			if (rightBoundIsPres)
				return cell.u_next.p - props_sk[getSkeletonIdx(cell)].p_out;
			else
				return cell.u_next.p - getCell(cur - cellsNum_z - 2).u_next.p;
		}

		inline double solve_eqRight_dp(int cur, int beta)
		{
			if (rightBoundIsPres)
				return 1.0;
			else
				return 1.0;
		}

		inline double solve_eqRight_dp_beta(int cur, int beta)
		{
			if (rightBoundIsPres)
				return 0.0;
			else
				return -1.0;
		}

		/*-------------- Top cells ------------------*/

		inline double solve_eqTop(int cur)
		{
			return getCell(cur).u_next.p - getCell(cur + 1).u_next.p;
		}

		inline double solve_eqTop_dp(int cur, int beta)
		{
			return 1.0;
		}

		inline double solve_eqTop_dp_beta(int cur, int beta)
		{
			return -1.0;
		}

		/*-------------- Bot cells ------------------*/

		inline double solve_eqBot(int cur)
		{
			return getCell(cur).u_next.p - getCell(cur - 1).u_next.p;
		}

		inline double solve_eqBot_dp(int cur, int beta)
		{
			return 1.0;
		}

		inline double solve_eqBot_dp_beta(int cur, int beta)
		{
			return -1.0;
		}

		/*-------------- Thermal funcs ------------------*/

		// Finds functional
		double solveH();

		FillFoo middleFoo;
		FillFoo rightFoo;
		FillFoo leftFoo;
		FillFoo topFoo;
		FillFoo botFoo;

	public:
		Oil_Perf_NIT();
		~Oil_Perf_NIT();
	
		void setPeriod(int period);
		double getRate(int cur);


		inline Iterator getMidIter()
		{
			return Iterator(*midIter);
		};
		inline const Iterator& getMidBegin()
		{
			return *midBegin;
		};
		inline const Iterator& getMidEnd()
		{
			return *midEnd;
		};

		inline Iterator getLeftIter()
		{
			return Iterator(*leftIter);
		};
		inline const Iterator& getLeftBegin()
		{
			return *leftBegin;
		};
		inline const Iterator& getLeftEnd()
		{
			return *leftEnd;
		};

		inline Iterator getRightIter()
		{
			return Iterator(*rightIter);
		};
		inline const Iterator& getRightBegin()
		{
			return *rightBegin;
		};
		inline const Iterator& getRightEnd()
		{
			return *rightEnd;
		};

	};
};

#endif /* OIL_PERF_NIT_H_ */