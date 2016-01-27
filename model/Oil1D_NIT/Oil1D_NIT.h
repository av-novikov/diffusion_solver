#ifndef OIL1D_NIT_H_
#define OIL1D_NIT_H_

#include <vector>
#include <map>

#include "util/utils.h"
#include "model/cells/Variables.hpp"
#include "model/cells/RadialCell.hpp"
#include "model/AbstractModel.hpp"

namespace oil1D_NIT
{
	typedef RadialCell<Var1phaseNIT> Cell;

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

		// Thermal diffusivity coefficient [m2/sec]
		double kappa_eff;
		// Mass heat capacity [J/kg/K]
		double c;
		// Thermal conductivity coefficient [W/m/K]
		double lambda;
	};

	struct Oil_Props
	{
		// Viscosity [cP]
		double visc;
		// Density of fluid in STC [kg/m3]
		double dens_stc;
		// Volume factor for well bore
		double b_bore;
		// Compessibility [1/Pa]
		double beta;

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

		// Perforated intervals
		std::vector<std::pair<int,int> > perfIntervals;
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
		// Number of cells in vertical direction
		int cellsNum_z;

		std::vector<Skeleton_Props> props_sk;
		Oil_Props props_oil;

		double depth_point;
	};

	class Oil1D_NIT : public AbstractModel<Var1phaseNIT, Properties, RadialCell, Oil1D_NIT>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class GRDECLSnapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class AbstractMethod;
		friend class Oil1DNITSolver;	
	protected:
		int skeletonsNum;
		std::vector<Skeleton_Props> props_sk;
		Oil_Props props_oil;

		// Number of cells in radial direction
		int cellsNum_r;
		// Number of cells in vertical direction
		int cellsNum_z;

		// BHP will be converted to the depth
		double depth_point;
		// During the time flow rate decreases 'e' times in well test [sec] 
		double alpha;

		void setProps(Properties& props);
		void makeDimLess();
		void buildGridLog();
		void setPerforated();
		void setInitialState();
		
		// Service functions
		inline double upwindIsCur(int cur, int beta) const
		{
			if(cells[cur].u_next.p < cells[beta].u_next.p)
				return 0.0;
			else
				return 1.0;
		};
		inline int getUpwindIdx(int cur, int beta) const
		{
			if(cells[cur].u_next.p < cells[beta].u_next.p)
				return beta;
			else
				return cur;
		};

		// Solving coefficients
		inline double getPoro(double p) const
		{
			return props_sk[0].m * (1.0 + props_sk[0].beta * (props_sk[0].p_init - p) );
		};
		inline double getRho(double p) const
		{
			return props_oil.dens_stc * (1.0 + props_oil.beta * p);
		};
		inline double getTrans(Cell& cell, Cell& beta) const
		{
			double k1, k2;
			k1 = (cell.r > props_sk[0].radius_eff ? props_sk[0].perm_r : props_sk[0].perm_eff);
			k2 = (beta.r > props_sk[0].radius_eff ? props_sk[0].perm_r : props_sk[0].perm_eff);
			return 4.0 * M_PI * k1 * k2 * (cell.r + sign(beta.num - cell.num) * cell.hr / 2.0) * props_sk[0].height / (k2 * cell.hr + k1 * beta.hr);
		};
		inline double getRho(Cell& cell, Cell& beta) const
		{
			return ( getRho(beta.u_next.p) * cell.hr + beta.hr * getRho(cell.u_next.p) ) / (beta.hr + cell.hr);
		};
		inline double getPerm_r(const Cell& cell) const
		{
			return (cell.r > props_sk[0].radius_eff ? props_sk[0].perm_r : props_sk[0].perm_eff);
		};
		inline double getNablaP(Cell& cell, int varNum)
		{
			Cell* nebr1;
			Cell* nebr2;
			double h;

			nebr1 = &cells[cell.num - 1];
			nebr2 = &cells[cell.num + 1];
			h = nebr2->r - nebr1->r;

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
		inline double getOilVelocity(Cell& cell, int varNum)
		{
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

			return -getPerm_r(cell) / props_oil.visc * getNablaP(cell, varNum);
		};
		inline double getCn(Cell& cell) const
		{
			return getPoro(cell.u_next.p) * getRho(cell.u_next.p) * props_oil.c + (1.0 - getPoro(cell.u_next.p)) * props_sk[0].dens_stc * props_sk[0].c;
		};
		inline double getAd(Cell& cell) const
		{
			return getPoro(cell.u_next.p) * getRho(cell.u_next.p) * props_oil.c * props_oil.ad;
		};
		inline double getLambda(Cell& cell) const
		{
			return getPoro(cell.u_next.p) * props_oil.lambda + (1.0-getPoro(cell.u_next.p)) * props_sk[0].lambda;
		};
		inline double getLambdaWorse(Cell& cell) const
		{
				return getPoro(cell.u_next.p) * props_oil.lambda + (1.0-getPoro(cell.u_next.p)) * props_sk[0].lambda * 1.0;
		};
		inline double getLambda(Cell& cell1, Cell& cell2) const
		{
			double cm = (cell1.r + sign(cell2.num - cell1.num) * cell1.hr / 2.0);
			if(cm > 0.1724 / R_dim)
				return ( cell1.hr * getLambda(cell2) + cell2.hr * getLambda(cell1) ) / (cell1.hr + cell2.hr);
			else 
				return ( cell1.hr * getLambdaWorse(cell2) + cell2.hr * getLambdaWorse(cell1) ) / (cell1.hr + cell2.hr);
		};
		inline double getJT(Cell& cell, int varNum)
		{
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
		
			return getRho(var->p) * props_oil.c * props_oil.jt * getOilVelocity(cell, varNum);
		};
		inline double getA(Cell& cell, int varNum)
		{
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
		
			return getRho(var->p) * props_oil.c * getOilVelocity(cell, varNum);
		};

		double solve_eq(int i);
		double solve_eq_dp(int i);
		double solve_eq_dp_beta(int i, int beta);
		
		double solve_eqLeft();
		double solve_eqLeft_dp();
		double solve_eqLeft_dp_beta();
		
		double solve_eqRight();
		double solve_eqRight_dp_beta();
		double solve_eqRight_dp();

	public:
		Oil1D_NIT();
		~Oil1D_NIT();

		void setPeriod(int period);
		double getRate();
	};

};

#endif OIL1D_NIT_H_

