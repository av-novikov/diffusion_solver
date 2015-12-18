#ifndef GASOIL_RZ_NIT_H_
#define GASOIL_RZ_NIT_H_

#include <vector>
#include <map>
#include <string>

#include "model/cells/Variables.hpp"
#include "model/cells/CylindricalCell.h"
#include "model/AbstractModel.hpp"
#include "util/Interpolate.h"
#include "util/utils.h"

namespace gasOil_rz_NIT
{
	typedef CylCell<Var2phaseNIT> Cell;

	struct Properties
	{
		// Vector of start times of periods [sec]
		std::vector<double> timePeriods;
		// Vector of rates [m3/day]
		std::vector<double> rates;
		// Vector of skins
		std::vector<double> skins;
		// Radius of damaged zone [m]
		std::vector<double> radius;
	
		// Time step limits
		// Initial time step [sec]
		double ht;
		// Minimal time step [sec]
		double ht_min;
		// Maximum time step [sec]
		double ht_max;
		// During the time flow rate decreases 'e' times in well test [sec] 
		double alpha;

		// Perforated intervals
		std::vector<std::pair<int,int> > perfIntervals;
		// Inner radius of well [m]
		double r_w;
		// Radius of formation [m]
		double r_e;
	
		// Number of cells in radial direction
		int cellsNum_r;
		// Number of cells in vertical direction
		int cellsNum_z;
	
		// Porosity
		double m; 
		// Height of formation [m]
		double h1, h2;
		// BHP will be converted to the depth
		double depth_point;

		// Compessibilities [1/Pa]
		double beta_oil;
		double beta_sk;

		// Permeability along radial direction [mD]
		std::vector<double> perm_r;
		// Vertical permeabilities of cells [mD]
		std::vector<double> perm_z;
		// Permeability along vertical direction [mD]
		//double perm_z;

		// Data set (saturation, relative oil permeability)
		std::vector< std::pair<double,double> > kr_oil;
		// Data set (saturation, relative gas permeability)
		std::vector< std::pair<double,double> > kr_gas;
	
		// Viscosity of oil [cP]
		double visc_oil;
		// Viscosity of gas [cP]
		double visc_gas;

		// Data set (pressure, oil volume factor) ([Pa], [m3/m3])
		std::vector< std::pair<double,double> > B_oil;
		// Data set (pressure, gas volume factor) ([Pa], [m3/m3])
		std::vector< std::pair<double,double> > B_gas;

		// Data set (pressure, gas content in oil) ([Pa], [m3/m3])
		std::vector< std::pair<double,double> > Rs;
	
		// Densities in STC [kg/m3]
		double dens_oil_stc;
		double dens_gas_stc;
		double dens_sk_stc;

		// Volume factor for well bore
		double b_oil_bore;

		// Initial formation pressure, right boundary pressure [Pa]
		double p_init;
		// Pressure of fully saturated oil [Pa]
		double p_sat;

		// Initial oil saturation, right boundary of saturation
		double s_init;
		// Initial temperature, right boundary temperature [K]
		double T_init;

		// Thermal properties
		// Joule-thompson coefficients [K/Pa]
		double jt_oil;
		double jt_gas;

		// Adiabatic coefficients [K/Pa]
		double ad_oil;
		double ad_gas;

		// Mass heat capacities [J/kg/K]
		double c_oil;
		double c_gas;
		double c_sk;

		// Heat of phase transition [J/kg]
		double L;

		// Thermal diffusivity coefficient [m2/sec]
		double kappa_eff;

		// Thermal conductivity coefficients [W/m/K]
		double lambda_sk_r;
		double lambda_sk_z;
		double lambda_oil;
		double lambda_gas;
	};

	struct Skeleton_Props
	{
		// Porosity in STC
		double m; 
		// Permeability along radial direction [mD]
		std::vector<double> perm_r;
		// Permeability along vertical direction [mD]
		std::vector<double> perm_z;
		// Permeability of colmatage zone [mD]
		std::vector<std::vector<double> > perm_eff;
		// Radius of colmatage zone [m]
		std::vector<std::vector<double> > r_eff;
		// Vector of skins
		std::vector<std::vector<double> > skin;
		// Top and bottom depth of perforation
		double h1, h2;

		// Height of formation [m]
		double height;
		// Density of skeleton matter in STC [kg/m3]
		double dens_stc;
		// Compessibility [1/Pa]
		double beta;

		// Thermal properties

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

	class GasOil_RZ_NIT : public AbstractModel<Var2phaseNIT, Properties, CylCell>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class GRDECLSnapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class AbstractMethod;
		friend class GasOil2DNITSolver;	

	protected:
		// Continuum properties
		Skeleton_Props props_sk;
		Fluid_Props props_oil;
		Fluid_Props props_gas;

		std::vector<double> r_eff;
		std::vector<double> Perm_eff;
		std::vector<double> skin;

		// Number of cells in radial direction
		int cellsNum_r;
		// Number of cells in vertical direction
		int cellsNum_z;

		// Gas content in oil
		Interpolate* Rs;
		Interpolate* Prs;

		// Heat of phase transition [J/kg]
		double L;
		// BHP will be converted to the depth
		double depth_point;
		// During the time flow rate decreases 'e' times in well test [sec] 
		double alpha;

		// Set all properties
		void setProps(Properties& props);
		// Make all properties dimensionless
		void makeDimLess();
		// Build grid
		void buildGridLog();
		// Set perforated cells
		void setPerforated();

		// Service functions
		inline double upwindIsCur(int cur, int beta)
		{
			if(cells[cur].u_next.p < cells[beta].u_next.p)
				return 0.0;
			else
				return 1.0;
			//return 0.0 ? cells[cur].u_next.p < cells[beta].u_next.p : 1.0;
		};
		inline int getUpwindIdx(int cur, int beta)
		{
			if(cells[cur].u_next.p < cells[beta].u_next.p)
				return beta;
			else
				return cur;
			//return beta ? cells[cur].u_next.p < cells[beta].u_next.p : cur;
		};
		inline void getNeighborIdx(int cur, int* const neighbor)
		{
			neighbor[0] = cur - cellsNum_z - 2; 
			neighbor[1] = cur + cellsNum_z + 2;
			neighbor[2] = cur - 1;
			neighbor[3] = cur + 1;
		};

		// Solving coefficients
		inline double getPoro(double p)
		{
			return props_sk.m * (1.0 + props_sk.beta * (p - varInit.p) );
		};
		inline double getPoro_dp(double p)
		{
			return props_sk.m * props_sk.beta;
		};
		inline double getTrans(int idx1, int idx2)
		{
			Cell& cell1 = cells[idx1];
			Cell& cell2 = cells[idx2];
			double k1, k2, S;

			if( abs(idx1 - idx2) == 1) {
				k1 = props_sk.perm_z[idx1 % (cellsNum_z+2)];
				k2 = props_sk.perm_z[idx2 % (cellsNum_z+2)];
				if(k1 == 0.0 && k2 == 0.0)
					return 0.0;
				S = 2.0 * M_PI * cell1.r * cell1.hr;
				return 2.0 * k1 * k2 * S / (k1 * cell2.hz + k2 * cell1.hz);
			} else {
				k1 = (cells[idx1].r > r_eff[idx1 % (cellsNum_z+2)] ? props_sk.perm_r[idx1 % (cellsNum_z+2)] : Perm_eff[idx1 % (cellsNum_z+2)]);
				k2 = (cells[idx2].r > r_eff[idx2 % (cellsNum_z+2)] ? props_sk.perm_r[idx2 % (cellsNum_z+2)] : Perm_eff[idx2 % (cellsNum_z+2)]);
				S = 2.0 * M_PI * cell1.hz * (cell1.r + sign(idx2 - idx1) * cell1.hr / 2.0);
				return 2.0 * k1 * k2 * S / (k1 * cell2.hr + k2 * cell1.hr);
			}
		};
		inline double getPerm_r(double r, int idx)
		{
			return (r > r_eff[idx % (cellsNum_z+2)] ? props_sk.perm_r[idx % (cellsNum_z+2)] : Perm_eff[idx % (cellsNum_z+2)]);
		};
		inline double getKr_oil(double sat_oil)
		{
			if(sat_oil > 1.0)
				return 1.0;
			else
				return props_oil.kr->Solve(sat_oil);
		};
		inline double getKr_oil_ds(double sat_oil)
		{
			if(sat_oil > 1.0)
				return props_oil.kr->DSolve(1.0);
			else
				return props_oil.kr->DSolve(sat_oil);
		};
		inline double getKr_gas(double sat_oil)
		{
			if(sat_oil > 1.0)
				return 0.0;
			else
				return props_gas.kr->Solve(sat_oil);
		};
		inline double getKr_gas_ds(double sat_oil)
		{
			if(sat_oil > 1.0)
				return props_gas.kr->DSolve(1.0);
			else
				return props_gas.kr->DSolve(sat_oil);
		};
		inline double getB_oil(double p, double p_bub, bool SATUR)
		{
			if(SATUR)
				return props_oil.b->Solve(p);
			else
				return props_oil.b->Solve(p_bub) * (1.0 + props_oil.beta * (p_bub - p));
		};
		inline double getBoreB_oil(double p, double p_bub, bool SATUR)
		{
			//return props_oil.b_bore;
			return getB_oil(p, p_bub, SATUR);
		};
		inline double getB_oil_dp(double p, double p_bub, bool SATUR)
		{
			if(SATUR)
				return props_oil.b->DSolve(p);
			else
				return -props_oil.b->Solve(p_bub) * props_oil.beta;
		};
		inline double getB_gas(double p)
		{
			return props_gas.b->Solve(p);
		};
		inline double getB_gas_dp(double p)
		{
			return props_gas.b->DSolve(p);
		};
		inline double getRs(double p, double p_bub, bool SATUR)
		{
			if(SATUR)
				return Rs->Solve(p);
			else
				return Rs->Solve(p_bub);
		};
		inline double getRs_dp(double p, double p_bub, bool SATUR)
		{
			if(SATUR)
				return Rs->DSolve(p);
			else
				return 0.0;
		};
		inline double getPresFromRs(double rs)
		{
			return Prs->Solve(rs);
		};
		inline void solveP_bub()
		{
			int idx;
			double factRs, dissGas;

			for(int i = 0; i < cellsNum_r+2; i++)
				for(int j = 0; j < cellsNum_z+2; j++)
				{
					idx = i * (cellsNum_z + 2) + j;

					Var2phaseNIT& next = cells[idx].u_next;
					Var2phaseNIT& prev = cells[idx].u_prev;

					if(next.s > 1.0)
						next.s = 1.0;

					dissGas = (1.0 - next.s) * getB_oil(next.p, next.p_bub, next.SATUR) / ( (1.0 - next.s) * getB_oil(next.p, next.p_bub, next.SATUR) + next.s * getB_gas(next.p));
					factRs = getRs(prev.p, prev.p_bub, next.SATUR) + dissGas;

					if(getRs(next.p, next.p, next.SATUR) > factRs)
					{
						next.p_bub = getPresFromRs(factRs);
						next.SATUR = false;
					} else {
						next.p_bub = next.p;
						next.SATUR = true;
					}
				}
		};

		// Thermal functions
		inline double getRho_oil(double p, double p_bub, bool SATUR)
		{
			return (props_oil.dens_stc + getRs(p, p_bub, SATUR) * props_gas.dens_stc) / getB_oil(p, p_bub, SATUR);
		};
		inline double getRho_gas(double p)
		{
			return props_gas.dens_stc / getB_gas(p);
		};
		inline double getNablaP(Cell& cell, int varNum, int axis)
		{
			Cell* nebr1;
			Cell* nebr2;
			double h;

			switch(axis)
			{
			case R_AXIS:
				nebr1 = &cells[cell.num - cellsNum_z - 2];
				nebr2 = &cells[cell.num + cellsNum_z + 2];
				h = nebr2->r - nebr1->r;
				break;
			case Z_AXIS:
				nebr1 = &cells[cell.num - 1];
				nebr2 = &cells[cell.num + 1];
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
			Var2phaseNIT* var;
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
				return -getPerm_r(cell.r, cell.num) * getKr_oil(var->s) / props_oil.visc * getNablaP(cell, varNum, axis);
			case Z_AXIS:
				return -props_sk.perm_z[cell.num % (cellsNum_z+2)] * getKr_oil(var->s) / props_oil.visc * getNablaP(cell, varNum, axis);
			}
		};
		inline double getGasVelocity(Cell& cell, int varNum, int axis)
		{
			Var2phaseNIT* var;
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
				return -getPerm_r(cell.r, cell.num) * getKr_gas(var->s) / props_gas.visc * getNablaP(cell, varNum, axis);
			case Z_AXIS:
				return -props_sk.perm_z[cell.num % (cellsNum_z+2)] * getKr_gas(var->s) / props_gas.visc * getNablaP(cell, varNum, axis);
			}
		};
		inline double getCn(Var2phaseNIT& var)
		{
			return getPoro(var.p) * (var.s * getRho_oil(var.p, var.p_bub, var.SATUR) * props_oil.c +
						(1.0 - var.s) * getRho_gas(var.p) * props_gas.c ) + 
						(1.0 - getPoro(var.p)) * props_sk.dens_stc * props_sk.c;
		};
		inline double getAd(Var2phaseNIT& var)
		{
			return getPoro(var.p) * (var.s * getRho_oil(var.p, var.p_bub, var.SATUR) * props_oil.c * props_oil.ad +
								(1.0 - var.s) * getRho_gas(var.p) * props_gas.c * props_gas.ad );
		};
		inline double getLambda(Var2phaseNIT& var, int axis)
		{
			switch(axis)
			{
			case R_AXIS:
				return getPoro(var.p) * (var.s * props_oil.lambda + (1.0-var.s) * props_gas.lambda) + 
					(1.0-getPoro(var.p)) * props_sk.lambda_r;
			case Z_AXIS:
				return getPoro(var.p) * (var.s * props_oil.lambda + (1.0-var.s) * props_gas.lambda) + 
					(1.0-getPoro(var.p)) * props_sk.lambda_z;
			}
		};
		inline double getLambdaWorse(Var2phaseNIT& var, int axis)
		{
			switch(axis)
			{
			case R_AXIS:
				return getPoro(var.p) * (var.s * props_oil.lambda + (1.0-var.s) * props_gas.lambda) + 
					(1.0-getPoro(var.p)) * props_sk.lambda_r * 4.0;
			case Z_AXIS:
				return getPoro(var.p) * (var.s * props_oil.lambda + (1.0-var.s) * props_gas.lambda) + 
					(1.0-getPoro(var.p)) * props_sk.lambda_z * 4.0;
			}
		};
		inline double getLambda(Cell& cell1, Cell& cell2)
		{
			double cm;
			if( abs(cell1.num - cell2.num) == 1) {
				cm = cell1.r;
				if(cm > 0.2524 / R_dim)
					return ( cell1.hz * getLambda(cell2.u_next, Z_AXIS) + cell2.hz * getLambda(cell1.u_next, Z_AXIS) ) / (cell1.hz + cell2.hz);
				else 
					return ( cell1.hz * getLambdaWorse(cell2.u_next, Z_AXIS) + cell2.hz * getLambdaWorse(cell1.u_next, Z_AXIS) ) / (cell1.hz + cell2.hz);
			} else {
				cm = (cell1.r + sign(cell2.num - cell1.num) * cell1.hr / 2.0);
				if(cm > 0.2524 / R_dim)
					return ( cell1.hr * getLambda(cell2.u_next, R_AXIS) + cell2.hr * getLambda(cell1.u_next, R_AXIS) ) / (cell1.hr + cell2.hr);
				else 
					return ( cell1.hr * getLambdaWorse(cell2.u_next, R_AXIS) + cell2.hr * getLambdaWorse(cell1.u_next, R_AXIS) ) / (cell1.hr + cell2.hr);
			}
		};
		inline double getJT(Cell& cell, int varNum, int axis)
		{
			Var2phaseNIT* var;
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
		
			return getRho_oil(var->p, var->p_bub, var->SATUR) * props_oil.c * props_oil.jt * getOilVelocity(cell, varNum, axis) + 
					getRho_gas(var->p) * props_gas.c * props_gas.jt * getGasVelocity(cell, varNum, axis);
		};
		inline double getA(Cell& cell, int varNum, int axis)
		{
			Var2phaseNIT* var;
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
		
			return getRho_oil(var->p, var->p_bub, var->SATUR) * props_oil.c * getOilVelocity(cell, varNum, axis) + 
					getRho_gas(var->p) * props_gas.c * getGasVelocity(cell, varNum, axis);
		};

		// First eqn
		double solve_eq1(int cur);
		double solve_eq1_dp(int cur);
		double solve_eq1_ds(int cur);
		double solve_eq1_dp_beta(int cur, int beta);
		double solve_eq1_ds_beta(int cur, int beta);

		// Second eqn
		double solve_eq2(int cur);
		double solve_eq2_dp(int cur);
		double solve_eq2_ds(int cur);
		double solve_eq2_dp_beta(int cur, int beta);
		double solve_eq2_ds_beta(int cur, int beta);

		// Left boundary condition
		double solve_eqLeft(int cur);
		double solve_eqLeft_dp(int cur);
		double solve_eqLeft_ds(int cur);
		double solve_eqLeft_dp_beta(int cur);
		double solve_eqLeft_ds_beta(int cur);

		// Solves intensity of phase transition
		double solve_eq3(int cur);

		Snapshotter<GasOil_RZ_NIT>* snapshotter;
		void snapshot(int i);
		void snapshot_all(int i);

	public:
		GasOil_RZ_NIT();
		~GasOil_RZ_NIT();
	
		void setSnapshotter(std::string type);
		void setPeriod(int period);
	};
};

#endif /* GASOIL_RZ_NIT_H_ */