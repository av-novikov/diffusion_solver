#ifndef GAS1D_H_
#define GAS1D_H_

#include <vector>
#include <map>

#include "util/utils.h"
#include "model/cells/Variables.hpp"
#include "model/cells/RadialCell.hpp"
#include "model/AbstractModel.hpp"

namespace gas1D
{
	typedef RadialCell<Var1phase> Cell;

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
	};

	struct Gas_Props
	{
		// Viscosity [cP]
		double visc;
		// Density of fluid in STC [kg/m3]
		double dens_stc;
		// Z-factor ( supercompressibility factor )
		Interpolate* z;
		// Volume factor for well pipe
		double b_bore;
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
		Gas_Props props_gas;

		double depth_point;

		// Data set (saturation, relative oil permeability)
		std::vector< std::pair<double,double> > z_factor;
	};

	class Gas1D : public AbstractModel<Var1phase, Properties, RadialCell>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class GRDECLSnapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class AbstractMethod;
		friend class Gas1DSolver;	

	protected:
		// Continuum properties
		int skeletonsNum;
		std::vector<Skeleton_Props> props_sk;
		Gas_Props props_gas;

		// Number of cells in radial direction
		int cellsNum_r;
		// Number of cells in vertical direction
		int cellsNum_z;

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
		// Set initial state
		void setInitialState();

		inline double getPdivZ(double p) const
		{
			return p / props_gas.z->Solve(p);
		};
		inline double getPdivZ(const Cell& cell1, const Cell& cell2) const
		{
			return ( getPdivZ(cell1.u_next.p) * cell2.hr + getPdivZ(cell2.u_next.p) * cell1.hr ) / (cell1.hr + cell2.hr);
		};
		inline double getPdivZ_dp(double p) const
		{
			const double z = props_gas.z->Solve(p);
			return (1.0 - p / z * props_gas.z->DSolve(p)) / z;
		};
		inline double getPdivZ_dp(const Cell& cell, const Cell& beta) const
		{
			return ( getPdivZ_dp(cell.u_next.p) * beta.hr ) / (cell.hr + beta.hr);
		};
		inline double getPdivZ_dp_beta(const Cell& cell, const Cell& beta) const
		{
			return ( getPdivZ_dp(beta.u_next.p) * cell.hr ) / (cell.hr + beta.hr);
		};
		inline double getTrans(Cell& cell, Cell& beta) const
		{
			double k1, k2;
			k1 = (cell.r > props_sk[0].radius_eff ? props_sk[0].perm_r : props_sk[0].perm_eff);
			k2 = (beta.r > props_sk[0].radius_eff ? props_sk[0].perm_r : props_sk[0].perm_eff);
			return 4.0 * M_PI * k1 * k2 * (cell.r + sign(beta.num - cell.num) * cell.hr / 2.0) * props_sk[0].height / (k2 * cell.hr + k1 * beta.hr);
		};
		inline double getBoreB_gas(double p) const
		{
			return props_gas.b_bore;
			//return 101325.0 / P_dim / p;
		};

		Snapshotter<Gas1D>* snapshotter;
		void snapshot(int i);
		void snapshot_all(int i);

		double solve_eq(int i);
		double solve_eq_dp(int i);
		double solve_eq_dp_beta(int i, int beta);
		
		double solve_eqLeft();
		double solve_eqLeft_dp();
		double solve_eqLeft_dp_beta();
		
		double solve_eqRight();
		double solve_eqRight_dp();
		double solve_eqRight_dp_beta();

	public:
		Gas1D();
		~Gas1D();
	
		void setSnapshotter(std::string type);
		void setPeriod(int period);

	};
};

#endif /* GAS1D_H_ */