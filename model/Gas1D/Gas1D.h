#ifndef GAS1D_H_
#define GAS1D_H_

#include <vector>
#include <map>

#include "util/utils.h"
#include "model/cells/Variables.hpp"
#include "model/cells/RadialCell.hpp"
#include "model/AbstractModel.hpp"

namespace gas1d
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
		int skeletonsNum = 1;
		std::vector<Skeleton_Props> props_sk;
		Gas_Props props_gas;

		// Number of cells in radial direction
		int cellsNum_r;
		// Number of cells in vertical direction
		int cellsNum_z = 1;

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

	public:
		Gas1D();
		~Gas1D();
	
		void setSnapshotter(std::string type);
		void setPeriod(int period);

	};
};

#endif /* GAS1D_H_ */