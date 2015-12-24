#ifndef GAS1D_SIMPLE_H_
#define GAS1D_SIMPLE_H_

#include <vector>
#include <map>

#include "util/utils.h"
#include "model/cells/Variables.hpp"
#include "model/cells/RadialCell.hpp"
#include "model/AbstractModel.hpp"
#include "model/Gas1D/Gas1D.h"

namespace gas1D
{
	class Gas1D_simple : public AbstractModel<Var1phase, Properties, RadialCell>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class GRDECLSnapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class AbstractMethod;
		template<typename> friend class Gas1DSolver;	

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

		inline double getZ(double p)
		{
			return props_gas.z->Solve(p);
		};
		inline double getZ_dp(double p)
		{
			return props_gas.z->DSolve(p);
		};
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

		Snapshotter<Gas1D_simple>* snapshotter;
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
		Gas1D_simple();
		~Gas1D_simple();
	
		void setSnapshotter(std::string type);
		void setPeriod(int period);

		double getRate();
	};
};

#endif /* GAS1D_SIMPLE_H_ */