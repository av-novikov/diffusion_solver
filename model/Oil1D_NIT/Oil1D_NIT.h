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

	struct Properties
	{
		std::vector<double> timePeriods;
		std::vector<double> rates;
		std::vector<double> skins;
		std::vector<double> radius;

		std::vector<std::pair<int,int> > perfIntervals;

		double ht;
		double ht_min;
		double ht_max;

		double alpha;

		int cellsNum_r;

		double r_w;
		double r_e;
		double height;
		double m;
		double perm;
		double dens_sk_stc;
		double beta_sk;

		double visc_oil;
		double dens_oil_stc;
		double beta_oil;

		double p_init;
	};

	struct Skeleton_Props
	{
		double m;
		double beta;
		double dens_stc;

		double perm;

		std::vector<double> perm_eff;
		std::vector<double> r_eff;
		std::vector<double> skin;

		double height;
	};

	struct Oil_Props
	{
		double visc;
		double dens_stc;
		double beta;
	};

	class Oil1D_NIT : public AbstractModel<Var1phaseNIT, Properties, RadialCell>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class GRDECLSnapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class AbstractMethod;
		friend class Oil1DNITSolver;	
	protected:
		Skeleton_Props props_sk;
		Oil_Props props_oil;

		double r_eff;
		double Perm_eff;
		double skin;
		
		int cellsNum_r;

		double alpha;

		void setProps(Properties& props);
		void makeDimLess();
		void buildGridLog();
		void setPerforated();
		
		// Service functions
		inline double upwindIsCur(int cur, int beta)
		{
			if(cells[cur].u_next.p < cells[beta].u_next.p)
				return 0.0;
			else
				return 1.0;
		};
		inline int getUpwindIdx(int cur, int beta)
		{
			if(cells[cur].u_next.p < cells[beta].u_next.p)
				return beta;
			else
				return cur;
		};

		// Solving coefficients
		inline double getPoro(double p)
		{
			return props_sk.m * (1.0 + props_sk.beta * (varInit.p - p) );
		};
		inline double getRho(double p)
		{
			return props_oil.dens_stc * (1.0 + props_oil.beta * p);
		};
		inline double getTrans(Cell& cell, Cell& beta)
		{
			double k1, k2;
			k1 = (cell.r > r_eff ? props_sk.perm : Perm_eff);
			k2 = (beta.r > r_eff ? props_sk.perm : Perm_eff);
			return 4.0 * M_PI * k1 * k2 * (cell.r + sign(beta.num - cell.num) * cell.hr / 2.0) * props_sk.height / (k2 * cell.hr + k1 * beta.hr);
		};
		inline double getRho(Cell& cell, Cell& beta)
		{
			return ( getRho(beta.u_next.p) * cell.hr + beta.hr * getRho(cell.u_next.p) ) / (beta.hr + cell.hr);
		};
		inline double getPwf()
		{
			return varInit.p - Qcell[0] * props_oil.visc * log(r_e / r_w) / 2.0 / M_PI / props_sk.height / props_sk.perm;
		};

		Snapshotter<Oil1D_NIT>* snapshotter;
		void snapshot(int i);
		void snapshot_all(int i);

		double solve_eq(int i);
		double solve_eq_dp(int i);
		double solve_eq_dp_beta(int i, int beta);
		
		double solve_left();
		double solve_left_dp();
		double solve_left_dp_beta();
		
		double solve_right();
		double solve_right_dp_beta();
		double solve_right_dp();

	public:
		Oil1D_NIT();
		~Oil1D_NIT();

		void setPeriod(int period);
		void setSnapshotter(std::string type);
	};

};

#endif OIL1D_NIT_H_

