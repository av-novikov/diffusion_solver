#ifndef GASOILELLIPTIC_HPP_
#define GASOILELLIPTIC_HPP_

#include <vector>
#include <map>
#include <string>

#include "model/cells/Variables.hpp"
#include "model/GasOil_Elliptic/Properties.hpp"
#include "model/cells/EllipticCell.hpp"
#include "model/AbstractModel.hpp"
#include "boost/math/special_functions/ellint_2.hpp"

namespace gasOil_elliptic
{
	static const int stencil = 7;
	static const int Lstencil = 3;
	static const int Rstencil = 2;
	static const int Vstencil = 2;

	typedef Var2phase Variable;
	typedef TapeVarGasOil TapeVariable;
	typedef EllipticCell<Variable, Skeleton_Props> Cell;
	typedef Cell::Point Point;
	template <typename TVariable> using TCell = EllipticCell<TVariable, Skeleton_Props>;

	class GasOil_Elliptic : public AbstractModel<Variable, Properties, TCell, GasOil_Elliptic>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class GRDECLSnapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class AbstractMethod;
		template<typename> friend class AbstractSolver;
		friend class GasOilEllipticSolver;
	protected:

		// Continuum properties
		int skeletonsNum;
		std::vector<Skeleton_Props> props_sk;
		std::vector<Skeleton_Props>::iterator sk_well;
		Oil_Props props_oil;
		Gas_Props props_gas;

		double l;
		int cellsNum_mu;
		int cellsNum_nu;
		int cellsNum_z;

		void setInitialState();
		void setProps(Properties& props);
		void makeDimLess();
		void buildGridLog();
		void setPerforated();
		void setUnused();
		void buildWellCells();

		std::vector<Cell> wellCells;
		std::map<int,int> wellNebrMap;

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
		inline double getTrans(Cell& cell, Cell& beta) const
		{
			double k1, k2, S;

			if (fabs(cell.z - beta.z) > EQUALITY_TOLERANCE) {
				const double dz1 = fabs(cell.z - sk_well->h_well);
				const double dz2 = fabs(beta.z - sk_well->h_well);
				k1 = (dz1 > cell.props->radius_eff_z ? cell.props->perm_z : cell.props->perm_eff_z);
				k2 = (dz2 > beta.props->radius_eff_z ? beta.props->perm_z : beta.props->perm_eff_z);
				if (k1 == 0.0 && k2 == 0.0)
					return 0.0;
				S = Cell::a * Cell::a * (sinh(cell.mu) * sinh(cell.mu) + sin(cell.nu) * sin(cell.nu)) * cell.hmu * cell.hnu;
				return 2.0 * k1 * k2 * S / (k1 * beta.hz + k2 * cell.hz);
			}
			else if (fabs(cell.mu - beta.mu) > EQUALITY_TOLERANCE) {
				k1 = (cell.mu > cell.props->radius_eff_mu ? cell.props->perm_mu : cell.props->perm_eff_mu);
				k2 = (beta.mu > beta.props->radius_eff_mu ? beta.props->perm_mu : beta.props->perm_eff_mu);
				if (k1 == 0.0 && k2 == 0.0)
					return 0.0;
				const double eccentr = 1.0 / cosh(cell.mu + sign(beta.mu - cell.mu) * cell.hmu / 2.0);
				S = Cell::a * cosh(cell.mu + sign(beta.mu - cell.mu) * cell.hmu / 2.0) * cell.hz * (boost::math::ellint_2(eccentr, cell.nu + cell.hnu / 2.0) - boost::math::ellint_2(eccentr, cell.nu - cell.hnu / 2.0));
				// S = Cell::a * cosh(cell.mu + sign(beta.mu - cell.mu) * cell.hmu / 2.0) * cell.hz * sqrt(1.0 - eccentr * eccentr * cos(cell.nu) * cos(cell.nu)) * cell.hnu;
				return 2.0 * k1 * k2 * S / Cell::a / 
					( k1 * beta.hmu * sqrt(sinh(beta.mu) * sinh(beta.mu) * cos(beta.nu) * cos(beta.nu) + cosh(beta.mu) * cosh(beta.mu) * sin(beta.nu) * sin(beta.nu))
					+ k2 * cell.hmu * sqrt(sinh(cell.mu) * sinh(cell.mu) * cos(cell.nu) * cos(cell.nu) + cosh(cell.mu) * cosh(cell.mu) * sin(cell.nu) * sin(cell.nu)) );
			}
			else if (fabs(cell.nu - beta.nu) > EQUALITY_TOLERANCE) {
				k1 = (cell.mu > cell.props->radius_eff_mu ? cell.props->perm_mu : cell.props->perm_eff_mu);
				if (k1 == 0.0)
					return 0.0;
				const double nu = cell.nu + sign(beta.nu - cell.nu) * cell.hnu / 2.0;
				S = cell.hz * Cell::a * cell.hmu * sqrt(sinh(cell.mu) * sinh(cell.mu) * cos(nu) * cos(nu) + cosh(cell.mu) * cosh(cell.mu) * sin(nu) * sin(nu));
				return 2.0 * k1 * S / Cell::a / 
					( beta.hnu * sqrt(sinh(beta.mu) * sinh(beta.mu) * cos(beta.nu) * cos(beta.nu) + cosh(beta.mu) * cosh(beta.mu) * sin(beta.nu) * sin(beta.nu)) + 
					+ cell.hnu * sqrt(sinh(cell.mu) * sinh(cell.mu) * cos(cell.nu) * cos(cell.nu) + cosh(cell.mu) * cosh(cell.mu) * sin(cell.nu) * sin(cell.nu)));
			}
		};
	public:
		GasOil_Elliptic();
		~GasOil_Elliptic();

		double* x;
		double* y;

		void setPeriod(int period);
	};
};


#endif /* GASOILELLIPTIC_HPP_ */
