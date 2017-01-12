#ifndef GASOILELLIPTIC_HPP_
#define GASOILELLIPTIC_HPP_

#include <vector>
#include <map>
#include <string>
#include <initializer_list>

#define ADOLC_ADVANCED_BRANCHING

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

		inline const Cell& getUpwindIdx(const Cell& cell, const Cell& nebr) const
		{
			assert(cell.isUsed);
			assert(nebr.isUsed);
			if (cell.u_next.p < nebr.u_next.p)
				return nebr;
			else
				return cell;
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
		inline Cell& getCell(const int num)
		{
			return cells[num];
		};
		inline Cell& getCell(const int num, const int beta)
		{
			Cell& nebr = cells[beta];
			if (nebr.isUsed)
				return nebr;
			else
				return wellCells[wellNebrMap.at(num)];
		};
		inline const Cell& getCell(const int num) const
		{
			return cells[num];
		};
		inline const Cell& getCell(const int num, const int beta) const
		{
			const Cell& nebr = cells[beta];
			if (nebr.isUsed)
				return nebr;
			else
				return wellCells[wellNebrMap.at(num)];
		};
		inline void getNeighbors(const Cell& cell, Cell** const neighbor)
		{
			// FIXME: Only for inner cells
			assert(cell.isUsed);

			neighbor[0] = &getCell(cell.num, cell.num - cellsNum_z - 2);
			neighbor[1] = &getCell(cell.num, cell.num + cellsNum_z + 2);
			neighbor[2] = &getCell(cell.num, cell.num - 1);
			neighbor[3] = &getCell(cell.num, cell.num + 1);

			if (cell.num < (cellsNum_mu + 2) * (cellsNum_z + 2))
				neighbor[4] = &getCell(cell.num, cell.num + (cellsNum_mu + 2) * (cellsNum_z + 2) * (cellsNum_nu - 1));
			else
				neighbor[4] = &getCell(cell.num, cell.num - (cellsNum_mu + 2) * (cellsNum_z + 2));
			if (cell.num < (cellsNum_mu + 2) * (cellsNum_z + 2) * (cellsNum_nu - 1))
				neighbor[5] = &getCell(cell.num, cell.num + (cellsNum_mu + 2) * (cellsNum_z + 2));
			else
				neighbor[5] = &getCell(cell.num, cell.num - (cellsNum_mu + 2) * (cellsNum_z + 2) * (cellsNum_nu - 1));
		};
		inline void getStencil(const Cell& cell, Cell** const neighbor)
		{
			// FIXME: Only for inner cells
			assert(cell.isUsed);

			neighbor[0] = &getCell(cell.num);
			neighbor[1] = &getCell(cell.num, cell.num - cellsNum_z - 2);
			neighbor[2] = &getCell(cell.num, cell.num + cellsNum_z + 2);
			neighbor[3] = &getCell(cell.num, cell.num - 1);
			neighbor[4] = &getCell(cell.num, cell.num + 1);

			if (cell.num < (cellsNum_mu + 2) * (cellsNum_z + 2))
				neighbor[5] = &getCell(cell.num, cell.num + (cellsNum_mu + 2) * (cellsNum_z + 2) * (cellsNum_nu - 1));
			else
				neighbor[5] = &getCell(cell.num, cell.num - (cellsNum_mu + 2) * (cellsNum_z + 2));
			if (cell.num < (cellsNum_mu + 2) * (cellsNum_z + 2) * (cellsNum_nu - 1))
				neighbor[6] = &getCell(cell.num, cell.num + (cellsNum_mu + 2) * (cellsNum_z + 2));
			else
				neighbor[6] = &getCell(cell.num, cell.num - (cellsNum_mu + 2) * (cellsNum_z + 2) * (cellsNum_nu - 1));
		};
		inline double getTrans(const Cell& cell, const Cell& beta) const
		{
			double k1, k2, S;

			if (fabs(cell.z - beta.z) > EQUALITY_TOLERANCE) {
				const double dz1 = fabs(cell.z - sk_well->h_well);
				const double dz2 = fabs(beta.z - sk_well->h_well);
				k1 = (dz1 > cell.props->radius_eff_z ? cell.props->perm_z : cell.props->perm_eff_z);
				k2 = (dz2 > beta.props->radius_eff_z ? beta.props->perm_z : beta.props->perm_eff_z);
				if (k1 == 0.0 && k2 == 0.0)
					return 0.0;
				S = Cell::getH(cell.mu, cell.nu) * Cell::getH(cell.mu, cell.nu) * cell.hmu * cell.hnu;
				return 2.0 * k1 * k2 * S / (k1 * beta.hz + k2 * cell.hz);
			}
			else if (fabs(cell.mu - beta.mu) > EQUALITY_TOLERANCE) {
				k1 = (cell.mu > cell.props->radius_eff_mu ? cell.props->perm_mu : cell.props->perm_eff_mu);
				k2 = (beta.mu > beta.props->radius_eff_mu ? beta.props->perm_mu : beta.props->perm_eff_mu);
				if (k1 == 0.0 && k2 == 0.0)
					return 0.0;

				const double mu = cell.mu + sign(beta.mu - cell.mu) * cell.hmu / 2.0;
				S = Cell::getH(mu, cell.nu) * cell.hnu * cell.hz;
				return 2.0 * k1 * k2 * S /
					( k1 * beta.hmu * Cell::getH(beta.mu, beta.nu)
					+ k2 * cell.hmu * Cell::getH(cell.mu, cell.nu) );
			}
			else if (fabs(cell.nu - beta.nu) > EQUALITY_TOLERANCE) {
				k1 = (cell.mu > cell.props->radius_eff_mu ? cell.props->perm_mu : cell.props->perm_eff_mu);
				if (k1 == 0.0)
					return 0.0;
				const double nu = cell.nu + sign(beta.nu - cell.nu) * cell.hnu / 2.0;
				S = cell.hz * cell.hmu * Cell::getH(cell.mu, nu);
				return 2.0 * k1 * S /
					( beta.hnu * Cell::getH(beta.mu, beta.nu) + 
					+ cell.hnu * Cell::getH(cell.mu, cell.mu) );
			}
		};

		//void setVariables(int cur);
		void solve_eqMiddle(int cur);
		//void solve_eqLeft(int cur);
		//void solve_eqRight(int cur);
		//void solve_eqVertical(int cur);

	public:
		GasOil_Elliptic();
		~GasOil_Elliptic();

		double* x;
		double* y;

		void setPeriod(int period);
	};
};


#endif /* GASOILELLIPTIC_HPP_ */