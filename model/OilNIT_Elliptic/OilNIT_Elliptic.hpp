#ifndef OILNITELLIPTIC_HPP_
#define OILNITELLIPTIC_HPP_

#include <vector>
#include <map>
#include <string>

#define ADOLC_ADVANCED_BRANCHING

#include "model/cells/Variables.hpp"
#include "model/OilNIT_Elliptic/Properties.hpp"
#include "model/cells/EllipticCell.hpp"
#include "model/AbstractModel.hpp"

#include <cassert>

namespace oilnit_elliptic
{
	static const int MU_AXIS = 0;
	static const int NU_AXIS = 1;
	static const int Z_AXIS = 2;

	static const int stencil = 7;
	static const int Lstencil = 2;
	static const int TLstencil = 3;
	static const int Rstencil = 2;
	static const int Vstencil = 2;

	typedef Var1phaseNIT Variable;
	typedef TapeVar1Phase TapeVariable;
	typedef TapeVar1PhaseNIT TapeVariableNIT;
	typedef EllipticCell<Variable, Skeleton_Props> Cell;
	typedef Cell::Point Point;
	typedef Cell::Type Type;
	template <typename TVariable> using TCell = EllipticCell<TVariable, Skeleton_Props>;

	class OilNIT_Elliptic : public AbstractModel<Variable, Properties, TCell, OilNIT_Elliptic>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class GRDECLSnapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class AbstractMethod;
		template<typename> friend class AbstractSolver;
		template<typename> friend class OilNITEllipticSolver;
	protected:

		// Continuum properties
		int skeletonsNum;
		std::vector<Skeleton_Props> props_sk;
		std::vector<Skeleton_Props>::iterator sk_well;
		Oil_Props props_oil;

		double l;
		double mu_init;
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
		void setRateDeviation(int num, double ratio);
		double solveH();

		double Q_sum_quater;
		std::map<int, double> Qcell_ellipse;
		std::vector<Cell> wellCells;
		std::map<int, int> wellNebrMap;
		std::map<int, std::pair<int, int> > nebrMap;

		inline const std::vector<int> getPerforationIndices(const int idx) const
		{
			const Cell& cell = wellCells[idx];
			if (cell.type == Type::WELL_LAT)
				if (cell.nu == 0.0 || cell.nu == M_PI)
					return{ idx, idx + cellsNum_nu, idx + 2 * cellsNum_nu };
				else
					return{ idx, cellsNum_nu - idx, cellsNum_nu + idx,
							2 * cellsNum_nu - idx, 2 * cellsNum_nu + idx, 3 * cellsNum_nu - idx };
			else if (cell.type == Type::WELL_TOP)
				return{ idx - cellsNum_nu, 2 * cellsNum_nu - idx,
						idx, 3 * cellsNum_nu - idx, cellsNum_nu + idx, 4 * cellsNum_nu - idx };
		};
		inline const std::vector<int> getPerforationIndicesSameType(const int idx) const 
		{
			const Cell& cell = wellCells[idx];
			
			if (cell.type == Type::WELL_LAT)
				if (cell.nu == 0.0 || cell.nu == M_PI)
					return{ 0, cellsNum_nu / 2};
				else
					return{ idx % (cellsNum_nu / 2), cellsNum_nu - (idx % (cellsNum_nu / 2))};
			else
			{
				const int startIdx = cellsNum_nu + (idx % (cellsNum_nu / 2));
				return{ startIdx, 3 * cellsNum_nu - startIdx, cellsNum_nu + startIdx, 4 * cellsNum_nu - startIdx };
			}
		};
		inline const std::vector<int> getSymmetricalWellIndices(const int idx) const
		{
			std::vector<int> indices;
			const Cell& cell = wellCells[idx];
			if (cell.type == Type::WELL_LAT)
			{
				if (cell.nu == 0.0 || cell.nu == M_PI_2)
				{
					indices.resize(2);
					indices[0] = idx;
					indices[1] = idx + cellsNum_nu / 2;
				}
				else
				{
					indices.resize(4);
					indices[0] = idx;
					indices[1] = cellsNum_nu - idx;
					indices[2] = cellsNum_nu / 2 + idx;
					indices[3] = cellsNum_nu / 2 - idx;
				}
			}
			else if (cell.type == Type::WELL_TOP)
			{
				const int topStartIdx = cellsNum_nu;
				const int botStartIdx = 2 * cellsNum_nu;

				if (cell.nu == 0.0 || cell.nu == M_PI_2)
				{
					indices.resize(4);
					indices[0] = topStartIdx + idx;
					indices[1] = topStartIdx + idx + cellsNum_nu / 2;

					indices[2] = botStartIdx + idx;
					indices[3] = botStartIdx + idx + cellsNum_nu / 2;
				}
				else
				{
					indices.resize(8);
					indices[0] = topStartIdx + idx;
					indices[1] = topStartIdx + cellsNum_nu - idx;
					indices[2] = topStartIdx + cellsNum_nu / 2 + idx;
					indices[3] = topStartIdx + cellsNum_nu / 2 - idx;

					indices[4] = botStartIdx + idx;
					indices[5] = botStartIdx + cellsNum_nu - idx;
					indices[6] = botStartIdx + cellsNum_nu / 2 + idx;
					indices[7] = botStartIdx + cellsNum_nu / 2 - idx;
				}
			}

			return indices;
		};
		inline const Cell& getUpwindCell(const Cell& cell, const Cell& nebr) const
		{
			assert(cell.isUsed);
			assert(nebr.isUsed);
			if (cell.u_next.p < nebr.u_next.p)
				return nebr;
			else
				return cell;
		};
		inline const int getSkeletonIdx(const Cell& cell) const
		{
			int idx = 0;
			/*while (idx < props_sk.size())
			{
				if (cell.z <= props_sk[idx].h2 + EQUALITY_TOLERANCE)
					return idx;
				idx++;
			}*/
			const double nu = (cell.nu >= M_PI) ? 2.0 * M_PI - cell.nu : cell.nu;
			return int(5.0 * (nu - EQUALITY_TOLERANCE) / M_PI);
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
		inline int getCellIdx(const int num, const int beta)
		{
			Cell& nebr = cells[beta];
			if (nebr.isUsed)
				return beta;
			else
				return cellsNum + wellNebrMap.at(num);
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
		inline double getDamagedPerm(const Skeleton_Props& props, const int i) const
		{
			const double R = log((cosh(M_PI * l / props.height) - cos(2.0 * M_PI * props.h_well / props.height)) / 
								(cosh(M_PI * r_w / props.height) - cos(2.0 * M_PI * props.h_well / props.height)) * 
								(cosh(M_PI * l / props.height) - 1.0) / (cosh(M_PI * r_w / props.height) - 1.0));
			return props.perm_mu / (props.skins[i] / 0.08 / R / log(props.radiuses_eff[i] / r_w) + 1.0);
		};
		inline void getNeighbors(const Cell& cell, Cell** const neighbor)
		{
			// FIXME: Only for inner cells
			assert(cell.isUsed);

			switch (cell.type)
			{
			case Type::MIDDLE:
				// Special neighbor search for center cells
				if (cell.num % ((cellsNum_mu + 2) * (cellsNum_z + 2)) > cellsNum_z + 1)
					neighbor[0] = &getCell(cell.num, cell.num - cellsNum_z - 2);
				else
				{
					int nu_idx = cell.num / ((cellsNum_z + 2) * (cellsNum_mu + 2));
					neighbor[0] = &getCell(cell.num, cellsNum + cell.num - 
						2 * nu_idx * (cellsNum_z + 2) * (cellsNum_mu + 2));
				}
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
				break;

			case Type::MIDDLE_SIDE:
				neighbor[0] = &getCell(cell.num, cell.num + cellsNum_z + 2);
				neighbor[1] = &getCell(cell.num, cell.num - 1);
				neighbor[2] = &getCell(cell.num, cell.num + 1);

				if (cell.num < (cellsNum_mu + 2) * (cellsNum_z + 2))
					neighbor[3] = &getCell(cell.num, cell.num + (cellsNum_mu + 2) * (cellsNum_z + 2) * (cellsNum_nu - 1));
				else
					neighbor[3] = &getCell(cell.num, cell.num - (cellsNum_mu + 2) * (cellsNum_z + 2));
				if (cell.num < (cellsNum_mu + 2) * (cellsNum_z + 2) * (cellsNum_nu - 1))
					neighbor[4] = &getCell(cell.num, cell.num + (cellsNum_mu + 2) * (cellsNum_z + 2));
				else
					neighbor[4] = &getCell(cell.num, cell.num - (cellsNum_mu + 2) * (cellsNum_z + 2) * (cellsNum_nu - 1));
				break;

			case Type::RIGHT:
				neighbor[0] = &getCell(cell.num - cellsNum_z - 2);
				break;

			case Type::TOP:
				neighbor[0] = &getCell(cell.num + 1);
				break;

			case Type::BOTTOM:
				neighbor[0] = &getCell(cell.num - 1);
				break;

			case Type::WELL_LAT:
				neighbor[0] = &cells[nebrMap[cell.num].first];
				neighbor[1] = &cells[nebrMap[cell.num].second];
				break;

			case Type::WELL_TOP:
				neighbor[0] = &cells[nebrMap[cell.num].first];
				neighbor[1] = &cells[nebrMap[cell.num].second];
				break;

			case Type::WELL_BOT:
				neighbor[0] = &cells[nebrMap[cell.num].first];
				neighbor[1] = &cells[nebrMap[cell.num].second];
				break;
			}
		};
		inline double getPerm_mu(const Cell& cell) const
		{
			const double dz = fabs(cell.z - sk_well->h_well);
			const double kmu = ((dz <= cell.props->radius_eff_z && cell.mu <= cell.props->radius_eff_mu)
									? cell.props->perm_eff_mu : cell.props->perm_mu);
			return kmu;
			/*if (cos(cell.nu) > 0.0)
				return kmu;
			else
				return kmu / 2.0;*/
		};
		inline double getPerm_z(const Cell& cell) const
		{
			const double dz = fabs(cell.z - sk_well->h_well);
			const double kz = ((dz <= cell.props->radius_eff_z && cell.mu <= cell.props->radius_eff_mu) 
									? cell.props->perm_eff_z : cell.props->perm_z);
			return kz;
			/*if (cos(cell.nu) > 0.0)
				return kz;
			else
				return kz / 2.0;*/
		};
		inline double getDistance(const Cell& cell1, const Cell& cell2) const
		{
			double dist;
			if (fabs(cell1.z - cell2.z) > EQUALITY_TOLERANCE)
				dist = cell1.z - cell2.z;
			else if (fabs(cell1.mu - cell2.mu) > EQUALITY_TOLERANCE ||
				(cell1.mu == mu_init && cell2.mu == mu_init && fabs(cell1.nu - cell2.nu) > EQUALITY_TOLERANCE))
			{
				if (cell1.mu == mu_init && cell2.mu == mu_init)
					dist = Cell::getH(cell1.mu, cell1.nu) * (cell1.mu + cell2.mu);
				else
					dist = sign(cell1.mu - cell2.mu) * (Cell::getH(cell1.mu, cell1.nu) * cell1.hmu + Cell::getH(cell2.mu, cell2.nu) * cell2.hmu) / 2.0;
			}
			else if (fabs(cell1.nu - cell2.nu) > EQUALITY_TOLERANCE)
				dist = sign(cell1.nu - cell2.nu) * (Cell::getH(cell1.mu, cell1.nu) * cell1.hnu + Cell::getH(cell2.mu, cell2.nu) * cell2.hnu) / 2.0;

			return dist;
		};
		inline double getTrans(const Cell& cell, const Cell& beta) const
		{
			double k1, k2, S;
			const double dz1 = fabs(cell.z - sk_well->h_well);
			const double dz2 = fabs(beta.z - sk_well->h_well);

			if (fabs(cell.z - beta.z) > EQUALITY_TOLERANCE) {
				k1 = getPerm_z(cell);
				k2 = getPerm_z(beta);
				if (k1 == 0.0 && k2 == 0.0)
					return 0.0;

				S = Cell::getH(cell.mu, cell.nu) * Cell::getH(cell.mu, cell.nu) * cell.hmu * cell.hnu;
				return 2.0 * k1 * k2 * S / (k1 * beta.hz + k2 * cell.hz);
			}
			else if (fabs(cell.mu - beta.mu) > EQUALITY_TOLERANCE ||
				(cell.mu == mu_init && beta.mu == mu_init && fabs(cell.nu - beta.nu) > EQUALITY_TOLERANCE)) {

				k1 = getPerm_mu(cell);
				k2 = getPerm_mu(beta);
				if (k1 == 0.0 && k2 == 0.0)
					return 0.0;

				double mu;
				if (cell.mu == mu_init && beta.mu == mu_init)
					mu = 0.0;
				else
					mu = cell.mu + sign(beta.mu - cell.mu) * cell.hmu / 2.0;
				S = Cell::getH(mu, cell.nu) * cell.hnu * cell.hz;
				return 2.0 * k1 * k2 * S /
					(k1 * beta.hmu * Cell::getH(beta.mu, beta.nu)
						+ k2 * cell.hmu * Cell::getH(cell.mu, cell.nu));
			}
			else if (fabs(cell.nu - beta.nu) > EQUALITY_TOLERANCE) {
				k1 = getPerm_mu(cell);
				k2 = getPerm_mu(beta);
				if (k1 == 0.0 && k2 == 0.0)
					return 0.0;

				double nu;
				if (fabs(beta.nu - cell.nu) > 2.0 * cell.hnu + EQUALITY_TOLERANCE)
					nu = (cell.nu > beta.nu ? beta.nu - beta.hnu / 2 : cell.nu - cell.hnu / 2.0);
				else
					nu = cell.nu + sign(beta.nu - cell.nu) * cell.hnu / 2.0;
				S = cell.hz * cell.hmu * Cell::getH(cell.mu, nu);

				return 2.0 * k1 * k2 * S /
					(k1 * beta.hnu * Cell::getH(beta.mu, beta.nu) +
						+ k2 * cell.hnu * Cell::getH(cell.mu, cell.nu));
			}
		};
		inline adouble getDensity(adouble p1, const Cell& cell1, adouble p2, const Cell& cell2) const
		{
			double r1, r2;
			if (fabs(cell1.z - cell2.z) > EQUALITY_TOLERANCE) {
				r1 = cell1.hz;		r2 = cell2.hz;
			} else if (fabs(cell1.mu - cell2.mu) > EQUALITY_TOLERANCE ||
				(cell1.mu == mu_init && cell2.mu == mu_init && fabs(cell1.nu - cell2.nu) > EQUALITY_TOLERANCE)) {
				r1 = cell1.hmu * Cell::getH(cell1.mu, cell1.nu);
				r2 = cell2.hmu * Cell::getH(cell2.mu, cell2.nu);
			} else if (fabs(cell1.nu - cell2.nu) > EQUALITY_TOLERANCE) {
				r1 = cell1.hnu * Cell::getH(cell1.mu, cell1.nu);
				r2 = cell2.hnu * Cell::getH(cell2.mu, cell2.nu);
			}

			return (props_oil.getDensity(p1) * (adouble)r2 + props_oil.getDensity(p2) * (adouble)r1) / (adouble)(r1 + r2);
		};

		inline double getPerm(const Cell& cell, const int axis) const
		{
			if (axis == MU_AXIS)
				return getPerm_mu(cell);
			else if (axis == NU_AXIS)
				return getPerm_mu(cell);
			else if (axis == Z_AXIS)
				return getPerm_z(cell);
		};
		inline double getNablaP(Cell& cell, Cell** neighbor, const int axis)
		{
			Cell* nebr1;
			Cell* nebr2;
			double h, r_eff;

			if(axis == MU_AXIS)
			{
				if (cell.type == Type::MIDDLE_SIDE)
				{
					nebr1 = &cell;	nebr2 = neighbor[0];
				}
				else
				{
					nebr1 = neighbor[0];	nebr2 = neighbor[1];
					r_eff = cell.props->radius_eff_mu;
					if ((nebr1->mu < r_eff) && (nebr2->mu > r_eff))
					{
						if (cell.mu > r_eff)
							nebr1 = &cell;
						else
							nebr2 = &cell;
					}
				}

				if(nebr1->mu == mu_init && nebr2->mu == mu_init)
					h = Cell::getH(cell.mu, cell.nu) * (nebr2->mu + nebr1->mu);
				else
					h =	Cell::getH(cell.mu, cell.nu) * (nebr2->mu - nebr1->mu);
			}
			else if (axis == NU_AXIS)
			{
				if (cell.type == Type::MIDDLE_SIDE)
				{
					nebr1 = neighbor[3];	nebr2 = neighbor[4];
				}
				else
				{
					nebr1 = neighbor[4];	nebr2 = neighbor[5];
				}

				if (fabs(nebr2->nu - nebr1->nu) > 2.0 * cell.hnu + EQUALITY_TOLERANCE)
					h = Cell::getH(cell.mu, cell.nu) * (nebr2->nu - nebr1->nu + 2.0 * M_PI);
				else
					h = Cell::getH(cell.mu, cell.nu) * (nebr2->nu - nebr1->nu);
			}
			else if (axis == Z_AXIS)
			{
				if (cell.type == Type::MIDDLE_SIDE)
				{
					nebr1 = neighbor[1];	nebr2 = neighbor[2];
				}
				else
				{
					nebr1 = neighbor[2];	nebr2 = neighbor[3];

					r_eff = cell.props->radius_eff_z;
					const double dz_cur = fabs(cell.z - sk_well->h_well);
					const double dz1 = fabs(nebr1->z - sk_well->h_well);
					const double dz2 = fabs(nebr2->z - sk_well->h_well);
					if ((dz1 < r_eff) && (dz2 > r_eff))
					{
						if (dz_cur > r_eff)
							nebr1 = &cell;
						else
							nebr2 = &cell;
					}
				}

				h = nebr2->z - nebr1->z;
			}

			return (nebr2->u_next.p - nebr1->u_next.p) / h;
		};
		inline double getVelocity(Cell& cell, Cell** neighbor, const int axis)
		{
			return -getPerm(cell, axis) / props_oil.getViscosity(cell.u_next.p).value() * getNablaP(cell, neighbor, axis);
		};
		inline double getCn(const Cell& cell) const
		{
			return cell.props->getPoro(cell.u_next.p).value() * props_oil.getDensity(cell.u_next.p).value() * props_oil.c +
					(1.0 - cell.props->getPoro(cell.u_next.p).value()) * cell.props->getDensity(cell.u_next.p).value() * 
								cell.props->c;
		};
		inline double getAd(const Cell& cell) const
		{
			return cell.props->getPoro(cell.u_next.p).value() * (props_oil.getDensity(cell.u_next.p).value() * props_oil.ad * props_oil.c);
		};
		inline double getLambda(const Cell& cell, const int axis) const
		{
			if (axis == Z_AXIS)
				return cell.props->getPoro(cell.u_next.p).value() * props_oil.lambda +
						(1.0 - cell.props->getPoro(cell.u_next.p).value()) * cell.props->lambda_z;
			else
				return cell.props->getPoro(cell.u_next.p).value() * props_oil.lambda +
					(1.0 - cell.props->getPoro(cell.u_next.p).value()) * cell.props->lambda_r;
		};
		inline double getJT(Cell& cell, Cell** neighbor, const int axis)
		{
			return props_oil.getDensity(cell.u_next.p).value() *
				props_oil.c * props_oil.jt * getVelocity(cell, neighbor, axis);
		};
		inline double getA(Cell& cell, Cell** neighbor, const int axis)
		{
			return props_oil.getDensity(cell.u_next.p).value() * 
						props_oil.c * getVelocity(cell, neighbor, axis);
		};
		struct DivIndices 
		{
			double ther;
			double pres;
			DivIndices(double _ther, double _pres) : ther(_ther), pres(_pres) {};
		};
		inline DivIndices getDivCoeff(Cell& cell, Cell& beta, Cell** neighbor)
		{
			double r1, r2, lambda, A, a;
			if (fabs(cell.z - beta.z) > EQUALITY_TOLERANCE)
			{
				A = getA(cell, neighbor, Z_AXIS);
				a = sign(A) * std::max(0.0, sign(cell.z - beta.z) * A);
				r1 = cell.hz;		r2 = beta.hz;
				lambda = (r1 * getLambda(beta, Z_AXIS) + r2 * getLambda(cell, Z_AXIS)) / (r1 + r2) / r1;
			}
			else if (fabs(cell.mu - beta.mu) > EQUALITY_TOLERANCE ||
				(cell.mu == mu_init && beta.mu == mu_init && fabs(cell.nu - beta.nu) > EQUALITY_TOLERANCE))
			{
				double dmu;
				if (cell.mu == mu_init && beta.mu == mu_init)
					dmu = 1.0;
				else
					dmu = sign(cell.mu - beta.mu);
				A = getA(cell, neighbor, MU_AXIS);
				a = sign(A) * std::max(0.0, dmu * A);
				r1 = cell.hmu * Cell::getH(cell.mu, cell.nu);
				r2 = beta.hmu * Cell::getH(beta.mu, beta.nu);
				lambda = (r1 * getLambda(beta, MU_AXIS) + r2 * getLambda(cell, MU_AXIS)) / (r1 + r2) / r1;
			}
			else if (fabs(cell.nu - beta.nu) > EQUALITY_TOLERANCE)
			{
				A = getA(cell, neighbor, NU_AXIS);
				a = sign(A) * std::max(0.0, sign(cell.nu - beta.nu) * A);
				r1 = cell.hnu * Cell::getH(cell.mu, cell.nu);
				r2 = beta.hnu * Cell::getH(beta.mu, beta.nu);
				lambda = (r1 * getLambda(beta, NU_AXIS) + r2 * getLambda(cell, NU_AXIS)) / (r1 + r2) / r1;
			}

			DivIndices coeff (a + lambda, a * props_oil.jt);
			return coeff;
		};

		void solve_eqMiddle(const Cell& cell, const int val);
		void solve_eqWell(const Cell& cell, const int val);
		void solve_eqRight(const Cell& cell, const int val);
		void solve_eqVertical(const Cell& cell, const int val);
		void setVariables(const Cell& cell, const int val);

	public:
		OilNIT_Elliptic();
		~OilNIT_Elliptic();

		double* x;
		double* y;
		double** jac;

		void setPeriod(int period);
		double getRate(int cur) const;
	};
};

#endif /* OILNITELLIPTIC_HPP_ */