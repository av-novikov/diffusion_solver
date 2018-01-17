#ifndef WAXNIT_HPP_
#define WAXNIT_HPP_

#include "model/WaxNIT/Properties.hpp"
#include "model/cells/Variables.hpp"
#include "model/cells/CylCell2D.h"
#include "model/Basic2d/Basic2d.hpp"
#include "model/Basic2d/Basic2dSolver.hpp"
#include "util/Interpolate.h"

namespace wax_nit
{
	static const int R_AXIS = 0;
	static const int Z_AXIS = 1;

	static const int stencil = basic2d::stencil;
	static const int Lstencil = basic2d::Lstencil;
	static const int Rstencil = basic2d::Rstencil;
	static const int Vstencil = basic2d::Vstencil;

	typedef VarWaxNIT Variable;
	typedef TapeVarWaxNIT TapeVariable;
	typedef NewCylCell2D<Variable, Skeleton_Props> Cell;
	template <typename TVariable> using TCell = NewCylCell2D<TVariable, Skeleton_Props>;

	class WaxNIT : public basic2d::Basic2d<Variable, Properties, Skeleton_Props, TCell, WaxNIT>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class basic2d::Basic2dSolver;
		friend class WaxNITSolver;
	protected:
		Water_Props props_wat;
		Oil_Props props_oil;
		Gas_Props props_gas;

		void setProps(Properties& props);
		void makeDimLess();
		void setInitialState();

		// Service functions
		inline adouble getOilVelocity(const Cell& cell, const int axis) const
		{
			const TapeVariable& next = var[0];
			adouble tmp;
			adouble is_R_axis = (axis == R_AXIS) ? true : false;
			adouble is_Z_axis = (axis == Z_AXIS) ? true : false;

			condassign(tmp, is_R_axis, 
				-cell.props->getPermCoseni_r(next.m) * props_oil.getKr(next.s_w, next.s_o, cell.props) / props_oil.getViscosity(next.p) * 
				(var[4].p - var[3].p) / (cells[cell.num + cellsNum_z + 2].r - cells[cell.num - cellsNum_z - 2].r));
			condassign(tmp, is_Z_axis,
				-cell.props->getPermCoseni_z(next.m) * props_oil.getKr(next.s_w, next.s_o, cell.props) / props_oil.getViscosity(next.p) *
				(var[2].p - var[1].p) / (cells[cell.num + 1].z - cells[cell.num - 1].z));
			return tmp;
		};
		inline adouble getWatVelocity(const Cell& cell, const int axis) const
		{
			const TapeVariable& next = var[0];
			adouble tmp;
			adouble is_R_axis = (axis == R_AXIS) ? true : false;
			adouble is_Z_axis = (axis == Z_AXIS) ? true : false;

			condassign(tmp, is_R_axis,
				-cell.props->getPermCoseni_r(next.m) * props_wat.getKr(next.s_w, next.s_o, cell.props) / props_wat.getViscosity(next.p) *
				(var[4].p - var[3].p) / (cells[cell.num + cellsNum_z + 2].r - cells[cell.num - cellsNum_z - 2].r));
			condassign(tmp, is_Z_axis,
				-cell.props->getPermCoseni_z(next.m) * props_wat.getKr(next.s_w, next.s_o, cell.props) / props_wat.getViscosity(next.p) *
				(var[2].p - var[1].p) / (cells[cell.num + 1].z - cells[cell.num - 1].z));
		};
		inline adouble getGasVelocity(const Cell& cell, const int axis) const
		{
			const TapeVariable& next = var[0];
			adouble tmp;
			adouble is_R_axis = (axis == R_AXIS) ? true : false;
			adouble is_Z_axis = (axis == Z_AXIS) ? true : false;

			condassign(tmp, is_R_axis,
				-cell.props->getPermCoseni_r(next.m) * props_gas.getKr(next.s_w, next.s_o, cell.props) / props_gas.getViscosity(next.p) *
				(var[4].p - var[3].p) / (cells[cell.num + cellsNum_z + 2].r - cells[cell.num - cellsNum_z - 2].r));
			condassign(tmp, is_Z_axis,
				-cell.props->getPermCoseni_z(next.m) * props_gas.getKr(next.s_w, next.s_o, cell.props) / props_gas.getViscosity(next.p) *
				(var[2].p - var[1].p) / (cells[cell.num + 1].z - cells[cell.num - 1].z));
		};
		inline adouble getOilVelocityAbs(const Cell& cell) const 
		{
			adouble vel_r = getOilVelocity(cell, R_AXIS);
			adouble vel_z = getOilVelocity(cell, Z_AXIS);
			return sqrt(vel_r * vel_r + vel_z * vel_z);
		};
		inline adouble getTrans(const Cell& cell, adouble m_cell, const Cell& beta, adouble m_beta) const
		{
			adouble k1, k2, S;

			if (abs(cell.num - beta.num) == 1) {
				k1 = cell.props->getPermCoseni_z(m_cell);
				k2 = beta.props->getPermCoseni_z(m_beta);
				if (k1 == 0.0 && k2 == 0.0)
					return 0.0;
				S = 2.0 * M_PI * cell.r * cell.hr;
				return 2.0 * k1 * k2 * S / (k1 * beta.hz + k2 * cell.hz);
			}
			else {
				k1 = cell.props->getPermCoseni_r(m_cell);
				k2 = beta.props->getPermCoseni_r(m_beta);
				//k1 = (cell.r > cell.props->radius_eff ? cell.props->perm_r : cell.props->perm_eff);
				//k2 = (beta.r > beta.props->radius_eff ? beta.props->perm_r : beta.props->perm_eff);
				if (k1 == 0.0 && k2 == 0.0)
					return 0.0;
				S = 2.0 * M_PI * cell.hz * (cell.r + sign(beta.num - cell.num) * cell.hr / 2.0);
				return 2.0 * k1 * k2 * S / (k1 * beta.hr + k2 * cell.hr);
			}
		};
		inline adouble getAverage(adouble p1, const Cell& cell1, adouble p2, const Cell& cell2) const
		{
			double r1, r2;
			if (abs(cell1.num - cell2.num) == 1)
			{
				r1 = cell1.hz;		r2 = cell2.hz;
			}
			else
			{
				r1 = cell1.hr;		r2 = cell2.hr;
			}

			return (p1 * (adouble)r2 + p2 * (adouble)r1) / (adouble)(r1 + r2);
		};

		TapeVariable* var;
		adouble* h;
		inline adouble getCn(const Cell& cell) const
		{
			const TapeVariable& next = var[0];
			return next.m * (next.s_o * props_oil.getRho(next.p, next.p_bub, next.SATUR) * props_oil.c +
					next.s_w * props_wat.getRho(next.p) * props_wat.c +
					(1.0 - next.s_w - next.s_o) * props_gas.getRho(next.p) * props_gas.c) +
					(1.0 - next.m) * cell.props->dens_stc * cell.props->c;
		};
		inline adouble getAd(const Cell& cell) const
		{
			const TapeVariable& next = var[0];
			return next.m * (next.s_o * props_oil.getRho(next.p, next.p_bub, next.SATUR) * props_oil.ad * props_oil.c +
					next.s_w * props_wat.getRho(next.p) * props_wat.ad * props_wat.c +
					(1.0 - next.s_w - next.s_o) * props_gas.getRho(next.p) * props_gas.ad * props_gas.c);
		};
		inline adouble getLambda(const Cell& cell, const int axis) const
		{
			const TapeVariable& next = var[0];
			adouble tmp;
			adouble is_R_axis = (axis == R_AXIS) ? true : false;
			adouble is_Z_axis = (axis == Z_AXIS) ? true : false;
			condassign(tmp, is_R_axis, next.m *	(next.s_o * props_oil.lambda + next.s_w * props_wat.lambda + (1.0 - next.s_o - next.s_w) * props_gas.lambda) +
				(1.0 - next.m) * cell.props->lambda_r);
			condassign(tmp, is_Z_axis, next.m * (next.s_o * props_oil.lambda + next.s_w * props_wat.lambda + (1.0 - next.s_o - next.s_w) * props_gas.lambda) +
				(1.0 - next.m) * cell.props->lambda_z);
			return tmp;
		};
		inline adouble getJT(const Cell& cell, const int axis) const
		{
			const TapeVariable& next = var[0];
			return	props_oil.getRho(next.p, next.p_bub, next.SATUR) * props_oil.c * props_oil.jt * getOilVelocity(cell, axis) +
					props_wat.getRho(next.p) * props_wat.c * props_wat.jt * getWatVelocity(cell, axis) +
					props_gas.getRho(next.p) * props_gas.c * props_gas.jt * getGasVelocity(cell, axis);
		};
		inline adouble getA(const Cell& cell, const int axis) const
		{
			const TapeVariable& next = var[0];
			return	props_oil.getRho(next.p, next.p_bub, next.SATUR) * props_oil.c * getOilVelocity(cell, axis) +
					props_wat.getRho(next.p) * props_wat.c * getWatVelocity(cell, axis) +
					props_gas.getRho(next.p) * props_gas.c * getGasVelocity(cell, axis);
		};
		adouble phaseTrans(const Cell& cell);
		struct DivIndices
		{
			adouble ther;
			adouble pres;
			DivIndices(adouble _ther, adouble _pres) : ther(_ther), pres(_pres) {};
		};
		inline DivIndices getDivCoeff(const Cell& cell, const Cell& beta) const
		{
			adouble r1, r2, lambda, A, a, jt, JT;
			if (abs(cell.num - beta.num) == 1)
			{
				A = getA(cell, Z_AXIS);
				a = fmax((adouble)0.0, signA(cell.z - beta.z) * A);
				JT = getJT(cell, Z_AXIS);
				jt = fmax(0.0, sign(cell.z - beta.z) * JT);
				r1 = cell.hz;		r2 = beta.hz;
				lambda = (r1 * getLambda(beta, Z_AXIS) + r2 * getLambda(cell, Z_AXIS)) / (r1 + r2) / r1;
			} else
			{
				A = getA(cell, R_AXIS);
				a = fmax((adouble)0.0, signA(cell.r - beta.r) * A);
				JT = getJT(cell, Z_AXIS);
				jt = fmax(0.0, sign(cell.z - beta.z) * JT);
				r1 = cell.hr;		r2 = beta.hr;
				lambda = (r1 * getLambda(beta, Z_AXIS) + r2 * getLambda(cell, Z_AXIS)) / (r1 + r2) / r1;
			}
			DivIndices coeff(a + lambda, jt);
			return coeff;
		};

		void solve_eqMiddle(const Cell& cell);
		void solve_eqLeft(const Cell& cell);
		void solve_eqRight(const Cell& cell);
		void solve_eqVertical(const Cell& cell);
		void setVariables(const Cell& cell);

	public:
		WaxNIT();
		~WaxNIT();

		double getRate(int cur) const;

		double* x;
		double* y;
		double** jac;

		static const int var_size = Variable::size - 1;
	};
};


#endif /* WAXNIT_HPP_ */