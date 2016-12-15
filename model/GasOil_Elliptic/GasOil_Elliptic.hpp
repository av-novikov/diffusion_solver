#ifndef GASOILELLIPTIC_HPP_
#define GASOILELLIPTIC_HPP_

#include <vector>
#include <map>
#include <string>

#include "model/cells/Variables.hpp"
#include "model/GasOil_Elliptic/Properties.hpp"
#include "model/cells/EllipticCell.hpp"
#include "model/AbstractModel.hpp"

namespace gasOil_elliptic
{
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
		friend class GasOilEllipticSolver;
	protected:

		// Continuum properties
		int skeletonsNum;
		std::vector<Skeleton_Props> props_sk;
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

	public:
		GasOil_Elliptic();
		~GasOil_Elliptic();

		void setPeriod(int period);
	};
};


#endif /* GASOILELLIPTIC_HPP_ */
