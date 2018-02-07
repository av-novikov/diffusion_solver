#ifndef FRACTUREMODEL_HPP_
#define FRACTUREMODEL_HPP_

#include "model/AbstractModel.hpp"
#include "model/cells/Variables.hpp"
#include "model/cells/HexCell.hpp"

namespace frac
{
	static const int stencil = 7;

	struct Properties
	{
		int cellsNum_x, cellsNum_y, cellsNum_z;
		double w2, l2, height;

		double p_init, c_init;
	};

	typedef VarFrac TapeVariable;
	typedef TapeVarFrac Variable;
	//template <typename TVariable> using TCell = HexCell<TVariable, EmptyStruct>;
	typedef HexCell<Variable> Cell;
	typedef Cell::Type Type;

	class FractureModel : public AbstractModel<Variable, Properties, HexCell, FractureModel>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class snapshotter::VTKSnapshotter;
	protected:
		int cellsNum_x, cellsNum_y, cellsNum_z;
		double l2, w2, height;
		double p_init, c_init;

		void setProps(const Properties& props);
		void makeDimLess();
		void setInitialState();
		void buildGrid();

		TapeVariable* var;
		adouble* h;

	public:
		FractureModel();
		~FractureModel();
		void setPeriod(int period);

		double* x;
		double* y;
		double** jac;

		static const int var_size = Variable::size;
	};
}

#endif /* FRACTUREMODEL_HPP_  */