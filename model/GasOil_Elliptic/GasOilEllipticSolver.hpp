#ifndef GASOILELLIPTICSOLVER_HPP_
#define GASOILELLIPTICSOLVER_HPP_

#include <iostream>
#include <cstdlib>
#include <map>

#include "model/AbstractSolver.hpp"
#include "method/ParalutionInterface.h"
#include "model/GasOil_Elliptic/GasOil_Elliptic.hpp"

namespace gasOil_elliptic
{
	class GasOilEllipticSolver : public AbstractSolver<GasOil_Elliptic>
	{
	protected:
		void control() {};
		void doNextStep() {};
		void solveStep() {};
		void writeData() {};

	public:
		GasOilEllipticSolver(GasOil_Elliptic* _model);
		~GasOilEllipticSolver();
	};
};

#endif /* GASOILELLIPTICSOLVER_HPP_ */
