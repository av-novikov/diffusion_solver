#ifndef ACIDFRACSOLVER_HPP_
#define ACIDFRACSOLVER_HPP_

#include "model/Acid/frac/AcidFracModel.hpp"

namespace acidfrac
{
	class AcidFracSolver
	{
	protected:
		AcidFrac* model;
	public:
		AcidFracSolver(AcidFrac* _model) : model(_model) {}
		~AcidFracSolver() {};

		void start() {};
	};
}

#endif /* ACIDFRACSOLVER_HPP_ */
