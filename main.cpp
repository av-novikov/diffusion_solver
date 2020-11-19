#include "config.hpp"

using namespace issues;

int main(int ac, char* av[])
{ 
	using namespace acidrecfrac;
	int res = run<Properties, AcidRecFrac, AcidRecFracSolver>(ac, av);
	//using namespace acid2drec;
	//int res = run<Properties, Acid2dRecModel, Acid2dRecSolver<ParSolver>>(ac, av);
	//using namespace acid2dnit;
	//int res = run<Properties, Acid2dNIT, Acid2dNITSolver>(ac, av);

	return res;
}