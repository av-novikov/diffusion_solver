#include "config.hpp"

using namespace issues;

int main(int ac, char* av[])
{ 
	using namespace acidrecfrac;
    int res = run<Properties,AcidRecFrac,AcidRecFracSolver>(ac, av);
	return res;
}