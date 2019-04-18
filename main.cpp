#include "config.hpp"

using namespace issues;

int main(int argc, char* argv[])
{ 
	using namespace acidrecfrac;
    run<Properties,AcidRecFrac,AcidRecFracSolver>();
	return 0;
}