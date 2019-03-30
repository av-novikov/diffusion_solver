#include "config.hpp"
#include "Scene.h"

using namespace issues;

int main(int argc, char* argv[])
{ 
	using namespace acidrecfracmov;
	run<Properties, AcidRecFracMov, AcidRecFracMovSolver>();
	return 0;
}