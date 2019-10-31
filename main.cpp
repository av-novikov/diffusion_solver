#include "config.hpp"

using namespace issues;

int main(int ac, char* av[])
{ 
	using namespace acid2drec;
	int res = run<Properties, Acid2dRecModel, Acid2dRecSolver>(ac, av);
	return res;
}