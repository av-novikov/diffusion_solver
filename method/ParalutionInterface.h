#ifndef PARALUTIONINTERFACE_H_
#define PARALUTIONINTERFACE_H_

#include "paralution.hpp"

class ParSolver
{
protected:
	paralution::LocalVector<double> x;
	paralution::LocalVector<double> Rhs;
	paralution::LocalMatrix<double> Mat;

	int* ind_i;
	int* ind_j;
	double* a;
	
	int* ind_rhs;
	double* rhs;

public:
	ParSolver();
	~ParSolver();
};

#endif /* PARALUTIONINTERFACE_H_ */