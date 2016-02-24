#ifndef PARALUTIONINTERFACE_H_
#define PARALUTIONINTERFACE_H_

#include "paralution.hpp"

class ParSolver
{
public:
	paralution::LocalVector<double> x;
	paralution::LocalVector<double> Rhs;
	paralution::LocalMatrix<double> Mat;
	paralution::GMRES<paralution::LocalMatrix<double>, paralution::LocalVector<double>, double > ls;

	int* ind_i;
	int* ind_j;
	double* a;
	
	int* ind_rhs;
	double* rhs;

	virtual void Solve();
	virtual void copySolution() = 0;

	ParSolver();
	~ParSolver();
};

#endif /* PARALUTIONINTERFACE_H_ */