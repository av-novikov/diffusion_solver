#ifndef PARALUTIONINTERFACE_H_
#define PARALUTIONINTERFACE_H_

#include "paralution.hpp"

class ParSolver
{
protected:
	paralution::LocalVector<double> x;
	paralution::LocalVector<double> Rhs;
	paralution::LocalMatrix<double> Mat;
	paralution::GMRES<paralution::LocalMatrix<double>, paralution::LocalVector<double>, double > ls;

	bool isAssembled;
	int matSize;

public:
	void Init(int vecSize);
	void Assemble(const int* ind_i, const int* ind_j, const double* a, const int counter, const int* ind_rhs, const double* rhs);
	void AssembleUpdate(const double* a, const int* ind_rhs, const double* rhs);
	void Solve();

	const paralution::LocalVector<double>& getSolution();

	ParSolver();
	~ParSolver();
};

#endif /* PARALUTIONINTERFACE_H_ */