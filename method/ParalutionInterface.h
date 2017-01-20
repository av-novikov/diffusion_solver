#ifndef PARALUTIONINTERFACE_H_
#define PARALUTIONINTERFACE_H_

#include <string>

#include "paralution.hpp"

class ParSolver
{
protected:
	paralution::LocalVector<double> x;
	paralution::LocalVector<double> Rhs;
	paralution::LocalMatrix<double> Mat;
	paralution::BiCGStab<paralution::LocalMatrix<double>, paralution::LocalVector<double>, double > bicgstab;
	void SolveBiCGStab();
	paralution::GMRES<paralution::LocalMatrix<double>, paralution::LocalVector<double>, double > gmres;
	void SolveGMRES();
	paralution::ILU<paralution::LocalMatrix<double>, paralution::LocalVector<double>, double> p;

	bool isAssembled;
	bool isPrecondBuilt;
	int matSize;

	inline void writeSystem()
	{
		Mat.WriteFileMTX("snaps/mat.mtx");
		Rhs.WriteFileASCII("snaps/rhs.dat");
		x.WriteFileASCII("snaps/x.dat");
	};

	double initRes, finalRes;
	int iterNum;
	const std::string resHistoryFile;
	void getResiduals();

public:
	void Init(int vecSize);
	void Assemble(const int* ind_i, const int* ind_j, const double* a, const int counter, const int* ind_rhs, const double* rhs);
	void Solve();

	const paralution::LocalVector<double>& getSolution();

	ParSolver();
	~ParSolver();
};

#endif /* PARALUTIONINTERFACE_H_ */