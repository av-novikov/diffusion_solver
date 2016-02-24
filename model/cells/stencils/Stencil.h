#ifndef STENCIL_H_
#define STENCIL_H_

#include "model\AbstractModel.hpp"
#include "util\utils.h"

#include <vector>
#include <functional>
#include <new>
#include <initializer_list>

/*--------------------MidStencil--------------------*/

template <class modelType>
class MidStencil
{
protected:
	modelType* model;
	FillFoo foo;

	double* matrix;
	int* ind_i;
	int* ind_j;
	double* rhs;

public:
	MidStencil(modelType* _model, FillFoo _foo);
	~MidStencil();

	void setStorage(double* _matrix, int* _ind_i, int* _ind_j, double* _rhs);

	void fill(int cellIdx, int* counter);
	void fillIndex(int cellIdx, int *counter);
};

/*--------------------LeftStencil--------------------*/

template <class modelType>
class LeftStencil
{
protected:
	modelType* model;
	FillFoo foo;

	double* matrix;
	int* ind_i;
	int* ind_j;
	double* rhs;

public:
	LeftStencil(modelType* _model, FillFoo _foo);
	~LeftStencil();

	void setStorage(double* _matrix, int* _ind_i, int* _ind_j, double* _rhs);

	void fill(int cellIdx, int *counter);
	void fillIndex(int cellIdx, int *counter);
};

/*--------------------RightStencil--------------------*/

template <class modelType>
class RightStencil
{
protected:
	modelType* model;
	FillFoo foo;

	double* matrix;
	int* ind_i;
	int* ind_j;
	double* rhs;

public:
	RightStencil(modelType* _model, FillFoo _foo);
	~RightStencil();

	void setStorage(double* _matrix, int* _ind_i, int* _ind_j, double* _rhs);

	void fill(int cellIdx, int *counter);
	void fillIndex(int cellIdx, int *counter);
};

/*--------------------TopStencil--------------------*/

template <class modelType>
class TopStencil
{
protected:
	modelType* model;
	FillFoo foo;

	double* matrix;
	int* ind_i;
	int* ind_j;
	double* rhs;

public:
	TopStencil(modelType* _model, FillFoo _foo);
	~TopStencil();

	void setStorage(double* _matrix, int* _ind_i, int* _ind_j, double* _rhs);

	void fill(int cellIdx, int *counter);
	void fillIndex(int cellIdx, int *counter);
};

/*--------------------BotStencil--------------------*/

template <class modelType>
class BotStencil
{
protected:
	modelType* model;
	FillFoo foo;

	double* matrix;
	int* ind_i;
	int* ind_j;
	double* rhs;

public:
	BotStencil(modelType* _model, FillFoo _foo);
	~BotStencil();

	void setStorage(double* _matrix, int* _ind_i, int* _ind_j, double* _rhs);

	void fill(int cellIdx, int *counter);
	void fillIndex(int cellIdx, int *counter);
};

/*--------------------UsedStencils--------------------*/

template <class modelType>
class UsedStencils
{
public:
	UsedStencils(modelType* _model);
	~UsedStencils();

	void setStorages(double* _matrix, int* _ind_i, int* _ind_j, double* _rhs);

	MidStencil<modelType>* middle;
	LeftStencil<modelType>* left;
	RightStencil<modelType>* right;
	TopStencil<modelType>* top;
	BotStencil<modelType>* bot;
};

#endif /* STENCIL_H_ */
