#ifndef UTILS_H_
#define UTILS_H_

#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include "util/Interpolate.h"

#define BAR_TO_PA 1.E5
#define P_ATM 1.01325
#define EQUALITY_TOLERANCE 1.E-6

using std::vector;
using std::pair;
using std::string;
using std::sort;

inline bool IsNan(double a)
{
	if (a!=a)  return true;
	return false;
};

inline double sign(int a)
{
	if (a > 0) return 1.0;
	else if (a < 0) return -1.0;
	else return 0.0;
};

inline double MilliDarcyToM2(double perm)
{
	return perm * 0.986923 * 1.E-15;
};

inline double M2toMilliDarcy(double perm)
{
	return perm * 1.E15 / 0.986923;
};

inline double cPToPaSec(double visc)
{
 return visc/1000.0;
};

template <typename T>
inline std::string to_string(T value)
{
	std::ostringstream os ;
	os << value;
	return os.str();
};

struct sort_pair_first {
    bool operator() (const std::pair<double,double> &left, const std::pair<double,double> &right) 
	{
        return left.first < right.first;
    }
};

struct sort_pair_second {
    bool operator() (const std::pair<double,double> &left, const std::pair<double,double> &right) 
	{
        return left.second < right.second;
    }
};

inline Interpolate* setDataset(vector< pair<double,double> >& vec, const double xDim, const double yDim)
{
	sort(vec.begin(), vec.end(), sort_pair_first());

	const int N = vec.size();
	double* x = new double [N];
	double* y  = new double [N];

	for (int i = 0; i < N; i++)
	{
		x[i] = vec[i].first / xDim;
		y[i] = vec[i].second / yDim;
	}

	return new Interpolate(x, y, N);
};

inline Interpolate* setInvDataset(vector< pair<double,double> >& vec, const double xDim, const double yDim)
{
	sort(vec.begin(), vec.end(), sort_pair_second());

	const int N = vec.size();
	double* x = new double [N];
	double* y  = new double [N];

	for (int i = 0; i < N; i++)
	{
		x[i] = vec[i].second / xDim;
		y[i] = vec[i].first / yDim;
	}

	return new Interpolate(x, y, N);
};

#endif /* UTILS_H_ */