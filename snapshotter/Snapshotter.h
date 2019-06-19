#ifndef SNAPSHOTTER_H_
#define SNAPSHOTTER_H_

#include <string>
#include <vector>
#include <algorithm>
#include <fstream>

#include "model/cells/CylCell2D.h"

template <class modelType>
class Snapshotter {
protected:
	modelType* model;
	std::string pattern;
	std::string name;

	double r_dim;
	double t_dim;
	double P_dim;
	double T_dim;

	int nx, ny, nz;

	std::string replace(std::string filename, std::string from, std::string to);
	std::string getFileName(int i);
	std::string getFileName(int i, const std::string name);
public:
	Snapshotter();
	virtual ~Snapshotter();

	void setModel(modelType* _model);

	virtual void dump(int i) = 0;
	virtual void dump_all(int i) = 0;
};

#endif /* SNAPSHOTTER_H_ */