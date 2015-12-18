#ifndef SNAPSHOTTER_H_
#define SNAPSHOTTER_H_

#include <string>
#include <vector>
#include <algorithm>
#include <fstream>

#include "model/cells/CylindricalCell.h"

template <class modelType>
class Snapshotter {
protected:
	modelType* model;
	static const std::string prefix;
	std::string pattern;

	double r_dim;
	double t_dim;
	double T_dim;

	int nx, ny;

	std::string replace(std::string filename, std::string from, std::string to);
	std::string getFileName(int i);
public:
	Snapshotter();
	virtual ~Snapshotter();

	void setModel(modelType* _model);

	virtual void dump(int i) = 0;
	virtual void dump_all(int i) = 0;
};

#endif /* SNAPSHOTTER_H_ */