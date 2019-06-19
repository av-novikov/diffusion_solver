#ifndef VTKSNAPSHOTTER_H_ 
#define VTKSNAPSHOTTER_H_

#include <vector>
#include <string>

#include "snapshotter/Snapshotter.h"

template <class modelType>
class VTKSnapshotter : public Snapshotter<modelType> {
protected:
    const std::string prefix;
public:
	VTKSnapshotter(const std::string _prefix = "snaps/");
	~VTKSnapshotter();

	void dump(int i);
	void dump_all(int i);
};

#endif /* VTKSNAPSHOTTER_H_ */