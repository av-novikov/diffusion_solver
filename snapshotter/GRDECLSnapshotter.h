#ifndef GRDECLSNAPSHOTTER_H_
#define GRDECLSNAPSHOTTER_H_

#include "snapshotter/Snapshotter.h"

template <class modelType>
class GRDECLSnapshotter : public Snapshotter<modelType>
{
protected:

public:
	GRDECLSnapshotter()
	{
		pattern = prefix + "snap_%{STEP}.grdecl";
	};
	~GRDECLSnapshotter()
	{

	};

	void dump(int i) {};
	void dump_all(int i) {};
};

#endif /* GRDECLSNAPSHOTTER_H_ */