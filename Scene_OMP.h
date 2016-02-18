#ifndef SCENE_OMP_H_
#define SCENE_OMP_H_

#include <new>
#include <string>

#include "model/Oil1D/Oil1D.h"
#include "model/Oil1D/Oil1DSolver.h"

#include "model/Gas1D/Gas1D.h"
#include "model/Gas1D/Gas1D_simple.h"
#include "model/Gas1D/Gas1DSolver.h"

#include "model/Oil1D_NIT/Oil1D_NIT.h"
#include "model/Oil1D_NIT/Oil1DNITSolver.h"

#include "model/Oil_RZ/Oil_RZ.h"
#include "model/Oil_RZ/OilRZSolver.h"

#include "model/GasOil_RZ/GasOil_RZ.h"
#include "model/GasOil_RZ/GasOil2DSolver.h"

#include "model/GasOil_RZ_NIT/GasOil_RZ_NIT.h"
#include "model/GasOil_RZ_NIT/GasOil2DNITSolver.h"

#include "model/3D/GasOil_3D/GasOil_3D.h"
#include "model/3D/GasOil_3D/GasOil3DSolver.h"

#include "model/3D/GasOil_3D_NIT/GasOil_3D_NIT.h"
#include "model/3D/GasOil_3D_NIT/GasOil3DNITSolver.h"

#include "tests/gas1D-test.h"

template <int n, class modelType, class methodType, typename propsType>
class Scene_OMP
{
protected:
	modelType* model [n];
	methodType* method [n];

public:
	Scene_OMP()
	{
		for(int i = 0; i < n; i++)
			model[i] = new modelType();
	};

	~Scene_OMP()
	{
		delete [] model;
		delete [] method;
	};
	
	void load(propsType& props)
	{
		for(int j = 0; j < n; j++)
		{
			model[j]->load(props);
			method[j] = new methodType(model[j]);
		}
	};

	void setSnapshotterType(std::string type)
	{
		for(int i = 0; i < n; i++)
		{
			model[i]->setSnapshotter(type, model[i]);
		}
	}

	void start()
	{
		method[0]->start();
	};

	modelType* getModel(int i) const
	{
		return model[i];
	};
};

#endif /* SCENE_OMP_H_ */