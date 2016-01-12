#include "Scene.h"

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

#include "tests/gas1D-test.h"

using namespace std;

template <class modelType, class methodType, typename propsType>
Scene<modelType, methodType, propsType>::Scene()
{
	model = new modelType();
}

template <class modelType, class methodType, typename propsType>
Scene<modelType, methodType, propsType>::~Scene()
{
	delete model;
	delete method;
}

template <class modelType, class methodType, typename propsType>
void Scene<modelType, methodType, propsType>::load(propsType& props)
{
	model->load(props);
	method = new methodType(model);
}

template <class modelType, class methodType, typename propsType>
void Scene<modelType, methodType, propsType>::load(propsType& props, int i)
{
}

template <>
void Scene<gas1D::Gas1D, gas1D::Gas1DSol, gas1D::Properties>::load(gas1D::Properties& props, int i)
{
	model->load(props);
	method = new gas1D::Gas1DSolver<gas1D::Gas1D>(model, i);
}

template <class modelType, class methodType, typename propsType>
void Scene<modelType, methodType, propsType>::setSnapshotterType(std::string type)
{
	model->setSnapshotter(type, model);
}

template <class modelType, class methodType, typename propsType>
void Scene<modelType, methodType, propsType>::start()
{
	method->start();
}

template <class modelType, class methodType, typename propsType>
modelType* Scene<modelType, methodType, propsType>::getModel() const
{
	return model;
}

template class Scene<oil1D::Oil1D, oil1D::Oil1DSolver, oil1D::Properties>;
template class Scene<gas1D::Gas1D, gas1D::Gas1DSol, gas1D::Properties>;
template class Scene<gas1D::Gas1D_simple, gas1D::Gas1DSolSimp, gas1D::Properties>;
template class Scene<oil1D_NIT::Oil1D_NIT, oil1D_NIT::Oil1DNITSolver, oil1D_NIT::Properties>;
template class Scene<oil_rz::Oil_RZ, oil_rz::OilRZSolver, oil_rz::Properties>;
template class Scene<gasOil_rz::GasOil_RZ, gasOil_rz::GasOil2DSolver, gasOil_rz::Properties>;
template class Scene<gasOil_rz_NIT::GasOil_RZ_NIT, gasOil_rz_NIT::GasOil2DNITSolver, gasOil_rz_NIT::Properties>;

template class Scene<Gas1D_Wrapped, gas1D::Gas1DSol, gas1D::Properties>;