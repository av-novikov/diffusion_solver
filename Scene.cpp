#include "Scene.h"

#include "model/GasOil_RZ/GasOil_RZ.h"
#include "model/GasOil_RZ/GasOil2DSolver.h"

#include "model/VPP2d/VPP2d.hpp"
#include "model/VPP2d/VPPSolver.hpp"

#include "paralution.hpp"

using namespace std;
using namespace paralution;

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

template class Scene<gasOil_rz::GasOil_RZ, gasOil_rz::GasOil2DSolver, gasOil_rz::Properties>;
template class Scene<vpp2d::VPP2d, vpp2d::VPPSolver, vpp2d::Properties>;