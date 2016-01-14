#include "tests/base-test.h"

#include "Scene.h"
//#include "model/Gas1D/Gas1D.h"
#include "tests/gas1D-test.h"
#include "model/Gas1D/Gas1D_simple.h"
#include "model/Gas1D/Gas1DSolver.h"

template <typename propsType, typename sceneType>
BaseTest<propsType, sceneType>::BaseTest()
{
}

template <typename propsType, typename sceneType>
BaseTest<propsType, sceneType>::~BaseTest()
{
}

template <typename propsType, typename sceneType>
void BaseTest<propsType, sceneType>::run()
{
	props = getProps();
	scene.load(*props);
	scene.setSnapshotterType("none");
	scene.start();
}

template class BaseTest<gas1D::Properties, Scene<gas1D::Gas1D_simple, gas1D::Gas1DSolSimp, gas1D::Properties> >;
template class BaseTest<gas1D::Properties, Scene<Gas1D_Wrapped, gas1D::Gas1DSol, gas1D::Properties> >;