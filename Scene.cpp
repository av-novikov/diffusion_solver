#include "Scene.h"
#include "paralution.hpp"
#include <mpi.h>

#include "method/ParalutionInterface.h"
#include "method/HypreInterface.hpp"

#include "model/GasOil_RZ/GasOil_RZ.h"
#include "model/GasOil_RZ/GasOil2DSolver.h"

#include "model/VPP2d/VPP2d.hpp"
#include "model/VPP2d/VPPSolver.hpp"

#include "model/Bingham1d/Bingham1d.hpp"
#include "model/Bingham1d/BingSolver.hpp"

#include "model/GasOil_Elliptic/GasOil_Elliptic.hpp"
#include "model/GasOil_Elliptic/GasOilEllipticSolver.hpp"

#include "model/BlackOilNIT_Elliptic/BlackOilNIT_Elliptic.hpp"
#include "model/BlackOilNIT_Elliptic/BlackOilNITEllipticSolver.hpp"

#include "model/OilNIT_Elliptic/OilNIT_Elliptic.hpp"
#include "model/OilNIT_Elliptic/OilNITEllipticSolver.hpp"

#include "model/BlackOil_RZ/BlackOil_RZ.hpp"
#include "model/BlackOil_RZ/BlackOil2DSolver.hpp"

#include "model/Acid/2d/Acid2d.hpp"
#include "model/Acid/2d/Acid2dSolver.hpp"

#include "model/Acid/2dnit/Acid2dNIT.hpp"
#include "model/Acid/2dnit/Acid2dNITSolver.hpp"

#include "model/Acid/1d/Acid1d.hpp"
#include "model/Acid/1d/Acid1dSolver.hpp"

#include "model/Acid/frac/AcidFracModel.hpp"
#include "model/Acid/frac/AcidFracSolver.hpp"

#include "model/Acid/ellfrac/AcidEllFracModel.hpp"
#include "model/Acid/ellfrac/AcidEllFracSolver.hpp"

#include "model/Acid/recfrac/AcidRecFracModel.hpp"
#include "model/Acid/recfrac/AcidRecFracSolver.hpp"

#include "model/Acid/recfrac/RecFracProd.hpp"
#include "model/Acid/recfrac/RecFracProdSolver.hpp"

#include "model/Acid/recfracmov/AcidRecFracMovModel.hpp"
#include "model/Acid/recfracmov/AcidRecFracMovSolver.hpp"

#include "model/Acid/2drec/Acid2dRecModel.hpp"
#include "model/Acid/2drec/Acid2dRecSolver.hpp"

#include "model/WaxNIT/1d/WaxNIT1d.hpp"
#include "model/WaxNIT/1d/WaxNIT1dSolver.hpp"

#include "model/WaxNIT/2d/WaxNIT.hpp"
#include "model/WaxNIT/2d/WaxNITSolver.hpp"

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
Scene<acid2d::Acid2d, acid2d::Acid2dSolver, acid2d::Properties>::~Scene()
{
	paralution::stop_paralution();

	delete model;
	delete method;
}
Scene<acid2dnit::Acid2dNIT, acid2dnit::Acid2dNITSolver, acid2dnit::Properties>::~Scene()
{
	paralution::stop_paralution();

	delete model;
	delete method;
}
Scene<acidfrac::AcidFrac, acidfrac::AcidFracSolver, acidfrac::Properties>::~Scene()
{
	paralution::stop_paralution();

	delete model;
	delete method;
}
Scene<acidellfrac::AcidEllFrac, acidellfrac::AcidEllFracSolver, acidellfrac::Properties>::~Scene()
{
	paralution::stop_paralution();

	delete model;
	delete method;
}
Scene<acidrecfrac::AcidRecFrac, acidrecfrac::AcidRecFracSolver, acidrecfrac::Properties>::~Scene()
{
	paralution::stop_paralution();

	delete model;
	delete method;
}
Scene<acidrecfracmov::AcidRecFracMov, acidrecfracmov::AcidRecFracMovSolver, acidrecfracmov::Properties>::~Scene()
{
	paralution::stop_paralution();

	delete model;
	delete method;
}
Scene<wax_nit::WaxNIT, wax_nit::WaxNITSolver, wax_nit::Properties>::~Scene()
{
	paralution::stop_paralution();

	delete model;
	delete method;
}
Scene<wax_nit1d::WaxNIT1d, wax_nit1d::WaxNIT1dSolver, wax_nit1d::Properties>::~Scene()
{
	paralution::stop_paralution();

	delete model;
	delete method;
}
Scene<gasOil_elliptic::GasOil_Elliptic, gasOil_elliptic::GasOilEllipticSolver, gasOil_elliptic::Properties>::~Scene()
{
	paralution::stop_paralution();

	delete model;
	delete method;
}
Scene<oilnit_elliptic::OilNIT_Elliptic, oilnit_elliptic::OilNITEllipticSolver<ParSolver>, oilnit_elliptic::Properties>::~Scene()
{
	paralution::stop_paralution();

	delete model;
	delete method;
}
template <class modelType, class methodType, typename propsType>
void Scene<modelType, methodType, propsType>::load(propsType& props)
{
	model->load(props);
	method = new methodType(model);
}
void Scene<gasOil_elliptic::GasOil_Elliptic, gasOil_elliptic::GasOilEllipticSolver, gasOil_elliptic::Properties>::load(gasOil_elliptic::Properties& props)
{
	model->load(props);
	paralution::init_paralution();
	method = new gasOil_elliptic::GasOilEllipticSolver(model);
}
void Scene<acid2d::Acid2d, acid2d::Acid2dSolver, acid2d::Properties>::load(acid2d::Properties& props)
{
	model->load(props);
	paralution::init_paralution();
	method = new acid2d::Acid2dSolver(model);
}
void Scene<acid2dnit::Acid2dNIT, acid2dnit::Acid2dNITSolver, acid2dnit::Properties>::load(acid2dnit::Properties& props)
{
	model->load(props);
	paralution::init_paralution();
	method = new acid2dnit::Acid2dNITSolver(model);
}
void Scene<acid2dnit::Acid2dNIT, acid2dnit::Acid2dNITSolver, acid2dnit::Properties>::load(acid2dnit::Properties& props, int argc, char* argv[])
{
	model->load(props);
	model->setSnapshotter("VTK", model);
	paralution::init_paralution();
	method = new acid2dnit::Acid2dNITSolver(model);
}
void Scene<acidfrac::AcidFrac, acidfrac::AcidFracSolver, acidfrac::Properties>::load(acidfrac::Properties& props)
{
	model->load(props);
	paralution::init_paralution();
	method = new acidfrac::AcidFracSolver(model);
}
void Scene<acidellfrac::AcidEllFrac, acidellfrac::AcidEllFracSolver, acidellfrac::Properties>::load(acidellfrac::Properties& props)
{
	model->load(props);
	paralution::init_paralution();
	method = new acidellfrac::AcidEllFracSolver(model);
}
void Scene<acidrecfrac::AcidRecFrac, acidrecfrac::AcidRecFracSolver, acidrecfrac::Properties>::load(acidrecfrac::Properties& props)
{
	model->load(props);
	paralution::init_paralution();
	method = new acidrecfrac::AcidRecFracSolver(model);
}
void Scene<acidrecfracmov::AcidRecFracMov, acidrecfracmov::AcidRecFracMovSolver, acidrecfracmov::Properties>::load(acidrecfracmov::Properties& props)
{
	model->load(props);
	paralution::init_paralution();
	method = new acidrecfracmov::AcidRecFracMovSolver(model);
}
void Scene<wax_nit::WaxNIT, wax_nit::WaxNITSolver, wax_nit::Properties>::load(wax_nit::Properties& props)
{
	model->load(props);
	paralution::init_paralution();
	method = new wax_nit::WaxNITSolver(model);
}
void Scene<wax_nit1d::WaxNIT1d, wax_nit1d::WaxNIT1dSolver, wax_nit1d::Properties>::load(wax_nit1d::Properties& props)
{
	model->load(props);
	paralution::init_paralution();
	method = new wax_nit1d::WaxNIT1dSolver(model);
}
void Scene<oilnit_elliptic::OilNIT_Elliptic, oilnit_elliptic::OilNITEllipticSolver<ParSolver>, oilnit_elliptic::Properties>::load(oilnit_elliptic::Properties& props)
{
	model->load(props);
	paralution::init_paralution();
	method = new oilnit_elliptic::OilNITEllipticSolver<ParSolver>(model);
}
void Scene<blackoilnit_elliptic::BlackOilNIT_Elliptic, blackoilnit_elliptic::BlackOilNITEllipticSolver<ParSolver>, blackoilnit_elliptic::Properties>::load(blackoilnit_elliptic::Properties& props)
{
	model->load(props);
	paralution::init_paralution();
	method = new blackoilnit_elliptic::BlackOilNITEllipticSolver<ParSolver>(model);
}
void Scene<acid2drec::Acid2dRecModel, acid2drec::Acid2dRecSolver<ParSolver>, acid2drec::Properties>::load(acid2drec::Properties& props)
{
	model->load(props);
	paralution::init_paralution();
	method = new acid2drec::Acid2dRecSolver<ParSolver>(model);
}
void Scene<acid2drec::Acid2dRecModel, acid2drec::Acid2dRecSolver<ParSolver>, acid2drec::Properties>::load(acid2drec::Properties& props, int argc, char* argv[])
{
	model->load(props);
	paralution::init_paralution();
	method = new acid2drec::Acid2dRecSolver<ParSolver>(model);
}

Scene<acid2drec::Acid2dRecModel, acid2drec::Acid2dRecSolver<HypreSolver>, acid2drec::Properties>::~Scene()
{
	MPI_Finalize();

	delete model;
	delete method;
}
void Scene<acid2drec::Acid2dRecModel, acid2drec::Acid2dRecSolver<HypreSolver>, acid2drec::Properties>::load(acid2drec::Properties& props, int argc, char* argv[])
{
	model->load(props);
	MPI_Init(&argc, &argv);
	method = new acid2drec::Acid2dRecSolver<HypreSolver>(model);
}

Scene<acid2drec::Acid2dRecModel, acid2drec::Acid2dRecSolver<ParSolver>, acid2drec::Properties>::~Scene()
{
	paralution::stop_paralution();

	delete model;
	delete method;
}
#ifdef USE_HYPRE
Scene<oilnit_elliptic::OilNIT_Elliptic, oilnit_elliptic::OilNITEllipticSolver<HypreSolver>, oilnit_elliptic::Properties>::~Scene()
{
	MPI_Finalize();

	delete model;
	delete method;
}
Scene<blackoilnit_elliptic::BlackOilNIT_Elliptic, blackoilnit_elliptic::BlackOilNITEllipticSolver<HypreSolver>, blackoilnit_elliptic::Properties>::~Scene()
{
	MPI_Finalize();

	delete model;
	delete method;
}
void Scene<blackoilnit_elliptic::BlackOilNIT_Elliptic, blackoilnit_elliptic::BlackOilNITEllipticSolver<HypreSolver>, blackoilnit_elliptic::Properties>::load(blackoilnit_elliptic::Properties& props, int argc, char* argv[])
{
	model->load(props);
	MPI_Init(&argc, &argv);
	method = new blackoilnit_elliptic::BlackOilNITEllipticSolver<HypreSolver>(model);
}
void Scene<oilnit_elliptic::OilNIT_Elliptic, oilnit_elliptic::OilNITEllipticSolver<HypreSolver>, oilnit_elliptic::Properties>::load(oilnit_elliptic::Properties& props, int argc, char* argv[])
{
	model->load(props);
	MPI_Init(&argc, &argv);
	method = new oilnit_elliptic::OilNITEllipticSolver<HypreSolver>(model);
}
#endif /* USE_HYPRE */

template <class modelType, class methodType, typename propsType>
void Scene<modelType, methodType, propsType>::load(propsType& props, int i) {}
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
template class Scene<bing1d::Bingham1d, bing1d::Bing1dSolver, bing1d::Properties>;
template class Scene<gasOil_elliptic::GasOil_Elliptic, gasOil_elliptic::GasOilEllipticSolver, gasOil_elliptic::Properties>;
template class Scene<oilnit_elliptic::OilNIT_Elliptic, oilnit_elliptic::OilNITEllipticSolver<ParSolver>, oilnit_elliptic::Properties>;
template class Scene<blackoilnit_elliptic::BlackOilNIT_Elliptic, blackoilnit_elliptic::BlackOilNITEllipticSolver<ParSolver>, blackoilnit_elliptic::Properties>;
#ifdef USE_HYPRE
template class Scene<oilnit_elliptic::OilNIT_Elliptic, oilnit_elliptic::OilNITEllipticSolver<HypreSolver>, oilnit_elliptic::Properties>;
template class Scene<blackoilnit_elliptic::BlackOilNIT_Elliptic, blackoilnit_elliptic::BlackOilNITEllipticSolver<HypreSolver>, blackoilnit_elliptic::Properties>;
#endif /* USE_HYPRE */
template class Scene<blackoil_rz::BlackOil_RZ, blackoil_rz::BlackOil2dSolver, blackoil_rz::Properties>;
template class Scene<acid2d::Acid2d, acid2d::Acid2dSolver, acid2d::Properties>;
template class Scene<acid2dnit::Acid2dNIT, acid2dnit::Acid2dNITSolver, acid2dnit::Properties>;
template class Scene<acid1d::Acid1d, acid1d::Acid1dSolver, acid1d::Properties>;
template class Scene<acidfrac::AcidFrac, acidfrac::AcidFracSolver, acidfrac::Properties>;
template class Scene<acidellfrac::AcidEllFrac, acidellfrac::AcidEllFracSolver, acidellfrac::Properties>;
template class Scene<acidrecfrac::AcidRecFrac, acidrecfrac::AcidRecFracSolver, acidrecfrac::Properties>;
template class Scene<acidrecfracmov::AcidRecFracMov, acidrecfracmov::AcidRecFracMovSolver, acidrecfracmov::Properties>;
template class Scene<acid2drec::Acid2dRecModel, acid2drec::Acid2dRecSolver<ParSolver>, acid2drec::Properties>;
template class Scene<acid2drec::Acid2dRecModel, acid2drec::Acid2dRecSolver<HypreSolver>, acid2drec::Properties>;
template class Scene<wax_nit::WaxNIT, wax_nit::WaxNITSolver, wax_nit::Properties>;
template class Scene<wax_nit1d::WaxNIT1d, wax_nit1d::WaxNIT1dSolver, wax_nit1d::Properties>;