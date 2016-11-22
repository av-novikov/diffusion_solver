#include "snapshotter/GRDECLSnapshotter.h"
#include "model/AbstractModel.hpp"

#include "model/GasOil_RZ/GasOil_RZ.h"

#include "model/Acid/2d/Acid2d.hpp"
#include "model/VPP2d/VPP2d.hpp"

using namespace std;

template <class modelType>
GRDECLSnapshotter<modelType>::GRDECLSnapshotter()
{
	pattern = prefix + "snap_%{STEP}.grdecl";
}

template <class modelType>
GRDECLSnapshotter<modelType>::~GRDECLSnapshotter()
{
}

template <class modelType>
void GRDECLSnapshotter<modelType>::dump(int i)
{
}

template <class modelType>
void GRDECLSnapshotter<modelType>::dump_all(int i)
{
}

template class GRDECLSnapshotter<gasOil_rz::GasOil_RZ>;

template class GRDECLSnapshotter<acid2d::Acid2d>;
template class GRDECLSnapshotter<vpp2d::VPP2d>;