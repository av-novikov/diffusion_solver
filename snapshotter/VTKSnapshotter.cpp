#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolygon.h>

#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkHexahedron.h>

#include "snapshotter/VTKSnapshotter.h"

#include "model/GasOil_RZ/GasOil_RZ.h"
#include "model/Acid/2dnit/Acid2dNIT.hpp"
#include "model/Acid/2d/Acid2d.hpp"
#include "model/Acid/1d/Acid1d.hpp"
#include "model/Acid/frac/AcidFracModel.hpp"
#include "model/Acid/ellfrac/AcidEllFracModel.hpp"
#include "model/Acid/recfrac/AcidRecFracModel.hpp"
#include "model/Acid/recfrac/RecFracProd.hpp"
#include "model/Acid/recfracmov/AcidRecFracMovModel.hpp"
#include "model/VPP2d/VPP2d.hpp"
#include "model/Bingham1d/Bingham1d.hpp"
#include "model/GasOil_Elliptic/GasOil_Elliptic.hpp"
#include "model/BlackOilNIT_Elliptic/BlackOilNIT_Elliptic.hpp"
#include "model/OilNIT_Elliptic/OilNIT_Elliptic.hpp"
#include "model/BlackOil_RZ/BlackOil_RZ.hpp"
#include "model/WaxNIT/2d/WaxNIT.hpp"
#include "model/WaxNIT/1d/WaxNIT1d.hpp"

#include <cmath>

#define VIEW_MULTIPLIER 1
#define FRAC_WIDTH_MULT 20.0

using namespace std;

template<>
VTKSnapshotter<acidfrac::AcidFrac>::VTKSnapshotter(const std::string _prefix) : prefix(_prefix)
{
	pattern = prefix + "AcidFrac_%{NAME}_%{STEP}.vtu";
}
template<>
VTKSnapshotter<acidellfrac::AcidEllFrac>::VTKSnapshotter(const std::string _prefix) : prefix(_prefix)
{
	pattern = prefix + "AcidEllFrac_%{NAME}_%{STEP}.vtu";
}
template<>
VTKSnapshotter<acidrecfrac::AcidRecFrac>::VTKSnapshotter(const std::string _prefix) : prefix(_prefix)
{
	pattern = prefix + "AcidRecFrac_%{NAME}_%{STEP}.vtu";
}
template<>
VTKSnapshotter<acidrecfrac_prod::RecFracProd>::VTKSnapshotter(const std::string _prefix) : prefix(_prefix)
{
    pattern = prefix + "AcidRecFrac_prod_%{STEP}.vtu";
}
template<>
VTKSnapshotter<acidrecfracmov::AcidRecFracMov>::VTKSnapshotter(const std::string _prefix) : prefix(_prefix)
{
	pattern = prefix + "AcidRecFracMov_%{NAME}_%{STEP}.vtu";
}
template <class modelType>
VTKSnapshotter<modelType>::VTKSnapshotter(const std::string _prefix) : prefix(_prefix)
{
	pattern = prefix + "snap_%{STEP}.vtp";
}
VTKSnapshotter<gasOil_rz::GasOil_RZ>::VTKSnapshotter(const std::string _prefix) : prefix(_prefix)
{
	pattern = prefix + "GasOil_RZ_%{STEP}.vtp";
}
VTKSnapshotter<acid2dnit::Acid2dNIT>::VTKSnapshotter(const std::string _prefix) : prefix(_prefix)
{
	pattern = prefix + "Acid2dNIT_%{STEP}.vtp";
}
VTKSnapshotter<acid2d::Acid2d>::VTKSnapshotter(const std::string _prefix) : prefix(_prefix)
{
	pattern = prefix + "Acid2d_%{STEP}.vtp";
}
VTKSnapshotter<acid1d::Acid1d>::VTKSnapshotter(const std::string _prefix) : prefix(_prefix)
{
	pattern = prefix + "Acid1d_%{STEP}.vtp";
}
VTKSnapshotter<vpp2d::VPP2d>::VTKSnapshotter(const std::string _prefix) : prefix(_prefix)
{
	pattern = prefix + "VPP2d_%{STEP}.vtp";
}
VTKSnapshotter<bing1d::Bingham1d>::VTKSnapshotter(const std::string _prefix) : prefix(_prefix)
{
	pattern = prefix + "Bing1d_%{STEP}.vtp";
}
VTKSnapshotter<gasOil_elliptic::GasOil_Elliptic>::VTKSnapshotter(const std::string _prefix) : prefix(_prefix)
{
	pattern = prefix + "GasOil_El_%{STEP}.vtu";
}
VTKSnapshotter<oilnit_elliptic::OilNIT_Elliptic>::VTKSnapshotter(const std::string _prefix) : prefix(_prefix)
{
	pattern = prefix + "OilNIT_El_%{STEP}.vtu";
}
VTKSnapshotter<blackoilnit_elliptic::BlackOilNIT_Elliptic>::VTKSnapshotter(const std::string _prefix) : prefix(_prefix)
{
	pattern = prefix + "BlackOilNIT_El_%{STEP}.vtu";
}
VTKSnapshotter<blackoil_rz::BlackOil_RZ>::VTKSnapshotter(const std::string _prefix) : prefix(_prefix)
{
	pattern = prefix + "BlackOil_RZ_%{STEP}.vtp";
}
VTKSnapshotter<wax_nit::WaxNIT>::VTKSnapshotter(const std::string _prefix) : prefix(_prefix)
{
	pattern = prefix + "WaxNIT_%{STEP}.vtp";
}
VTKSnapshotter<wax_nit1d::WaxNIT1d>::VTKSnapshotter(const std::string _prefix) : prefix(_prefix)
{
	pattern = prefix + "WaxNIT1d_%{STEP}.vtp";
}
template <class modelType>
VTKSnapshotter<modelType>::~VTKSnapshotter()
{
}
template <class modelType>
void VTKSnapshotter<modelType>::dump(int i)
{
}
template <class modelType>
void VTKSnapshotter<modelType>::dump_all(int i)
{
}

void VTKSnapshotter<gasOil_rz::GasOil_RZ>::dump_all(int i)
{
	using namespace gasOil_rz;

	auto grid =	vtkSmartPointer<vtkPolyData>::New();
	auto points = vtkSmartPointer<vtkPoints>::New();
	for(int j = 1; j < ny; j++)
	{
		Cell& cell = model->cells[ j ];
		points->InsertNextPoint(r_dim * (0.9 * cell.r), -r_dim * (cell.z-cell.hz/2.0), 0.0);
	}
	for(int k = 1; k < nx; k++)
	{
		for(int j = 1; j < ny; j++)
		{
			Cell& cell = model->cells[ k * ny + j ];
			points->InsertNextPoint(r_dim * (cell.r-cell.hr/2.0), -r_dim * (cell.z-cell.hz/2.0), 0.0);
		}
	}
	grid->SetPoints(points);

	// Data
	auto polygons =	vtkSmartPointer<vtkCellArray>::New();
	auto polygon = vtkSmartPointer<vtkPolygon>::New();
	polygon->GetPointIds()->SetNumberOfIds(4);

	auto pres = vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");
	auto p_bub = vtkSmartPointer<vtkDoubleArray>::New();
	p_bub->SetName("buble_point");
	auto sat_oil = vtkSmartPointer<vtkDoubleArray>::New();
	sat_oil->SetName("oilSaturation");
	auto sat_gas = vtkSmartPointer<vtkDoubleArray>::New();
	sat_gas->SetName("gasSaturation");
	auto satur = vtkSmartPointer<vtkIntArray>::New();
	satur->SetName("SATUR");
	auto vel_oil = vtkSmartPointer<vtkDoubleArray>::New();
	vel_oil->SetName("oilVelocity");
	vel_oil->SetNumberOfComponents(3);
	auto vel_gas = vtkSmartPointer<vtkDoubleArray>::New();
	vel_gas->SetName("gasVelocity");
	vel_gas->SetNumberOfComponents(3);

	int k, j, idx, idx1;
	double vel [3];

	for(j = 0; j < ny-2; j++)
	{
		idx = j;
		idx1 = j + 1;
		Cell& cell = model->cells[idx1];

		polygon->GetPointIds()->SetId(0, idx);
		polygon->GetPointIds()->SetId(1, idx+ny-1);
		polygon->GetPointIds()->SetId(2, idx+ny);
		polygon->GetPointIds()->SetId(3, idx+1);
		polygons->InsertNextCell(polygon);

		pres->InsertNextValue(cell.u_next.p * P_dim / BAR_TO_PA);
		p_bub->InsertNextValue(cell.u_next.p_bub * P_dim / BAR_TO_PA);
		satur->InsertNextValue(cell.u_next.SATUR);
		sat_oil->InsertNextValue(cell.u_next.s);
		sat_gas->InsertNextValue(1.0 - cell.u_next.s);
		vel[0] = 0.0;
		vel[1] = 0.0;	
		vel[2] = 0.0;
		vel_oil->InsertNextTuple(vel);
		vel[0] = 0.0;
		vel[1] = 0.0;	
		vel[2] = 0.0;	
		vel_gas->InsertNextTuple(vel);
		//vel[0] = model->getOilVelocity(cell, NEXT, R_AXIS);	vel[1] = model->getOilVelocity(cell, NEXT, Z_AXIS);	
		//vel_oil->InsertNextValue(vel);
	}

	// Middle cells
	for(k = 1; k < nx-1; k++)
	{
		for(j = 0; j < ny-2; j++)
		{
			idx = k * (ny-1) + j;
			idx1 = k * ny + j + 1;
			Cell& cell = model->cells[idx1];

			polygon->GetPointIds()->SetId(0, idx);
			polygon->GetPointIds()->SetId(1, idx+ny-1);
			polygon->GetPointIds()->SetId(2, idx+ny);
			polygon->GetPointIds()->SetId(3, idx+1);
			polygons->InsertNextCell(polygon);

			pres->InsertNextValue(cell.u_next.p * P_dim / BAR_TO_PA);
			p_bub->InsertNextValue(cell.u_next.p_bub * P_dim / BAR_TO_PA);
			satur->InsertNextValue(cell.u_next.SATUR);
			sat_oil->InsertNextValue(cell.u_next.s);
			sat_gas->InsertNextValue(1.0 - cell.u_next.s);
			vel[0] = r_dim / t_dim * model->getOilVelocity(cell, NEXT, R_AXIS);
			vel[1] = r_dim / t_dim * model->getOilVelocity(cell, NEXT, Z_AXIS);	
			vel[2] = 0.0;
			vel_oil->InsertNextTuple(vel);
			vel[0] = r_dim / t_dim * model->getGasVelocity(cell, NEXT, R_AXIS);
			vel[1] = r_dim / t_dim * model->getGasVelocity(cell, NEXT, Z_AXIS);	
			vel[2] = 0.0;	
			vel_gas->InsertNextTuple(vel);
			//vel[0] = model->getOilVelocity(cell, NEXT, R_AXIS);	vel[1] = model->getOilVelocity(cell, NEXT, Z_AXIS);	
			//vel_oil->InsertNextValue(vel);
		}
	}

	grid->SetPolys(polygons);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(pres);
	fd->AddArray(p_bub);
	fd->AddArray(satur);
	fd->AddArray(sat_oil);
	fd->AddArray(sat_gas);
	fd->AddArray(vel_oil);
	fd->AddArray(vel_gas);

	// Writing
	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	
	writer->SetFileName( getFileName(i).c_str() );
	writer->SetInputData(grid);
	writer->Write();
}
void VTKSnapshotter<acid2d::Acid2d>::dump_all(int i)
{
	using namespace acid2d;

	// Grid
	auto grid = vtkSmartPointer<vtkPolyData>::New();

	// Points
	auto points = vtkSmartPointer<vtkPoints>::New();
	for (int j = 1; j < ny; j++)
	{
		Cell& cell = model->cells[j];
		points->InsertNextPoint(r_dim * (0.99 * cell.r), -r_dim * (cell.z - cell.hz / 2.0), 0.0);
	}
	for (int k = 1; k < nx; k++)
	{
		for (int j = 1; j < ny; j++)
		{
			Cell& cell = model->cells[k * ny + j];
			points->InsertNextPoint(r_dim * (cell.r - cell.hr / 2.0), -r_dim * (cell.z - cell.hz / 2.0), 0.0);
		}
	}
	grid->SetPoints(points);

	// Data
	auto polygons =	vtkSmartPointer<vtkCellArray>::New();
	auto polygon = vtkSmartPointer<vtkPolygon>::New();
	polygon->GetPointIds()->SetNumberOfIds(4);
	auto poro = vtkSmartPointer<vtkDoubleArray>::New();
	poro->SetName("porosity");
	auto perm = vtkSmartPointer<vtkDoubleArray>::New();
	perm->SetName("permeability");
	auto pres =	vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");
	auto p_bub = vtkSmartPointer<vtkDoubleArray>::New();
	p_bub->SetName("bublePoint");
	auto satur = vtkSmartPointer<vtkIntArray>::New();
	satur->SetName("SATUR");
	auto sat_w = vtkSmartPointer<vtkDoubleArray>::New();
	sat_w->SetName("WaterSaturation");
	auto sat_o = vtkSmartPointer<vtkDoubleArray>::New();
	sat_o->SetName("OilSaturation");
	auto sat_g = vtkSmartPointer<vtkDoubleArray>::New();
	sat_g->SetName("GasSaturation");
	auto conc_a = vtkSmartPointer<vtkDoubleArray>::New();
	conc_a->SetName("AcidConcentration");
	auto conc_w = vtkSmartPointer<vtkDoubleArray>::New();
	conc_w->SetName("WaterConcentration");
	auto conc_s = vtkSmartPointer<vtkDoubleArray>::New();
	conc_s->SetName("SaltConcentration");
	auto vel_w = vtkSmartPointer<vtkDoubleArray>::New();
	vel_w->SetName("WaterVelocity");		vel_w->SetNumberOfComponents(3);
	auto vel_o = vtkSmartPointer<vtkDoubleArray>::New();
	vel_o->SetName("OilVelocity");			vel_o->SetNumberOfComponents(3);
	auto vel_g = vtkSmartPointer<vtkDoubleArray>::New();
	vel_g->SetName("GasVelocity");			vel_g->SetNumberOfComponents(3);

	int k, j, idx, idx1;

	double vel[3];

	for (j = 0; j < ny - 2; j++)
	{
		idx = j;
		idx1 = j + 1;
		Cell& cell = model->cells[idx1];
		Variable& next = cell.u_next;

		polygon->GetPointIds()->SetId(0, idx);
		polygon->GetPointIds()->SetId(1, idx + ny - 1);
		polygon->GetPointIds()->SetId(2, idx + ny);
		polygon->GetPointIds()->SetId(3, idx + 1);
		polygons->InsertNextCell(polygon);

		poro->InsertNextValue(next.m);
		perm->InsertNextValue(M2toMilliDarcy(cell.props->getPermCoseni_r(next.m).value() * r_dim * r_dim));
		pres->InsertNextValue(next.p * P_dim / BAR_TO_PA);
		p_bub->InsertNextValue(next.p_bub * P_dim / BAR_TO_PA);
		satur->InsertNextValue(next.SATUR);
		sat_w->InsertNextValue(next.sw);
		sat_o->InsertNextValue(next.so);
		sat_g->InsertNextValue(1.0 - next.sw - next.so);
		conc_a->InsertNextValue(next.xa);
		conc_w->InsertNextValue(next.xw);
		conc_s->InsertNextValue(1.0 - next.xw - next.xa);
		vel[0] = 0.0;
		vel[1] = 0.0;
		vel[2] = 0.0;
		vel_w->InsertNextTuple(vel);
		vel_o->InsertNextTuple(vel);
		vel_g->InsertNextTuple(vel);
	}

	// Middle cells
	for (k = 1; k < nx - 1; k++)
	{
		for (j = 0; j < ny - 2; j++)
		{
			idx = k * (ny - 1) + j;
			idx1 = k * ny + j + 1;
			Cell& cell = model->cells[idx1];
			Variable& next = cell.u_next;

			polygon->GetPointIds()->SetId(0, idx);
			polygon->GetPointIds()->SetId(1, idx + ny - 1);
			polygon->GetPointIds()->SetId(2, idx + ny);
			polygon->GetPointIds()->SetId(3, idx + 1);
			polygons->InsertNextCell(polygon);

			poro->InsertNextValue(next.m);
			perm->InsertNextValue(M2toMilliDarcy(cell.props->getPermCoseni_r(next.m).value() * r_dim * r_dim));
			pres->InsertNextValue(next.p * P_dim / BAR_TO_PA);
			p_bub->InsertNextValue(next.p_bub * P_dim / BAR_TO_PA);
			satur->InsertNextValue(next.SATUR);
			sat_w->InsertNextValue(next.sw);
			sat_o->InsertNextValue(next.so);
			sat_g->InsertNextValue(1.0 - next.sw - next.so);
			conc_a->InsertNextValue(next.xa);
			conc_w->InsertNextValue(next.xw);
			conc_s->InsertNextValue(1.0 - next.xw - next.xa);
			vel[0] = r_dim / t_dim * model->getWaterVelocity(cell, R_AXIS);
			vel[1] = r_dim / t_dim * model->getWaterVelocity(cell, Z_AXIS);
			vel[2] = 0.0;
			vel_w->InsertNextTuple(vel);
			vel[0] = r_dim / t_dim * model->getOilVelocity(cell, R_AXIS);
			vel[1] = r_dim / t_dim * model->getOilVelocity(cell, Z_AXIS);
			vel[2] = 0.0;
			vel_o->InsertNextTuple(vel);
			vel[0] = r_dim / t_dim * model->getGasVelocity(cell, R_AXIS);
			vel[1] = r_dim / t_dim * model->getGasVelocity(cell, Z_AXIS);
			vel[2] = 0.0;
			vel_g->InsertNextTuple(vel);
		}
	}

	grid->SetPolys(polygons);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(poro);
	fd->AddArray(perm);
	fd->AddArray(pres);
	fd->AddArray(p_bub);
	fd->AddArray(satur);
	fd->AddArray(sat_w);
	fd->AddArray(sat_o);
	fd->AddArray(sat_g);
	fd->AddArray(conc_a);
	fd->AddArray(conc_w);
	fd->AddArray(conc_s);
	fd->AddArray(vel_w);
	fd->AddArray(vel_o);
	fd->AddArray(vel_g);

	// Writing
	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

	writer->SetFileName(getFileName(i).c_str());
	writer->SetInputData(grid);
	writer->Write();
}
void VTKSnapshotter<acid2dnit::Acid2dNIT>::dump_all(int i)
{
	using namespace acid2dnit;

	// Grid
	auto grid = vtkSmartPointer<vtkPolyData>::New();

	// Points
	auto points = vtkSmartPointer<vtkPoints>::New();
	for (int j = 1; j < ny; j++)
	{
		Cell& cell = model->cells[j];
		points->InsertNextPoint(r_dim * (0.99 * cell.r), -r_dim * (cell.z - cell.hz / 2.0), 0.0);
	}
	for (int k = 1; k < nx; k++)
	{
		for (int j = 1; j < ny; j++)
		{
			Cell& cell = model->cells[k * ny + j];
			points->InsertNextPoint(r_dim * (cell.r - cell.hr / 2.0), -r_dim * (cell.z - cell.hz / 2.0), 0.0);
		}
	}
	grid->SetPoints(points);

	// Data
	auto polygons = vtkSmartPointer<vtkCellArray>::New();
	auto polygon = vtkSmartPointer<vtkPolygon>::New();
	polygon->GetPointIds()->SetNumberOfIds(4);
	auto poro = vtkSmartPointer<vtkDoubleArray>::New();
	poro->SetName("porosity");
	auto perm = vtkSmartPointer<vtkDoubleArray>::New();
	perm->SetName("permeability");
	auto pres = vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");
	auto sat_w = vtkSmartPointer<vtkDoubleArray>::New();
	sat_w->SetName("WaterSaturation");
	auto sat_o = vtkSmartPointer<vtkDoubleArray>::New();
	sat_o->SetName("OilSaturation");
	auto conc_a = vtkSmartPointer<vtkDoubleArray>::New();
	conc_a->SetName("AcidConcentration");
	auto conc_w = vtkSmartPointer<vtkDoubleArray>::New();
	conc_w->SetName("WaterConcentration");
	auto conc_co2 = vtkSmartPointer<vtkDoubleArray>::New();
	conc_co2->SetName("CO2Concentration");
	auto conc_s = vtkSmartPointer<vtkDoubleArray>::New();
	conc_s->SetName("SaltConcentration");
	auto temp = vtkSmartPointer<vtkDoubleArray>::New();
	temp->SetName("temperature");
	auto vel_w = vtkSmartPointer<vtkDoubleArray>::New();
	vel_w->SetName("WaterVelocity");		vel_w->SetNumberOfComponents(3);
	auto vel_o = vtkSmartPointer<vtkDoubleArray>::New();
	vel_o->SetName("OilVelocity");			vel_o->SetNumberOfComponents(3);

	int k, j, idx, idx1;

	double vel[3];

	for (j = 0; j < ny - 2; j++)
	{
		idx = j;
		idx1 = j + 1;
		Cell& cell = model->cells[idx1];
		Variable& next = cell.u_next;

		polygon->GetPointIds()->SetId(0, idx);
		polygon->GetPointIds()->SetId(1, idx + ny - 1);
		polygon->GetPointIds()->SetId(2, idx + ny);
		polygon->GetPointIds()->SetId(3, idx + 1);
		polygons->InsertNextCell(polygon);

		poro->InsertNextValue(next.m);
		perm->InsertNextValue(M2toMilliDarcy(cell.props->getPermCoseni_r(next.m, next.p).value() * r_dim * r_dim));
		pres->InsertNextValue(next.p * P_dim / BAR_TO_PA);
		sat_w->InsertNextValue(next.sw);
		sat_o->InsertNextValue(1.0 - next.sw);
		conc_a->InsertNextValue(next.xa);
		conc_w->InsertNextValue(next.xw);
		conc_s->InsertNextValue(next.xs);
		conc_co2->InsertNextValue(1.0 - next.xw - next.xa - next.xs);
		temp->InsertNextValue(next.t * T_dim);
		vel[0] = 0.0;
		vel[1] = 0.0;
		vel[2] = 0.0;
		vel_w->InsertNextTuple(vel);
		vel_o->InsertNextTuple(vel);
	}

	// Middle cells
	for (k = 1; k < nx - 1; k++)
	{
		for (j = 0; j < ny - 2; j++)
		{
			idx = k * (ny - 1) + j;
			idx1 = k * ny + j + 1;
			Cell& cell = model->cells[idx1];
			Variable& next = cell.u_next;

			polygon->GetPointIds()->SetId(0, idx);
			polygon->GetPointIds()->SetId(1, idx + ny - 1);
			polygon->GetPointIds()->SetId(2, idx + ny);
			polygon->GetPointIds()->SetId(3, idx + 1);
			polygons->InsertNextCell(polygon);

			poro->InsertNextValue(next.m);
			perm->InsertNextValue(M2toMilliDarcy(cell.props->getPermCoseni_r(next.m, next.p).value() * r_dim * r_dim));
			pres->InsertNextValue(next.p * P_dim / BAR_TO_PA);
			sat_w->InsertNextValue(next.sw);
			sat_o->InsertNextValue(1.0 - next.sw);
			conc_a->InsertNextValue(next.xa);
			conc_w->InsertNextValue(next.xw);
			conc_s->InsertNextValue(next.xs);
			conc_co2->InsertNextValue(1.0 - next.xw - next.xa - next.xs);
			temp->InsertNextValue(next.t * T_dim);
			vel[0] = r_dim / t_dim * model->getWatVel(cell, R_AXIS);
			vel[1] = r_dim / t_dim * model->getWatVel(cell, Z_AXIS);
			vel[2] = 0.0;
			vel_w->InsertNextTuple(vel);
			vel[0] = r_dim / t_dim * model->getOilVel(cell, R_AXIS);
			vel[1] = r_dim / t_dim * model->getOilVel(cell, Z_AXIS);
			vel[2] = 0.0;
			vel_o->InsertNextTuple(vel);
		}
	}

	grid->SetPolys(polygons);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(poro);
	fd->AddArray(perm);
	fd->AddArray(pres);
	fd->AddArray(sat_w);
	fd->AddArray(sat_o);
	fd->AddArray(conc_a);
	fd->AddArray(conc_w);
	fd->AddArray(conc_s);
	fd->AddArray(conc_co2);
	fd->AddArray(vel_w);
	fd->AddArray(vel_o);
	fd->AddArray(temp);

	// Writing
	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

	writer->SetFileName(getFileName(i).c_str());
	writer->SetInputData(grid);
	writer->Write();
}
void VTKSnapshotter<acid1d::Acid1d>::dump_all(int i)
{
	using namespace acid1d;

	// Grid
	auto grid = vtkSmartPointer<vtkPolyData>::New();

	// Points
	auto points = vtkSmartPointer<vtkPoints>::New();

	points->InsertNextPoint(r_dim * (0.99 * model->cells[0].x), 0.0, 0.0);
	points->InsertNextPoint(r_dim * (0.99 * model->cells[0].x), -r_dim * model->props_sk.height, 0.0);
	for (int k = 1; k < nx; k++)
	{
		Cell& cell = model->cells[k];
		points->InsertNextPoint(r_dim * (cell.x - cell.hx / 2.0), 0.0, 0.0);
		points->InsertNextPoint(r_dim * (cell.x - cell.hx / 2.0), -r_dim * model->props_sk.height, 0.0);
	}
	grid->SetPoints(points);

	// Data
	auto polygons = vtkSmartPointer<vtkCellArray>::New();
	auto polygon = vtkSmartPointer<vtkPolygon>::New();
	polygon->GetPointIds()->SetNumberOfIds(4);
	auto poro = vtkSmartPointer<vtkDoubleArray>::New();
	poro->SetName("porosity");
	auto perm = vtkSmartPointer<vtkDoubleArray>::New();
	perm->SetName("permeability");
	auto pres = vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");
	auto sat_w = vtkSmartPointer<vtkDoubleArray>::New();
	sat_w->SetName("WaterSaturation");
	auto sat_o = vtkSmartPointer<vtkDoubleArray>::New();
	sat_o->SetName("OilSaturation");
	auto conc_a = vtkSmartPointer<vtkDoubleArray>::New();
	conc_a->SetName("AcidConcentration");
	auto conc_co2 = vtkSmartPointer<vtkDoubleArray>::New();
	conc_co2->SetName("Co2Concentration");
	auto conc_w = vtkSmartPointer<vtkDoubleArray>::New();
	conc_w->SetName("WaterConcentration");
	auto conc_s = vtkSmartPointer<vtkDoubleArray>::New();
	conc_s->SetName("SaltConcentration");
	auto vel_w = vtkSmartPointer<vtkDoubleArray>::New();
	vel_w->SetName("WaterVelocity");		vel_w->SetNumberOfComponents(3);
	auto vel_o = vtkSmartPointer<vtkDoubleArray>::New();
	vel_o->SetName("OilVelocity");			vel_o->SetNumberOfComponents(3);

	int k, j, idx, idx1;

	double vel[3];

	idx = 0;
	Variable& next = model->cells[0].u_next;

	polygon->GetPointIds()->SetId(0, idx);
	polygon->GetPointIds()->SetId(1, idx + ny - 1);
	polygon->GetPointIds()->SetId(2, idx + ny);
	polygon->GetPointIds()->SetId(3, idx + 1);
	polygons->InsertNextCell(polygon);

	poro->InsertNextValue(next.m);
	perm->InsertNextValue(M2toMilliDarcy(model->props_sk.getPermCoseni(next.m).value() * r_dim * r_dim));
	pres->InsertNextValue(next.p * P_dim / BAR_TO_PA);
	sat_w->InsertNextValue(next.sw);
	sat_o->InsertNextValue(1.0 - next.sw);
	conc_co2->InsertNextValue(1.0 - next.xa - next.xw - next.xs);
	conc_a->InsertNextValue(next.xa);
	conc_w->InsertNextValue(next.xw);
	conc_s->InsertNextValue(1.0 - next.xw - next.xa);
	vel[0] = 0.0;
	vel[1] = 0.0;
	vel[2] = 0.0;
	vel_w->InsertNextTuple(vel);
	vel_o->InsertNextTuple(vel);

	// Middle cells
	for (k = 1; k < nx - 1; k++)
	{
			idx = k * (ny - 1);
			idx1 = k;
			Cell& cell = model->cells[idx1];
			Variable& next = cell.u_next;

			polygon->GetPointIds()->SetId(0, idx);
			polygon->GetPointIds()->SetId(1, idx + ny - 1);
			polygon->GetPointIds()->SetId(2, idx + ny);
			polygon->GetPointIds()->SetId(3, idx + 1);
			polygons->InsertNextCell(polygon);

			poro->InsertNextValue(next.m);
			perm->InsertNextValue(M2toMilliDarcy(model->props_sk.getPermCoseni(next.m).value() * r_dim * r_dim));
			pres->InsertNextValue(next.p * P_dim / BAR_TO_PA);
			sat_w->InsertNextValue(next.sw);
			sat_o->InsertNextValue(1.0 - next.sw);
			conc_a->InsertNextValue(next.xa);
			conc_w->InsertNextValue(next.xw);
			conc_co2->InsertNextValue(1.0 - next.xa - next.xw - next.xs);
			conc_s->InsertNextValue(1.0 - next.xw - next.xa);
			vel[0] = r_dim / t_dim * model->getWaterVelocity(cell);
			vel[1] = 0.0;
			vel[2] = 0.0;
			vel_w->InsertNextTuple(vel);
			vel[0] = r_dim / t_dim * model->getOilVelocity(cell);
			vel[1] = 0.0;
			vel[2] = 0.0;
			vel_o->InsertNextTuple(vel);
	}

	grid->SetPolys(polygons);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(poro);
	fd->AddArray(perm);
	fd->AddArray(pres);
	fd->AddArray(sat_w);
	fd->AddArray(sat_o);
	fd->AddArray(conc_a);
	fd->AddArray(conc_w);
	fd->AddArray(conc_s);
	fd->AddArray(conc_co2);
	fd->AddArray(vel_w);
	fd->AddArray(vel_o);

	// Writing
	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

	writer->SetFileName(getFileName(i).c_str());
	writer->SetInputData(grid);
	writer->Write();
}
void VTKSnapshotter<acidfrac::AcidFrac>::dump_all(int i)
{
	using namespace acidfrac;
	const double& w2 = model->props_frac.w2;

	// Grid
	auto grid_frac = vtkSmartPointer<vtkUnstructuredGrid>::New();
	auto grid_poro = vtkSmartPointer<vtkUnstructuredGrid>::New();

	// Points
	const double PORO_MULT_Y = 10;
	const double MULT_X = 0.1;

	auto points_frac = vtkSmartPointer<vtkPoints>::New();
	auto points_poro = vtkSmartPointer<vtkPoints>::New();
	for (int k = 0; k < nz - 1; k++)
	{
		for (int j = 0; j < ny; j++)
		{
			const FracCell& cell = model->cells_frac[j + k * ny];
			const FracCell& xnebr = model->cells_frac[j + k * ny + ny * nz];
			points_frac->InsertNextPoint(-MULT_X * r_dim * xnebr.hx / 10.0, r_dim * (cell.y + cell.hy / 2.0) * FRAC_WIDTH_MULT, r_dim * (cell.z + cell.hz / 2.0));
		}
	}
	for(int i = 0; i < nx - 1; i++)
		for (int k = 0; k < nz - 1; k++)
		{
			for (int j = 0; j < ny; j++)
			{
				FracCell& cell = model->cells_frac[j + k * ny + i * ny * nz];
				points_frac->InsertNextPoint(MULT_X * r_dim * (cell.x + cell.hx / 2.0), r_dim * (cell.y + cell.hy / 2.0) * FRAC_WIDTH_MULT, r_dim * (cell.z + cell.hz / 2.0));
			}
		}
	grid_frac->SetPoints(points_frac);

	const int stop_idx = model->poro_grids[0].cellsNum / 3;
	// Poro points
	for (int i = 0; i < nx - 1; i++)
		for (int k = 0; k < nz - 1; k++)
		{
			if (i != 0 && k != 0)
			{
				auto& frac_cell = model->cells_frac[ny - 1 + k * ny + i * ny * nz];
				const auto& poro_grid = model->poro_grids[model->frac2poro[frac_cell.num]];

				points_poro->InsertNextPoint(MULT_X * r_dim * (frac_cell.x - frac_cell.hx / 2.0), r_dim * ((FRAC_WIDTH_MULT) * w2), r_dim * (frac_cell.z - frac_cell.hz / 2.0));
				points_poro->InsertNextPoint(MULT_X * r_dim * (frac_cell.x + frac_cell.hx / 2.0), r_dim * ((FRAC_WIDTH_MULT) * w2), r_dim * (frac_cell.z - frac_cell.hz / 2.0));
				points_poro->InsertNextPoint(MULT_X * r_dim * (frac_cell.x + frac_cell.hx / 2.0), r_dim * ((FRAC_WIDTH_MULT) * w2), r_dim * (frac_cell.z + frac_cell.hz / 2.0));
				points_poro->InsertNextPoint(MULT_X * r_dim * (frac_cell.x - frac_cell.hx / 2.0), r_dim * ((FRAC_WIDTH_MULT) * w2), r_dim * (frac_cell.z + frac_cell.hz / 2.0));

				double width = poro_grid.cells[1].hx / 10.0 * PORO_MULT_Y;
				for (int j = 0; j < stop_idx + 1; j++)
				{
					const PoroCell& cell = poro_grid.cells[j];
					points_poro->InsertNextPoint(MULT_X * r_dim * (frac_cell.x - frac_cell.hx / 2.0), r_dim * (width + (FRAC_WIDTH_MULT) * w2 + PORO_MULT_Y * cell.hx), r_dim * (frac_cell.z - frac_cell.hz / 2.0));
					points_poro->InsertNextPoint(MULT_X * r_dim * (frac_cell.x + frac_cell.hx / 2.0), r_dim * (width + (FRAC_WIDTH_MULT) * w2 + PORO_MULT_Y * cell.hx), r_dim * (frac_cell.z - frac_cell.hz / 2.0));
					points_poro->InsertNextPoint(MULT_X * r_dim * (frac_cell.x + frac_cell.hx / 2.0), r_dim * (width + (FRAC_WIDTH_MULT) * w2 + PORO_MULT_Y * cell.hx), r_dim * (frac_cell.z + frac_cell.hz / 2.0));
					points_poro->InsertNextPoint(MULT_X * r_dim * (frac_cell.x - frac_cell.hx / 2.0), r_dim * (width + (FRAC_WIDTH_MULT) * w2 + PORO_MULT_Y * cell.hx), r_dim * (frac_cell.z + frac_cell.hz / 2.0));
					width += PORO_MULT_Y * cell.hx;
				}
			}
		}
	grid_poro->SetPoints(points_poro);

	// Data
	auto hexs_frac = vtkSmartPointer<vtkCellArray>::New();
	auto poro_frac = vtkSmartPointer<vtkDoubleArray>::New();
	poro_frac->SetName("porosity");
	auto perm_frac = vtkSmartPointer<vtkDoubleArray>::New();
	perm_frac->SetName("permeability");
	auto pres_frac = vtkSmartPointer<vtkDoubleArray>::New();
	pres_frac->SetName("pressure");
	auto sat_w_frac = vtkSmartPointer<vtkDoubleArray>::New();
	sat_w_frac->SetName("WaterSaturation");
	auto sat_o_frac = vtkSmartPointer<vtkDoubleArray>::New();
	sat_o_frac->SetName("OilSaturation");
	auto conc_a_frac = vtkSmartPointer<vtkDoubleArray>::New();
	conc_a_frac->SetName("AcidConcentration");
	auto conc_w_frac = vtkSmartPointer<vtkDoubleArray>::New();
	conc_w_frac->SetName("WaterConcentration");
	auto conc_co2_frac = vtkSmartPointer<vtkDoubleArray>::New();
	conc_co2_frac->SetName("CO2Concentration");
	auto conc_s_frac = vtkSmartPointer<vtkDoubleArray>::New();
	conc_s_frac->SetName("SaltConcentration");
	auto trans = vtkSmartPointer<vtkDoubleArray>::New();
	trans->SetName("Transmissibility");

	int np = ny * (nz - 1);
	int np_frac = np * nx;
	double vel[3];
	for (int k = 0; k < nz - 2; k++)
	{
		for (int j = 0; j < ny - 1; j++)
		{
			const FracCell& cell = model->cells_frac[j + 1 + (k + 1) * ny];
			const auto& grid0 = model->poro_grids[model->frac2poro[model->cells_frac[model->getRowOuter(cell.num + ny * nz)].num]];
			const auto& props = model->props_sk[0];
			const auto& next = cell.u_next;
			vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();
			poro_frac->InsertNextValue(props.m_init);
			perm_frac->InsertNextValue(M2toMilliDarcy(props.perm * r_dim * r_dim));
			pres_frac->InsertNextValue(next.p * P_dim / BAR_TO_PA);
			sat_w_frac->InsertNextValue(1.0);
			sat_o_frac->InsertNextValue(0.0);
			conc_a_frac->InsertNextValue(next.c);
			conc_w_frac->InsertNextValue(1.0 - next.c);
			conc_s_frac->InsertNextValue(0.0);
			conc_co2_frac->InsertNextValue(0.0);
			trans->InsertNextValue(grid0.trans);

			hex->GetPointIds()->SetId(0, j + k * ny);
			hex->GetPointIds()->SetId(1, j + 1 + k * ny);
			hex->GetPointIds()->SetId(2, j + 1 + k * ny + np);
			hex->GetPointIds()->SetId(3, j + k * ny + np);

			hex->GetPointIds()->SetId(4, j + (k + 1) * ny);
			hex->GetPointIds()->SetId(5, j + 1 + (k + 1) * ny);
			hex->GetPointIds()->SetId(6, j + 1 + (k + 1) * ny + np);
			hex->GetPointIds()->SetId(7, j + (k + 1) * ny + np);

			hexs_frac->InsertNextCell(hex);
		}
	}
	for (int i = 1; i < nx - 1; i++)
	{
		for (int k = 0; k < nz - 2; k++)
		{
			for (int j = 0; j < ny - 1; j++)
			{
				const FracCell& cell = model->cells_frac[j + 1 + (k + 1) * ny + i * nz * ny];
				const auto& next = cell.u_next;
				const auto& props = model->props_sk[0];
				const auto& grid = model->poro_grids[model->frac2poro[model->cells_frac[model->getRowOuter(cell.num)].num]];
				vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();
				poro_frac->InsertNextValue(props.m_init);
				perm_frac->InsertNextValue(M2toMilliDarcy(props.perm * r_dim * r_dim));
				pres_frac->InsertNextValue(next.p * P_dim / BAR_TO_PA);
				sat_w_frac->InsertNextValue(1.0);
				sat_o_frac->InsertNextValue(0.0);
				conc_a_frac->InsertNextValue(next.c);
				conc_w_frac->InsertNextValue(1.0 - next.c);
				conc_s_frac->InsertNextValue(0.0);
				conc_co2_frac->InsertNextValue(0.0);
				trans->InsertNextValue(grid.trans);

				hex->GetPointIds()->SetId(0, j + k * ny + i * np);
				hex->GetPointIds()->SetId(1, j + 1 + k * ny + i * np);
				hex->GetPointIds()->SetId(2, j + 1 + k * ny + (i + 1) * np);
				hex->GetPointIds()->SetId(3, j + k * ny + (i + 1) * np);

				hex->GetPointIds()->SetId(4, j + (k + 1) * ny + i * np);
				hex->GetPointIds()->SetId(5, j + 1 + (k + 1) * ny + i * np);
				hex->GetPointIds()->SetId(6, j + 1 + (k + 1) * ny + (i + 1) * np);
				hex->GetPointIds()->SetId(7, j + (k + 1) * ny + (i + 1) * np);

				hexs_frac->InsertNextCell(hex);
			}
		}
	}

	// Data
	auto hexs_poro = vtkSmartPointer<vtkCellArray>::New();
	auto poro_poro = vtkSmartPointer<vtkDoubleArray>::New();
	poro_poro->SetName("porosity");
	auto perm_poro = vtkSmartPointer<vtkDoubleArray>::New();
	perm_poro->SetName("permeability");
	auto pres_poro = vtkSmartPointer<vtkDoubleArray>::New();
	pres_poro->SetName("pressure");
	auto sat_w_poro = vtkSmartPointer<vtkDoubleArray>::New();
	sat_w_poro->SetName("WaterSaturation");
	auto sat_o_poro = vtkSmartPointer<vtkDoubleArray>::New();
	sat_o_poro->SetName("OilSaturation");
	auto conc_a_poro = vtkSmartPointer<vtkDoubleArray>::New();
	conc_a_poro->SetName("AcidConcentration");
	auto conc_w_poro = vtkSmartPointer<vtkDoubleArray>::New();
	conc_w_poro->SetName("WaterConcentration");
	auto conc_co2_poro = vtkSmartPointer<vtkDoubleArray>::New();
	conc_co2_poro->SetName("CO2Concentration");
	auto conc_s_poro = vtkSmartPointer<vtkDoubleArray>::New();
	conc_s_poro->SetName("SaltConcentration");

	int counter = 0;
	for (int i = 1; i < nx - 1; i++)
	{
		for (int k = 0; k < nz - 2; k++)
		{
			FracCell& frac_cell = model->cells_frac[ny - 1 + (k + 1) * ny + i * nz * ny];
			assert(frac_cell.type == FracType::FRAC_OUT);
			const auto& poro_grid = model->poro_grids[model->frac2poro[frac_cell.num]];
			
			/*const PoroCell& cell = poro_grid.cells[0];
			const auto& next = cell.u_next;
			poro_poro->InsertNextValue(next.m);
			perm_poro->InsertNextValue(M2toMilliDarcy(poro_grid.props_sk->getPermCoseni(next.m, next.p).value() * r_dim * r_dim));
			pres_poro->InsertNextValue(next.p * P_dim / BAR_TO_PA);
			sat_w_poro->InsertNextValue(next.sw);
			sat_o_poro->InsertNextValue(1.0 - next.sw);
			conc_a_poro->InsertNextValue(next.xa);
			conc_w_poro->InsertNextValue(next.xw);
			conc_s_poro->InsertNextValue(next.xs);
			conc_co2_poro->InsertNextValue(1.0 - next.xw - next.xa - next.xs);

			vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();
			hex->GetPointIds()->SetId(0, ny - 1 + k * ny + i * np);
			hex->GetPointIds()->SetId(1, np_frac + counter);
			hex->GetPointIds()->SetId(2, np_frac + counter + 1);
			hex->GetPointIds()->SetId(3, ny - 1 + k * ny + (i + 1) * np);

			hex->GetPointIds()->SetId(4, ny - 1 + (k + 1) * ny + i * np);
			hex->GetPointIds()->SetId(5, np_frac + counter + 3);
			hex->GetPointIds()->SetId(6, np_frac + counter + 2);
			hex->GetPointIds()->SetId(7, ny - 1 + (k + 1) * ny + (i + 1) * np);
			hexs_poro->InsertNextCell(hex);*/

			for (int j = 1; j < stop_idx + 1; j++)
			{
				const PoroCell& cell = poro_grid.cells[j-1];
				const auto& next = cell.u_next;
				poro_poro->InsertNextValue(next.m);
				perm_poro->InsertNextValue(M2toMilliDarcy(poro_grid.props_sk->getPermCoseni(next.m, next.p).value() * r_dim * r_dim));
				pres_poro->InsertNextValue(next.p * P_dim / BAR_TO_PA);
				sat_w_poro->InsertNextValue(next.sw);
				sat_o_poro->InsertNextValue(1.0 - next.sw);
				conc_a_poro->InsertNextValue(next.xa);
				conc_w_poro->InsertNextValue(next.xw);
				conc_s_poro->InsertNextValue(next.xs);
				conc_co2_poro->InsertNextValue(1.0 - next.xw - next.xa - next.xs);

				vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();
				hex->GetPointIds()->SetId(0, counter + 4 * (j - 1));
				hex->GetPointIds()->SetId(1, counter + 4 * (j - 1) + 1);
				hex->GetPointIds()->SetId(2, counter + 4 * j + 1);
				hex->GetPointIds()->SetId(3, counter + 4 * j);

				hex->GetPointIds()->SetId(4, counter + 4 * (j - 1) + 3);
				hex->GetPointIds()->SetId(5, counter + 4 * (j - 1) + 2);
				hex->GetPointIds()->SetId(6, counter + 4 * j + 2);
				hex->GetPointIds()->SetId(7, counter + 4 * j + 3);
				hexs_poro->InsertNextCell(hex);
			}

			counter += 4 * (stop_idx + 2);
		}
	}

	grid_frac->SetCells(VTK_HEXAHEDRON, hexs_frac);
	vtkCellData* fd_frac = grid_frac->GetCellData();
	fd_frac->AddArray(poro_frac);
	fd_frac->AddArray(perm_frac);
	fd_frac->AddArray(pres_frac);
	fd_frac->AddArray(sat_w_frac);
	fd_frac->AddArray(sat_o_frac);
	fd_frac->AddArray(conc_a_frac);
	fd_frac->AddArray(conc_w_frac);
	fd_frac->AddArray(conc_s_frac);
	fd_frac->AddArray(conc_co2_frac);
	fd_frac->AddArray(trans);

	grid_poro->SetCells(VTK_HEXAHEDRON, hexs_poro);
	vtkCellData* fd_poro = grid_poro->GetCellData();
	fd_poro->AddArray(poro_poro);
	fd_poro->AddArray(perm_poro);
	fd_poro->AddArray(pres_poro);
	fd_poro->AddArray(sat_w_poro);
	fd_poro->AddArray(sat_o_poro);
	fd_poro->AddArray(conc_a_poro);
	fd_poro->AddArray(conc_w_poro);
	fd_poro->AddArray(conc_s_poro);
	fd_poro->AddArray(conc_co2_poro);

	// Writing
	auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(getFileName(i, "frac").c_str());
	writer->SetInputData(grid_frac);
	writer->Write();

	writer->SetFileName(getFileName(i, "poro").c_str());
	writer->SetInputData(grid_poro);
	writer->Write();
}
void VTKSnapshotter<acidellfrac::AcidEllFrac>::dump_all(int i)
{
	using namespace acidellfrac;
	const double& w2 = model->props_frac.w2;
	const double MU_FRAC_MULT = 1.0;
	const double X_MULT = 0.05;
	const double MAX_MU = model->cells_poro.back().c.mu / 15.0;
	const double MU_0 = model->cells_poro[0].c.mu;
	// Grid
	auto grid_frac = vtkSmartPointer<vtkUnstructuredGrid>::New();
	auto grid_poro = vtkSmartPointer<vtkUnstructuredGrid>::New();
	// Frac points
	auto points_frac = vtkSmartPointer<vtkPoints>::New();
	for (int k = 0; k < nz - 1; k++)
	{
		for (int j = 0; j < ny; j++)
		{
			const FracCell& cell = model->cells_frac[j + k * ny];
			const FracCell& xnebr = model->cells_frac[j + k * ny + ny * nz];
			const auto cart_point = point::getCartesianFromElliptic(MU_FRAC_MULT * (cell.c.mu + cell.h.mu / 2.0), cell.c.nu + xnebr.h.nu / 20.0, cell.c.z + cell.h.z / 2.0);
			points_frac->InsertNextPoint(r_dim * X_MULT * cart_point.x, r_dim * cart_point.y, r_dim * cart_point.z);
		}
	}
	for (int i = 0; i < nx - 1; i++)
		for (int k = 0; k < nz - 1; k++)
		{
			for (int j = 0; j < ny; j++)
			{
				const FracCell& cell = model->cells_frac[j + k * ny + i * ny * nz];
				const auto cart_point = point::getCartesianFromElliptic(MU_FRAC_MULT * (cell.c.mu + cell.h.mu / 2.0), cell.c.nu - cell.h.nu / 2.0, cell.c.z + cell.h.z / 2.0);
				points_frac->InsertNextPoint(r_dim * X_MULT * cart_point.x, r_dim * cart_point.y, r_dim * cart_point.z);
			}
		}
	grid_frac->SetPoints(points_frac);
	// Poro points
	auto points_poro = vtkSmartPointer<vtkPoints>::New();
	const int ny_poro = model->cellsNum_mu_poro + 2;
	point::CartesianPoint3d point;
	for (int i = 0; i < nx - 1; i++)
		for (int k = 0; k < nz - 1; k++)
		{
			for (int j = 0; j < ny_poro; j++)
			{
				const PoroCell& cell = model->cells_poro[j + k * ny_poro + i * ny_poro * nz];
				if (j != 1)
					point = point::getCartesianFromElliptic(MU_FRAC_MULT * MU_0 + (cell.c.mu - cell.h.mu / 2.0) - MU_0, cell.c.nu - cell.h.nu / 2.0, cell.c.z + cell.h.z / 2.0);
				else
					point = point::getCartesianFromElliptic(MU_FRAC_MULT * MU_0 + (cell.c.mu - 2.0 * cell.h.mu / 5.0) - MU_0, cell.c.nu - cell.h.nu / 2.0, cell.c.z + cell.h.z / 2.0);
				points_poro->InsertNextPoint(r_dim * X_MULT * point.x, r_dim * point.y, r_dim * point.z);
			}
		}
	grid_poro->SetPoints(points_poro);
	/*for (int i = 0; i < nx - 1; i++)
		for (int k = 0; k < nz - 1; k++)
		{
			if (i != 0 && k != 0)
			{
				auto& frac_cell = model->cells_frac[ny - 1 + k * ny + i * ny * nz];
				const auto& poro_grid = model->poro_grids[model->frac2poro[frac_cell.num]];

				points_poro->InsertNextPoint(MULT_X * r_dim * (frac_cell.x - frac_cell.hx / 2.0), r_dim * ((FRAC_WIDTH_MULT)* w2), r_dim * (frac_cell.z - frac_cell.hz / 2.0));
				points_poro->InsertNextPoint(MULT_X * r_dim * (frac_cell.x + frac_cell.hx / 2.0), r_dim * ((FRAC_WIDTH_MULT)* w2), r_dim * (frac_cell.z - frac_cell.hz / 2.0));
				points_poro->InsertNextPoint(MULT_X * r_dim * (frac_cell.x + frac_cell.hx / 2.0), r_dim * ((FRAC_WIDTH_MULT)* w2), r_dim * (frac_cell.z + frac_cell.hz / 2.0));
				points_poro->InsertNextPoint(MULT_X * r_dim * (frac_cell.x - frac_cell.hx / 2.0), r_dim * ((FRAC_WIDTH_MULT)* w2), r_dim * (frac_cell.z + frac_cell.hz / 2.0));

				double width = poro_grid.cells[1].hx / 10.0 * PORO_MULT_Y;
				for (int j = 0; j < stop_idx + 1; j++)
				{
					const PoroCell& cell = poro_grid.cells[j];
					points_poro->InsertNextPoint(MULT_X * r_dim * (frac_cell.x - frac_cell.hx / 2.0), r_dim * (width + (FRAC_WIDTH_MULT)* w2 + PORO_MULT_Y * cell.hx), r_dim * (frac_cell.z - frac_cell.hz / 2.0));
					points_poro->InsertNextPoint(MULT_X * r_dim * (frac_cell.x + frac_cell.hx / 2.0), r_dim * (width + (FRAC_WIDTH_MULT)* w2 + PORO_MULT_Y * cell.hx), r_dim * (frac_cell.z - frac_cell.hz / 2.0));
					points_poro->InsertNextPoint(MULT_X * r_dim * (frac_cell.x + frac_cell.hx / 2.0), r_dim * (width + (FRAC_WIDTH_MULT)* w2 + PORO_MULT_Y * cell.hx), r_dim * (frac_cell.z + frac_cell.hz / 2.0));
					points_poro->InsertNextPoint(MULT_X * r_dim * (frac_cell.x - frac_cell.hx / 2.0), r_dim * (width + (FRAC_WIDTH_MULT)* w2 + PORO_MULT_Y * cell.hx), r_dim * (frac_cell.z + frac_cell.hz / 2.0));
					width += PORO_MULT_Y * cell.hx;
				}
			}
		}*/

	// Data
	auto hexs_frac = vtkSmartPointer<vtkCellArray>::New();
	auto poro_frac = vtkSmartPointer<vtkDoubleArray>::New();
	poro_frac->SetName("porosity");
	auto perm_frac = vtkSmartPointer<vtkDoubleArray>::New();
	perm_frac->SetName("permeability");
	auto pres_frac = vtkSmartPointer<vtkDoubleArray>::New();
	pres_frac->SetName("pressure");
	auto sat_w_frac = vtkSmartPointer<vtkDoubleArray>::New();
	sat_w_frac->SetName("WaterSaturation");
	auto sat_o_frac = vtkSmartPointer<vtkDoubleArray>::New();
	sat_o_frac->SetName("OilSaturation");
	auto conc_a_frac = vtkSmartPointer<vtkDoubleArray>::New();
	conc_a_frac->SetName("AcidConcentration");
	auto conc_w_frac = vtkSmartPointer<vtkDoubleArray>::New();
	conc_w_frac->SetName("WaterConcentration");
	auto conc_co2_frac = vtkSmartPointer<vtkDoubleArray>::New();
	conc_co2_frac->SetName("CO2Concentration");
	auto conc_s_frac = vtkSmartPointer<vtkDoubleArray>::New();
	conc_s_frac->SetName("SaltConcentration");
	//auto trans = vtkSmartPointer<vtkDoubleArray>::New();
	//trans->SetName("Transmissibility");
	auto type = vtkSmartPointer<vtkDoubleArray>::New();
	type->SetName("Type");

	int np = ny * (nz - 1);
	double vel[3];
	for (int k = 0; k < nz - 2; k++)
	{
		for (int j = 0; j < ny - 1; j++)
		{
			const FracCell& cell = model->cells_frac[j + 1 + (k + 1) * ny];
			assert(cell.type == FracType::FRAC_IN);
			const auto& props = model->props_sk[0];
			const auto& next = cell.u_next;
			vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();
			poro_frac->InsertNextValue(props.m_init);
			perm_frac->InsertNextValue(M2toMilliDarcy(props.perm * r_dim * r_dim));
			pres_frac->InsertNextValue(next.p * P_dim / BAR_TO_PA);
			sat_w_frac->InsertNextValue(1.0);
			sat_o_frac->InsertNextValue(0.0);
			conc_a_frac->InsertNextValue(next.c);
			conc_w_frac->InsertNextValue(1.0 - next.c);
			conc_s_frac->InsertNextValue(0.0);
			conc_co2_frac->InsertNextValue(0.0);
			//trans->InsertNextValue(grid0.trans);
			type->InsertNextValue(cell.type);

			hex->GetPointIds()->SetId(0, j + k * ny);
			hex->GetPointIds()->SetId(1, j + 1 + k * ny);
			hex->GetPointIds()->SetId(2, j + 1 + k * ny + np);
			hex->GetPointIds()->SetId(3, j + k * ny + np);

			hex->GetPointIds()->SetId(4, j + (k + 1) * ny);
			hex->GetPointIds()->SetId(5, j + 1 + (k + 1) * ny);
			hex->GetPointIds()->SetId(6, j + 1 + (k + 1) * ny + np);
			hex->GetPointIds()->SetId(7, j + (k + 1) * ny + np);

			hexs_frac->InsertNextCell(hex);
		}
	}
	for (int i = 1; i < nx - 1; i++)
	{
		for (int k = 0; k < nz - 2; k++)
		{
			for (int j = 0; j < ny - 1; j++)
			{
				const FracCell& cell = model->cells_frac[j + 1 + (k + 1) * ny + i * nz * ny];
				assert(cell.type == FracType::FRAC_MID || cell.type == FracType::FRAC_OUT);
				const auto& next = cell.u_next;
				const auto& props = model->props_sk[0];
				vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();
				poro_frac->InsertNextValue(props.m_init);
				perm_frac->InsertNextValue(M2toMilliDarcy(props.perm * r_dim * r_dim));
				pres_frac->InsertNextValue(next.p * P_dim / BAR_TO_PA);
				sat_w_frac->InsertNextValue(1.0);
				sat_o_frac->InsertNextValue(0.0);
				conc_a_frac->InsertNextValue(next.c);
				conc_w_frac->InsertNextValue(1.0 - next.c);
				conc_s_frac->InsertNextValue(0.0);
				conc_co2_frac->InsertNextValue(0.0);
				//trans->InsertNextValue(grid.trans);
				type->InsertNextValue(cell.type);

				hex->GetPointIds()->SetId(0, j + k * ny + i * np);
				hex->GetPointIds()->SetId(1, j + 1 + k * ny + i * np);
				hex->GetPointIds()->SetId(2, j + 1 + k * ny + (i + 1) * np);
				hex->GetPointIds()->SetId(3, j + k * ny + (i + 1) * np);

				hex->GetPointIds()->SetId(4, j + (k + 1) * ny + i * np);
				hex->GetPointIds()->SetId(5, j + 1 + (k + 1) * ny + i * np);
				hex->GetPointIds()->SetId(6, j + 1 + (k + 1) * ny + (i + 1) * np);
				hex->GetPointIds()->SetId(7, j + (k + 1) * ny + (i + 1) * np);

				hexs_frac->InsertNextCell(hex);
			}
		}
	}

	// Data
	auto hexs_poro = vtkSmartPointer<vtkCellArray>::New();
	auto poro_poro = vtkSmartPointer<vtkDoubleArray>::New();
	poro_poro->SetName("porosity");
	auto perm_poro = vtkSmartPointer<vtkDoubleArray>::New();
	perm_poro->SetName("permeability");
	auto pres_poro = vtkSmartPointer<vtkDoubleArray>::New();
	pres_poro->SetName("pressure");
	auto sat_w_poro = vtkSmartPointer<vtkDoubleArray>::New();
	sat_w_poro->SetName("WaterSaturation");
	auto sat_o_poro = vtkSmartPointer<vtkDoubleArray>::New();
	sat_o_poro->SetName("OilSaturation");
	auto conc_a_poro = vtkSmartPointer<vtkDoubleArray>::New();
	conc_a_poro->SetName("AcidConcentration");
	auto conc_w_poro = vtkSmartPointer<vtkDoubleArray>::New();
	conc_w_poro->SetName("WaterConcentration");
	auto conc_co2_poro = vtkSmartPointer<vtkDoubleArray>::New();
	conc_co2_poro->SetName("CO2Concentration");
	auto conc_s_poro = vtkSmartPointer<vtkDoubleArray>::New();
	conc_s_poro->SetName("SaltConcentration");

	int np_poro = ny_poro * (nz - 1);
	const auto& props_sk = model->props_sk[0];
	for (int i = 0; i < nx - 2; i++)
	{
		for (int k = 0; k < nz - 2; k++)
		{
			for (int j = 0; j < ny_poro - 1; j++)
			{
				const PoroCell& cell = model->cells_poro[j + (k + 1) * ny_poro + (i + 1) * nz * ny_poro];
				if (cell.c.mu < MAX_MU)
				{
					assert(cell.type == PoroType::MIDDLE || (cell.type == PoroType::WELL_LAT && j == 0));
					const auto& next = cell.u_next;
					vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();

					poro_poro->InsertNextValue(next.m);
					perm_poro->InsertNextValue(M2toMilliDarcy(props_sk.getPermCoseni(next.m, next.p).value() * r_dim * r_dim));
					pres_poro->InsertNextValue(next.p * P_dim / BAR_TO_PA);
					sat_w_poro->InsertNextValue(next.sw);
					sat_o_poro->InsertNextValue(1.0 - next.sw);
					conc_a_poro->InsertNextValue(next.xa);
					conc_w_poro->InsertNextValue(next.xw);
					conc_s_poro->InsertNextValue(next.xs);
					conc_co2_poro->InsertNextValue(1.0 - next.xw - next.xa - next.xs);

					hex->GetPointIds()->SetId(0, j + k * ny_poro + i * np_poro);
					hex->GetPointIds()->SetId(1, j + 1 + k * ny_poro + i * np_poro);
					hex->GetPointIds()->SetId(2, j + 1 + k * ny_poro + (i + 1) * np_poro);
					hex->GetPointIds()->SetId(3, j + k * ny_poro + (i + 1) * np_poro);

					hex->GetPointIds()->SetId(4, j + (k + 1) * ny_poro + i * np_poro);
					hex->GetPointIds()->SetId(5, j + 1 + (k + 1) * ny_poro + i * np_poro);
					hex->GetPointIds()->SetId(6, j + 1 + (k + 1) * ny_poro + (i + 1) * np_poro);
					hex->GetPointIds()->SetId(7, j + (k + 1) * ny_poro + (i + 1) * np_poro);

					hexs_poro->InsertNextCell(hex);
				}
			}
		}
	}

	grid_frac->SetCells(VTK_HEXAHEDRON, hexs_frac);
	vtkCellData* fd_frac = grid_frac->GetCellData();
	fd_frac->AddArray(poro_frac);
	fd_frac->AddArray(perm_frac);
	fd_frac->AddArray(pres_frac);
	fd_frac->AddArray(sat_w_frac);
	fd_frac->AddArray(sat_o_frac);
	fd_frac->AddArray(conc_a_frac);
	fd_frac->AddArray(conc_w_frac);
	fd_frac->AddArray(conc_s_frac);
	fd_frac->AddArray(conc_co2_frac);
	//fd_frac->AddArray(trans);
	fd_frac->AddArray(type);

	grid_poro->SetCells(VTK_HEXAHEDRON, hexs_poro);
	vtkCellData* fd_poro = grid_poro->GetCellData();
	fd_poro->AddArray(poro_poro);
	fd_poro->AddArray(perm_poro);
	fd_poro->AddArray(pres_poro);
	fd_poro->AddArray(sat_w_poro);
	fd_poro->AddArray(sat_o_poro);
	fd_poro->AddArray(conc_a_poro);
	fd_poro->AddArray(conc_w_poro);
	fd_poro->AddArray(conc_s_poro);
	fd_poro->AddArray(conc_co2_poro);

	// Writing
	auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(getFileName(i, "frac").c_str());
	writer->SetInputData(grid_frac);
	writer->Write();
	writer->SetFileName(getFileName(i, "poro").c_str());
	writer->SetInputData(grid_poro);
	writer->Write();
}
void VTKSnapshotter<acidrecfrac::AcidRecFrac>::dump_all(int i)
{
	using namespace acidrecfrac;
	const double& w2 = model->props_frac.w2;
	const double Y_MULT = 1.0;
	// Grid
	auto grid_frac = vtkSmartPointer<vtkUnstructuredGrid>::New();
	auto grid_poro = vtkSmartPointer<vtkUnstructuredGrid>::New();
	// Frac points
	auto points_frac = vtkSmartPointer<vtkPoints>::New();
	for (int k = 0; k < nz - 1; k++)
	{
		for (int j = 0; j < ny; j++)
		{
			const FracCell& cell = model->cells_frac[j + k * ny];
			const FracCell& xnebr = model->cells_frac[j + k * ny + ny * nz];
			points_frac->InsertNextPoint(r_dim * (cell.x - xnebr.hx / 30.0), Y_MULT * r_dim * (cell.y + cell.hy / 2.0), r_dim * (cell.z + cell.hz / 2.0));
		}
	}
	for (int i = 0; i < nx - 1; i++)
		for (int k = 0; k < nz - 1; k++)
		{
			for (int j = 0; j < ny; j++)
			{
				const FracCell& cell = model->cells_frac[j + k * ny + i * ny * nz];
				points_frac->InsertNextPoint(r_dim * (cell.x + cell.hx / 2.0), Y_MULT * r_dim * (cell.y + cell.hy / 2.0), r_dim * (cell.z + cell.hz / 2.0));
			}
		}
	grid_frac->SetPoints(points_frac);
	// Poro points
	auto points_poro = vtkSmartPointer<vtkPoints>::New();
	const int ny_poro = model->cellsNum_y_poro + 2;
	for (int i = 0; i < nx - 1; i++)
		for (int k = 0; k < nz - 1; k++)
		{
			for (int j = 0; j < ny_poro; j++)
			{
				const PoroCell& cell = model->cells_poro[j + k * ny_poro + i * ny_poro * nz];
				if (j != 1)
					points_poro->InsertNextPoint(r_dim * (cell.x + cell.hx / 2.0), Y_MULT * r_dim * (cell.y - cell.hy / 2.0), r_dim * (cell.z + cell.hz / 2.0));
				else
					points_poro->InsertNextPoint(r_dim * (cell.x + cell.hx / 2.0), Y_MULT * r_dim * (cell.y - 2.0 * cell.hy / 200.0), r_dim * (cell.z + cell.hz / 2.0));
			}
		}
	grid_poro->SetPoints(points_poro);
	// Data
	auto hexs_frac = vtkSmartPointer<vtkCellArray>::New();
	auto poro_frac = vtkSmartPointer<vtkDoubleArray>::New();
	poro_frac->SetName("porosity");
	auto perm_frac = vtkSmartPointer<vtkDoubleArray>::New();
	perm_frac->SetName("permeability");
	auto pres_frac = vtkSmartPointer<vtkDoubleArray>::New();
	pres_frac->SetName("pressure");
	auto sat_w_frac = vtkSmartPointer<vtkDoubleArray>::New();
	sat_w_frac->SetName("WaterSaturation");
	auto sat_o_frac = vtkSmartPointer<vtkDoubleArray>::New();
	sat_o_frac->SetName("OilSaturation");
	auto conc_a_frac = vtkSmartPointer<vtkDoubleArray>::New();
	conc_a_frac->SetName("AcidConcentration");
	auto conc_w_frac = vtkSmartPointer<vtkDoubleArray>::New();
	conc_w_frac->SetName("WaterConcentration");
	auto conc_co2_frac = vtkSmartPointer<vtkDoubleArray>::New();
	conc_co2_frac->SetName("CO2Concentration");
	auto conc_s_frac = vtkSmartPointer<vtkDoubleArray>::New();
	conc_s_frac->SetName("SaltConcentration");
	auto trans = vtkSmartPointer<vtkDoubleArray>::New();
	trans->SetName("Transmissibility");
	auto v_leak = vtkSmartPointer<vtkDoubleArray>::New();
	v_leak->SetName("Leakoff_Velocity");
	auto width = vtkSmartPointer<vtkDoubleArray>::New();
	width->SetName("Width");
	auto type = vtkSmartPointer<vtkDoubleArray>::New();
	type->SetName("Type");
    //auto reaction_frac = vtkSmartPointer<vtkDoubleArray>::New();
   // reaction_frac->SetName("ReactionSpeed");

	int np = ny * (nz - 1);
	std::array<double,3> vel;
	double buf;
	for (int k = 0; k < nz - 2; k++)
	{
		for (int j = 0; j < ny - 1; j++)
		{
			const FracCell& cell = model->cells_frac[j + 1 + (k + 1) * ny];
			assert(cell.type == FracType::FRAC_IN);
			const auto& props = model->props_sk[0];
			const auto& next = cell.u_next;
			vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();
			poro_frac->InsertNextValue(props.m_init);
			perm_frac->InsertNextValue(M2toMilliDarcy(props.perm * r_dim * r_dim));
			pres_frac->InsertNextValue(next.p * P_dim / BAR_TO_PA);
			sat_w_frac->InsertNextValue(1.0);
			sat_o_frac->InsertNextValue(0.0);
			conc_a_frac->InsertNextValue(next.c);
			conc_w_frac->InsertNextValue(1.0 - next.c);
			conc_s_frac->InsertNextValue(0.0);
			conc_co2_frac->InsertNextValue(0.0);
            trans->InsertNextValue(model->trans[k] * model->widths[k] * model->R_dim *
                                M2toMilliDarcy(props.perm) * model->R_dim * model->R_dim);
			width->InsertNextValue(model->widths[k] * r_dim);
			v_leak->InsertNextValue(0.0);
			type->InsertNextValue(cell.type);
            //reaction_frac->InsertNextValue(0.0);

			hex->GetPointIds()->SetId(0, j + k * ny);
			hex->GetPointIds()->SetId(1, j + 1 + k * ny);
			hex->GetPointIds()->SetId(2, j + 1 + k * ny + np);
			hex->GetPointIds()->SetId(3, j + k * ny + np);

			hex->GetPointIds()->SetId(4, j + (k + 1) * ny);
			hex->GetPointIds()->SetId(5, j + 1 + (k + 1) * ny);
			hex->GetPointIds()->SetId(6, j + 1 + (k + 1) * ny + np);
			hex->GetPointIds()->SetId(7, j + (k + 1) * ny + np);

			hexs_frac->InsertNextCell(hex);
		}
	}
	for (int i = 1; i < nx - 1; i++)
	{
		for (int k = 0; k < nz - 2; k++)
		{
			for (int j = 0; j < ny - 1; j++)
			{
				const FracCell& cell = model->cells_frac[j + 1 + (k + 1) * ny + i * nz * ny];
				assert(cell.type == FracType::FRAC_MID || cell.type == FracType::FRAC_OUT);
				const auto& next = cell.u_next;
				const auto& props = model->props_sk[0];
				vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();
				poro_frac->InsertNextValue(props.m_init);
				perm_frac->InsertNextValue(M2toMilliDarcy(props.perm * r_dim * r_dim));
				pres_frac->InsertNextValue(next.p * P_dim / BAR_TO_PA);
				sat_w_frac->InsertNextValue(1.0);
				sat_o_frac->InsertNextValue(0.0);
				conc_a_frac->InsertNextValue(next.c);
				conc_w_frac->InsertNextValue(1.0 - next.c);
				conc_s_frac->InsertNextValue(0.0);
				conc_co2_frac->InsertNextValue(0.0);
				trans->InsertNextValue(model->trans[k + (i - 1) * model->cellsNum_z] * model->widths[k + (i - 1) * model->cellsNum_z] * model->R_dim *
                                            M2toMilliDarcy(props.perm) * model->R_dim * model->R_dim);
				width->InsertNextValue(model->widths[k + (i - 1) * model->cellsNum_z] * r_dim);
				v_leak->InsertNextValue(model->getFlowLeak(cell).value() * r_dim / t_dim);
				type->InsertNextValue(cell.type);

				hex->GetPointIds()->SetId(0, j + k * ny + i * np);
				hex->GetPointIds()->SetId(1, j + 1 + k * ny + i * np);
				hex->GetPointIds()->SetId(2, j + 1 + k * ny + (i + 1) * np);
				hex->GetPointIds()->SetId(3, j + k * ny + (i + 1) * np);

				hex->GetPointIds()->SetId(4, j + (k + 1) * ny + i * np);
				hex->GetPointIds()->SetId(5, j + 1 + (k + 1) * ny + i * np);
				hex->GetPointIds()->SetId(6, j + 1 + (k + 1) * ny + (i + 1) * np);
				hex->GetPointIds()->SetId(7, j + (k + 1) * ny + (i + 1) * np);
				hexs_frac->InsertNextCell(hex);
			}
		}
	}

	// Data
	auto hexs_poro = vtkSmartPointer<vtkCellArray>::New();
	auto poro_poro = vtkSmartPointer<vtkDoubleArray>::New();
	poro_poro->SetName("porosity");
	auto perm_poro = vtkSmartPointer<vtkDoubleArray>::New();
	perm_poro->SetName("permeability");
	auto pres_poro = vtkSmartPointer<vtkDoubleArray>::New();
	pres_poro->SetName("pressure");
	auto sat_w_poro = vtkSmartPointer<vtkDoubleArray>::New();
	sat_w_poro->SetName("WaterSaturation");
	auto sat_o_poro = vtkSmartPointer<vtkDoubleArray>::New();
	sat_o_poro->SetName("OilSaturation");
	auto conc_a_poro = vtkSmartPointer<vtkDoubleArray>::New();
	conc_a_poro->SetName("AcidConcentration");
	auto conc_w_poro = vtkSmartPointer<vtkDoubleArray>::New();
	conc_w_poro->SetName("WaterConcentration");
	auto conc_co2_poro = vtkSmartPointer<vtkDoubleArray>::New();
	conc_co2_poro->SetName("CO2Concentration");
	auto conc_s_poro = vtkSmartPointer<vtkDoubleArray>::New();
	conc_s_poro->SetName("SaltConcentration");
    auto trans_poro = vtkSmartPointer<vtkDoubleArray>::New();
    trans_poro->SetName("Transmissibility");
    auto reaction_poro = vtkSmartPointer<vtkDoubleArray>::New();
    reaction_poro->SetName("ReactionSpeed");
    auto vel_poro = vtkSmartPointer<vtkDoubleArray>::New();
    vel_poro->SetName("Velocity");
    vel_poro->SetNumberOfComponents(3);
	auto darmkoller = vtkSmartPointer<vtkDoubleArray>::New();
	darmkoller->SetName("Darmkoller");

    double sum_width, sum_trans;
	int np_poro = ny_poro * (nz - 1);
	for (int i = 0; i < nx - 2; i++)
	{
		for (int k = 0; k < nz - 2; k++)
		{
            sum_width = sum_trans = 0.0;
			for (int j = 0; j < ny_poro - 1; j++)
			{
				const PoroCell& cell = model->cells_poro[j + (k + 1) * ny_poro + (i + 1) * nz * ny_poro];
				if (cell.y * r_dim < 0.5)
				{
					assert(cell.type == PoroType::MIDDLE || (cell.type == PoroType::WELL_LAT && j == 0));
					const auto& next = cell.u_next;
					vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();

					poro_poro->InsertNextValue(next.m);
					perm_poro->InsertNextValue(M2toMilliDarcy(cell.props->getPermCoseni(next.m, next.p).value() * r_dim * r_dim));
					pres_poro->InsertNextValue(next.p * P_dim / BAR_TO_PA);
					sat_w_poro->InsertNextValue(next.sw);
					sat_o_poro->InsertNextValue(1.0 - next.sw);
					conc_a_poro->InsertNextValue(next.xa);
					conc_w_poro->InsertNextValue(next.xw);
					conc_s_poro->InsertNextValue(next.xs);
					conc_co2_poro->InsertNextValue(1.0 - next.xw - next.xa - next.xs);

                    sum_width += cell.hy * r_dim;
                    sum_trans += M2toMilliDarcy(cell.props->getPermCoseni(next.m, next.p).value() * r_dim * r_dim) * cell.hy * r_dim;
                    trans_poro->InsertNextValue(sum_trans);

                    if (cell.type != PoroType::WELL_LAT)
                    {
                        buf = model->getReactionRateOutput(next, *cell.props).value();
                        reaction_poro->InsertNextValue(buf / r_dim / r_dim / r_dim / t_dim);
                        vel = model->getPoroWaterVelocity(cell);
                        vel[0] *= r_dim / t_dim;    vel[1] *= r_dim / t_dim;    vel[2] *= r_dim / t_dim;
                        vel_poro->InsertNextTuple(&vel[0]);
						darmkoller->InsertNextValue(model->getDarmkoller(cell, next, *cell.props));
                    }
                    else
                    {
						buf = 0.0;
                        reaction_poro->InsertNextValue(buf / r_dim / r_dim / r_dim / t_dim);
                        vel[0] = 0.0;   vel[2] = 0.0;
                        vel[1] = model->getFlowLeak(model->cells_frac[model->getFracNebr(cell.num)]).value() * r_dim / t_dim;
                        vel_poro->InsertNextTuple(&vel[0]);
						darmkoller->InsertNextValue(0.0);
                    }

					hex->GetPointIds()->SetId(0, j + k * ny_poro + i * np_poro);
					hex->GetPointIds()->SetId(1, j + 1 + k * ny_poro + i * np_poro);
					hex->GetPointIds()->SetId(2, j + 1 + k * ny_poro + (i + 1) * np_poro);
					hex->GetPointIds()->SetId(3, j + k * ny_poro + (i + 1) * np_poro);

					hex->GetPointIds()->SetId(4, j + (k + 1) * ny_poro + i * np_poro);
					hex->GetPointIds()->SetId(5, j + 1 + (k + 1) * ny_poro + i * np_poro);
					hex->GetPointIds()->SetId(6, j + 1 + (k + 1) * ny_poro + (i + 1) * np_poro);
					hex->GetPointIds()->SetId(7, j + (k + 1) * ny_poro + (i + 1) * np_poro);
					hexs_poro->InsertNextCell(hex);
				}
			}
		}
	 }

	grid_frac->SetCells(VTK_HEXAHEDRON, hexs_frac);
	vtkCellData* fd_frac = grid_frac->GetCellData();
	fd_frac->AddArray(poro_frac);
	fd_frac->AddArray(perm_frac);
	fd_frac->AddArray(pres_frac);
	fd_frac->AddArray(sat_w_frac);
	fd_frac->AddArray(sat_o_frac);
	fd_frac->AddArray(conc_a_frac);
	fd_frac->AddArray(conc_w_frac);
	fd_frac->AddArray(conc_s_frac);
	fd_frac->AddArray(conc_co2_frac);
	fd_frac->AddArray(trans);
	fd_frac->AddArray(width);
	fd_frac->AddArray(v_leak);
	fd_frac->AddArray(type);
    //fd_frac->AddArray(reaction_frac);

	grid_poro->SetCells(VTK_HEXAHEDRON, hexs_poro);
	vtkCellData* fd_poro = grid_poro->GetCellData();
	fd_poro->AddArray(poro_poro);
	fd_poro->AddArray(perm_poro);
	fd_poro->AddArray(pres_poro);
	fd_poro->AddArray(sat_w_poro);
	fd_poro->AddArray(sat_o_poro);
	fd_poro->AddArray(conc_a_poro);
	fd_poro->AddArray(conc_w_poro);
	fd_poro->AddArray(conc_s_poro);
	fd_poro->AddArray(conc_co2_poro);
    fd_poro->AddArray(trans_poro);
    fd_poro->AddArray(reaction_poro);
    fd_poro->AddArray(vel_poro);
	fd_poro->AddArray(darmkoller);

	// Writing
	auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(getFileName(i, "frac").c_str());
	writer->SetInputData(grid_frac);
	writer->Write();
	writer->SetFileName(getFileName(i, "poro").c_str());
	writer->SetInputData(grid_poro);
	writer->Write();
}
void VTKSnapshotter<acidrecfrac_prod::RecFracProd>::dump_all(int i)
{
	double dy;
    using namespace acidrecfrac_prod;
    const double Y_MULT = 1.0;
    // Grid
    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    // Frac points
    auto points = vtkSmartPointer<vtkPoints>::New();
    for (int k = 0; k < nz - 1; k++)
    {
        for (int j = 0; j < ny - 1; j++)
        {
            const auto& cell = model->cells[j + k * ny];
            const auto& xnebr = model->cells[j + k * ny + ny * nz];
            points->InsertNextPoint(r_dim * (cell.x - xnebr.hx / 30.0), Y_MULT * r_dim * (cell.y + cell.hy / 2.0), r_dim * (cell.z + cell.hz / 2.0));
        }
    }
    for (int i = 0; i < nx - 1; i++)
        for (int k = 0; k < nz - 1; k++)
        { 
            for (int j = 0; j < ny - 1; j++)
            {
                const auto& cell = model->cells[j + k * ny + i * ny * nz];
				const auto& right_nebr_cell = model->cells[j + k * ny + (i + 1) * ny * nz];
				dy = (right_nebr_cell.hx * cell.hy + cell.hx * right_nebr_cell.hy) / (cell.hx + right_nebr_cell.hx) / 2.0;
                points->InsertNextPoint(r_dim * (cell.x + cell.hx / 2.0), Y_MULT * r_dim * (cell.y + dy), r_dim * (cell.z + cell.hz / 2.0));
            }
        }
    grid->SetPoints(points);
    // Data
    auto hexs = vtkSmartPointer<vtkCellArray>::New();
	auto poro = vtkSmartPointer<vtkDoubleArray>::New();
	poro->SetName("porosity");
	auto perm_x = vtkSmartPointer<vtkDoubleArray>::New();
	perm_x->SetName("perm_x");
    auto perm_y = vtkSmartPointer<vtkDoubleArray>::New();
    perm_y->SetName("perm_y");
    auto perm_z = vtkSmartPointer<vtkDoubleArray>::New();
    perm_z->SetName("perm_z");
	auto pres = vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");

    int np = (ny - 1) * (nz - 1);
    std::array<double, 3> vel;
    double buf;
    for (int k = 0; k < nz - 2; k++)
    {
        for (int j = 0; j < ny - 2; j++)
        {
			if (j < 1)
			{
				const auto cell = model->cells[j + 1 + (k + 1) * ny];
				//assert(cell.type == FracType::FRAC_IN);
				const auto& props = model->props_sk[0];
				const auto& next = cell.u_next;
				vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();

				hex->GetPointIds()->SetId(0, j + k * (ny - 1));
				hex->GetPointIds()->SetId(1, j + 1 + k * (ny - 1));
				hex->GetPointIds()->SetId(2, j + 1 + k * (ny - 1) + np);
				hex->GetPointIds()->SetId(3, j + k * (ny - 1) + np);

				hex->GetPointIds()->SetId(4, j + (k + 1) * (ny - 1));
				hex->GetPointIds()->SetId(5, j + 1 + (k + 1) * (ny - 1));
				hex->GetPointIds()->SetId(6, j + 1 + (k + 1) * (ny - 1) + np);
				hex->GetPointIds()->SetId(7, j + (k + 1) * (ny - 1) + np);

				hexs->InsertNextCell(hex);
				poro->InsertNextValue(next.m);
				perm_x->InsertNextValue(M2toMilliDarcy(cell.u_next.kx) * r_dim * r_dim);
                perm_y->InsertNextValue(M2toMilliDarcy(cell.u_next.ky) * r_dim * r_dim);
                perm_z->InsertNextValue(M2toMilliDarcy(cell.u_next.kx) * r_dim * r_dim);
				pres->InsertNextValue(next.p * P_dim / BAR_TO_PA);
			}
        }
    }
    for (int i = 1; i < nx - 1; i++)
    {
        for (int k = 0; k < nz - 2; k++)
        {
            for (int j = 0; j < ny - 2; j++)
            {
                const auto& cell = model->cells[j + 1 + (k + 1) * ny + i * nz * ny];
                //assert(cell.type == FracType::FRAC_MID || cell.type == FracType::FRAC_OUT);
                const auto& next = cell.u_next;
                const auto& props = model->props_sk[0];
                vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();

                hex->GetPointIds()->SetId(0, j + k * (ny - 1) + i * np);
                hex->GetPointIds()->SetId(1, j + 1 + k * (ny - 1) + i * np);
                hex->GetPointIds()->SetId(2, j + 1 + k * (ny - 1) + (i + 1) * np);
                hex->GetPointIds()->SetId(3, j + k * (ny - 1) + (i + 1) * np);

                hex->GetPointIds()->SetId(4, j + (k + 1) * (ny - 1) + i * np);
                hex->GetPointIds()->SetId(5, j + 1 + (k + 1) * (ny - 1) + i * np);
                hex->GetPointIds()->SetId(6, j + 1 + (k + 1) * (ny - 1) + (i + 1) * np);
                hex->GetPointIds()->SetId(7, j + (k + 1) * (ny - 1) + (i + 1) * np);

                hexs->InsertNextCell(hex);
				poro->InsertNextValue(next.m);
				perm_x->InsertNextValue(M2toMilliDarcy(cell.u_next.kx)* r_dim * r_dim);
                perm_y->InsertNextValue(M2toMilliDarcy(cell.u_next.ky)* r_dim * r_dim);
                perm_z->InsertNextValue(M2toMilliDarcy(cell.u_next.kx)* r_dim * r_dim);
				pres->InsertNextValue(next.p * P_dim / BAR_TO_PA);
            }
        }
    }
    
    grid->SetCells(VTK_HEXAHEDRON, hexs);
    vtkCellData* fd = grid->GetCellData();
	fd->AddArray(poro);
	fd->AddArray(perm_x);
    fd->AddArray(perm_y);
    fd->AddArray(perm_z);
	fd->AddArray(pres);
    // Writing
    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(getFileName(i).c_str());
    writer->SetInputData(grid);
    writer->Write();
}
void VTKSnapshotter<acidrecfracmov::AcidRecFracMov>::dump_all(int i)
{
	using namespace acidrecfracmov;
	const double w2 = model->props_frac.w2_init;
	const double Y_MULT = 1.0;
	// Grid
	auto grid_frac = vtkSmartPointer<vtkUnstructuredGrid>::New();
	auto grid_poro = vtkSmartPointer<vtkUnstructuredGrid>::New();
	// Frac points
	auto points_frac = vtkSmartPointer<vtkPoints>::New();
	for (int k = 0; k < nz - 1; k++)
	{
		for (int j = 0; j < ny; j++)
		{
			const FracCell& cell = model->cells_frac[j + k * ny];
			const FracCell& xnebr = model->cells_frac[j + k * ny + ny * nz];
			points_frac->InsertNextPoint(r_dim * (cell.x - xnebr.hx / 30.0), Y_MULT * r_dim * (cell.y + cell.hy / 2.0), r_dim * (cell.z + cell.hz / 2.0));
		}
	}
	for (int i = 0; i < nx - 1; i++)
		for (int k = 0; k < nz - 1; k++)
		{
			for (int j = 0; j < ny; j++)
			{
				const FracCell& cell = model->cells_frac[j + k * ny + i * ny * nz];
				points_frac->InsertNextPoint(r_dim * (cell.x + cell.hx / 2.0), Y_MULT * r_dim * (cell.y + cell.hy / 2.0), r_dim * (cell.z + cell.hz / 2.0));
			}
		}
	grid_frac->SetPoints(points_frac);
	// Poro points
	auto points_poro = vtkSmartPointer<vtkPoints>::New();
	const int ny_poro = model->cellsNum_y_poro + 2;
	for (int i = 0; i < nx - 1; i++)
		for (int k = 0; k < nz - 1; k++)
		{
			for (int j = 0; j < ny_poro; j++)
			{
				const PoroCell& cell = model->cells_poro[j + k * ny_poro + i * ny_poro * nz];
				if (j != 1)
					points_poro->InsertNextPoint(r_dim * (cell.x + cell.hx / 2.0), Y_MULT * r_dim * (cell.y - cell.hy / 2.0), r_dim * (cell.z + cell.hz / 2.0));
				else
					points_poro->InsertNextPoint(r_dim * (cell.x + cell.hx / 2.0), Y_MULT * r_dim * (cell.y - 2.0 * cell.hy / 200.0), r_dim * (cell.z + cell.hz / 2.0));
			}
		}
	grid_poro->SetPoints(points_poro);
	// Data
	auto hexs_frac = vtkSmartPointer<vtkCellArray>::New();
	auto poro_frac = vtkSmartPointer<vtkDoubleArray>::New();
	poro_frac->SetName("porosity");
	auto perm_frac = vtkSmartPointer<vtkDoubleArray>::New();
	perm_frac->SetName("permeability");
	auto pres_frac = vtkSmartPointer<vtkDoubleArray>::New();
	pres_frac->SetName("pressure");
	auto sat_w_frac = vtkSmartPointer<vtkDoubleArray>::New();
	sat_w_frac->SetName("WaterSaturation");
	auto sat_o_frac = vtkSmartPointer<vtkDoubleArray>::New();
	sat_o_frac->SetName("OilSaturation");
	auto conc_a_frac = vtkSmartPointer<vtkDoubleArray>::New();
	conc_a_frac->SetName("AcidConcentration");
	auto conc_w_frac = vtkSmartPointer<vtkDoubleArray>::New();
	conc_w_frac->SetName("WaterConcentration");
	auto conc_co2_frac = vtkSmartPointer<vtkDoubleArray>::New();
	conc_co2_frac->SetName("CO2Concentration");
	auto conc_s_frac = vtkSmartPointer<vtkDoubleArray>::New();
	conc_s_frac->SetName("SaltConcentration");
	auto trans = vtkSmartPointer<vtkDoubleArray>::New();
	trans->SetName("Transmissibility");
	auto v_leak = vtkSmartPointer<vtkDoubleArray>::New();
	v_leak->SetName("Leakoff_Velocity");
	auto width = vtkSmartPointer<vtkDoubleArray>::New();
	width->SetName("Width");
	auto type = vtkSmartPointer<vtkDoubleArray>::New();
	type->SetName("Type");
	//auto reaction_frac = vtkSmartPointer<vtkDoubleArray>::New();
	// reaction_frac->SetName("ReactionSpeed");

	int np = ny * (nz - 1);
	std::array<double, 3> vel;
	double buf;
	for (int k = 0; k < nz - 2; k++)
	{
		for (int j = 0; j < ny - 1; j++)
		{
			const FracCell& cell = model->cells_frac[j + 1 + (k + 1) * ny];
			assert(cell.type == FracType::FRAC_IN);
			const auto& props = model->props_sk[0];
			const auto& next = cell.u_next;
			vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();
			poro_frac->InsertNextValue(props.m_init);
			perm_frac->InsertNextValue(M2toMilliDarcy(props.perm * r_dim * r_dim));
			pres_frac->InsertNextValue(next.p * P_dim / BAR_TO_PA);
			sat_w_frac->InsertNextValue(1.0);
			sat_o_frac->InsertNextValue(0.0);
			conc_a_frac->InsertNextValue(next.c);
			conc_w_frac->InsertNextValue(1.0 - next.c);
			conc_s_frac->InsertNextValue(0.0);
			conc_co2_frac->InsertNextValue(0.0);
			trans->InsertNextValue(model->trans[k] * model->widths[k] * model->R_dim *
				M2toMilliDarcy(props.perm) * model->R_dim * model->R_dim);
			width->InsertNextValue(model->widths[k] * r_dim);
			v_leak->InsertNextValue(0.0);
			type->InsertNextValue(cell.type);
			//reaction_frac->InsertNextValue(0.0);

			hex->GetPointIds()->SetId(0, j + k * ny);
			hex->GetPointIds()->SetId(1, j + 1 + k * ny);
			hex->GetPointIds()->SetId(2, j + 1 + k * ny + np);
			hex->GetPointIds()->SetId(3, j + k * ny + np);

			hex->GetPointIds()->SetId(4, j + (k + 1) * ny);
			hex->GetPointIds()->SetId(5, j + 1 + (k + 1) * ny);
			hex->GetPointIds()->SetId(6, j + 1 + (k + 1) * ny + np);
			hex->GetPointIds()->SetId(7, j + (k + 1) * ny + np);

			hexs_frac->InsertNextCell(hex);
		}
	}
	for (int i = 1; i < nx - 1; i++)
	{
		for (int k = 0; k < nz - 2; k++)
		{
			for (int j = 0; j < ny - 1; j++)
			{
				const FracCell& cell = model->cells_frac[j + 1 + (k + 1) * ny + i * nz * ny];
				assert(cell.type == FracType::FRAC_MID || cell.type == FracType::FRAC_OUT);
				const auto& next = cell.u_next;
				const auto& props = model->props_sk[0];
				vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();
				poro_frac->InsertNextValue(props.m_init);
				perm_frac->InsertNextValue(M2toMilliDarcy(props.perm * r_dim * r_dim));
				pres_frac->InsertNextValue(next.p * P_dim / BAR_TO_PA);
				sat_w_frac->InsertNextValue(1.0);
				sat_o_frac->InsertNextValue(0.0);
				conc_a_frac->InsertNextValue(next.c);
				conc_w_frac->InsertNextValue(1.0 - next.c);
				conc_s_frac->InsertNextValue(0.0);
				conc_co2_frac->InsertNextValue(0.0);
				trans->InsertNextValue(model->trans[k + (i - 1) * model->cellsNum_z] * model->widths[k + (i - 1) * model->cellsNum_z] * model->R_dim *
					M2toMilliDarcy(props.perm) * model->R_dim * model->R_dim);
				width->InsertNextValue(model->widths[k + (i - 1) * model->cellsNum_z] * r_dim);
				v_leak->InsertNextValue(model->getFlowLeak(cell).value() * r_dim / t_dim);
				type->InsertNextValue(cell.type);

				hex->GetPointIds()->SetId(0, j + k * ny + i * np);
				hex->GetPointIds()->SetId(1, j + 1 + k * ny + i * np);
				hex->GetPointIds()->SetId(2, j + 1 + k * ny + (i + 1) * np);
				hex->GetPointIds()->SetId(3, j + k * ny + (i + 1) * np);

				hex->GetPointIds()->SetId(4, j + (k + 1) * ny + i * np);
				hex->GetPointIds()->SetId(5, j + 1 + (k + 1) * ny + i * np);
				hex->GetPointIds()->SetId(6, j + 1 + (k + 1) * ny + (i + 1) * np);
				hex->GetPointIds()->SetId(7, j + (k + 1) * ny + (i + 1) * np);
				hexs_frac->InsertNextCell(hex);
			}
		}
	}

	// Data
	auto hexs_poro = vtkSmartPointer<vtkCellArray>::New();
	auto poro_poro = vtkSmartPointer<vtkDoubleArray>::New();
	poro_poro->SetName("porosity");
	auto perm_poro = vtkSmartPointer<vtkDoubleArray>::New();
	perm_poro->SetName("permeability");
	auto pres_poro = vtkSmartPointer<vtkDoubleArray>::New();
	pres_poro->SetName("pressure");
	auto sat_w_poro = vtkSmartPointer<vtkDoubleArray>::New();
	sat_w_poro->SetName("WaterSaturation");
	auto sat_o_poro = vtkSmartPointer<vtkDoubleArray>::New();
	sat_o_poro->SetName("OilSaturation");
	auto conc_a_poro = vtkSmartPointer<vtkDoubleArray>::New();
	conc_a_poro->SetName("AcidConcentration");
	auto conc_w_poro = vtkSmartPointer<vtkDoubleArray>::New();
	conc_w_poro->SetName("WaterConcentration");
	auto conc_co2_poro = vtkSmartPointer<vtkDoubleArray>::New();
	conc_co2_poro->SetName("CO2Concentration");
	auto conc_s_poro = vtkSmartPointer<vtkDoubleArray>::New();
	conc_s_poro->SetName("SaltConcentration");
	auto trans_poro = vtkSmartPointer<vtkDoubleArray>::New();
	trans_poro->SetName("Transmissibility");
	auto reaction_poro = vtkSmartPointer<vtkDoubleArray>::New();
	reaction_poro->SetName("ReactionSpeed");
	auto vel_poro = vtkSmartPointer<vtkDoubleArray>::New();
	vel_poro->SetName("Velocity");
	vel_poro->SetNumberOfComponents(3);
	auto darmkoller = vtkSmartPointer<vtkDoubleArray>::New();
	darmkoller->SetName("Darmkoller");

	double sum_width, sum_trans;
	int np_poro = ny_poro * (nz - 1);
	for (int i = 0; i < nx - 2; i++)
	{
		for (int k = 0; k < nz - 2; k++)
		{
			sum_width = sum_trans = 0.0;
			for (int j = 0; j < ny_poro - 1; j++)
			{
				const PoroCell& cell = model->cells_poro[j + (k + 1) * ny_poro + (i + 1) * nz * ny_poro];
				if (cell.y * r_dim < 0.1)
				{
					assert(cell.type == PoroType::MIDDLE || (cell.type == PoroType::WELL_LAT && j == 0));
					const auto& next = cell.u_next;
					vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();

					poro_poro->InsertNextValue(next.m);
					perm_poro->InsertNextValue(M2toMilliDarcy(cell.props->getPermCoseni(next.m, next.p).value() * r_dim * r_dim));
					pres_poro->InsertNextValue(next.p * P_dim / BAR_TO_PA);
					sat_w_poro->InsertNextValue(next.sw);
					sat_o_poro->InsertNextValue(1.0 - next.sw);
					conc_a_poro->InsertNextValue(next.xa);
					conc_w_poro->InsertNextValue(next.xw);
					conc_s_poro->InsertNextValue(next.xs);
					conc_co2_poro->InsertNextValue(1.0 - next.xw - next.xa - next.xs);

					sum_width += cell.hy * r_dim;
					sum_trans += M2toMilliDarcy(cell.props->getPermCoseni(next.m, next.p).value() * r_dim * r_dim) * cell.hy * r_dim;
					trans_poro->InsertNextValue(sum_trans);

					if (cell.type != PoroType::WELL_LAT)
					{
						buf = model->getReactionRateOutput(next, *cell.props).value();
						reaction_poro->InsertNextValue(buf / r_dim / r_dim / r_dim / t_dim);
						vel = model->getPoroWaterVelocity(cell);
						vel[0] *= r_dim / t_dim;    vel[1] *= r_dim / t_dim;    vel[2] *= r_dim / t_dim;
						vel_poro->InsertNextTuple(&vel[0]);
						darmkoller->InsertNextValue(model->getDarmkoller(cell, next, *cell.props));
					}
					else
					{
						buf = 0.0;
						reaction_poro->InsertNextValue(buf / r_dim / r_dim / r_dim / t_dim);
						vel[0] = 0.0;   vel[2] = 0.0;
						vel[1] = model->getFlowLeak(model->cells_frac[model->getFracNebr(cell.num)]).value() * r_dim / t_dim;
						vel_poro->InsertNextTuple(&vel[0]);
						darmkoller->InsertNextValue(0.0);
					}

					hex->GetPointIds()->SetId(0, j + k * ny_poro + i * np_poro);
					hex->GetPointIds()->SetId(1, j + 1 + k * ny_poro + i * np_poro);
					hex->GetPointIds()->SetId(2, j + 1 + k * ny_poro + (i + 1) * np_poro);
					hex->GetPointIds()->SetId(3, j + k * ny_poro + (i + 1) * np_poro);

					hex->GetPointIds()->SetId(4, j + (k + 1) * ny_poro + i * np_poro);
					hex->GetPointIds()->SetId(5, j + 1 + (k + 1) * ny_poro + i * np_poro);
					hex->GetPointIds()->SetId(6, j + 1 + (k + 1) * ny_poro + (i + 1) * np_poro);
					hex->GetPointIds()->SetId(7, j + (k + 1) * ny_poro + (i + 1) * np_poro);
					hexs_poro->InsertNextCell(hex);
				}
			}
		}
	}

	grid_frac->SetCells(VTK_HEXAHEDRON, hexs_frac);
	vtkCellData* fd_frac = grid_frac->GetCellData();
	fd_frac->AddArray(poro_frac);
	fd_frac->AddArray(perm_frac);
	fd_frac->AddArray(pres_frac);
	fd_frac->AddArray(sat_w_frac);
	fd_frac->AddArray(sat_o_frac);
	fd_frac->AddArray(conc_a_frac);
	fd_frac->AddArray(conc_w_frac);
	fd_frac->AddArray(conc_s_frac);
	fd_frac->AddArray(conc_co2_frac);
	fd_frac->AddArray(trans);
	fd_frac->AddArray(width);
	fd_frac->AddArray(v_leak);
	fd_frac->AddArray(type);
	//fd_frac->AddArray(reaction_frac);

	grid_poro->SetCells(VTK_HEXAHEDRON, hexs_poro);
	vtkCellData* fd_poro = grid_poro->GetCellData();
	fd_poro->AddArray(poro_poro);
	fd_poro->AddArray(perm_poro);
	fd_poro->AddArray(pres_poro);
	fd_poro->AddArray(sat_w_poro);
	fd_poro->AddArray(sat_o_poro);
	fd_poro->AddArray(conc_a_poro);
	fd_poro->AddArray(conc_w_poro);
	fd_poro->AddArray(conc_s_poro);
	fd_poro->AddArray(conc_co2_poro);
	fd_poro->AddArray(trans_poro);
	fd_poro->AddArray(reaction_poro);
	fd_poro->AddArray(vel_poro);
	fd_poro->AddArray(darmkoller);

	// Writing
	auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(getFileName(i, "frac").c_str());
	writer->SetInputData(grid_frac);
	writer->Write();
	writer->SetFileName(getFileName(i, "poro").c_str());
	writer->SetInputData(grid_poro);
	writer->Write();
}
void VTKSnapshotter<vpp2d::VPP2d>::dump_all(int i)
{
	using namespace vpp2d;

	// Grid
	vtkSmartPointer<vtkPolyData> grid =
		vtkSmartPointer<vtkPolyData>::New();

	// Points
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();
	for (int j = 1; j < ny; j++)
	{
		Cell& cell = model->cells[j];
		points->InsertNextPoint(r_dim * (0.9 * cell.r), -r_dim * (cell.z - cell.hz / 2.0), 0.0);
	}
	for (int k = 1; k < nx; k++)
	{
		for (int j = 1; j < ny; j++)
		{
			Cell& cell = model->cells[k * ny + j];
			points->InsertNextPoint(r_dim * (cell.r - cell.hr / 2.0), -r_dim * (cell.z - cell.hz / 2.0), 0.0);
		}
	}
	grid->SetPoints(points);

	// Data
	vtkSmartPointer<vtkCellArray> polygons =
		vtkSmartPointer<vtkCellArray>::New();

	vtkSmartPointer<vtkPolygon> polygon =
		vtkSmartPointer<vtkPolygon>::New();
	polygon->GetPointIds()->SetNumberOfIds(4);

	vtkSmartPointer<vtkDoubleArray> pres =
		vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");

	vtkSmartPointer<vtkDoubleArray> sat_w =
		vtkSmartPointer<vtkDoubleArray>::New();
	sat_w->SetName("WaterSaturation");

	vtkSmartPointer<vtkDoubleArray> sat_o =
		vtkSmartPointer<vtkDoubleArray>::New();
	sat_o->SetName("OilSaturation");

	vtkSmartPointer<vtkDoubleArray> conc_a =
		vtkSmartPointer<vtkDoubleArray>::New();
	conc_a->SetName("ActiveConcentration");

	vtkSmartPointer<vtkDoubleArray> vel_w =
		vtkSmartPointer<vtkDoubleArray>::New();
	vel_w->SetName("WaterVelocity");
	vel_w->SetNumberOfComponents(3);

	vtkSmartPointer<vtkDoubleArray> vel_o =
		vtkSmartPointer<vtkDoubleArray>::New();
	vel_o->SetName("OilVelocity");
	vel_o->SetNumberOfComponents(3);

	int k, j, idx, idx1;

	double vel[3];

	for (j = 0; j < ny - 2; j++)
	{
		idx = j;
		idx1 = j + 1;
		Cell& cell = model->cells[idx1];

		polygon->GetPointIds()->SetId(0, idx);
		polygon->GetPointIds()->SetId(1, idx + ny - 1);
		polygon->GetPointIds()->SetId(2, idx + ny);
		polygon->GetPointIds()->SetId(3, idx + 1);
		polygons->InsertNextCell(polygon);

		pres->InsertNextValue(cell.u_next.p * P_dim);
		sat_w->InsertNextValue(cell.u_next.s);
		sat_o->InsertNextValue(1.0 - cell.u_next.s);
		conc_a->InsertNextValue(cell.u_next.c);

		vel[0] = 0.0;
		vel[1] = 0.0;
		vel[2] = 0.0;
		vel_w->InsertNextTuple(vel);
		vel[0] = 0.0;
		vel[1] = 0.0;
		vel[2] = 0.0;
		vel_o->InsertNextTuple(vel);
	}

	// Middle cells
	for (k = 1; k < nx - 1; k++)
	{
		for (j = 0; j < ny - 2; j++)
		{
			idx = k * (ny - 1) + j;
			idx1 = k * ny + j + 1;
			Cell& cell = model->cells[idx1];

			polygon->GetPointIds()->SetId(0, idx);
			polygon->GetPointIds()->SetId(1, idx + ny - 1);
			polygon->GetPointIds()->SetId(2, idx + ny);
			polygon->GetPointIds()->SetId(3, idx + 1);
			polygons->InsertNextCell(polygon);

			pres->InsertNextValue(cell.u_next.p * P_dim);
			sat_w->InsertNextValue(cell.u_next.s);
			sat_o->InsertNextValue(1.0 - cell.u_next.s);
			conc_a->InsertNextValue(cell.u_next.c);
			vel[0] = r_dim / t_dim * model->getWaterVelocity(cell, NEXT, R_AXIS);
			vel[1] = r_dim / t_dim * model->getWaterVelocity(cell, NEXT, Z_AXIS);
			vel[2] = 0.0;
			vel_w->InsertNextTuple(vel);
			vel[0] = r_dim / t_dim * model->getOilVelocity(cell, NEXT, R_AXIS);
			vel[1] = r_dim / t_dim * model->getOilVelocity(cell, NEXT, Z_AXIS);
			vel[2] = 0.0;
			vel_o->InsertNextTuple(vel);
		}
	}

	grid->SetPolys(polygons);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(pres);
	fd->AddArray(sat_w);
	fd->AddArray(sat_o);
	fd->AddArray(conc_a);
	fd->AddArray(vel_w);
	fd->AddArray(vel_o);

	// Writing
	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();

	writer->SetFileName(getFileName(i).c_str());
	writer->SetInputData(grid);
	writer->Write();
}
void VTKSnapshotter<bing1d::Bingham1d>::dump_all(int i)
{
	using namespace bing1d;

	// Grid
	auto grid = vtkSmartPointer<vtkPolyData>::New();

	// Points
	auto points = vtkSmartPointer<vtkPoints>::New();

	Cell& cell = model->cells[0];
	points->InsertNextPoint(r_dim * (0.9 * cell.r), 0.0, 0.0);
	points->InsertNextPoint(r_dim * (0.9 * cell.r), -r_dim * model->props_sk.height, 0.0);

	for (int k = 1; k < nx; k++)
	{
		Cell& cell = model->cells[k];
		points->InsertNextPoint(r_dim * (cell.r - cell.hr / 2.0), 0.0, 0.0);
		points->InsertNextPoint(r_dim * (cell.r - cell.hr / 2.0), -r_dim * model->props_sk.height, 0.0);
	}
	grid->SetPoints(points);

	// Data
	auto polygons = vtkSmartPointer<vtkCellArray>::New();
	auto polygon = vtkSmartPointer<vtkPolygon>::New();
	polygon->GetPointIds()->SetNumberOfIds(4);

	auto pres = vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");
	auto vel_oil = vtkSmartPointer<vtkDoubleArray>::New();
	vel_oil->SetName("oilVelocity");

	int k, j, idx, idx1;

	cell = model->cells[0];
	polygon->GetPointIds()->SetId(0, 0);
	polygon->GetPointIds()->SetId(1, 2);
	polygon->GetPointIds()->SetId(2, 3);
	polygon->GetPointIds()->SetId(3, 1);
	polygons->InsertNextCell(polygon);
	pres->InsertNextValue(cell.u_next.p * P_dim / BAR_TO_PA);
	vel_oil->InsertNextValue(0.0);


	// Middle cells
	for (k = 1; k < nx - 1; k++)
	{
		Cell& cell = model->cells[k];

		polygon->GetPointIds()->SetId(0, 2 * k);
		polygon->GetPointIds()->SetId(1, 2 * k + 2);
		polygon->GetPointIds()->SetId(2, 2 * k + 3);
		polygon->GetPointIds()->SetId(3, 2 * k + 1);
		polygons->InsertNextCell(polygon);

		pres->InsertNextValue(cell.u_next.p * P_dim / BAR_TO_PA);
		vel_oil->InsertNextValue( r_dim / t_dim * model->getOilVelocity(cell) );
	}

	grid->SetPolys(polygons);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(pres);
	fd->AddArray(vel_oil);

	// Writing
	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

	writer->SetFileName(getFileName(i).c_str());
	writer->SetInputData(grid);
	writer->Write();
}
void VTKSnapshotter<gasOil_elliptic::GasOil_Elliptic>::dump_all(int snap_idx)
{
	using namespace gasOil_elliptic;

	// Grid
	auto grid =	vtkSmartPointer<vtkUnstructuredGrid>::New();

	// Points
	auto points = vtkSmartPointer<vtkPoints>::New();

	for (int m = 0; m < ny; m++)
	{
		for (int j = 1; j < nz; j++)
		{
			Cell& cell = model->cells[j + m * nx * nz];
			Point point = getCartesian<Cell>(0, cell.nu - cell.hnu / 2.0, cell.z - cell.hz / 2.0);
			points->InsertNextPoint(r_dim * point[0] * VIEW_MULTIPLIER,
									r_dim * point[1] * VIEW_MULTIPLIER,
									-r_dim * point[2]);
		}

		for (int k = 1; k < nx; k++)
		{
			for (int j = 1; j < nz; j++)
			{
				Cell& cell = model->cells[j + k * nz + m * nx * nz];
				Point point = getCartesian<Cell>(cell.mu - cell.hmu / 2.0, cell.nu - cell.hnu / 2.0, cell.z - cell.hz / 2.0);
				points->InsertNextPoint(r_dim * point[0] * VIEW_MULTIPLIER,
										r_dim * point[1] * VIEW_MULTIPLIER,
										-r_dim * point[2]);
			}
		}
	}
	grid->SetPoints(points);

	// Data
	auto hexs = vtkSmartPointer<vtkCellArray>::New();
	auto pres = vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");
	auto p_bub = vtkSmartPointer<vtkDoubleArray>::New();
	p_bub->SetName("buble_point");
	auto sat_oil = vtkSmartPointer<vtkDoubleArray>::New();
	sat_oil->SetName("oilSaturation");
	auto sat_gas = vtkSmartPointer<vtkDoubleArray>::New();
	sat_gas->SetName("gasSaturation");
	auto satur = vtkSmartPointer<vtkIntArray>::New();
	satur->SetName("SATUR");
	auto vel_oil = vtkSmartPointer<vtkDoubleArray>::New();
	vel_oil->SetName("oilVelocity");
	vel_oil->SetNumberOfComponents(3);
	auto vel_gas = vtkSmartPointer<vtkDoubleArray>::New();
	vel_gas->SetName("gasVelocity");
	vel_gas->SetNumberOfComponents(3);

	double vel[3];
	int l, j, k;

	for (l = 0; l < ny - 1; l++)
	{
		// Left
		for (j = 0; j < nz - 2; j++)
		{
			const Cell& cell = model->cells[l*nx*nz + j + 1];

			if (cell.isUsed)
			{
				vtkSmartPointer<vtkHexahedron> hex =
					vtkSmartPointer<vtkHexahedron>::New();
				hex->GetPointIds()->SetId(0, j + l * nx * (nz - 1));
				hex->GetPointIds()->SetId(1, j + (l + 1) * nx * (nz - 1));
				hex->GetPointIds()->SetId(2, j + (l + 1) * nx * (nz - 1) + 1);
				hex->GetPointIds()->SetId(3, j + l * nx * (nz - 1) + 1);

				hex->GetPointIds()->SetId(4, j + l * nx * (nz - 1) + (nz - 1));
				hex->GetPointIds()->SetId(5, j + (l + 1) * nx * (nz - 1) + (nz - 1));
				hex->GetPointIds()->SetId(6, j + (l + 1) * nx * (nz - 1) + 1 + (nz - 1));
				hex->GetPointIds()->SetId(7, j + l * nx * (nz - 1) + 1 + (nz - 1));

				hexs->InsertNextCell(hex);

				pres->InsertNextValue(cell.u_next.p * P_dim);
				p_bub->InsertNextValue(cell.u_next.p_bub * P_dim);
				sat_oil->InsertNextValue(cell.u_next.s);
				sat_gas->InsertNextValue(1.0 - cell.u_next.s);
				satur->InsertNextValue(cell.u_next.SATUR);
				vel[0] = 0.0;
				vel[1] = 0.0;
				vel[2] = 0.0;
				vel_oil->InsertNextTuple(vel);
				vel[0] = 0.0;
				vel[1] = 0.0;
				vel[2] = 0.0;
				vel_gas->InsertNextTuple(vel);
			}
		}

		// Middle cells
		for (k = 1; k < nx - 1; k++)
		{
			for (j = 0; j < nz - 2; j++)
			{
				const Cell& cell = model->cells[l*nx*nz + k*nz + j + 1];

				vtkSmartPointer<vtkHexahedron> hex =
					vtkSmartPointer<vtkHexahedron>::New();
				hex->GetPointIds()->SetId(0, j + k*(nz - 1) + l * nx * (nz - 1));
				hex->GetPointIds()->SetId(1, j + k*(nz - 1) + (l + 1) * nx * (nz - 1));
				hex->GetPointIds()->SetId(2, j + k*(nz - 1) + (l + 1) * nx * (nz - 1) + 1);
				hex->GetPointIds()->SetId(3, j + k*(nz - 1) + l * nx * (nz - 1) + 1);

				hex->GetPointIds()->SetId(4, j + (k + 1)*(nz - 1) + l * nx * (nz - 1));
				hex->GetPointIds()->SetId(5, j + (k + 1)*(nz - 1) + (l + 1) * nx * (nz - 1));
				hex->GetPointIds()->SetId(6, j + (k + 1)*(nz - 1) + (l + 1) * nx * (nz - 1) + 1);
				hex->GetPointIds()->SetId(7, j + (k + 1)*(nz - 1) + l * nx * (nz - 1) + 1);

				hexs->InsertNextCell(hex);

				pres->InsertNextValue(cell.u_next.p * P_dim);
				p_bub->InsertNextValue(cell.u_next.p_bub * P_dim);
				sat_oil->InsertNextValue(cell.u_next.s);
				sat_gas->InsertNextValue(1.0 - cell.u_next.s);
				satur->InsertNextValue(cell.u_next.SATUR);
				vel[0] = 0.0;// r_dim / t_dim * (cos(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) - sin(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
				vel[1] = 0.0;// r_dim / t_dim * (sin(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) + cos(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
				vel[2] = 0.0;// r_dim / t_dim * model->getOilVelocity(cell, NEXT, Z_AXIS);
				vel_oil->InsertNextTuple(vel);
				vel[0] = 0.0;// r_dim / t_dim * (cos(cell.phi) * model->getGasVelocity(cell, NEXT, R_AXIS) - sin(cell.phi) * model->getGasVelocity(cell, NEXT, PHI_AXIS));
				vel[1] = 0.0;// r_dim / t_dim * (sin(cell.phi) * model->getGasVelocity(cell, NEXT, R_AXIS) + cos(cell.phi) * model->getGasVelocity(cell, NEXT, PHI_AXIS));
				vel[2] = 0.0;// r_dim / t_dim * model->getGasVelocity(cell, NEXT, Z_AXIS);
				vel_gas->InsertNextTuple(vel);
			}
		}
	}

	// Last segment - left
	for (j = 0; j < nz - 2; j++)
	{
		const Cell& cell = model->cells[l*nx*nz + j + 1];

		if (cell.isUsed)
		{
			vtkSmartPointer<vtkHexahedron> hex =
				vtkSmartPointer<vtkHexahedron>::New();
			hex->GetPointIds()->SetId(0, j + l * nx * (nz - 1));
			hex->GetPointIds()->SetId(1, j);
			hex->GetPointIds()->SetId(2, j + 1);
			hex->GetPointIds()->SetId(3, j + l * nx * (nz - 1) + 1);

			hex->GetPointIds()->SetId(4, j + l * nx * (nz - 1) + (nz - 1));
			hex->GetPointIds()->SetId(5, j + (nz - 1));
			hex->GetPointIds()->SetId(6, j + 1 + (nz - 1));
			hex->GetPointIds()->SetId(7, j + l * nx * (nz - 1) + 1 + (nz - 1));

			hexs->InsertNextCell(hex);

			pres->InsertNextValue(cell.u_next.p * P_dim);
			p_bub->InsertNextValue(cell.u_next.p_bub * P_dim);
			sat_oil->InsertNextValue(cell.u_next.s);
			sat_gas->InsertNextValue(1.0 - cell.u_next.s);
			satur->InsertNextValue(cell.u_next.SATUR);
			vel[0] = 0.0;
			vel[1] = 0.0;
			vel[2] = 0.0;
			vel_oil->InsertNextTuple(vel);
			vel[0] = 0.0;
			vel[1] = 0.0;
			vel[2] = 0.0;
			vel_gas->InsertNextTuple(vel);
		}
	}

	// Last segment - middle
	for (k = 1; k < nx - 1; k++)
	{
		for (j = 0; j < nz - 2; j++)
		{
			const Cell& cell = model->cells[l*nx*nz + k*nz + j + 1];

			vtkSmartPointer<vtkHexahedron> hex =
				vtkSmartPointer<vtkHexahedron>::New();
			hex->GetPointIds()->SetId(0, j + k*(nz - 1) + l * nx * (nz - 1));
			hex->GetPointIds()->SetId(1, j + k*(nz - 1));
			hex->GetPointIds()->SetId(2, j + k*(nz - 1) + 1);
			hex->GetPointIds()->SetId(3, j + k*(nz - 1) + l * nx * (nz - 1) + 1);

			hex->GetPointIds()->SetId(4, j + (k + 1)*(nz - 1) + l * nx * (nz - 1));
			hex->GetPointIds()->SetId(5, j + (k + 1)*(nz - 1));
			hex->GetPointIds()->SetId(6, j + (k + 1)*(nz - 1) + 1);
			hex->GetPointIds()->SetId(7, j + (k + 1)*(nz - 1) + l * nx * (nz - 1) + 1);

			hexs->InsertNextCell(hex);

			pres->InsertNextValue(cell.u_next.p * P_dim);
			p_bub->InsertNextValue(cell.u_next.p_bub * P_dim);
			sat_oil->InsertNextValue(cell.u_next.s);
			sat_gas->InsertNextValue(1.0 - cell.u_next.s);
			satur->InsertNextValue(cell.u_next.SATUR);
			vel[0] = 0.0;// r_dim / t_dim * (cos(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) - sin(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
			vel[1] = 0.0;// r_dim / t_dim * (sin(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) + cos(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
			vel[2] = 0.0;// r_dim / t_dim * model->getOilVelocity(cell, NEXT, Z_AXIS);
			vel_oil->InsertNextTuple(vel);
			vel[0] = 0.0;// r_dim / t_dim * (cos(cell.phi) * model->getGasVelocity(cell, NEXT, R_AXIS) - sin(cell.phi) * model->getGasVelocity(cell, NEXT, PHI_AXIS));
			vel[1] = 0.0;// r_dim / t_dim * (sin(cell.phi) * model->getGasVelocity(cell, NEXT, R_AXIS) + cos(cell.phi) * model->getGasVelocity(cell, NEXT, PHI_AXIS));
			vel[2] = 0.0;// r_dim / t_dim * model->getGasVelocity(cell, NEXT, Z_AXIS);
			vel_gas->InsertNextTuple(vel);
		}
	}

	grid->SetCells(VTK_HEXAHEDRON, hexs);
	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(pres);
	fd->AddArray(p_bub);
	fd->AddArray(satur);
	fd->AddArray(sat_oil);
	fd->AddArray(sat_gas);
	fd->AddArray(vel_oil);
	fd->AddArray(vel_gas);

	// Writing
	auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(getFileName(snap_idx).c_str());
	writer->SetInputData(grid);
	writer->Write();
}
void VTKSnapshotter<oilnit_elliptic::OilNIT_Elliptic>::dump_all(int snap_idx)
{
	using namespace oilnit_elliptic;

	// Grid
	auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	// Points
	auto points = vtkSmartPointer<vtkPoints>::New();

	for (int m = 0; m < ny; m++)
	{
		for (int j = 1; j < nz; j++)
		{
			Cell& cell = model->cells[j + m * nx * nz];
			Point point = getCartesian<Cell>(0, cell.nu - cell.hnu / 2.0, cell.z - cell.hz / 2.0);
			points->InsertNextPoint(r_dim * point[0] * VIEW_MULTIPLIER,
				r_dim * point[1] * VIEW_MULTIPLIER,
				-r_dim * point[2]);
		}

		for (int k = 1; k < nx; k++)
		{
			for (int j = 1; j < nz; j++)
			{
				Cell& cell = model->cells[j + k * nz + m * nx * nz];
				Point point = getCartesian<Cell>(cell.mu - cell.hmu / 2.0, cell.nu - cell.hnu / 2.0, cell.z - cell.hz / 2.0);
				points->InsertNextPoint(r_dim * point[0] * VIEW_MULTIPLIER,
					r_dim * point[1] * VIEW_MULTIPLIER,
					-r_dim * point[2]);
			}
		}
	}
	grid->SetPoints(points);

	// Data
	auto hexs = vtkSmartPointer<vtkCellArray>::New();
	auto pres = vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");
	auto temp = vtkSmartPointer<vtkDoubleArray>::New();
	temp->SetName("temperature");
	auto perm_xy = vtkSmartPointer<vtkDoubleArray>::New();
	perm_xy->SetName("perm_xy");
	auto vel_oil = vtkSmartPointer<vtkDoubleArray>::New();
	vel_oil->SetName("oilVelocity");
	vel_oil->SetNumberOfComponents(3);

	double vel[3];
	int l, j, k;
	Cell* neighbor[6];

	for (l = 0; l < ny - 1; l++)
	{
		// Left
		for (j = 0; j < nz - 2; j++)
		{
			Cell& cell = model->cells[l*nx*nz + j + 1];

			if (cell.isUsed)
			{
				vtkSmartPointer<vtkHexahedron> hex =
					vtkSmartPointer<vtkHexahedron>::New();
				hex->GetPointIds()->SetId(0, j + l * nx * (nz - 1));
				hex->GetPointIds()->SetId(1, j + (l + 1) * nx * (nz - 1));
				hex->GetPointIds()->SetId(2, j + (l + 1) * nx * (nz - 1) + 1);
				hex->GetPointIds()->SetId(3, j + l * nx * (nz - 1) + 1);

				hex->GetPointIds()->SetId(4, j + l * nx * (nz - 1) + (nz - 1));
				hex->GetPointIds()->SetId(5, j + (l + 1) * nx * (nz - 1) + (nz - 1));
				hex->GetPointIds()->SetId(6, j + (l + 1) * nx * (nz - 1) + 1 + (nz - 1));
				hex->GetPointIds()->SetId(7, j + l * nx * (nz - 1) + 1 + (nz - 1));

				hexs->InsertNextCell(hex);

				temp->InsertNextValue(cell.u_next.t * T_dim);
				pres->InsertNextValue(cell.u_next.p * P_dim);
				perm_xy->InsertNextValue(M2toMilliDarcy(model->getPerm_mu(cell) * r_dim * r_dim));

				model->getNeighbors(cell, neighbor);
				auto vel_cart = cell.getVectorCartesian(model->getVelocity(cell, neighbor, MU_AXIS),
					model->getVelocity(cell, neighbor, NU_AXIS),
					model->getVelocity(cell, neighbor, Z_AXIS));
				vel[0] = vel_cart[0] * r_dim / t_dim;
				vel[1] = vel_cart[1] * r_dim / t_dim;
				vel[2] = -vel_cart[2] * r_dim / t_dim;
				vel_oil->InsertNextTuple(vel);
			}
		}

		// Middle cells
		for (k = 1; k < nx - 1; k++)
		{
			for (j = 0; j < nz - 2; j++)
			{
				Cell& cell = model->cells[l*nx*nz + k*nz + j + 1];

				vtkSmartPointer<vtkHexahedron> hex =
					vtkSmartPointer<vtkHexahedron>::New();
				hex->GetPointIds()->SetId(0, j + k*(nz - 1) + l * nx * (nz - 1));
				hex->GetPointIds()->SetId(1, j + k*(nz - 1) + (l + 1) * nx * (nz - 1));
				hex->GetPointIds()->SetId(2, j + k*(nz - 1) + (l + 1) * nx * (nz - 1) + 1);
				hex->GetPointIds()->SetId(3, j + k*(nz - 1) + l * nx * (nz - 1) + 1);

				hex->GetPointIds()->SetId(4, j + (k + 1)*(nz - 1) + l * nx * (nz - 1));
				hex->GetPointIds()->SetId(5, j + (k + 1)*(nz - 1) + (l + 1) * nx * (nz - 1));
				hex->GetPointIds()->SetId(6, j + (k + 1)*(nz - 1) + (l + 1) * nx * (nz - 1) + 1);
				hex->GetPointIds()->SetId(7, j + (k + 1)*(nz - 1) + l * nx * (nz - 1) + 1);

				hexs->InsertNextCell(hex);

				temp->InsertNextValue(cell.u_next.t * T_dim);
				pres->InsertNextValue(cell.u_next.p * P_dim);
				perm_xy->InsertNextValue(M2toMilliDarcy(model->getPerm_mu(cell) * r_dim * r_dim));

				model->getNeighbors(cell, neighbor);
				auto vel_cart = cell.getVectorCartesian(model->getVelocity(cell, neighbor, MU_AXIS),
					model->getVelocity(cell, neighbor, NU_AXIS),
					model->getVelocity(cell, neighbor, Z_AXIS));
				vel[0] = vel_cart[0] * r_dim / t_dim;
				vel[1] = vel_cart[1] * r_dim / t_dim;
				vel[2] = -vel_cart[2] * r_dim / t_dim;
				vel_oil->InsertNextTuple(vel);
			}
		}
	}

	// Last segment - left
	for (j = 0; j < nz - 2; j++)
	{
		Cell& cell = model->cells[l*nx*nz + j + 1];

		if (cell.isUsed)
		{
			vtkSmartPointer<vtkHexahedron> hex =
				vtkSmartPointer<vtkHexahedron>::New();
			hex->GetPointIds()->SetId(0, j + l * nx * (nz - 1));
			hex->GetPointIds()->SetId(1, j);
			hex->GetPointIds()->SetId(2, j + 1);
			hex->GetPointIds()->SetId(3, j + l * nx * (nz - 1) + 1);

			hex->GetPointIds()->SetId(4, j + l * nx * (nz - 1) + (nz - 1));
			hex->GetPointIds()->SetId(5, j + (nz - 1));
			hex->GetPointIds()->SetId(6, j + 1 + (nz - 1));
			hex->GetPointIds()->SetId(7, j + l * nx * (nz - 1) + 1 + (nz - 1));

			hexs->InsertNextCell(hex);

			temp->InsertNextValue(cell.u_next.t * T_dim);
			pres->InsertNextValue(cell.u_next.p * P_dim);
			perm_xy->InsertNextValue(M2toMilliDarcy(model->getPerm_mu(cell) * r_dim * r_dim));

			model->getNeighbors(cell, neighbor);
			auto vel_cart = cell.getVectorCartesian(model->getVelocity(cell, neighbor, MU_AXIS),
				model->getVelocity(cell, neighbor, NU_AXIS),
				model->getVelocity(cell, neighbor, Z_AXIS));
			vel[0] = vel_cart[0] * r_dim / t_dim;
			vel[1] = vel_cart[1] * r_dim / t_dim;
			vel[2] = -vel_cart[2] * r_dim / t_dim;
			vel_oil->InsertNextTuple(vel);
		}
	}

	// Last segment - middle
	for (k = 1; k < nx - 1; k++)
	{
		for (j = 0; j < nz - 2; j++)
		{
			Cell& cell = model->cells[l*nx*nz + k*nz + j + 1];

			vtkSmartPointer<vtkHexahedron> hex =
				vtkSmartPointer<vtkHexahedron>::New();
			hex->GetPointIds()->SetId(0, j + k*(nz - 1) + l * nx * (nz - 1));
			hex->GetPointIds()->SetId(1, j + k*(nz - 1));
			hex->GetPointIds()->SetId(2, j + k*(nz - 1) + 1);
			hex->GetPointIds()->SetId(3, j + k*(nz - 1) + l * nx * (nz - 1) + 1);

			hex->GetPointIds()->SetId(4, j + (k + 1)*(nz - 1) + l * nx * (nz - 1));
			hex->GetPointIds()->SetId(5, j + (k + 1)*(nz - 1));
			hex->GetPointIds()->SetId(6, j + (k + 1)*(nz - 1) + 1);
			hex->GetPointIds()->SetId(7, j + (k + 1)*(nz - 1) + l * nx * (nz - 1) + 1);

			hexs->InsertNextCell(hex);

			pres->InsertNextValue(cell.u_next.p * P_dim);
			temp->InsertNextValue(cell.u_next.t * T_dim);
			perm_xy->InsertNextValue(M2toMilliDarcy(model->getPerm_mu(cell) * r_dim * r_dim));
			
			model->getNeighbors(cell, neighbor);
			auto vel_cart = cell.getVectorCartesian(model->getVelocity(cell, neighbor, MU_AXIS),
													model->getVelocity(cell, neighbor, NU_AXIS), 
													model->getVelocity(cell, neighbor, Z_AXIS));
			vel[0] = vel_cart[0] * r_dim / t_dim;
			vel[1] = vel_cart[1] * r_dim / t_dim;
			vel[2] = -vel_cart[2] * r_dim / t_dim;
			vel_oil->InsertNextTuple(vel);
		}
	}

	grid->SetCells(VTK_HEXAHEDRON, hexs);
	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(pres);
	fd->AddArray(temp);
	fd->AddArray(perm_xy);
	fd->AddArray(vel_oil);

	// Writing
	auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(getFileName(snap_idx).c_str());
	writer->SetInputData(grid);
	writer->Write();
}
void VTKSnapshotter<blackoilnit_elliptic::BlackOilNIT_Elliptic>::dump_all(int snap_idx)
{
	using namespace blackoilnit_elliptic;

	// Grid
	auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	// Points
	auto points = vtkSmartPointer<vtkPoints>::New();

	for (int m = 0; m < ny; m++)
	{
		for (int j = 1; j < nz; j++)
		{
			Cell& cell = model->cells[j + m * nx * nz];
			Point point = getCartesian<Cell>(0, cell.nu - cell.hnu / 2.0, cell.z - cell.hz / 2.0);
			points->InsertNextPoint(r_dim * point[0] * VIEW_MULTIPLIER,
				r_dim * point[1] * VIEW_MULTIPLIER,
				-r_dim * point[2]);
		}

		for (int k = 1; k < nx; k++)
		{
			for (int j = 1; j < nz; j++)
			{
				Cell& cell = model->cells[j + k * nz + m * nx * nz];
				Point point = getCartesian<Cell>(cell.mu - cell.hmu / 2.0, cell.nu - cell.hnu / 2.0, cell.z - cell.hz / 2.0);
				points->InsertNextPoint(r_dim * point[0] * VIEW_MULTIPLIER,
					r_dim * point[1] * VIEW_MULTIPLIER,
					-r_dim * point[2]);
			}
		}
	}
	grid->SetPoints(points);

	// Data
	auto hexs = vtkSmartPointer<vtkCellArray>::New();
	auto pres = vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");
	auto p_bub = vtkSmartPointer<vtkDoubleArray>::New();
	p_bub->SetName("buble_point");
	auto sat_wat = vtkSmartPointer<vtkDoubleArray>::New();
	sat_wat->SetName("waterSaturation");
	auto sat_oil = vtkSmartPointer<vtkDoubleArray>::New();
	sat_oil->SetName("oilSaturation");
	auto sat_gas = vtkSmartPointer<vtkDoubleArray>::New();
	sat_gas->SetName("gasSaturation");
	auto satur = vtkSmartPointer<vtkIntArray>::New();
	satur->SetName("SATUR");
	auto temp = vtkSmartPointer<vtkDoubleArray>::New();
	temp->SetName("temperature");
	auto perm_xy = vtkSmartPointer<vtkDoubleArray>::New();
	perm_xy->SetName("perm_xy");
	auto vel_oil = vtkSmartPointer<vtkDoubleArray>::New();
	vel_oil->SetName("oilVelocity");
	vel_oil->SetNumberOfComponents(3);

	double vel[3];
	int l, j, k;
	Cell* neighbor[6];

	for (l = 0; l < ny - 1; l++)
	{
		// Left
		for (j = 0; j < nz - 2; j++)
		{
			Cell& cell = model->cells[l*nx*nz + j + 1];

			if (cell.isUsed)
			{
				vtkSmartPointer<vtkHexahedron> hex =
					vtkSmartPointer<vtkHexahedron>::New();
				hex->GetPointIds()->SetId(0, j + l * nx * (nz - 1));
				hex->GetPointIds()->SetId(1, j + (l + 1) * nx * (nz - 1));
				hex->GetPointIds()->SetId(2, j + (l + 1) * nx * (nz - 1) + 1);
				hex->GetPointIds()->SetId(3, j + l * nx * (nz - 1) + 1);

				hex->GetPointIds()->SetId(4, j + l * nx * (nz - 1) + (nz - 1));
				hex->GetPointIds()->SetId(5, j + (l + 1) * nx * (nz - 1) + (nz - 1));
				hex->GetPointIds()->SetId(6, j + (l + 1) * nx * (nz - 1) + 1 + (nz - 1));
				hex->GetPointIds()->SetId(7, j + l * nx * (nz - 1) + 1 + (nz - 1));

				hexs->InsertNextCell(hex);

				temp->InsertNextValue(cell.u_next.t * T_dim);
				pres->InsertNextValue(cell.u_next.p * P_dim);
				p_bub->InsertNextValue(cell.u_next.p_bub * P_dim / BAR_TO_PA);
				satur->InsertNextValue(cell.u_next.SATUR);
				sat_wat->InsertNextValue(cell.u_next.s_w);
				sat_oil->InsertNextValue(cell.u_next.s_o);
				sat_gas->InsertNextValue(1.0 - cell.u_next.s_w - cell.u_next.s_o);
				perm_xy->InsertNextValue(M2toMilliDarcy(model->getPerm_mu(cell) * r_dim * r_dim));

				model->getNeighbors(cell, neighbor);
				auto vel_cart = cell.getVectorCartesian(model->getOilVelocity(cell, neighbor, MU_AXIS),
					model->getOilVelocity(cell, neighbor, NU_AXIS),
					model->getOilVelocity(cell, neighbor, Z_AXIS));
				vel[0] = vel_cart[0] * r_dim / t_dim;
				vel[1] = vel_cart[1] * r_dim / t_dim;
				vel[2] = -vel_cart[2] * r_dim / t_dim;
				vel_oil->InsertNextTuple(vel);
			}
		}

		// Middle cells
		for (k = 1; k < nx - 1; k++)
		{
			for (j = 0; j < nz - 2; j++)
			{
				Cell& cell = model->cells[l*nx*nz + k*nz + j + 1];

				vtkSmartPointer<vtkHexahedron> hex =
					vtkSmartPointer<vtkHexahedron>::New();
				hex->GetPointIds()->SetId(0, j + k*(nz - 1) + l * nx * (nz - 1));
				hex->GetPointIds()->SetId(1, j + k*(nz - 1) + (l + 1) * nx * (nz - 1));
				hex->GetPointIds()->SetId(2, j + k*(nz - 1) + (l + 1) * nx * (nz - 1) + 1);
				hex->GetPointIds()->SetId(3, j + k*(nz - 1) + l * nx * (nz - 1) + 1);

				hex->GetPointIds()->SetId(4, j + (k + 1)*(nz - 1) + l * nx * (nz - 1));
				hex->GetPointIds()->SetId(5, j + (k + 1)*(nz - 1) + (l + 1) * nx * (nz - 1));
				hex->GetPointIds()->SetId(6, j + (k + 1)*(nz - 1) + (l + 1) * nx * (nz - 1) + 1);
				hex->GetPointIds()->SetId(7, j + (k + 1)*(nz - 1) + l * nx * (nz - 1) + 1);

				hexs->InsertNextCell(hex);

				temp->InsertNextValue(cell.u_next.t * T_dim);
				pres->InsertNextValue(cell.u_next.p * P_dim);
				p_bub->InsertNextValue(cell.u_next.p_bub * P_dim / BAR_TO_PA);
				satur->InsertNextValue(cell.u_next.SATUR);
				sat_wat->InsertNextValue(cell.u_next.s_w);
				sat_oil->InsertNextValue(cell.u_next.s_o);
				sat_gas->InsertNextValue(1.0 - cell.u_next.s_w - cell.u_next.s_o);
				perm_xy->InsertNextValue(M2toMilliDarcy(model->getPerm_mu(cell) * r_dim * r_dim));

				model->getNeighbors(cell, neighbor);
				auto vel_cart = cell.getVectorCartesian(model->getOilVelocity(cell, neighbor, MU_AXIS),
					model->getOilVelocity(cell, neighbor, NU_AXIS),
					model->getOilVelocity(cell, neighbor, Z_AXIS));
				vel[0] = vel_cart[0] * r_dim / t_dim;
				vel[1] = vel_cart[1] * r_dim / t_dim;
				vel[2] = -vel_cart[2] * r_dim / t_dim;
				vel_oil->InsertNextTuple(vel);
			}
		}
	}

	// Last segment - left
	for (j = 0; j < nz - 2; j++)
	{
		Cell& cell = model->cells[l*nx*nz + j + 1];

		if (cell.isUsed)
		{
			vtkSmartPointer<vtkHexahedron> hex =
				vtkSmartPointer<vtkHexahedron>::New();
			hex->GetPointIds()->SetId(0, j + l * nx * (nz - 1));
			hex->GetPointIds()->SetId(1, j);
			hex->GetPointIds()->SetId(2, j + 1);
			hex->GetPointIds()->SetId(3, j + l * nx * (nz - 1) + 1);

			hex->GetPointIds()->SetId(4, j + l * nx * (nz - 1) + (nz - 1));
			hex->GetPointIds()->SetId(5, j + (nz - 1));
			hex->GetPointIds()->SetId(6, j + 1 + (nz - 1));
			hex->GetPointIds()->SetId(7, j + l * nx * (nz - 1) + 1 + (nz - 1));

			hexs->InsertNextCell(hex);

			temp->InsertNextValue(cell.u_next.t * T_dim);
			pres->InsertNextValue(cell.u_next.p * P_dim);
			p_bub->InsertNextValue(cell.u_next.p_bub * P_dim / BAR_TO_PA);
			satur->InsertNextValue(cell.u_next.SATUR);
			sat_wat->InsertNextValue(cell.u_next.s_w);
			sat_oil->InsertNextValue(cell.u_next.s_o);
			sat_gas->InsertNextValue(1.0 - cell.u_next.s_w - cell.u_next.s_o);
			perm_xy->InsertNextValue(M2toMilliDarcy(model->getPerm_mu(cell) * r_dim * r_dim));

			model->getNeighbors(cell, neighbor);
			auto vel_cart = cell.getVectorCartesian(model->getOilVelocity(cell, neighbor, MU_AXIS),
				model->getOilVelocity(cell, neighbor, NU_AXIS),
				model->getOilVelocity(cell, neighbor, Z_AXIS));
			vel[0] = vel_cart[0] * r_dim / t_dim;
			vel[1] = vel_cart[1] * r_dim / t_dim;
			vel[2] = -vel_cart[2] * r_dim / t_dim;
			vel_oil->InsertNextTuple(vel);
		}
	}

	// Last segment - middle
	for (k = 1; k < nx - 1; k++)
	{
		for (j = 0; j < nz - 2; j++)
		{
			Cell& cell = model->cells[l*nx*nz + k*nz + j + 1];

			vtkSmartPointer<vtkHexahedron> hex =
				vtkSmartPointer<vtkHexahedron>::New();
			hex->GetPointIds()->SetId(0, j + k*(nz - 1) + l * nx * (nz - 1));
			hex->GetPointIds()->SetId(1, j + k*(nz - 1));
			hex->GetPointIds()->SetId(2, j + k*(nz - 1) + 1);
			hex->GetPointIds()->SetId(3, j + k*(nz - 1) + l * nx * (nz - 1) + 1);

			hex->GetPointIds()->SetId(4, j + (k + 1)*(nz - 1) + l * nx * (nz - 1));
			hex->GetPointIds()->SetId(5, j + (k + 1)*(nz - 1));
			hex->GetPointIds()->SetId(6, j + (k + 1)*(nz - 1) + 1);
			hex->GetPointIds()->SetId(7, j + (k + 1)*(nz - 1) + l * nx * (nz - 1) + 1);

			hexs->InsertNextCell(hex);

			pres->InsertNextValue(cell.u_next.p * P_dim);
			temp->InsertNextValue(cell.u_next.t * T_dim);
			p_bub->InsertNextValue(cell.u_next.p_bub * P_dim / BAR_TO_PA);
			satur->InsertNextValue(cell.u_next.SATUR);
			sat_wat->InsertNextValue(cell.u_next.s_w);
			sat_oil->InsertNextValue(cell.u_next.s_o);
			sat_gas->InsertNextValue(1.0 - cell.u_next.s_w - cell.u_next.s_o);
			perm_xy->InsertNextValue(M2toMilliDarcy(model->getPerm_mu(cell) * r_dim * r_dim));

			model->getNeighbors(cell, neighbor);
			auto vel_cart = cell.getVectorCartesian(model->getOilVelocity(cell, neighbor, MU_AXIS),
				model->getOilVelocity(cell, neighbor, NU_AXIS),
				model->getOilVelocity(cell, neighbor, Z_AXIS));
			vel[0] = vel_cart[0] * r_dim / t_dim;
			vel[1] = vel_cart[1] * r_dim / t_dim;
			vel[2] = -vel_cart[2] * r_dim / t_dim;
			vel_oil->InsertNextTuple(vel);
		}
	}

	grid->SetCells(VTK_HEXAHEDRON, hexs);
	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(pres);
	fd->AddArray(p_bub);
	fd->AddArray(satur);
	fd->AddArray(sat_wat);
	fd->AddArray(sat_oil);
	fd->AddArray(sat_gas);
	fd->AddArray(temp);
	fd->AddArray(perm_xy);
	fd->AddArray(vel_oil);

	// Writing
	auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(getFileName(snap_idx).c_str());
	writer->SetInputData(grid);
	writer->Write();
}
void VTKSnapshotter<blackoil_rz::BlackOil_RZ>::dump_all(int snap_idx)
{
	using namespace blackoil_rz;

	auto grid = vtkSmartPointer<vtkPolyData>::New();
	auto points = vtkSmartPointer<vtkPoints>::New();
	for (int j = 1; j < ny; j++)
	{
		Cell& cell = model->cells[j];
		points->InsertNextPoint(r_dim * (0.9 * cell.r), -r_dim * (cell.z - cell.hz / 2.0), 0.0);
	}
	for (int k = 1; k < nx; k++)
	{
		for (int j = 1; j < ny; j++)
		{
			Cell& cell = model->cells[k * ny + j];
			points->InsertNextPoint(r_dim * (cell.r - cell.hr / 2.0), -r_dim * (cell.z - cell.hz / 2.0), 0.0);
		}
	}
	grid->SetPoints(points);

	// Data
	auto polygons = vtkSmartPointer<vtkCellArray>::New();
	auto polygon = vtkSmartPointer<vtkPolygon>::New();
	polygon->GetPointIds()->SetNumberOfIds(4);

	auto pres = vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");
	auto p_bub = vtkSmartPointer<vtkDoubleArray>::New();
	p_bub->SetName("buble_point");
	auto sat_wat = vtkSmartPointer<vtkDoubleArray>::New();
	sat_wat->SetName("waterSaturation");
	auto sat_oil = vtkSmartPointer<vtkDoubleArray>::New();
	sat_oil->SetName("oilSaturation");
	auto sat_gas = vtkSmartPointer<vtkDoubleArray>::New();
	sat_gas->SetName("gasSaturation");
	auto satur = vtkSmartPointer<vtkIntArray>::New();
	satur->SetName("SATUR");
	auto vel_oil = vtkSmartPointer<vtkDoubleArray>::New();
	vel_oil->SetName("oilVelocity");
	vel_oil->SetNumberOfComponents(3);
	auto vel_gas = vtkSmartPointer<vtkDoubleArray>::New();
	vel_gas->SetName("gasVelocity");
	vel_gas->SetNumberOfComponents(3);

	int k, j, idx, idx1;
	double vel[3];

	for (j = 0; j < ny - 2; j++)
	{
		idx = j;
		idx1 = j + 1;
		Cell& cell = model->cells[idx1];

		polygon->GetPointIds()->SetId(0, idx);
		polygon->GetPointIds()->SetId(1, idx + ny - 1);
		polygon->GetPointIds()->SetId(2, idx + ny);
		polygon->GetPointIds()->SetId(3, idx + 1);
		polygons->InsertNextCell(polygon);

		pres->InsertNextValue(cell.u_next.p * P_dim / BAR_TO_PA);
		p_bub->InsertNextValue(cell.u_next.p_bub * P_dim / BAR_TO_PA);
		satur->InsertNextValue(cell.u_next.SATUR);
		sat_wat->InsertNextValue(cell.u_next.s_w);
		sat_oil->InsertNextValue(cell.u_next.s_o);
		sat_gas->InsertNextValue(1.0 - cell.u_next.s_w - cell.u_next.s_o);
		vel[0] = 0.0;
		vel[1] = 0.0;
		vel[2] = 0.0;
		vel_oil->InsertNextTuple(vel);
		vel[0] = 0.0;
		vel[1] = 0.0;
		vel[2] = 0.0;
		vel_gas->InsertNextTuple(vel);
		//vel[0] = model->getOilVelocity(cell, NEXT, R_AXIS);	vel[1] = model->getOilVelocity(cell, NEXT, Z_AXIS);	
		//vel_oil->InsertNextValue(vel);
	}

	// Middle cells
	for (k = 1; k < nx - 1; k++)
	{
		for (j = 0; j < ny - 2; j++)
		{
			idx = k * (ny - 1) + j;
			idx1 = k * ny + j + 1;
			Cell& cell = model->cells[idx1];

			polygon->GetPointIds()->SetId(0, idx);
			polygon->GetPointIds()->SetId(1, idx + ny - 1);
			polygon->GetPointIds()->SetId(2, idx + ny);
			polygon->GetPointIds()->SetId(3, idx + 1);
			polygons->InsertNextCell(polygon);

			pres->InsertNextValue(cell.u_next.p * P_dim / BAR_TO_PA);
			p_bub->InsertNextValue(cell.u_next.p_bub * P_dim / BAR_TO_PA);
			satur->InsertNextValue(cell.u_next.SATUR);
			sat_wat->InsertNextValue(cell.u_next.s_w);
			sat_oil->InsertNextValue(cell.u_next.s_o);
			sat_gas->InsertNextValue(1.0 - cell.u_next.s_w - cell.u_next.s_o);
			vel[0] = 0.0;
			vel[1] = 0.0;
			vel[2] = 0.0;
			vel_oil->InsertNextTuple(vel);
			vel[0] = 0.0;
			vel[1] = 0.0;
			vel[2] = 0.0;
			vel_gas->InsertNextTuple(vel);
			//vel[0] = model->getOilVelocity(cell, NEXT, R_AXIS);	vel[1] = model->getOilVelocity(cell, NEXT, Z_AXIS);	
			//vel_oil->InsertNextValue(vel);
		}
	}

	grid->SetPolys(polygons);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(pres);
	fd->AddArray(p_bub);
	fd->AddArray(satur);
	fd->AddArray(sat_wat);
	fd->AddArray(sat_oil);
	fd->AddArray(sat_gas);
	fd->AddArray(vel_oil);
	fd->AddArray(vel_gas);

	// Writing
	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

	writer->SetFileName(getFileName(snap_idx).c_str());
	writer->SetInputData(grid);
	writer->Write();
}
void VTKSnapshotter<wax_nit::WaxNIT>::dump_all(int snap_idx)
{
	using namespace wax_nit;

	auto grid = vtkSmartPointer<vtkPolyData>::New();
	auto points = vtkSmartPointer<vtkPoints>::New();
	for (int j = 1; j < ny; j++)
	{
		Cell& cell = model->cells[j];
		points->InsertNextPoint(r_dim * (0.9 * cell.r), -r_dim * (cell.z - cell.hz / 2.0), 0.0);
	}
	for (int k = 1; k < nx; k++)
	{
		for (int j = 1; j < ny; j++)
		{
			Cell& cell = model->cells[k * ny + j];
			points->InsertNextPoint(r_dim * (cell.r - cell.hr / 2.0), -r_dim * (cell.z - cell.hz / 2.0), 0.0);
		}
	}
	grid->SetPoints(points);

	// Data
	auto polygons = vtkSmartPointer<vtkCellArray>::New();
	auto polygon = vtkSmartPointer<vtkPolygon>::New();
	polygon->GetPointIds()->SetNumberOfIds(4);

	auto poro = vtkSmartPointer<vtkDoubleArray>::New();
	poro->SetName("porosity");
	auto perm = vtkSmartPointer<vtkDoubleArray>::New();
	perm->SetName("permeability");
	auto perm_dim = vtkSmartPointer<vtkDoubleArray>::New();
	perm_dim->SetName("permeability_dimless");
	auto temp = vtkSmartPointer<vtkDoubleArray>::New();
	temp->SetName("temperature");
	auto temp_sat = vtkSmartPointer<vtkDoubleArray>::New();
	temp_sat->SetName("wax_saturation_temp");
	auto pres = vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");
	auto p_bub = vtkSmartPointer<vtkDoubleArray>::New();
	p_bub->SetName("buble_point");
	auto sat_wat = vtkSmartPointer<vtkDoubleArray>::New();
	sat_wat->SetName("waterSaturation");
	auto sat_oil = vtkSmartPointer<vtkDoubleArray>::New();
	sat_oil->SetName("oilSaturation");
	auto sat_gas = vtkSmartPointer<vtkDoubleArray>::New();
	sat_gas->SetName("gasSaturation");
	auto sat_wax = vtkSmartPointer<vtkDoubleArray>::New();
	sat_wax->SetName("waxSaturation");
	auto satur_gas = vtkSmartPointer<vtkIntArray>::New();
	satur_gas->SetName("satur_gas");
	auto satur_wax = vtkSmartPointer<vtkIntArray>::New();
	satur_wax->SetName("satur_wax");
	auto vel_oil = vtkSmartPointer<vtkDoubleArray>::New();
	vel_oil->SetName("oilVelocity");
	vel_oil->SetNumberOfComponents(3);
	auto vel_gas = vtkSmartPointer<vtkDoubleArray>::New();
	vel_gas->SetName("gasVelocity");
	vel_gas->SetNumberOfComponents(3);
	auto vel_wat = vtkSmartPointer<vtkDoubleArray>::New();
	vel_wat->SetName("waterVelocity");
	vel_wat->SetNumberOfComponents(3);

	int k, j, idx, idx1;
	double vel[3];

	for (j = 0; j < ny - 2; j++)
	{
		idx = j;
		idx1 = j + 1;
		Cell& cell = model->cells[idx1];

		polygon->GetPointIds()->SetId(0, idx);
		polygon->GetPointIds()->SetId(1, idx + ny - 1);
		polygon->GetPointIds()->SetId(2, idx + ny);
		polygon->GetPointIds()->SetId(3, idx + 1);
		polygons->InsertNextCell(polygon);

		poro->InsertNextValue(cell.u_next.m);
		perm->InsertNextValue(M2toMilliDarcy(cell.props->getPermCoseni_r(cell.u_next.m).value() * r_dim * r_dim));
		perm_dim->InsertNextValue(cell.props->getPermCoseni_r(cell.u_next.m).value() / cell.props->getPermCoseni_r(cell.props->m_init).value());
		temp->InsertNextValue((cell.u_next.t - cell.props->t_init) * T_dim);
		temp_sat->InsertNextValue((cell.u_next.t_bub - cell.props->t_init) * T_dim);
		pres->InsertNextValue(cell.u_next.p * P_dim / BAR_TO_PA);
		p_bub->InsertNextValue(cell.u_next.p_bub * P_dim / BAR_TO_PA);
		satur_gas->InsertNextValue(cell.u_next.satur_gas);
		satur_wax->InsertNextValue(cell.u_next.satur_wax);
		sat_wat->InsertNextValue(cell.u_next.s_w);
		sat_oil->InsertNextValue(cell.u_next.s_o);
		sat_gas->InsertNextValue(cell.u_next.s_g);
		sat_wax->InsertNextValue(1.0 - cell.u_next.s_w - cell.u_next.s_o - cell.u_next.s_g);
		vel[0] = 0.0;	vel[1] = 0.0;	vel[2] = 0.0;
		vel_oil->InsertNextTuple(vel);
		vel_gas->InsertNextTuple(vel);
		vel_wat->InsertNextTuple(vel);
	}

	// Middle cells
	for (k = 1; k < nx - 1; k++)
	{
		for (j = 0; j < ny - 2; j++)
		{
			idx = k * (ny - 1) + j;
			idx1 = k * ny + j + 1;
			Cell& cell = model->cells[idx1];

			polygon->GetPointIds()->SetId(0, idx);
			polygon->GetPointIds()->SetId(1, idx + ny - 1);
			polygon->GetPointIds()->SetId(2, idx + ny);
			polygon->GetPointIds()->SetId(3, idx + 1);
			polygons->InsertNextCell(polygon);

			poro->InsertNextValue(cell.u_next.m);
			perm->InsertNextValue(M2toMilliDarcy(cell.props->getPermCoseni_r(cell.u_next.m).value() * r_dim * r_dim));
			perm_dim->InsertNextValue(cell.props->getPermCoseni_r(cell.u_next.m).value() / cell.props->getPermCoseni_r(cell.props->m_init).value());
			temp->InsertNextValue((cell.u_next.t - cell.props->t_init) * T_dim);
			temp_sat->InsertNextValue((cell.u_next.t_bub - cell.props->t_init) * T_dim);
			pres->InsertNextValue(cell.u_next.p * P_dim / BAR_TO_PA);
			p_bub->InsertNextValue(cell.u_next.p_bub * P_dim / BAR_TO_PA);
			satur_gas->InsertNextValue(cell.u_next.satur_gas);
			satur_wax->InsertNextValue(cell.u_next.satur_wax);
			sat_wat->InsertNextValue(cell.u_next.s_w);
			sat_oil->InsertNextValue(cell.u_next.s_o);
			sat_gas->InsertNextValue(cell.u_next.s_g);
			sat_wax->InsertNextValue(1.0 - cell.u_next.s_w - cell.u_next.s_o - cell.u_next.s_g);
			vel[0] = r_dim / t_dim * model->getOilVel(cell, R_AXIS);
			vel[1] = r_dim / t_dim * model->getOilVel(cell, Z_AXIS);
			vel[2] = 0.0;
			vel_oil->InsertNextTuple(vel);
			vel[0] = r_dim / t_dim * model->getGasVel(cell, R_AXIS);
			vel[1] = r_dim / t_dim * model->getGasVel(cell, Z_AXIS);
			vel[2] = 0.0;
			vel_gas->InsertNextTuple(vel);
			vel[0] = r_dim / t_dim * model->getWatVel(cell, R_AXIS);
			vel[1] = r_dim / t_dim * model->getWatVel(cell, Z_AXIS);
			vel[2] = 0.0;
			vel_wat->InsertNextTuple(vel);
		}
	}

	grid->SetPolys(polygons);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(poro);
	fd->AddArray(perm);
	fd->AddArray(perm_dim);
	fd->AddArray(temp);
	fd->AddArray(temp_sat);
	fd->AddArray(pres);
	fd->AddArray(p_bub);
	fd->AddArray(satur_gas);
	fd->AddArray(satur_wax);
	fd->AddArray(sat_wat);
	fd->AddArray(sat_oil);
	fd->AddArray(sat_gas);
	fd->AddArray(sat_wax);
	fd->AddArray(vel_oil);
	fd->AddArray(vel_gas);
	fd->AddArray(vel_wat);

	// Writing
	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

	writer->SetFileName(getFileName(snap_idx).c_str());
	writer->SetInputData(grid);
	writer->Write();
}
void VTKSnapshotter<wax_nit1d::WaxNIT1d>::dump_all(int i)
{
	using namespace wax_nit1d;

	// Grid
	auto grid = vtkSmartPointer<vtkPolyData>::New();
	// Points
	auto points = vtkSmartPointer<vtkPoints>::New();

	points->InsertNextPoint(r_dim * (0.99 * model->cells[0].x), 0.0, 0.0);
	points->InsertNextPoint(r_dim * (0.99 * model->cells[0].x), -r_dim * model->props_sk.height, 0.0);
	for (int k = 1; k < nx; k++)
	{
		Cell& cell = model->cells[k];
		points->InsertNextPoint(r_dim * (cell.x - cell.hx / 2.0), 0.0, 0.0);
		points->InsertNextPoint(r_dim * (cell.x - cell.hx / 2.0), -r_dim * model->props_sk.height, 0.0);
	}
	grid->SetPoints(points);

	// Data
	auto polygons = vtkSmartPointer<vtkCellArray>::New();
	auto polygon = vtkSmartPointer<vtkPolygon>::New();
	polygon->GetPointIds()->SetNumberOfIds(4);
	auto poro = vtkSmartPointer<vtkDoubleArray>::New();
	poro->SetName("porosity");
	auto perm = vtkSmartPointer<vtkDoubleArray>::New();
	perm->SetName("permeability");
	auto perm_dim = vtkSmartPointer<vtkDoubleArray>::New();
	perm_dim->SetName("permeability_dimless");
	auto pres = vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");
	auto sat_oil = vtkSmartPointer<vtkDoubleArray>::New();
	sat_oil->SetName("oilSaturation");
	auto sat_wax = vtkSmartPointer<vtkDoubleArray>::New();
	sat_wax->SetName("waxSaturation");
	auto vel_oil = vtkSmartPointer<vtkDoubleArray>::New();
	vel_oil->SetName("oilVelocity");
	vel_oil->SetNumberOfComponents(3);

	int k, j, idx, idx1;

	double vel[3];

	idx = 0;
	Variable& next = model->cells[0].u_next;
	const auto& props = model->props_sk;

	polygon->GetPointIds()->SetId(0, idx);
	polygon->GetPointIds()->SetId(1, idx + ny - 1);
	polygon->GetPointIds()->SetId(2, idx + ny);
	polygon->GetPointIds()->SetId(3, idx + 1);
	polygons->InsertNextCell(polygon);

	poro->InsertNextValue(next.m);
	perm->InsertNextValue(M2toMilliDarcy(props.getPermCoseni(next.m).value() * r_dim * r_dim));
	perm_dim->InsertNextValue(props.getPermCoseni(next.m).value() / props.getPermCoseni(props.m_init).value());
	pres->InsertNextValue(next.p * P_dim / BAR_TO_PA);
	sat_oil->InsertNextValue(1.0 - next.s_p);
	sat_wax->InsertNextValue(next.s_p);
	vel[0] = 0.0;	vel[1] = 0.0;	vel[2] = 0.0;
	vel_oil->InsertNextTuple(vel);

	// Middle cells
	for (k = 1; k < nx - 1; k++)
	{
		idx = k * (ny - 1);
		idx1 = k;
		Cell& cell = model->cells[idx1];
		Variable& next = cell.u_next;

		polygon->GetPointIds()->SetId(0, idx);
		polygon->GetPointIds()->SetId(1, idx + ny - 1);
		polygon->GetPointIds()->SetId(2, idx + ny);
		polygon->GetPointIds()->SetId(3, idx + 1);
		polygons->InsertNextCell(polygon);

		poro->InsertNextValue(cell.u_next.m);
		perm->InsertNextValue(M2toMilliDarcy(props.getPermCoseni(cell.u_next.m).value() * r_dim * r_dim));
		perm_dim->InsertNextValue(props.getPermCoseni(cell.u_next.m).value() / props.getPermCoseni(props.m_init).value());
		pres->InsertNextValue(cell.u_next.p * P_dim / BAR_TO_PA);
		sat_oil->InsertNextValue(1.0 - next.s_p);
		sat_wax->InsertNextValue(next.s_p);
		vel[0] = r_dim / t_dim * model->getOilVel(cell);
		vel[1] = 0.0;
		vel[2] = 0.0;
		vel_oil->InsertNextTuple(vel);
	}

	grid->SetPolys(polygons);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(poro);
	fd->AddArray(perm);
	fd->AddArray(perm_dim);
	fd->AddArray(pres);
	fd->AddArray(sat_oil);
	fd->AddArray(sat_wax);
	fd->AddArray(vel_oil);
	// Writing
	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(getFileName(i).c_str());
	writer->SetInputData(grid);
	writer->Write();
}

template class VTKSnapshotter<gasOil_rz::GasOil_RZ>;
template class VTKSnapshotter<acid2d::Acid2d>;
template class VTKSnapshotter<acid2dnit::Acid2dNIT>;
template class VTKSnapshotter<acid1d::Acid1d>;
template class VTKSnapshotter<wax_nit::WaxNIT>;
template class VTKSnapshotter<wax_nit1d::WaxNIT1d>;
template class VTKSnapshotter<acidfrac::AcidFrac>;
template class VTKSnapshotter<acidellfrac::AcidEllFrac>;
template class VTKSnapshotter<acidrecfrac::AcidRecFrac>;
template class VTKSnapshotter<acidrecfrac_prod::RecFracProd>;
template class VTKSnapshotter<vpp2d::VPP2d>;
template class VTKSnapshotter<bing1d::Bingham1d>;
template class VTKSnapshotter<gasOil_elliptic::GasOil_Elliptic>;
template class Snapshotter<oilnit_elliptic::OilNIT_Elliptic>;
template class Snapshotter<blackoilnit_elliptic::BlackOilNIT_Elliptic>;
template class Snapshotter<blackoil_rz::BlackOil_RZ>;