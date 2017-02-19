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
#include "model/Acid/2d/Acid2d.hpp"
#include "model/VPP2d/VPP2d.hpp"
#include "model/Bingham1d/Bingham1d.hpp"
#include "model/GasOil_Elliptic/GasOil_Elliptic.hpp"
#include "model/GasOilNIT_Elliptic/GasOilNIT_Elliptic.hpp"
#include "model/OilNIT_Elliptic/OilNIT_Elliptic.hpp"

#include <cmath>

#define VIEW_MULTIPLIER 0.01

using namespace std;

template <class modelType>
VTKSnapshotter<modelType>::VTKSnapshotter()
{
	pattern = prefix + "snap_%{STEP}.vtp";
}
template <>
VTKSnapshotter<gasOil_rz::GasOil_RZ>::VTKSnapshotter()
{
	pattern = prefix + "GasOil_RZ_%{STEP}.vtp";
}
template <>
VTKSnapshotter<acid2d::Acid2d>::VTKSnapshotter()
{
	pattern = prefix + "Acid2d_%{STEP}.vtp";
}
template <>
VTKSnapshotter<vpp2d::VPP2d>::VTKSnapshotter()
{
	pattern = prefix + "VPP2d_%{STEP}.vtp";
}
VTKSnapshotter<bing1d::Bingham1d>::VTKSnapshotter()
{
	pattern = prefix + "Bing1d_%{STEP}.vtp";
}
VTKSnapshotter<gasOil_elliptic::GasOil_Elliptic>::VTKSnapshotter()
{
	pattern = prefix + "GasOil_El_%{STEP}.vtu";
}
VTKSnapshotter<oilnit_elliptic::OilNIT_Elliptic>::VTKSnapshotter()
{
	pattern = prefix + "OilNIT_El_%{STEP}.vtu";
}
VTKSnapshotter<gasOilnit_elliptic::GasOilNIT_Elliptic>::VTKSnapshotter()
{
	pattern = prefix + "GasOilNIT_El_%{STEP}.vtu";
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

template <>
void VTKSnapshotter<gasOil_rz::GasOil_RZ>::dump_all(int i)
{
	using gasOil_rz::Cell;

	// Grid
	vtkSmartPointer<vtkPolyData> grid =
		vtkSmartPointer<vtkPolyData>::New();

	// Points
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();
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
	vtkSmartPointer<vtkCellArray> polygons =
		vtkSmartPointer<vtkCellArray>::New();

	vtkSmartPointer<vtkPolygon> polygon = 
		vtkSmartPointer<vtkPolygon>::New();
	polygon->GetPointIds()->SetNumberOfIds(4);

	vtkSmartPointer<vtkDoubleArray> pres = 
		vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");

	vtkSmartPointer<vtkDoubleArray> p_bub =
		vtkSmartPointer<vtkDoubleArray>::New();
	p_bub->SetName("buble_point");

	vtkSmartPointer<vtkDoubleArray> sat_oil = 
		vtkSmartPointer<vtkDoubleArray>::New();
	sat_oil->SetName("oilSaturation");

	vtkSmartPointer<vtkDoubleArray> sat_gas = 
		vtkSmartPointer<vtkDoubleArray>::New();
	sat_gas->SetName("gasSaturation");

	vtkSmartPointer<vtkIntArray> satur =
		vtkSmartPointer<vtkIntArray>::New();
	satur->SetName("SATUR");

	vtkSmartPointer<vtkDoubleArray> vel_oil = 
		vtkSmartPointer<vtkDoubleArray>::New();
	vel_oil->SetName("oilVelocity");
	vel_oil->SetNumberOfComponents(3);

	vtkSmartPointer<vtkDoubleArray> vel_gas = 
		vtkSmartPointer<vtkDoubleArray>::New();
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
	for(k = 1; k < nx-2; k++)
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
	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	
	writer->SetFileName( getFileName(i).c_str() );
	writer->SetInputData(grid);
	writer->Write();
}
template <>
void VTKSnapshotter<acid2d::Acid2d>::dump_all(int i)
{
	using acid2d::Cell;

	// Grid
	vtkSmartPointer<vtkPolyData> grid =
		vtkSmartPointer<vtkPolyData>::New();

	// Points
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();
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
	vtkSmartPointer<vtkCellArray> polygons =
		vtkSmartPointer<vtkCellArray>::New();

	vtkSmartPointer<vtkPolygon> polygon =
		vtkSmartPointer<vtkPolygon>::New();
	polygon->GetPointIds()->SetNumberOfIds(4);

	vtkSmartPointer<vtkDoubleArray> poro =
		vtkSmartPointer<vtkDoubleArray>::New();
	poro->SetName("porosity");

	vtkSmartPointer<vtkDoubleArray> pres =
		vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");

	vtkSmartPointer<vtkDoubleArray> sat_l =
		vtkSmartPointer<vtkDoubleArray>::New();
	sat_l->SetName("LiquidSaturation");

	vtkSmartPointer<vtkDoubleArray> sat_g =
		vtkSmartPointer<vtkDoubleArray>::New();
	sat_g->SetName("GasSaturation");

	vtkSmartPointer<vtkDoubleArray> conc_a =
		vtkSmartPointer<vtkDoubleArray>::New();
	conc_a->SetName("AcidConcentration");

	vtkSmartPointer<vtkDoubleArray> conc_s =
		vtkSmartPointer<vtkDoubleArray>::New();
	conc_s->SetName("SaltConcentration");

	vtkSmartPointer<vtkDoubleArray> vel_l =
		vtkSmartPointer<vtkDoubleArray>::New();
	vel_l->SetName("LiquidVelocity");
	vel_l->SetNumberOfComponents(3);

	vtkSmartPointer<vtkDoubleArray> vel_g =
		vtkSmartPointer<vtkDoubleArray>::New();
	vel_g->SetName("gasVelocity");
	vel_g->SetNumberOfComponents(3);

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
		pres->InsertNextValue(cell.u_next.p);
		sat_l->InsertNextValue(cell.u_next.s);
		sat_g->InsertNextValue(1.0 - cell.u_next.s);
		conc_a->InsertNextValue(1.0 - cell.u_next.Ya);
		conc_s->InsertNextValue(1.0 - cell.u_next.Ys);
		vel[0] = 0.0;
		vel[1] = 0.0;
		vel[2] = 0.0;
		vel_l->InsertNextTuple(vel);
		vel[0] = 0.0;
		vel[1] = 0.0;
		vel[2] = 0.0;
		vel_g->InsertNextTuple(vel);
	}

	// Middle cells
	for (k = 1; k < nx - 2; k++)
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
			pres->InsertNextValue(cell.u_next.p);
			sat_l->InsertNextValue(cell.u_next.s);
			sat_g->InsertNextValue(1.0 - cell.u_next.s);
			conc_a->InsertNextValue(1.0 - cell.u_next.Ya);
			conc_s->InsertNextValue(1.0 - cell.u_next.Ys);
			vel[0] = r_dim / t_dim * model->getLiquidVelocity(cell, NEXT, R_AXIS);
			vel[1] = r_dim / t_dim * model->getLiquidVelocity(cell, NEXT, Z_AXIS);
			vel[2] = 0.0;
			vel_l->InsertNextTuple(vel);
			vel[0] = r_dim / t_dim * model->getGasVelocity(cell, NEXT, R_AXIS);
			vel[1] = r_dim / t_dim * model->getGasVelocity(cell, NEXT, Z_AXIS);
			vel[2] = 0.0;
			vel_g->InsertNextTuple(vel);
		}
	}

	grid->SetPolys(polygons);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(poro);
	fd->AddArray(pres);
	fd->AddArray(sat_l);
	fd->AddArray(sat_g);
	fd->AddArray(conc_a);
	fd->AddArray(conc_s);
	fd->AddArray(vel_l);
	fd->AddArray(vel_g);

	// Writing
	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();

	writer->SetFileName(getFileName(i).c_str());
	writer->SetInputData(grid);
	writer->Write();
}
template <>
void VTKSnapshotter<vpp2d::VPP2d>::dump_all(int i)
{
	using vpp2d::Cell;

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
	for (k = 1; k < nx - 2; k++)
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
template <>
void VTKSnapshotter<bing1d::Bingham1d>::dump_all(int i)
{
	using bing1d::Cell;

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
template <>
void VTKSnapshotter<gasOil_elliptic::GasOil_Elliptic>::dump_all(int snap_idx)
{
	using gasOil_elliptic::Cell;
	using gasOil_elliptic::Point;

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
	using oilnit_elliptic::Cell;
	using oilnit_elliptic::Point;

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
	auto polygon = vtkSmartPointer<vtkPolygon>::New();
	polygon->GetPointIds()->SetNumberOfIds(4);
	auto pres = vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");
	auto temp = vtkSmartPointer<vtkDoubleArray>::New();
	temp->SetName("temperature");
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

				temp->InsertNextValue(cell.u_next.t * T_dim);
				pres->InsertNextValue(cell.u_next.p * P_dim);
				vel[0] = 0.0;
				vel[1] = 0.0;
				vel[2] = 0.0;
				vel_oil->InsertNextTuple(vel);
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

				temp->InsertNextValue(cell.u_next.t * T_dim);
				pres->InsertNextValue(cell.u_next.p * P_dim);
				vel[0] = 0.0;// r_dim / t_dim * (cos(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) - sin(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
				vel[1] = 0.0;// r_dim / t_dim * (sin(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) + cos(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
				vel[2] = 0.0;// r_dim / t_dim * model->getOilVelocity(cell, NEXT, Z_AXIS);
				vel_oil->InsertNextTuple(vel);
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

			temp->InsertNextValue(cell.u_next.t * T_dim);
			pres->InsertNextValue(cell.u_next.p * P_dim);

			model->getNeighbors(cell, neighbor);
			vel[0] = 0.0;
			vel[1] = 0.0;
			vel[2] = 0.0;
			vel_oil->InsertNextTuple(vel);
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
			temp->InsertNextValue(cell.u_next.t * T_dim);
			
			model->getNeighbors(cell, neighbor);
			vel[0] = 0.0;// r_dim / t_dim * (cos(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) - sin(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
			vel[1] = 0.0;// r_dim / t_dim * (sin(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) + cos(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
			vel[2] = 0.0;// r_dim / t_dim * model->getOilVelocity(cell, NEXT, Z_AXIS);
			vel_oil->InsertNextTuple(vel);
		}
	}

	grid->SetCells(VTK_HEXAHEDRON, hexs);
	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(pres);
	fd->AddArray(temp);
	fd->AddArray(vel_oil);

	// Writing
	auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(getFileName(snap_idx).c_str());
	writer->SetInputData(grid);
	writer->Write();
}

template class VTKSnapshotter<gasOil_rz::GasOil_RZ>;
template class VTKSnapshotter<acid2d::Acid2d>;
template class VTKSnapshotter<vpp2d::VPP2d>;
template class VTKSnapshotter<bing1d::Bingham1d>;
template class VTKSnapshotter<gasOil_elliptic::GasOil_Elliptic>;
template class Snapshotter<oilnit_elliptic::OilNIT_Elliptic>;
template class Snapshotter<gasOilnit_elliptic::GasOilNIT_Elliptic>;