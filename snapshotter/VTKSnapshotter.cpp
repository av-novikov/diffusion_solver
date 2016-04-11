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

#include "model/Oil1D/Oil1D.h"
#include "model/Gas1D/Gas1D.h"
#include "model/Gas1D/Gas1D_simple.h"
#include "model/Oil1D_NIT/Oil1D_NIT.h"
#include "model/Oil_RZ/Oil_RZ.h"
#include "model/Oil_RZ_NIT/Oil_RZ_NIT.h"
#include "model/GasOil_RZ/GasOil_RZ.h"
#include "model/GasOil_RZ_NIT/GasOil_RZ_NIT.h"

#include "model/3D/GasOil_3D/GasOil_3D.h"
#include "model/3D/GasOil_3D_NIT/GasOil_3D_NIT.h"

#include "model/3D/Perforation/GasOil_Perf.h"
#include "model/3D/Perforation/GasOil_Perf_NIT.h"

#include <cmath>

using namespace std;

template <class modelType>
VTKSnapshotter<modelType>::VTKSnapshotter()
{
	pattern = prefix + "snap_%{STEP}.vtp";
}

template <>
VTKSnapshotter<oil_rz::Oil_RZ>::VTKSnapshotter()
{
	pattern = prefix + "Oil_RZ_%{STEP}.vtp";
}

template <>
VTKSnapshotter<oil_rz_nit::Oil_RZ_NIT>::VTKSnapshotter()
{
	pattern = prefix + "Oil_RZ_NIT_%{STEP}.vtp";
}

template <>
VTKSnapshotter<gasOil_rz::GasOil_RZ>::VTKSnapshotter()
{
	pattern = prefix + "GasOil_RZ_%{STEP}.vtp";
}

template <>
VTKSnapshotter<gasOil_rz_NIT::GasOil_RZ_NIT>::VTKSnapshotter()
{
	pattern = prefix + "GasOil_RZ_NIT_%{STEP}.vtp";
}

template <>
VTKSnapshotter<gasOil_3d::GasOil_3D>::VTKSnapshotter()
{
	pattern = prefix + "GasOil_3D_%{STEP}.vtu";
}

template <>
VTKSnapshotter<gasOil_perf::GasOil_Perf>::VTKSnapshotter()
{
	pattern = prefix + "GasOil_Perf_%{STEP}.vtu";
}

template <>
VTKSnapshotter<gasOil_perf_nit::GasOil_Perf_NIT>::VTKSnapshotter()
{
	pattern = prefix + "GasOil_Perf_NIT_%{STEP}.vtu";
}

template <>
VTKSnapshotter<gasOil_3d_NIT::GasOil_3D_NIT>::VTKSnapshotter()
{
	pattern = prefix + "GasOil_3D_NIT_%{STEP}.vtu";
}

template <class modelType>
VTKSnapshotter<modelType>::~VTKSnapshotter()
{
}

template <class modelType>
void VTKSnapshotter<modelType>::dump(int i)
{
}

template <>
void VTKSnapshotter<gasOil_rz_NIT::GasOil_RZ_NIT>::dump(int i)
{
	using gasOil_rz_NIT::Cell;

	// Grid
	vtkSmartPointer<vtkPolyData> grid =
		vtkSmartPointer<vtkPolyData>::New();

	// Points
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();
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

	vtkSmartPointer<vtkDoubleArray> temp = 
		vtkSmartPointer<vtkDoubleArray>::New();
	temp->SetName("temperature");

	vtkSmartPointer<vtkDoubleArray> pres = 
		vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");

	vtkSmartPointer<vtkDoubleArray> sat_oil = 
		vtkSmartPointer<vtkDoubleArray>::New();
	sat_oil->SetName("oilSaturation");

	vtkSmartPointer<vtkDoubleArray> sat_gas = 
		vtkSmartPointer<vtkDoubleArray>::New();
	sat_gas->SetName("gasSaturation");

	vtkSmartPointer<vtkDoubleArray> vel_oil = 
		vtkSmartPointer<vtkDoubleArray>::New();
	vel_oil->SetName("oilVelocity");
	vel_oil->SetNumberOfComponents(3);

	vtkSmartPointer<vtkDoubleArray> vel_gas = 
		vtkSmartPointer<vtkDoubleArray>::New();
	vel_gas->SetName("gasVelocity");
	vel_gas->SetNumberOfComponents(3);

	int k, j, idx, idx1;

	// Left border
	/*polygon->GetPointIds()->SetId(0, 0);
	polygon->GetPointIds()->SetId(1, 0);
	polygon->GetPointIds()->SetId(2, 0);
	polygon->GetPointIds()->SetId(3, 0);
	grid->InsertNextCell(polygon->GetCellType(), polygon->GetPointIds());
	for(j = 0; j < ny-2; j++)
	{
		polygon->GetPointIds()->SetId(0, j);
		polygon->GetPointIds()->SetId(1, j);
		polygon->GetPointIds()->SetId(2, j+1);
		polygon->GetPointIds()->SetId(3, j+1);
		grid->InsertNextCell(polygon->GetCellType(), polygon->GetPointIds());
	}
	polygon->GetPointIds()->SetId(0, j);
	polygon->GetPointIds()->SetId(1, j);
	polygon->GetPointIds()->SetId(2, j);
	polygon->GetPointIds()->SetId(3, j);
	grid->InsertNextCell(polygon->GetCellType(), polygon->GetPointIds());*/

	double vel [3];

	// Middle cells
	for(k = 0; k < nx-2; k++)
	{
/*		idx = k * (ny-1);
		polygon->GetPointIds()->SetId(0, idx);
		polygon->GetPointIds()->SetId(1, idx+ny-1);
		polygon->GetPointIds()->SetId(2, idx+ny-1);
		polygon->GetPointIds()->SetId(3, idx);
		grid->InsertNextCell(polygon->GetCellType(), polygon->GetPointIds());*/
		for(j = 0; j < ny-2; j++)
		{
			idx = k * (ny-1) + j;
			idx1 = (k+1) * ny + j + 1;
			Cell& cell = model->cells[idx1];

			polygon->GetPointIds()->SetId(0, idx);
			polygon->GetPointIds()->SetId(1, idx+ny-1);
			polygon->GetPointIds()->SetId(2, idx+ny);
			polygon->GetPointIds()->SetId(3, idx+1);
			polygons->InsertNextCell(polygon);

			temp->InsertNextValue(cell.u_next.t * T_dim);
			pres->InsertNextValue(cell.u_next.p);
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
/*		idx++;
		polygon->GetPointIds()->SetId(0, idx);
		polygon->GetPointIds()->SetId(1, idx+ny-1);
		polygon->GetPointIds()->SetId(2, idx+ny-1);
		polygon->GetPointIds()->SetId(3, idx);
		grid->InsertNextCell(polygon->GetCellType(), polygon->GetPointIds());*/
	}
	
	// Right border
/*	idx = (nx-2) * (ny-1);
	polygon->GetPointIds()->SetId(0, idx);
	polygon->GetPointIds()->SetId(1, idx);
	polygon->GetPointIds()->SetId(2, idx);
	polygon->GetPointIds()->SetId(3, idx);
	grid->InsertNextCell(polygon->GetCellType(), polygon->GetPointIds());
	for(j = 0; j < ny-2; j++)
	{
		idx = (nx-2) * (ny-1) + j;
		polygon->GetPointIds()->SetId(0, idx);
		polygon->GetPointIds()->SetId(1, idx);
		polygon->GetPointIds()->SetId(2, idx+1);
		polygon->GetPointIds()->SetId(3, idx+1);
		grid->InsertNextCell(polygon->GetCellType(), polygon->GetPointIds());
	}
	idx++;
	polygon->GetPointIds()->SetId(0, idx);
	polygon->GetPointIds()->SetId(1, idx);
	polygon->GetPointIds()->SetId(2, idx);
	polygon->GetPointIds()->SetId(3, idx);
	grid->InsertNextCell(polygon->GetCellType(), polygon->GetPointIds());*/

	grid->SetPolys(polygons);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(temp);
	fd->AddArray(pres);
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
void VTKSnapshotter<gasOil_rz::GasOil_RZ>::dump(int i)
{
	using gasOil_rz::Cell;

	// Grid
	vtkSmartPointer<vtkPolyData> grid =
		vtkSmartPointer<vtkPolyData>::New();

	// Points
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();
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

	vtkSmartPointer<vtkDoubleArray> sat_oil = 
		vtkSmartPointer<vtkDoubleArray>::New();
	sat_oil->SetName("oilSaturation");

	vtkSmartPointer<vtkDoubleArray> sat_gas = 
		vtkSmartPointer<vtkDoubleArray>::New();
	sat_gas->SetName("gasSaturation");

	vtkSmartPointer<vtkDoubleArray> vel_oil = 
		vtkSmartPointer<vtkDoubleArray>::New();
	vel_oil->SetName("oilVelocity");
	vel_oil->SetNumberOfComponents(3);

	vtkSmartPointer<vtkDoubleArray> vel_gas = 
		vtkSmartPointer<vtkDoubleArray>::New();
	vel_gas->SetName("gasVelocity");
	vel_gas->SetNumberOfComponents(3);

	int k, j, idx, idx1;

	// Left border
	/*polygon->GetPointIds()->SetId(0, 0);
	polygon->GetPointIds()->SetId(1, 0);
	polygon->GetPointIds()->SetId(2, 0);
	polygon->GetPointIds()->SetId(3, 0);
	grid->InsertNextCell(polygon->GetCellType(), polygon->GetPointIds());
	for(j = 0; j < ny-2; j++)
	{
		polygon->GetPointIds()->SetId(0, j);
		polygon->GetPointIds()->SetId(1, j);
		polygon->GetPointIds()->SetId(2, j+1);
		polygon->GetPointIds()->SetId(3, j+1);
		grid->InsertNextCell(polygon->GetCellType(), polygon->GetPointIds());
	}
	polygon->GetPointIds()->SetId(0, j);
	polygon->GetPointIds()->SetId(1, j);
	polygon->GetPointIds()->SetId(2, j);
	polygon->GetPointIds()->SetId(3, j);
	grid->InsertNextCell(polygon->GetCellType(), polygon->GetPointIds());*/

	double vel [3];

	// Middle cells
	for(k = 0; k < nx-2; k++)
	{
/*		idx = k * (ny-1);
		polygon->GetPointIds()->SetId(0, idx);
		polygon->GetPointIds()->SetId(1, idx+ny-1);
		polygon->GetPointIds()->SetId(2, idx+ny-1);
		polygon->GetPointIds()->SetId(3, idx);
		grid->InsertNextCell(polygon->GetCellType(), polygon->GetPointIds());*/
		for(j = 0; j < ny-2; j++)
		{
			idx = k * (ny-1) + j;
			idx1 = (k+1) * ny + j + 1;
			Cell& cell = model->cells[idx1];

			polygon->GetPointIds()->SetId(0, idx);
			polygon->GetPointIds()->SetId(1, idx+ny-1);
			polygon->GetPointIds()->SetId(2, idx+ny);
			polygon->GetPointIds()->SetId(3, idx+1);
			polygons->InsertNextCell(polygon);

			pres->InsertNextValue(cell.u_next.p);
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
/*		idx++;
		polygon->GetPointIds()->SetId(0, idx);
		polygon->GetPointIds()->SetId(1, idx+ny-1);
		polygon->GetPointIds()->SetId(2, idx+ny-1);
		polygon->GetPointIds()->SetId(3, idx);
		grid->InsertNextCell(polygon->GetCellType(), polygon->GetPointIds());*/
	}
	
	// Right border
/*	idx = (nx-2) * (ny-1);
	polygon->GetPointIds()->SetId(0, idx);
	polygon->GetPointIds()->SetId(1, idx);
	polygon->GetPointIds()->SetId(2, idx);
	polygon->GetPointIds()->SetId(3, idx);
	grid->InsertNextCell(polygon->GetCellType(), polygon->GetPointIds());
	for(j = 0; j < ny-2; j++)
	{
		idx = (nx-2) * (ny-1) + j;
		polygon->GetPointIds()->SetId(0, idx);
		polygon->GetPointIds()->SetId(1, idx);
		polygon->GetPointIds()->SetId(2, idx+1);
		polygon->GetPointIds()->SetId(3, idx+1);
		grid->InsertNextCell(polygon->GetCellType(), polygon->GetPointIds());
	}
	idx++;
	polygon->GetPointIds()->SetId(0, idx);
	polygon->GetPointIds()->SetId(1, idx);
	polygon->GetPointIds()->SetId(2, idx);
	polygon->GetPointIds()->SetId(3, idx);
	grid->InsertNextCell(polygon->GetCellType(), polygon->GetPointIds());*/

	grid->SetPolys(polygons);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(pres);
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
void VTKSnapshotter<oil1D::Oil1D>::dump(int i)
{
	using oil1D::Cell;

	// Grid
	vtkSmartPointer<vtkPolyData> grid =
		vtkSmartPointer<vtkPolyData>::New();

	// Points
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();
	for(int k = 1; k < nx; k++)
	{
		Cell& cell = model->cells[k];
		points->InsertNextPoint(r_dim * (cell.r-cell.hr/2.0), 0.0, 0.0);
	}
	for(int k = 1; k < nx; k++)
	{
		Cell& cell = model->cells[k];
		points->InsertNextPoint(r_dim * (cell.r-cell.hr/2.0), -r_dim * model->props_sk.height, 0.0);
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

	/*vtkSmartPointer<vtkDoubleArray> vel_oil = 
		vtkSmartPointer<vtkDoubleArray>::New();
	vel_oil->SetName("oilVelocity");
	vel_oil->SetNumberOfComponents(3);*/

	//double vel [3];

	// Middle cells
	for(int k = 0; k < nx-2; k++)
	{
			Cell& cell = model->cells[k+1];

			polygon->GetPointIds()->SetId(0, k);
			polygon->GetPointIds()->SetId(1, k+1);
			polygon->GetPointIds()->SetId(2, k+nx);
			polygon->GetPointIds()->SetId(3, k+nx-1);
			polygons->InsertNextCell(polygon);

			pres->InsertNextValue(cell.u_next.p);
			/*vel[0] = r_dim / t_dim * model->getOilVelocity(cell, NEXT, R_AXIS);
			vel[1] = r_dim / t_dim * model->getOilVelocity(cell, NEXT, Z_AXIS);	
			vel[2] = 0.0;
			vel_oil->InsertNextTuple(vel);*/
	}

	grid->SetPolys(polygons);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(pres);
	//fd->AddArray(vel_oil);

	// Writing
	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	
	writer->SetFileName( getFileName(i).c_str() );
	writer->SetInputData(grid);
	writer->Write();
}

template <>
void VTKSnapshotter<gas1D::Gas1D>::dump(int i)
{
	using gas1D::Cell;

	// Grid
	vtkSmartPointer<vtkPolyData> grid =
		vtkSmartPointer<vtkPolyData>::New();

	// Points
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();
	for(int k = 1; k < nx; k++)
	{
		Cell& cell = model->cells[k];
		points->InsertNextPoint(r_dim * (cell.r-cell.hr/2.0), 0.0, 0.0);
	}
	for(int k = 1; k < nx; k++)
	{
		Cell& cell = model->cells[k];
		points->InsertNextPoint(r_dim * (cell.r-cell.hr/2.0), -r_dim * model->props_sk[0].height, 0.0);
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

	/*vtkSmartPointer<vtkDoubleArray> vel_oil = 
		vtkSmartPointer<vtkDoubleArray>::New();
	vel_oil->SetName("oilVelocity");
	vel_oil->SetNumberOfComponents(3);*/

	//double vel [3];

	// Middle cells
	for(int k = 0; k < nx-2; k++)
	{
			Cell& cell = model->cells[k+1];

			polygon->GetPointIds()->SetId(0, k);
			polygon->GetPointIds()->SetId(1, k+1);
			polygon->GetPointIds()->SetId(2, k+nx);
			polygon->GetPointIds()->SetId(3, k+nx-1);
			polygons->InsertNextCell(polygon);

			pres->InsertNextValue(cell.u_next.p);
	}

	grid->SetPolys(polygons);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(pres);
	//fd->AddArray(vel_oil);

	// Writing
	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	
	writer->SetFileName( getFileName(i).c_str() );
	writer->SetInputData(grid);
	writer->Write();
}

template <>
void VTKSnapshotter<oil1D_NIT::Oil1D_NIT>::dump(int i)
{
	using oil1D_NIT::Cell;

	// Grid
	vtkSmartPointer<vtkPolyData> grid =
		vtkSmartPointer<vtkPolyData>::New();

	// Points
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();
	for(int k = 1; k < nx; k++)
	{
		Cell& cell = model->cells[k];
		points->InsertNextPoint(r_dim * (cell.r-cell.hr/2.0), 0.0, 0.0);
	}
	for(int k = 1; k < nx; k++)
	{
		Cell& cell = model->cells[k];
		points->InsertNextPoint(r_dim * (cell.r-cell.hr/2.0), -r_dim * model->props_sk[0].height, 0.0);
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

	vtkSmartPointer<vtkDoubleArray> temp = 
		vtkSmartPointer<vtkDoubleArray>::New();
	temp->SetName("temperature");

	vtkSmartPointer<vtkDoubleArray> vel_oil = 
		vtkSmartPointer<vtkDoubleArray>::New();
	vel_oil->SetName("oilVelocity");

	// Middle cells
	for(int k = 0; k < nx-2; k++)
	{
			Cell& cell = model->cells[k+1];

			polygon->GetPointIds()->SetId(0, k);
			polygon->GetPointIds()->SetId(1, k+1);
			polygon->GetPointIds()->SetId(2, k+nx);
			polygon->GetPointIds()->SetId(3, k+nx-1);
			polygons->InsertNextCell(polygon);

			temp->InsertNextValue(cell.u_next.t * T_dim);
			pres->InsertNextValue(cell.u_next.p * model->P_dim / BAR_TO_PA);
			vel_oil->InsertNextValue(model->getOilVelocity(cell, NEXT));
	}

	grid->SetPolys(polygons);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(pres);
	fd->AddArray(temp);
	fd->AddArray(vel_oil);

	// Writing
	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	
	writer->SetFileName( getFileName(i).c_str() );
	writer->SetInputData(grid);
	writer->Write();
}

template <>
void VTKSnapshotter<gas1D::Gas1D_simple>::dump(int i)
{
	using gas1D::Cell;

	// Grid
	vtkSmartPointer<vtkPolyData> grid =
		vtkSmartPointer<vtkPolyData>::New();

	// Points
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();
	for(int k = 1; k < nx; k++)
	{
		Cell& cell = model->cells[k];
		points->InsertNextPoint(r_dim * (cell.r-cell.hr/2.0), 0.0, 0.0);
	}
	for(int k = 1; k < nx; k++)
	{
		Cell& cell = model->cells[k];
		points->InsertNextPoint(r_dim * (cell.r-cell.hr/2.0), -r_dim * model->props_sk[0].height, 0.0);
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

	/*vtkSmartPointer<vtkDoubleArray> vel_oil = 
		vtkSmartPointer<vtkDoubleArray>::New();
	vel_oil->SetName("oilVelocity");
	vel_oil->SetNumberOfComponents(3);*/

	//double vel [3];

	// Middle cells
	for(int k = 0; k < nx-2; k++)
	{
			Cell& cell = model->cells[k+1];

			polygon->GetPointIds()->SetId(0, k);
			polygon->GetPointIds()->SetId(1, k+1);
			polygon->GetPointIds()->SetId(2, k+nx);
			polygon->GetPointIds()->SetId(3, k+nx-1);
			polygons->InsertNextCell(polygon);

			pres->InsertNextValue(cell.u_next.p);
			/*vel[0] = r_dim / t_dim * model->getOilVelocity(cell, NEXT, R_AXIS);
			vel[1] = r_dim / t_dim * model->getOilVelocity(cell, NEXT, Z_AXIS);	
			vel[2] = 0.0;
			vel_oil->InsertNextTuple(vel);*/
	}

	grid->SetPolys(polygons);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(pres);
	//fd->AddArray(vel_oil);

	// Writing
	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	
	writer->SetFileName( getFileName(i).c_str() );
	writer->SetInputData(grid);
	writer->Write();
}

template <class modelType>
void VTKSnapshotter<modelType>::dump_all(int i)
{
}

template <>
void VTKSnapshotter<gasOil_rz_NIT::GasOil_RZ_NIT>::dump_all(int i)
{
	using gasOil_rz_NIT::Cell;

	// Grid
	vtkSmartPointer<vtkPolyData> grid =
		vtkSmartPointer<vtkPolyData>::New();

	// Points
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();
	for(int j = 1; j < ny; j++)
	{
		Cell& cell = model->cells[ j ];
		points->InsertNextPoint(r_dim * (0.99 * cell.r), -r_dim * (cell.z-cell.hz/2.0), 0.0);
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

	vtkSmartPointer<vtkDoubleArray> temp = 
		vtkSmartPointer<vtkDoubleArray>::New();
	temp->SetName("temperature");

	vtkSmartPointer<vtkDoubleArray> pres = 
		vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");

	vtkSmartPointer<vtkDoubleArray> sat_oil = 
		vtkSmartPointer<vtkDoubleArray>::New();
	sat_oil->SetName("oilSaturation");

	vtkSmartPointer<vtkDoubleArray> sat_gas = 
		vtkSmartPointer<vtkDoubleArray>::New();
	sat_gas->SetName("gasSaturation");

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

		temp->InsertNextValue(cell.u_next.t * T_dim);
		pres->InsertNextValue(cell.u_next.p);
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

			temp->InsertNextValue(cell.u_next.t * T_dim);
			pres->InsertNextValue(cell.u_next.p);
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
		}
	}
	grid->SetPolys(polygons);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(temp);
	fd->AddArray(pres);
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
		points->InsertNextPoint(r_dim * (0.99 * cell.r), -r_dim * (cell.z-cell.hz/2.0), 0.0);
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

	vtkSmartPointer<vtkDoubleArray> sat_oil = 
		vtkSmartPointer<vtkDoubleArray>::New();
	sat_oil->SetName("oilSaturation");

	vtkSmartPointer<vtkDoubleArray> sat_gas = 
		vtkSmartPointer<vtkDoubleArray>::New();
	sat_gas->SetName("gasSaturation");

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

		pres->InsertNextValue(cell.u_next.p);
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

			pres->InsertNextValue(cell.u_next.p);
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
void VTKSnapshotter<gasOil_3d::GasOil_3D>::dump_all(int i)
{
	using gasOil_3d::Cell;
	
	// Grid
	vtkSmartPointer<vtkUnstructuredGrid> grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	// Points
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();

	for(int m = 0; m < ny; m++)
	{
		for(int j = 1; j < nz; j++)
		{
			Cell& cell = model->cells[ j + m * nx * nz ];
			points->InsertNextPoint(r_dim * (0.99 * cell.r) * cos(cell.phi - cell.hphi / 2.0), r_dim * (0.99 * cell.r) * sin(cell.phi - cell.hphi / 2.0), -r_dim * (cell.z-cell.hz/2.0));
		}
	
		for(int k = 1; k < nx; k++)
		{
			for(int j = 1; j < nz; j++)
			{
				Cell& cell = model->cells[ j + k * nz + m * nx * nz];
				points->InsertNextPoint(r_dim * (cell.r-cell.hr/2.0) * cos(cell.phi - cell.hphi / 2.0), r_dim * (cell.r-cell.hr/2.0) * sin(cell.phi - cell.hphi / 2.0), -r_dim * (cell.z-cell.hz/2.0));
			}
		}
	}
	grid->SetPoints(points);

	// Data
	vtkSmartPointer<vtkDoubleArray> pres = 
		vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");

	vtkSmartPointer<vtkDoubleArray> sat_oil = 
		vtkSmartPointer<vtkDoubleArray>::New();
	sat_oil->SetName("oilSaturation");

	vtkSmartPointer<vtkDoubleArray> sat_gas = 
		vtkSmartPointer<vtkDoubleArray>::New();
	sat_gas->SetName("gasSaturation");

	vtkSmartPointer<vtkDoubleArray> vel_oil = 
		vtkSmartPointer<vtkDoubleArray>::New();
	vel_oil->SetName("oilVelocity");
	vel_oil->SetNumberOfComponents(3);

	vtkSmartPointer<vtkDoubleArray> vel_gas = 
		vtkSmartPointer<vtkDoubleArray>::New();
	vel_gas->SetName("gasVelocity");
	vel_gas->SetNumberOfComponents(3);

	vtkSmartPointer<vtkIntArray> satur =
		vtkSmartPointer<vtkIntArray>::New();
	satur->SetName("SATUR");

	double vel [3];

	vtkSmartPointer<vtkCellArray> hexs = 
		vtkSmartPointer<vtkCellArray>::New();
	
	int l, j, k;
	for(l = 0; l < ny-1; l++)
	{
		for(j = 0; j < nz-2; j++)
		{	
			Cell& cell = model->cells[ l*nx*nz + j + 1 ];

			vtkSmartPointer<vtkHexahedron> hex =
				vtkSmartPointer<vtkHexahedron>::New();
			hex->GetPointIds()->SetId(0, j + l * nx * (nz-1));
			hex->GetPointIds()->SetId(1, j + (l+1) * nx * (nz-1) );
			hex->GetPointIds()->SetId(2, j + (l+1) * nx * (nz-1) + 1);
			hex->GetPointIds()->SetId(3, j + l * nx * (nz-1) + 1);

			hex->GetPointIds()->SetId(4, j + l * nx * (nz-1) + (nz-1) );
			hex->GetPointIds()->SetId(5, j + (l+1) * nx * (nz-1) + (nz-1) );
			hex->GetPointIds()->SetId(6, j + (l+1) * nx * (nz-1) + 1 + (nz-1) );
			hex->GetPointIds()->SetId(7, j + l * nx * (nz-1) + 1 + (nz-1) );

			hexs->InsertNextCell(hex);

			pres->InsertNextValue(cell.u_next.p * P_dim);
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

		// Middle cells
		for(k = 1; k < nx-1; k++)
		{
			for(j = 0; j < nz-2; j++)
			{		
				Cell& cell = model->cells[ l*nx*nz + k*nz + j + 1 ];

				vtkSmartPointer<vtkHexahedron> hex =
					vtkSmartPointer<vtkHexahedron>::New();
				hex->GetPointIds()->SetId(0, j + k*(nz-1) + l * nx * (nz-1));
				hex->GetPointIds()->SetId(1, j + k*(nz-1) + (l+1) * nx * (nz-1) );
				hex->GetPointIds()->SetId(2, j + k*(nz-1) + (l+1) * nx * (nz-1) + 1);
				hex->GetPointIds()->SetId(3, j + k*(nz-1) + l * nx * (nz-1) + 1);

				hex->GetPointIds()->SetId(4, j + (k+1)*(nz-1) + l * nx * (nz-1) );
				hex->GetPointIds()->SetId(5, j + (k+1)*(nz-1) + (l+1) * nx * (nz-1) );
				hex->GetPointIds()->SetId(6, j + (k+1)*(nz-1) + (l+1) * nx * (nz-1) + 1);
				hex->GetPointIds()->SetId(7, j + (k+1)*(nz-1) + l * nx * (nz-1) + 1);

				hexs->InsertNextCell(hex);

				pres->InsertNextValue(cell.u_next.p * P_dim);
				sat_oil->InsertNextValue(cell.u_next.s);
				sat_gas->InsertNextValue(1.0 - cell.u_next.s);
				satur->InsertNextValue(cell.u_next.SATUR);
				vel[0] = r_dim / t_dim * (cos(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) - sin(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
				vel[1] = r_dim / t_dim * (sin(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) + cos(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
				vel[2] = r_dim / t_dim * model->getOilVelocity(cell, NEXT, Z_AXIS);
				vel_oil->InsertNextTuple(vel);
				vel[0] = r_dim / t_dim * (cos(cell.phi) * model->getGasVelocity(cell, NEXT, R_AXIS) - sin(cell.phi) * model->getGasVelocity(cell, NEXT, PHI_AXIS));
				vel[1] = r_dim / t_dim * (sin(cell.phi) * model->getGasVelocity(cell, NEXT, R_AXIS) + cos(cell.phi) * model->getGasVelocity(cell, NEXT, PHI_AXIS));
				vel[2] = r_dim / t_dim * model->getGasVelocity(cell, NEXT, Z_AXIS);
				vel_gas->InsertNextTuple(vel);
			}
		}
	}

	for(j = 0; j < nz-2; j++)
	{	
		Cell& cell = model->cells[ l*nx*nz + j + 1 ];

		vtkSmartPointer<vtkHexahedron> hex =
			vtkSmartPointer<vtkHexahedron>::New();
		hex->GetPointIds()->SetId(0, j + l * nx * (nz-1));
		hex->GetPointIds()->SetId(1, j);
		hex->GetPointIds()->SetId(2, j + 1);
		hex->GetPointIds()->SetId(3, j + l * nx * (nz-1) + 1);

		hex->GetPointIds()->SetId(4, j + l * nx * (nz-1) + (nz-1) );
		hex->GetPointIds()->SetId(5, j + (nz-1) );
		hex->GetPointIds()->SetId(6, j + 1 + (nz-1) );
		hex->GetPointIds()->SetId(7, j + l * nx * (nz-1) + 1 + (nz-1) );

		hexs->InsertNextCell(hex);

		pres->InsertNextValue(cell.u_next.p * P_dim);
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
	
	for(k = 1; k < nx-1; k++)
	{
		for(j = 0; j < nz-2; j++)
		{		
			Cell& cell = model->cells[ l*nx*nz + k*nz + j + 1 ];

			vtkSmartPointer<vtkHexahedron> hex =
				vtkSmartPointer<vtkHexahedron>::New();
			hex->GetPointIds()->SetId(0, j + k*(nz-1) + l * nx * (nz-1));
			hex->GetPointIds()->SetId(1, j + k*(nz-1) );
			hex->GetPointIds()->SetId(2, j + k*(nz-1) + 1);
			hex->GetPointIds()->SetId(3, j + k*(nz-1) + l * nx * (nz-1) + 1);

			hex->GetPointIds()->SetId(4, j + (k+1)*(nz-1) + l * nx * (nz-1) );
			hex->GetPointIds()->SetId(5, j + (k+1)*(nz-1) );
			hex->GetPointIds()->SetId(6, j + (k+1)*(nz-1) + 1 );
			hex->GetPointIds()->SetId(7, j + (k+1)*(nz-1) + l * nx * (nz-1) + 1 );

			hexs->InsertNextCell(hex);

			pres->InsertNextValue(cell.u_next.p * P_dim);
			sat_oil->InsertNextValue(cell.u_next.s);
			sat_gas->InsertNextValue(1.0 - cell.u_next.s);
			satur->InsertNextValue(cell.u_next.SATUR);
			vel[0] = r_dim / t_dim * (cos(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) - sin(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
			vel[1] = r_dim / t_dim * (sin(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) + cos(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
			vel[2] = r_dim / t_dim * model->getOilVelocity(cell, NEXT, Z_AXIS);
			vel_oil->InsertNextTuple(vel);
			vel[0] = r_dim / t_dim * (cos(cell.phi) * model->getGasVelocity(cell, NEXT, R_AXIS) - sin(cell.phi) * model->getGasVelocity(cell, NEXT, PHI_AXIS));
			vel[1] = r_dim / t_dim * (sin(cell.phi) * model->getGasVelocity(cell, NEXT, R_AXIS) + cos(cell.phi) * model->getGasVelocity(cell, NEXT, PHI_AXIS));
			vel[2] = r_dim / t_dim * model->getGasVelocity(cell, NEXT, Z_AXIS);
			vel_gas->InsertNextTuple(vel);
		}
	}


	grid->SetCells(VTK_HEXAHEDRON, hexs);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(pres);
	fd->AddArray(sat_oil);
	fd->AddArray(sat_gas);
	fd->AddArray(satur);
	fd->AddArray(vel_oil);
	fd->AddArray(vel_gas);

	// Writing
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
		vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	
	writer->SetFileName( getFileName(i).c_str() );
	writer->SetInputData(grid);
	writer->Write();
}

template <>
void VTKSnapshotter<gasOil_perf::GasOil_Perf>::dump_all(int i)
{
	using gasOil_perf::Cell;

	// Grid
	vtkSmartPointer<vtkUnstructuredGrid> grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	// Points
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();

	for (int m = 0; m < ny; m++)
	{
		for (int j = 1; j < nz; j++)
		{
			Cell& cell = model->cells[j + m * nx * nz];
			points->InsertNextPoint(r_dim * (0.99 * cell.r) * cos(cell.phi - cell.hphi / 2.0), r_dim * (0.99 * cell.r) * sin(cell.phi - cell.hphi / 2.0), -r_dim * (cell.z - cell.hz / 2.0));
		}

		for (int k = 1; k < nx; k++)
		{
			for (int j = 1; j < nz; j++)
			{
				Cell& cell = model->cells[j + k * nz + m * nx * nz];
				points->InsertNextPoint(r_dim * (cell.r - cell.hr / 2.0) * cos(cell.phi - cell.hphi / 2.0), r_dim * (cell.r - cell.hr / 2.0) * sin(cell.phi - cell.hphi / 2.0), -r_dim * (cell.z - cell.hz / 2.0));
			}
		}
	}

	vector<std::pair<int,int> >::iterator it;
	int counter = 0;
	for (it = model->perfTunnels.begin(); it != model->perfTunnels.end(); ++it)
	{
		counter += (4 * it->second + 1);

		Cell& cell = model->tunnelCells[counter-1];
		points->InsertNextPoint(r_dim * (0.99 * cell.r) * cos(cell.phi - 0.9 * cell.hphi / 2.0), r_dim * (0.99 * cell.r) * sin(cell.phi - 0.9 * cell.hphi / 2.0), -r_dim * (cell.z - 0.9 * cell.hz / 2.0));
		points->InsertNextPoint(r_dim * (0.99 * cell.r) * cos(cell.phi + 0.9 * cell.hphi / 2.0), r_dim * (0.99 * cell.r) * sin(cell.phi + 0.9 * cell.hphi / 2.0), -r_dim * (cell.z - 0.9 * cell.hz / 2.0));
		points->InsertNextPoint(r_dim * (0.99 * cell.r) * cos(cell.phi + 0.9 * cell.hphi / 2.0), r_dim * (0.99 * cell.r) * sin(cell.phi + 0.9 * cell.hphi / 2.0), -r_dim * (cell.z + 0.9 * cell.hz / 2.0));
		points->InsertNextPoint(r_dim * (0.99 * cell.r) * cos(cell.phi - 0.9 * cell.hphi / 2.0), r_dim * (0.99 * cell.r) * sin(cell.phi - 0.9 * cell.hphi / 2.0), -r_dim * (cell.z + 0.9 * cell.hz / 2.0));
		
		for (int row = 0; row < it->second; row++)
		{
			Cell& cell = model->cells[it->first + (it->second-row-1) * nz];
			points->InsertNextPoint(r_dim * (cell.r + cell.hr / 2.0) * cos(cell.phi - 0.9 * cell.hphi / 2.0), r_dim * (cell.r + cell.hr / 2.0) * sin(cell.phi - 0.9 * cell.hphi / 2.0), -r_dim * (cell.z - 0.9 * cell.hz / 2.0));
			points->InsertNextPoint(r_dim * (cell.r + cell.hr / 2.0) * cos(cell.phi + 0.9 * cell.hphi / 2.0), r_dim * (cell.r + cell.hr / 2.0) * sin(cell.phi + 0.9 * cell.hphi / 2.0), -r_dim * (cell.z - 0.9 * cell.hz / 2.0));
			points->InsertNextPoint(r_dim * (cell.r + cell.hr / 2.0) * cos(cell.phi + 0.9 * cell.hphi / 2.0), r_dim * (cell.r + cell.hr / 2.0) * sin(cell.phi + 0.9 * cell.hphi / 2.0), -r_dim * (cell.z + 0.9 * cell.hz / 2.0));
			points->InsertNextPoint(r_dim * (cell.r + cell.hr / 2.0) * cos(cell.phi - 0.9 * cell.hphi / 2.0), r_dim * (cell.r + cell.hr / 2.0) * sin(cell.phi - 0.9 * cell.hphi / 2.0), -r_dim * (cell.z + 0.9 * cell.hz / 2.0));		
		}
	}

	grid->SetPoints(points);

	// Data
	vtkSmartPointer<vtkDoubleArray> pres =
		vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");

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

	double vel[3];

	vtkSmartPointer<vtkCellArray> hexs =
		vtkSmartPointer<vtkCellArray>::New();

	int l, j, k;
	for (l = 0; l < ny - 1; l++)
	{
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

				pres->InsertNextValue(cell.u_next.p * P_dim);
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
				Cell& cell = model->cells[l*nx*nz + k*nz + j + 1];
				if (cell.isUsed)
				{
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
					sat_oil->InsertNextValue(cell.u_next.s);
					sat_gas->InsertNextValue(1.0 - cell.u_next.s);
					satur->InsertNextValue(cell.u_next.SATUR);
					vel[0] = r_dim / t_dim * (cos(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) - sin(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
					vel[1] = r_dim / t_dim * (sin(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) + cos(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
					vel[2] = r_dim / t_dim * model->getOilVelocity(cell, NEXT, Z_AXIS);
					vel_oil->InsertNextTuple(vel);
					vel[0] = r_dim / t_dim * (cos(cell.phi) * model->getGasVelocity(cell, NEXT, R_AXIS) - sin(cell.phi) * model->getGasVelocity(cell, NEXT, PHI_AXIS));
					vel[1] = r_dim / t_dim * (sin(cell.phi) * model->getGasVelocity(cell, NEXT, R_AXIS) + cos(cell.phi) * model->getGasVelocity(cell, NEXT, PHI_AXIS));
					vel[2] = r_dim / t_dim * model->getGasVelocity(cell, NEXT, Z_AXIS);
					vel_gas->InsertNextTuple(vel);
				}
			}
		}
	}

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

			pres->InsertNextValue(cell.u_next.p * P_dim);
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

	for (k = 1; k < nx - 1; k++)
	{
		for (j = 0; j < nz - 2; j++)
		{
			Cell& cell = model->cells[l*nx*nz + k*nz + j + 1];
			if (cell.isUsed)
			{
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
				sat_oil->InsertNextValue(cell.u_next.s);
				sat_gas->InsertNextValue(1.0 - cell.u_next.s);
				satur->InsertNextValue(cell.u_next.SATUR);
				vel[0] = r_dim / t_dim * (cos(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) - sin(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
				vel[1] = r_dim / t_dim * (sin(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) + cos(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
				vel[2] = r_dim / t_dim * model->getOilVelocity(cell, NEXT, Z_AXIS);
				vel_oil->InsertNextTuple(vel);
				vel[0] = r_dim / t_dim * (cos(cell.phi) * model->getGasVelocity(cell, NEXT, R_AXIS) - sin(cell.phi) * model->getGasVelocity(cell, NEXT, PHI_AXIS));
				vel[1] = r_dim / t_dim * (sin(cell.phi) * model->getGasVelocity(cell, NEXT, R_AXIS) + cos(cell.phi) * model->getGasVelocity(cell, NEXT, PHI_AXIS));
				vel[2] = r_dim / t_dim * model->getGasVelocity(cell, NEXT, Z_AXIS);
				vel_gas->InsertNextTuple(vel);
			}
		}
	}

	const int pointsNum = nx * ny * (nz - 1);
	int tunNum = 0;
	int ind;
	counter = 0;
	for (it = model->perfTunnels.begin(); it != model->perfTunnels.end(); ++it)
	{
		ind = it->first + (it->second + 1) * nz;
		l = int(ind / (nz * nx));//phi
		k = int((ind % (nz * nx)) / nz);//r
		j = int((ind % (nz * nx)) % nz)-1;//z

		vtkSmartPointer<vtkHexahedron> hex =
			vtkSmartPointer<vtkHexahedron>::New();	
		
		hex->GetPointIds()->SetId(0, j + k*(nz - 1) + l * nx * (nz - 1));
		hex->GetPointIds()->SetId(1, j + k*(nz - 1) + (l + 1) * nx * (nz - 1));
		hex->GetPointIds()->SetId(2, j + k*(nz - 1) + (l + 1) * nx * (nz - 1) + 1);
		hex->GetPointIds()->SetId(3, j + k*(nz - 1) + l * nx * (nz - 1) + 1);

		hex->GetPointIds()->SetId(4, pointsNum + 4 * (it->second + 1) * tunNum);
		hex->GetPointIds()->SetId(5, pointsNum + 4 * (it->second + 1) * tunNum + 1);
		hex->GetPointIds()->SetId(6, pointsNum + 4 * (it->second + 1) * tunNum + 2);
		hex->GetPointIds()->SetId(7, pointsNum + 4 * (it->second + 1) * tunNum + 3);
		hexs->InsertNextCell(hex);

		counter += 4 * it->second + 1;
		Cell& cell = model->tunnelCells[counter-1];

		pres->InsertNextValue(cell.u_next.p * P_dim);
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

		for (int row = 0; row < it->second; row++)
		{
			ind = it->first + (it->second + 1 - row) * nz;
			l = int(ind / (nz * nx));//phi
			k = int((ind % (nz * nx)) / nz);//r
			j = int((ind % (nz * nx)) % nz) - 1;//z

			// Right
			vtkSmartPointer<vtkHexahedron> hex_right =
				vtkSmartPointer<vtkHexahedron>::New();
			hex_right->GetPointIds()->SetId(0, j + k*(nz - 1) + l * nx * (nz - 1));
			hex_right->GetPointIds()->SetId(1, pointsNum + 4 * (it->second + 1) * tunNum + 4 * row);
			hex_right->GetPointIds()->SetId(2, pointsNum + 4 * (it->second + 1) * tunNum + 4 * row + 3);
			hex_right->GetPointIds()->SetId(3, j + k*(nz - 1) + l * nx * (nz - 1) + 1);
			
			hex_right->GetPointIds()->SetId(4, j + (k-1)*(nz - 1) + l * nx * (nz - 1));
			hex_right->GetPointIds()->SetId(5, pointsNum + 4 * (it->second + 1) * tunNum + 4 * (row + 1));
			hex_right->GetPointIds()->SetId(6, pointsNum + 4 * (it->second + 1) * tunNum + 4 * (row + 1) + 3);
			hex_right->GetPointIds()->SetId(7, j + (k-1)*(nz - 1) + l * nx * (nz - 1) + 1);
			hexs->InsertNextCell(hex_right);

			Cell& cell_right = model->tunnelCells[counter - 4*(row+1) - 1];

			pres->InsertNextValue(cell_right.u_next.p * P_dim);
			sat_oil->InsertNextValue(cell_right.u_next.s);
			sat_gas->InsertNextValue(1.0 - cell_right.u_next.s);
			satur->InsertNextValue(cell.u_next.SATUR);
			vel[0] = 0.0;
			vel[1] = 0.0;
			vel[2] = 0.0;
			vel_oil->InsertNextTuple(vel);
			vel[0] = 0.0;
			vel[1] = 0.0;
			vel[2] = 0.0;
			vel_gas->InsertNextTuple(vel);

			// Top
			vtkSmartPointer<vtkHexahedron> hex_top =
				vtkSmartPointer<vtkHexahedron>::New();
			hex_top->GetPointIds()->SetId(0, pointsNum + 4 * (it->second + 1) * tunNum + 4 * row + 3);
			hex_top->GetPointIds()->SetId(1, pointsNum + 4 * (it->second + 1) * tunNum + 4 * row + 2);
			hex_top->GetPointIds()->SetId(2, j + k*(nz - 1) + (l+1) * nx * (nz - 1) + 1);
			hex_top->GetPointIds()->SetId(3, j + k*(nz - 1) + l * nx * (nz - 1) + 1);
			
			hex_top->GetPointIds()->SetId(4, pointsNum + 4 * (it->second + 1) * tunNum + 4 * (row + 1) + 3);
			hex_top->GetPointIds()->SetId(5, pointsNum + 4 * (it->second + 1) * tunNum + 4 * (row + 1) + 2);
			hex_top->GetPointIds()->SetId(6, j + (k-1)*(nz - 1) + (l + 1) * nx * (nz - 1) + 1);
			hex_top->GetPointIds()->SetId(7, j + (k-1)*(nz - 1) + l * nx * (nz - 1) + 1);
			hexs->InsertNextCell(hex_top);

			Cell& cell_top = model->tunnelCells[counter - 4 * (row + 1)];

			pres->InsertNextValue(cell_top.u_next.p * P_dim);
			sat_oil->InsertNextValue(cell_top.u_next.s);
			sat_gas->InsertNextValue(1.0 - cell_top.u_next.s);
			satur->InsertNextValue(cell.u_next.SATUR);
			vel[0] = 0.0;
			vel[1] = 0.0;
			vel[2] = 0.0;
			vel_oil->InsertNextTuple(vel);
			vel[0] = 0.0;
			vel[1] = 0.0;
			vel[2] = 0.0;
			vel_gas->InsertNextTuple(vel);

			// Left
			vtkSmartPointer<vtkHexahedron> hex_left =
				vtkSmartPointer<vtkHexahedron>::New();
			hex_left->GetPointIds()->SetId(0, pointsNum + 4 * (it->second + 1) * tunNum + 4 * row + 1);
			hex_left->GetPointIds()->SetId(1, j + k*(nz - 1) + (l + 1) * nx * (nz - 1));
			hex_left->GetPointIds()->SetId(2, j + k*(nz - 1) + (l + 1) * nx * (nz - 1) + 1);
			hex_left->GetPointIds()->SetId(3, pointsNum + 4 * (it->second + 1) * tunNum + 4 * row + 2);
			
			hex_left->GetPointIds()->SetId(4, pointsNum + 4 * (it->second + 1) * tunNum + 4 * (row+1) + 1);
			hex_left->GetPointIds()->SetId(5, j + (k-1)*(nz - 1) + (l + 1) * nx * (nz - 1));
			hex_left->GetPointIds()->SetId(6, j + (k-1)*(nz - 1) + (l + 1) * nx * (nz - 1) + 1);
			hex_left->GetPointIds()->SetId(7, pointsNum + 4 * (it->second + 1) * tunNum + 4 * (row+1) + 2);
			hexs->InsertNextCell(hex_left);

			Cell& cell_left = model->tunnelCells[counter - 4 * (row + 1) + 1];

			pres->InsertNextValue(cell_left.u_next.p * P_dim);
			sat_oil->InsertNextValue(cell_left.u_next.s);
			sat_gas->InsertNextValue(1.0 - cell_left.u_next.s);
			satur->InsertNextValue(cell.u_next.SATUR);
			vel[0] = 0.0;
			vel[1] = 0.0;
			vel[2] = 0.0;
			vel_oil->InsertNextTuple(vel);
			vel[0] = 0.0;
			vel[1] = 0.0;
			vel[2] = 0.0;
			vel_gas->InsertNextTuple(vel);

			// Bottom
			vtkSmartPointer<vtkHexahedron> hex_bot =
				vtkSmartPointer<vtkHexahedron>::New();
			hex_bot->GetPointIds()->SetId(0, j + k*(nz - 1) + l * nx * (nz - 1));
			hex_bot->GetPointIds()->SetId(1, j + k*(nz - 1) + (l + 1) * nx * (nz - 1));
			hex_bot->GetPointIds()->SetId(2, pointsNum + 4 * (it->second + 1) * tunNum + 4 * row + 1);
			hex_bot->GetPointIds()->SetId(3, pointsNum + 4 * (it->second + 1) * tunNum + 4 * row);

			hex_bot->GetPointIds()->SetId(4, j + (k - 1)*(nz - 1) + l * nx * (nz - 1));
			hex_bot->GetPointIds()->SetId(5, j + (k - 1)*(nz - 1) + (l + 1) * nx * (nz - 1));
			hex_bot->GetPointIds()->SetId(6, pointsNum + 4 * (it->second + 1) * tunNum + 4 * (row + 1) + 1);
			hex_bot->GetPointIds()->SetId(7, pointsNum + 4 * (it->second + 1) * tunNum + 4 * (row + 1));
			hexs->InsertNextCell(hex_bot);

			Cell& cell_bot = model->tunnelCells[counter - 4 * (row + 1) + 2];

			pres->InsertNextValue(cell_bot.u_next.p * P_dim);
			sat_oil->InsertNextValue(cell_bot.u_next.s);
			sat_gas->InsertNextValue(1.0 - cell_bot.u_next.s);
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

		tunNum++;
	}

	grid->SetCells(VTK_HEXAHEDRON, hexs);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(pres);
	fd->AddArray(sat_oil);
	fd->AddArray(sat_gas);
	fd->AddArray(satur);
	fd->AddArray(vel_oil);
	fd->AddArray(vel_gas);

	// Writing
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
		vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

	writer->SetFileName(getFileName(i).c_str());
	writer->SetInputData(grid);
	writer->Write();
}

template <>
void VTKSnapshotter<gasOil_perf_nit::GasOil_Perf_NIT>::dump_all(int i)
{
	using gasOil_perf_nit::Cell;

	// Grid
	vtkSmartPointer<vtkUnstructuredGrid> grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	// Points
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();

	for (int m = 0; m < ny; m++)
	{
		for (int j = 1; j < nz; j++)
		{
			Cell& cell = model->cells[j + m * nx * nz];
			points->InsertNextPoint(r_dim * (0.99 * cell.r) * cos(cell.phi - cell.hphi / 2.0), r_dim * (0.99 * cell.r) * sin(cell.phi - cell.hphi / 2.0), -r_dim * (cell.z - cell.hz / 2.0));
		}

		for (int k = 1; k < nx; k++)
		{
			for (int j = 1; j < nz; j++)
			{
				Cell& cell = model->cells[j + k * nz + m * nx * nz];
				points->InsertNextPoint(r_dim * (cell.r - cell.hr / 2.0) * cos(cell.phi - cell.hphi / 2.0), r_dim * (cell.r - cell.hr / 2.0) * sin(cell.phi - cell.hphi / 2.0), -r_dim * (cell.z - cell.hz / 2.0));
			}
		}
	}

	vector<std::pair<int, int> >::iterator it;
	int counter = 0;
	for (it = model->perfTunnels.begin(); it != model->perfTunnels.end(); ++it)
	{
		counter += (4 * it->second + 1);

		Cell& cell = model->tunnelCells[counter - 1];
		points->InsertNextPoint(r_dim * (0.99 * cell.r) * cos(cell.phi - 0.9 * cell.hphi / 2.0), r_dim * (0.99 * cell.r) * sin(cell.phi - 0.9 * cell.hphi / 2.0), -r_dim * (cell.z - 0.9 * cell.hz / 2.0));
		points->InsertNextPoint(r_dim * (0.99 * cell.r) * cos(cell.phi + 0.9 * cell.hphi / 2.0), r_dim * (0.99 * cell.r) * sin(cell.phi + 0.9 * cell.hphi / 2.0), -r_dim * (cell.z - 0.9 * cell.hz / 2.0));
		points->InsertNextPoint(r_dim * (0.99 * cell.r) * cos(cell.phi + 0.9 * cell.hphi / 2.0), r_dim * (0.99 * cell.r) * sin(cell.phi + 0.9 * cell.hphi / 2.0), -r_dim * (cell.z + 0.9 * cell.hz / 2.0));
		points->InsertNextPoint(r_dim * (0.99 * cell.r) * cos(cell.phi - 0.9 * cell.hphi / 2.0), r_dim * (0.99 * cell.r) * sin(cell.phi - 0.9 * cell.hphi / 2.0), -r_dim * (cell.z + 0.9 * cell.hz / 2.0));

		for (int row = 0; row < it->second; row++)
		{
			Cell& cell = model->cells[it->first + (it->second - row - 1) * nz];
			points->InsertNextPoint(r_dim * (cell.r + cell.hr / 2.0) * cos(cell.phi - 0.9 * cell.hphi / 2.0), r_dim * (cell.r + cell.hr / 2.0) * sin(cell.phi - 0.9 * cell.hphi / 2.0), -r_dim * (cell.z - 0.9 * cell.hz / 2.0));
			points->InsertNextPoint(r_dim * (cell.r + cell.hr / 2.0) * cos(cell.phi + 0.9 * cell.hphi / 2.0), r_dim * (cell.r + cell.hr / 2.0) * sin(cell.phi + 0.9 * cell.hphi / 2.0), -r_dim * (cell.z - 0.9 * cell.hz / 2.0));
			points->InsertNextPoint(r_dim * (cell.r + cell.hr / 2.0) * cos(cell.phi + 0.9 * cell.hphi / 2.0), r_dim * (cell.r + cell.hr / 2.0) * sin(cell.phi + 0.9 * cell.hphi / 2.0), -r_dim * (cell.z + 0.9 * cell.hz / 2.0));
			points->InsertNextPoint(r_dim * (cell.r + cell.hr / 2.0) * cos(cell.phi - 0.9 * cell.hphi / 2.0), r_dim * (cell.r + cell.hr / 2.0) * sin(cell.phi - 0.9 * cell.hphi / 2.0), -r_dim * (cell.z + 0.9 * cell.hz / 2.0));
		}
	}

	grid->SetPoints(points);

	// Data
	vtkSmartPointer<vtkDoubleArray> temp =
		vtkSmartPointer<vtkDoubleArray>::New();
	temp->SetName("temperature");

	vtkSmartPointer<vtkDoubleArray> pres =
		vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");

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

	double vel[3];

	vtkSmartPointer<vtkCellArray> hexs =
		vtkSmartPointer<vtkCellArray>::New();

	int l, j, k;
	for (l = 0; l < ny - 1; l++)
	{
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
				Cell& cell = model->cells[l*nx*nz + k*nz + j + 1];
				if (cell.isUsed)
				{
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
					sat_oil->InsertNextValue(cell.u_next.s);
					sat_gas->InsertNextValue(1.0 - cell.u_next.s);
					satur->InsertNextValue(cell.u_next.SATUR);
					vel[0] = r_dim / t_dim * (cos(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) - sin(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
					vel[1] = r_dim / t_dim * (sin(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) + cos(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
					vel[2] = r_dim / t_dim * model->getOilVelocity(cell, NEXT, Z_AXIS);
					vel_oil->InsertNextTuple(vel);
					vel[0] = r_dim / t_dim * (cos(cell.phi) * model->getGasVelocity(cell, NEXT, R_AXIS) - sin(cell.phi) * model->getGasVelocity(cell, NEXT, PHI_AXIS));
					vel[1] = r_dim / t_dim * (sin(cell.phi) * model->getGasVelocity(cell, NEXT, R_AXIS) + cos(cell.phi) * model->getGasVelocity(cell, NEXT, PHI_AXIS));
					vel[2] = r_dim / t_dim * model->getGasVelocity(cell, NEXT, Z_AXIS);
					vel_gas->InsertNextTuple(vel);
				}
			}
		}
	}

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

	for (k = 1; k < nx - 1; k++)
	{
		for (j = 0; j < nz - 2; j++)
		{
			Cell& cell = model->cells[l*nx*nz + k*nz + j + 1];
			if (cell.isUsed)
			{
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

				temp->InsertNextValue(cell.u_next.t * T_dim);
				pres->InsertNextValue(cell.u_next.p * P_dim);
				sat_oil->InsertNextValue(cell.u_next.s);
				sat_gas->InsertNextValue(1.0 - cell.u_next.s);
				satur->InsertNextValue(cell.u_next.SATUR);
				vel[0] = r_dim / t_dim * (cos(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) - sin(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
				vel[1] = r_dim / t_dim * (sin(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) + cos(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
				vel[2] = r_dim / t_dim * model->getOilVelocity(cell, NEXT, Z_AXIS);
				vel_oil->InsertNextTuple(vel);
				vel[0] = r_dim / t_dim * (cos(cell.phi) * model->getGasVelocity(cell, NEXT, R_AXIS) - sin(cell.phi) * model->getGasVelocity(cell, NEXT, PHI_AXIS));
				vel[1] = r_dim / t_dim * (sin(cell.phi) * model->getGasVelocity(cell, NEXT, R_AXIS) + cos(cell.phi) * model->getGasVelocity(cell, NEXT, PHI_AXIS));
				vel[2] = r_dim / t_dim * model->getGasVelocity(cell, NEXT, Z_AXIS);
				vel_gas->InsertNextTuple(vel);
			}
		}
	}

	const int pointsNum = nx * ny * (nz - 1);
	int tunNum = 0;
	int ind;
	counter = 0;
	for (it = model->perfTunnels.begin(); it != model->perfTunnels.end(); ++it)
	{
		ind = it->first + (it->second + 1) * nz;
		l = int(ind / (nz * nx));//phi
		k = int((ind % (nz * nx)) / nz);//r
		j = int((ind % (nz * nx)) % nz) - 1;//z

		vtkSmartPointer<vtkHexahedron> hex =
			vtkSmartPointer<vtkHexahedron>::New();

		hex->GetPointIds()->SetId(0, j + k*(nz - 1) + l * nx * (nz - 1));
		hex->GetPointIds()->SetId(1, j + k*(nz - 1) + (l + 1) * nx * (nz - 1));
		hex->GetPointIds()->SetId(2, j + k*(nz - 1) + (l + 1) * nx * (nz - 1) + 1);
		hex->GetPointIds()->SetId(3, j + k*(nz - 1) + l * nx * (nz - 1) + 1);

		hex->GetPointIds()->SetId(4, pointsNum + 4 * (it->second + 1) * tunNum);
		hex->GetPointIds()->SetId(5, pointsNum + 4 * (it->second + 1) * tunNum + 1);
		hex->GetPointIds()->SetId(6, pointsNum + 4 * (it->second + 1) * tunNum + 2);
		hex->GetPointIds()->SetId(7, pointsNum + 4 * (it->second + 1) * tunNum + 3);
		hexs->InsertNextCell(hex);

		counter += 4 * it->second + 1;
		Cell& cell = model->tunnelCells[counter - 1];

		temp->InsertNextValue(cell.u_next.t * T_dim);
		pres->InsertNextValue(cell.u_next.p * P_dim);
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

		for (int row = 0; row < it->second; row++)
		{
			ind = it->first + (it->second + 1 - row) * nz;
			l = int(ind / (nz * nx));//phi
			k = int((ind % (nz * nx)) / nz);//r
			j = int((ind % (nz * nx)) % nz) - 1;//z

												// Right
			vtkSmartPointer<vtkHexahedron> hex_right =
				vtkSmartPointer<vtkHexahedron>::New();
			hex_right->GetPointIds()->SetId(0, j + k*(nz - 1) + l * nx * (nz - 1));
			hex_right->GetPointIds()->SetId(1, pointsNum + 4 * (it->second + 1) * tunNum + 4 * row);
			hex_right->GetPointIds()->SetId(2, pointsNum + 4 * (it->second + 1) * tunNum + 4 * row + 3);
			hex_right->GetPointIds()->SetId(3, j + k*(nz - 1) + l * nx * (nz - 1) + 1);

			hex_right->GetPointIds()->SetId(4, j + (k - 1)*(nz - 1) + l * nx * (nz - 1));
			hex_right->GetPointIds()->SetId(5, pointsNum + 4 * (it->second + 1) * tunNum + 4 * (row + 1));
			hex_right->GetPointIds()->SetId(6, pointsNum + 4 * (it->second + 1) * tunNum + 4 * (row + 1) + 3);
			hex_right->GetPointIds()->SetId(7, j + (k - 1)*(nz - 1) + l * nx * (nz - 1) + 1);
			hexs->InsertNextCell(hex_right);

			Cell& cell_right = model->tunnelCells[counter - 4 * (row + 1) - 1];

			temp->InsertNextValue(cell_right.u_next.t * T_dim);
			pres->InsertNextValue(cell_right.u_next.p * P_dim);
			sat_oil->InsertNextValue(cell_right.u_next.s);
			sat_gas->InsertNextValue(1.0 - cell_right.u_next.s);
			satur->InsertNextValue(cell.u_next.SATUR);
			vel[0] = 0.0;
			vel[1] = 0.0;
			vel[2] = 0.0;
			vel_oil->InsertNextTuple(vel);
			vel[0] = 0.0;
			vel[1] = 0.0;
			vel[2] = 0.0;
			vel_gas->InsertNextTuple(vel);

			// Top
			vtkSmartPointer<vtkHexahedron> hex_top =
				vtkSmartPointer<vtkHexahedron>::New();
			hex_top->GetPointIds()->SetId(0, pointsNum + 4 * (it->second + 1) * tunNum + 4 * row + 3);
			hex_top->GetPointIds()->SetId(1, pointsNum + 4 * (it->second + 1) * tunNum + 4 * row + 2);
			hex_top->GetPointIds()->SetId(2, j + k*(nz - 1) + (l + 1) * nx * (nz - 1) + 1);
			hex_top->GetPointIds()->SetId(3, j + k*(nz - 1) + l * nx * (nz - 1) + 1);

			hex_top->GetPointIds()->SetId(4, pointsNum + 4 * (it->second + 1) * tunNum + 4 * (row + 1) + 3);
			hex_top->GetPointIds()->SetId(5, pointsNum + 4 * (it->second + 1) * tunNum + 4 * (row + 1) + 2);
			hex_top->GetPointIds()->SetId(6, j + (k - 1)*(nz - 1) + (l + 1) * nx * (nz - 1) + 1);
			hex_top->GetPointIds()->SetId(7, j + (k - 1)*(nz - 1) + l * nx * (nz - 1) + 1);
			hexs->InsertNextCell(hex_top);

			Cell& cell_top = model->tunnelCells[counter - 4 * (row + 1)];

			temp->InsertNextValue(cell_top.u_next.t * T_dim);
			pres->InsertNextValue(cell_top.u_next.p * P_dim);
			sat_oil->InsertNextValue(cell_top.u_next.s);
			sat_gas->InsertNextValue(1.0 - cell_top.u_next.s);
			satur->InsertNextValue(cell.u_next.SATUR);
			vel[0] = 0.0;
			vel[1] = 0.0;
			vel[2] = 0.0;
			vel_oil->InsertNextTuple(vel);
			vel[0] = 0.0;
			vel[1] = 0.0;
			vel[2] = 0.0;
			vel_gas->InsertNextTuple(vel);

			// Left
			vtkSmartPointer<vtkHexahedron> hex_left =
				vtkSmartPointer<vtkHexahedron>::New();
			hex_left->GetPointIds()->SetId(0, pointsNum + 4 * (it->second + 1) * tunNum + 4 * row + 1);
			hex_left->GetPointIds()->SetId(1, j + k*(nz - 1) + (l + 1) * nx * (nz - 1));
			hex_left->GetPointIds()->SetId(2, j + k*(nz - 1) + (l + 1) * nx * (nz - 1) + 1);
			hex_left->GetPointIds()->SetId(3, pointsNum + 4 * (it->second + 1) * tunNum + 4 * row + 2);

			hex_left->GetPointIds()->SetId(4, pointsNum + 4 * (it->second + 1) * tunNum + 4 * (row + 1) + 1);
			hex_left->GetPointIds()->SetId(5, j + (k - 1)*(nz - 1) + (l + 1) * nx * (nz - 1));
			hex_left->GetPointIds()->SetId(6, j + (k - 1)*(nz - 1) + (l + 1) * nx * (nz - 1) + 1);
			hex_left->GetPointIds()->SetId(7, pointsNum + 4 * (it->second + 1) * tunNum + 4 * (row + 1) + 2);
			hexs->InsertNextCell(hex_left);

			Cell& cell_left = model->tunnelCells[counter - 4 * (row + 1) + 1];

			temp->InsertNextValue(cell_left.u_next.t * T_dim);
			pres->InsertNextValue(cell_left.u_next.p * P_dim);
			sat_oil->InsertNextValue(cell_left.u_next.s);
			sat_gas->InsertNextValue(1.0 - cell_left.u_next.s);
			satur->InsertNextValue(cell.u_next.SATUR);
			vel[0] = 0.0;
			vel[1] = 0.0;
			vel[2] = 0.0;
			vel_oil->InsertNextTuple(vel);
			vel[0] = 0.0;
			vel[1] = 0.0;
			vel[2] = 0.0;
			vel_gas->InsertNextTuple(vel);

			// Bottom
			vtkSmartPointer<vtkHexahedron> hex_bot =
				vtkSmartPointer<vtkHexahedron>::New();
			hex_bot->GetPointIds()->SetId(0, j + k*(nz - 1) + l * nx * (nz - 1));
			hex_bot->GetPointIds()->SetId(1, j + k*(nz - 1) + (l + 1) * nx * (nz - 1));
			hex_bot->GetPointIds()->SetId(2, pointsNum + 4 * (it->second + 1) * tunNum + 4 * row + 1);
			hex_bot->GetPointIds()->SetId(3, pointsNum + 4 * (it->second + 1) * tunNum + 4 * row);

			hex_bot->GetPointIds()->SetId(4, j + (k - 1)*(nz - 1) + l * nx * (nz - 1));
			hex_bot->GetPointIds()->SetId(5, j + (k - 1)*(nz - 1) + (l + 1) * nx * (nz - 1));
			hex_bot->GetPointIds()->SetId(6, pointsNum + 4 * (it->second + 1) * tunNum + 4 * (row + 1) + 1);
			hex_bot->GetPointIds()->SetId(7, pointsNum + 4 * (it->second + 1) * tunNum + 4 * (row + 1));
			hexs->InsertNextCell(hex_bot);

			Cell& cell_bot = model->tunnelCells[counter - 4 * (row + 1) + 2];

			temp->InsertNextValue(cell_bot.u_next.t * T_dim);
			pres->InsertNextValue(cell_bot.u_next.p * P_dim);
			sat_oil->InsertNextValue(cell_bot.u_next.s);
			sat_gas->InsertNextValue(1.0 - cell_bot.u_next.s);
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

		tunNum++;
	}

	grid->SetCells(VTK_HEXAHEDRON, hexs);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(temp);
	fd->AddArray(pres);
	fd->AddArray(sat_oil);
	fd->AddArray(sat_gas);
	fd->AddArray(satur);
	fd->AddArray(vel_oil);
	fd->AddArray(vel_gas);

	// Writing
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
		vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

	writer->SetFileName(getFileName(i).c_str());
	writer->SetInputData(grid);
	writer->Write();
}

template <>
void VTKSnapshotter<gasOil_3d_NIT::GasOil_3D_NIT>::dump_all(int i)
{
	using gasOil_3d_NIT::Cell;
	
	// Grid
	vtkSmartPointer<vtkUnstructuredGrid> grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	// Points
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();

	for(int m = 0; m < ny; m++)
	{
		for(int j = 1; j < nz; j++)
		{
			Cell& cell = model->cells[ j + m * nx * nz ];
			points->InsertNextPoint(r_dim * (0.99 * cell.r) * cos(cell.phi - cell.hphi / 2.0), r_dim * (0.99 * cell.r) * sin(cell.phi - cell.hphi / 2.0), -r_dim * (cell.z-cell.hz/2.0));
		}
	
		for(int k = 1; k < nx; k++)
		{
			for(int j = 1; j < nz; j++)
			{
				Cell& cell = model->cells[ j + k * nz + m * nx * nz];
				points->InsertNextPoint(r_dim * (cell.r-cell.hr/2.0) * cos(cell.phi - cell.hphi / 2.0), r_dim * (cell.r-cell.hr/2.0) * sin(cell.phi - cell.hphi / 2.0), -r_dim * (cell.z-cell.hz/2.0));
			}
		}
	}
	grid->SetPoints(points);

	// Data
	vtkSmartPointer<vtkDoubleArray> temp = 
		vtkSmartPointer<vtkDoubleArray>::New();
	temp->SetName("temperature");

	vtkSmartPointer<vtkDoubleArray> pres = 
		vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");

	vtkSmartPointer<vtkDoubleArray> sat_oil = 
		vtkSmartPointer<vtkDoubleArray>::New();
	sat_oil->SetName("oilSaturation");

	vtkSmartPointer<vtkDoubleArray> sat_gas = 
		vtkSmartPointer<vtkDoubleArray>::New();
	sat_gas->SetName("gasSaturation");

	vtkSmartPointer<vtkDoubleArray> vel_oil = 
		vtkSmartPointer<vtkDoubleArray>::New();
	vel_oil->SetName("oilVelocity");
	vel_oil->SetNumberOfComponents(3);

	vtkSmartPointer<vtkDoubleArray> vel_gas = 
		vtkSmartPointer<vtkDoubleArray>::New();
	vel_gas->SetName("gasVelocity");
	vel_gas->SetNumberOfComponents(3);

	double vel [3];

	vtkSmartPointer<vtkCellArray> hexs = 
		vtkSmartPointer<vtkCellArray>::New();
	
	int l, j, k;
	for(l = 0; l < ny-1; l++)
	{
		for(j = 0; j < nz-2; j++)
		{	
			Cell& cell = model->cells[ l*nx*nz + j + 1 ];

			vtkSmartPointer<vtkHexahedron> hex =
				vtkSmartPointer<vtkHexahedron>::New();
			hex->GetPointIds()->SetId(0, j + l * nx * (nz-1));
			hex->GetPointIds()->SetId(1, j + (l+1) * nx * (nz-1) );
			hex->GetPointIds()->SetId(2, j + (l+1) * nx * (nz-1) + 1);
			hex->GetPointIds()->SetId(3, j + l * nx * (nz-1) + 1);

			hex->GetPointIds()->SetId(4, j + l * nx * (nz-1) + (nz-1) );
			hex->GetPointIds()->SetId(5, j + (l+1) * nx * (nz-1) + (nz-1) );
			hex->GetPointIds()->SetId(6, j + (l+1) * nx * (nz-1) + 1 + (nz-1) );
			hex->GetPointIds()->SetId(7, j + l * nx * (nz-1) + 1 + (nz-1) );

			hexs->InsertNextCell(hex);

			temp->InsertNextValue(cell.u_next.t * T_dim);
			pres->InsertNextValue(cell.u_next.p * P_dim);
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
		}

		// Middle cells
		for(k = 1; k < nx-1; k++)
		{
			for(j = 0; j < nz-2; j++)
			{	
				Cell& cell = model->cells[ l*nx*nz + k*nz + j + 1 ];

				vtkSmartPointer<vtkHexahedron> hex =
					vtkSmartPointer<vtkHexahedron>::New();
				hex->GetPointIds()->SetId(0, j + k*(nz-1) + l * nx * (nz-1));
				hex->GetPointIds()->SetId(1, j + k*(nz-1) + (l+1) * nx * (nz-1) );
				hex->GetPointIds()->SetId(2, j + k*(nz-1) + (l+1) * nx * (nz-1) + 1);
				hex->GetPointIds()->SetId(3, j + k*(nz-1) + l * nx * (nz-1) + 1);

				hex->GetPointIds()->SetId(4, j + (k+1)*(nz-1) + l * nx * (nz-1) );
				hex->GetPointIds()->SetId(5, j + (k+1)*(nz-1) + (l+1) * nx * (nz-1) );
				hex->GetPointIds()->SetId(6, j + (k+1)*(nz-1) + (l+1) * nx * (nz-1) + 1);
				hex->GetPointIds()->SetId(7, j + (k+1)*(nz-1) + l * nx * (nz-1) + 1);

				hexs->InsertNextCell(hex);

				temp->InsertNextValue(cell.u_next.t * T_dim);
				pres->InsertNextValue(cell.u_next.p * P_dim);
				sat_oil->InsertNextValue(cell.u_next.s);
				sat_gas->InsertNextValue(1.0 - cell.u_next.s);
				vel[0] = r_dim / t_dim * (cos(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) - sin(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
				vel[1] = r_dim / t_dim * (sin(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) + cos(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
				vel[2] = r_dim / t_dim * model->getOilVelocity(cell, NEXT, Z_AXIS);
				vel_oil->InsertNextTuple(vel);
				vel[0] = r_dim / t_dim * (cos(cell.phi) * model->getGasVelocity(cell, NEXT, R_AXIS) - sin(cell.phi) * model->getGasVelocity(cell, NEXT, PHI_AXIS));
				vel[1] = r_dim / t_dim * (sin(cell.phi) * model->getGasVelocity(cell, NEXT, R_AXIS) + cos(cell.phi) * model->getGasVelocity(cell, NEXT, PHI_AXIS));
				vel[2] = r_dim / t_dim * model->getGasVelocity(cell, NEXT, Z_AXIS);
				vel_gas->InsertNextTuple(vel);
			}
		}
	}

	for(j = 0; j < nz-2; j++)
	{	
		Cell& cell = model->cells[ l*nx*nz + j + 1 ];

		vtkSmartPointer<vtkHexahedron> hex =
			vtkSmartPointer<vtkHexahedron>::New();
		hex->GetPointIds()->SetId(0, j + l * nx * (nz-1));
		hex->GetPointIds()->SetId(1, j);
		hex->GetPointIds()->SetId(2, j + 1);
		hex->GetPointIds()->SetId(3, j + l * nx * (nz-1) + 1);

		hex->GetPointIds()->SetId(4, j + l * nx * (nz-1) + (nz-1) );
		hex->GetPointIds()->SetId(5, j + (nz-1) );
		hex->GetPointIds()->SetId(6, j + 1 + (nz-1) );
		hex->GetPointIds()->SetId(7, j + l * nx * (nz-1) + 1 + (nz-1) );

		hexs->InsertNextCell(hex);

		temp->InsertNextValue(cell.u_next.t * T_dim);
		pres->InsertNextValue(cell.u_next.p * P_dim);
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
	}
	
	for(k = 1; k < nx-1; k++)
	{
		for(j = 0; j < nz-2; j++)
		{		
			Cell& cell = model->cells[ l*nx*nz + k*nz + j + 1 ];

			vtkSmartPointer<vtkHexahedron> hex =
				vtkSmartPointer<vtkHexahedron>::New();
			hex->GetPointIds()->SetId(0, j + k*(nz-1) + l * nx * (nz-1));
			hex->GetPointIds()->SetId(1, j + k*(nz-1) );
			hex->GetPointIds()->SetId(2, j + k*(nz-1) + 1);
			hex->GetPointIds()->SetId(3, j + k*(nz-1) + l * nx * (nz-1) + 1);

			hex->GetPointIds()->SetId(4, j + (k+1)*(nz-1) + l * nx * (nz-1) );
			hex->GetPointIds()->SetId(5, j + (k+1)*(nz-1) );
			hex->GetPointIds()->SetId(6, j + (k+1)*(nz-1) + 1 );
			hex->GetPointIds()->SetId(7, j + (k+1)*(nz-1) + l * nx * (nz-1) + 1 );

			hexs->InsertNextCell(hex);

			temp->InsertNextValue(cell.u_next.t * T_dim);
			pres->InsertNextValue(cell.u_next.p * P_dim);
			sat_oil->InsertNextValue(cell.u_next.s);
			sat_gas->InsertNextValue(1.0 - cell.u_next.s);
			vel[0] = r_dim / t_dim * (cos(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) - sin(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
			vel[1] = r_dim / t_dim * (sin(cell.phi) * model->getOilVelocity(cell, NEXT, R_AXIS) + cos(cell.phi) * model->getOilVelocity(cell, NEXT, PHI_AXIS));
			vel[2] = r_dim / t_dim * model->getOilVelocity(cell, NEXT, Z_AXIS);
			vel_oil->InsertNextTuple(vel);
			vel[0] = r_dim / t_dim * (cos(cell.phi) * model->getGasVelocity(cell, NEXT, R_AXIS) - sin(cell.phi) * model->getGasVelocity(cell, NEXT, PHI_AXIS));
			vel[1] = r_dim / t_dim * (sin(cell.phi) * model->getGasVelocity(cell, NEXT, R_AXIS) + cos(cell.phi) * model->getGasVelocity(cell, NEXT, PHI_AXIS));
			vel[2] = r_dim / t_dim * model->getGasVelocity(cell, NEXT, Z_AXIS);
			vel_gas->InsertNextTuple(vel);
		}
	}


	grid->SetCells(VTK_HEXAHEDRON, hexs);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(temp);
	fd->AddArray(pres);
	fd->AddArray(sat_oil);
	fd->AddArray(sat_gas);
	fd->AddArray(vel_oil);
	fd->AddArray(vel_gas);

	// Writing
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
		vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	
	writer->SetFileName( getFileName(i).c_str() );
	writer->SetInputData(grid);
	writer->Write();
}

template <>
void VTKSnapshotter<oil_rz::Oil_RZ>::dump_all(int i)
{
	using oil_rz::Cell;

	// Grid
	vtkSmartPointer<vtkPolyData> grid =
		vtkSmartPointer<vtkPolyData>::New();

	// Points
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();
	for(int j = 1; j < ny; j++)
	{
		Cell& cell = model->cells[ j ];
		points->InsertNextPoint(r_dim * (0.99 * cell.r), -r_dim * (cell.z-cell.hz/2.0), 0.0);
	}
	for(int k = 1; k < nx; k++)
	{
		for(int j = 1; j < ny; j++)
		{
			Cell& cell = model->cells[ k * ny + j ];
			double ttemp = r_dim * (cell.r - cell.hr / 2.0);
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

	vtkSmartPointer<vtkDoubleArray> vel_oil = 
		vtkSmartPointer<vtkDoubleArray>::New();
	vel_oil->SetName("oilVelocity");
	vel_oil->SetNumberOfComponents(3);

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

		pres->InsertNextValue(cell.u_next.p);
		vel[0] = 0.0;
		vel[1] = 0.0;	
		vel[2] = 0.0;
		vel_oil->InsertNextTuple(vel);
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

			pres->InsertNextValue(cell.u_next.p);
			vel[0] = r_dim / t_dim * model->getOilVelocity(cell, NEXT, R_AXIS);
			vel[1] = -r_dim / t_dim * model->getOilVelocity(cell, NEXT, Z_AXIS);	
			vel[2] = 0.0;
			vel_oil->InsertNextTuple(vel);
		}
	}

	grid->SetPolys(polygons);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(pres);
	fd->AddArray(vel_oil);

	// Writing
	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	
	writer->SetFileName( getFileName(i).c_str() );
	writer->SetInputData(grid);
	writer->Write();
}

template <>
void VTKSnapshotter<oil_rz_nit::Oil_RZ_NIT>::dump_all(int i)
{
	using oil_rz_nit::Cell;

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
			double ttemp = r_dim * (cell.r - cell.hr / 2.0);
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

	vtkSmartPointer<vtkDoubleArray> temp =
		vtkSmartPointer<vtkDoubleArray>::New();
	temp->SetName("temperature");

	vtkSmartPointer<vtkDoubleArray> pres =
		vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");

	vtkSmartPointer<vtkDoubleArray> vel_oil =
		vtkSmartPointer<vtkDoubleArray>::New();
	vel_oil->SetName("oilVelocity");
	vel_oil->SetNumberOfComponents(3);

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

		temp->InsertNextValue(cell.u_next.t * T_dim);
		pres->InsertNextValue(cell.u_next.p);
		vel[0] = 0.0;
		vel[1] = 0.0;
		vel[2] = 0.0;
		vel_oil->InsertNextTuple(vel);
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

			temp->InsertNextValue(cell.u_next.t * T_dim);
			pres->InsertNextValue(cell.u_next.p);
			vel[0] = r_dim / t_dim * model->getOilVelocity(cell, NEXT, R_AXIS);
			vel[1] = r_dim / t_dim * model->getOilVelocity(cell, NEXT, Z_AXIS);
			vel[2] = 0.0;
			vel_oil->InsertNextTuple(vel);
		}
	}

	grid->SetPolys(polygons);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(temp);
	fd->AddArray(pres);
	fd->AddArray(vel_oil);

	// Writing
	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();

	writer->SetFileName(getFileName(i).c_str());
	writer->SetInputData(grid);
	writer->Write();
}



template class VTKSnapshotter<oil1D::Oil1D>;
template class VTKSnapshotter<gas1D::Gas1D>;
template class VTKSnapshotter<gas1D::Gas1D_simple>;
template class VTKSnapshotter<oil1D_NIT::Oil1D_NIT>;
template class VTKSnapshotter<oil_rz::Oil_RZ>;
template class VTKSnapshotter<oil_rz_nit::Oil_RZ_NIT>;
template class VTKSnapshotter<gasOil_rz::GasOil_RZ>;
template class VTKSnapshotter<gasOil_rz_NIT::GasOil_RZ_NIT>;

template class VTKSnapshotter<gasOil_3d::GasOil_3D>;
template class VTKSnapshotter<gasOil_3d_NIT::GasOil_3D_NIT>;

template class VTKSnapshotter<gasOil_perf::GasOil_Perf>;
template class VTKSnapshotter<gasOil_perf_nit::GasOil_Perf_NIT>;