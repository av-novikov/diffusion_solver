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

#include "snapshotter/VTKSnapshotter.h"

#include "model/Oil1D/Oil1D.h"
#include "model/Oil1D_NIT/Oil1D_NIT.h"
#include "model/Oil_RZ/Oil_RZ.h"
#include "model/GasOil_RZ/GasOil_RZ.h"
#include "model/GasOil_RZ_NIT/GasOil_RZ_NIT.h"

using namespace std;

template <class modelType>
VTKSnapshotter<modelType>::VTKSnapshotter()
{
	pattern = prefix + "snap_%{STEP}.vtp";
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
	for(k = 0; k < nx-2; k++)
	{
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
	for(k = 0; k < nx-2; k++)
	{
		for(j = 0; j < ny-2; j++)
		{
			idx = (k+1) * (ny-1) + j;
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
	for(k = 0; k < nx-2; k++)
	{
		for(j = 0; j < ny-2; j++)
		{
			idx = (k+1) * (ny-1) + j;
			idx1 = (k+1) * ny + j + 1;
			Cell& cell = model->cells[idx1];

			polygon->GetPointIds()->SetId(0, idx);
			polygon->GetPointIds()->SetId(1, idx+ny-1);
			polygon->GetPointIds()->SetId(2, idx+ny);
			polygon->GetPointIds()->SetId(3, idx+1);
			polygons->InsertNextCell(polygon);

			pres->InsertNextValue(cell.u_next.p);
			vel[0] = r_dim / t_dim * model->getOilVelocity(cell, NEXT, R_AXIS);
			vel[1] = r_dim / t_dim * model->getOilVelocity(cell, NEXT, Z_AXIS);	
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

template class VTKSnapshotter<oil1D::Oil1D>;
template class VTKSnapshotter<oil1D_NIT::Oil1D_NIT>;
template class VTKSnapshotter<oil_rz::Oil_RZ>;
template class VTKSnapshotter<gasOil_rz::GasOil_RZ>;
template class VTKSnapshotter<gasOil_rz_NIT::GasOil_RZ_NIT>;