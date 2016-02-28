#include "snapshotter/GRDECLSnapshotter.h"
#include "model/AbstractModel.hpp"

#include "model/Oil1D/Oil1D.h"
#include "model/Gas1D/Gas1D.h"
#include "model/Gas1D/Gas1D_simple.h"
#include "model/Oil1D_NIT/Oil1D_NIT.h"
#include "model/Oil_RZ/Oil_RZ.h"
#include "model/GasOil_RZ/GasOil_RZ.h"
#include "model/GasOil_RZ_NIT/GasOil_RZ_NIT.h"

#include "model/3D/GasOil_3D/GasOil_3D.h"
#include "model/3D/GasOil_3D_NIT/GasOil_3D_NIT.h"

#include "model/3D/Perforation/GasOil_Perf.h"

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

template <>
void GRDECLSnapshotter<gasOil_rz_NIT::GasOil_RZ_NIT>::dump(int i)
{
	string gridName = getFileName(i);

	ofstream grid;
	grid.open(gridName.c_str(), ofstream::out);

	int n = model->cellsNum_r * model->cellsNum_z;

	grid << "SPECGRID" << endl;
	grid << "  " << model->cellsNum_r << "  " << model->cellsNum_z << "  " << 1;
	grid << "  " << 1 << "  " << "F /" << endl;

	grid << "COORDSYS" << endl;
	grid << "  " << "1  34  'INCOMP  '  /" << endl;

	grid << "COORD" << endl;
	int idx = 0;
	double R_dim = model->R_dim;
	double z;
	int counter = 0;
	for(int k = 0; k < model->cellsNum_z+1; k++)
	{
		idx = k+1;
		z = R_dim * (model->cells[idx].z - model->cells[idx].hz / 2.0);
		for(int j = 0; j < model->cellsNum_r+1; j++)
		{
			idx = (j+1)*(model->cellsNum_z+2) + (k+1);
			grid << R_dim * (model->cells[idx].r - model->cells[idx].hr / 2.0) << "  " \
					<< z << " " \
					<< 1.0;
			counter++;
			if(counter % 2 == 0)
				grid << endl;
			else
				grid << "  ";
		}
	}

	if(counter % 2 == 0)
		grid << endl;
	else
		grid << endl << endl;

	grid << "ACTNUM" << endl;
	counter = 0;
	while(counter < n)
	{
		grid << counter << "  ";
		if(counter % 6 == 0)
			grid << endl;
		counter++;
	}

	grid << endl;

	// Pressure
	grid << "PRES" << endl;
	counter = 0;
	for(int k = 0; k < model->cellsNum_z; k++)
	{
		for(int j = 0; j < model->cellsNum_r; j++)
		{
			idx = (j+1)*(model->cellsNum_z+2) + (k+1);
			grid << model->cells[idx].u_next.p;
			
			counter++;
			if(counter % 2 == 0)
				grid << endl;
			else
				grid << "  ";
		}
	}

	grid << endl;

	// Saturation
	grid << "SATUR" << endl;
	counter = 0;
	for(int k = 0; k < model->cellsNum_z; k++)
	{
		for(int j = 0; j < model->cellsNum_r; j++)
		{
			idx = (j+1)*(model->cellsNum_z+2) + (k+1);
			grid << model->cells[idx].u_next.s;
			
			counter++;
			if(counter % 2 == 0)
				grid << endl;
			else
				grid << "  ";
		}
	}

	grid.close();
}

template <>
void GRDECLSnapshotter<gasOil_rz_NIT::GasOil_RZ_NIT>::dump_all(int i)
{
	string gridName = getFileName(i);

	ofstream grid;
	grid.open(gridName.c_str(), ofstream::out);

	int n = (model->cellsNum_r + 2) * (model->cellsNum_z + 2);

	grid << "SPECGRID" << endl;
	grid << "  " << model->cellsNum_r+2 << "  " << model->cellsNum_z+2 << "  " << 1;
	grid << "  " << 1 << "  " << "F /" << endl;

	grid << "COORDSYS" << endl;
	grid << "  " << "1  34  'INCOMP  '  /" << endl;

	grid << "COORD" << endl;
	int idx = 0;
	double R_dim = model->R_dim;
	double z;
	int counter = 0;

	idx = 0;
	z = R_dim * (model->cells[idx].z - model->cells[idx].hz / 2.0);
	grid << R_dim * (model->cells[idx].r - model->cells[idx].hr / 2.0) << "  " \
			<< z << " " \
			<< 1.0 << " ";
	counter++;
	for(int k = 0; k < model->cellsNum_r+2; k++)
	{
		idx = k * (model->cellsNum_z+2);
		grid << R_dim * (model->cells[idx].r + model->cells[idx].hr / 2.0) << "  " \
			<< z << " " \
			<< 1.0;
		counter++;

		if(counter % 2 == 0)
			grid << endl;
		else
			grid << "  ";
	}

	for(int k = 0; k < model->cellsNum_z+2; k++)
	{
		idx = k;
		z = R_dim * (model->cells[idx].z + model->cells[idx].hz / 2.0);
		grid << R_dim * (model->cells[idx].r - model->cells[idx].hr / 2.0) << "  " \
			<< z << " " \
			<< 1.0;
		counter++;
		
		if(counter % 2 == 0)
			grid << endl;
		else
			grid << "  ";
		
		for(int j = 0; j < model->cellsNum_r+2; j++)
		{
			idx = j*(model->cellsNum_z+2) + k;
			grid << R_dim * (model->cells[idx].r + model->cells[idx].hr / 2.0) << "  " \
					<< z << " " \
					<< 1.0;
			counter++;

			if(counter % 2 == 0)
				grid << endl;
			else
				grid << "  ";
		}
	}

	grid << endl;

	grid << "ACTNUM" << endl;
	counter = 0;
	while(counter < n)
	{
		grid << counter << "  ";
		if(counter % 6 == 0)
			grid << endl;
		counter++;
	}

	grid << endl;

	// Pressure
	grid << "PRES" << endl;
	counter = 0;
	for(int k = 0; k < model->cellsNum_z+2; k++)
	{
		for(int j = 0; j < model->cellsNum_r+2; j++)
		{
			idx = j * (model->cellsNum_z+2) + k;
			grid << model->cells[idx].u_next.p;
			
			counter++;
			if(counter % 6 == 0)
				grid << endl;
			else
				grid << "  ";
		}
	}

	grid << endl;

	// Saturation
	grid << "SATUR" << endl;
	counter = 0;
	for(int k = 0; k < model->cellsNum_z+2; k++)
	{
		for(int j = 0; j < model->cellsNum_r+2; j++)
		{
			idx = j * (model->cellsNum_z+2) + k;
			grid << model->cells[idx].u_next.s;
			
			counter++;
			if(counter % 6 == 0)
				grid << endl;
			else
				grid << "  ";
		}
	}

	grid.close();
}

template class GRDECLSnapshotter<oil1D::Oil1D>;
template class GRDECLSnapshotter<gas1D::Gas1D>;
template class GRDECLSnapshotter<gas1D::Gas1D_simple>;
template class GRDECLSnapshotter<oil1D_NIT::Oil1D_NIT>;
template class GRDECLSnapshotter<oil_rz::Oil_RZ>;
template class GRDECLSnapshotter<gasOil_rz::GasOil_RZ>;
template class GRDECLSnapshotter<gasOil_rz_NIT::GasOil_RZ_NIT>;

template class GRDECLSnapshotter<gasOil_3d::GasOil_3D>;
template class GRDECLSnapshotter<gasOil_3d_NIT::GasOil_3D_NIT>;

template class GRDECLSnapshotter<gasOil_perf::GasOil_Perf>;