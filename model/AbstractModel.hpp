#ifndef ABSTRACTMODEL_HPP_
#define ABSTRACTMODEL_HPP_

#include <vector>
#include <string>
#include <map>

#include "model\cells\stencils\Stencil.h"
#include "snapshotter/VTKSnapshotter.h"
#include "snapshotter/GRDECLSnapshotter.h"
#include "snapshotter/Snapshotter.h"

template <typename varType, typename propsType,
template <typename varType> class cellType, class modelType>
class AbstractModel
{
	template<typename> friend class AbstractSolver;
	template<typename> friend class GRDECLSnapshotter;
	template<typename> friend class VTKSnapshotter;
	template<typename> friend class Snapshotter;
	
	protected:
		varType varInit;
		std::vector<cellType<varType> > cells;
		
		// Spacial properties
		double height_perf;
		double r_w;
		double r_e;
		double Volume;
		int cellsNum;
		
		// Rate of the well
		double Q_sum;
		double Pwf;
		// Ranges of perforated cells numbers
		std::vector<std::pair<int,int> > perfIntervals;
		// Vector of <cell number, rate in the cell> for left border cells
		std::map<int,double> Qcell;

		// Temporary properties
		double ht;
		double ht_min;
		double ht_max;
		
		// Number of periods
		int periodsNum;
		// End times of periods [sec]
		std::vector<double> period;
		// Oil rates [m3/day]
		std::vector<double> rate;
		// Vector of BHPs [bar]
		std::vector<double> pwf;
		// If left boundary condition would be 2nd type
		bool leftBoundIsRate;
		// If right boundary condition would be 1st type
		bool rightBoundIsPres;

		virtual void buildGridLog() = 0;
		virtual void setProps(propsType& props) = 0;
		virtual void makeDimLess() = 0;
		virtual void setPerforated() = 0;
		virtual void setInitialState();
		
		Snapshotter<modelType>* snapshotter;
		bool isWriteSnaps;

	public:
		AbstractModel();
		virtual ~AbstractModel();		
	
		// Dimensions
		double t_dim;
		double R_dim;
		double Z_dim;
		double P_dim;
		double T_dim;
		double Q_dim;

		void setSnapshotter(std::string type, modelType* model);

		void load(propsType& props);
		virtual void setPeriod(int period) = 0;
		virtual void setWellborePeriod(int period, double cut_t);
		int getCellsNum();

		void snapshot(int i);
		void snapshot_all(int i);
};

#endif /* ABSTRACTMODEL_HPP_ */
