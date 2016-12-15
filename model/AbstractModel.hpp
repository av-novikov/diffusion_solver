#ifndef ABSTRACTMODEL_HPP_
#define ABSTRACTMODEL_HPP_

#include <vector>
#include <string>
#include <map>

#include "snapshotter/VTKSnapshotter.h"
#include "snapshotter/GRDECLSnapshotter.h"
#include "snapshotter/Snapshotter.h"
#include "util/utils.h"

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
		// BHP will be converted to the depth
		double depth_point;
		// During the time flow rate decreases 'e' times in well test [sec] 
		double alpha;

		virtual void buildGridLog() = 0;
		virtual void setProps(propsType& props) = 0;
		virtual void makeDimLess() = 0;
		virtual void setPerforated() = 0;
		virtual void setInitialState()
		{
			std::vector<cellType<varType> >::iterator it;
			for (it = cells.begin(); it != cells.end(); ++it)
				it->u_prev = it->u_iter = it->u_next = varInit;
		};
		
		Snapshotter<modelType>* snapshotter;
		bool isWriteSnaps;

	public:
		AbstractModel() { isWriteSnaps = true; };
		virtual ~AbstractModel() {};
	
		// Dimensions
		double t_dim;
		double R_dim;
		double Z_dim;
		double P_dim;
		double T_dim;
		double Q_dim;

		void setSnapshotter(std::string type, modelType* model)
		{
			if (type == "VTK") {
				snapshotter = new VTKSnapshotter<modelType>();
				snapshotter->setModel(model);
				isWriteSnaps = true;
			}
			else if (type == "GRDECL") {
				snapshotter = new GRDECLSnapshotter<modelType>();
				snapshotter->setModel(model);
				isWriteSnaps = true;
			}
			else {
				isWriteSnaps = false;
			}
		};

		void load(propsType& props)
		{
			setProps(props);

			buildGridLog();
			setPerforated();
			setInitialState();
		};
		virtual void setPeriod(int period) = 0;
		virtual void setWellborePeriod(int period, double cut_t) {};
		int getCellsNum() {	return cellsNum; };

		void snapshot(int i) { snapshotter->dump(i); };
		void snapshot_all(int i) { snapshotter->dump_all(i); }
};

#endif /* ABSTRACTMODEL_HPP_ */
