#ifndef BASIC2DSOLVER_HPP_
#define BASIC2DSOLVER_HPP_

#include <iostream>
#include <vector>
#include <map>

#include "model/AbstractSolver.hpp"
#include "method/sweep.h"
#include "model/Basic2d/Basic2d.hpp"

namespace basic2d
{
	template <class modelType>
	class Basic2dSolver : public AbstractSolver<modelType>, public Sweep
	{
	protected:
		void control()
		{
			writeData();

			if (cur_t >= model->period[curTimePeriod])
			{
				curTimePeriod++;
				model->ht = model->ht_min;
				model->setPeriod(curTimePeriod);
			}

			if (model->ht <= model->ht_max && iterations < 6)
				model->ht = model->ht * 1.5;
			else if (iterations > 6 && model->ht > model->ht_min)
				model->ht = model->ht / 1.5;

			if (cur_t + model->ht > model->period[curTimePeriod])
				model->ht = model->period[curTimePeriod] - cur_t;

			cur_t += model->ht;
			model->time = cur_t;
		}
		void doNextStep()
		{
			solveStep();

			if (n > 1 && model->Q_sum > EQUALITY_TOLERANCE)
			{
				double H0 = fabs(model->solveH());
				if (H0 > 0.1)
				{
					printWellRates();
					fillq();

					double mult = 0.9;
					double H = H0;

					while (H > H0 / 50.0 || H > 0.05)
					{
						solveDq(mult);

						int i = 0;
						std::map<int, double>::iterator it = model->Qcell.begin();

						while (it != model->Qcell.end())
						{
							q[i] += mult * dq[i];
							it->second = q[i];
							i++;	++it;
						}

						solveStep();
						printWellRates();

						H = fabs(model->solveH());
					}
				}
			}
		}
		void fillq()
		{
			int i = 0;
			std::map<int, double>::iterator it = model->Qcell.begin();
			while (it != model->Qcell.end())
			{
				q[i++] = it->second;
				++it;
			}
		}
		void fillDq()
		{
			for (int i = 0; i < n; i++)
				dq[i] = 0.0;
		}
		void solveDq(double mult)
		{
			fillDq();
			filldPdQ(mult);
			solveStep();
			solveSystem();

			int i = 0;
			std::map<int, double>::iterator it;
			for (it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
			{
				std::cout << "dq[" << it->first << "] = " << dq[i++] * model->Q_dim * 86400.0 << std::endl;
			}
			std::cout << std::endl;
		}
		void solveSystem()
		{
			double s = 0.0, p1, p2;
			std::map<int, double>::iterator it;

			for (int i = 0; i < n - 1; i++)
			{
				for (int j = 0; j < n - 1; j++)
				{
					s = 0.0;
					for (int k = 0; k < n - 1; k++)
						s += (dpdq[k + 1][j] - dpdq[k][j]) * (dpdq[k + 1][i] - dpdq[k][i]);
					mat[i][j] = s;
				}

				s = 0.0;
				it = model->Qcell.begin();
				for (int k = 0; k < n - 1; k++)
				{
					p1 = model->cells[it->first].u_next.p;
					p2 = model->cells[(++it)->first].u_next.p;
					s += (p2 - p1) * (dpdq[k + 1][i] - dpdq[k][i]);
				}
				b[i] = -s;
			}

			MC_LU rateSystem(mat, b);
			rateSystem.LU_Solve();

			s = 0.0;
			for (int i = 0; i < n - 1; i++)
			{
				dq[i + 1] = rateSystem.ptResult[i];
				s += rateSystem.ptResult[i];
			}
			dq[0] = -s;
		}
		void filldPdQ(double mult)
		{
			double p1, p2, ratio;
			ratio = mult * 0.001 / (double)(n);

			int i = 0, j = 0;
			std::map<int, double>::iterator it0 = model->Qcell.begin();
			std::map<int, double>::iterator it1 = model->Qcell.begin();
			std::map<int, double>::iterator it2 = it1;	++it2;
			while (it1 != model->Qcell.end())
			{
				j = 0;
				it2 = model->Qcell.begin();		++it2;
				while (it2 != model->Qcell.end())
				{
					model->setRateDeviation(it2->first, -ratio);
					model->setRateDeviation(it0->first, ratio);
					solveStep();
					p1 = model->cells[it1->first].u_next.p;

					model->setRateDeviation(it2->first, 2.0 * ratio);
					model->setRateDeviation(it0->first, -2.0 * ratio);
					solveStep();
					p2 = model->cells[it1->first].u_next.p;

					model->setRateDeviation(it2->first, -ratio);
					model->setRateDeviation(it0->first, ratio);

					dpdq[i][j++] = (p2 - p1) / (2.0 * ratio * model->Q_sum);

					++it2;
				}
				i++;
				++it1;
			}
		}

		int n;
		bool isChange;
		TVector<double> dq;
		TVector<double> q;
		TMatrix<double> dpdq;
		TMatrix<double> mat;
		TVector<double> b;

		inline void printWellRates()
		{
			double DQ = model->Q_sum;
			std::map<int, double>::iterator it;
			int k = 0;
			for (it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
			{
				std::cout << "Rate in " << it->first << " = " << it->second * model->Q_dim * 86400.0 << "\t";
				std::cout << "Press in " << it->first << " = " << model->cells[it->first].u_next.p << std::endl;
				DQ -= it->second;
				k++;
			}
			std::cout << "Summary rate deviation = " << DQ * model->Q_dim * 86400.0 << std::endl;
			std::cout << std::endl;
		};
	public:
		Basic2dSolver(modelType* _model) : AbstractSolver<modelType>(_model) 
		{
			Initialize(model->cellsNum_r + 2, model->var_size * (model->cellsNum_z + 2));

			n = model->Qcell.size();
			dq.Initialize(n);
			q.Initialize(n);

			dpdq.Initialize(n, n - 1);
			mat.Initialize(n - 1, n - 1);
			b.Initialize(n - 1);

			Tt = model->period[model->period.size() - 1];
		};
		~Basic2dSolver() {};
	};

};

#endif /* BASIC2DSOLVER_HPP_ */
