#include "model/3D/GasOil_3D/Par3DSolver.h"

#include <algorithm>

using namespace std;
using namespace gasOil_3d;

Par3DSolver::Par3DSolver(GasOil_3D* _model) : AbstractSolver<GasOil_3D>(_model)
{
	// Output streams
	plot_Pdyn.open("snaps/P_dyn.dat", ofstream::out);
	plot_Sdyn.open("snaps/S_dyn.dat", ofstream::out);
	plot_qcells.open("snaps/q_cells.dat", ofstream::out);

	// Flow rate optimization structures
	n = model->Qcell.size();
	dq.Initialize(n);
	q.Initialize(n);
	dpdq.Initialize(n, n - 1);
	mat.Initialize(n - 1, n - 1);
	b.Initialize(n - 1);

	// Time settings
	t_dim = model->t_dim;
	Tt = model->period[model->period.size() - 1];

	// Memory allocating
	ind_i = new int[7 * 4 * model->cellsNum];
	ind_j = new int[7 * 4 * model->cellsNum];
	a = new double[7 * 4 * model->cellsNum];
	ind_rhs = new int[2 * model->cellsNum];
	rhs = new double[2 * model->cellsNum];

	// Stencils allocating
	stencils = new UsedStencils<GasOil_3D>(model);
	stencils->setStorages(a, ind_i, ind_j, rhs);
}

Par3DSolver::~Par3DSolver()
{
	// Closing streams
	plot_Pdyn.close();
	plot_Sdyn.close();
	plot_qcells.close();
}

void Par3DSolver::writeData()
{
	double p = 0.0, s = 0.0;

	plot_qcells << cur_t * t_dim / 3600.0;

	map<int, double>::iterator it;
	for (it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
	{
		p += model->cells[it->first].u_next.p * model->P_dim;
		s += model->cells[it->first].u_next.s;
		if (model->leftBoundIsRate)
			plot_qcells << "\t" << it->second * model->Q_dim * 86400.0;
		else
			plot_qcells << "\t" << model->getRate(it->first) * model->Q_dim * 86400.0;
	}

	plot_Pdyn << cur_t * t_dim / 3600.0 << "\t" << p / (double)(model->Qcell.size()) << endl;
	plot_Sdyn << cur_t * t_dim / 3600.0 << "\t" << s / (double)(model->Qcell.size()) << endl;

	plot_qcells << endl;
}

void Par3DSolver::control()
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
}

void Par3DSolver::start()
{
	int counter = 0;
	iterations = 8;

	fillIndices();
	solver.Init(2 * model->cellsNum);

	model->setPeriod(curTimePeriod);
	while (cur_t < Tt)
	{
		control();
		if (model->isWriteSnaps)
			model->snapshot_all(counter++);
		doNextStep();
		copyTimeLayer();
		cout << "---------------------NEW TIME STEP---------------------" << endl;
	}
	if (model->isWriteSnaps)
		model->snapshot_all(counter++);
	writeData();
}

void Par3DSolver::doNextStep()
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
				map<int, double>::iterator it = model->Qcell.begin();

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

void Par3DSolver::copySolution(const paralution::LocalVector<double>& sol)
{
	for (int i = 0; i < model->cellsNum; i++)
	{
		model->cells[i].u_next.p += sol[2 * i];
		model->cells[i].u_next.s += sol[2 * i + 1];
	}
}

void Par3DSolver::solveStep()
{
	int cellIdx, varIdx;
	double err_newton = 1.0;
	double averPresPrev = averValue(1);
	double averSatPrev = averValue(2);
	double averPres, averSat;
	double dAverPres = 1.0, dAverSat = 1.0;

	iterations = 0;
	while (err_newton > 1.e-4 && (dAverSat > 1.e-8 || dAverPres > 1.e-4) && iterations < 8)
	{
		copyIterLayer();

		fill();
		solver.Assemble(ind_i, ind_j, a, elemNum, ind_rhs, rhs);
		solver.Solve();
		copySolution( solver.getSolution() );

		model->solveP_bub();

		err_newton = convergance(cellIdx, varIdx);

		averPres = averValue(1);					averSat = averValue(2);
		dAverPres = fabs(averPres - averPresPrev);	dAverSat = fabs(averSat - averSatPrev);
		averPresPrev = averPres;					averSatPrev = averSat;

		//		if(varIdx == PRES)
		//			cout << "BadPresValue[" << cellIdx  << "]: " << model->cells[cellIdx].u_next.p << endl;
		//		else if(varIdx == SAT)
		//			cout << "BadSatValue[" << cellIdx  << "]: " << model->cells[cellIdx].u_next.s << endl;

		iterations++;
	}

	cout << "Newton Iterations = " << iterations << endl;
}

void Par3DSolver::fillIndices()
{
	int idx, nebr;
	int counter = 0;
	Iterator it;
	map<int, double>::iterator itPerf;

	// Left
	for (it = model->getLeftBegin(); it != model->getLeftEnd(); ++it)
	{
		idx = it.getIdx();
		itPerf = model->Qcell.find(idx);
		if (itPerf == model->Qcell.end())
		{
			nebr = idx + model->cellsNum_z + 2;
			ind_i[counter] = 2 * idx;
			ind_j[counter++] = 2 * idx;

			ind_i[counter] = 2 * idx;
			ind_j[counter++] = 2 * idx + 1;

			ind_i[counter] = 2 * idx;
			ind_j[counter++] = 2 * nebr;

			ind_i[counter] = 2 * idx;
			ind_j[counter++] = 2 * nebr + 1;

			ind_i[counter] = 2 * idx + 1;
			ind_j[counter++] = 2 * idx;

			ind_i[counter] = 2 * idx + 1;
			ind_j[counter++] = 2 * idx + 1;

			ind_i[counter] = 2 * idx + 1;
			ind_j[counter++] = 2 * nebr;

			ind_i[counter] = 2 * idx + 1;
			ind_j[counter++] = 2 * nebr + 1;
		}
		else
			stencils->left->fillIndex(idx, &counter);
	}

	// Middle
	int res;
	for (it = model->getMidBegin(); it != model->getMidEnd(); ++it)
	{
		idx = it.getIdx();
		res = idx % (model->cellsNum_z + 2);
		if( res == 0 )
			stencils->top->fillIndex(idx, &counter);
		else if( res == model->cellsNum_z + 1)
			stencils->bot->fillIndex(idx, &counter);
		else
			stencils->middle->fillIndex(idx, &counter);	
	}

	// Right
	for (it = model->getRightBegin(); it != model->getRightEnd(); ++it)
	{
		idx = it.getIdx();
		stencils->right->fillIndex(idx, &counter);
	}

	/*// Left
	for (int i = 0; i < model->cellsNum_phi; i++)
	{
		idx = i * (model->cellsNum_r + 2) * (model->cellsNum_z + 2);
		for (int j = 0; j < model->cellsNum_z + 2; j++)
		{
			it = model->Qcell.find(idx);
			if (it == model->Qcell.end())
			{
				nebr = idx + model->cellsNum_z + 2;
				ind_i[counter] = 2 * idx;
				ind_j[counter++] = 2 * idx;

				ind_i[counter] = 2 * idx;
				ind_j[counter++] = 2 * idx + 1;

				ind_i[counter] = 2 * idx;
				ind_j[counter++] = 2 * nebr;

				ind_i[counter] = 2 * idx;
				ind_j[counter++] = 2 * nebr + 1;

				ind_i[counter] = 2 * idx + 1;
				ind_j[counter++] = 2 * idx;

				ind_i[counter] = 2 * idx + 1;
				ind_j[counter++] = 2 * idx + 1;

				ind_i[counter] = 2 * idx + 1;
				ind_j[counter++] = 2 * nebr;

				ind_i[counter] = 2 * idx + 1;
				ind_j[counter++] = 2 * nebr + 1;
			}
			else
				stencils->left->fillIndex(idx, &counter);

			idx++;
		}
	}

	// Middle
	for (int i = 0; i < model->cellsNum_phi; i++)
	{
		idx = i * (model->cellsNum_r + 2) * (model->cellsNum_z + 2) + model->cellsNum_z + 2;
		for (int j = 1; j < model->cellsNum_r + 1; j++)
		{
			stencils->top->fillIndex(idx, &counter);
			idx++;
			for (int k = 1; k < model->cellsNum_z + 1; k++)
			{
				stencils->middle->fillIndex(idx, &counter);
				idx++;
			}
			stencils->bot->fillIndex(idx, &counter);
			idx++;
		}
	}

	// Right
	for (int i = 0; i < model->cellsNum_phi; i++)
	{
		idx = i * (model->cellsNum_r + 2) * (model->cellsNum_z + 2) + (model->cellsNum_r + 1) * (model->cellsNum_z + 2);
		for (int j = 0; j < model->cellsNum_z + 2; j++)
		{
			stencils->right->fillIndex(idx, &counter);
			idx++;
		}
	}*/

	elemNum = counter;

	for (int i = 0; i < 2*model->cellsNum; i++)
		ind_rhs[i] = i;
}

void Par3DSolver::fill()
{
	int idx;
	int counter = 0;
	Iterator it;
	map<int, double>::iterator itPerf;

	// Left
	for (it = model->getLeftBegin(); it != model->getLeftEnd(); ++it)
	{
		idx = it.getIdx();
		itPerf = model->Qcell.find(idx);
		if (itPerf == model->Qcell.end())
		{
			a[counter++] = 1.0;
			a[counter++] = 0.0;
			a[counter++] = -1.0;
			a[counter++] = 0.0;

			a[counter++] = 0.0;
			a[counter++] = 1.0;
			a[counter++] = 0.0;
			a[counter++] = -1.0;

			rhs[2 * idx] = 0.0;
			rhs[2 * idx + 1] = 0.0;
		}
		else
			stencils->left->fill(idx, &counter);
	}

	// Middle
	int res;
	for (it = model->getMidBegin(); it != model->getMidEnd(); ++it)
	{
		idx = it.getIdx();
		res = idx % (model->cellsNum_z + 2);
		if (res == 0)
			stencils->top->fill(idx, &counter);
		else if (res == model->cellsNum_z + 1)
			stencils->bot->fill(idx, &counter);
		else
			stencils->middle->fill(idx, &counter);
	}

	// Right
	for (it = model->getRightBegin(); it != model->getRightEnd(); ++it)
	{
		idx = it.getIdx();
		stencils->right->fill(idx, &counter);
	}

	/*int idx;
	int counter = 0;
	map<int, double>::iterator it;

	// Left
	for (int i = 0; i < model->cellsNum_phi; i++)
	{
		idx = i * (model->cellsNum_r + 2) * (model->cellsNum_z + 2);
		for (int j = 0; j < model->cellsNum_z + 2; j++)
		{
			it = model->Qcell.find(idx);
			if (it == model->Qcell.end())
			{
				a[counter++] = 1.0;
				a[counter++] = 0.0;
				a[counter++] = -1.0;
				a[counter++] = 0.0;

				a[counter++] = 0.0;
				a[counter++] = 1.0;
				a[counter++] = 0.0;
				a[counter++] = -1.0;

				rhs[2 * idx] = 0.0;
				rhs[2 * idx + 1] = 0.0;
			}
			else
				stencils->left->fill(idx, &counter);

			idx++;
		}
	}

	// Middle
	for (int i = 0; i < model->cellsNum_phi; i++)
	{
		idx = i * (model->cellsNum_r + 2) * (model->cellsNum_z + 2) + model->cellsNum_z + 2;
		for (int j = 1; j < model->cellsNum_r + 1; j++)
		{
			stencils->top->fill(idx, &counter);
			idx++;
			for (int k = 1; k < model->cellsNum_z + 1; k++)
			{	
				stencils->middle->fill(idx, &counter);
				idx++;
			}
			stencils->bot->fill(idx, &counter);
			idx++;
		}
	}

	// Right
	for (int i = 0; i < model->cellsNum_phi; i++)
	{
		idx = i * (model->cellsNum_r + 2) * (model->cellsNum_z + 2) + (model->cellsNum_r + 1) * (model->cellsNum_z + 2);
		for (int j = 0; j < model->cellsNum_z + 2; j++)
		{
			stencils->right->fill(idx, &counter);
			idx++;
		}
	}*/
}

void Par3DSolver::fillq()
{
	int i = 0;
	map<int, double>::iterator it = model->Qcell.begin();
	while (it != model->Qcell.end())
	{
		q[i++] = it->second;
		++it;
	}
}

void Par3DSolver::fillDq()
{
	for (int i = 0; i < n; i++)
		dq[i] = 0.0;
}

void Par3DSolver::solveDq(double mult)
{
	fillDq();
	filldPdQ(mult);
	solveStep();
	solveSystem();

	int i = 0;
	map<int, double>::iterator it;
	for (it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
	{
		cout << "dq[" << it->first << "] = " << dq[i++] * model->Q_dim * 86400.0 << endl;
	}
	cout << endl;
}

void Par3DSolver::solveSystem()
{
	double s = 0.0, p1, p2;
	map<int, double>::iterator it;

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

void Par3DSolver::filldPdQ(double mult)
{
	double p1, p2, ratio;
	ratio = mult * 0.001 / (double)(n);

	int i = 0, j = 0;
	map<int, double>::iterator it0 = model->Qcell.begin();
	map<int, double>::iterator it1 = model->Qcell.begin();
	map<int, double>::iterator it2 = it1;	++it2;
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