#include "model/3D/Perforation/ParPerfSolver.h"

#include <algorithm>
#include <map>

using namespace std;
using namespace gasOil_perf;

ParPerfSolver::ParPerfSolver(GasOil_Perf* _model) : AbstractSolver<GasOil_Perf>(_model)
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
	ind_rhs = new int[2 * (model->cellsNum + model->tunnelCells.size())];
	rhs = new double[2 * (model->cellsNum + model->tunnelCells.size())];

	// Stencils allocating
	stencils = new UsedStencils<GasOil_Perf>(model);
	stencils->setStorages(a, ind_i, ind_j, rhs);
}

ParPerfSolver::~ParPerfSolver()
{
	// Closing streams
	plot_Pdyn.close();
	plot_Sdyn.close();
	plot_qcells.close();
}

void ParPerfSolver::writeData()
{
	double p = 0.0, s = 0.0, q = 0;

	plot_qcells << cur_t * t_dim / 3600.0;

	map<int, double>::iterator it;
	for (it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
	{
		p += model->tunnelCells[it->first].u_next.p * model->P_dim;
		s += model->tunnelCells[it->first].u_next.s;
		if (model->leftBoundIsRate)
			plot_qcells << "\t" << it->second * model->Q_dim * 86400.0;
		else
		{
			plot_qcells << "\t" << model->getRate(it->first) * model->Q_dim * 86400.0;
			q += model->getRate(it->first);
		}
	}
	
	plot_Pdyn << cur_t * t_dim / 3600.0 << "\t" << p / (double)(model->Qcell.size()) << endl;
	plot_Sdyn << cur_t * t_dim / 3600.0 << "\t" << s / (double)(model->Qcell.size()) << endl;

	if (model->leftBoundIsRate)
		plot_qcells << "\t" << model->Q_sum * model->Q_dim * 86400.0 << endl;
	else
		plot_qcells << "\t" << q * model->Q_dim * 86400.0 << endl;
}

void ParPerfSolver::control()
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

void ParPerfSolver::start()
{
	int counter = 0;
	iterations = 8;

	fillIndices();
	solver.Init(2 * (model->cellsNum + model->tunnelCells.size()));

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

void ParPerfSolver::doNextStep()
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

void ParPerfSolver::copySolution(const paralution::LocalVector<double>& sol)
{
	for (int i = 0; i < model->cellsNum; i++)
	{
		model->cells[i].u_next.p += sol[2 * i];
		model->cells[i].u_next.s += sol[2 * i + 1];
	}

	for (int i = model->cellsNum; i < model->cellsNum + model->tunnelCells.size(); i++)
	{
		model->tunnelCells[i-model->cellsNum].u_next.p += sol[2 * i];
		model->tunnelCells[i-model->cellsNum].u_next.s += sol[2 * i + 1];
	}
}

void ParPerfSolver::solveStep()
{
	int cellIdx, varIdx;
	double err_newton = 1.0;
	double averPresPrev = averValue(0);
	double averSatPrev = averValue(1);
	double averPres, averSat;
	double dAverPres = 1.0, dAverSat = 1.0;

	iterations = 0;
	while (err_newton > 1.e-4 /*&& (dAverSat > 1.e-9 || dAverPres > 1.e-7)*/ && iterations < 10)
	{
		copyIterLayer();

		fill();
		solver.Assemble(ind_i, ind_j, a, elemNum, ind_rhs, rhs);
		solver.Solve();
		copySolution( solver.getSolution() );

		model->solveP_bub();

		err_newton = convergance(cellIdx, varIdx);

		averPres = averValue(0);					averSat = averValue(1);
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

void ParPerfSolver::fillIndices()
{
	int idx, nebr;
	int counter = 0;
	Iterator it;
	map<int, double>::iterator itPerf;

	// Left
	for (it = model->getLeftBegin(); it != model->getLeftEnd(); ++it)
	{
		idx = it.getIdx();

		if (it->isUsed)
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
		else {
			ind_i[counter] = 2 * idx;
			ind_j[counter++] = 2 * idx;

			ind_i[counter] = 2 * idx + 1;
			ind_j[counter++] = 2 * idx + 1;
		}
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

	// Tunnel cells
	vector<Cell>::iterator itr;
	for (itr = model->tunnelCells.begin(); itr != model->tunnelCells.end(); ++itr)
	{
		stencils->left->fillIndex(itr->num, &counter);
	}

	elemNum = counter;

	for (int i = 0; i < 2*(model->cellsNum + model->tunnelCells.size()); i++)
		ind_rhs[i] = i;
}

void ParPerfSolver::fill()
{
	int idx;
	int counter = 0;
	Iterator it;
	map<int, double>::iterator itPerf;

	// Left
	for (it = model->getLeftBegin(); it != model->getLeftEnd(); ++it)
	{
		idx = it.getIdx();

		if (it->isUsed)
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
		else {
			a[counter++] = 1.0;
			a[counter++] = 1.0;

			rhs[2 * idx] = 0.0;
			rhs[2 * idx + 1] = 0.0;
		}
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

	// Tunnel cells
	vector<Cell>::iterator itr;
	for (itr = model->tunnelCells.begin(); itr != model->tunnelCells.end(); ++itr)
	{
		stencils->left->fill(itr->num, &counter);
	}
}

void ParPerfSolver::fillq()
{
	int i = 0;
	map<int, double>::iterator it = model->Qcell.begin();
	while (it != model->Qcell.end())
	{
		q[i++] = it->second;
		++it;
	}
}

void ParPerfSolver::fillDq()
{
	for (int i = 0; i < n; i++)
		dq[i] = 0.0;
}

void ParPerfSolver::solveDq(double mult)
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

void ParPerfSolver::solveSystem()
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

void ParPerfSolver::filldPdQ(double mult)
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