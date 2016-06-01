#include "model/3D/Perforation/OilPerfNITSolver.h"

#include <algorithm>
#include <map>

using namespace std;
using namespace oil_perf_nit;

OilPerfNITSolver::OilPerfNITSolver(Oil_Perf_NIT* _model) : AbstractSolver<Oil_Perf_NIT>(_model)
{
	// Output streams
	plot_Tdyn.open("snaps/T_dyn.dat", ofstream::out);
	plot_Pdyn.open("snaps/P_dyn.dat", ofstream::out);
	plot_qcells.open("snaps/q_cells.dat", ofstream::out);
	
	// Flow rate optimization structures
	n = model->Qcell.size();
	dq.Initialize(n);
	q.Initialize(n);
	dpdq.Initialize(n, n - 1);
	mat.Initialize(n - 1, n - 1);
	b.Initialize(n - 1);

	// Time settings
	T_dim = model->T_dim;
	t_dim = model->t_dim;
	Tt = model->period[model->period.size() - 1];

	// Memory allocating
	ind_i = new int[7 * model->cellsNum];
	ind_j = new int[7 * model->cellsNum];
	a = new double[7 * model->cellsNum];
	ind_rhs = new int[model->cellsNum + model->tunnelCells.size()];
	rhs = new double[model->cellsNum + model->tunnelCells.size()];

	tind_i = new int[7 * (model->cellsNum + model->tunnelCells.size())];
	tind_j = new int[7 * (model->cellsNum + model->tunnelCells.size())];
	ta = new double[7 * (model->cellsNum + model->tunnelCells.size())];
	tind_rhs = new int[model->cellsNum + model->tunnelCells.size()];
	trhs = new double[model->cellsNum + model->tunnelCells.size()];

	// Stencils allocating
	stencils = new UsedStencils<Oil_Perf_NIT>(model);
	stencils->setStorages(a, ind_i, ind_j, rhs);
}

OilPerfNITSolver::~OilPerfNITSolver()
{
	// Closing streams
	plot_Tdyn.close();
	plot_Pdyn.close();
	plot_qcells.close();
}

void OilPerfNITSolver::writeData()
{
	double t = 0.0, p = 0.0, s = 0.0, q = 0;

	plot_qcells << cur_t * t_dim / 3600.0;

	map<int, double>::iterator it;
	for (it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
	{
	//	t += model->tunnelCells[it->first].u_next.t;
	//	p += model->tunnelCells[it->first].u_next.p * model->P_dim;
	//	s += model->tunnelCells[it->first].u_next.s;
		if (model->leftBoundIsRate)
			plot_qcells << "\t" << it->second * model->Q_dim * 86400.0;
		else
		{
			plot_qcells << "\t" << model->getRate(it->first) * model->Q_dim * 86400.0;
			q += model->getRate(it->first);
		}
	}

	int counter = 0;
	int sum = 0;
	for (int i = 0; i < model->perfTunnels.size(); i++)
	{
		if (model->perfTunnels[i].second != 0)
		{
			t += model->tunnelCells[counter].u_next.t;
			p += model->tunnelCells[counter].u_next.p * model->P_dim;

			t += model->tunnelCells[counter].u_next.t;
			p += model->tunnelCells[counter].u_next.p * model->P_dim;

			t += model->tunnelCells[counter].u_next.t;
			p += model->tunnelCells[counter].u_next.p * model->P_dim;

			t += model->tunnelCells[counter].u_next.t;
			p += model->tunnelCells[counter].u_next.p * model->P_dim;

			sum += 4;
			counter += ((model->perfTunnels[i].second - 1) * 4 + 1);
		}
		else
		{
			t += model->tunnelCells[counter].u_next.t;
			p += model->tunnelCells[counter].u_next.p * model->P_dim;

			sum++;
		}
	}
	
	plot_Tdyn << cur_t * t_dim / 3600.0 << "\t" << t / (double)(sum) * T_dim << endl;
	plot_Pdyn << cur_t * t_dim / 3600.0 << "\t" << p / (double)(sum) / 100000.0 << endl;

	if (model->leftBoundIsRate)
		plot_qcells << "\t" << model->Q_sum * model->Q_dim * 86400.0 << endl;
	else
		plot_qcells << "\t" << q * model->Q_dim * 86400.0 << endl;
}

void OilPerfNITSolver::control()
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

void OilPerfNITSolver::start()
{
	int counter = 0;
	iterations = 8;

	fillIndices(PRES);
	fillIndices(TEMP);
	pres_solver.Init( model->cellsNum + model->tunnelCells.size() );
	temp_solver.Init( model->cellsNum + model->tunnelCells.size() );

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

void OilPerfNITSolver::doNextStep()
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

void OilPerfNITSolver::copySolution(const paralution::LocalVector<double>& sol, int key)
{
	if (key == PRES)
	{
		for (int i = 0; i < model->cellsNum; i++)
			model->cells[i].u_next.p += sol[i];

		for (int i = model->cellsNum; i < model->cellsNum + model->tunnelCells.size(); i++)
			model->tunnelCells[i - model->cellsNum].u_next.p += sol[i];
	}
	else if (key == TEMP)
	{
		for (int i = 0; i < model->cellsNum; i++)
			model->cells[i].u_next.t = sol[i];

		for (int i = model->cellsNum; i < model->cellsNum + model->tunnelCells.size(); i++)
			model->tunnelCells[i - model->cellsNum].u_next.t = sol[i];
	}
}

void OilPerfNITSolver::solveStep()
{
	int cellIdx, varIdx;
	double err_newton = 1.0;
	double averPresPrev = averValue(1);
	double averPres;
	double dAverPres = 1.0;

	iterations = 0;
	while (err_newton > 1.e-4 && /*(dAverSat > 1.e-8 || dAverPres > 1.e-4) &&*/ iterations < 10)
	{
		copyIterLayer();

		fill(PRES);
		pres_solver.Assemble(ind_i, ind_j, a, presElemNum, ind_rhs, rhs);
		pres_solver.Solve();
		copySolution( pres_solver.getSolution(), PRES );

		err_newton = convergance(cellIdx, varIdx);

		averPres = averValue(1);
		dAverPres = fabs(averPres - averPresPrev);
		averPresPrev = averPres;

		//		if(varIdx == PRES)
		//			cout << "BadPresValue[" << cellIdx  << "]: " << model->cells[cellIdx].u_next.p << endl;
		//		else if(varIdx == SAT)
		//			cout << "BadSatValue[" << cellIdx  << "]: " << model->cells[cellIdx].u_next.s << endl;

		iterations++;
	}

	fill(TEMP);
	temp_solver.Assemble(tind_i, tind_j, ta, tempElemNum, tind_rhs, trhs);
	temp_solver.Solve();
	copySolution(temp_solver.getSolution(), TEMP);

	cout << "Newton Iterations = " << iterations << endl;
}

void OilPerfNITSolver::fillIndices(int key)
{
	int idx, nebr;
	int counter = 0;
	Iterator it;
	map<int, double>::iterator itPerf;

	if (key == PRES)
	{
		// Left
		for (it = model->getLeftBegin(); it != model->getLeftEnd(); ++it)
		{
			idx = it.getIdx();

			if (it->isUsed)
			{
				nebr = idx + model->cellsNum_z + 2;
				ind_i[counter] = idx;
				ind_j[counter++] = idx;

				ind_i[counter] = idx;
				ind_j[counter++] = nebr;
			}
			else {
				ind_i[counter] = idx;
				ind_j[counter++] = idx;
			}
		}

		// Middle
		int res;
		for (it = model->getMidBegin(); it != model->getMidEnd(); ++it)
		{
			idx = it.getIdx();
			res = idx % (model->cellsNum_z + 2);
			if (res == 0)
				stencils->top->fillIndex(idx, &counter);
			else if (res == model->cellsNum_z + 1)
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

		presElemNum = counter;

		for (int i = 0; i < model->cellsNum + model->tunnelCells.size(); i++)
			ind_rhs[i] = i;
	}
	else if (key == TEMP)
	{
		// Left
		for (it = model->getLeftBegin(); it != model->getLeftEnd(); ++it)
		{
			idx = it.getIdx();

			if (it->isUsed)
			{
				nebr = idx + model->cellsNum_z + 2;
				tind_i[counter] = idx;
				tind_j[counter++] = idx;

				tind_i[counter] = idx;
				tind_j[counter++] = nebr;
			}
			else {
				tind_i[counter] = idx;
				tind_j[counter++] = idx;
			}
		}

		// Middle
		int res;
		for (it = model->getMidBegin(); it != model->getMidEnd(); ++it)
		{
			idx = it.getIdx();
			res = idx % (model->cellsNum_z + 2);
			if (res == 0)
			{
				tind_i[counter] = idx;
				tind_j[counter++] = idx;

				tind_i[counter] = idx;
				tind_j[counter++] = idx + 1;
			}
			else if (res == model->cellsNum_z + 1)
			{
				tind_i[counter] = idx;
				tind_j[counter++] = idx;

				tind_i[counter] = idx;
				tind_j[counter++] = idx - 1;
			} 
			else
			{
				Cell* nebr[7];
				model->getStencilIdx(idx, nebr);

				if (nebr[0]->isUsed)
				{
					for (int j = 0; j < 7; j++)
					{
						tind_i[counter] = idx;
						if (nebr[j]->isUsed)
							tind_j[counter++] = nebr[j]->num;
						else
							tind_j[counter++] = model->cellsNum + model->getCell(nebr[0]->num, nebr[j]->num).num;
					}
				}
				else
				{
					tind_i[counter] = idx;
					tind_j[counter++] = idx;
				}
			}
		}

		// Right
		for (it = model->getRightBegin(); it != model->getRightEnd(); ++it)
		{
			idx = it.getIdx();

			tind_i[counter] = idx;
			tind_j[counter++] = idx;
		}

		// Tunnel cells
		vector<Cell>::iterator itr;
		for (itr = model->tunnelCells.begin(); itr != model->tunnelCells.end(); ++itr)
		{
			tind_i[counter] = itr->num + model->cellsNum;
			tind_j[counter++] = itr->num + model->cellsNum;

			tind_i[counter] = itr->num + model->cellsNum;
			tind_j[counter++] = model->nebrMap[itr->num].first;

			tind_i[counter] = itr->num + model->cellsNum;
			tind_j[counter++] = model->nebrMap[itr->num].second;
		}

		tempElemNum = counter;

		for (int i = 0; i < model->cellsNum + model->tunnelCells.size(); i++)
			tind_rhs[i] = i;
	}
}

void OilPerfNITSolver::fill(int key)
{
	int idx;
	int counter = 0;
	Iterator it;
	map<int, double>::iterator itPerf;

	if (key == PRES)
	{
		// Left
		for (it = model->getLeftBegin(); it != model->getLeftEnd(); ++it)
		{
			idx = it.getIdx();

			if (it->isUsed)
			{
				a[counter++] = 1.0;
				a[counter++] = -1.0;

				rhs[idx] = 0.0;
			}
			else {
				a[counter++] = 1.0;

				rhs[idx] = 0.0;
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
	else if (key == TEMP)
	{
		// Left
		for (it = model->getLeftBegin(); it != model->getLeftEnd(); ++it)
		{
			idx = it.getIdx();

			if (it->isUsed)
			{
				ta[counter++] = 1.0;
				ta[counter++] = -1.0;

				trhs[idx] = 0.0;
			}
			else {
				ta[counter++] = 1.0;
				trhs[idx] = 0.0;
			}
		}

		// Middle
		int res;
		for (it = model->getMidBegin(); it != model->getMidEnd(); ++it)
		{
			idx = it.getIdx();
			res = idx % (model->cellsNum_z + 2);
			if (res == 0)
			{
				ta[counter++] = 1.0;
				ta[counter++] = -1.0;
				trhs[idx] = 0.0;
			}
			else if (res == model->cellsNum_z + 1)
			{
				ta[counter++] = 1.0;
				ta[counter++] = -1.0;
				trhs[idx] = 0.0;
			}
			else
			{
				Cell* nebr[7];
				model->getStencilIdx(idx, nebr);

				if (nebr[0]->isUsed)
				{
					ta[counter+1] = -2.0 * (max(model->getA(*nebr[0], NEXT, R_AXIS), 0.0) +
						model->getLambda(*nebr[0], *nebr[1]) * (nebr[0]->r - nebr[0]->hr / 2.0) / nebr[0]->r / nebr[0]->hr) / (nebr[0]->hr + nebr[1]->hr);
					ta[counter+2] = 2.0 * (min(model->getA(*nebr[0], NEXT, R_AXIS), 0.0) -
						model->getLambda(*nebr[0], *nebr[2]) * (nebr[0]->r + nebr[0]->hr / 2.0) / nebr[0]->r / nebr[0]->hr) / (nebr[0]->hr + nebr[2]->hr);
					ta[counter+3] = -2.0 * (max(model->getA(*nebr[0], NEXT, Z_AXIS), 0.0) +
						model->getLambda(*nebr[0], *nebr[3]) / nebr[0]->hz) / (nebr[0]->hz + nebr[3]->hz);
					ta[counter+4] = 2.0 * (min(model->getA(*nebr[0], NEXT, Z_AXIS), 0.0) -
						model->getLambda(*nebr[0], *nebr[4]) / nebr[0]->hz) / (nebr[0]->hz + nebr[4]->hz);
					ta[counter+5] = -2.0 * (max(model->getA(*nebr[0], NEXT, PHI_AXIS), 0.0) +
						model->getLambda(*nebr[0], *nebr[5]) / nebr[0]->r / nebr[0]->hphi) / nebr[0]->r / (nebr[0]->hphi + nebr[5]->hphi);
					ta[counter+6] = 2.0 * (min(model->getA(*nebr[0], NEXT, PHI_AXIS), 0.0) -
						model->getLambda(*nebr[0], *nebr[6]) / nebr[0]->r / nebr[0]->hphi) / nebr[0]->r / (nebr[0]->hphi + nebr[6]->hphi);
					ta[counter] = model->getCn(*nebr[0]) / model->ht 
						- ta[counter + 1]
						- ta[counter + 2]
						- ta[counter + 3]
						- ta[counter + 4]
						- ta[counter + 5]
						- ta[counter + 6];

					trhs[idx] = model->getCn(*nebr[0]) * nebr[0]->u_prev.t / model->ht +
						model->getAd(*nebr[0]) * (nebr[0]->u_next.p - nebr[0]->u_prev.p) / model->ht -
						model->getJT(*nebr[0], NEXT, R_AXIS) * model->getNablaP(*nebr[0], NEXT, R_AXIS) -
						model->getJT(*nebr[0], NEXT, PHI_AXIS) * model->getNablaP(*nebr[0], NEXT, PHI_AXIS) -
						model->getJT(*nebr[0], NEXT, Z_AXIS) * model->getNablaP(*nebr[0], NEXT, Z_AXIS);

					counter += 7;
				}
				else
				{
					ta[counter++] = 1.0;
					trhs[idx] = 0.0;
				}
			}
		}

		// Right
		for (it = model->getRightBegin(); it != model->getRightEnd(); ++it)
		{
			idx = it.getIdx();

			ta[counter++] = 1.0;

			trhs[idx] = model->props_sk[model->getSkeletonIdx(model->cells[idx])].t_init;
		}

		// Tunnel cells
		vector<Cell>::iterator itr;
		for (itr = model->tunnelCells.begin(); itr != model->tunnelCells.end(); ++itr)
		{
			Cell& cell = *itr;
			Cell& nebr1 = model->getCell(model->nebrMap[itr->num].first);
			Cell& nebr2 = model->getCell(model->nebrMap[itr->num].second);

			if (fabs(nebr2.r - nebr1.r) > EQUALITY_TOLERANCE)
			{
				ta[counter++] = 1.0 / (nebr1.r - cell.r);
				ta[counter++] = -1.0 / (nebr2.r - nebr1.r) - 1.0 / (nebr1.r - cell.r);
				ta[counter++] = 1.0 / (nebr2.r - nebr1.r);
			}
			else if (fabs(nebr2.z - nebr1.z) > EQUALITY_TOLERANCE)
			{
				ta[counter++] = 1.0 / (nebr1.z - cell.z);
				ta[counter++] = -1.0 / (nebr2.z - nebr1.z) - 1.0 / (nebr1.z - cell.z);
				ta[counter++] = 1.0 / (nebr2.z - nebr1.z);
			}
			else if (fabs(nebr2.phi - nebr1.phi) > EQUALITY_TOLERANCE)
			{
				ta[counter++] = 1.0 / (nebr1.phi - cell.phi) / nebr1.r;
				ta[counter++] = -1.0 / (nebr2.phi - nebr1.phi) / nebr1.r - 1.0 / (nebr1.phi - cell.phi) / nebr1.r;
				ta[counter++] = 1.0 / (nebr2.phi - nebr1.phi) / nebr1.r;
			}

			trhs[itr->num + model->cellsNum] = 0.0;
		}
	}
}

void OilPerfNITSolver::fillq()
{
	int i = 0;
	map<int, double>::iterator it = model->Qcell.begin();
	while (it != model->Qcell.end())
	{
		q[i++] = it->second;
		++it;
	}
}

void OilPerfNITSolver::fillDq()
{
	for (int i = 0; i < n; i++)
		dq[i] = 0.0;
}

void OilPerfNITSolver::solveDq(double mult)
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

void OilPerfNITSolver::solveSystem()
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

void OilPerfNITSolver::filldPdQ(double mult)
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