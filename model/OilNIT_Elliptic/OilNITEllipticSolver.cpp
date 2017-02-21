#include "model/OilNIT_Elliptic/OilNITEllipticSolver.hpp"

using namespace std;
using namespace oilnit_elliptic;

OilNITEllipticSolver::OilNITEllipticSolver(OilNIT_Elliptic* _model) : AbstractSolver<OilNIT_Elliptic>(_model)
{
	// Output streams
	plot_Pdyn.open("snaps/P_dyn.dat", ofstream::out);
	plot_Tdyn.open("snaps/T_dyn.dat", ofstream::out);
	plot_qcells.open("snaps/q_cells.dat", ofstream::out);

	// Flow rate optimization structures
	n = model->Qcell_ellipse.size();
	
	dq.Initialize(n);
	q.Initialize(n);
	dpdq.Initialize(n, n - 1);
	mat.Initialize(n - 1, n - 1);
	b.Initialize(n - 1);

	// Time settings
	t_dim = model->t_dim;
	Tt = model->period[model->period.size() - 1];

	// Memory allocating
	ind_i = new int[stencil * (model->cellsNum + model->wellCells.size())];
	ind_j = new int[stencil * (model->cellsNum + model->wellCells.size())];
	tind_i = new int[stencil * (model->cellsNum + model->wellCells.size())];
	tind_j = new int[stencil * (model->cellsNum + model->wellCells.size())];
	a = new double[stencil * (model->cellsNum + model->wellCells.size())];
	ind_rhs = new int[model->cellsNum + model->wellCells.size()];
	rhs = new double[model->cellsNum + model->wellCells.size()];
}
OilNITEllipticSolver::~OilNITEllipticSolver()
{
	plot_Pdyn.close();
	plot_Tdyn.close();
	plot_qcells.close();

	delete[] ind_i, ind_j, tind_i, tind_j, ind_rhs;
	delete[] a, rhs;

}
void OilNITEllipticSolver::writeData()
{
	double p = 0.0, t = 0.0, q = 0.0;

	plot_qcells << cur_t * t_dim / 3600.0;
	plot_Tdyn << cur_t * t_dim / 3600.0;
	plot_Pdyn << cur_t * t_dim / 3600.0;

	map<int, double>::iterator it;
	for (it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
	{
		p += model->wellCells[it->first].u_next.p * model->P_dim;
		t += model->wellCells[it->first].u_next.t * model->T_dim;
		if (model->leftBoundIsRate) {
			plot_Pdyn << "\t" << model->wellCells[it->first].u_next.p * model->P_dim / BAR_TO_PA;
			plot_qcells << "\t" << it->second * model->Q_dim * 86400.0;
		} else
		{
			plot_qcells << "\t" << model->getRate(it->first) * model->Q_dim * 86400.0;
			q += model->getRate(it->first);
		}
		plot_Tdyn << "\t" << model->wellCells[it->first].u_next.t * model->T_dim;
	}

	plot_Pdyn << "\t" << p / (double)(model->Qcell.size()) / BAR_TO_PA << endl;

	if (model->leftBoundIsRate)
		plot_qcells << "\t" << model->Q_sum * model->Q_dim * 86400.0 << endl;
	else
		plot_qcells << "\t" << q * model->Q_dim * 86400.0 << endl;

	plot_Tdyn << "\t" << t / (double)(model->Qcell.size()) << endl;
}
void OilNITEllipticSolver::control()
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
void OilNITEllipticSolver::start()
{
	int counter = 0;
	iterations = 8;

	fillIndices();
	pres_solver.Init(model->cellsNum + model->wellCells.size());
	temp_solver.Init(model->cellsNum + model->wellCells.size());

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
void OilNITEllipticSolver::copySolution(const paralution::LocalVector<double>& sol, const int val)
{
	for (int i = 0; i < model->cellsNum; i++)
		model->cells[i].u_next.values[val] += sol[i];

	for (int i = model->cellsNum; i < model->cellsNum + model->wellCells.size(); i++)
		model->wellCells[i - model->cellsNum].u_next.values[val] += sol[i];
}
void OilNITEllipticSolver::solveStep()
{
	int cellIdx, varIdx;
	double err_newton = 1.0;
	double averPresPrev = averValue(0);
	double averPres;
	double dAverPres = 1.0;

	iterations = 0;
	while (err_newton > 1.e-4 /*&& (dAverSat > 1.e-9 || dAverPres > 1.e-7)*/ && iterations < 10)
	{
		copyIterLayer();

		fill(PRES);
		pres_solver.Assemble(ind_i, ind_j, a, elemNum, ind_rhs, rhs);
		pres_solver.Solve();
		copySolution(pres_solver.getSolution(), PRES);

		err_newton = convergance(cellIdx, varIdx);

		averPres = averValue(0);				
		dAverPres = fabs(averPres - averPresPrev);	
		averPresPrev = averPres;

		iterations++;
	}

	fill(TEMP);
	temp_solver.Assemble(tind_i, tind_j, a, telemNum, ind_rhs, rhs);
	temp_solver.Solve();
	copySolution(temp_solver.getSolution(), TEMP);

	cout << "Newton Iterations = " << iterations << endl;
}
void OilNITEllipticSolver::doNextStep()
{
	solveStep();

	if (n > 1 && model->Q_sum != -model->Q_sum)
	{
		double H0 = model->P_dim * model->P_dim * fabs(model->solveH()) / BAR_TO_PA / BAR_TO_PA;
		if (H0 > 0.1)
		{
			printWellRates();
			fillq();

			double mult = 0.9;
			double H = H0;

			while (H > H0 / 1000.0 || H > 0.05)
			{
				solveDq(mult);

				int i = 0;
				auto it = model->Qcell_ellipse.begin();

				while (it != model->Qcell_ellipse.end())
				{
					q[i] += mult * dq[i];
					it->second = q[i];
					auto indices = model->getPerforationIndicesSameType(it->first);
					for (auto& idx : indices)
						model->Qcell[idx] = it->second;

					i++;	++it;
				}

				solveStep();
				printWellRates();

				H = fabs(model->solveH());
			}
		}
	}
}
void OilNITEllipticSolver::fill(const int val)
{
	int counter = 0;
	int i;
	Cell* neighbor[7];

	for (const auto& cell : model->cells)
	{
		if (cell.isUsed)
		{
			model->setVariables(cell, val);

			i = 0;
			const auto mat_idx = getMatrixStencil(cell, val);
			for (const int idx : mat_idx)
			{
				a[counter++] = model->jac[0][i];	i++;
			}

			rhs[cell.num] = -model->y[0];
		}
		else
		{
			a[counter++] = 1.0;		rhs[cell.num] = 0.0;
		}
	}

	for (const auto& cell : model->wellCells)
	{
		model->setVariables(cell, val);

		i = 0;
		const auto mat_idx = getMatrixStencil(cell, val);
		for (const int idx : mat_idx)
		{
			a[counter++] = model->jac[0][i];
			i++;
		}

		rhs[cell.num + model->cellsNum] = -model->y[0];
	}
}
void OilNITEllipticSolver::fillIndices()
{
	int pres_counter = 0, temp_counter = 0;

	for (const auto& cell : model->cells)
	{
		if (cell.isUsed)
		{
			const auto pres_idx = getMatrixStencil(cell, PRES);
			for (const int idx : pres_idx)
			{
				ind_i[pres_counter] = cell.num;			ind_j[pres_counter++] = idx;
			}

			const auto temp_idx = getMatrixStencil(cell, TEMP);
			for (const int idx : temp_idx)
			{
				tind_i[temp_counter] = cell.num;			tind_j[temp_counter++] = idx;
			}
		}
		else
		{
			ind_i[pres_counter] = cell.num;			ind_j[pres_counter++] = cell.num;
			tind_i[temp_counter] = cell.num;		tind_j[temp_counter++] = cell.num;
		}
	}

	for (const auto& cell : model->wellCells)
	{
		const auto pres_idx = getMatrixStencil(cell, PRES);
		for (const int idx : pres_idx)
		{
			ind_i[pres_counter] = cell.num + model->cellsNum;			ind_j[pres_counter++] = idx;
		}

		const auto temp_idx = getMatrixStencil(cell, TEMP);
		for (const int idx : temp_idx)
		{
			tind_i[temp_counter] = cell.num + model->cellsNum;			tind_j[temp_counter++] = idx;
		}
	}

	elemNum = pres_counter;		telemNum = temp_counter;

	for (int i = 0; i < model->cellsNum + model->wellCells.size(); i++)
		ind_rhs[i] = i;
}

void OilNITEllipticSolver::fillq()
{
	int i = 0;
	map<int, double>::iterator it = model->Qcell_ellipse.begin();
	while (it != model->Qcell_ellipse.end())
	{
		q[i++] = it->second;
		++it;
	}
}
void OilNITEllipticSolver::fillDq()
{
	for (int i = 0; i < n; i++)
		dq[i] = 0.0;
}
void OilNITEllipticSolver::solveDq(double mult)
{
	fillDq();
	filldPdQ(mult);
	solveStep();
	solveSystem();

	int i = 0;
	map<int, double>::iterator it;
	for (it = model->Qcell_ellipse.begin(); it != model->Qcell_ellipse.end(); ++it)
	{
		cout << "dq[" << it->first << "] = " << dq[i++] * model->Q_dim * 86400.0 << endl;
	}
	cout << endl;
}
void OilNITEllipticSolver::solveSystem()
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
		it = model->Qcell_ellipse.begin();
		for (int k = 0; k < n - 1; k++)
		{
			p1 = model->wellCells[it->first].u_next.p;
			p2 = model->wellCells[(++it)->first].u_next.p;
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
void OilNITEllipticSolver::filldPdQ(double mult)
{
	double p1, p2, ratio;
	ratio = mult * 0.001 / (double)(n);

	int i = 0, j = 0;
	map<int, double>::iterator it0 = model->Qcell_ellipse.begin();
	map<int, double>::iterator it1 = model->Qcell_ellipse.begin();
	map<int, double>::iterator it2 = it1;	++it2;
	while (it1 != model->Qcell_ellipse.end())
	{
		j = 0;
		it2 = model->Qcell_ellipse.begin();		++it2;
		while (it2 != model->Qcell_ellipse.end())
		{
			model->setRateDeviation(it2->first, -ratio);
			model->setRateDeviation(it0->first, ratio);
			solveStep();
			p1 = model->wellCells[it1->first].u_next.p;

			model->setRateDeviation(it2->first, 2.0 * ratio);
			model->setRateDeviation(it0->first, -2.0 * ratio);
			solveStep();
			p2 = model->wellCells[it1->first].u_next.p;

			model->setRateDeviation(it2->first, -ratio);
			model->setRateDeviation(it0->first, ratio);

			dpdq[i][j++] = (p2 - p1) / (2.0 * ratio * model->Q_sum_quater);

			++it2;
		}
		i++;
		++it1;
	}
}