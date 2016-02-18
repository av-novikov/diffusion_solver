#include "model/3D/GasOil_3D/Par3DSolver.h"

using namespace std;
using namespace gasOil_3d;

Par3DSolver::Par3DSolver(GasOil_3D* _model) : AbstractSolver<GasOil_3D>(_model)
{
	plot_Pdyn.open("snaps/P_dyn.dat", ofstream::out);
	plot_Sdyn.open("snaps/S_dyn.dat", ofstream::out);
	plot_qcells.open("snaps/q_cells.dat", ofstream::out);

	t_dim = model->t_dim;

	n = model->Qcell.size();
	dq.Initialize(n);
	q.Initialize(n);

	dpdq.Initialize(n, n - 1);
	mat.Initialize(n - 1, n - 1);
	b.Initialize(n - 1);

	Tt = model->period[model->period.size() - 1];

	ind_i = new int[ 7 * 2 * model->cellsNum ];
	ind_j = new int[ 7 * 2 * model->cellsNum ];
	a = new double[ 7 * 2 * model->cellsNum ];

	ind_rhs = new int[ 2 * model->cellsNum ];
	rhs = new double[ 2 * model->cellsNum ];
}

Par3DSolver::~Par3DSolver()
{
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
		p += model->cells[it->first].u_next.p;
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

		Solve();

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

void Par3DSolver::Solve()
{
	Rhs.Zeros();

	int counter = 0;

	// Left cells
	fillLeft();

	// Middle cells
	fillMiddle();

	// Right cells
	fillRight();
	
	Mat.Assemble(ind_i, ind_j, a, counter, "A", 2 * model->cellsNum, 2 * model->cellsNum);

}

void Par3DSolver::fillLeft()
{
	for (auto it : model->Qcell)
	{
		Cell& curr = model->cells[it.first];
		Cell& nebr = model->cells[it.first + model->cellsNum_z + 2];
		int idx = 2 * ((it.first % ((model->cellsNum_z + 2) * (model->cellsNum_r + 2))) + int(it.first / ((model->cellsNum_z + 2) * (model->cellsNum_r + 2))) * (model->cellsNum_z + 2));

		// First eqn
		ind_i[counter] = 2 * it.first;
		ind_j[counter] = 2 * it.first;
		a[counter++] = model->solve_eqLeft_dp(it.first);

		ind_i[counter] = 2 * it.first;
		ind_j[counter] = 2 * it.first + 1;
		a[counter++] = model->solve_eqLeft_ds(it.first);

		ind_i[counter] = 2 * it.first;
		ind_j[counter] = 2 * (it.first + model->cellsNum_z + 2);
		a[counter++] = model->solve_eqLeft_dp_beta(it.first);

		ind_i[counter] = 2 * it.first;
		ind_j[counter] = 2 * (it.first + model->cellsNum_z + 2) + 1;
		a[counter++] = model->solve_eqLeft_ds_beta(it.first);

		rhs[2 * it.first] = -model->solve_eqLeft(it.first) +
			a[counter - 4] * curr.u_next.p + a[counter - 3] * curr.u_next.s +
			a[counter - 2] * nebr.u_next.p + a[counter - 1] * nebr.u_next.s;


		// Second eqn
		ind_i[counter] = 2 * it.first + 1;
		ind_j[counter] = 2 * it.first + 1;
		a[counter++] = 1.0;

		ind_i[counter] = 2 * it.first + 1;
		ind_j[counter] = 2 * (it.first + model->cellsNum_z + 2) + 1;
		a[counter++] = (curr.r - nebr.r) / (model->cells[it.first + 2 * model->cellsNum_z + 4].r - nebr.r) - 1.0;

		ind_i[counter] = 2 * it.first + 1;
		ind_j[counter] = 2 * (it.first + 2 * (model->cellsNum_z + 2)) + 1;
		a[counter++] = -(curr.r - nebr.r) / (model->cells[it.first + 2 * model->cellsNum_z + 4].r - nebr.r);

		rhs[2 * it.first + 1] = 0.0;
	}
}

void Par3DSolver::fillMiddle()
{

}

void Par3DSolver::fillRight()
{

}

void Par3DSolver::fillTop(int i)
{

}

void Par3DSolver::fillBottom(int i)
{

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