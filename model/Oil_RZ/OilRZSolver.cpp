#include "model/Oil_RZ/OilRZSolver.h"

using namespace std;
using namespace oil_rz;

OilRZSolver::OilRZSolver(Oil_RZ* _model) : AbstractSolver<Oil_RZ>(_model)
{
	Initialize(model->cellsNum_r+2, model->cellsNum_z+2);

	plot_Pdyn.open("snaps/P_dyn.dat", ofstream::out);
	plot_Pavg.open("snaps/Pavg.dat", ofstream::out);
	plot_qcells.open("snaps/q_cells.dat", ofstream::out);

	t_dim = model->t_dim;

	isChange = false;
	n = model->Qcell.size();
	dq.Initialize(n);
	q.Initialize(n);

	dpdq.Initialize(n, n-1);
	mat.Initialize(n-1, n-1);
	b.Initialize(n-1);

	Tt = model->period[model->period.size()-1];
}

OilRZSolver::~OilRZSolver()
{
	plot_Pdyn.close();
	plot_Pavg.close();
	plot_qcells.close();
}

void OilRZSolver::writeData()
{
	double p = 0.0, q = 0.0, p_avg = 0.0;

	plot_qcells << cur_t * t_dim / 3600.0;

	map<int, double>::iterator it;
	for (it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
	{
		p += model->cells[it->first].u_next.p * model->P_dim;
		if (model->leftBoundIsRate)
			plot_qcells << "\t" << it->second * model->Q_dim * 86400.0;
		else
		{
			plot_qcells << "\t" << model->getRate(it->first) * model->Q_dim * 86400.0;
			q += model->getRate(it->first);
		}
	}

	for (int i = 0; i < model->cellsNum; i++)
		p_avg += model->cells[i].u_next.p * model->cells[i].V * model->P_dim;
	p_avg /= model->Volume;

	plot_Pdyn << cur_t * t_dim / 3600.0 << "\t" << p / (double)(model->Qcell.size()) << endl;
	plot_Pavg << cur_t * t_dim / 3600.0 << "\t" << p_avg << endl;

	if (model->leftBoundIsRate)
		plot_qcells << "\t" << model->Q_sum * model->Q_dim * 86400.0 << endl;
	else
		plot_qcells << "\t" << q * model->Q_dim * 86400.0 << endl;
}

void OilRZSolver::control()
{
	writeData();

	if(cur_t >= model->period[curTimePeriod])
	{
		curTimePeriod++;
		model->ht = model->ht_min;
		model->setPeriod(curTimePeriod);
		isChange = true;
	}
	else
		isChange = false;

	if(model->ht <= model->ht_max && iterations < 6)
		model->ht = model->ht * 1.5;
	else if(iterations > 6 && model->ht > model->ht_min)
		model->ht = model->ht / 1.5;

	if(cur_t + model->ht > model->period[curTimePeriod])
		model->ht = model->period[curTimePeriod] - cur_t;

	cur_t += model->ht;
}

void OilRZSolver::doNextStep()
{
	solveStep();

	if(n > 1 && model->Q_sum > EQUALITY_TOLERANCE)
	{
		double H0 = fabs( model->solveH() );
		if( H0 > 0.1 )
		{
			printWellRates();
			fillq();

			double mult = 0.9;
			double H = H0;			

			while(H > H0 / 50.0 || H > 0.05)
			{
				solveDq(mult);

				int i = 0;
				map<int,double>::iterator it = model->Qcell.begin();

				while(it != model->Qcell.end())
				{
					//if(isChange)
					//	q[i] += 5.0 * mult * H * dq[i];
					//else
						q[i] += mult * dq[i];
					it->second = q[i];
					i++;	++it;
				}
				
				solveStep();
				printWellRates();

				H = fabs( model->solveH() );
			}
		}
	}
}

void OilRZSolver::fillq()
{
	int i = 0;
	map<int,double>::iterator it = model->Qcell.begin();
	while(it != model->Qcell.end())
	{
		q[i++] = it->second;
		++it;
	}
}

void OilRZSolver::fillDq()
{
	for(int i = 0; i < n; i++)
		dq[i] = 0.0;
}

void OilRZSolver::solveDq(double mult)
{
	fillDq();
	filldPdQ(mult);
	solveStep();
	solveSystem();

	int i = 0;
	map<int,double>::iterator it;
	for(it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
	{
		cout << "dq[" << it->first << "] = " << dq[i++] * model->Q_dim * 86400.0 << endl;
	}
	cout << endl;
}

void OilRZSolver::solveSystem()
{
	double s = 0.0, p1, p2;
	map<int,double>::iterator it;

	for(int i = 0; i < n-1; i++)
	{
		for(int j = 0; j < n-1; j++)
		{
			s = 0.0;
			for(int k = 0; k < n-1; k++)
				s += ( dpdq[k+1][j] - dpdq[k][j] ) * ( dpdq[k+1][i] - dpdq[k][i] );
			mat[i][j] = s;
		}
		
		s = 0.0;
		it = model->Qcell.begin();
		for(int k = 0; k < n-1; k++)
		{
			p1 = model->cells[ it->first ].u_next.p;
			p2 = model->cells[ (++it)->first ].u_next.p;
			s += ( p2 - p1 ) * ( dpdq[k+1][i] - dpdq[k][i] );
		}
		b[i] = -s;
	}

	MC_LU rateSystem(mat, b);
	rateSystem.LU_Solve();

	s = 0.0;
	for(int i = 0; i < n-1; i++)
	{
		dq[i+1] = rateSystem.ptResult[i];
		s += rateSystem.ptResult[i];
	}
	dq[0] = - s;
}

void OilRZSolver::filldPdQ(double mult)
{
	double p1, p2, ratio;
	ratio = mult * 0.001 / (double)(n);
	
	int i = 0, j = 0;
	map<int,double>::iterator it0 = model->Qcell.begin();
	map<int,double>::iterator it1 = model->Qcell.begin();
	map<int,double>::iterator it2 = it1;	++it2;
	while(it1 != model->Qcell.end())
	{
		j = 0;
		it2 = model->Qcell.begin();		++it2;
		while(it2 != model->Qcell.end())
		{
			model->setRateDeviation(it2->first, -ratio);
			model->setRateDeviation(it0->first, ratio);
			solveStep();
			p1 = model->cells[ it1->first ].u_next.p;

			model->setRateDeviation(it2->first, 2.0 * ratio);
			model->setRateDeviation(it0->first, -2.0 * ratio);
			solveStep();
			p2 = model->cells[ it1->first ].u_next.p;

			model->setRateDeviation(it2->first, -ratio);
			model->setRateDeviation(it0->first, ratio);

			dpdq[i][j++] = (p2 - p1) / ( 2.0 * ratio * model->Q_sum);

			++it2;
		}
		i++;
		++it1;
	}
}

void OilRZSolver::solveStep()
{
	int cellIdx, varIdx;
	double err_newton = 1.0;
	double averPresPrev = averValue(1);
	double averPres;
	double dAverPres = 1.0;
	
	iterations = 0;
	while( err_newton > 1.e-4 && dAverPres > 1.e-4 && iterations < 8 )
	{	
		copyIterLayer();

		Solve(model->cellsNum_r+1, model->cellsNum_z+2, PRES);
		construction_from_fz(model->cellsNum_r+2, model->cellsNum_z+2, PRES);
 
		err_newton = convergance(cellIdx, varIdx);

		averPres = averValue(0);
		dAverPres = fabs(averPres - averPresPrev);
		averPresPrev = averPres;

		iterations++;
	}

	cout << "Newton Iterations = " << iterations << endl;
}

void OilRZSolver::construction_from_fz(int N, int n, int key)
{
	vector<Cell>::iterator it;
	if (key == PRES)
	{
		for(int i = 0; i < N; i++)
		{
			for(int j = 0; j < model->cellsNum_z+2; j++)
			{
				model->cells[i*(model->cellsNum_z+2) + j].u_next.p = fz[i][j+1];
			}
		}
	}
}

void OilRZSolver::LeftBoundAppr(int MZ, int key)
{
	for(int i = 0; i < MZ; i++)
	{
		for (int j = 0 ; j < MZ; j++)
		{
			C[i][j] = 0.0;
			B[i][j] = 0.0;
			A[i][j] = 0.0;
		}
		C[i][i] = 1.0;
		B[i][i] = -1.0;
		RightSide[i][0] = 0.0;
	}

	map<int,double>::iterator it;
	int idx = 0;
	if(key == PRES)
	{
		for(it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
		{
			Cell& curr = model->cells[it->first];
			Cell& nebr = model->cells[it->first+model->cellsNum_z+2];
			
			idx = it->first;

			// First eqn
			C[idx][idx] = model->solve_eqLeft_dp(it->first);
			B[idx][idx] = model->solve_eqLeft_dp_beta(it->first);
			RightSide[idx][0] = -model->solve_eqLeft(it->first) + 
								C[idx][idx] * curr.u_next.p +
								B[idx][idx] * nebr.u_next.p;
		}
	}

	construction_bz(MZ, 2);
}

void OilRZSolver::RightBoundAppr(int MZ, int key)
{
	for(int i = 0; i < MZ; i++)
	{
		for (int j = 0 ; j < MZ; j++)
		{
			C[i][j] = 0.0;
			B[i][j] = 0.0;
			A[i][j] = 0.0;
		}
		RightSide[i][0] = 0.0;
	}

	if(key == PRES)
	{
		int idx = 0;
		
		for(int i = (model->cellsNum_r+1)*(model->cellsNum_z+2); i < (model->cellsNum_r+2)*(model->cellsNum_z+2); i++)
		{
			Cell& curr = model->cells[i];
			Cell& nebr = model->cells[i-model->cellsNum_z-2];

			// First eqn
			A[idx][idx] = model->solve_eqRight_dp(i);
			B[idx][idx] = model->solve_eqRight_dp_beta(i);
			RightSide[idx][0] = -model->solve_eqRight(i) + 
								A[idx][idx] * curr.u_next.p +
								B[idx][idx] * nebr.u_next.p;

			idx++;
		}
	}

	construction_bz(MZ,1);
}

void OilRZSolver::MiddleAppr(int current, int MZ, int key)
{
	for(int i = 0; i < MZ; i++)
	{
		for(int j = 0; j < MZ; j++)
		{
			A[i][j] = 0.0;
			B[i][j] = 0.0;
			C[i][j] = 0.0;
		}
		RightSide[i][0] = 0.0;
	}

	if(key == PRES)
	{
		int idx = 0;
		int i = current * (model->cellsNum_z+2);
		
		// Top cell
		TopAppr(i, key);
		idx++;

		// Middle cells
		for(i = current * (model->cellsNum_z+2) + 1; i < (current+1) * (model->cellsNum_z+2) - 1; i++)
		{
			Var1phase& next = model->cells[i].u_next;
			
			C[idx][idx] = model->solve_eq_dp_beta(i, i - model->cellsNum_z - 2);
			B[idx][idx-1] = model->solve_eq_dp_beta(i, i-1);
			B[idx][idx] = model->solve_eq_dp(i);
			B[idx][idx+1] = model->solve_eq_dp_beta(i, i+1);
			A[idx][idx] = model->solve_eq_dp_beta(i, i + model->cellsNum_z + 2);
			RightSide[idx][0] = -model->solve_eq(i) + 
								C[idx][idx] * model->cells[i-model->cellsNum_z-2].u_next.p + 
								B[idx][idx-1] * model->cells[i-1].u_next.p +
								B[idx][idx] * model->cells[i].u_next.p + 
								B[idx][idx+1] * model->cells[i+1].u_next.p + 
								A[idx][idx] * model->cells[i+model->cellsNum_z+2].u_next.p;	
			idx++;
		}

		// Bottom cell
		BottomAppr(i, key);
	}

	construction_bz(MZ,2);
}

void OilRZSolver::TopAppr(int i, int key)
{
	if(key == PRES)
	{
		// First eqn
		//B[0][0] = 1.0;
		//B[0][1] = -1.0;

		B[0][0] = model->solve_eqTop_dp(i);
		B[0][1] = model->solve_eqTop_dp_beta(i);
		RightSide[0][0] = -model->solve_eqTop(i) +
			B[0][0] * model->cells[i].u_next.p +
			B[0][1] * model->cells[i + 1].u_next.p;
	}
}

void OilRZSolver::BottomAppr(int i, int key)
{
	if(key == PRES)
	{
		int idx = model->cellsNum_z + 1;

		// First eqn
		/*B[idx][idx] = 1.0;
		B[idx][idx-1] = -1.0;*/

		B[idx][idx] = model->solve_eqBot_dp(i);
		B[idx][idx-1] = model->solve_eqBot_dp_beta(i);
		RightSide[idx][0] = -model->solve_eqBot(i) +
			B[idx][idx] * model->cells[i].u_next.p +
			B[idx][idx-1] * model->cells[i - 1].u_next.p;
	}
}