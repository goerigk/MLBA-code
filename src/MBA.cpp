#include "MBA.h"

#include "ilcplex/ilocplex.h"
#include <iostream>
#include <algorithm>

ILOSTLBEGIN

//double TIMELIMIT=300;
double TIMELIMIT=600;

using namespace std;

void MBA::generate(int _n, int _m, double dense)
{
	n = _n;
	m = _m;
	
	//generate node costs
	c.resize(n);
	for (int i=0; i<n; ++i)
	{
		c[i].resize(m);
		for (int j=0; j<m; ++j)
			c[i][j] = rand()%100 + 1;
	}
	
	//generate graph
	//sample n^2 paths
	suc.resize(n);
	pre.resize(n);
	
	e_suc.resize(n);
	e_pre.resize(n);
	
	int e_counter = 0;
	
	for (int i=0; i<n; ++i)
	{
		suc[i].resize(m);
		pre[i].resize(m);
		
		e_suc[i].resize(m);
		e_pre[i].resize(m);
		
		for (int j=1; j<m; ++j)
		{
			suc[i][j-1].insert(i);
			pre[i][j].insert(i);
			
			Edge e;
			e.fromi = i;
			e.fromj = j-1;
			e.toi = i;
			e.toj = j;
			e.c = c[i][j];
			e.between = -1;
			edges.push_back(e);
			
			e_suc[i][j-1].insert(e_counter);
			e_pre[i][j].insert(e_counter);
			++e_counter;
		}
		
	}
	
	//int maxp = 0;
	//if (dense == 0)
	int maxp = dense*n;
	//else
		//maxp = 0.5*n*n;
	
	for (int k=0; k<maxp; ++k)
	{
		int curi = rand()%n;
			
		for (int l=1; l<m; ++l)
		{
			int nexti = rand()%n;
				
			//avoid existing edges!
			if (suc[curi][l-1].count(nexti) < 0.5)
			{
				suc[curi][l-1].insert(nexti);
				pre[nexti][l].insert(curi);
				
				Edge e;
				e.fromi = curi;
				e.fromj = l-1;
				e.toi = nexti;
				e.toj = l;
				e.c = c[nexti][l];
				e.between = -1;
				edges.push_back(e);
				
				e_suc[curi][l-1].insert(e_counter);
				e_pre[nexti][l].insert(e_counter);
				++e_counter;
			}
			
			curi = nexti;
		}
	}
	
	paths_are_generated = false;
	para_method0_L = 1;
	para_post = false;
}

void MBA::print_problem()
{
	for (int i=0; i<n; ++i)	
	{
		for (int j=0; j<m; ++j)
			cout<<c[i][j]<<" ";
		cout<<"\n";
	}
	
}


Solution MBA::solve_form1b()
{
	IloEnv env;
	IloModel model(env);
			
	vector<vector<vector<IloNumVar> > > cplexx(n);
	for (int i=0; i<n; ++i)
	{
		cplexx[i].resize(m);
		for (int j=0; j<m; ++j)
		{
			cplexx[i][j].resize(n);
			for (int k=0; k<n; ++k)
				cplexx[i][j][k] = IloNumVar(env,0,1,ILOBOOL);
		}
	}
		
	IloNumVar cplexD(env,0,IloInfinity,ILOFLOAT);
	
	for (int k=0; k<n; ++k)
		model.add(cplexx[k][0][k] == 1);
	
	for (int j=0; j<m; ++j)
		for (int k=0; k<n; ++k)
		{
			IloExpr con(env);
			for (int i=0; i<n; ++i)
				con += cplexx[i][j][k];
			model.add(con == 1);
		}
		
	for (int i=0; i<n; ++i)
		for (int j=0; j<m; ++j)
		{
			IloExpr con(env);
			for (int k=0; k<n; ++k)
				con += cplexx[i][j][k];
			model.add(con == 1);
		}
		
	for (int k=0; k<n; ++k)
	{
		IloExpr con(env);
		for (int i=0; i<n; ++i)
			for (int j=0; j<m; ++j)
				con += c[i][j] * cplexx[i][j][k];
		model.add(con <= cplexD);
	}
	
	for (int i=0; i<n; ++i)
		for (int j=0; j<m-1; ++j)
			for (int k=0; k<n; ++k)
			{
				IloExpr con(env);
				for (set<int>::iterator it=suc[i][j].begin(); it!=suc[i][j].end(); ++it)
					con += cplexx[*it][j+1][k];
				model.add(con >= cplexx[i][j][k]);
			}
			
	
	model.add(IloMinimize(env, cplexD));

	IloCplex cplex(model);	
	//warmstart
	{
		vector<vector<int> > greedy = gen_paths0_L(c,1);
		
		
		IloNumVarArray startVar(env);
		IloNumArray startVal(env);
		for (int k=0; k<n; ++k)		
			for (int j=0; j<m; ++j)
			{
				for (int i=0; i<n; ++i)
				{
					startVar.add(cplexx[i][j][k]);
					if (greedy[k][j] == i)
						startVal.add(1);
					else
						startVal.add(0);
						
				}
				
					//cout<<k<<","<<j<<","<<greedy[k][j]<<"\n";
			}
		cplex.addMIPStart(startVar, startVal);
		startVal.end();
		startVar.end();
	}

			
	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Threads, 1);
	cplex.setParam(IloCplex::TiLim, TIMELIMIT);
    //cplex.setParam(IloCplex::WorkMem, 128);
    //cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, 128);
    //cplex.setParam(IloCplex::Param::MIP::Strategy::File, 0);
    
    
	//cplex.setParam(IloCplex::Param::Emphasis::MIP,CPX_MIPEMPHASIS_FEASIBILITY);
	
	Solution sol;
	
	double start = clock();	
	//long nodelim = 	cplex.getParam(IloCplex::NodeLim);
	//cplex.setParam(IloCplex::NodeLim, 0);
	bool result = cplex.solve();	
	//if ( result )
	//{
		//sol.root = cplex.getBestObjValue();
		//cplex.setParam(IloCplex::NodeLim, nodelim);
		//cplex.setParam(IloCplex::TiLim, TIMELIMIT - (clock() - start)/CLOCKS_PER_SEC);
		//result = cplex.solve();
	//}
	double time = (clock() - start)/CLOCKS_PER_SEC;
	
	//cout<<n<<";"<<m<<flush;
	//cout<<";"<<time<<";"<<cplex.getObjValue()<<"\n";

	if (result)
	{
		sol.time = time;
		sol.obj = cplex.getObjValue();
		sol.status = cplex.getStatus();
		sol.gap = cplex.getMIPRelativeGap();
		sol.vars = cplex.getNcols();
		sol.cons = cplex.getNrows();
		sol.lb = cplex.getBestObjValue();
	}
	else
	{
		sol.time = TIMELIMIT;
		sol.obj = -1;
		sol.status = -1;
		sol.gap = -1;
	}
	
	//print solution
	//for (int k=0; k<n; ++k)
	//{
		//cout<<"Tour "<<k<<": ";
		//for (int i=0; i<n; ++i)
			//for (int j=0; j<m; ++j)	
				//if (cplex.getValue(cplexx[i][j][k]) > 0.5)
					//cout<<"("<<i<<","<<j<<") ";
		//cout<<"\n";
		//cout<<"Cost: ";
		//for (int i=0; i<n; ++i)
			//for (int j=0; j<m; ++j)	
				//if (cplex.getValue(cplexx[i][j][k]) > 0.5)
					//cout<<c[i][j]<<" + ";
		//cout<<"\n";
	//}
	
	env.end();
	
	return sol;
}

Solution MBA::solve_form4_old()
{
	IloEnv env;
	IloModel model(env);
	
	vector<vector<IloNumVar> > cplexx(n);
	for (int i=0; i<n; ++i)
	{
		cplexx[i].resize(m);
		for (int j=0; j<m; ++j)
			cplexx[i][j] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
	}
	
	vector<vector<vector<vector<IloNumVar> > > > cplexz(n);
	for (int i=0; i<n; ++i)
	{
		cplexz[i].resize(m);
		for (int j=0; j<m; ++j)
		{
			cplexz[i][j].resize(n);
			for (int k=0; k<n; ++k)
			{
				cplexz[i][j][k].resize(n);
				for (int l=0; l<m; ++l)
					cplexz[i][j][k][l] = IloNumVar(env, 0, 1, ILOBOOL);
			}
		}
	}
	
	IloNumVar cplexD(env,0,IloInfinity,ILOFLOAT);
	
	for (int j=0; j<m; ++j)
		for (int l=0; l<m; ++l)
			if (j != l)
				for (int k=0; k<n; ++k)
				{
					IloExpr con(env);
					for (int i=0; i<n; ++i)
						con += cplexz[i][j][k][l];
					model.add(con == 1);
				}
				
	for (int i=0; i<n; ++i)
		for (int j=0; j<m; ++j)
		{
			IloExpr con(env);
			for (int k=0; k<n; ++k)
				for (int l=0; l<m; ++l)
					if (l!=j)
						con += c[k][l] * cplexz[i][j][k][l];
			model.add(cplexx[i][j] >= c[i][j] + con);
		}
	
	for (int i=0; i<n; ++i)
		for (int j=0; j<m; ++j)
			model.add(cplexD >= cplexx[i][j]);


	model.add(IloMinimize(env, cplexD));
			
	IloCplex cplex(model);

	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Threads, 1);
	cplex.setParam(IloCplex::TiLim, TIMELIMIT);

	double start = clock();
	bool result = cplex.solve();			
	double time = (clock() - start)/CLOCKS_PER_SEC;

	Solution sol;
	sol.time = time;
	sol.obj = cplex.getObjValue();
	
	env.end();
	
	return sol;
	
}

Solution MBA::solve_form4()
{
	double bigM = 0;
	for (int j=0; j<m; ++j)
	{
		double tmax = 0;
		for (int i=0; i<n; ++i)
			if (tmax < c[i][j])
				tmax = c[i][j];
		bigM += tmax;
	}
	
	IloEnv env;
	IloModel model(env);
	
	vector<vector<IloNumVar> > cplexD(n);
	for (int i=0; i<n; ++i)
	{
		cplexD[i].resize(m);
		for (int j=0; j<m; ++j)
			cplexD[i][j] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
	}
	
	IloNumVar cplexobjD(env,0,IloInfinity,ILOFLOAT);
	
	vector<vector<vector<IloNumVar> > > cplexalpha(n);
	for (int k=0; k<n; ++k)
	{
		cplexalpha[k].resize(m-1);
		for (int j=0; j<m-1; ++j)
		{
			cplexalpha[k][j].resize(n);
			for (int i=0; i<n; ++i)
				cplexalpha[k][j][i] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
		}
	}
	
	vector<vector<vector<IloNumVar> > > cplexz(n);
	for (int i=0; i<n; ++i)
	{
		cplexz[i].resize(m-1);
		for (int j=0; j<m-1; ++j)
		{
			cplexz[i][j].resize(n);
			for (int k=0; k<n; ++k)			
				cplexz[i][j][k] = IloNumVar(env, 0, 1, ILOBOOL);
		}
	}
	
	
	for (int i=0; i<n; ++i)
		for (int j=0; j<m-1; ++j)
		{
			IloExpr con(env);
			for (int k=0; k<n; ++k)
				con += cplexz[i][j][k];
			model.add(con == 1);
		}
	
	
	for (int k=0; k<n; ++k)
		for (int j=0; j<m-1; ++j)
		{
			IloExpr con(env);
			for (int i=0; i<n; ++i)
				con += cplexz[i][j][k];
			model.add(con == 1);
		}

	for (int i=0; i<n; ++i)
		model.add(cplexD[i][0] == c[i][0]);
	
	for (int i=0; i<n; ++i)
		for (int j=1; j<m; ++j)
		{
			IloExpr con(env);
			for (int k=0; k<n; ++k)
				con += cplexalpha[k][j-1][i];
			model.add(cplexD[i][j] == c[i][j] + con);
		}
		
	for (int i=0; i<n; ++i)
		for (int j=0; j<m-1; ++j)
			for (int k=0; k<n; ++k)
				model.add(cplexalpha[k][j][i] >= cplexD[k][j] - bigM*(1-cplexz[k][j][i]));
				
	for (int i=0; i<n; ++i)
		model.add(cplexobjD >= cplexD[i][m-1]);

	model.add(IloMinimize(env, cplexobjD));
			
	IloCplex cplex(model);

	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Threads, 1);
	cplex.setParam(IloCplex::TiLim, TIMELIMIT);

	double start = clock();
	bool result = cplex.solve();			
	double time = (clock() - start)/CLOCKS_PER_SEC;

	Solution sol;
	if (result)
	{
		sol.time = time;
		sol.obj = cplex.getObjValue();
		sol.status = cplex.getStatus();
	}
	else
	{
		sol.time = TIMELIMIT;
		sol.obj = -1;
		sol.status = -1;
	}
	
	env.end();
	
	return sol;
	
}


Solution MBA::solve_form3()
{
	IloEnv env;
	IloModel model(env);
	
	IloNumVar cplexD(env,0,IloInfinity,ILOFLOAT);
	
	vector<vector<vector<IloNumVar> > > cplexz(n);
	for (int i=0; i<n; ++i)
	{
		cplexz[i].resize(m-1);
		for (int j=0; j<m-1; ++j)
		{
			cplexz[i][j].resize(n);
			for (int k=0; k<n; ++k)			
				cplexz[i][j][k] = IloNumVar(env, 0, 1, ILOBOOL);
		}
	}
	
	vector<vector<vector<IloNumVar> > > cplexx(n);
	for (int i=0; i<n; ++i)
	{
		cplexx[i].resize(n);
		for (int k=0; k<n; ++k)			
		{
			cplexx[i][k].resize(m-1);
			//always add + 1 for j index to get 2,...,m
			for (int j=0; j<m-1; ++j)
				cplexx[i][k][j] = IloNumVar(env, 0, 1, ILOBOOL);
		}
	}
	
	
	for (int i=0; i<n; ++i)
		for (int j=0; j<m-1; ++j)
		{
			IloExpr con(env);
			for (int k=0; k<n; ++k)
				con += cplexz[i][j][k];
			model.add(con == 1);
		}
	
	
	for (int k=0; k<n; ++k)
		for (int j=0; j<m-1; ++j)
		{
			IloExpr con(env);
			for (int i=0; i<n; ++i)
				con += cplexz[i][j][k];
			model.add(con == 1);
		}

	for (int i=0; i<n; ++i)
		for (int k=0; k<n; ++k)
			for (int l=0; l<n; ++l)
				model.add(cplexz[i][0][k] + cplexz[k][1][l] <= cplexx[i][l][0] + 1);
	
	for (int i=0; i<n; ++i)
		for (int k=0; k<n; ++k)
			for (int l=0; l<n; ++l)
				for (int j=1; j<m-1; ++j)
					model.add(cplexx[i][k][j-1] + cplexz[k][j][l] <= cplexx[i][l][j] + 1);
					
	for (int i=0; i<n; ++i)
	{
		IloExpr con(env);
		for (int k=0; k<n; ++k)
			for (int j=0; j<m-1; ++j)
				con += c[k][j+1]*cplexx[i][k][j];
		model.add(cplexD >= c[i][0] + con);
	}				
				
	
	model.add(IloMinimize(env, cplexD));
			
	IloCplex cplex(model);

	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Threads, 1);
	cplex.setParam(IloCplex::TiLim, TIMELIMIT);

	double start = clock();
	bool result = cplex.solve();			
	double time = (clock() - start)/CLOCKS_PER_SEC;

	Solution sol;
	if (result)
	{
		sol.time = time;
		sol.obj = cplex.getObjValue();
		sol.status = cplex.getStatus();
	}
	else
	{
		sol.time = TIMELIMIT;
		sol.obj = -1;
		sol.status = -1;
	}
	
	env.end();
	
	return sol;
	
}

vector<vector<int> > MBA::gen_paths0_L(std::vector<std::vector<double> > pcost, int L)
{
	//greedy pricing with L lookahead
	// L>=1 (nr of layers in greedy + 1)
	
	vector<double> costs(n);
	vector<vector<int> > npaths(n);
	vector<int> curi(n);
	
	for (int i=0; i<n; ++i)
	{
		costs[i] = pcost[i][0];
		curi[i] = i;
		npaths[i].push_back(i);
	}	
	
	for (int j=0; j<m-1; ++j)
	{		
		if (L > m-0.5-j)
			L = m-1-j;
		
		IloEnv env;
		IloModel model(env);
	
		IloNumVar cplexD(env,0,IloInfinity,ILOFLOAT);
	
		vector<vector<vector<IloNumVar> > > cplexx(n);
		for (int i=0; i<n; ++i)
		{
			cplexx[i].resize(L+1);
			for (int l=0; l<L+1; ++l)
			{
				cplexx[i][l].resize(n);
				for (int k=0; k<n; ++k)
					cplexx[i][l][k] = IloNumVar(env, 0, 1, ILOBOOL);
			}
		}
		
		
		for (int k=0; k<n; ++k)
			model.add(cplexx[k][0][k] == 1);
			
	
		for (int l=1; l<L+1; ++l)
			for (int k=0; k<n; ++k)
			{
				IloExpr con(env);
				for (int i=0; i<n; ++i)
					con += cplexx[i][l][k];
				model.add(con == 1);
			}
		
		for (int i=0; i<n; ++i)
			for (int l=1; l<L+1; ++l)
			{
				IloExpr con(env);
				for (int k=0; k<n; ++k)
					con += cplexx[i][l][k];
				model.add(con == 1);
			}
		
		for (int k=0; k<n; ++k)
		{
			IloExpr con(env);
			for (int i=0; i<n; ++i)
				for (int l=1; l<L+1; ++l)
					con += pcost[i][j+l] * cplexx[i][l][k];
			model.add(costs[k] + con <= cplexD);
		}
		
		
		for (int i=0; i<n; ++i)
			for (int l=0; l<L; ++l)
				for (int k=0; k<n; ++k)
				{
					IloExpr con(env);
					for (set<int>::iterator it=suc[i][j+l].begin(); it!=suc[i][j+l].end(); ++it)
						con += cplexx[*it][l+1][k];
					model.add(con >= cplexx[i][l][k]);
				}

		
		model.add(IloMinimize(env, cplexD));
		
		IloCplex cplex(model);

		cplex.setOut(env.getNullStream());
		cplex.setParam(IloCplex::Threads, 1);

		bool result = cplex.solve();			

		vector<bool> moved(n,false);
		vector<double> ncosts(n);

		for (int i=0; i<n; ++i)
			for (int k=0; k<n; ++k)
				if (cplex.getValue(cplexx[i][1][k]) > 0.5)
				{
					costs[k] += pcost[i][j+1];
					ncosts[i] = costs[k];
					for (int l=0; l<n; ++l)
						if (curi[l] == k && !moved[l])
						{
							moved[l] = true;
							curi[l] = i;
							npaths[l].push_back(i);
						}
				}
	
		for (int i=0; i<n; ++i)
			costs[i] = ncosts[i];
			
		env.end();	
	}
    
    double maxcost = 0;
    for (int i=0; i<n; ++i)
        maxcost = max(maxcost, costs[i]);
    cout<<"Greedy;"<<maxcost<<"\n";
	
	if (para_post)
		return post_improve(pcost,npaths);
	
	return npaths;	
}


vector<vector<int> > MBA::gen_paths0(vector<vector<double> > pcost)
{
	cout<<"WARNING: USING gen_paths0\n"<<flush;
	
	//greedy pricing
	
	vector<double> costs(n);
	vector<vector<int> > npaths(n);
	vector<int> curi(n);
	
	for (int i=0; i<n; ++i)
	{
		costs[i] = pcost[i][0];
		curi[i] = i;
		npaths[i].push_back(i);
	}	
	
		
	for (int j=0; j<m-1; ++j)
	{		
		IloEnv env;
		IloModel model(env);
	
		IloNumVar cplexD(env,0,IloInfinity,ILOFLOAT);
	
		vector<vector<IloNumVar> > cplexx(n);
		for (int i=0; i<n; ++i)
		{
			cplexx[i].resize(n);
			for (int k=0; k<n; ++k)
				if (suc[i][j].count(k) > 0.5)
					cplexx[i][k] = IloNumVar(env, 0, 1, ILOBOOL);
				else
					cplexx[i][k] = IloNumVar(env, 0, 0, ILOFLOAT);
		}
		
		for (int i=0; i<n; ++i)
		{
			IloExpr con(env);
			for (int k=0; k<n; ++k)
				con += cplexx[i][k];
			model.add(con == 1);
		}
		
		for (int k=0; k<n; ++k)
		{
			IloExpr con(env);
			for (int i=0; i<n; ++i)
				con += cplexx[i][k];
			model.add(con == 1);
		}
		
		for (int i=0; i<n; ++i)
		{
			IloExpr con(env);
			for (int k=0; k<n; ++k)
				con += pcost[k][j+1]*cplexx[i][k];
			model.add(cplexD >= costs[i] + con);
		}
		
		model.add(IloMinimize(env, cplexD));
		
		IloCplex cplex(model);

		cplex.setOut(env.getNullStream());
		cplex.setParam(IloCplex::Threads, 1);

		bool result = cplex.solve();			

		vector<bool> moved(n,false);
		vector<double> ncosts(n);

		for (int i=0; i<n; ++i)
			for (int k=0; k<n; ++k)
				if (cplex.getValue(cplexx[i][k]) > 0.5)
				{
					costs[i] += pcost[k][j+1];
					ncosts[k] = costs[i];
					for (int l=0; l<n; ++l)
						if (curi[l] == i && !moved[l])
						{
							moved[l] = true;
							curi[l] = k;
							npaths[l].push_back(k);
						}
				}
	
		for (int i=0; i<n; ++i)
			costs[i] = ncosts[i];
			
		env.end();	
	}
	
	return npaths;	
}

Solution MBA::solve_form0(int method)
{
	//start with greedy solution
	double start = clock();	
	vector<vector<int> > paths = gen_paths0_L(c,para_method0_L);
		
	bool imp = true;
	int oldtime = 0;
	
	do
	{
		imp = false;
		int P = paths.size();

		//calculate path lenghts
		vector<double> plen(P,0);
		for (int p=0; p<P; ++p)
			for (int j=0; j<m; ++j)
				plen[p] += c[paths[p][j]][j];
		
		//solve dual LP
		IloEnv env;
		IloModel model(env);
		
		vector<vector<IloNumVar> > cplexu(n);
		for (int i=0; i<n; ++i)
		{
			cplexu[i].resize(m);
			for (int j=0; j<m; ++j)
				cplexu[i][j] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
		}
		
		vector<IloNumVar> cplexr(n);
		for (int k=0; k<n; ++k)
			cplexr[k] = IloNumVar(env, 0, 1, ILOFLOAT);
			
		for (int k=0; k<n; ++k)
			for (int p=0; p<P; ++p)
			{
				IloExpr con(env);
				for (int j=0; j<m; ++j)
					con += cplexu[paths[p][j]][j];
				model.add(con <= plen[p]*cplexr[k]);
			}
			
		{
			IloExpr con(env);
			for (int k=0; k<n; ++k)
				con += cplexr[k];
			model.add(con <= 1);
		}

		IloExpr obj(env);
		for (int i=0; i<n; ++i)
			for (int j=0; j<m; ++j)
				obj += cplexu[i][j];
		
		model.add(IloMaximize(env, obj));		
		
		IloCplex cplex(model);

		cplex.setOut(env.getNullStream());
		cplex.setParam(IloCplex::Threads, 1);
		//cplex.setParam(IloCplex::TiLim, TIMELIMIT);
		
		bool result = cplex.solve();			
				
		//solve pricing problem
	
		vector<vector<double> > optu(n);
		for (int i=0; i<n; ++i)
		{
			optu[i].resize(m);
			for (int j=0; j<m; ++j)
				optu[i][j] = cplex.getValue(cplexu[i][j]);
		}
		vector<double> optr(n);
		for (int k=0; k<n; ++k)
			optr[k] = cplex.getValue(cplexr[k]);
		
		env.end();
		
		//find min k
		int mink=0;
		for (int k=0; k<n; ++k)
			if (optr[mink] > optr[k])
				mink = k;
	
		vector<vector<double> > redcost(n);
		for (int i=0; i<n; ++i)
		{
			redcost[i].resize(m);
			for (int j=0; j<m; ++j)
				redcost[i][j] = optr[mink]*c[i][j]-optu[i][j];
		}
	
		//find pricing paths
		vector<vector<int> > npaths;
		if (method == 0)
			npaths = gen_paths0_L(redcost,para_method0_L);
		else if (method == 1)
			npaths = gen_paths1(redcost);
		else if (method == 2)
			npaths = gen_paths2(redcost);
		else if (method == 3)
			npaths = gen_paths3(redcost);
		
		//add paths
		set<vector<int> > spaths;
		for (int p=0; p<P; ++p)
			spaths.insert(paths[p]);
		for (int p=0; p<npaths.size(); ++p)	
			spaths.insert(npaths[p]);
		paths.clear();
		for (set<vector<int> >::iterator it=spaths.begin(); it!=spaths.end(); ++it)
			paths.push_back(*it);
		if (paths.size() > P)
			imp = true;
	
		double ctime = (clock() - start)/CLOCKS_PER_SEC;
		if (((int)ctime)/60 > oldtime)
		{
			oldtime = ((int)ctime)/60;
			cout<<"P;"<<ctime<<";"<<paths.size()<<"\n"<<flush;
		}
	
	
	}while(imp && (clock() - start)/CLOCKS_PER_SEC < TIMELIMIT);
	
	if (method < 2)	
		cout<<"pre;"<<method<<";"<<(clock() - start)/CLOCKS_PER_SEC<<";";
	else if (method == 2)
		cout<<"pre;"<<method<<";"<<para_method2_nr<<";"<<(clock() - start)/CLOCKS_PER_SEC<<";";
	else if (method == 3)
		cout<<"pre;"<<method<<";"<<para_method3_thresh<<";"<<(clock() - start)/CLOCKS_PER_SEC<<";";
		
	cout<<paths.size()<<"\n"<<flush;
	
	
	return solve_form0_master(paths);
}

Solution MBA::solve_form0_master(vector<vector<int> > paths)
{
	//paths have been generated
	
	int P = paths.size();
	
	//calculate path lenghts
	vector<double> plen(P,0);
	for (int p=0; p<P; ++p)
		for (int j=0; j<m; ++j)
			plen[p] += c[paths[p][j]][j];
			
	IloEnv env;
	IloModel model(env);
	
	IloNumVar cplexD(env,0,IloInfinity,ILOFLOAT);
	
	vector<IloNumVar> cplexx(P);
	for (int p=0; p<P; ++p)
		cplexx[p] = IloNumVar(env, 0, 1, ILOBOOL);
	
	for (int i=0; i<n; ++i)
		for (int j=0; j<m; ++j)
		{
			IloExpr con(env);
			for (int p=0; p<P; ++p)
				if (paths[p][j] == i)
					con += cplexx[p];
			model.add(con == 1);
		}
		
	for (int p=0; p<P; ++p)
		model.add(plen[p]*cplexx[p] <= cplexD);
	
	model.add(IloMinimize(env, cplexD));
		
	IloCplex cplex(model);

	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Threads, 1);
	//cplex.setParam(IloCplex::TiLim, TIMELIMIT);

	double start = clock();
	bool result = cplex.solve();			
	double time = (clock() - start)/CLOCKS_PER_SEC;

	Solution sol;
	if (result)
	{
		sol.time = time;
		sol.obj = cplex.getObjValue();
		sol.status = cplex.getStatus();
	}
	else
	{
		sol.time = TIMELIMIT;
		sol.obj = -1;
		sol.status = -1;
	}

	env.end();
	
	return sol;
	
}

Solution MBA::solve_greedy()
{
	double start = clock();	
	
	Solution sol;
	sol.status = -1;
		
	vector<vector<int> > paths = gen_paths0(c);
	vector<int> costs(n,0);
	for (int i=0; i<n; ++i)
		for (int j=0; j<m; ++j)
			costs[i] += c[paths[i][j]][j];
			
	sort(costs.rbegin(),costs.rend());
		
	sol.obj = costs[0];
	
	sol.time = (clock() - start)/CLOCKS_PER_SEC;
	
	return sol;
}

Solution MBA::solve_greedy_L(int L)
{
	double start = clock();	
	
	Solution sol;
	sol.status = -1;
	
	vector<vector<int> > paths = gen_paths0_L(c, L);
	vector<int> costs(n,0);
	for (int i=0; i<n; ++i)
		for (int j=0; j<m; ++j)
			costs[i] += c[paths[i][j]][j];
			
	sort(costs.rbegin(),costs.rend());
		
	sol.obj = costs[0];
	
	sol.time = (clock() - start)/CLOCKS_PER_SEC;
	
	return sol;
}

void MBA::create_all_paths()
{
	allpaths.clear();
	
	//enumerate all paths
	set<vector<int> > partials;
	for (int i=0; i<n; ++i)
	{
		vector<int> path;
		path.push_back(i);
		partials.insert(path);
	}
	
	for (int j=0; j<m-1; ++j)
	{
		set<vector<int> > npartials;
		
		for (set<vector<int> >::iterator it=partials.begin(); it!=partials.end(); ++it)
			for (set<int>::iterator itsuc=suc[(*it)[j]][j].begin(); itsuc!=suc[(*it)[j]][j].end(); ++itsuc)
			{
				vector<int> npath = *it;
				npath.push_back(*itsuc);
				npartials.insert(npath);
			}
			
		partials = npartials;
	}
	
	for (set<vector<int> >::iterator it=partials.begin(); it!=partials.end(); ++it)
		allpaths.push_back(*it);
	
	cout<<"Numer of possible paths: "<<partials.size()<<"\n";
}

vector<vector<int> > MBA::gen_paths1(vector<vector<double> > pcost)
{
	//full enumeration pricing
	
	if (!paths_are_generated)
	{
		create_all_paths();
		paths_are_generated = true;
	}
	
	vector<vector<int> > paths;
	int P = allpaths.size();
	for (int p=0; p<P; ++p)
	{
		double curcost = 0;
		for (int j=0; j<m; ++j)
			curcost += pcost[allpaths[p][j]][j];
		if (curcost < -0.001)
			paths.push_back(allpaths[p]);
	}
	
	return paths;
	
}

vector<vector<int> > MBA::gen_paths2(vector<vector<double> > pcost)
{
	//restricted enumeration pricing
	//add only nr many paths
	
	if (!paths_are_generated)
	{
		create_all_paths();
		paths_are_generated = true;
	}
	
	vector<vector<int> > paths;
	random_shuffle(allpaths.begin(), allpaths.end());
	int P = allpaths.size();
	int found = 0;
	for (int p=0; p<P && found < para_method2_nr; ++p)
	{
		double curcost = 0;
		for (int j=0; j<m; ++j)
			curcost += pcost[allpaths[p][j]][j];
		if (curcost < -0.001)
		{
			paths.push_back(allpaths[p]);
			++found;
		}
	}
	
	return paths;
	
}

void MBA::set_para_method2(int _nr)
{
	para_method2_nr = _nr;
}

void MBA::reset()
{
	paths_are_generated = false;
	allpaths.clear();
}


vector<vector<int> > MBA::gen_paths3(vector<vector<double> > pcost)
{
	//solve relaxation with current costs
	//then decompose into paths
	IloEnv env;
	IloModel model(env);
			
	vector<vector<vector<IloNumVar> > > cplexx(n);
	for (int i=0; i<n; ++i)
	{
		cplexx[i].resize(m);
		for (int j=0; j<m; ++j)
		{
			cplexx[i][j].resize(n);
			for (int k=0; k<n; ++k)
				cplexx[i][j][k] = IloNumVar(env,0,1,ILOFLOAT);
		}
	}
		
	IloNumVar cplexD(env,0,IloInfinity,ILOFLOAT);
	
	
	for (int j=0; j<m; ++j)
		for (int k=0; k<n; ++k)
		{
			IloExpr con(env);
			for (int i=0; i<n; ++i)
				con += cplexx[i][j][k];
			model.add(con == 1);
		}
		
	for (int i=0; i<n; ++i)
		for (int j=0; j<m; ++j)
		{
			IloExpr con(env);
			for (int k=0; k<n; ++k)
				con += cplexx[i][j][k];
			model.add(con == 1);
		}
		
	for (int k=0; k<n; ++k)
	{
		IloExpr con(env);
		for (int i=0; i<n; ++i)
			for (int j=0; j<m; ++j)
				con += pcost[i][j] * cplexx[i][j][k];
		model.add(con <= cplexD);
	}
	
	for (int i=0; i<n; ++i)
		for (int j=0; j<m-1; ++j)
			for (int k=0; k<n; ++k)
			{
				IloExpr con(env);
				for (set<int>::iterator it=suc[i][j].begin(); it!=suc[i][j].end(); ++it)
					con += cplexx[*it][j+1][k];
				model.add(con >= cplexx[i][j][k]);
			}
			
	for (int i=0; i<n; ++i)
		for (int j=1; j<m; ++j)
			for (int k=0; k<n; ++k)
			{
				IloExpr con(env);
				for (set<int>::iterator it=pre[i][j].begin(); it!=pre[i][j].end(); ++it)
					con += cplexx[*it][j-1][k];
				model.add(con >= cplexx[i][j][k]);
			}
	
	model.add(IloMinimize(env, cplexD));
			
	IloCplex cplex(model);

	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Threads, 1);

	bool result = cplex.solve();			


	//decompose into paths
	vector<vector<int> > paths;
	for (int k=0; k<n; ++k)
	{
		vector<vector<double> > flow(n);
		for (int i=0; i<n; ++i)
		{
			flow[i].resize(m);
			for (int j=0; j<m; ++j)
				if (cplex.getValue(cplexx[i][j][k]) > para_method3_thresh)
					flow[i][j] = cplex.getValue(cplexx[i][j][k]);
		}
		
		//for (int i=0; i<n; ++i)
			//for (int j=0; j<m; ++j)
				//if (flow[i][j] > 0)
					//cout<<i<<","<<j<<","<<flow[i][j]<<"\n";
		
		vector<int> order(n);
		for (int i=0; i<n; ++i)
			order[i] = i;
		random_shuffle(order.begin(), order.end());
		
		bool usedup = false;
		
		do
		{
			vector<int> path;
			bool lookout = true;
			int cur = -1;
			
			for (int i=0; lookout && i<n; ++i)
				if (flow[order[i]][0] > 0)
				{
					lookout = false;
					path.push_back(order[i]);
					cur = order[i];
					//cout<<cur<<"\n"<<flush;
				}
			
			if (lookout)
				break;
			
			for (int j=0; j<m-1; ++j)
			{
				//assert(!lookout);
				if (lookout)
					break;
				lookout = true;
				for (int i=0; lookout && i<n; ++i)
					if (suc[cur][j].count(order[i]) > 0.5 && flow[order[i]][j+1] > 0)				
					{
						lookout = false;
						path.push_back(order[i]);
						cur = order[i];
						//cout<<cur<<"\n"<<flush;
					}					
			}
			
			//assert(path.size() == m);
			if (path.size() < m)
				break;
				
			paths.push_back(path);
			
			//reduce flow along path
			double minval = 1;
			for (int j=0; j<m; ++j)
				if (flow[path[j]][j] < minval)
					minval = flow[path[j]][j];
			for (int j=0; j<m; ++j)
				flow[path[j]][j] -= minval;
			
			double gmaxval = 0;
			for (int i=0; i<n; ++i)
				for (int j=0; j<m; ++j)
					if (gmaxval < flow[i][j])
						gmaxval = flow[i][j];
			if (gmaxval < 0.01)
				usedup = true;
				
		}
		while(!usedup);
	}
	
	env.end();

	return paths;
	
}

Solution MBA::solve_form1c()
{
	IloEnv env;
	IloModel model(env);
			
	vector<vector<vector<IloNumVar> > > cplexx(n);
	for (int i=0; i<n; ++i)
	{
		cplexx[i].resize(m);
		for (int j=0; j<m; ++j)
		{
			cplexx[i][j].resize(n);
			for (int k=0; k<n; ++k)
				cplexx[i][j][k] = IloNumVar(env,0,1,ILOBOOL);
		}
	}
		
	IloNumVar cplexD(env,0,IloInfinity,ILOFLOAT);
	
	
	for (int j=0; j<m; ++j)
		for (int k=0; k<n; ++k)
		{
			IloExpr con(env);
			for (int i=0; i<n; ++i)
				con += cplexx[i][j][k];
			model.add(con == 1);
		}
		
	for (int i=0; i<n; ++i)
		for (int j=0; j<m; ++j)
		{
			IloExpr con(env);
			for (int k=0; k<n; ++k)
				con += cplexx[i][j][k];
			model.add(con == 1);
		}
		
	for (int k=0; k<n; ++k)
	{
		IloExpr con(env);
		for (int i=0; i<n; ++i)
			for (int j=0; j<m; ++j)
				con += c[i][j] * cplexx[i][j][k];
		model.add(con <= cplexD);
	}
	
	for (int i=0; i<n; ++i)
		for (int j=0; j<m-1; ++j)
			for (int k=0; k<n; ++k)
			{
				IloExpr con(env);
				for (set<int>::iterator it=suc[i][j].begin(); it!=suc[i][j].end(); ++it)
					con += cplexx[*it][j+1][k];
				model.add(con >= cplexx[i][j][k]);
			}
			
	
	for (int i=0; i<n; ++i)
		for (int j=1; j<m; ++j)
			for (int k=0; k<n; ++k)
			{
				IloExpr con(env);
				for (set<int>::iterator it=pre[i][j].begin(); it!=pre[i][j].end(); ++it)
					con += cplexx[*it][j-1][k];
				model.add(con >= cplexx[i][j][k]);
			}
	
	model.add(IloMinimize(env, cplexD));
			
	IloCplex cplex(model);

	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Threads, 1);
	cplex.setParam(IloCplex::TiLim, TIMELIMIT);

	double start = clock();
	bool result = cplex.solve();			
	double time = (clock() - start)/CLOCKS_PER_SEC;
	
	//cout<<n<<";"<<m<<flush;
	//cout<<";"<<time<<";"<<cplex.getObjValue()<<"\n";

	Solution sol;
	if (result)
	{
		sol.time = time;
		sol.obj = cplex.getObjValue();
		sol.status = cplex.getStatus();
	}
	else
	{
		sol.time = TIMELIMIT;
		sol.obj = -1;
		sol.status = -1;
	}
	
	//print solution
	//for (int k=0; k<n; ++k)
	//{
		//cout<<"Tour "<<k<<": ";
		//for (int i=0; i<n; ++i)
			//for (int j=0; j<m; ++j)	
				//if (cplex.getValue(cplexx[i][j][k]) > 0.5)
					//cout<<"("<<i<<","<<j<<") ";
		//cout<<"\n";
		//cout<<"Cost: ";
		//for (int i=0; i<n; ++i)
			//for (int j=0; j<m; ++j)	
				//if (cplex.getValue(cplexx[i][j][k]) > 0.5)
					//cout<<c[i][j]<<" + ";
		//cout<<"\n";
	//}
	
	env.end();
	
	return sol;
}


void MBA::set_para_method3(double _thresh)
{
	para_method3_thresh = _thresh;
}

void MBA::set_para_method0(int _L)
{
	para_method0_L = _L;
}

void MBA::set_para_post(bool _post)
{
	para_post = _post;
}

vector<vector<int> > MBA::post_improve(vector<vector<double> > pcost, vector<vector<int> > current)
{
	//try to improve current solution locally
	
	bool improved;
	
	do
	{
		//cout<<"post round\n"<<flush;
		
		improved = false;
	
		for (int j=0; j<m-1; ++j)
		{	
			//calculate the costs in each node
			vector<double> lcost(n);
			for (int i=0; i<n; ++i)
			{
				int k = current[i][j];
				for (int jp=0; jp<j+1; ++jp)
					lcost[k] += pcost[current[i][jp]][jp];
			}
			
			vector<double> rcost(n);
			for (int i=0; i<n; ++i)
			{
				int k = current[i][j+1];
				for (int jp=j+1; jp<m; ++jp)
					rcost[k] += pcost[current[i][jp]][jp];
			}
			
			IloEnv env;
			IloModel model(env);

			IloNumVar cplexD(env,0,IloInfinity,ILOFLOAT);

			vector<vector<IloNumVar> > cplexx(n);
			for (int i=0; i<n; ++i)
			{
				cplexx[i].resize(n);
				for (int k=0; k<n; ++k)
					if (suc[i][j].count(k) > 0.5)
						cplexx[i][k] = IloNumVar(env, 0, 1, ILOBOOL);
					else
						cplexx[i][k] = IloNumVar(env, 0, 0, ILOFLOAT);
			}
			
			for (int i=0; i<n; ++i)
			{
				IloExpr con(env);
				for (int k=0; k<n; ++k)
					con += cplexx[i][k];
				model.add(con == 1);
			}
			
			for (int k=0; k<n; ++k)
			{
				IloExpr con(env);
				for (int i=0; i<n; ++i)
					con += cplexx[i][k];
				model.add(con == 1);
			}
			
			for (int i=0; i<n; ++i)
			{
				IloExpr con(env);
				for (int k=0; k<n; ++k)
					con += (lcost[i] + rcost[k])*cplexx[i][k];
				model.add(cplexD >= con);
			}
			
			model.add(IloMinimize(env, cplexD));
			
			IloCplex cplex(model);

			cplex.setOut(env.getNullStream());
			cplex.setParam(IloCplex::Threads, 1);

			bool result = cplex.solve();			

			double oldcost = 0;
			for (int i=0; i<n; ++i)
			{
				double tmax = 0;
				for (int jp=0; jp<m; ++jp)
					tmax += pcost[current[i][jp]][jp];
				if (tmax > oldcost)
					oldcost = tmax;
			}
			
			if (cplex.getObjValue() < oldcost - 0.01)
			{
				//cout<<oldcost<<";"<<cplex.getObjValue()<<"\n"<<flush;
				
				improved = true;
				//update current paths
				vector<vector<int> > npaths;
				
				vector<int> left(n);
				vector<int> right(n);
				
				for (int i=0; i<n; ++i)
				{
					left[current[i][j]] = i;
					right[current[i][j+1]] = i;
				}
				
				for (int i=0; i<n; ++i)
					for (int k=0; k<n; ++k)
						if (cplex.getValue(cplexx[i][k]) > 0.5)
						{
							vector<int> npath;
							for (int jp=0; jp<j+1; ++jp)
								npath.push_back(current[left[i]][jp]);
							for (int jp=j+1; jp<m; ++jp)
								npath.push_back(current[right[k]][jp]);
							npaths.push_back(npath);							
						}
				current = npaths;
			}


			env.end();	
		}
		
	}while(improved);
	
	return current;		
}




Solution MBA::solve_sparse()
{
	int E = edges.size();
	
	IloEnv env;
	IloModel model(env);
			
	vector<vector<IloNumVar> > cplexx(E);
	for (int e=0; e<E; ++e)
	{
		cplexx[e].resize(n);
		for (int k=0; k<n; ++k)
			cplexx[e][k] = IloNumVar(env,0,1,ILOBOOL);
	}
		
	IloNumVar cplexD(env,0,IloInfinity,ILOFLOAT);
	
	for (int k=0; k<n; ++k)
	{
		IloExpr con(env);
		for (set<int>::iterator it=e_suc[k][0].begin(); it!=e_suc[k][0].end(); ++it)
			con += cplexx[*it][k];
		model.add(con == 1);
	}
	
	for (int k=0; k<n; ++k)
	{
		IloExpr con(env);
		for (int e=0; e<E; ++e)
			con += edges[e].c * cplexx[e][k];
		con += c[k][0];
		model.add(con <= cplexD);
	}
	
	for (int i=0; i<n; ++i)
		for (int j=1; j<m-1; ++j)
			for (int k=0; k<n; ++k)
			{
				IloExpr con(env);
				for (set<int>::iterator it=e_suc[i][j].begin(); it!=e_suc[i][j].end(); ++it)
					con += cplexx[*it][k];
				for (set<int>::iterator it=e_pre[i][j].begin(); it!=e_pre[i][j].end(); ++it)
					con -= cplexx[*it][k];
				model.add(con == 0);
			}
		
	for (int i=0; i<n; ++i)
		for (int j=1; j<m; ++j)
		{
			IloExpr con(env);
			for (int k=0; k<n; ++k)
				for (set<int>::iterator it=e_pre[i][j].begin(); it!=e_pre[i][j].end(); ++it)
					con += cplexx[*it][k];
			model.add(con == 1);
		}
		
	
	model.add(IloMinimize(env, cplexD));

	IloCplex cplex(model);	
	//warmstart
	{
		vector<vector<int> > greedy = gen_paths0_L(c,1);
		
		IloNumVarArray startVar(env);
		IloNumArray startVal(env);
		for (int k=0; k<n; ++k)		
			for (int e=0; e<E; ++e)
				{
					startVar.add(cplexx[e][k]);
					if (greedy[k][edges[e].fromj] == edges[e].fromi && greedy[k][edges[e].toj] == edges[e].toi)
					{
						startVal.add(1);
						//cout<<"setting "<<k<<","<<e<<"="<<edges[e].fromj<<","<<edges[e].fromi<<","<<edges[e].toj<<","<<edges[e].toi<<"\n";
					}
					else
						startVal.add(0);
				}
		cplex.addMIPStart(startVar, startVal);
		startVal.end();
		startVar.end();
	}

			
	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Threads, 1);
	cplex.setParam(IloCplex::TiLim, TIMELIMIT);
	//cplex.setParam(IloCplex::Param::Emphasis::MIP,CPX_MIPEMPHASIS_FEASIBILITY);
	
	double start = clock();
	bool result = cplex.solve();			
	double time = (clock() - start)/CLOCKS_PER_SEC;
	
	//cout<<n<<";"<<m<<flush;
	//cout<<";"<<time<<";"<<cplex.getObjValue()<<"\n";

	Solution sol;
	if (result)
	{
		sol.time = time;
		sol.obj = cplex.getObjValue();
		sol.status = cplex.getStatus();
		sol.gap = cplex.getMIPRelativeGap();
	}
	else
	{
		sol.time = TIMELIMIT;
		sol.obj = -1;
		sol.status = -1;
		sol.gap = -1;
	}
	
	//print solution
	//for (int k=0; k<n; ++k)
	//{
		//cout<<"Tour "<<k<<": ";
		//for (int i=0; i<n; ++i)
			//for (int j=0; j<m; ++j)	
				//if (cplex.getValue(cplexx[i][j][k]) > 0.5)
					//cout<<"("<<i<<","<<j<<") ";
		//cout<<"\n";
		//cout<<"Cost: ";
		//for (int i=0; i<n; ++i)
			//for (int j=0; j<m; ++j)	
				//if (cplex.getValue(cplexx[i][j][k]) > 0.5)
					//cout<<c[i][j]<<" + ";
		//cout<<"\n";
	//}
	
	env.end();
	
	return sol;
}



Solution MBA::solve_sparse_var()
{
	int E = edges.size();
	
	IloEnv env;
	IloModel model(env);
			
	vector<vector<IloNumVar> > cplexx(E);
	for (int e=0; e<E; ++e)
	{
		cplexx[e].resize(n);
		for (int k=0; k<n; ++k)
			cplexx[e][k] = IloNumVar(env,0,1,ILOBOOL);
	}
		
	IloNumVar cplexD(env,0,IloInfinity,ILOFLOAT);
	
	
	//preprocess reachability
	//for (int k=0; k<n; ++k)
	//{
		//set<int> reachable;
		//reachable.insert(k);
		
		//for (int j=1; j<m; ++j)
		//{
			//set<int> nreachable;
			//for (set<int>::iterator it=reachable.begin(); it!= reachable.end(); ++it)
				//for (set<int>::iterator it2=suc[*it][j-1].begin(); it2!=suc[*it][j-1].end(); ++it2)
					//nreachable.insert(*it2);
					
			//reachable = nreachable;
			//for (int i=0; i<n; ++i)
				//if (reachable.count(i) < 0.5)
					//for (set<int>::iterator it=e_pre[i][j].begin(); it!=e_pre[i][j].end(); ++it)
						//model.add(cplexx[*it][k] == 0);
		//}
	//}
	
	
	for (int k=0; k<n; ++k)
	{
		IloExpr con(env);
		for (set<int>::iterator it=e_suc[k][0].begin(); it!=e_suc[k][0].end(); ++it)
			con += cplexx[*it][k];
		model.add(con == 1);
	}
	
	for (int k=0; k<n; ++k)
	{
		IloExpr con(env);
		for (int e=0; e<E; ++e)
			con += edges[e].c * cplexx[e][k];
		con += c[k][0];
		model.add(con <= cplexD);
	}
	
	for (int i=0; i<n; ++i)
		for (int j=1; j<m-1; ++j)
			for (int k=0; k<n; ++k)
			{
				IloExpr con(env);
				for (set<int>::iterator it=e_suc[i][j].begin(); it!=e_suc[i][j].end(); ++it)
					con += cplexx[*it][k];
				for (set<int>::iterator it=e_pre[i][j].begin(); it!=e_pre[i][j].end(); ++it)
					con -= cplexx[*it][k];
				model.add(con == 0);
			}
		
	for (int i=0; i<n; ++i)
		for (int j=1; j<m; ++j)
		{
			IloExpr con(env);
			for (int k=0; k<n; ++k)
				for (set<int>::iterator it=e_pre[i][j].begin(); it!=e_pre[i][j].end(); ++it)
					con += cplexx[*it][k];
			model.add(con == 1);
		}		
		
	for (int i=0; i<n; ++i)
		for (int j=1; j<m-1; ++j)
		{
			IloExpr con(env);
			for (int k=0; k<n; ++k)
				for (set<int>::iterator it=e_suc[i][j].begin(); it!=e_suc[i][j].end(); ++it)
					con += cplexx[*it][k];
			model.add(con == 1);
		}
		
	//aggregated case constraints
	//if (edges[0].between > -0.5)
	//{
		//for (int i=0; i<n; ++i)
			//for (int j=0; j<m-1; ++j)
			//{
				//IloExpr con(env);
				//for (int e=0; e<E; ++e)
					//if (edges[e].fromj == j && edges[e].between == i)
						//for (int k=0; k<n; ++k)
							//con += cplexx[e][k];
				//model.add(con == 1);
			//}
				
	//}
	
	model.add(IloMinimize(env, cplexD));

	IloCplex cplex(model);	
	//warmstart
	
	if (edges[0].between < -0.5)
	{
		vector<vector<int> > greedy = gen_paths0_L(c,1);
		
		IloNumVarArray startVar(env);
		IloNumArray startVal(env);
		for (int k=0; k<n; ++k)		
			for (int e=0; e<E; ++e)
				{
					startVar.add(cplexx[e][k]);
					if (greedy[k][edges[e].fromj] == edges[e].fromi && greedy[k][edges[e].toj] == edges[e].toi)
					{
						startVal.add(1);
						//cout<<"setting "<<k<<","<<e<<"="<<edges[e].fromj<<","<<edges[e].fromi<<","<<edges[e].toj<<","<<edges[e].toi<<"\n";
					}
					else
						startVal.add(0);
				}
		cplex.addMIPStart(startVar, startVal);
		startVal.end();
		startVar.end();
	}
	else
	{
		m = 2*m-1;
		vector<vector<int> > greedy = gen_paths0_L(c,1);
		m = (m+1)/2;
		
		IloNumVarArray startVar(env);
		IloNumArray startVal(env);
		for (int k=0; k<n; ++k)		
			for (int e=0; e<E; ++e)
				{
					startVar.add(cplexx[e][k]);
					if (greedy[k][2*edges[e].fromj] == edges[e].fromi && greedy[k][2*edges[e].toj] == edges[e].toi && greedy[k][2*edges[e].fromj+1] == edges[e].between)
					{
						startVal.add(1);
						//cout<<"setting "<<k<<","<<e<<"="<<edges[e].fromj<<","<<edges[e].fromi<<","<<edges[e].between<<","<<edges[e].toj<<","<<edges[e].toi<<"\n";
					}
					else
						startVal.add(0);
				}
		cplex.addMIPStart(startVar, startVal);
		startVal.end();
		startVar.end();
	}

			
	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Threads, 1);
	cplex.setParam(IloCplex::TiLim, TIMELIMIT);
    //cplex.setParam(IloCplex::WorkMem, 128);
    //cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, 128);
    //cplex.setParam(IloCplex::Param::MIP::Strategy::File, 0);
	//cplex.setParam(IloCplex::Param::Emphasis::MIP,CPX_MIPEMPHASIS_FEASIBILITY);
	    
	Solution sol;
	
	double start = clock();	
	//long nodelim = 	cplex.getParam(IloCplex::NodeLim);
	//cplex.setParam(IloCplex::NodeLim, 0);
	bool result = cplex.solve();	
	//if ( result )
	//{
		//sol.root = cplex.getBestObjValue();
		//cplex.setParam(IloCplex::NodeLim, nodelim);
		//cplex.setParam(IloCplex::TiLim, TIMELIMIT - (clock() - start)/CLOCKS_PER_SEC);
		//result = cplex.solve();
	//}
	double time = (clock() - start)/CLOCKS_PER_SEC;
	
	//cout<<n<<";"<<m<<flush;
	//cout<<";"<<time<<";"<<cplex.getObjValue()<<"\n";

	if (result)
	{
		sol.time = time;
		sol.obj = cplex.getObjValue();
		sol.status = cplex.getStatus();
		sol.gap = cplex.getMIPRelativeGap();
		sol.vars = cplex.getNcols();
		sol.cons = cplex.getNrows();
		sol.lb = cplex.getBestObjValue();
	}
	else
	{
		sol.time = TIMELIMIT;
		sol.obj = -1;
		sol.status = -1;
		sol.gap = -1;
	}
	
	//print solution
	//for (int k=0; k<n; ++k)
	//{
		//cout<<"Tour "<<k<<": ";
		//for (int e=0; e<E; ++e)
			//if (cplex.getValue(cplexx[e][k]) > 0.5)
					//cout<<"[("<<edges[e].fromi<<","<<edges[e].fromj<<")-("<<edges[e].toi<<","<<edges[e].toj<<")] ";
		//cout<<"\n";
		//cout<<"Cost: ";
		//for (int e=0; e<E; ++e)
			//if (cplex.getValue(cplexx[e][k]) > 0.5)
					//cout<<edges[e].c<<" + ";
		//cout<<"\n";
	//}
	
	env.end();
	
	return sol;
}

Solution MBA::solve_triple()
{
	assert((m-1)%2 == 0);
	
	
	//update edge set and run sparse_var formulation
	vector<Edge> n_edges;
	
	int curj=0;
	for (int j=0; j<m-2; j+=2)
	{
		for (int i1=0; i1<n; ++i1)
			for (int i2=0; i2<n; ++i2)
				for (int i3=0; i3<n; ++i3)
					if (suc[i1][j].count(i2) > 0.5 && suc[i2][j+1].count(i3) > 0.5)
					{
						Edge e;
						e.fromi=i1;
						e.toi = i3;
						e.fromj = curj;
						e.toj = curj+1;
						e.c = c[i2][j+1] + c[i3][j+2];
						e.between = i2;
						
						n_edges.push_back(e);
					}
		++curj;
	}
	
	m = curj + 1;
	
	edges = n_edges;
	int E = edges.size();
	e_suc.clear();
	e_pre.clear();
	e_suc.resize(n);
	e_pre.resize(n);
	for (int i=0; i<n; ++i)
	{
		e_suc[i].resize(m);
		e_pre[i].resize(m);
	}
	
	
	for (int e=0; e<E; ++e)
	{
		e_suc[edges[e].fromi][edges[e].fromj].insert(e);
		e_pre[edges[e].toi][edges[e].toj].insert(e);
	}
	
	return solve_sparse_var();
}

void MBA::print_stats()
{
	//print nr of nodes and edges
	cout<<"Stats;"<<n*m<<";"<<edges.size()<<"\n";
}
