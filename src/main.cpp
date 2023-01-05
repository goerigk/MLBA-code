#include "MBA.h"

#include <iostream>
#include <ctime>
#include <cstdlib>

using namespace std;

int main (int argc, char* argv[])
{
	//usage: ./MBA n m densefactor seed
	
	srand(10*atoi(argv[4]));
	
	int n = atoi(argv[1]);
	int m = atoi(argv[2]);

	MBA algo;
	algo.generate(n, m, atof(argv[3]));
	
	algo.print_stats();
	
	Solution sol;
	
	sol = algo.solve_form1b();
	cout<<"Form1b;"<<sol.vars<<";"<<sol.cons<<";"<<sol.time<<";"<<sol.lb<<";"<<sol.obj<<";"<<sol.gap<<";"<<sol.status<<"\n"<<flush;
	
	sol = algo.solve_sparse_var();
	cout<<"Sparsevar;"<<sol.vars<<";"<<sol.cons<<";"<<sol.time<<";"<<sol.lb<<";"<<sol.obj<<";"<<sol.gap<<";"<<sol.status<<"\n"<<flush;
	

	return 0;
}
