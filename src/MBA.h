#ifndef _H_MBA_H_
#define _H_MBA_H_

#include <vector>
#include <set>

struct Solution
{
	double time;
	double obj;
	int status;
	double gap;
	
	int vars;
	int cons;
	double root;
	double lb;
};

struct Edge
{
	int fromi, fromj;
	int toi, toj;
	double c;
	
	int between;
};

class MBA
{
	public:
		void generate(int _n, int _m, double dense);
		
		Solution solve_form0(int method);
		Solution solve_form1b();
		Solution solve_form1c();
		Solution solve_form3();
		Solution solve_form4();
		
		Solution solve_sparse();
		Solution solve_sparse_var();
		Solution solve_triple();
		
		Solution solve_greedy();
		Solution solve_greedy_L(int L);
		
		Solution solve_form4_old();
	
		void print_problem();
		void print_stats();
		
		void set_para_method0(int _L);
		void set_para_method2(int _nr);
		void set_para_method3(double _thresh);
		
		void reset();
		
		void set_para_post(bool _post);
		
	private:
		bool paths_are_generated;
		void create_all_paths();
	
		std::vector<std::vector<int> > gen_paths0(std::vector<std::vector<double> > pcost);		
		std::vector<std::vector<int> > gen_paths1(std::vector<std::vector<double> > pcost);
		std::vector<std::vector<int> > gen_paths2(std::vector<std::vector<double> > pcost);
		std::vector<std::vector<int> > gen_paths3(std::vector<std::vector<double> > pcost);
		
		std::vector<std::vector<int> > gen_paths0_L(std::vector<std::vector<double> > pcost, int L);
		
		std::vector<std::vector<int> > post_improve(std::vector<std::vector<double> > pcost, std::vector<std::vector<int> > current);
		
		Solution solve_form0_master(std::vector<std::vector<int> > paths);
	
		std::vector<std::vector<double> > c;
		
		int n,m;
		
		std::vector<std::vector<std::set<int> > > suc;
		std::vector<std::vector<std::set<int> > > pre;
		
		std::vector<Edge> edges;
		std::vector<std::vector<std::set<int> > > e_suc;
		std::vector<std::vector<std::set<int> > > e_pre;
		
		std::vector<std::vector<int> > allpaths;
		
		int para_method0_L;
		int para_method2_nr;
		double para_method3_thresh;
		bool para_post;
};

#endif
