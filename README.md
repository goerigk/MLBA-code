# MLBA-code
Code to generate and solve multi-level bottleneck assignment problems.

It accompanies the paper on "Multi-level bottleneck assignment problems: complexity and sparsity-exploiting formulations" written by Trivikram Dokka and Marc Goerigk (see our preprint https://arxiv.org/abs/1910.12504).

Instances are generated through the "generate" function. The two main solution methods that are discussed in the paper are "solve_form1b" and "solve_sparse_var", but many more variantes were tested and are still provided in the code.

Code is written in C++ and requires Cplex to run. A typical compilation call might look as follows:

c++ src/main.cpp src/MBA.cpp -o MBA -I/opt/ibm/ILOG/CPLEX_Studio128/cplex/include/ -I/opt/ibm/ILOG/CPLEX_Studio128/concert/include/ -L/opt/ibm/ILOG/CPLEX_Studio128/cplex/lib/x86-64_linux/static_pic/ -lilocplex -lcplex -L/opt/ibm/ILOG/CPLEX_Studio128/concert/lib/x86-64_linux/static_pic/ -lconcert -lpthread -O3 -DIL_STD -Wno-ignored-attributes -ldl
