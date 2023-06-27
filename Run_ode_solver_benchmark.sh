# Run ODE solver benchmarks with 1 command line input 
# $1 - analysis to perform can be:
#   Test_all - test all ODE solvers (Fig. 2 and 3 (linear solvers) in the paper)
#   Test_random_parameters - test random parameters for subset of models (Fig. 2c)
#   Test_random_parameters_big_models - test random parameters for big models (Fig. 3d-i)

if [ $1 == "Test_all" ];then
    cd Master-Thesis 
    bash Benchmarks/Run_ode_solvers.sh Test_all
    exit 0
fi


if [ $1 == "Test_random_parameters" ];then
    cd Master-Thesis 
    bash Benchmarks/Run_ode_solvers.sh Test_random_parameters
    exit 0
fi


if [ $1 == "Test_random_parameters_big_models" ];then
    cd Master-Thesis 
    bash Benchmarks/Run_ode_solvers.sh Test_random_parameters_big_models
    exit 0
fi
