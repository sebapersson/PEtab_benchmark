# Run ODE solver benchmarks with 1 command line input (for accepted arguments see README)


# Neede to load Conda environment to have AMICI 
eval "$(conda shell.bash hook)"
conda activate PeTab


if [ $1 == "Nominal_values" ];then
    # Run the benchmark in Julia
    cd Master-Thesis
    bash Benchmarks/Run_cost_grad_hessian.sh Gradient_cost_small_models
    cd ..

    # Download AMICI and Benchmark collection
    git clone https://github.com/AMICI-dev/AMICI.git

    # Make a directory where to store results
    mkdir Intermediate/
    mkdir Intermediate/AMICI_cost_grad/

    cd AMICI
    git checkout develop
    # Download benchmark collection and associated benchmark scripts
    git clone --depth 1 https://github.com/benchmarking-initiative/Benchmark-Models-PEtab.git 
    export BENCHMARK_COLLECTION="$(pwd)/Benchmark-Models-PEtab/Benchmark-Models/"

    for i in {1..10}
    do
	AMICI_PARALLEL_COMPILE=2 tests/benchmark-models/test_benchmark_collection.sh
	pathSave="../Intermediate/AMICI_cost_grad/computation_times${i}.csv"
	mv tests/benchmark-models/computation_times.csv ${pathSave}
    done

    cd ..
    rm -rf AMICI
fi 


if [ $1 == "Fix_parameters" ];then
    cd Master-Thesis
    bash Benchmarks/Run_cost_grad_hessian.sh Fix_parameters
fi 


if [ $1 == "Test_chunks" ];then
    cd Master-Thesis
    bash Benchmarks/Run_cost_grad_hessian.sh Test_chunks
fi 


if [ $1 == "Test_adjoint_random_julia" ];then
    cd Master-Thesis
    bash Benchmarks/Run_cost_grad_hessian.sh Hessian_nominal_values
fi 


if [ $1 == "Test_adjoint_random_julia" ];then
    cd Master-Thesis
    bash Benchmarks/Run_cost_grad_hessian.sh Hessian_cost_small_models
fi 


if [ $1 == "Test_adjoint_random" ];then

    # Julia
    cd Julia_adjoint
    runJulia="/home/sebpe/julia-1.8.5-linux-x86_64/julia-1.8.5/bin/julia --project=. --threads=1"
    ${runJulia} Test_adjoint_random_p Boehm_JProteomeRes2014
    ${runJulia} Test_adjoint_random_p Bachmann_MSB2011
    ${runJulia} Test_adjoint_random_p Lucarelli_CellSystems2018
    ${runJulia} Test_adjoint_random_p Smith_BMCSystBiol2013
    cd ..

    # Python (AMICI)
    python amici_gradient_benchmark/gradient_benchmark_amici.py Boehm_JProteomeRes2014
    python amici_gradient_benchmark/gradient_benchmark_amici.py Bachmann_MSB2011
    python amici_gradient_benchmark/gradient_benchmark_amici.py Lucarelli_CellSystems2018
    python amici_gradient_benchmark/gradient_benchmark_amici.py Smith_BMCSystBiol2013 
fi 
