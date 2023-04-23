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

if [ $1 == "Adjoint" ];then
    python amici_gradient_benchmark/gradient_benchmark_amici.py Boehm_JProteomeRes2014
    python amici_gradient_benchmark/gradient_benchmark_amici.py Bachmann_MSB2011.csv
    python amici_gradient_benchmark/gradient_benchmark_amici.py Lucarelli_CellSystems2018
    python amici_gradient_benchmark/gradient_benchmark_amici.py Isensee_JCB2018 
    python amici_gradient_benchmark/gradient_benchmark_amici.py Smith_BMCSystBiol2013 
fi 
