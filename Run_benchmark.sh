# Ensure we run the script from the correct directory (in order for paths to set correctly)
currentDir=${PWD##*/}
if [ ! $currentDir == "PEtab_benchmark" ]; then
    >&2 echo "Error : Script must be run from directory PEtab_benchmark to set it paths correctly"
    exit 1
fi

# To enable Python packages in Julia
eval "$(conda shell.bash hook)"
conda activate PeTab

nMultiStarts="1000"

if [ $1 == "Boehm_JProteomeRes2014" ];then
    echo "Running benchmark for Boehm_JProteomeRes2014"
    juliaOptimizers="OptimIPNewtonAutoHess OptimIPNewtonGN FidesBFGS FidesGN FidesAutoHess"
    cd Master-Thesis 
    bash Benchmarks/Run_parameter_estimation.sh Boehm_JProteomeRes2014 ${nMultiStarts} "${juliaOptimizers}"
    cd ../pypesto_benchmark
    python benchmark.py Boehm_JProteomeRes2014 fides.hessian=FIM ${nMultiStarts}
    python benchmark.py Boehm_JProteomeRes2014 fides.hessian=BFGS ${nMultiStarts}
    python benchmark.py Boehm_JProteomeRes2014 fides ${nMultiStarts}
fi


if [ $1 == "Fiedler_BMC2016" ];then
    echo "Running benchmark for Fiedler_BMC2016"
    juliaOptimizers="OptimIPNewtonAutoHess OptimIPNewtonGN FidesBFGS FidesGN FidesAutoHess"
    cd Master-Thesis 
    bash Benchmarks/Run_parameter_estimation.sh Fiedler_BMC2016 ${nMultiStarts} "${juliaOptimizers}"
    cd ../pypesto_benchmark
    python benchmark.py Fiedler_BMC2016 fides.hessian=FIM ${nMultiStarts}
    python benchmark.py Fiedler_BMC2016 fides.hessian=BFGS ${nMultiStarts}
    python benchmark.py Fiedler_BMC2016 fides ${nMultiStarts}
fi


if [ $1 == "Fujita_SciSignal2010" ];then
    echo "Running benchmark for Fujita_SciSignal2010"
    juliaOptimizers="OptimIPNewtonAutoHess OptimIPNewtonGN FidesBFGS FidesGN FidesAutoHess"
    cd Master-Thesis 
    bash Benchmarks/Run_parameter_estimation.sh Fujita_SciSignal2010 ${nMultiStarts} "${juliaOptimizers}"
    cd ../pypesto_benchmark
    python benchmark.py Fujita_SciSignal2010 fides.hessian=FIM ${nMultiStarts}
    python benchmark.py Fujita_SciSignal2010 fides.hessian=BFGS ${nMultiStarts}
    python benchmark.py Fujita_SciSignal2010 fides ${nMultiStarts}
fi


if [ $1 == "Bachmann_MSB2011" ];then
    echo "Running benchmark for Bachmann_MSB2011"
    juliaOptimizers="OptimIPNewtonGN FidesBFGS FidesGN"
    cd Master-Thesis 
    bash Benchmarks/Run_parameter_estimation.sh Bachmann_MSB2011 ${nMultiStarts} "${juliaOptimizers}"
    cd ../pypesto_benchmark
    python benchmark.py Bachmann_MSB2011 fides.hessian=FIM ${nMultiStarts}
    python benchmark.py Bachmann_MSB2011 fides.hessian=BFGS ${nMultiStarts}
    python benchmark.py Bachmann_MSB2011 fides ${nMultiStarts}
fi
