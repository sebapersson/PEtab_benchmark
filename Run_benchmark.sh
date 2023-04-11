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


if [ $1 == "Beer_MolBioSystems2014" ];then
    echo "Running benchmark for Beer_MolBioSystems2014"
    juliaOptimizers="OptimIPNewtonAutoHess OptimIPNewtonGN FidesBFGS FidesGN FidesAutoHess"
    cd Master-Thesis 
    bash Benchmarks/Run_parameter_estimation.sh Beer_MolBioSystems2014 ${nMultiStarts} "${juliaOptimizers}"
    cd ../pypesto_benchmark
    python benchmark.py Beer_MolBioSystems2014 fides.hessian=FIM ${nMultiStarts}
    python benchmark.py Beer_MolBioSystems2014 fides.hessian=BFGS ${nMultiStarts}
    python benchmark.py Beer_MolBioSystems2014 fides ${nMultiStarts}
fi


if [ $1 == "Brannmark_JBC2010" ];then
    echo "Running benchmark for Brannmark_JBC2010"
    juliaOptimizers="OptimIPNewtonAutoHess OptimIPNewtonGN FidesBFGS FidesGN FidesAutoHess"
    cd Master-Thesis 
    bash Benchmarks/Run_parameter_estimation.sh Brannmark_JBC2010 ${nMultiStarts} "${juliaOptimizers}"
    cd ../pypesto_benchmark
    python benchmark.py Brannmark_JBC2010 fides.hessian=FIM ${nMultiStarts}
    python benchmark.py Brannmark_JBC2010 fides.hessian=BFGS ${nMultiStarts}
    python benchmark.py Brannmark_JBC2010 fides ${nMultiStarts}
fi


if [ $1 == "Bruno_JExpBot2016" ];then
    echo "Running benchmark for Bruno_JExpBot2016"
    juliaOptimizers="OptimIPNewtonAutoHess OptimIPNewtonGN FidesBFGS FidesGN FidesAutoHess"
    cd Master-Thesis 
    bash Benchmarks/Run_parameter_estimation.sh Bruno_JExpBot2016 ${nMultiStarts} "${juliaOptimizers}"
    cd ../pypesto_benchmark
    python benchmark.py Bruno_JExpBot2016 fides.hessian=FIM ${nMultiStarts}
    python benchmark.py Bruno_JExpBot2016 fides.hessian=BFGS ${nMultiStarts}
    python benchmark.py Bruno_JExpBot2016 fides ${nMultiStarts}
fi


if [ $1 == "Weber_BMC2015" ];then
    echo "Running benchmark for Weber_BMC2015"
    juliaOptimizers="OptimIPNewtonBlockAutoDiff OptimIPNewtonGN FidesBFGS FidesGN FidesBlockAutoHess"
    cd Master-Thesis 
    bash Benchmarks/Run_parameter_estimation.sh Weber_BMC2015 ${nMultiStarts} "${juliaOptimizers}"
    cd ../pypesto_benchmark
    python benchmark.py Weber_BMC2015 fides.hessian=FIM ${nMultiStarts}
    python benchmark.py Weber_BMC2015 fides.hessian=BFGS ${nMultiStarts}
    python benchmark.py Weber_BMC2015 fides ${nMultiStarts}
fi


if [ $1 == "Sneyd_PNAS2002" ];then
    echo "Running benchmark for Sneyd_PNAS2002"
    juliaOptimizers="OptimIPNewtonAutoHess OptimIPNewtonGN FidesBFGS FidesGN FidesAutoHess"
    cd Master-Thesis 
    bash Benchmarks/Run_parameter_estimation.sh Sneyd_PNAS2002 ${nMultiStarts} "${juliaOptimizers}"
    cd ../pypesto_benchmark
    python benchmark.py Sneyd_PNAS2002 fides.hessian=FIM ${nMultiStarts}
    python benchmark.py Sneyd_PNAS2002 fides.hessian=BFGS ${nMultiStarts}
    python benchmark.py Sneyd_PNAS2002 fides ${nMultiStarts}
fi


if [ $1 == "Lucarelli_CellSystems2018" ];then
    echo "Running benchmark for Lucarelli_CellSystems2018"
    juliaOptimizers="FidesBFGS FidesGN"
    cd Master-Thesis 
    bash Benchmarks/Run_parameter_estimation.sh Lucarelli_CellSystems2018 ${nMultiStarts} "${juliaOptimizers}"
    cd ../pypesto_benchmark
    python benchmark.py Lucarelli_CellSystems2018 fides.hessian=FIM ${nMultiStarts}
    python benchmark.py Lucarelli_CellSystems2018 fides.hessian=BFGS ${nMultiStarts}
    python benchmark.py Lucarelli_CellSystems2018 fides ${nMultiStarts}
fi


if [ $1 == "Schwen_PONE2014" ];then
    echo "Running benchmark for Schwen_PONE2014"
    juliaOptimizers="OptimIPNewtonBlockAutoDiff OptimIPNewtonGN FidesBFGS FidesGN FidesBlockAutoHess"
    cd Master-Thesis 
    bash Benchmarks/Run_parameter_estimation.sh Schwen_PONE2014 ${nMultiStarts} "${juliaOptimizers}"
    cd ../pypesto_benchmark
    python benchmark.py Schwen_PONE2014 fides.hessian=FIM ${nMultiStarts}
    python benchmark.py Schwen_PONE2014 fides.hessian=BFGS ${nMultiStarts}
    python benchmark.py Schwen_PONE2014 fides ${nMultiStarts}
fi


if [ $1 == "Elowitz_Nature2000" ];then
    echo "Running benchmark for Elowitz_Nature2000"
    juliaOptimizers="OptimIPNewtonAutoHess OptimIPNewtonGN FidesBFGS FidesGN FidesAutoHess"
    cd Master-Thesis 
    bash Benchmarks/Run_parameter_estimation.sh Elowitz_Nature2000 ${nMultiStarts} "${juliaOptimizers}"
    cd ../pypesto_benchmark
    python benchmark.py Elowitz_Nature2000 fides.hessian=FIM ${nMultiStarts}
    python benchmark.py Elowitz_Nature2000 fides.hessian=BFGS ${nMultiStarts}
    python benchmark.py Elowitz_Nature2000 fides ${nMultiStarts}
fi


if [ $1 == "Crauste_CellSystems2017" ];then
    echo "Running benchmark for Crauste_CellSystems2017"
    juliaOptimizers="OptimIPNewtonAutoHess OptimIPNewtonGN FidesBFGS FidesGN FidesAutoHess"
    cd Master-Thesis 
    bash Benchmarks/Run_parameter_estimation.sh Crauste_CellSystems2017 ${nMultiStarts} "${juliaOptimizers}"
    cd ../pypesto_benchmark
    python benchmark.py Crauste_CellSystems2017 fides.hessian=FIM ${nMultiStarts}
    python benchmark.py Crauste_CellSystems2017 fides.hessian=BFGS ${nMultiStarts}
    python benchmark.py Crauste_CellSystems2017 fides ${nMultiStarts}
fi


if [ $1 == "Isensee_JCB2018" ];then
    echo "Running benchmark for Isensee_JCB2018"
    juliaOptimizers="FidesBFGS FidesGN"
    cd Master-Thesis 
    bash Benchmarks/Run_parameter_estimation.sh Isensee_JCB2018 ${nMultiStarts} "${juliaOptimizers}"
    cd ../pypesto_benchmark
    python benchmark.py Isensee_JCB2018 fides.hessian=FIM ${nMultiStarts}
    python benchmark.py Isensee_JCB2018 fides.hessian=BFGS ${nMultiStarts}
    python benchmark.py Isensee_JCB2018 fides ${nMultiStarts}
fi
