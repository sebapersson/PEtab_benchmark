# PEtab Benchmark

This repository contains scripts for running the PEtab parameter estimation benchmark using pyPESTO (via the AMICI interface) and the precursor to PEtab.jl (before the formal package repository was created).

**Note**: The benchmark takes a substantial amount of time to run. Therefore, we provide all the intermediate result files on GitHub, which can be found under *Intermediate/* and *Master-thesis/Intermediate/Benchmarks*.

## Setup

To set up the benchmark, follow these steps:

1. Make sure you have a Julia 1.8.5 executable (or a later version) located at *path\_to\_julia*.
2. Install the necessary dependencies for AMICI and pyPESTO.
3. Run the command ```bash setup.sh``` from the project's root directory to install all the required packages.

**Note**: You will need to manually set the Julia executable path in the *setup.sh* file.

**Note**: The compilation time for Julia can be significant ([xkcd reference](https://xkcd.com/303/)).

To process the results, you will need R (version $\geq$ 4.0).

## ODE Solver Benchmark

### Running the Benchmark

After successfully setting up the environment, you can run the ODE solver benchmark by executing ```bash Run_ode_solver_benchmark.sh benchmark```, where ```benchmark``` refers to the specific benchmark you want to run. There are three different options available:

* ```Test_all``` : Test 29 benchmark models using all the ODE solvers listed in Table 2. For a subset of stiff solvers, different linear solvers are also tested.
* ```Test_random_parameters``` : Benchmark ODE solver performance for random parameter vectors on a subset of models (Figure 2c).
* ```Test_random_parameters_big_models``` : Benchmark ODE solver performance for random parameter vectors on a subset of larger models (Figure 3).

**Note**: You will need to manually set the Julia executable path in the *Master_thesis/Benchmark/Run_ode_solvers.sh* file.

### Processing Results

The results can be processed using the *Process_ode_solvers.R* script (make sure to run it from the script's location to set the path correctly). The results can be found in the *Results/ODE_solvers* folder.

## Cost, Gradient, and Hessian Benchmark

After successfully setting up the environment, you can run the cost, gradient, and Hessian benchmark by executing ```bash Run_gradient_benchmark.sh benchmark```, where ```benchmark``` refers to the specific benchmark you want to run. There are several options available:

* ```Nominal_values``` : Evaluate cost and gradient runtime for AMICI and PEtab.jl using reported parameter values on a subset of benchmark models (Figure 4).
* ```Hessian_nominal_values``` : Evaluate cost and Hessian runtime for PEtab.jl using reported parameter values on a subset of benchmark models (Figure 4).
* ```Fix_parameters``` : Evaluate gradient runtime by increasing the number of parameters to compute the gradient, while comparing default chunk-size versus not using chunking (Figure 3).
* ```Test_chunks``` : Evaluate gradient runtime for each possible chunk-size for the Bachmann, Iseense, and Lucarelli models (Figure 3).
* ```Test_adjoint_random_julia``` : Benchmark adjoint sensitivity analysis runtime for random parameter vectors using PEtab.jl on a subset of models (Figure 5).
* ```Test_adjoint_random_amici``` : Benchmark adjoint sensitivity analysis runtime for random parameter vectors using AMICI on a subset of models (Figure 5).

**Note**: You will need to manually set the Julia executable path in the *Master_thesis/Benchmark/Run_gradient_benchmark.sh* file.

### Processing Results

The results can be processed using the *Process_cost_gradient_hessian.R* script (make sure to run it from the script's location to set the path correctly). The results can be found in the *Results/Cost_gradient_hessian* folder.

## Parameter Estimation Benchmark

### Running the Benchmark

After successfully setting up the environment, you can run the benchmark by executing ```bash Run_benchmark.sh modelName```. For example, ```bash Run_benchmark.sh Boehm_JProteomeRes2014``` will run the benchmark for the Boehm model with 1000 multi-starts, using the optimizers specified in the configuration file.

**Note**: You need to manually set the Julia executable path in the *Master_thesis/Benchmark/Run_benchmark.sh* file.

The available optimization algorithms that can be provided for Julia are:

| Optimizer      | Description |
| ----------- | ----------- |
| IpoptAutoHess      | Ipopt using full Hessian (via autodiff)|
| IpoptBlockAutoDiff   | Ipopt using block-approximated Hessian (via autodiff)        |
| IpoptLBFGS   | Ipopt using L-BFGS Hessian approximation        |
| OptimIPNewtonAutoHess   | Optim.jl interior point Newton with full Hessian (via autodiff)   |
| OptimIPNewtonBlockAutoDiff   | Optim.jl interior point Newton with block-approximated Hessian (via autodiff)        |
| OptimIPNewtonGN   | Optim.jl interior point Newton with Gauss-Newton Hessian approximation        |
| OptimLBFGS   | Optim.jl L-BFGS method|
| FidesAutoHess      | Fides using full Hessian (via autodiff)|
| FidesBlockAutoHess   | Fides using block-approximated Hessian (via autodiff)        |
| FidesGN   | Fides using Gauss-Newton Hessian approximation

For pyPESTO, Fides with different Hessian approximations is used.

### Running the Bruno parameter estimation benchmark in Julia

To run the Bruno benchmark a different approach is used, because due to a bug in an earlier PEtab version the benchmark was performed using a later version of PEtab.jl (2.14.0). To perform the Bruno parameter estimation benchmark, enter the *Julia_new_pest* directory. Here, initialize the Julia environment by executing the following command in the Julia REPL:

```julia
] instantiate
```

Afterward, the benchmark can be initiated from within the *Julia_new_pest* directory using the command:

```bash
path_to_julia --threads=1 --project=. Run_bruno.jl
```

### Processing Results

The parameter estimation results can be processed using the *Process_results.R* script (make sure to run it from the script's location to set the path correctly). This script will generate individual convergence plots for each model, as well as summary visuals and tables. The results can be found in the *Results/Parameter_estimation* folder.

## Stochastic simulators benchmark

### Running the Benchmark

The stochastic simulators benchmark is located in the *Stochastic_benchmark* folder, which is download by the `setup.sh` script. Within the *Stochastic_benchmark* directory, there's another `setup.sh` script dedicated to installing the necessary packages for conducting the benchmark using SBMLImporter, RoadRunner, ReactionNetworkImporters, and PySB. In the `setup.sh` we recommend to use, as in the benchmarks, Julia version 1.10.

To execute the benchmark, navigate to the specific software directory and initiate one of the provided bash scripts. For instance, to run the benchmark with SBMLImporter using the `RSSACR` simulator, enter the following commands in a bash terminal:

```bash
cd SBMLImporter
bash run_benchmarks_sbmlimporter_ssa_RSSACR_1.sh
```

This process will execute the benchmark for several models; multistate, mulisite3, egfr_net, and fceri_gamma2. Note that the BCR model requires a separate script. The benchmarks are divided across multiple bash scripts to facilitate workload parallelization.

### Processing Results

The results can be processed using the *Plot.R* script (make sure to run it from the script's location to set the path correctly) in the *Plots* folder. This script will generate individual all the plots found in the paper.
