# PEtab benchmark

This GitHub repository contains scripts for running PEtab parameter estimation benchmark using PyPesto (via the AMICI interface) and the Julia PEtab importer (the package does not have a cool name yet). 

## Setup 

Given i) a Julia 1.8.5 executable (or later version) at *path\_to\_julia*, ii) and the AMICI dependencies PyPesto and PEtab.jl can be installed by running ```bash setup.sh``` from the project root directory. 

1. **Note** : the Julia executable path must be set manually in *setup.sh*.
2. **Note** : the compilation time for Julia can be [significant](https://xkcd.com/303/). 

## Running the benchmark

Given that the setup was successful a benchmark can be run by executing ```bash Run_benchmark.sh modelName```, so for example ```bash Run_benchmark.sh Boehm_JProteomeRes2014``` will run the benchmark for the Boehm model.

1 **Note** : In the *Master_thesis/Benchmark/Run_benchmark.sh* file the Julia executable path must be set manually.

The available optimization algorithms that can be provided for the Julia are:

| Optimizer      | Description |
| ----------- | ----------- |
| IpoptAutoHess      | Ipopt using full hessian (via autodiff)|
| IpoptBlockAutoDiff   | Ipopt using block approximated hessian (via autodiff)        |
| IpoptLBFGS   | Ipopt using L-BFGS hessian approximation        |
| OptimIPNewtonAutoHess   | Optim.jl interior point Newton full hessian (via autodiff)   |
| OptimIPNewtonBlockAutoDiff   | Optim.jl interior point Newton via bloack approximated hessian (via autodiff)        |
| OptimIPNewtonGN   | Optim.jl interior point Newton with Gauss-Newton hessian approxmiation        |
| OptimLBFGS   | Optim.jl L-BFGS method|
| FidesAutoHess      | Fides using full hessian (via autodiff)|
| FidesBlockAutoHess   | Fides using block approximated hessian (via autodiff)        |
| FidesGN   | Fides using Gauss Newton hessian approximation.