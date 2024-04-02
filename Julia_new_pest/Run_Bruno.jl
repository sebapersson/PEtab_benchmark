using PEtab, OrdinaryDiffEq, DataFrames, CSV, Optim, PyCall, LinearAlgebra

BLAS.set_num_threads(1)

# For Fides
ENV["PYTHON"] = "/home/sebpe/anaconda3/envs/PeTab/bin/python"
import Pkg
Pkg.build("PyCall")

function writeFile(dirSave::String,
                   θ_opt::Vector{Float64},
                   parameterNames::Vector{String},
                   finalCost,
                   runTime,
                   retCode,
                   nIter,
                   startGuess,
                   alg::String,
                   solver::String,
                   absTol::String,
                   relTol::String;
                   tag::String="")

    # Save after each iteration (do not loose data)
    pathFile = joinpath(dirSave, "Estimation_statistics" * tag * ".csv")
    dataSave = [alg finalCost runTime retCode nIter startGuess solver absTol relTol]
    dataSave = DataFrame(dataSave, ["Alg", "Cost", "Run_time", "Ret_code", "N_iter", "Start_guess", "Solver", "absTol", "relTol"])
    shouldAppend = isfile(pathFile) ? true : false
    CSV.write(pathFile, dataSave, append=shouldAppend)

    # Save optimal parameter vector
    pathFile = joinpath(dirSave, "Minimizer" * tag * ".csv")
    _dataSave = Matrix{Any}(undef, (1, length(θ_opt)+5))
    _dataSave[:] .= vcat(θ_opt, startGuess, alg, solver, absTol, relTol)
    dataSaveθ = DataFrame(_dataSave, vcat(parameterNames, "Start_guess", "Alg", "Solver", "absTol", "relTol"))
    shouldAppend = isfile(pathFile) ? true : false
    CSV.write(pathFile, dataSaveθ, append=shouldAppend)

end

nmultistarts = 1000

path_yaml = joinpath(@__DIR__, "..", "Master-Thesis", "Intermediate", "PeTab_models", "model_Bruno_JExpBot2016", "Bruno_JExpBot2016.yaml")
petab_model = PEtabModel(path_yaml, build_julia_files=true)
prob1 = PEtabODEProblem(petab_model, ode_solver=ODESolver(Rodas5P(), abstol=1e-8, reltol=1e-8),
                        gradient_method=:ForwardDiff, hessian_method=:ForwardDiff)
prob2 = PEtabODEProblem(petab_model, ode_solver=ODESolver(Rodas5P(), abstol=1e-8, reltol=1e-8),
                        gradient_method=:ForwardEquations, hessian_method=:GaussNewton,
                        reuse_sensitivities=true)
prob3 = PEtabODEProblem(petab_model, ode_solver=ODESolver(Rodas5P(), abstol=1e-8, reltol=1e-8),
                        gradient_method=:ForwardDiff, hessian_method=:GaussNewton,
                        reuse_sensitivities=false)

θ_names = prob1.θ_names

pathCube = joinpath(petab_model.dir_julia, "Cube_benchmark.csv")
cube = Matrix(CSV.read(pathCube, DataFrame))

dir_result = joinpath(@__DIR__, "..", "Master-Thesis", "Intermediate", "Benchmarks", "Parameter_estimation", petab_model.model_name)
if !isdir(dir_result)
    mkpath(dir_result)
end

# Precompile Julia code
@info "Precompiling the code ..."
x = prob1.θ_nominalT
_ = prob1.compute_cost(x)
_ = prob2.compute_cost(x)
_ = prob3.compute_cost(x)
g1 = prob1.compute_gradient(x)
g2 = prob2.compute_gradient(x)
g3 = prob3.compute_gradient(x)
h1 = zeros(length(x), length(x))
h2, h3 = similar(h1), similar(h1)
prob1.compute_hessian!(h1, x)
prob2.compute_hessian!(h2, x)
prob3.compute_hessian!(h3, x)
println("Done")

for i in 1:nmultistarts

    @info "i = $i"
    p0 = cube[i, :]

    res = calibrate_model(prob1, p0, IPNewton(),
                          options=Optim.Options(iterations = 1000, show_trace = false, allow_f_increases=true,
                                                successive_f_tol = 3, f_tol=1e-8, g_tol=1e-6, x_tol=0.0))
    writeFile(dir_result, res.xmin, string.(θ_names), res.fmin, res.runtime, string(res.converged),
              res.n_iterations, i, "optimIPNewtonAutoHess", "Rodas5P", string(1e-8), string(1e-8))

    res = calibrate_model(prob3, p0, IPNewton(),
                          options=Optim.Options(iterations = 1000, show_trace = false, allow_f_increases=true,
                                                successive_f_tol = 3, f_tol=1e-8, g_tol=1e-6, x_tol=0.0))
    writeFile(dir_result, res.xmin, string.(θ_names), res.fmin, res.runtime, string(res.converged),
              res.n_iterations, i, "OptimIPNewtonGN", "Rodas5P", string(1e-8), string(1e-8))

    res = calibrate_model(prob1, p0, Fides(nothing; verbose=false),
                          options=py"{'maxiter' : 1000, 'fatol' : 0.0, 'frtol' : 1e-8, 'xtol' : 0.0, 'gatol' : 1e-6, 'grtol' : 0.0}"o)
    writeFile(dir_result, res.xmin, string.(θ_names), res.fmin, res.runtime, string(res.converged),
              res.n_iterations, i, "FidesAutoHess", "Rodas5P", string(1e-8), string(1e-8))

    res = calibrate_model(prob2, p0, Fides(nothing; verbose=false),
                          options=py"{'maxiter' : 1000, 'fatol' : 0.0, 'frtol' : 1e-8, 'xtol' : 0.0, 'gatol' : 1e-6, 'grtol' : 0.0}"o)
    writeFile(dir_result, res.xmin, string.(θ_names), res.fmin, res.runtime, string(res.converged),
              res.n_iterations, i, "FidesGN", "Rodas5P", string(1e-8), string(1e-8))

    res = calibrate_model(prob1, p0, Fides(:BFGS; verbose=false),
                          options=py"{'maxiter' : 1000, 'fatol' : 0.0, 'frtol' : 1e-8, 'xtol' : 0.0, 'gatol' : 1e-6, 'grtol' : 0.0}"o)
    writeFile(dir_result, res.xmin, string.(θ_names), res.fmin, res.runtime, string(res.converged),
              res.n_iterations, i, "FidesBFGS", "Rodas5P", string(1e-8), string(1e-8))
end
