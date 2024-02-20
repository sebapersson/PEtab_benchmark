using DataFrames
using CSV
using Random
using LinearAlgebra
using Printf
using SciMLSensitivity
using Zygote
using OrdinaryDiffEq
using Sundials
using FiniteDifferences
using YAML
using BenchmarkTools
using PEtab

BLAS.set_num_threads(1)


function get_petab_problem(petab_model::PEtabModel,
                           gradient_method::Symbol,
                           sensealg,
                           ode_solver::ODESolver,
                           ode_solver_gradient::ODESolver,
                           sparse_jacobian::Bool)::PEtabODEProblem

    petab_problem = PEtabODEProblem(petab_model,
                                    gradient_method=gradient_method,
                                    ode_solver=ode_solver,
                                    ode_solver_gradient=ode_solver_gradient,
                                    sparse_jacobian=sparse_jacobian,
                                    sensealg=sensealg,
                                    verbose=false)
    return petab_problem
end


function get_runtime(compute_gradient::F, run_time, n_repeat, gradient, θ_est) where F

    compute_gradient(gradient, θ_est)
    b_grad =  @benchmark $compute_gradient($gradient, $θ_est) samples=(n_repeat+1) seconds=100000 evals=1
    println("b_grad1 = ", b_grad)
    b_grad =  @benchmark $compute_gradient($gradient, $θ_est) samples=(n_repeat+1) seconds=100000 evals=1
    println("b_grad2 = ", b_grad)
    run_time .= b_grad.times[2:end] .* 1e-9
end


function benchmarkCostGrad(petab_model::PEtabModel,
                           path_file_save,
                           gradientInfo,
                           ode_solver::ODESolver,
                           ode_solver_gradient::ODESolver,
                           ode_solver_name,
                           sparse_jacobian::Bool,
                           θ_est;
                           nParamFixed=nothing,
                           n_repeat=5,
                           iParameter=0,
                           path_save_gradient=nothing,
                           shouldSave::Bool=true)

    println("Running model ", petab_model.model_name)

    run_time = Vector{Float64}(undef, n_repeat)

    what_compute = "Gradient"
    gradient_method, sensealg, method_info = gradientInfo
    petab_problem = get_petab_problem(petab_model, gradient_method, sensealg, ode_solver, ode_solver_gradient, sparse_jacobian)
    gradient = zeros(length(θ_est))

    # Enzyme does not currently handle callbacks
    if petab_model.model_name == "model_Isensee_JCB2018" && (method_info == "InterpolatingAdjoint(autojacvec=EnzymeVJP())" || method_info == "QuadratureAdjoint(autojacvec=EnzymeVJP())")
        return
    end
    if petab_model.model_name == "Smith_BMCSystBiol2013" && (method_info == "InterpolatingAdjoint(autojacvec=EnzymeVJP())" || method_info == "QuadratureAdjoint(autojacvec=EnzymeVJP())")
        return
    end
    if petab_model.model_name == "model_Chen_MSB2009" && (method_info == "InterpolatingAdjoint(autojacvec=EnzymeVJP())" || method_info == "QuadratureAdjoint(autojacvec=EnzymeVJP())")
        return
    end

    println("Precompiling the code for method ", method_info)
    if gradient_method != "FiniteDifferences"
        compute_gradient! = petab_problem.compute_gradient!
    else
        compute_gradient! = (gradient, x) -> begin gradient .= FiniteDifferences.grad(central_fdm(4, 1), petab_problem.compute_cost, x)[1] end
    end
    local can_eval = true
    try
        compute_gradient!(gradient, θ_est)
    catch
        can_eval = false
    end
    if all(gradient .== 0.0) || can_eval == false
        run_time .= Inf
    else
        GC.gc(), GC.gc(), GC.gc()
        sleep(0.2)
        get_runtime(compute_gradient!, run_time, n_repeat, gradient, θ_est)
        println("BGrad = ", run_time)
    end

    writeParamFixed = isnothing(nParamFixed) ? 0 : nParamFixed
    data_save = DataFrame(Time = run_time,
                         What_calc=what_compute,
                         Method_info=method_info,
                         Model = petab_model.model_name,
                         abstol = ode_solver_gradient.abstol,
                         reltol = ode_solver_gradient.reltol,
                         N_param_fixed=writeParamFixed,
                         I_parameter=iParameter,
                         chunk_size = "0",
                         solver = ode_solver_name)

    if isfile(path_file_save) && shouldSave == true
        CSV.write(path_file_save, data_save, append = true)
    elseif shouldSave == true
        CSV.write(path_file_save, data_save)
    end

    if !isnothing(path_save_gradient)
        data_save = DataFrame(reshape(gradient, 1, length(gradient)), :auto)
        rename!(data_save, petab_problem.θ_names)
        data_save[!, :Method_info] .= method_info
        data_save[!, :Model] .= petab_model.model_name
        data_save[!, :abstol] .= ode_solver_gradient.abstol
        data_save[!, :reltol] .= ode_solver_gradient.reltol
        data_save[!, :I_parameter] .= iParameter
        data_save[!, :solver] .= ode_solver_name
        if isfile(path_save_gradient)
            CSV.write(path_save_gradient, data_save, append = true)
        else
            CSV.write(path_save_gradient, data_save)
        end
    end
end


model_test = "Smith_BMCSystBiol2013"
model_test = "Lucarelli_CellSystems2018"
model_test = "Bachmann_MSB2011"
model_test = "Boehm_JProteomeRes2014"

if ARGS[1] == "Test_adjoint_random_p"

    Random.seed!(123)

    dir_save = joinpath(@__DIR__, "..", "Master-Thesis", "Intermediate", "Benchmarks", "Adjoint_new")
    path_save = joinpath(dir_save, "Test_adjoint_random_new.csv")
    if !isdir(dir_save)
        mkpath(dir_save)
    end
    model_test = ARGS[2]
    path_save_gradient = joinpath(dir_save, "Test_adjoint_random_gradient_new" * model_test * ".csv")

    ode_solvers = [CVODE_BDF()]
    ode_solvers_name = ["CVODE_BDF"]
    sensealgsCheck = [[:Adjoint, InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true)), "InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true))"],
                      [:Adjoint, InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false)), "InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false))"],
                      [:Adjoint, InterpolatingAdjoint(autojacvec=EnzymeVJP()), "InterpolatingAdjoint(autojacvec=EnzymeVJP())"],
                      [:Adjoint, QuadratureAdjoint(abstol=1e-8, reltol=1e-8, autojacvec=ReverseDiffVJP(true)), "QuadratureAdjoint(autojacvec=ReverseDiffVJP(true))"],
                      [:Adjoint, QuadratureAdjoint(abstol=1e-8, reltol=1e-8, autojacvec=ReverseDiffVJP(false)), "QuadratureAdjoint(autojacvec=ReverseDiffVJP(false))"],
                      [:Adjoint, QuadratureAdjoint(abstol=1e-8, reltol=1e-8, autojacvec=EnzymeVJP()), "QuadratureAdjoint(autojacvec=EnzymeVJP())"],
                      [:Adjoint, GaussAdjoint(autojacvec=ReverseDiffVJP(true)), "GaussAdjoint(autojacvec=ReverseDiffVJP(true))"],
                      [:Adjoint, GaussAdjoint(autojacvec=ReverseDiffVJP(false)), "GaussAdjoint(autojacvec=ReverseDiffVJP(false))"],
                      [:Adjoint, GaussAdjoint(autojacvec=EnzymeVJP()), "GaussAdjoint(autojacvec=EnzymeVJP())"]]

    tolerances = [[1e-8, 1e-8], [1e-8, 1e-6]]
    dtmin = 1e-14
    if model_test != "Smith_BMCSystBiol2013"
        dir_model = joinpath(@__DIR__, "..", "Master-Thesis", "Intermediate", "PeTab_models", "model_" * model_test)
    else
        dir_model = joinpath(@__DIR__, "..", "Master-Thesis", "Intermediate", "PeTab_models", model_test)
    end
    path_yaml = joinpath(dir_model, model_test * ".yaml")
    petab_model = PEtabModel(path_yaml, build_julia_files=true)

    path_cube = joinpath(@__DIR__, "..", "Intermediate", "Parameters_test_gradient", model_test * ".csv")
    cube = Matrix(CSV.read(path_cube, DataFrame))
    p0 = cube[5, :]

    for i in 1:50
        @info "i = $i"
        if i == 1
            _petab_problem = PEtabODEProblem(petab_model)
            θ_est = _petab_problem.θ_nominalT
        else
            θ_est = cube[i, :]
        end
        for tols in tolerances
            reltol, abstol= tols
            for j in eachindex(ode_solvers)
                # Check Gradient
                ode_solver = ODESolver(ode_solvers[j], abstol=abstol, reltol=reltol, dtmin=dtmin)
                ode_solver_gradient = ODESolver(ode_solvers[j], abstol=abstol, reltol=reltol, dtmin=dtmin)
                for sensealgInfo in sensealgsCheck
                    benchmarkCostGrad(petab_model, path_save, sensealgInfo, ode_solver, ode_solver_gradient, ode_solvers_name[j], false, θ_est, path_save_gradient=path_save_gradient, n_repeat=5)
                end
            end
        end

        ode_solver = ODESolver(Rodas4P(), abstol=1e-8/100, reltol=1e-8/100, dtmin=dtmin)
        ode_solver_gradient = ODESolver(Rodas4P(), abstol=1e-8/100, reltol=1e-8/100, dtmin=dtmin)
        benchmarkCostGrad(petab_model, path_save, [:ForwardDiff, nothing, "ForwardDiff"], ode_solver, ode_solver_gradient, "Rodas4P", false, θ_est, path_save_gradient=path_save_gradient, n_repeat=5)
    end
end
