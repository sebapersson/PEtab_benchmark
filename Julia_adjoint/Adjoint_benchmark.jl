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


function getPEtabProblem(petabModel::PEtabModel, 
                         gradientMethod::Symbol, 
                         sensealg, 
                         odeSolverOptions::ODESolverOptions, 
                         odeSolverGradientOptions::ODESolverOptions,
                         sparseJacobian::Bool)::PEtabODEProblem

    petabProblem = createPEtabODEProblem(petabModel, 
                                         gradientMethod=gradientMethod,
                                         odeSolverOptions=odeSolverOptions, 
                                         odeSolverGradientOptions=odeSolverGradientOptions, 
                                         sparseJacobian=sparseJacobian, 
                                         sensealg=sensealg, 
                                         verbose=false)                         

    return petabProblem
end


function getRunTime(computeGradient::F, runTime, nRepeat, gradient, θ_est) where F

    computeGradient(gradient, θ_est)
    bGrad =  @benchmark $computeGradient($gradient, $θ_est) samples=(nRepeat+1) seconds=100000 evals=1
    println("bGrad1 = ", bGrad)
    bGrad =  @benchmark $computeGradient($gradient, $θ_est) samples=(nRepeat+1) seconds=100000 evals=1
    println("bGrad2 = ", bGrad)
    runTime .= bGrad.times[2:end] .* 1e-9

end


function benchmarkCostGrad(petabModel::PEtabModel, 
                           pathFileSave,
                           gradientInfo,
                           odeSolverOptions::ODESolverOptions, 
                           odeSolverGradientOptions::ODESolverOptions,
                           odeSolverName,
                           sparseJacobian::Bool, 
                           θ_est;
                           nParamFixed=nothing, 
                           nRepeat=5, 
                           iParameter=0,
                           pathSaveGradient=nothing,
                           shouldSave::Bool=true)

    println("Running model ", petabModel.modelName)

    
    runTime = Vector{Float64}(undef, nRepeat)

    whatCompute = "Gradient"
    gradientMethod, sensealg, methodInfo = gradientInfo
    petabProblem = getPEtabProblem(petabModel, gradientMethod, sensealg, odeSolverOptions, odeSolverGradientOptions, sparseJacobian)
    gradient = zeros(length(θ_est))
        
    # Enzyme does not currently handle callbacks 
    if petabModel.modelName == "model_Isensee_JCB2018" && (methodInfo == "InterpolatingAdjoint(autojacvec=EnzymeVJP())" || methodInfo == "QuadratureAdjoint(autojacvec=EnzymeVJP())")
        return 
    end
    if petabModel.modelName == "Smith_BMCSystBiol2013" && (methodInfo == "InterpolatingAdjoint(autojacvec=EnzymeVJP())" || methodInfo == "QuadratureAdjoint(autojacvec=EnzymeVJP())")
        return 
    end
    if petabModel.modelName == "model_Chen_MSB2009" && (methodInfo == "InterpolatingAdjoint(autojacvec=EnzymeVJP())" || methodInfo == "QuadratureAdjoint(autojacvec=EnzymeVJP())")
        return 
    end

    println("Precompiling the code for method ", methodInfo)
    if gradientMethod != "FiniteDifferences"
        computeGradient! = petabProblem.computeGradient!
    else
        computeGradient! = (gradient, x) -> begin gradient .= FiniteDifferences.grad(central_fdm(4, 1), petabProblem.computeCost, x)[1] end
    end
    local canEval = true
    try 
        computeGradient!(gradient, θ_est)
    catch 
        canEval = false
    end
    if all(gradient .== 0.0) || canEval == false
        runTime .= Inf
    else
        GC.gc(), GC.gc(), GC.gc()
        sleep(0.2)
        getRunTime(computeGradient!, runTime, nRepeat, gradient, θ_est) 
        println("BGrad = ", runTime)
    end

    writeParamFixed = isnothing(nParamFixed) ? 0 : nParamFixed
    dataSave = DataFrame(Time = runTime, 
                         What_calc=whatCompute,
                         Method_info=methodInfo,
                         Model = petabModel.modelName, 
                         absTol = absTol, 
                         relTol = relTol,
                         N_param_fixed=writeParamFixed,
                         I_parameter=iParameter,
                         chunk_size = "0",
                         solver = odeSolverName)

    if isfile(pathFileSave) && shouldSave == true
        CSV.write(pathFileSave, dataSave, append = true)
    elseif shouldSave == true
        CSV.write(pathFileSave, dataSave)
    end

    if !isnothing(pathSaveGradient)
        dataSave = DataFrame(reshape(gradient, 1, length(gradient)), :auto)
        rename!(dataSave, petabProblem.θ_estNames)
        dataSave[!, :Method_info] .= methodInfo
        dataSave[!, :Model] .= petabModel.modelName
        dataSave[!, :absTol] .= absTol
        dataSave[!, :relTol] .= relTol
        dataSave[!, :I_parameter] .= iParameter
        dataSave[!, :solver] .= odeSolverName
        if isfile(pathSaveGradient)
            CSV.write(pathSaveGradient, dataSave, append = true)
        else
            CSV.write(pathSaveGradient, dataSave)
        end
    end
end


modelTest = "Smith_BMCSystBiol2013"

if ARGS[1] == "Test_adjoint_random_p"

    Random.seed!(123)

    dirSave = joinpath(@__DIR__, "..", "Master-Thesis", "Intermediate", "Benchmarks", "Cost_grad_hess_fix")
    pathSave = joinpath(dirSave, "Test_adjoint_random_fix.csv")
    if !isdir(dirSave)
        mkpath(dirSave)
    end
    modelTest = ARGS[2]
    pathSaveGradient = joinpath(dirSave, "Test_adjoint_random_gradient_fix_" * modelTest * ".csv")
                            
    odeSolvers = [CVODE_BDF()]
    odeSolversName = ["CVODE_BDF"]                 
    sensealgsCheck = [[:Adjoint, InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true)), "InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true))"],
                      [:Adjoint, InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false)), "InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false))"],
                      [:Adjoint, InterpolatingAdjoint(autojacvec=EnzymeVJP()), "InterpolatingAdjoint(autojacvec=EnzymeVJP())"], 
                      [:Adjoint, QuadratureAdjoint(abstol=1e-8, reltol=1e-8, autojacvec=ReverseDiffVJP(true)), "QuadratureAdjoint(autojacvec=ReverseDiffVJP(true))"],
                      [:Adjoint, QuadratureAdjoint(abstol=1e-8, reltol=1e-8, autojacvec=ReverseDiffVJP(false)), "QuadratureAdjoint(autojacvec=ReverseDiffVJP(false))"], 
                      [:Adjoint, QuadratureAdjoint(abstol=1e-8, reltol=1e-8, autojacvec=EnzymeVJP()), "QuadratureAdjoint(autojacvec=EnzymeVJP())"]]

    absTol, relTol = 1e-8, 1e-8
    dtmin = 1e-14
    if modelTest != "Smith_BMCSystBiol2013"
        dirModel = joinpath(@__DIR__, "..", "Master-Thesis", "Intermediate", "PeTab_models", "model_" * modelTest)
    else
        dirModel = joinpath(@__DIR__, "..", "Master-Thesis", "Intermediate", "PeTab_models", modelTest)
    end
    pathYML = joinpath(dirModel, modelTest * ".yaml")
    petabModel = readPEtabModel(pathYML, forceBuildJuliaFiles=true)

    pathCube = joinpath(@__DIR__, "..", "Intermediate", "Parameters_test_gradient", modelTest * ".csv")
    cube = Matrix(CSV.read(pathCube, DataFrame))

    for i in 1:50
        if i == 1
            _petabProblem = createPEtabODEProblem(petabModel)
            θ_est = _petabProblem.θ_nominalT
        else
            θ_est = cube[i, :]
        end
        j = 1
        for j in eachindex(odeSolvers)
            # Check Gradient 
            odeSolverOptions = ODESolverOptions(odeSolvers[j], abstol=absTol, reltol=relTol, dtmin=dtmin)
            odeSolverGradientOptions = ODESolverOptions(odeSolvers[j], abstol=absTol, reltol=relTol, dtmin=dtmin)
            for sensealgInfo in sensealgsCheck
                benchmarkCostGrad(petabModel, pathSave, sensealgInfo, odeSolverOptions, odeSolverGradientOptions, odeSolversName[j], false, θ_est, pathSaveGradient=pathSaveGradient, nRepeat=5)
            end            
        end
         
        odeSolverOptions = ODESolverOptions(Rodas4P(), abstol=absTol/100, reltol=relTol/100, dtmin=dtmin)
        odeSolverGradientOptions = ODESolverOptions(Rodas4P(), abstol=absTol/100, reltol=relTol/100, dtmin=dtmin)
        benchmarkCostGrad(petabModel, pathSave, [:ForwardDiff, nothing, "ForwardDiff"], odeSolverOptions, odeSolverGradientOptions, "Rodas4P", false, θ_est, pathSaveGradient=pathSaveGradient, nRepeat=5)
        benchmarkCostGrad(petabModel, pathSave, [:ForwardDiff, nothing, "FiniteDifferences"], odeSolverOptions, odeSolverGradientOptions, "Rodas4P", false, θ_est, pathSaveGradient=pathSaveGradient, nRepeat=5)
    end
end
