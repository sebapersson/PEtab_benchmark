#!/usr/bin/env python3

"""
Simulate a PEtab problem and compare results to reference values
"""
import os

import numpy as np
import pandas as pd
import sys
from pathlib import Path

import petab
import amici
from amici.petab_objective import simulate_petab, RDATAS, SLLH
from amici.petab_import import import_petab_problem

"""Simulate the model specified on the command line"""

if __name__ == "__main__":
    model = sys.argv[1]

    benchmark_dir = Path(__file__).parents[1] / 'pypesto_benchmark' / \
                    'Benchmark-Models-PEtab' / 'Benchmark-Models'
    petab_yaml = benchmark_dir / model / f'{model}.yaml'

    # load PEtab files
    problem = petab.Problem.from_yaml(petab_yaml)
    petab.flatten_timepoint_specific_output_overrides(problem)

    # load model
    amici_model = import_petab_problem(problem)
    amici_solver = amici_model.getSolver()

    amici_solver.setAbsoluteTolerance(1e-8)
    amici_solver.setRelativeTolerance(1e-8)
    amici_solver.setMaxSteps(int(1e4))
    if model in ('Brannmark_JBC2010', 'Isensee_JCB2018'):
        amici_model.setSteadyStateSensitivityMode(
            amici.SteadyStateSensitivityMode.integrationOnly
        )

    sim_times = list()
    preeq_times = list()
    failures = list()
    grads = list()

    parameter_dir = Path(__file__).parents[1] / 'Intermediate' / \
                    'Parameters_test_gradient'
    parameters = pd.read_csv(parameter_dir / f'{model}.csv')

    for ir, parameter in parameters.iterrows():
        amici_solver.setSensitivityMethod(amici.SensitivityMethod.adjoint)
        amici_solver.setSensitivityOrder(amici.SensitivityOrder.first)

        problem_parameters = dict(parameter)

        res = simulate_petab(
            petab_problem=problem, amici_model=amici_model,
            solver=amici_solver, problem_parameters=problem_parameters,
            scaled_parameters=True
        )
        sim_times.append(sum(r.cpu_time + r.cpu_timeB for r in res[RDATAS]) / 1000)
        preeq_times.append(sum(r.preeq_cpu_time + r.preeq_cpu_timeB for r in res[RDATAS]) / 1000)
        failures.append(not all(r.status == amici.AMICI_SUCCESS for r in res[RDATAS]))
        grads.append([r.sllh for r in res[RDATAS]])
        
        msg = f'{ir + 1}/{len(parameters)} done'

        if failures[-1]:
            msg += ' (failed)'
        print(msg)

    sim_times = np.asarray(sim_times)
    preeq_times = np.asarray(preeq_times)
    failures = np.asarray(failures)
    avg_time = sim_times[np.logical_not(failures)].mean()
    df_grads = pd.DataFrame(np.concatenate(grads), columns=list(parameters.columns))

    print(f'finished {model}, average execution time was {avg_time} s with {failures.sum()} failures.')

    df = pd.DataFrame({
        ('simulation time',): sim_times,
        ('preequilibraiton time',): preeq_times,
        ('failed',): failures
    }).to_csv(
        parameter_dir / f'{model}_results.csv'
    )
    df_grads.to_csv(parameter_dir / f'{model}_results_grad.csv')
