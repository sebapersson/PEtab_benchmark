#!/usr/bin/env python3

"""
Simulate a PEtab problem and compare results to reference values
"""
import numpy as np
import pandas as pd
import sys
from pathlib import Path

import petab
import amici
from amici.petab_objective import simulate_petab, RDATAS
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

    times = list()

    failures = list()

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
            scaled_parameters=True,
            scaled_gradients=True,
        )
        times.append(sum(r.cpu_time + r.cpu_timeB for r in res[RDATAS]))
        failures.append(not all(r.status == amici.AMICI_SUCCESS for r in res[RDATAS]))
        print(f'{ir + 1}/{len(parameters)} done')

    times = np.asarray(times)
    failures = np.asarray(failures)
    avg_time = times[np.logical_not(failures)].mean()

    print(f'finished {model}, average execution time was {avg_time} with {failures.sum()} failures.')

    df = pd.DataFrame(dict(time=times, failed=failures)).to_csv(
        parameter_dir / f'{model}_results.csv'
    )
