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
from amici.petab_objective import simulate_petab
def main():
    """Simulate the model specified on the command line"""

    model = sys.argv[1]

    benchmark_dir = Path(__file__).parent / 'pypesto_benchmark' / 'Benchmark-Models-PEtab'
    petab_yaml = 'Benchmark-Models-PEtab'

    # load PEtab files
    problem = petab.Problem.from_yaml(petab_yaml)
    petab.flatten_timepoint_specific_output_overrides(problem)
    problem.import_model()

    # load model
    amici_model = problem.create_model()
    amici_solver = problem.create_solver()

    amici_solver.setAbsoluteTolerance(1e-8)
    amici_solver.setRelativeTolerance(1e-8)
    amici_solver.setMaxSteps(int(1e4))
    if model in ('Brannmark_JBC2010', 'Isensee_JCB2018'):
        amici_model.setSteadyStateSensitivityMode(
            amici.SteadyStateSensitivityMode.integrationOnly
        )

    times = dict()

    failures = 0

    parameters = pd.read_csv(parameter_file)

    for parameter in parameters:
        amici_solver.setSensitivityMethod(amici.SensitivityMethod.adjoint)
        amici_solver.setSensitivityOrder(amici.SensitivityOrder.first)

        problem_parameters = dict(zip(parameters.index, parameter))

        res_repeats = [
            simulate_petab(
                petab_problem=problem, amici_model=amici_model,
                solver=amici_solver, problem_parameters=problem_parameters,
                scaled_gradients=True,
                scaled_parameters=True,
            )
            for _ in range(3)  # repeat to get more stable timings
        ]

        times.append(np.mean([
            sum(r.cpu_time + r.cpu_timeB + r.preeq_cpu_time + r.preeq_cpu_timeB for r in res[RDATAS]) / 1000
            # only forwards/backwards simulation
            for res in res_repeats
        ]))
        failures += not all(rdata.status == amici.AMICI_SUCCESS for rdata in rdatas)

    pd.Series(times).to_csv(
        f'./amici_gradient_benchmark/{model}_benchmark.csv'
    )


if __name__ == "__main__":
    main()
