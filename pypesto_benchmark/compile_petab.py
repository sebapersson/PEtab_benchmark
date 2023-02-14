import pypesto.petab
import petab
import os
import sys

import pandas as pd
import numpy as np

folder_base = os.path.join(os.path.dirname(__file__),
                           'Benchmark-Models-PEtab',
                           'Benchmark-Models')


def preprocess_problem(problem, model, extend_bounds):
    if model in ['Brannmark_JBC2010', 'Fiedler_BMC2016']:
        petab.flatten_timepoint_specific_output_overrides(problem)

    if np.isfinite(extend_bounds):
        problem.parameter_df[petab.LOWER_BOUND] /= extend_bounds
        problem.parameter_df[petab.UPPER_BOUND] *= extend_bounds
    else:
        problem.parameter_df.loc[
            problem.parameter_df[petab.PARAMETER_SCALE] == petab.LIN,
            petab.LOWER_BOUND
        ] = - np.inf
        problem.parameter_df.loc[
            problem.parameter_df[petab.PARAMETER_SCALE] != petab.LIN,
            petab.LOWER_BOUND
        ] = 0
        problem.parameter_df[petab.UPPER_BOUND] = np.inf


def load_problem(model, force_compile=False, extend_bounds=1.0):
    yaml_config = os.path.join(folder_base, model, model + '.yaml')
    petab_problem = petab.Problem.from_yaml(yaml_config)
    preprocess_problem(petab_problem, model, extend_bounds)
    importer = pypesto.petab.PetabImporter(petab_problem, validate_petab=False)
    problem = importer.create_problem(force_compile=force_compile)

    pnames = problem.get_reduced_vector(problem.x_names, problem.x_free_indices)

    pstart = pd.read_csv(os.path.join(os.path.dirname(__file__), 'starting_points', f'{model}.csv'))

    problem.x_guesses_full = problem.get_full_vector(pstart[pnames].values, problem.x_fixed_vals)

    return petab_problem, problem


if __name__ == '__main__':
    MODEL_NAME = sys.argv[1]

    load_problem(MODEL_NAME, force_compile=True)
