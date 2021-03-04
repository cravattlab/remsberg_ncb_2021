from .analysis import (
    analyze_tmt,
    analyze_silac
)


import yaml

def analyze(analysis_params_filepath: str, user: str):
    with open(analysis_params_filepath, 'r') as f:
        analysis_params = yaml.safe_load(f)

    analysis_params['user'] = user

    analysis_fn_map = {
        # 'tmt_isotop': analyze_tmt_isotop,
        'tmt': analyze_tmt,
        'silac': analyze_silac
    }

    analysis_fn = analysis_fn_map[analysis_params['type']]

    return analysis_fn(analysis_params), analysis_params
