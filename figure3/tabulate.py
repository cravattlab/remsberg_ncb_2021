from utils import get_timestamped_report_path, get_newest_file
from db import UniprotDatabase, UniprotFlatfile
from constants import (
    MIN_INTENSITY_PER_CONTROL_CHANNEL,
    STANDARD_DEVIATION_FILTER_THRESHOLD,
    MIN_NUM_UNIQUE_SEQUENCES,
    MIN_NUM_DATASETS,
    HYDROXYLAMINE_REDUCTION_THRESHOLD,
    CACHE_DIR,
    DYNAMIC_FOLD_CHANGE_THRESHOLD,
    CHASE_FOLD_CHANGE_THRESHOLD
)
import pandas as pd
import numpy as np

def tabulate():
    time_course_path = get_newest_file('filtered_timecourse_*.xlsx')
    hydroxylamine_path = get_newest_file('filtered_hydroxylamine_*.xlsx')

    hydroxylamine = pd.read_excel(hydroxylamine_path, index_col='uniprot')
    hydroxylamine = hydroxylamine[hydroxylamine.mean_reduction >= HYDROXYLAMINE_REDUCTION_THRESHOLD]
    hydroxylamine = hydroxylamine[['gene_name', 'description', 'mean_reduction']]

    timecourse = pd.read_excel(time_course_path, index_col='uniprot')

    chase_col = 'DMSO t0/t1'
    palm_protect_col = 'Palm M t1/DMSO t1'
    abd_protect_col = '957 t1/DMSO t1'

    timecourse[chase_col] = timecourse['DMSO t0 (mean)']/timecourse['DMSO t1 (mean)']
    timecourse[palm_protect_col] = timecourse['Palm M t1 (mean)']/timecourse['DMSO t1 (mean)']
    timecourse[abd_protect_col] = timecourse['957 t1 (mean)']/timecourse['DMSO t1 (mean)']

    dynamic = timecourse[chase_col] >= DYNAMIC_FOLD_CHANGE_THRESHOLD
    palm_protected = timecourse[palm_protect_col] >= CHASE_FOLD_CHANGE_THRESHOLD
    abd_protected = timecourse[abd_protect_col] >= CHASE_FOLD_CHANGE_THRESHOLD

    timecourse['dynamic'] = 'no'
    timecourse.loc[dynamic, 'dynamic'] = 'yes'

    desired_cols = [
        'gene_name', 'description', 'dynamic', 'DMSO t0 (mean)', 'Palm M t0 (mean)',
        '957 t0 (mean)', 'DMSO t1 (mean)', 'Palm M t1 (mean)', '957 t1 (mean)'
    ]

    output_path = get_timestamped_report_path('supp_table_extra_sheets_{}.xlsx')

    with pd.ExcelWriter(output_path) as w:
        hydroxylamine.to_excel(w, sheet_name='Hydroxylamine Sensitive')
        timecourse[dynamic & ~palm_protected & ~abd_protected][desired_cols].drop(columns=['dynamic']).to_excel(w, sheet_name='Dynamic and not regulated')
        timecourse[palm_protected][desired_cols].to_excel(w, sheet_name='Palm M Regulated')
        timecourse[abd_protected][desired_cols].to_excel(w, sheet_name='ABD957 M Regulated')

def main():
    tabulate()

if __name__ == '__main__':
    main()
