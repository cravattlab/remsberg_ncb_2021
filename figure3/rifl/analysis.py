
from rifl.models import (
    database,
    create_tables,
)
from rifl.dataset import (
    TmtDataset,
    SilacDataset
)
import rifl.filters.isobaric as isobaric_filters
import rifl.filters.isotopic as isotopic_filters
import rifl.config as config
import rifl.models as m
import rifl.utils as utils

# installed modules
import pandas as pd

# base modules
from dataclasses import dataclass
from typing import Optional
from functools import partial
from collections import defaultdict
from collections.abc import Iterable
from abc import ABC, abstractmethod
import itertools
import dataclasses
import pathlib
import operator



@dataclass
class AnalysisParams:
    name: Optional[str] = None
    datasets: Optional[dict] = None
    filters: Optional[dict] = None
    fasta_db_path: Optional[str] = None
    datasets_base_path: Optional[str] = None
    output_folder: Optional[str] = None

    @classmethod
    def from_dict(cls, dct):
        return cls(**{f.name: dct[f.name] for f in dataclasses.fields(cls)})


@dataclass
class SilacAnalysisParams(AnalysisParams):
    parser: Optional[str] = config.DEFAULT_CIMAGE_PARSER
    combine: Optional[str] = config.DEFAULT_SILAC_COMBINE_LEVEL
    datasets_base_url: Optional[str] = ''
    dta_folder_name: Optional[str] = config.DEFAULT_DTA_FOLDER_NAME


class _Analysis(ABC):
    def __init__(self, params: AnalysisParams) -> None:
        # create all tables if they do not already exist
        create_tables()

        self.params = params
        self.name = params.name
        self.filters = self.init_filters(params.filters)
        self.fasta_db_path = pathlib.Path(params.fasta_db_path)
        self.output_path = pathlib.Path(params.output_folder)

        self.datasets = self.init_datasets(
            params.datasets,
            base_path=pathlib.Path(params.datasets_base_path)
        )

        self.fasta_db = m.FastaDatabaseIndex.get_or_create(self.fasta_db_path)
        self.results = None

    @property
    def filtered_out(self) -> set:
        filtered = set()

        for filtered_from_dataset in self._analysis.filters.values():
            for filters in filtered_from_dataset.values():
                for filtered_ids in filters.values():
                    filtered.update(filtered_ids)

        return filtered

    def _new_analysis(self):
        return m.AnalysisModel(
            name=self.name,
            datasets=[d.as_dict() for d in self.datasets]
        )

    def run(self):
        analysis = self._new_analysis()
        results, filtered_out = self.apply_filters()
        analysis.filters = filtered_out
        analysis.save()
        self._analysis = analysis
        database.close()

    def init_filters(self, filters_dict: dict):
        filter_objects = defaultdict(list)

        for level, fs in filters_dict.items():
            if fs is None:
                continue

            for f in fs:
                filter_objects[level].append(self.filter_module.Filter.from_dict(f, self.filter_module))

        return filter_objects

    def apply_filters(self):
        result = dict()

        for dataset in self.datasets:
            _, filtered_out = dataset.apply_filters()
            result[dataset.name] = filtered_out

        return (None, result)

    @property
    def experiments_included(self):
        return list(itertools.chain.from_iterable(d.experiments for d in self.datasets))

    @property
    def experiment_ids_included(self):
        return [e.id for e in self.experiments_included]

class TmtAnalysis(_Analysis):
    data_model = m.TmtData
    filter_module = isobaric_filters
    dataset_class = TmtDataset

    def init_datasets(self, datasets: dict, base_path=''):
        return [TmtDataset(
            name=dataset['name'],
            filters=self.filters,
            replicate_paths=[pathlib.Path(base_path).joinpath(r) for r in dataset['replicates']],
            channel_layout=dataset['channel_layout'],
            control_channels=dataset['control_channels'],
            fasta_db_path=self.fasta_db_path
        ) for dataset in datasets]

    def report(self, filename: str = ''):
        pass

    def filter_report(self, filename: str = ''):
        pass


class SilacAnalysis(_Analysis):
    data_model = m.SilacData
    filter_module = isotopic_filters
    dataset_class = SilacDataset

    def init_datasets(self, datasets: dict, base_path=''):
        return [self.dataset_class(
            name=dataset['name'],
            filters=self.filters,
            replicate_paths=[pathlib.Path(base_path).joinpath(r) for r in dataset['replicates']],
            fasta_db_path=self.fasta_db_path,
            parser=self.params.parser,
            base_url=self.params.datasets_base_url,
            combine_level=self.params.combine,
            dta_folder_name=self.params.dta_folder_name,
        ) for dataset in datasets]

    def filter_report(self):
        datasets_to_filter = self.experiment_ids_included

        m = self.data_model

        query = (m
            .select(
                m.id,
                m.experiment,
                m.uniprot,
                m.symbol,
                m.description,
                m.sequence,
                m.clean_sequence,
                m.ratio,
                m.num_ms2,
                m.rsquared,
                m.charge,
                m.meta
            )
            .where(m.experiment_id.in_(datasets_to_filter))
        )

        df = pd.DataFrame.from_records(list(query.dicts()))
        df = self.dataset_class.generate_id(df)

        for dataset in self.datasets:
            df.loc[df.experiment.isin(dataset.experiment_ids_included), 'condition'] = dataset.name

        df = df.set_index('id')

        df = df[[
            '_id',
            'experiment',
            'condition',
            'uniprot',
            'symbol',
            'sequence',
            'meta',
            'num_ms2',
            'rsquared',
            'ratio',
        ]]

        for filtered_out in self._analysis.filters.values():
            for cat, f in filtered_out.items():
                for filter_name, filtered_ids in f.items():
                    df.loc[filtered_ids, f'{cat}.{filter_name}'] = False

        # this should probably be just done in SQL
        # experiment_ids = df.experiment.unique().tolist()
        # q2 = Experiment.select(Experiment.id, Experiment.source_url).where(Experiment.id.in_(experiment_ids))
        # experiments = dict(list(q2.tuples()))
        # df.link = df.apply(lambda x: experiments[x.experiment] + x.link.split('"')[1], axis=1)

        report_output_name_template = 'filter_report_{}_{}_{{}}.xlsx'.format(
            self.name,
            self._analysis.id
        )

        report_output_path = utils.get_timestamped_report_path(
            report_output_name_template,
            self.output_path
        )

        # # df.set_index(['seq_id', 'condition', 'experiment']).sort_index(level=0).to_excel(report_output_path)
        df.set_index(['_id', 'experiment', 'condition']).sort_index(level=0).to_excel(report_output_path)
        return df


def analyze_tmt(analysis_params: dict):
    database.init(f'{analysis_params.get("user", config.DEFAULT_DB_NAME)}.db', pragmas=config.pragmas)
    params = AnalysisParams.from_dict(analysis_params)

    analysis = TmtAnalysis(params=params)
    analysis.run()
    return analysis

def analyze_silac(analysis_params: dict):
    database.init(f'{analysis_params.get("user", config.DEFAULT_DB_NAME)}.db', pragmas=config.pragmas)
    params = SilacAnalysisParams.from_dict(analysis_params)

    analysis = SilacAnalysis(params=params)
    analysis.run()
    return analysis
