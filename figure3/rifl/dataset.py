from rifl.filters.shared import Filter, FilterLevel
import rifl.utils as utils
import rifl.models as m
import rifl.config as config

from pandas.core.common import flatten
import pandas as pd
import numpy as np

from typing import List, Dict
from abc import ABC, abstractmethod
from collections import defaultdict
import pathlib


class _Dataset(ABC):
    
    @classmethod
    @abstractmethod
    def generate_id(cls, df):
        '''Generate a unique identifier per item.'''

    @abstractmethod
    def as_df(self, required_columns):
        '''Method for creating a pandas dataframe from the dataset.'''

    @abstractmethod
    def apply_filters(self, filters):
        '''Apply all filters'''

    @property
    def experiment_ids_included(self):
        return [e.id for e in self.experiments]

    def get_required_columns(self):
        filter_requirements = set()

        for level, fs in self.filters.items():
            filter_requirements.update(flatten(f.requires.database_fields for f in fs))

        requirements = self.base_requirements.union(filter_requirements)

        return [getattr(self.data_model, r) for r in requirements]

    def as_dict(self):
        return {
            'name': self.name,
            'experiments': self.experiment_ids_included
        }

class SilacDataset(_Dataset):
    base_requirements = {'id', 'experiment', 'uniprot', 'ratio', 'clean_sequence'}
    data_model = m.SilacData
    experiment_model = m.SilacExperiment

    def __init__(self,
        name: str,
        filters: Dict[FilterLevel, List[Filter]],
        replicate_paths: List[pathlib.Path],
        fasta_db_path: str,
        parser: str,
        base_url: str,
        combine_level: str = config.DEFAULT_SILAC_COMBINE_LEVEL,
        dta_folder_name: str = config.DEFAULT_DTA_FOLDER_NAME,
        meta = None
    ) -> None:
        self.name = name
        self.filters = filters
        self.replicate_paths = replicate_paths
        self.fasta_db_path = fasta_db_path
        self.meta = meta
        self.parser = parser
        self.base_url = base_url
        self.combine_level = combine_level
        self.dta_folder_name = dta_folder_name
        self.experiments = self._get_or_create_experiments()

    @classmethod
    def generate_id(cls, df):
        df['_id'] = df['uniprot']
        return df

    def as_df(self, required_columns):
        Data = self.data_model

        query = (Data
            .select(*required_columns)
            .where(Data.experiment_id.in_(self.as_dict()['experiments']))
        )

        df = pd.DataFrame.from_records(list(query.dicts()))
        df = df.set_index('id')

        return df

    def apply_filters(self):
        df = self.as_df(self.get_required_columns())
        df = self.generate_id(df)

        # replace zero with NaN to avoid counting in mean and stats calc
        df = df.replace(0, np.NaN)

        filtered_out = {}

        if 'individual' in self.filters:
            df, filtered_out['individual'] = self.apply_individual_filters(df)

        if 'experiment' in self.filters:
            df, filtered_out['experiment'] = self.apply_experiment_filters(df)

        if 'aggregate' in self.filters:
            df, filtered_out['aggregate'] = self.apply_aggregate_filters(df)

        df = df[~df.exclude]

        return df, filtered_out

    def apply_individual_filters(self, df):
        filter_category = 'individual'
        filtered_out = {}

        for f in self.filters.get('individual'):
            reqs = f.requires.dataframe_columns

            if isinstance(reqs, str):
                reqs = [reqs]

            f_colname = f'{filter_category}.{f.name}'

            df[f_colname] = df[reqs].apply(
                lambda x: f.fn(*(x[r] for r in reqs)),
                axis=1,
            )

            filtered_ids = set(df[~df[f_colname]].index.values.tolist())
            filtered_out[f.name] = filtered_ids

        df['exclude'] = ~df.filter(regex=f'^{filter_category}\.').all(axis=1)

        return df, filtered_out

    def apply_experiment_filters(self, df):
        filter_category = 'experiment'
        filtered_out = {}

        df = df[~df.exclude].copy()

        for f in self.filters.get('experiment'):
            f_colname = f'{filter_category}.{f.name}'

            grouped = df[~df.exclude].groupby(['_id', 'experiment'])
        
            df[f_colname] = grouped[f.requires.dataframe_columns].apply(
                lambda x: pd.Series(f.fn(x.tolist()), index=x.index)
            )

            filtered_out[f.name] = set(df[df[f_colname] == False].index.values.tolist())

            # make sure to exclude items that were filtered out for next filter
            df.loc[df[f_colname] == False, 'exclude'] = True

        return df, filtered_out


    def apply_aggregate_filters(self, df):
        filter_category = 'aggregate'
        filtered_out = defaultdict(set)

        df = df[~df.exclude].copy()

        for f in self.filters.get('aggregate'):
            f_colname = f'{filter_category}.{f.name}'

            # note that this implies our filters can only access the ratio column
            # instead could apply here and avoid two groupbys
            grouped = df[~df.exclude].reset_index().groupby(['_id', 'experiment']).agg({
                'ratio': 'median',
                'clean_sequence': list,
                'id': list
            })

            for _, g in grouped.groupby(level='_id'):
                data = g[f.requires.dataframe_columns].values.tolist()
                ids = g.id.values
                ids_to_filter_out = [x[1] for x in zip(f.fn(data), ids) if not x[0]]

                for y in ids_to_filter_out:
                    filtered_out[f.name].update(y)

            df.loc[filtered_out[f.name], f_colname] = False
            df.loc[filtered_out[f.name], 'exclude'] = True

        return df, filtered_out

    def _get_or_create_experiments(self):
        return [self.experiment_model.get_or_create(
            name=self.name,
            cimage_folder_path=path,
            fasta_db_path=self.fasta_db_path,
            base_url=self.base_url,
            parser=self.parser,
            dta_folder_name=self.dta_folder_name,
        ) for path in self.replicate_paths]

class _TmtDataset(_Dataset):
    base_requirements = {'id', 'experiment', 'uniprot', 'channel_intensities'}

    def as_df(self, required_columns):
        Data = self.data_model

        query = (Data
            .select(*required_columns)
            .where(Data.experiment_id.in_(self.as_dict()['experiments']))
        )

        df = pd.DataFrame.from_records(list(query.dicts()))
        df = df.set_index('id')

        for channel in self.channel_layout:
            channel_num, channel_name = next(iter(channel.items()))
            new_col_name = f'{channel_name}.intensity'
            # append number to column to distinguish duplicates
            new_col_name = f'{new_col_name}_{len(df.filter(regex=new_col_name).columns)}'
            df[new_col_name] = df.channel_intensities.str.get(channel_num - 1)
        df = df.drop(columns=['channel_intensities'])

        return df


class TmtDataset(_TmtDataset):
    data_model = m.TmtData
    experiment_model = m.TmtExperiment

    def __init__(self,
        name: str,
        filters: Dict[FilterLevel, List[Filter]],
        replicate_paths: List[str],
        channel_layout: dict,
        control_channels: str,
        fasta_db_path: str,
        meta = None
    ) -> None:
        self.name = name
        self.filters = filters
        self.replicate_paths = replicate_paths
        self.channel_layout = channel_layout
        self.control_channels = control_channels
        self.fasta_db_path = fasta_db_path
        self.meta = meta
        self.experiments = self._get_or_create_experiments()

    @classmethod
    def generate_id(cls, df):
        df['_id'] = df['uniprot']
        return df

    def _get_or_create_experiments(self):
        return [self.experiment_model.get_by_file_hash(path) or self.experiment_model.create_from_census(
            name=self.name,
            census_file_path=path,
            fasta_db_path=self.fasta_db_path
        ) for path in self.replicate_paths]

    def apply_filters(self):
        df = self.as_df(self.get_required_columns())
        df = self.generate_id(df)

        control_cols = df.filter(regex=f'{self.control_channels}.intensity')
        data_cols = df.filter(regex=f'.+\.intensity')

        percentages = data_cols.div(control_cols.mean(axis=1), axis=0).mul(100, axis=0)
        percentages = percentages.rename(columns={
            c: c.replace('intensity', 'percent_of_control')
            for c in df.filter(regex='\.intensity').columns
        })

        df = df.join(percentages)

        # note that this is done here since attrs get overwritten when two
        # dataframes (df + percentages) are joined
        df.attrs['control_intensities'] = list(control_cols)
        df.attrs['channel_intensities'] = list(data_cols)

        # replace zero with NaN to avoid counting in mean and stats calc
        # df = df.replace(0, np.NaN)

        filtered_out = {}

        if 'individual' in self.filters:
            df, filtered_out['individual'] = self.apply_individual_filters(df)

        if 'experiment' in self.filters:
            df, filtered_out['experiment'] = self.apply_experiment_filters(df)

        if 'aggregate' in self.filters:
            df, filtered_out['aggregate'] = self.apply_aggregate_filters(df)

        df = df[~df.exclude]

        return df, filtered_out

    def apply_individual_filters(self, df):
        filtered_out = {}

        for f in self.filters.get('individual'):
            reqs = f.requires.dataframe_columns
            filter_colname = f'individual.{f.name}'

            try:
                df[filter_colname] = f.fn(df[df.attrs.get(reqs, reqs)].values.T)
            except:
                df[filter_colname] = df[df.attrs.get(reqs, reqs)].apply(f.fn)

            filtered_out[f.name] = set(df[~df[filter_colname]].index.values.tolist())

        df['exclude'] = ~df.filter(regex=f'^individual\.').all(axis=1)

        return df, filtered_out

    def apply_experiment_filters(self, df):
        filtered_out = {}

        df = df[~df.exclude].copy()

        grouped = df.groupby(['_id', 'experiment']).agg({
            'uniprot': 'first',
            'experiment': list,
            **{c: 'mean' for c in df.filter(regex='percent_of_').columns}
        })

        for f in self.filters.get('experiment'):
            reqs = f.requires.dataframe_columns
            filter_colname = f'experiment.{f.name}'

            try:
                df[filter_colname] = f.fn(df[df.attrs.get(reqs, reqs)].values.T)
            except:
                df[filter_colname] = df[df.attrs.get(reqs, reqs)].apply(f.fn)

            filtered_out[f.name] = set(df[~df[filter_colname]].index.values.tolist())

        return df, filtered_out

    def apply_aggregate_filters(self, df):
        filtered_out = {}

        df  = df[~df.exclude].copy()

        for f in self.filters.get('aggregate'):
            reqs = f.requires.dataframe_columns
            filter_colname = f'aggregate.{f.name}'

            grouped = df[~df.exclude].reset_index()

        return df, filtered_out


# class TmtIsotopDataset(_Dataset):
#     def __init__(self,
#         name: str,
#         replicate_paths: List[str],
#         channel_layout: dict,
#         control_channels: str,
#         fasta_db_path: str,
#         reactive_residue: str,
#         meta = None
#     ) -> None:
#         self.name = name
#         self.replicate_paths = replicate_paths
#         self.channel_layout = channel_layout
#         self.control_channels = control_channels
#         self.fasta_db_path = fasta_db_path
#         self.reactive_residue = reactive_residue
#         self.meta = meta
#         self.experiments = self._get_or_create_experiments()

#     def _get_or_create_experiments(self):
#         return [m.TmtExperiment(
#             self.name,
#             path,
#             self.channel_layout,
#             reactive_residue=self.reactive_residue,
#             fasta_db_path=self.fasta_db_path
#         ) for path in self.replicate_paths]
