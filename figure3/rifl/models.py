from rifl.parse import (
    parse_census,
    parse_dtaselect,
    parse_cimage_flatfile
)
import rifl.utils as utils
import rifl.config as config
import rifl.exceptions as exceptions

from Bio import SwissProt
from playhouse.sqlite_ext import (
    SqliteExtDatabase,
    JSONField,
)
from peewee import *

import datetime
import hashlib
import pathlib
import json


database = SqliteExtDatabase(None)

def create_tables():
    with database:
        database.create_tables([
            AnalysisModel,
            Blacklist,
            Whitelist,
            FastaDatabase,
            FastaDatabaseIndex,
            SilacExperiment,
            IsotopExperiment,
            TmtExperiment,
            TmtIsotopExperiment,
            SilacData,
            IsotopData,
            TmtData,
            TmtIsotopData,
        ])


class BaseModel(Model):
    class Meta:
        database = database
        legacy_table_names = False

class AnalysisModel(BaseModel):
    """Data on analyses of multiple datasets."""
    name = TextField()
    # analysis_hash = TextField(default=utils.get_git_revision_hash)
    analysis_hash = TextField(default='')
    date = DateField(default=datetime.datetime.now)
    datasets = JSONField()
    filters = JSONField(json_dumps=lambda x: json.dumps(x, cls=utils.JSONSetEncoder), null=True)

class BaseExperiment(BaseModel):
    """Holds experimental metadata."""
    name = TextField()
    date = DateField(default=datetime.datetime.now)
    meta = JSONField(null=True)

class SilacExperiment(BaseExperiment):
    """Holds experimental metadata."""
    file_hash = TextField(unique=True)
    dta_file_hash = TextField(unique=True)
    cimage_folder_name = TextField()
    database = DeferredForeignKey('FastaDatabaseIndex')

    @classmethod
    def get_dta_hash(cls, cimage_folder_path: pathlib.Path, dta_folder_name: str = config.DEFAULT_DTA_FOLDER_NAME) -> str:
        dta_file_path = cls.get_dta_file_in_folder(cimage_folder_path, dta_folder_name)
        return utils.get_file_hash(dta_file_path)

    @classmethod
    def get_file_hash(cls, cimage_folder_path: pathlib.Path, parser: str, dta_folder_name: str = config.DEFAULT_DTA_FOLDER_NAME) -> str:
        if parser == 'flatfile':
            file_path = cls.get_flatfile_in_folder(cimage_folder_path, dta_folder_name)
        elif parser in ('protein', 'peptide'):
            file_path = cimage_folder_path.joinpath(f'combined_{dta_folder_name}.txt')

        return utils.get_file_hash(file_path)

    @classmethod
    def get_by_path_hash(cls, cimage_folder_path: pathlib.Path, parser: str, dta_folder_name: str = config.DEFAULT_DTA_FOLDER_NAME):
        dta_file_hash = cls.get_dta_hash(cimage_folder_path, dta_folder_name)
        file_hash = cls.get_file_hash(cimage_folder_path, parser, dta_folder_name)

        return cls.get_or_none(
            (cls.file_hash == file_hash) |
            (cls.dta_file_hash == dta_file_hash)
        )

    @classmethod
    def get_flatfile_in_folder(cls, cimage_folder_path: pathlib.Path, dta_folder_name: str = config.DEFAULT_DTA_FOLDER_NAME):
        search_glob = 'output_*.to_excel.txt'
        flatfile_path = next(cimage_folder_path.glob(search_glob), None)
    
        if flatfile_path is None:
            raise exceptions.CimageFlatFileNotFound(
                f'Cimage flatfile not found in folder {cimage_folder_path}, using pattern {search_glob}'
            )

        return flatfile_path

    @classmethod
    def get_dta_file_in_folder(cls, cimage_folder_path: pathlib.Path, dta_folder_name: str = config.DEFAULT_DTA_FOLDER_NAME):
        search_glob = 'DTASelect-filter*.txt'
        dta_file_path = next(cimage_folder_path.glob(search_glob), None)

        if dta_file_path is None:
            raise exceptions.DtaFileNotFound(
                f'DTASelect file not found in folder {cimage_folder_path}, using pattern {search_glob}'
            )

        return dta_file_path

    @classmethod
    def get_or_create(cls,
        name: str,
        cimage_folder_path: pathlib.Path,
        fasta_db_path: pathlib.Path,
        base_url: str,
        parser: str = config.DEFAULT_CIMAGE_PARSER,
        dta_folder_name: str = config.DEFAULT_DTA_FOLDER_NAME,
    ):
        create_methods = {
            'flatfile': cls.create_from_flatfile
        }

        existing = None

        try:
            existing = cls.get_by_path_hash(cimage_folder_path, parser, dta_folder_name)
        except (exceptions.CimageFlatFileNotFound, exceptions.DtaFileNotFound):
            utils.download_cimage(
                base_url,
                cimage_folder_path,
                cimage_folder_path.name,
                dta_folder_name
            )

        if existing is not None:
            return existing

        return create_methods[parser](
            name=name,
            cimage_folder_path=cimage_folder_path,
            fasta_db_path=fasta_db_path,
            dta_folder_name=dta_folder_name
        )

    @classmethod
    def create_from_flatfile(cls,
        name: str,
        cimage_folder_path: pathlib.Path,
        fasta_db_path: pathlib.Path,
        dta_folder_name: str = config.DEFAULT_DTA_FOLDER_NAME
    ):
        cimage_folder_path = pathlib.Path(cimage_folder_path)

        flatfile_path = cls.get_flatfile_in_folder(cimage_folder_path, dta_folder_name)
        flatfile_hash = utils.get_file_hash(flatfile_path)

        dta_file_path = cls.get_dta_file_in_folder(cimage_folder_path, dta_folder_name)
        dta_file_hash = utils.get_file_hash(dta_file_path)

        fasta_db = FastaDatabaseIndex.get_or_create(fasta_db_path)

        data = parse_cimage_flatfile(
            path=flatfile_path,
            dta_file_path=dta_file_path,
            fasta_db_path=fasta_db_path
        )

        with database.atomic():
            experiment = SilacExperiment.create(
                name=name,
                file_hash=flatfile_hash,
                dta_file_hash=dta_file_hash,
                cimage_folder_name=cimage_folder_path.name,
                database=fasta_db
            )
            SilacData.bulk_create([SilacData(**d, experiment=experiment) for d in data], batch_size=1000)
            return experiment

class IsotopExperiment(BaseExperiment):
    """Hold information about individual experiments."""
    file_hash = TextField()
    reactive_residue = FixedCharField(default=config.DEFAULT_REACTIVE_RESIDUE, max_length=1)

class TmtExperiment(BaseExperiment):
    """Hold information about individual experiments."""
    file_hash = TextField(unique=True)
    file_path = TextField()

    @classmethod
    def get_by_file_hash(cls, path):
        file_hash = utils.get_file_hash(path)
        return cls.get_or_none(cls.file_hash == file_hash)

    @classmethod
    def create_from_census(cls,
        name: str,
        census_file_path: pathlib.Path,
        fasta_db_path: pathlib.Path
    ):
        census_file_path = pathlib.Path(census_file_path)
        file_hash = utils.get_file_hash(census_file_path)

        data = parse_census(census_file_path, fasta_db_path)

        with database.atomic():
            experiment = TmtExperiment.create(
                name=name,
                file_hash=file_hash,
                file_path=census_file_path.name,
            )
            TmtData.bulk_create([TmtData(**d, experiment=experiment) for d in data], batch_size=1000)
            return experiment


class TmtIsotopExperiment(TmtExperiment):
    """Hold information about individual experiments."""
    reactive_residue = FixedCharField(default=config.DEFAULT_REACTIVE_RESIDUE, max_length=1)

class BaseData(BaseModel):
    """Holds actual experimental data."""
    experiment = ForeignKeyField(BaseExperiment) # this needs to be overriden if using as base class
    uniprot = TextField(index=True)
    description = TextField()
    symbol = TextField(index=True)
    sequence = TextField(index=True)
    clean_sequence = TextField(index=True)
    unique_sequence = BooleanField(index=True)
    meta = JSONField(null=True)

class IsotopData(BaseData):
    """Holds actual experimental data."""
    cimage_sequence = TextField(index=True)
    residue = IntegerField(index=True)
    possible_residues = JSONField(index=True)
    mass = FloatField()
    ratio = FloatField(index=True)
    stats = TextField()
    run = IntegerField()
    charge = IntegerField()
    segment = IntegerField()
    num_ms2 = IntegerField(null=True)
    link = TextField()

class SilacData(BaseData):
    """Holds actual experimental data."""
    mass = FloatField()
    ratio = FloatField(index=True)
    rsquared = FloatField(index=True)
    charge = IntegerField(index=True)
    segment = IntegerField(index=True)
    num_ms2 = IntegerField(index=True)

class TmtData(BaseData):
    experiment = ForeignKeyField(TmtExperiment)
    channel_intensities = JSONField()
    scan_num = IntegerField(index=True)
    filename = TextField(index=True)

    @property
    def control_intensities(self):
        '''Note this is actually set in during the analysis. This is a dummy property to document this behavior.'''

class TmtIsotopData(TmtData):
    """Holds actual experimental data."""
    experiment = ForeignKeyField(TmtIsotopExperiment)
    residue = IntegerField(index=True)
    possible_residues = JSONField(index=True)

class Blacklist(BaseModel):
    """List of entries to not consider in analyses."""
    data = ForeignKeyField(BaseData)
    description = TextField(null=True)

    class Meta:
        indexes = ((('data_id', 'description'), True),)

class Whitelist(BaseModel):
    """List of entries to consider in analyses irrespective of filters."""
    data = ForeignKeyField(BaseData)
    description = TextField(null=True)
    filter_subset = TextField(null=True)

    class Meta:
        indexes = ((('data_id', 'description'), True),)

class FastaDatabaseIndex(BaseModel):
    """FASTA database index."""
    file_name = TextField(unique=True, index=True)
    file_hash = TextField(unique=True)
    description = TextField(null=True)

    @classmethod
    def get_or_create(cls, path: pathlib.Path, description: str = ''):
        db_path = pathlib.Path(path)
        file_hash = utils.get_file_hash(db_path)

        try:
            with database.atomic():
                db = FastaDatabaseIndex.create(
                    file_name=db_path.name,
                    file_hash=utils.get_file_hash(path),
                    description=description
                )

                FastaDatabase.from_file(db_path, db)

                return db

        except IntegrityError:
            db = FastaDatabaseIndex.get(
                file_name=db_path.name,
                file_hash=file_hash
            )

            return db

    def is_sequence_unique(self, sequence: str) -> bool:
        count = (FastaDatabase
            .select(fn.Count(FastaDatabase.id))
            .where(FastaDatabase.sequence.contains(sequence))
        ).get()

        assert count != 0, 'Sequence not found in FASTA database!'
        return count == 1


class FastaDatabase(BaseModel):
    """Contains a uniprot database used for searching mass spec data."""
    database = ForeignKeyField(FastaDatabaseIndex)
    uniprot = TextField(unique=True, index=True)
    description = TextField(unique=True, index=True)
    sequence = TextField(index=True)

    @classmethod
    def from_file(cls, path: pathlib.Path, database: FastaDatabaseIndex):
        uniprot_db = utils.uniprot_db_from_fasta(path)

        db_models = [
            FastaDatabase(
                database=database,
                uniprot=uniprot,
                description=entry.description,
                sequence=str(entry.seq),
            )
            for uniprot, entry
            in uniprot_db.items()
        ]

        FastaDatabase.bulk_create(db_models, batch_size=250)

class UniprotFlatfile(Model):
    """Contains data from uniprot flatfile."""
    primary_accession = TextField(unique=True, index=True)
    primary_tax_id = IntegerField(index=True)

    # these are dumped directly from data
    gene_name = TextField(index=True)
    organism = TextField(index=True)
    description = TextField(index=True)
    sequence = TextField()
    comments = JSONField()
    cross_references = JSONField()

    # these too but they are more secondary right now
    accessions = JSONField()
    annotation_update = JSONField()
    created = JSONField()
    data_class = TextField()
    entry_name = TextField()
    features = JSONField()
    host_organism = JSONField()
    host_taxonomy_id = JSONField()
    keywords = JSONField()
    molecule_type = JSONField(null=True)
    organelle = JSONField()
    organism_classification = JSONField()
    protein_existence = IntegerField()
    references = JSONField()
    seqinfo = JSONField()
    sequence_length = IntegerField()
    sequence_update = JSONField()
    taxonomy_id = JSONField()

    @classmethod
    def from_file(cls, path):
        db_records = []

        with open(path) as f:
            db = SwissProt.parse(f)

            for record in db:
                db_records.append(cls(
                    primary_accession=record.accessions[0],
                    primary_tax_id=int(record.taxonomy_id[0]),
                    gene_name=record.gene_name,
                    organism=record.organism,
                    description=record.description,
                    sequence=record.sequence,
                    comments=record.comments,
                    cross_references=record.cross_references,
                    accessions=record.accessions,
                    annotation_update=record.annotation_update,
                    created=record.created,
                    data_class=record.data_class,
                    entry_name=record.entry_name,
                    features=record.features,
                    host_organism=record.host_organism,
                    host_taxonomy_id=record.host_taxonomy_id,
                    keywords=record.keywords,
                    molecule_type=record.molecule_type,
                    organelle=record.organelle,
                    organism_classification=record.organism_classification,
                    protein_existence=record.protein_existence,
                    references=[x.__dict__ for x in record.references],
                    seqinfo=record.seqinfo,
                    sequence_length=record.sequence_length,
                    sequence_update=record.sequence_update,
                    taxonomy_id=record.taxonomy_id,
                ))

        with database:
            database.create_tables([cls])

        with database.atomic():
            cls.bulk_create(db_records, batch_size=250)
