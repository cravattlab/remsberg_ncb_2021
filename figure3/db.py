from peewee import *
import datetime
from playhouse.sqlite_ext import (
    SqliteExtDatabase,
    FTSModel,
    SearchField,
    RowIDField,
    JSONField,
)

_pragmas = [
    ('journal_mode', 'wal'),
    ('cache_size', -1000 * 32)
]

ext_database = SqliteExtDatabase('ext.db', pragmas=_pragmas)

def create_tables():
    with ext_database:
        ext_database.create_tables([
            UniprotFlatfile,
            UniprotDatabase
        ])

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

    class Meta:
        database = ext_database

class UniprotDatabase(Model):
    """Contains a uniprot database used for searching mass spec data."""
    uniprot = TextField(unique=True, index=True)
    description = TextField(unique=True, index=True)
    sequence = TextField(index=True)

    class Meta:
        database = ext_database


create_tables()
