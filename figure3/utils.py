from constants import (
    DATA_OUTPUT_PATH,
    DATA_INPUT_PATH,
    FIG_OUTPUT_PATH,
    UNIPROT_SEARCH_DB_PATH
)
from db import (
    ext_database,
    UniprotFlatfile,
    UniprotDatabase
)
from Bio import SwissProt, SeqIO
import datetime


def get_timestamped_report_path(template):
    datestamp = datetime.datetime.now().strftime('%Y%m%d%H%M')
    return DATA_OUTPUT_PATH.joinpath(template.format(datestamp))

def get_path_for_fig(name, datestamp=None):
    if not datestamp:
        datestamp = datetime.datetime.now().strftime('%Y%m%d%H%M')

    return FIG_OUTPUT_PATH.joinpath('{}_{}.svg'.format(name, datestamp))

def get_newest_file(glob='*', base_path=DATA_OUTPUT_PATH):
    newest_file_path = sorted(
        (f for f in base_path.glob(glob)),
        key=lambda x: x.stat().st_ctime, reverse=True
    )[0]

    return newest_file_path

def load_uniprot_flatfiles():
    dbs = [
        'uniprot_sprot_human.dat',
    ]

    for db in dbs:
        uniprot_flatfile_to_sqlite(DATA_INPUT_PATH.joinpath(db))

def create_uniprot_db():
    uniprot_db_to_sqlite(UNIPROT_SEARCH_DB_PATH)

def uniprot_flatfile_to_sqlite(path):
    db_records = []

    with open(path) as f:
        db = SwissProt.parse(f)

        for record in db:
            db_records.append(UniprotFlatfile(
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

    with ext_database:
        ext_database.create_tables([UniprotFlatfile])

    with ext_database.atomic():
        UniprotFlatfile.bulk_create(db_records, batch_size=250)

def uniprot_db_from_fasta(path):
    db = list(SeqIO.parse(str(path), 'fasta'))
    uniprot_db = {}

    for item in db:
        try:
            if not 'Reverse_' in item.id:
                # Reversed entries can have the same id as the parent entry, and thus would
                # lead to overwriting the correct sequence for a given uniprot
                # with the reverse sequence
                _id = item.id.split('|')[1]
                uniprot_db[_id] = item
        except IndexError:
            pass
    return uniprot_db

def uniprot_db_to_sqlite(path):
    uniprot_db = uniprot_db_from_fasta(path)

    db_models = [
        UniprotDatabase(
            uniprot=uniprot,
            description=entry.description,
            sequence=str(entry.seq),
        )
        for uniprot, entry
        in uniprot_db.items()
    ]

    with ext_database:
        ext_database.create_tables([UniprotDatabase])

    with ext_database.atomic():
        UniprotDatabase.bulk_create(db_models, batch_size=250)
