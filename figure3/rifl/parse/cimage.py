# modules defined in this project
from rifl.utils import uniprot_db_from_fasta, find_all
import rifl.models as m

# installed modules
from bs4 import BeautifulSoup
import requests

# base modules
from collections import defaultdict
from typing import Literal
import functools
import itertools
import pathlib
import re
import io
import csv
import urllib.parse as urlparse


delchars = {ord(c): None for c in map(chr, range(256)) if not c.isalpha()}
delchars_except_mod = {ord(c): None for c in map(chr, range(256)) if not c.isalpha() and not c == '*'}

# fasta_db_path = 'UniProt_Human_Cravattlab_nonredundant2_98id_11-05-2012_reversed.fasta'
# uniprot_db = uniprot_db_from_fasta(fasta_db_path)


def get_dataset_from_url(combined_dta_url: str, dta_select_url: str):
    """Get raw dataset from  given url."""
    res: Response = requests.get(combined_dta_url)
    res.raise_for_status() # raise exception if the url is not found

    dta_res: Response = requests.get(dta_select_url)
    dta_res.raise_for_status()
    ms2_counts, unique_sequences = parse_dtaselect(dta_res.text)

    return parse_peptide(res.text, ms2_counts, unique_sequences)

def get_dta_select_url(combined_dta_url):
    dta_folder_url = urlparse.urljoin(combined_dta_url, 'dta/')
    res = requests.get(dta_folder_url)
    soup = BeautifulSoup(res.text, 'html.parser')
    dta_select_link = soup.find_all('a', string=re.compile('DTASelect'))[0].attrs['href']
    return dta_folder_url + dta_select_link

def parse_dtaselect(dta_file_path: pathlib.Path):
    with dta_file_path.open('r') as f:
        raw_dta = f.read()

    data = raw_dta.splitlines()

    data_start_row = 0
    sequence_column = 0

    # first we find the row where data starts and where the sequence column is
    for index, row in enumerate(data[:100]):
        if row.startswith('Unique'):
            data_start_row = index + 1
            split_row =  row.strip().split('\t')
            sequence_column = split_row.index('Sequence')
            break

    counts = defaultdict(int)
    unique_sequences = set()

    UNIQUE_SYMBOL = '*'

    for row in data[data_start_row:]:
        if row:
            first_character = row[0]
        else:
            continue

        if first_character in ('*', '\t'):
            try:
                sequence = row.split('\t')[sequence_column]

                counts[sequence] = counts[sequence] + 1
                if first_character == UNIQUE_SYMBOL:
                    unique_sequences.add(sequence)
            # reached the end of the file where summary stats lie
            except:
                break

    return (counts, unique_sequences)

class _ParseCombined():
    header_pattern = re.compile('^\d.+')

    def __init__(self):
        pass

    def parse_file(self, filename):
        with open(filename) as combined_dta:
            return self.parse_io(combined_dta)

    def parse_from_string(self, data):
        return self.parse_io(io.StringIO(data))

    def parse_io(self, raw):
        # ditch headings
        headings = raw.readline()
        data = raw.readlines()
        return self.parse_data(data)
        
    def parse_data(self, data):
        groups = []

        for line in data:
            line_info = [ x.strip() for x in re.split('\t+', line.strip()) if len(x.strip()) ]
            # if we've reached a group header
            if self.header_pattern.match(line):
                groups.append(self._extract_header(line_info))
            # if we've reached a subgroup
            else:
                # append peptide to last group
                groups[-1]['peptides'].append(self._extract_line_data(line_info))

        return groups

    def _extract_header(self, line):
        pass

    def _extract_line_data(self, line):
        pass


class ParseProteinCombined(_ParseCombined):
    def _extract_header(self, line):
        return {
            'uniprot': line[1],
            'description': line[2],
            'symbol': line[3],
            'mean_ratio': float(line[4]),
            'peptides': []
        }

    def _extract_line_data(self, line):
        peptide_keys = [
            'sequence', 'mass', 'ratio', 'stats', 'mean', 'noqp', 'run',
            'charge', 'segment', 'link'
        ]
        return dict(zip(peptide_keys, line))

class ParsePeptideCombined(_ParseCombined):
    def _extract_header(self, line):
        return {
            'sequence': line[1],
            'mean_ratio': line[2],
            'peptides': []
        }

    def _extract_line_data(self, line):
        peptide_keys = [
            'uniprot', 'description', 'symbol', 'sequence',
            'mass', 'ratio', 'stats', 'run', 'charge', 'segment', 'link'
        ]

        return dict(zip(peptide_keys, line))


def parse_flatfile(path: pathlib.Path, dta_file_path: pathlib.Path, fasta_db_path: pathlib.Path):
    data = []

    ms2_counts, unique_sequences = parse_dtaselect(dta_file_path)
    fasta_db = m.FastaDatabaseIndex.get_or_create(fasta_db_path)

    delchars_except_mod = {ord(c): None for c in map(chr, range(256)) if not c.isalpha() and not c == '*'}

    with open(path, 'r') as f:
        # skip first line
        f.readline()

        for line in csv.reader(f, delimiter='\t'):
            sequence = line[4]
            clean_sequence = sequence.split('.')[1].translate(delchars)

            data.append({
                'uniprot': line[1],
                'description': line[2],
                'symbol': line[3],
                'sequence': sequence,
                'clean_sequence': clean_sequence,
                'unique_sequence': sequence in unique_sequences,
                # 'residue_number': get_residue_number(line[1], line[4]) if get_residue_number else None,
                'mass': float(line[5]),
                'charge': int(line[6]),
                'segment': int(line[7]),
                'ratio': float(line[8]),
                'num_ms2': ms2_counts.get(sequence, 0),
                'meta': {
                    'intensity': line[9],
                    'peptide_index': line[0],
                    'num_ms2_peaks': line[10].split('/')[0],
                    'num_candidate_peaks': line[10].split('/')[1],
                    'max_light_intensity': line[10].split('/')[2],
                    'light_noise': line[10].split('/')[3],
                    'max_heavy_intensity': line[10].split('/')[4],
                    'heavy_noise': line[10].split('/')[5],
                    'entry': line[12],
                    'link': line[13],
                },
                'rsquared': float(line[11]),
            })

    return data


def parse_protein(raw_dataset: str, ms2_counts: dict, unique_sequences: set):
    parser = ParseProteinCombined()
    data = []

    for group in parser.parse_io(raw_dataset):
        for peptide in group['peptides']:
            full_sequence = peptide['sequence']
            clean_sequence = full_sequence.split('.')[1].translate(delchars)

            data.append({
                'cimage_sequence': '',
                'clean_sequence': clean_sequence,
                'residue': 0,
                'possible_residues': [],
                'num_ms2': ms2_counts.get(full_sequence),
                'unique_sequence': full_sequence in unique_sequences,
                **peptide
            })

    return data

def parse_peptide(raw_dataset: str, ms2_counts: dict, unique_sequences: set):
    parser = ParsePeptideCombined()
    data = []
    
    delchars_except_mod = {ord(c): None for c in map(chr, range(256)) if not c.isalpha() and not c == '*'}

    for group in parser.parse_from_string(raw_dataset):
        cimage_sequence = group['sequence']

        for peptide in group['peptides']:
            full_sequence = peptide['sequence']
            clean_sequence = full_sequence.split('.')[1].translate(delchars)

            data.append({
                'cimage_sequence': cimage_sequence,
                'clean_sequence': clean_sequence,
                'residue': annotate_residue(peptide['uniprot'], full_sequence),
                'possible_residues': get_all_possible_residues(peptide['uniprot'], clean_sequence),
                'num_ms2': ms2_counts.get(full_sequence),
                'unique_sequence': full_sequence in unique_sequences,
                **peptide
            })

    return data


def annotate_residue(uniprot, sequence):
    if uniprot.startswith('Reverse_'):
        return 0

    protein_sequence = uniprot_db[uniprot].seq

    # take sequence between tryptic terminii and remove all modifications
    # except *, which indicates that the residue is labeled
    sequence = sequence.split('.')[1].translate(delchars_except_mod)
    position_in_sequence = sequence.index('*')
    sequence = sequence.replace('*', '')
    sequence_position_in_protein = protein_sequence.find(sequence)
    residue = sequence_position_in_protein + position_in_sequence

    return residue

def get_all_possible_residues(uniprot, sequence, residue='K'):
    if uniprot.startswith('Reverse_'):
        return [0]

    protein_sequence = uniprot_db[uniprot].seq
    sequence_position_in_protein = protein_sequence.find(sequence)

    positions_in_sequence = find_all(sequence, residue)
    positions_in_protein = [sequence_position_in_protein + x + 1 for x in positions_in_sequence]

    return positions_in_protein
