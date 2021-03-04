from rifl.utils import (
    ResidueLocator,
    get_first_float_in_string as find_float
)

import csv
import io
import re

delchars = {ord(c): None for c in map(chr, range(256)) if not c.isalpha()}
delchars_except_mod = {ord(c): None for c in map(chr, range(256)) if not c.isalpha() and not c == '*'}


class ParseCensus():
    def __init__(self):
        self.protein_headers = []
        self.peptide_headers = []

    def parse_file(self, filename):
        headers = []

        with open(filename) as f:
            where_was_i = 0
            line = f.readline()
            # collect the headers first
            while line:
                if line.startswith('H'):
                    headers.append(line.strip().split('\t'))
                    where_was_i = f.tell()
                    line = f.readline()
                else:
                    # roll it back so we can pick up where we left
                    # off and parse the data
                    f.seek(where_was_i)
                    break

            self.protein_headers = headers[-2][2:]
            self.peptide_headers = headers[-1][2:]

            # and now parse the data with a real parser
            reader = csv.reader(f, delimiter='\t')
            return self.parse_io(reader)

    def parse_io(self, reader):
        data = []

        for row in reader:
            if row[0] == 'P':
                protein = dict(zip(self.protein_headers, row[1:]))
                protein['peptides'] = []
                data.append(protein)
            elif row[0] == 'S':
                peptide = dict(zip(self.peptide_headers, row[1:]))
                protein['peptides'].append(peptide)

        return data

def parse_census(filename, fasta_db_path):
    parser = ParseCensus()
    raw = parser.parse_file(filename)

    data = []

    oxidized_methionine_regex = re.compile(f'(M)\([\d\.]+\)')
    symbol_regex = re.compile(r'GN=(\w+)')

    for item in raw:
        symbol = symbol_regex.search(item['DESCRIPTION'])
        protein = {
            'uniprot': item['LOCUS'],
            'description': ' '.join(item['DESCRIPTION'].split('=')[0].split()[:-1]),
            'symbol': symbol[1] if symbol else ''
        }

        for peptide in item['peptides']:
            sequence = peptide['SEQUENCE']
            sequence = oxidized_methionine_regex.sub(r'\1#', sequence)
            # sequence = re.sub(r'\([\d\.]+\)', '*', peptide['SEQUENCE'])
            intensities = {find_float(k): int(v) for k, v in peptide.items() if k.startswith('m/z')}

            data.append({
                'sequence': sequence,
                'clean_sequence': sequence.split('.')[1].translate(delchars),
                'unique_sequence': peptide['UNIQUE'] == 'U',
                'channel_intensities': list(intensities.values()),
                'channel_masses': list(intensities.keys()),
                'scan_num': peptide['ScanNum'],
                'filename': peptide['Filename'],
                # 'meta': {
                #     'raw': peptide
                # },
                **protein
            })

    return data

def parse_census_residue(filename, fasta_db_path, residue):
    parser = ParseCensus()
    raw = parser.parse_file(filename)

    data = []

    locator = ResidueLocator(fasta_db_path, residue=residue)

    reactive_residue_regex = re.compile(f'({residue})\([\d\.]+\)')
    oxidized_methionine_regex = re.compile(f'(M)\([\d\.]+\)')

    for item in raw:
        protein = {
            'uniprot': item['LOCUS'],
            'description': '',
            'symbol': ''
        }

        for peptide in item['peptides']:
            sequence = reactive_residue_regex.sub(r'\1*', peptide['SEQUENCE'])
            sequence = oxidized_methionine_regex.sub(r'\1#', sequence)
            # sequence = re.sub(r'\([\d\.]+\)', '*', peptide['SEQUENCE'])
            intensities = {find_float(k): int(v) for k, v in peptide.items() if k.startswith('m/z')}

            data.append({
                'sequence': sequence,
                'clean_sequence': sequence.split('.')[1].translate(delchars),
                'unique_sequence': peptide['UNIQUE'] == 'U',
                'residue': locator.locate_residue(protein['uniprot'], sequence),
                'possible_residues': locator.get_all_possible_residues(protein['uniprot'], sequence),
                'channel_intensities': list(intensities.values()),
                'channel_masses': list(intensities.keys()),
                # 'meta': {
                #     'raw': peptide
                # },
                **protein
            })

    return data
