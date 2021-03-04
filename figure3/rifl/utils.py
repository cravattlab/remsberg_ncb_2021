from Bio import SeqIO
from bs4 import BeautifulSoup
import pandas as pd
import numpy as np
import requests

from typing import List, Iterable
import urllib.parse as urlparse
import datetime
import dataclasses
import itertools
import json
import hashlib
import inspect
import subprocess
import pathlib
import base64
import difflib
import re

_words_only = re.compile('[\W_]+', re.UNICODE)

class ResidueLocator:
    delchars = {ord(c): None for c in map(chr, range(256)) if not c.isalpha()}

    def __init__(self, fasta_db_path, residue, mod_symbol='*'):
        self.uniprot_db = uniprot_db_from_fasta(fasta_db_path)
        self.symbol = mod_symbol
        self.delchars_except_mod = {ord(c): None for c in map(chr, range(256)) if not c.isalpha() and not c == self.symbol}
        self.residue = residue
    
    def locate_residue(self, uniprot: str, sequence: str) -> int:
        if uniprot.startswith('Reverse_'):
            return 0

        protein_sequence = self.uniprot_db[uniprot].seq

        # take sequence between tryptic terminii and remove all modifications
        # except *, which indicates that the residue is labeled
        try:
            sequence = sequence.split('.')[1].translate(self.delchars_except_mod)
            position_in_sequence = sequence.index(self.symbol)
            sequence = sequence.replace(self.symbol, '')
            sequence_position_in_protein = protein_sequence.find(sequence)
            residue = sequence_position_in_protein + position_in_sequence
        except ValueError:
            # handle unmodified peptides
            residue = 0

        return residue

    def get_all_possible_residues(self, uniprot: str, sequence: str) -> int:
        if uniprot.startswith('Reverse_'):
            return [0]

        sequence = sequence.translate(self.delchars)
        protein_sequence = self.uniprot_db[uniprot].seq
        sequence_position_in_protein = protein_sequence.find(sequence)

        positions_in_sequence = find_all(sequence, self.residue)
        positions_in_protein = [sequence_position_in_protein + x + 1 for x in positions_in_sequence]

        return positions_in_protein


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

def get_git_revision_hash():
    # credit to Yuji Tomita: https://stackoverflow.com/a/21901260/383744
    return subprocess.check_output(['git', 'rev-parse', 'HEAD']).strip()

def hash_fn_source(fn):
    return hashlib.md5(inspect.getsource(fn).encode()).hexdigest()

def get_short_dataframe_hash(df, hash_length=4):
    h = hashlib.sha1(pd.util.hash_pandas_object(df, index=True).values)
    le_hash = _words_only.sub('', base64.b64encode(h.digest()).decode())
    return le_hash[:hash_length]

class JSONSetEncoder(json.JSONEncoder):
    def default(self, obj):
       if isinstance(obj, set):
          return list(obj)
       return json.JSONEncoder.default(self, obj)

def get_datestamp():
    return datetime.datetime.now().strftime('%Y%m%d%H%M')

def get_timestamped_report_path(template, data_output_path, datestamp=None):
    if not datestamp:
        datestamp = get_datestamp()
    return data_output_path.joinpath(template.format(datestamp))

def get_newest_file(base_path: pathlib.Path, glob: str ='*') -> pathlib.Path:
    newest_file_path = next(iter(sorted(
        (f for f in base_path.glob(glob)),
        key=lambda x: x.stat().st_ctime, reverse=True
    )), None)

    return newest_file_path

def pretty_residues(residue_list):
    return ', '.join(map(
        'C{}'.format,
        sorted(set(residue_list))
    ))

def get_longest_overlapping_sequence(sequences: List[str]) -> str:
    sequences = list(set(sequences))

    if len(sequences) == 1:
        return sequences[0]

    d = difflib.SequenceMatcher(a=sequences[0], b=sequences[1])
    match = d.find_longest_match(0, len(sequences[0]), 0, len(sequences[1]))

    overlap_seq = sequences[0]
    d = difflib.SequenceMatcher(a=overlap_seq)

    for seq in sequences[1:]:
        d.set_seq2(seq)
        match = d.find_longest_match(0, len(overlap_seq), 0, len(seq))
        overlap_seq = seq[match[1]:match[1] + match[2]]
        d.set_seq1(overlap_seq)

    return overlap_seq

# https://stackoverflow.com/a/4665027/383744
def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub) # use start += 1 to find overlapping matches

def partition(pred, iterable):
    'Use a predicate to partition entries into true entries and false entries'
    # partition(is_odd, range(10)) --> 0 2 4 6 8   and  1 3 5 7 9
    t1, t2 = itertools.tee(iterable)
    return list(filter(pred, t2)), list(itertools.filterfalse(pred, t1))

def compare_dfs(before, after):
    combined  = pd.concat([before, after])
    diffs = combined.drop_duplicates(keep=False)
    return diffs

def compare_output_excels(before, after):
    return(compare_dfs(
        pd.read_excel(before, index_col=0),
        pd.read_excel(after, index_col=0),
    ))


find_float_regex = re.compile('\d+(?:\.\d+)?')

def get_first_float_in_string(string: str) -> float:
    # https://stackoverflow.com/a/26192945/383744
    return float(find_float_regex.findall(string)[0])

def get_file_hash(path: str, chunk_size: int = 8192) -> str:
    with open(path, 'rb') as f:
        file_hash = hashlib.sha256()
        while chunk := f.read(8192):
            file_hash.update(chunk)

    return file_hash.hexdigest()

def combine_hashes(hashes: Iterable[str]) -> str:
    return hashlib.sha1('|'.join(hashes))

def slash_join(*args):
    '''
    Joins strings to form an URL

    From: https://codereview.stackexchange.com/a/175423/2397
    '''
    return '/'.join(arg.strip('/') for arg in args)


def flatten(list_of_lists: Iterable[Iterable]) -> list:
    '''Yield items from any nested iterable; see Reference.
    
    https://stackoverflow.com/a/40857703/383744
    '''
    for x in list_of_lists:
        if isinstance(x, Iterable) and not isinstance(x, (str, bytes)):
            for sub_x in flatten(x):
                yield sub_x
        else:
            yield x

def get_dta_select_url(dta_folder_url: str):
    res = requests.get(dta_folder_url)
    soup = BeautifulSoup(res.text, 'html.parser')
    dta_select_regex = 'DTASelect-filter_.+\.txt'
    dta_select_link = soup.find_all('a', string=re.compile(dta_select_regex))[0].attrs['href']
    return f'{dta_folder_url}/{dta_select_link}'

def download_cimage(
    base_url: str,
    download_path: pathlib.Path,
    cimage_folder_name: str,
    dta_folder_name: str
) -> bool:
    flatfile_name = 'output_rt_10_sn_2.5.to_excel.txt'
    combined_dta_name =  f'combined_{dta_folder_name}.txt'

    base = slash_join(base_url, cimage_folder_name)
    dta_folder_url = slash_join(base, dta_folder_name)
    dta_select_url = get_dta_select_url(dta_folder_url)
    combined_dta_url = slash_join(base, combined_dta_name)
    flatfile_url = slash_join(dta_folder_url, 'output', flatfile_name)

    download_path.mkdir(exist_ok=True, parents=True)

    with download_path.joinpath(combined_dta_name).open('w') as f:
        res: requests.Response = requests.get(combined_dta_url)
        res.raise_for_status() # raise exception if the url is not found
        f.write(res.text)

    with download_path.joinpath(dta_select_url.split('/')[-1]).open('w') as f:
        res: requests.Response = requests.get(dta_select_url)
        res.raise_for_status() # raise exception if the url is not found
        f.write(res.text)

    with download_path.joinpath(flatfile_name).open('w') as f:
        res: requests.Response = requests.get(flatfile_url)
        res.raise_for_status() # raise exception if the url is not found
        f.write(res.text)

    return True
