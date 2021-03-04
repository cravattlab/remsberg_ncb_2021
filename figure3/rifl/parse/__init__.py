from .census import parse_census
from .cimage import (
    parse_dtaselect,
    parse_flatfile as parse_cimage_flatfile,
    parse_peptide as parse_cimage_peptide,
    parse_protein as parse_cimage_protein,
)


tmt_parsers = {
    'census': parse_census
}

cimage_parsers = {
    'flatfile': parse_cimage_flatfile,
    'peptide': parse_cimage_peptide,
    'protein': parse_cimage_protein
}
