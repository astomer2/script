import os
import random
import re
from collections import defaultdict
from functools import lru_cache
from typing import Sequence, Tuple

from icecream import ic

ic.configureOutput(includeContext=True)

from utils_comm.log_util import logger

SEED = 1
# np.random.seed(SEED)
# random.seed(SEED)
SEQUENCE = "Sequence"
CREATED_LENGTH_NAME = "len"
# 20 natrual
amino_acids = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "L",
    "M",
    "N",
    "P",
    "K",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
]
aminoacids_plus_BX = amino_acids + ["B", "X"]
### Non-containing methionine sequences were preferred. Methionine，简写M，Met, 甲硫氨酸, for diversity, not used in codes.
# often_AMP_aminoacids = aminoacids.copy()
# often_AMP_aminoacids.remove('M')
C_TERMINUS = "cTerminus"
N_TERMINUS = "nTerminus"
LEAST_SEQ_LENGTH = 2
MAX_SEQ_LEN = 50
NOT_HEMOLYTIC_KEY = "isNotHemolytic"
ACTIVITY_KEY = "activity"

not_digits_at_head_end_pat = re.compile(r"^\D*|\D*$")
valid_units = ("µg/ml", "µM")
NAs = ["NA", "na", "Na", "nA", "N/A", "n/a"]
NA = "NA"
full_replace_pat = re.compile(r"\s+|[>=<]")
non_az_at_head_pat = re.compile(r"^[^a-zA-Z]+")
hemolysis_names = ["hemolisis", "hemolysis"]
cytotoxicity_names = ["cytotoxicity", "cell death"]
# However, C-Terminal Amidation is common enough that we make an exception
CTERM_AMIDATION_TERMS = [
    "C-Terminal amidation",
    "C-Terminus: AMD",
    "C-Terminal",
    "C-termianal amidation",
]


def is_cterminal_amidation(mod):
    for term in CTERM_AMIDATION_TERMS:
        if term in mod:
            return True
    return False


def is_seq_len_valid(seq: str):
    if LEAST_SEQ_LENGTH <= len(seq) <= MAX_SEQ_LEN:
        return True
    return False


def is_80percent_natural(seq):
    """ """
    count = 0
    for char in seq:
        if char in amino_acids:
            count += 1
    if count / len(seq) >= 0.8:
        return True
    return False


def is_natural_and_variable_length(seq, min_len=5, max_len=15):
    """Includes the end."""
    if (
        seq
        and len(seq) >= min_len
        and len(seq) <= max_len
        and is_natural_only_upper(seq)
    ):
        return True
    return False


def is_natural_and_length_5_15(seq):
    if seq and len(seq) >= 5 and len(seq) <= 15 and is_natural_only_upper(seq):
        return True
    return False


def is_natural_and_length_5_10(seq):
    if seq and len(seq) >= 5 and len(seq) <= 10 and is_natural_only_upper(seq):
        return True
    return False


def is_natural_and_max_length_30(seq):
    if (
        seq
        and len(seq) >= LEAST_SEQ_LENGTH
        and len(seq) <= 30
        and is_natural_only_upper(seq)
    ):
        return True
    return False


def is_natural_and_2_50_length(seq):
    """2~50, inclusive, is_natural_only_upper"""
    if (
        seq
        and len(seq) >= LEAST_SEQ_LENGTH
        and len(seq) <= MAX_SEQ_LEN
        and is_natural_only_upper(seq)
    ):
        return True
    return False


def is_natural_only_upper(seq):
    """If the char is not in upper case, treat as not natural"""
    if isinstance(seq, str) and seq:
        for aa in seq:
            if aa not in amino_acids:
                return False
        return True
    return False


def is_natural_includ_lower_case(seq):
    """includ_lower_case, that's D aas in DBAASP"""
    if isinstance(seq, str) and seq:
        for aa in seq:
            if aa.upper() not in amino_acids:
                return False
        return True
    return False


def randomChoice(l):
    return l[random.randint(0, len(l) - 1)]


def cal_novelty(generated_seqs, all_pos_neg_seqs, verbose=False):
    novel_seq = []
    if len(generated_seqs) == 0:
        novelty_ratio = 0
        logger.info("cal_novelty input len(generated_seqs) is zero!")
    else:
        for s in generated_seqs:
            if s not in all_pos_neg_seqs:
                novel_seq.append(s)
        novelty_ratio = len(novel_seq) / len(generated_seqs) * 100
    if verbose:
        logger.info(
            "len(novel_seq) %s, novelty_ratio %.1f%%, available pos and neg seq num %s",
            len(novel_seq),
            novelty_ratio,
            len(all_pos_neg_seqs),
        )
    return novel_seq, novelty_ratio


def cal_uniqueness(seqs, verbose=False):
    """Use dict to store the repetition num which is used"""
    unique_seqs_dict = defaultdict(int)
    for s in seqs:
        if len(s) > 0:
            unique_seqs_dict[s] += 1
    uniqueness_ratio = (len(unique_seqs_dict) / len(seqs)) * 100
    if verbose:
        logger.info(
            "len(unique_seqs_dict) %s, uniqueness_ratio %.1f%%, all input seq num %s.",
            len(unique_seqs_dict),
            uniqueness_ratio,
            len(seqs),
        )
    return unique_seqs_dict, uniqueness_ratio


def convert_plus_minus_as_only_plus(mynumber):
    """10±2 as 10+2, returns 12.0"""
    try:
        return sum(map(float, mynumber.split("±")))
    except:
        return float("inf")


def check_active(unit, concentration):
    if (
        (unit == "µM" and concentration < 10)
        or (unit == "nM" and concentration < 10000)
        or (unit == "µg/ml" and concentration < 32)
    ):
        return True


def check_inactive(unit, concentration):
    if (
        (unit == "µM" and concentration > 10)
        or (unit == "nM" and concentration > 10000)
        or (unit == "µg/ml" and concentration > 32)
    ):
        return True


def get_terminus_names(dbaasp_peptide):
    n_terminus = dbaasp_peptide.get(N_TERMINUS)
    if n_terminus:
        n_name = n_terminus["name"]
    else:
        n_name = "nan"
    c_terminus = dbaasp_peptide.get(C_TERMINUS)
    if c_terminus:
        c_name = c_terminus["name"]
    else:
        c_name = "nan"
    return n_name, c_name


def is_valid_terminus(n_name, c_name):
    if (n_name == "nan" or n_name == "ACT") and (c_name == "nan" or c_name == "AMD"):
        return True


def sum_plus_minus(s):
    return sum(map(float, s.split("±")))


def get_concentration(species):
    concentration_str = species["concentration"].replace(" ", "")
    concentration_str = concentration_str.replace("–", "-")
    concentration_str = concentration_str.replace("->", "-")
    concentration_str = concentration_str.replace("0,", "0.")
    concentration_str = concentration_str.replace(",", "")
    # try:
    #     if concentration_str[0] == '<':
    #         pass
    # except Exception as identifier:
    #     logger.info(f'concentration_str {type(concentration_str)} {len(concentration_str)}: {concentration_str}')
    #     raise Exception(identifier)
    if len(concentration_str) == 0:
        concentration = 0
    elif concentration_str[0] == "<":
        if concentration_str[1] == "=":
            concentration = convert_plus_minus_as_only_plus(concentration_str[2:])
        else:
            concentration = convert_plus_minus_as_only_plus(concentration_str[1:])
    elif concentration_str[0] == ">" or concentration_str in NAs:
        concentration = float("inf")
    elif "-" in concentration_str:
        concentrations = concentration_str.split("-")
        concentration = convert_plus_minus_as_only_plus(
            concentrations[0]
        ) + convert_plus_minus_as_only_plus(concentrations[1])
        concentration /= 2
    else:
        concentration = convert_plus_minus_as_only_plus(concentration_str)
    return concentration, concentration_str


@lru_cache(maxsize=None)
def get_lysis_ratio(lysis_value: str):
    """activityMeasureForLysisValue in hemoliticCytotoxicActivities, e.g.
    0% Hemolysis
    11.07-12.70
    LC50
    """
    s = not_digits_at_head_end_pat.sub("", lysis_value.split("%")[0])
    value = convert_str_to_float(s)
    return value


def has_hemolysis(lysis_value):
    """ """
    for hemolysis_name in hemolysis_names:
        if hemolysis_name in lysis_value:
            return True
    return False


@lru_cache(maxsize=None)
def convert_str_to_float(concentration_str):
    if "-" in concentration_str:
        # 12.5-25.0 => 12.5 prefer to keep the lower value, lower value prefer to be toxic
        concentrations = concentration_str.split("-")
        concentration = float(concentrations[0])
    elif "±" in concentration_str:
        # 10.7±4.6 => 10.7
        concentrations = concentration_str.split("±")
        concentration = float(concentrations[0])
    else:
        concentration_str = concentration_str.replace("upto", "")
        concentration = float(concentration_str)
    return concentration


def filter_short_amino_acids(df, shortest_length=5):
    df = df[df.Sequence.map(lambda x: len(x) >= shortest_length)]
    return df


def read_fasta_file(file):
    """Parses a file with FASTA formatted sequences

    Returns:
        A list of tuples with the description and sequence of each sequence in the file,
        that's [(description, seq)]
    """
    if not os.path.isfile(file):
        raise IOError("File not found/readable: {}".format(file))

    sequences = []
    description, cur_seq = None, []
    with open(file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cur_seq:
                    sequences.append((description, "".join(cur_seq)))
                description = line[1:]
                cur_seq = []
            elif line:
                cur_seq.append(line)
    if description and cur_seq:
        sequences.append((description, "".join(cur_seq)))  # last seq

    return sequences


def parse_fasta(fasta_string: str) -> Tuple[Sequence[str], Sequence[str]]:
    """Parses FASTA string and returns list of strings with amino-acid sequences, from alphafold.

    Arguments:
      fasta_string: The string contents of a FASTA file.

    Returns:
      A tuple of two lists:
      * A list of sequences.
      * A list of sequence descriptions taken from the comment lines. In the
        same order as the sequences.
    """
    sequences = []
    descriptions = []
    index = -1
    for line in fasta_string.splitlines():
        line = line.strip()
        if line.startswith(">"):
            index += 1
            descriptions.append(line[1:])  # Remove the '>' at the beginning.
            sequences.append("")
            continue
        elif not line:
            continue  # Skip blank lines.
        sequences[index] += line

    return sequences, descriptions


def a_mutation(seq: str):
    """replace each amino acid with A 丙氨酸 to test its performance, if the peformance
    reduces large, this original amino acid is important
    """
    new_seqs = []
    for i in range(len(seq)):
        if seq[i] != "A":
            new_seq = seq[:i] + "A" + seq[i + 1 :]
            new_seqs.append(new_seq)
    return new_seqs


def all_aa_mutation(seq: str):
    """replace each amino acid with all other amino acids to test its performance."""
    new_seqs = []
    for i in range(len(seq)):
        for aa in amino_acids:
            if seq[i] != aa:
                new_seq = seq[:i] + aa + seq[i + 1 :]
                new_seqs.append(new_seq)
    return new_seqs


def get_specify_aa_position(seq:str,aa:str) -> list:
    """return the position of the specify amino acid in the sequence"""
    positions = [index for index, char in enumerate(seq) if char == aa]
    return positions