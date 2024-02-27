""" 通过穷举的方法构建4-7长度的肽库，肽库的规模大概是20^4+20^5+20^6+20^7 = 1,347,360,000
然后可以用stability筛掉stability不好的肽
剩下的来做大规模的screening
"""


import random
from pathlib import Path
import os, sys
sys.path.append(os.path.abspath('.'))
from icecream import ic
ic.configureOutput(includeContext=True, argToStringFunction=str)
from utils_comm.log_util import logger
from functools import cache
from utils_comm.file_util import FileUtil
from utils_peptide.aa_util import get_each_aa_count, get_fraction_values, amino_acids


out_dir = Path('/mnt/sda/bio_drug_corpus/AIpep/random_generate')
out_dir.mkdir(exist_ok=1)
random_5_7_file = out_dir / 'random_5_7.txt'
temp_out_dir = Path(__file__).parent / 'temp_data'
temp_out_dir.mkdir(exist_ok=1, parents=1)
MAX_LINE_NUM_PER_TXT_FILE = 400_000_000


def generate_random_natural_seqs(min_len=5, max_len=7):
    """ the current max_len 8 is not successful, maybe the memory limits. """
    out_file = out_dir / f'random_{min_len}_{max_len}.txt'
    all_aas = []
    for length in range(min_len, max_len+1):
        logger.info(f'length {length} starts')
        aas = generate_all_random_sequences_in_fixed_len(length=length)
        logger.info(len(aas))
        all_aas.extend(aas)
    logger.info(len(all_aas))
    with open(out_file, 'w', encoding='utf-8') as f:
        for aas in all_aas:
            f.write(f'{aas}\n')


@cache
def generate_all_random_sequences_in_fixed_len(chars=amino_acids, length=4):
    """  """
    if length == 1:
        return chars
    seqs = []
    last_seqs = generate_all_random_sequences_in_fixed_len(chars, length-1)
    for char in chars:
        for last_seq in last_seqs:
            seqs.append(char+last_seq)
    return seqs


def generate_with_fixed_ratio_for_each_aa(num=1000, length=15, use_average_frequency=True, fixed_frequency=None,
                                          out_file=None, calc_generated_frequency=False, seed=0):
    """  """
    if fixed_frequency is None and use_average_frequency:
        frequency = [0.05] * len(amino_acids)
    else:
        frequency = fixed_frequency
    ic(frequency)
    random.seed(seed)
    seqs = set()
    while len(seqs) < num:
        sequence = random.choices(amino_acids, k=length, weights=frequency)
        seqs.add(''.join(sequence))
    if out_file:
        FileUtil.write_lines_to_txt(seqs, out_file)
    if calc_generated_frequency:
        aa_count = get_each_aa_count(seqs)
        ratios = get_fraction_values(aa_count, aa_count.keys())
        logger.info('%s', ratios)
    seqs = sorted(seqs)
    return seqs


def generate_8_from_5_7():
    """ len(5-6): 67,200,000, len(seqs_5_7): 1,347,200,000 """
    min_len = 5
    max_len = 7
    seqs_5_7 = FileUtil.read_lines_from_txt(out_dir / f'random_{min_len}_{max_len}.txt')
    ic(len(seqs_5_7))
    # len(seqs_5_7): 1,347,200,000
    seqs_7 = [seq for seq in seqs_5_7 if len(seq) == max_len]
    ic(len(seqs_5_7))
    ic(seqs_5_7[:5], seqs_5_7[-5:])
    seqs_8 = []
    for aa in amino_acids:
        ic(aa)
        for seq in seqs_7:
            seqs_8.append(aa+seq)
        FileUtil.write_lines_to_txt(seqs_8,
                                out_dir / f'random_seq_8_startswith_{aa}.txt')
        seqs_8.clear()


def split_huge_txt_files(input_dir_or_file):
    """  """
    input_dir_or_file = Path(input_dir_or_file)
    if input_dir_or_file.is_file():
        split_dir = input_dir_or_file.parent / f'{input_dir_or_file.stem}-split'
        split_huge_txt_file(input_dir_or_file, split_dir)
    else:
        split_dir = input_dir_or_file.with_name(f'{input_dir_or_file.stem}-split')
        for file in input_dir_or_file.iterdir():
            split_huge_txt_file(file, split_dir)


def split_huge_txt_file(file, split_dir):
    """  """
    ic(file)
    split_dir.mkdir(exist_ok=True)
    seqs = FileUtil.read_lines_from_txt(file)
    s_i = 0
    end_i = MAX_LINE_NUM_PER_TXT_FILE
    while s_i < len(seqs):
        sub_seqs = seqs[s_i: end_i]
        FileUtil.write_lines_to_txt(sub_seqs,
                                split_dir / f'{file.stem}-{s_i}-{end_i}.txt')
        s_i = end_i
        end_i += MAX_LINE_NUM_PER_TXT_FILE

    
if __name__ == "__main__":
    # aa_seqs = generate_random_sequences_in_fixed_len(len=4)
    # generate_random_natural_seqs(min_len=4, max_len=4)
    generate_8_from_5_7()
    # split_huge_txt_files(random_5_7_file)

    # fixed_ratio_for_each_aa_file = temp_out_dir / 'fixed_ratio_for_each_aa.txt'
    # generate_with_fixed_ratio_for_each_aa(out_file=fixed_ratio_for_each_aa_file, calc_generated_frequency=1)
    logger.info('end')