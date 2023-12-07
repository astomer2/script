from pathlib import Path
import pandas as pd
import numpy as np
from pandas import DataFrame
import os, sys
import re, random
import hashlib
from typing import Tuple, List, Sequence
import json
import logging


logger = logging.getLogger()
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(filename)s %(lineno)d: %(message)s',
                    datefmt='%m-%d %H:%M:%S')
RANK_FILENAME = 'ranking_debug.json'


def convert_multimer_CF_input_to_AF_input(CF_input_dir, AF_input_dir):
    """
    CFmultimer input fasta
    >APLKRQ
    AVKFPQLCKFCDVRFSTCDNQKSCMSNCSITSICEKPQEVCVAVWRKNDENITLETVCCDPKLPYCDFILEDAASPKCIMKEKKKPGETFFMCSCSSDECNDNIIFSE:
    APLKRQ

    AFmultimer input fasta
    >3e6t_A
    LDKLHVTSTRPQYVRIKNWGSGEILHDTLHHKATSDFTCKSKSCLGSIMNPKSLTRGPRDKPTPLEELLPHAIEFINQYYGSFKEAKIEEHLARLEAVTKEIET
    >3e6t_X
    ADLAHLPF
    """
    CF_input_dir = Path(CF_input_dir)
    AF_input_dir = Path(AF_input_dir)
    for file in CF_input_dir.iterdir():
        if file.is_file() and file.suffix == '.fasta':
            sequences, descriptions = parse_fasta(file)
            new_seqs = []
            new_seqs_descriptions = []
            for seq, desc in zip(sequences, descriptions):
                seqs = seq.split(':')
                descs = [desc]*2
                assert len(seqs) == len(descs)
                new_seqs.append(seqs)
                new_seqs_descriptions.append(descs)
            AF_file = AF_input_dir / file.name
            with open(AF_file, 'w', encoding='utf-8') as f:
                for seqs, descs in zip(new_seqs, new_seqs_descriptions):
                    for seq, desc in zip(seqs, descs):
                        AF_line = f'>{desc}\n{seq}\n'
                        f.write(AF_line)


def get_hash_sha1(x):
  return hashlib.sha1(x.encode()).hexdigest()


def get_hash(x):
  return hashlib.sha3_256(x.encode()).hexdigest()


def parse_fasta(fasta_string: str) -> Tuple[Sequence[str], Sequence[str]]:
  """Parses FASTA string and returns list of strings with amino-acid sequences.

  Arguments:
    fasta_string: The string contents of a FASTA file.

  Returns:
    A tuple of two lists:
    * A list of sequences.
    * A list of sequence descriptions taken from the comment lines. In the
      same order as the sequences.
  """
  if isinstance(fasta_string, Path) or Path(fasta_string).is_file():
    with open(fasta_string, 'r', encoding='utf-8') as f:
      lines = f.readlines()
  elif isinstance(fasta_string, str) and not Path(fasta_string).is_file():
    lines = fasta_string.splitlines()
  sequences = []
  descriptions = []
  index = -1
  for line in lines:
    line = line.strip()
    if line.startswith('>'):
      index += 1
      descriptions.append(line[1:])  # Remove the '>' at the beginning.
      sequences.append('')
      continue
    elif not line:
      continue  # Skip blank lines.
    sequences[index] += line

  return sequences, descriptions


def parse_pdb_id_chains(pdb_id_chains):
    """
    pdb_id_chains is '5qtu_A_H' or '5qtu_AH', or '3E7G_A'
    """
    items = pdb_id_chains.split('_')
    pdb_id = items[0]
    if len(items) == 3:
        prot_chain, pep_chain = items[1], items[2]
    elif len(items) == 2:
      if len(items[1]) == 2:
        prot_chain, pep_chain = list(items[1])
      elif len(items[1]) == 1:
        prot_chain = items[1]
        logger.info('Mock pep_chain: X')
        pep_chain = 'X'
    return pdb_id, prot_chain, pep_chain


def get_af_best_ranking_score(af_out_dir: Path, score_type):
  """  """

  file = af_out_dir / RANK_FILENAME
  with open(file, 'r', encoding='utf-8') as f:
    rank_data = json.load(f)
    best_rank_score = rank_data[score_type][rank_data['order'][0]]
  return best_rank_score


def get_afm_best_ranking_score(af_out_dir: Path):
  """  """
  file = af_out_dir / RANK_FILENAME
  with open(file, 'r', encoding='utf-8') as f:
    rank_data = json.load(f)
    best_rank_score = rank_data['iptm+ptm'][rank_data['order'][0]]
  return best_rank_score


def get_afm_top_ranking_scores(af_out_dir: Path, top_num=5):
  """  """
  file = af_out_dir / RANK_FILENAME
  if not file.is_file():
    logger.info('file %s is not file, exists %s', file, file.exists())
    return 0, 0, [0]
  with open(file, 'r', encoding='utf-8') as f:
    rank_data = json.load(f)
    afm_top_ranking_scores = []
    for i in range(top_num):
        rank_score = rank_data['iptm+ptm'][rank_data['order'][i]]
        if i == 0:
          afm_top_score = rank_score
        afm_top_ranking_scores.append(rank_score)
  afm_mean_top5_score = sum(afm_top_ranking_scores) / top_num
  return afm_top_score, afm_mean_top5_score, afm_top_ranking_scores


def filename_is_standard_positive(file: Path):
  """  """
  items = file.stem.split('_')
  if len(items) == 3 and len(items[0]) == 4 and len(items[1]) == 1 and len(items[2]) == 1:
    return items
  return False


if __name__ == "__main__":
    pass
