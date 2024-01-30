import sys
import os
import pickle
import csv
import argparse
import numpy as np
import pandas as pd
from scipy import spatial
from collections import defaultdict
from itertools import product, combinations, permutations
from Bio.PDB import PDBParser, parse_pdb_header
import Bio
import vmd
from vmd import molecule, atomsel
import logging
# Read and return motifs from the cyclic peptide database.

#Find matches between PPI motifs and cyclic peptide database motifs based on RMSD.
def find_match(database_name, mtfs , motif_type, frmsd_threshold):
    if motif_type:
        motif_type = 'consecutive'
    else:
        motif_type = 'not_consecutive'
    cyclo_mtfs = pickle.load(open(database_name, 'rb'))
    database_length = len(cyclo_mtfs)
    logging.info("Number of motifs in the database: {}".format(database_length))
    
    match = dict()
    all_match = []
    matches = []
    d_ppi = []
    
    for pp in mtfs:
        m_count = len(mtfs[pp]["m"])
        d_ppi = []
        for j in range(0, m_count):
            d_ppi_single = mtfs[pp]["m"][j]
            d_ppi.append(d_ppi_single)
        for i in range(0,database_length):
            for c in cyclo_mtfs[i]:
                peptide = c
                for d in (cyclo_mtfs[i][c]):
                    d_cyclo = []
                    for k in range(0, m_count):
                        d_cyclo_single = cyclo_mtfs[i][c][d]["m"][k]
                        d_cyclo.append(d_cyclo_single)
                    sum = 0
                    for r in range(0, m_count):
                        dd = (d_ppi[r]-d_cyclo[r])**2
                        sum = sum + dd
                    rmsd = np.sqrt(sum)
                    if rmsd < frmsd_threshold:
                        rmsd = rmsd.round(decimals=4)
                        match = ( peptide , rmsd,  cyclo_mtfs[i][c][d]["resid"], mtfs[pp]["resid"] )
                        all_match.append(match) 
    return(all_match)