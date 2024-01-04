import os, sys
sys.path.append(os.path.abspath('.'))
import math
import pandas as pd
# from pandarallel import pandarallel
# pandarallel.initialize(progress_bar=False)
from rdkit.Chem import Lipinski
from rdkit.Chem.rdmolfiles import MolFromFASTA, MolToSmiles, MolFromSmiles
import numpy as np
from rapidfuzz.string_metric import levenshtein as lev_dist
from pandas import DataFrame, Series
from pathlib import Path
from modlamp.descriptors import GlobalDescriptor, PeptideDescriptor
from icecream import ic
ic.configureOutput(includeContext=True, argToStringFunction=str)
from utils_comm.log_util import get_logger
import logging


logger = get_logger(name=__name__, log_file=None, log_level=logging.INFO)
amino_acids = (
    "A","C","D","E","F","G","H","I","L","M","N","P","K","Q","R","S","T","V","W","Y")
# in dbaasp dataset, d aa is in lower cases and l aa in upper cases.
d_aminoacids = [
    "a","c","d","e","f","g","h","i","l","m","n","p","k","q","r","s","t","v","w","y"]
# there are slightlydifferent hydrophobic_aa, use wiki version in table, not in pic. In wiki, pic and table are diff.
# https://en.wikipedia.org/wiki/Amino_acid
hydrophobic_aa = ["A","I", "L", "V", "M", "W", "F", "Y"]
polar_aa = set(('S', 'T', 'N', 'H', 'Q', 'G'))
special_aa = set(('P', 'C'))
apolar_aa = set(('A', 'L', 'V', 'I', 'M'))
charged_aa = set(('E', 'D', 'K', 'R'))
aromatic_aa = set(('W', 'Y', 'F'))
pos_aa = ['K', 'R']
neg_aa = ['D', 'E']
SEQUENCE = "Sequence"


def find_seqNN(seq, dataframe):
    best_dist = float("inf")
    dists = dataframe["Sequence"].map(lambda seq2 : lev_dist(seq, seq2))
    NNi = np.argmin(dists)
    best_dist = dists.iloc[NNi]
    NN = dataframe["Sequence"].iloc[NNi]
    return best_dist, NN


def seq_to_smiles(seq):
    """
    Args:
        seq: AAHAEINEAGRE
    Returns:
        CC[C@H](C)[C@H](NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](Cc1c[nH]cn1)NC(=O)[C@H](C)NC(=O)[C@H](C)N)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](C)C(=O)NCC(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](CCC(=O)O)C(=O)O
    """
    mol = MolFromFASTA(seq, flavor=True, sanitize = True)
    smiles = MolToSmiles(mol, isomericSmiles=True)
    return smiles


def jaccard_distance(a, b):
    """Estimates the Jaccard distance of two binary arrays based on their hashes.

    Arguments:
      a {numpy.ndarray} -- An array containing hash values.
      b {numpy.ndarray} -- An array containing hash values.

    Returns:
      float -- The estimated Jaccard distance.
    """
    # The Jaccard distance of Minhashed values is estimated by
    return 1.0 - np.float(np.count_nonzero(a == b)) / np.float(len(a))


def calc_neg_count(seq):
    neg = (seq.count('D') + seq.count('E'))
    return neg


def calc_pos_count(seq):
    pos = (seq.count('K') + seq.count('R'))
    return pos


def calc_aa_ratio(seq, aa):
    aa_f = seq.count(aa)/len(seq)
    return aa_f


def calc_hac(smiles):
    mol = MolFromSmiles(smiles)
    hac = Lipinski.HeavyAtomCount(mol)
    return hac


def calc_hydr(seq):
    hydr_count = 0
    for aa in hydrophobic_aa:
        hydr_count += seq.count(aa)
    return hydr_count


def calc_hydrophobic_ratio(seq):
    hydr_count = 0
    for aa in hydrophobic_aa:
        hydr_count += seq.count(aa)
    hydrophobic_ratio = hydr_count / len(seq)
    return hydrophobic_ratio


def hydropatch(seq):
    patch = ""
    patches = []
    for aa in seq:
        if aa in hydrophobic_aa:
            patch+=aa
        else:
            if patch != "":
                patches.append(len(patch))
            patch=""
    if patch != "":
        patches.append(len(patch))
    return np.array(patches)


def calc_hba(smiles):
    mol = MolFromSmiles(smiles)
    hba = Lipinski.NumHAcceptors(mol)
    return hba


def calc_hbd(smiles):
    mol = MolFromSmiles(smiles)
    hbd = Lipinski.NumHDonors(mol)
    return hbd


def mean(patches):
    if len(patches) == 0:
        return 0
    return round(patches.mean(), 2)


def d_aa(seq):
    for aa in d_aminoacids:
        if aa in seq:
            return True
    return False


""" calculates GRAVY
GRAVY (grand average of hydropathy) value for the protein sequences you enter. The GRAVY value is calculated by adding the hydropathy value for each residue and dividing by the length of the sequence (Kyte and Doolittle; 1982).
https://pubmed.ncbi.nlm.nih.gov/7108955/
https://www.bioinformatics.org/sms2/protein_gravy.html
Actually available in PeptideDescriptor
    sequences = ['MQKSPLEKASFISKLA']
    print(cal_GRAVY(sequences[0]))
    peptide_desc = PeptideDescriptor(sequences, 'gravy')  # gravy
    peptide_desc.calculate_global()
    print(peptide_desc.descriptor)
GRAVY (grand average of hydropathy), the index is obtained by dividing the sum of hydropathy values of all parts of a protein or peptide by its sequence length, which determines the hydrophobicity/hydrophilicity of the protein [17]. GRAVY index score (mean hydrophilicity and hydrophobicity) below zero indicates that the protein is globular and hydrophilic, while scores above zero are related to likely membrane proteins and hydrophobicity
亲水性太强的不透皮, Maybe filter too hydrophilic peptides
"""
GRAVY_scales = {
    'A':  1.8, 'R': -4.5, 'N': -3.5,
    'D': -3.5, 'C':  2.5, 'Q': -3.5,
    'E': -3.5, 'G':  -0.4, 'H':  -3.2,
    'I':  4.5, 'L':  3.8, 'K': -3.9,
    'M':  1.9, 'F':  2.8, 'P':  -1.6,
    'S': -0.8, 'T':  -0.7, 'W':  -0.9,
    'Y':  -1.3, 'V':  4.2
}


def cal_GRAVY(seq):
    """  """
    value = 0
    for char in seq:
        value += GRAVY_scales[char]
    value = value / len(seq)
    return value


"""
Calculates a set of properties from a protein sequence:
    - hydrophobicity (according to a particular scale)
    - mean hydrophobic dipole moment assuming it is an alpha-helix.
    - total charge (at pH 7.4)
    - amino acid composition
    - discimination factor according to Rob Keller (IJMS, 2011)
Essentially the same as HeliQuest (reproduces the same values).
Author:
  Joao Rodrigues
  j.p.g.l.m.rodrigues@gmail.com
"""
scales = {'Fauchere-Pliska': {'A':  0.31, 'R': -1.01, 'N': -0.60,
                              'D': -0.77, 'C':  1.54, 'Q': -0.22,
                              'E': -0.64, 'G':  0.00, 'H':  0.13,
                              'I':  1.80, 'L':  1.70, 'K': -0.99,
                              'M':  1.23, 'F':  1.79, 'P':  0.72,
                              'S': -0.04, 'T':  0.26, 'W':  2.25,
                              'Y':  0.96, 'V':  1.22},

          'Eisenberg': {'A':  0.25, 'R': -1.80, 'N': -0.64,
                        'D': -0.72, 'C':  0.04, 'Q': -0.69,
                        'E': -0.62, 'G':  0.16, 'H': -0.40,
                        'I':  0.73, 'L':  0.53, 'K': -1.10,
                        'M':  0.26, 'F':  0.61, 'P': -0.07,
                        'S': -0.26, 'T': -0.18, 'W':  0.37,
                        'Y':  0.02, 'V':  0.54},
          }
# https://proteinstructures.com/sequence/amino-acids/
aa_charge = {'E': -1, 'D': -1, 'K': 1, 'R': 1}


def assign_hydrophobicity(sequence, scale='Fauchere-Pliska'):  # noqa: E302
    """Assigns a hydrophobicity value to each amino acid in the sequence"""
    hscale = scales.get(scale, None)
    if not hscale:
        raise KeyError('{} is not a supported scale. '.format(scale))

    hvalues = []
    for aa in sequence:
        sc_hydrophobicity = hscale.get(aa, None)
        if sc_hydrophobicity is None:
            raise KeyError('Amino acid not defined in scale: {}'.format(aa))
        hvalues.append(sc_hydrophobicity)

    return hvalues


def calculate_moment(array, angle=100):
    """Calculates the hydrophobic dipole moment from an array of hydrophobicity
    values. Formula defined by Eisenberg, 1982 (Nature). Returns the average
    moment (normalized by sequence length)
    uH = sqrt(sum(Hi cos(i*d))**2 + sum(Hi sin(i*d))**2),
    where i is the amino acid index and d (delta) is an angular value in
    degrees (100 for alpha-helix, 180 for beta-sheet).
    """
    sum_cos, sum_sin = 0.0, 0.0
    for i, hv in enumerate(array):
        rad_inc = ((i*angle)*math.pi)/180.0
        sum_cos += hv * math.cos(rad_inc)
        sum_sin += hv * math.sin(rad_inc)
    if len(array) != 0:
        return math.sqrt(sum_cos**2 + sum_sin**2) / len(array)
    else:
        print(array)
        return 0


def calc_charge(sequence, charge_dict=aa_charge):
    """Calculates the charge of the peptide sequence at pH 7.4
    https://proteinstructures.com/sequence/amino-acids/
    Charged amino acids
    The charged amino acids at neutral pH (around 7.4) carry a single charge in the side chain. There are four of them, two basic, lysine (Lys, K) and arginine (Arg, R) with a positive charge at neutral pH, and two acidic, aspartate (Asp, D) and glutamate (Glu, E) carrying a negative charge at neutral pH. A so-called salt bridge is formed by the interaction between positively and negatively charged amino acid side chains. It has been been found to be important for the stabilization of protein three-dimensional structure
    """
    sc_charges = [charge_dict.get(aa, 0) for aa in sequence]
    return sum(sc_charges)


def calc_discrimination(mean_uH, total_charge):
    """Returns a discrimination factor according to Rob Keller (IJMS, 2011)
    A sequence with d>0.68 can be considered a potential lipid-binding region.
    """
    d = 0.944*mean_uH + 0.33*total_charge
    return d


def calculate_composition(sequence):
    """Returns a dictionary with percentages per classes"""
    # Residue character table
    polar_aa = set(('S', 'T', 'N', 'H', 'Q', 'G'))
    special_aa = set(('P', 'C'))
    apolar_aa = set(('A', 'L', 'V', 'I', 'M'))
    charged_aa = set(('E', 'D', 'K', 'R'))
    aromatic_aa = set(('W', 'Y', 'F'))

    n_p, n_s, n_a, n_ar, n_c = 0, 0, 0, 0, 0
    for aa in sequence:
        if aa in polar_aa:
            n_p += 1
        elif aa in special_aa:
            n_s += 1
        elif aa in apolar_aa:
            n_a += 1
        elif aa in charged_aa:
            n_c += 1
        elif aa in aromatic_aa:
            n_ar += 1

    return {'polar': n_p, 'special': n_s,
            'apolar': n_a, 'charged': n_c, 'aromatic': n_ar}


def analyze_sequence(name=None, sequence=None, window=18, verbose=False):
    """Runs all the above on a sequence. Pretty prints the results"""
    w = window
    outdata = []  # for csv writing

    # Processing...
    seq_len = len(sequence)
    print('[+] Analysing sequence {} ({} aa.)'.format(name, seq_len))
    print('[+] Using a window of {} aa.'.format(w))
    for seq_range in range(0, seq_len):

        seq_w = sequence[seq_range:seq_range+w]
        if seq_range and len(seq_w) < w:
            break

        # Numerical values
        z = calc_charge(seq_w)
        seq_h = assign_hydrophobicity(seq_w)
        av_h = sum(seq_h)/len(seq_h)
        av_uH = calculate_moment(seq_h)
        d = calc_discrimination(av_uH, z)

        # AA composition
        aa_comp = calculate_composition(seq_w)
        n_tot_pol = aa_comp['polar'] + aa_comp['charged']
        n_tot_apol = aa_comp['apolar'] + aa_comp['aromatic'] + aa_comp['special']  # noqa: E501
        n_charged = aa_comp['charged']  # noqa: E501
        n_aromatic = aa_comp['aromatic']  # noqa: E501

        _t = [name, sequence, seq_range+1, w, seq_w, z, av_h, av_uH, d,
              n_tot_pol, n_tot_apol, n_charged, n_aromatic]
        outdata.append(_t)

        if verbose:
            print('  Window {}: {}-{}-{}'.format(seq_range+1, seq_range,
                                                 seq_w, seq_range+w))
            print('    z={:<3d} <H>={:4.3f} <uH>={:4.3f} D={:4.3f}'.format(z, av_h,  # noqa: E501
                                                                           av_uH, d))  # noqa: E501
            print('    Amino acid composition')
            print('      Polar    : {:3d} / {:3.2f}%'.format(n_tot_pol, n_tot_pol*100/w))  # noqa: E501
            print('      Non-Polar: {:3d} / {:3.2f}%'.format(n_tot_apol, n_tot_apol*100/w))  # noqa: E501
            print('      Charged  : {:3d} / {:3.2f}%'.format(n_charged, n_charged*100/w))  # noqa: E501
            print('      Aromatic : {:3d} / {:3.2f}%'.format(n_aromatic, n_aromatic*100/w))  # noqa: E501
            print()

    return outdata


def calc_hydr_moment(seq):
    """ NB: only supports upper case!
    https://pubs.rsc.org/en/content/articlelanding/1982/fs/fs9821700109
    The structure of a protein can be analysed in terms of what may be called the “hydrophobic moments” of (1) the entire molecule and (2) of the segments of secondary structure that make up the polypeptide chain. The zeroth moment is defined as the sum of the hydrophobicities of the amino-acid residues of the structure under consideration; it is the analogue of the net charge of a cluster of point charges. The first moment, or hydrophobic dipole moment, is the analogue of the electric dipole moment of a cluster of charges. Just as the electric dipole moment measures the asymmetry of the charge distribution, the hydrophobic dipole moment measures the amphiphilicity (asymmetry of hydrophobicity) of the structure. A large hydrophobic dipole moment indicates that a structure is predominantly hydrophobic on one side and predominantly hydrophilic on the other. A quadrupole hydrophobic moment may be similarly defined. It indicates whether a protein is more hydrophobic in its interior (as for a globular protein in aqueous solution) or at its surface (as for a membrane protein).
    """
    hdr = assign_hydrophobicity(seq, "Eisenberg")
    return calculate_moment(hdr)


def calc_other_length_ratio(row):
    l = row.length
    pos = row.positive_count
    hy = row.hydrophobic
    return (l - (pos+hy)) / l


def calc_pos_length_ratio(row):
    l = row.length
    pos = row.positive_count
    return pos/l


def calc_pos_vs_hydrophobic_ratio(row):
    pos = row.positive_count
    hy = row.hydrophobic
    if hy>0:
        return pos/hy
    else:
        return pos


def write_fastafile(row, folder):
    fasta_seq = row["Sequence"]
    fasta = ">{}\n{}".format(fasta_seq, fasta_seq)
    ID = str(row["ID"])
    name = folder / f'{ID}+.seq'
    with open(name, "w") as output:
        output.write(fasta)


def get_filename(row):
    ID = str(row["ID"])
    name = ID+".seq"
    return name


def fileloc(row, folder):
    ID = str(row["ID"])
    name = str((folder / f'{ID}+.seq').absolute())
    return name


def count_ss(ss, pred = "H"):
    return ss.count(pred)


def fract_ss(ss, pred = "H"):
    if len(ss)!=0:
        return ss.count(pred)/len(ss)
    else:
        return 0


def read_spider(row, folder):
    ss = []
    ID = str(row["ID"])
    name = ID+".seq.i2"
    with open(folder / name) as infile:
        for line in infile.readlines()[1:]:
            line = line.strip()
            line = line.split(" ")
            ss.append(line[2])
    return ss


def count_ss(ss, pred = "H"):
    return ss.count(pred)


def fract_ss(ss, pred = "H"):
    if len(ss)!=0:
        return ss.count(pred)/len(ss)
    else:
        return 0


def is_natural_convert_to_upper(seq):
    """ Not diff lower and upper, auto convert to upper """
    if isinstance(seq, str) and seq:
        for aa in seq:
            if aa.upper() not in amino_acids:
                return False
        return True
    return False


def is_natural_only_supper(seq):
    """ If the char is not in upper case, treat as not natural """
    if isinstance(seq, str) and seq:
        for aa in seq:
            if aa not in amino_acids:
                return False
        return True
    return False


def net_charge(sequence):

    '''
    Utility function to calculate net charge of a sequence
    Reference: http://www.chem.ucalgary.ca/courses/351/Carey5th/Ch27/ch27-1-4-2.html

    Parameters
    ----------
    sequence:   str
                peptide sequence

    Returns
    -------
    charge: float
            net charge of sequence
    '''

    acidic = [sequence.count('D'), sequence.count('E'), sequence.count('C'), sequence.count('Y')]
    basic = [sequence.count('R'), sequence.count('K'), sequence.count('H')]

    acidic_pKa = [math.pow(10, 3.65), math.pow(10, 4.25), math.pow(10, 8.37), math.pow(10, 10.46)]
    basic_pKa = [math.pow(10, 10.76), math.pow(10, 9.74), math.pow(10, 7.59)]

    basic_coeff = [x*(1/(x+math.pow(10, 7))) for x in basic_pKa]
    acidic_coeff = [math.pow(10, 7)/(x+math.pow(10, 7)) for x in acidic_pKa]

    charge = - sum(np.multiply(acidic_coeff, acidic)) + sum(np.multiply(basic_coeff, basic))
    return charge


def log_df_basic_info(df: DataFrame):
    logger.info(f'df.shape {df.shape}')
    logger.info(f'df.columns {df.columns.to_list()}')
    logger.info(f'df.head()\n{df.head()}')
    logger.info(f'df.tail()\n{df.tail()}')
    logger.info('df.describe\n%s', df.describe())


def cal_properties_on_df_step1(
        df_all: DataFrame,
        out_file: Path,
        spider3_input_files_dir: Path,
        shortest_length=4,
        only_natrual=True,
        only_upper=True,
    ):
    """  """
    df_all["length"] = df_all.Sequence.map(len)
    logger.info(f'len(df) {len(df_all)}')
    df_all = df_all[df_all['length'] >= shortest_length]
    logger.info(f'len(df) {len(df_all)} only selecting length >= {shortest_length}')

    logger.info('Calculate D_AA')
    df_all["D_AA"] = df_all.Sequence.parallel_map(d_aa)

    logger.info('Calculate is_natural')
    df_all['natural'] = df_all.Sequence.parallel_map(is_natural_convert_to_upper)

    if only_natrual:
        df_all = df_all[df_all['natural']]
    logger.info(f'len(df) {len(df_all)} only selecting only natural')
    # In dbaasp corpus, only upper is only select l aa
    if only_upper:
        df_all = df_all[df_all['Sequence'].parallel_map(lambda x: x.isupper())]
    logger.info(f'len(df) {len(df_all)} only selecting only upper')

    logger.info('Calculate aa_ratio')
    for aa in amino_acids:
        df_all[f"{aa}_ratio"] = df_all.Sequence.parallel_map(lambda x: calc_aa_ratio(x, aa))

    logger.info('Calculate positive_count')
    df_all["positive_count"] = df_all.Sequence.parallel_map(calc_pos_count)
    df_all["negative_count"] = df_all.Sequence.parallel_map(calc_neg_count)
    df_all["pos_len_ratio"] =  df_all.parallel_apply(calc_pos_length_ratio, axis=1)

    logger.info('Calculate SMILES series')
    df_all["SMILES"] = df_all.Sequence.parallel_map(seq_to_smiles)
    # df_all["MAP4"] = df_all.SMILES.parallel_map(calc_map4)  # which is not used actually.
    df_all["HAC"] = df_all.SMILES.parallel_map(calc_hac)
    df_all["HBA"] = df_all.SMILES.parallel_map(calc_hba)
    df_all["HBD"] = df_all.SMILES.parallel_map(calc_hbd)

    logger.info('Calculate hydrophobic series')
    df_all["hydrophobic"] = df_all.Sequence.parallel_map(calc_hydr)
    df_all["hydrophobic_patches"] = df_all.Sequence.parallel_map(hydropatch)
    df_all["hydrophobic_patches_num"] = df_all.hydrophobic_patches.parallel_map(len)
    df_all["hydrophobic_patches_len"] = df_all.hydrophobic_patches.parallel_map(mean)
    df_all["hydro_len_ratio"] = df_all.parallel_apply(lambda x: x.hydrophobic / x.length, axis=1)
    df_all["HydroMoment"] = df_all.Sequence.parallel_map(calc_hydr_moment)
    df_all["charge"] = df_all.Sequence.parallel_map(calc_charge)
    df_all["hydrophobicity"] = df_all.Sequence.parallel_map(assign_hydrophobicity)
    df_all["av_hydrophobicity"] = df_all["hydrophobicity"].parallel_map(lambda x: sum(x)/len(x))
    df_all["discrimination"] = df_all.parallel_apply(lambda x: calc_discrimination(x.HydroMoment, x.charge), axis=1)
    df_all["other_length_ratio"] =  df_all.parallel_apply(calc_other_length_ratio, axis=1)
    df_all["pos_hy_ratio"] =  df_all.parallel_apply(calc_pos_vs_hydrophobic_ratio, axis=1)

    logger.info('Calculate fasta series for Spider')
    df_all.parallel_apply(write_fastafile, axis=1)
    df_all["SpiderFilename"] = df_all.parallel_apply(get_filename, axis=1)
    df_all["SpiderFileloc"] = df_all.parallel_apply(fileloc, axis=1)
    df_all[["SpiderFilename", "SpiderFileloc"]].to_csv(
        spider3_input_files_dir / "file_list-hem", header=False, index=False, sep=' ')

    log_df_basic_info(df_all)
    logger.info('Calculation ends')
    df_all.to_pickle(out_file)
    return df_all


"""
step2, run SPIDER3-Single_np
http://zhouyq-lab.szbl.ac.cn/download/ SPIDER3-Single_np
FIRST activate p27 environment, with python2.7, numpy and scipy!!
conda activate p27
cd /home/qcdong/sw/bio_chemical/SPIDER3-Single_np
sh run.sh
tail -f log.log file

run.sh
    file_list=/home/qcdong/codes/MLpeptide/cache/calculate_and_plot_properties/spider3_input_files_dir/file_list-hem
    out_dir=/home/qcdong/codes/MLpeptide/cache/calculate_and_plot_properties/spider3_output_dir
    nohup ./impute_script_np_with_arg.sh $file_list $out_dir > log.log 2>&1 &
"""


def cal_properties_from_spider3_step3(df_all: DataFrame, out_file: Path = None):
    """ after step2 run step3, add fraction_PredHelical etc
    """
    logger.info('cal_properties_from_spider3_step3, add fraction_PredHelical etc')
    df_all["SS"] = df_all.apply(read_spider, axis=1)
    df_all["countH"] = df_all.SS.map(count_ss)
    df_all["fraction_PredHelical"] = df_all.SS.map(fract_ss)
    df_all["fraction_PredBetaSheet"] = df_all.SS.map(lambda x : fract_ss(x, "E"))
    df_all["fraction_PredCoil"] = df_all.SS.map(lambda x : fract_ss(x, "C"))
    log_df_basic_info(df_all)
    if out_file:
        df_all.to_pickle(out_file)
    return df_all


def get_physico_chemical_properties(sequences = ['KKLQRSDLLRTK','KKLASCNNIPPR'], debug_print=0,
                                    enable_calculate_more=True, amide=False):
    """
    param amide: {boolean} whether the sequences have an amidated C-terminus. Actually this value is True as default in GlobalDescriptor.

    Comparison with local calculators for AMP projects.
        1, the local HydroMoment calculation use different eisenberg scale with the modlamp tool, and it was used in the final selection feature in the successful AMP paper and codes.
        The modlamp tool calculation maybe more accurate and which is about 2 times of local HydroMoment values, as the Eisenberg scale is about 2 times of the local.
        In real comparison, there two are very similar in KDE plot. So keeping only one is rational.
        2, the "charge" calculation is simpler but also rational, the "Charge_modlamp" is from modlamp

    GlobalDescriptor descriptor featurenames:
        ['Length', 'MW', 'Charge', 'ChargeDensity', 'Isoelectric point', 'InstabilityInd', 'Aromaticity', 'AliphaticInd',
        'BomanInd', 'HydrophRatio'] + ['HydroMoment', 'Gravy]
        len(featurenames) = 10
    Returns:
        numpy.ndarray, 2d shape: (len(sequences), len(featurenames))
    """
    desc = GlobalDescriptor(sequences)
    desc.calculate_all(amide=amide)
    if enable_calculate_more:
        peptide_desc = PeptideDescriptor(sequences, 'eisenberg')
        peptide_desc.calculate_moment()
        peptide_desc_gravy = PeptideDescriptor(sequences, 'gravy')
        peptide_desc_gravy.calculate_global()
        descriptor = np.hstack((desc.descriptor, peptide_desc.descriptor, peptide_desc_gravy.descriptor))
    else:
        descriptor = desc.descriptor

    if debug_print:
        print(desc.featurenames)
        print(descriptor)
        print(descriptor.shape)

    return descriptor


def convert_2d_values_to_row_col(input_array):
    """ There is an existent function in DataFrame to convert 2d ndarray to df,  DataFrame.from_records(input_array)
    Inputs: 
        input_array is a 2d np array or 3d python list. if 2d np array, will convert it to 3d list.
            np.array([[2, 3], [1, 3], [1, 1]]), e.g. from descriptor nparray from modlamp, array shape 3, 2
            python list: [[[2, 3], [1, 3]], [[2, 4], [5, 9]]], array shape 2, 2, 2

    Returns:
        dataFrame
            0  1
        0   2  3
        1   1  3
        2   2  4
        3   5  9
    """
    if isinstance(input_array, np.ndarray):
        # Series cannot read >1 dimension shape numpey, so convert np to list
        input_array = np.expand_dims(input_array, axis=0).tolist()
    input_array = Series(input_array)
    s = input_array.apply(Series, 1).stack()
    s.index = s.index.droplevel(-1)
    out = s.apply(lambda x: Series(x))
    out.reset_index(drop=True, inplace=True)
    return out


def cal_properties_by_modlamp(df_all: DataFrame, out_file: Path = None):
    """ use modlamp to calculate global descriptors
    """
    if isinstance(df_all, (str, Path)):
        df_all = pd.read_csv(df_all)
    logger.info('cal_properties by modlamp to calculate global descriptors')
    sequences = df_all['Sequence'].to_list()
    descriptor = get_physico_chemical_properties(sequences)
    logger.info('descriptor calculation finished and starts to conver to dataframe')
    out = DataFrame.from_records(descriptor)
    out.columns = ['Length', 'Molecular weight(g/mol)', 'Charge', 'ChargeDensity(charge/MW)', 'Isoelectric point', 
                   'InstabilityInd', 'Aromaticity', 'AliphaticInd', 'BomanInd', 'HydrophRatio', 'HydroMoment', 'GRAVY']
    # del out['Length']
    # hydro_len_ratio is the same as HydrophRatio
    # del df_all['hydro_len_ratio']
    # logger.debug(f' {df_all.shape}')
    df_all.reset_index(drop=True, inplace=True)
    df_all = pd.concat([df_all, out], axis=1)
    # logger.debug(f' {df_all.shape}')
    log_df_basic_info(df_all)
    if out_file:
        out_file = Path(out_file)
        if out_file.suffix == '.csv':
            df_all.to_csv(out_file, index=False, sep=',')
        else:
            df_all.to_pickle(out_file)
    return df_all


def select_key_properties_step5(
    df_all: DataFrame, out_file: Path, task_name,
    key_properties = [
        'Sequence', 'length', 'HydroMoment', 'HydrophRatio', 'pos_hy_ratio',
        'charge', 'Charge_modlamp', 'InstabilityInd', 'BomanInd', 'fraction_PredHelical',
        ]
):
    """  Only Left AA, only Natural
    Notes:
        1, the local HydroMoment calculation use different eisenberg scale with the modlamp tool, and it was used in the final selection feature in the successful AMP paper and codes.
        The modlamp tool calculation maybe more accurate and which is about 2 times of local HydroMoment values, as the Eisenberg scale is about 2 times of the local.
        In real comparison, there two are very similar in KDE plot. So keeping only one is rational.
        2, the "charge" calculation is simpler but also rational, the "Charge_modlamp" is from modlamp
        3, fraction_PredHelical from SPIDER3-Single_np results
    """
    logger.info(f'select_key_properties_step5, key_properties num: {len(key_properties)}')
    columns = ['ID', task_name] + key_properties
    df_all = df_all[columns]
    log_df_basic_info(df_all)
    df_all.to_pickle(out_file)
    return df_all


def calc_peptide_screening_score(df:DataFrame):
    """
    Inputs:
        df with 'pep_seq' column
    Returns:
        df with 8 extra columns, 7 separate scores and 1 composition score.

    Consider solubility_screening (2 scores) and stability_screening factors(5 scores).
    solubility_screening is less important than stability screening, mostly for references.
    solubility_ratio = 0.3
    stability_ratio = 0.7

    solubility_score_avg only have 3 values: 0, 0.5, 1, because the 2 sub scores only have 0/1 values.
    stability_multiple_P_adjacent_S have many values.

    https://www.sigmaaldrich.com/HK/zh/technical-documents/technical-article/protein-biology/protein-labeling-and-modification/peptide-stability

    https://www.thermofisher.com/hk/en/home/life-science/protein-biology/protein-biology-learning-center/protein-biology-resource-library/pierce-protein-methods/peptide-design.html
    """

    df['solubility_hydrophobic_ratio'] = df['pep_seq'].map(calc_hydrophobic_ratio)
    df['solubility_low_hydrophobic'] = df['pep_seq'].map(hydrophobic_ratio_below_50_percent_score)
    df['solubility_least_1_charged_aa_for_every_5'] = df['pep_seq'].map(least_1_charged_aa_for_every_5_score)

    df['stability_multiple_C_M_W'] = df['pep_seq'].map(multiple_C_M_W_aas_score)
    df['stability_n_terminal_Q'] = df['pep_seq'].map(n_terminal_Q_score)
    df['stability_n_terminal_N'] = df['pep_seq'].map(n_terminal_N_score)
    df['stability_multiple_P_adjacent_S'] = df['pep_seq'].map(multiple_P_adjacent_S_score)
    df['stability_beta_sheet'] = df['pep_seq'].map(beta_sheet_score)

    df['solubility_score_avg'] = (
        df['solubility_low_hydrophobic'] + df['solubility_least_1_charged_aa_for_every_5']) / 2
    df['stability_score_avg'] = (
        df['stability_multiple_C_M_W'] + df['stability_n_terminal_Q'] + df['stability_n_terminal_N'] +
        df['stability_multiple_P_adjacent_S'] + df['stability_beta_sheet']
    ) / 5

    solubility_ratio = 0.3
    stability_ratio = 0.7
    df['composition_score'] = (
        solubility_ratio * df['solubility_score_avg'] + stability_ratio * df['stability_score_avg'])

    return df


def select_peptide_screening_full_score(df):
    """  """
    stability_selected_df = df[df['stability_score_avg'] == 1]
    logger.info(f'len(stability_selected_df) {len(stability_selected_df)}')
    solubility_selected_df = df[df['solubility_score_avg'] == 1]
    logger.info(f'len(solubility_selected_df) {len(solubility_selected_df)}')
    selected_df = stability_selected_df[stability_selected_df['solubility_score_avg'] == 1]
    logger.info(f'len solubility + stability: {len(selected_df)}')
    return selected_df


def hydrophobic_ratio_below_50_percent_score(seq):
    """ solubility_screening 1. Keep hydrophobic amino acid content below 50% of the total sequence length """
    hydrophobic_ratio = calc_hydrophobic_ratio(seq)
    if hydrophobic_ratio < 0.5:
        logger.debug(f'seq {seq} hydrophobic_ratio < 0.5, return 1')
        return 1
    return 0


def least_1_charged_aa_for_every_5_score(seq):
    """ solubility_screening 2. A rule of thumb in designing soluble peptides.
    at least one charged amino acid for every five amino acids.
    """
    limit = 5
    start = 0
    for aa in seq:
        if aa in charged_aa:
            start = 0
        else:
            start += 1
        if start == limit:
            logger.debug(f'seq {seq} Not one charged amino acid for every five amino acids, return 0')
            return 0
    return 1


def multiple_C_M_W_aas_score(seq):
    """ stability_screening 1.
    Multiple Cys (C), Met (M) or Trp (W) amino acids may be difficult to obtain in high purity partly due to the susceptibility of oxidation and/or side reactions. Choose sequences which minimize these residues or choose conservative replacements for these amino acids.

    Notes: the total count of either C M W are sumed, not calculate the separate aa count sum.

    Returns:
        3 aas of CMW, score 0, e.g. ACDWM.
        2 aas of CWM, score 0.5
        0/1 aa of CWM, score 1
    """
    special_multiple_aas = ['C', 'M', 'W']
    score = multiple_aa_score(seq, special_multiple_aas, consider_limit_num_one=False)
    return score


def n_terminal_Q_score(seq):
    """ stability_screening 2.
    N-terminal Gln (Q) is unstable and may cyclize to pyroglutamate when exposed to the acidic conditions of cleavage.
    TIP: Amidate the N-terminus of the sequence or substitute this amino acid.
    """
    if seq.endswith('Q'):
        return 0
    return 1


def n_terminal_N_score(seq):
    """ stability_screening 3.
    Asparagine (N) has a protecting group that is difficult to remove when placed at the N-terminus of a peptide sequence.
    TIP: Remove the Asn at this location, substitute with another amino acid, or lengthen the peptide by one amino acid residue.
    """
    if seq.endswith('N'):
        return 0
    return 1


def multiple_P_adjacent_S_score(seq):
    """ stability_screening 4.
    Multiple prolines (P) or adjacent serines (S) in a sequence can result in a product that is lower in purity or contains many deletion products. Multiple prolines can also undergo a cis/trans isomerization, resulting in an apparent lower purity product.
    """
    score = multiple_aa_score(seq, 'P')

    # when score == 0, worest; if score > 0, there is possible score becomes 0 when considering adjacent serines (S).
    if score == 0:
        return score
    if 'SS' in seq:
        return 0
    return score


def multiple_aa_score(seq:str, target_aa, consider_limit_num_one:bool = False):
    """ calculate the separate aa count sum. 
    consider_limit_num_one: 
        if True, consider low_limit(1) score.
        if False, no consider low_limit(1) but treat is the same as count 0, all the same score 1.

    uncertain: whether to keep mid_limit or low_limit filter.
    """
    multiple_limit = 3
    mid_limit = 2
    low_limit = 1
    special_count = 0
    for aa in seq:
        if isinstance(target_aa, str):
            if aa == target_aa:
                special_count += 1
        elif isinstance(target_aa, list):
            if aa in target_aa:
                special_count += 1

    if special_count >= multiple_limit:
        score = 0
    elif special_count >= mid_limit:
        score = 0.5
    elif consider_limit_num_one and special_count >= low_limit:
        score = 0.75
    else:
        score = 1
    return score


def beta_sheet_score(seq):
    """ stability_screening 5.
    Beta sheet formation is a concern as it causes incomplete solvation of the growing peptide chain and will result in a higher incidence of deletion sequences in the final product. Avoid sequences that contain multiple or adjacent Val (V), Ile (I), Tyr (Y), Phe (F), Trp (W), Leu (L), Gln (Q), and Thr (T).
    Break the pattern by making conservative replacements, for example, inserting a Gly (G) or Pro (P) at every third residue, replacing Gln (Q) with Asn (N), or replacing Thr (T) with Ser (S).

    Notes: excludes Trp (W) because it already exists in multiple_C_M_W_aas_score(), and another article has no this aa.

    Notes: calculate the separate aa count sum, which is different with multiple_C_M_W_aas_score() which uses total sum.
    """
    beta_sheet_aas = ['V', 'I', 'Y', 'F', 'L', 'Q', 'T']
    for aa in beta_sheet_aas:
        score = multiple_aa_score(seq, aa)
        if score == 0:
            return score

    for aa in beta_sheet_aas:
        aas = aa + aa
        if aas in seq:
            return 0
    return score


def create_full_names_for_final_delivery(df:DataFrame):
    """ Replaces pep_seq/prot_seq as full names """
    if 'prot_seq' in df:
        df.rename(columns={'prot_seq': 'protein_seq'}, inplace=True)
    if 'pep_seq' in df:
        df.rename(columns={'pep_seq': 'peptide_seq'}, inplace=True)
    return df


def calc_3same_aas_count(seq: str):
    """  """
    count = 0
    i = 0
    n_chars = len(seq)
    while i < n_chars:
        char = seq[i]
        if i + 2 < n_chars and seq[i+1] == char and seq[i+2] == char:
            count += 1
            start_i = i + 3
            while start_i < n_chars:
                if seq[start_i] == char:
                    start_i += 1
                else:
                    break
            i = start_i
        else:
            i += 1
    return count


def calc_4same_aas_count(seq: str):
    """  """
    count = 0
    i = 0
    n_chars = len(seq)
    while i < n_chars:
        char = seq[i]
        if i + 3 < len(seq) and seq[i+1] == char and seq[i+2] == char and seq[i+3] == char:
            count += 1
            start_i = i + 4
            while start_i < n_chars:
                if seq[start_i] == char:
                    start_i += 1
                else:
                    break
            i = start_i
        else:
            i += 1
    return count


def has_4same_aas(seq: str):
    """  """
    for i, char in enumerate(seq):
        if i + 3 < len(seq) and seq[i+1] == char and seq[i+2] == char and seq[i+3] == char:
            return True
    return False


def calc_3_same_aas_ratio_seqs(seqs):
    """  """
    count = 0
    for seq in seqs:
        count += calc_3same_aas_count(seq)
    return count / len(seqs)


def calc_4_same_aas_ratio_seqs(seqs):
    """  """
    count = 0
    for seq in seqs:
        count += calc_4same_aas_count(seq)
    return count / len(seqs)


def tiny_test_PeptideDescriptor(tmp_out_file='peptide_utils/peptide_properties_tmp.csv'):
    """  """
    sequences = ['MQKSPLEKASFISKLA', 'AFDGHLKI','KKLQRSDLLRTK','KKLASCNNIPPR', 'IANNEADAISLDGGQVFEAG']
    test_only_gravy = 0
    if test_only_gravy:
        ic(cal_GRAVY(sequences[0]))
        peptide_desc = PeptideDescriptor(sequences, 'gravy')  # gravy
        peptide_desc.calculate_global()
        ic(peptide_desc.descriptor)
        ic(peptide_desc.descriptor.shape)
        return
    df = DataFrame({SEQUENCE: sequences})
    cal_properties_by_modlamp(df, tmp_out_file)
    hydr_moment_values = [calc_hydr_moment(sequence) for sequence in sequences]
    ic(hydr_moment_values)


def test_calc_peptide_quality_score():
    """  """
    ss = ['AEDKRPPPPP', 'AEDKRPPPP']
    input_ss = [[s] for s in ss]
    df = DataFrame(input_ss, columns=['pep_seq'])
    logger.info(f'{df}')

    df = calc_peptide_screening_score(df)
    df.to_csv('1.csv', index=False, sep=',')


if __name__ == "__main__":
    try:
        # test_calc_peptide_quality_score()
        tiny_test_PeptideDescriptor()
        # descriptor = get_physico_chemical_properties(sequences = ['AFDGHLKI','KKLQRSDLLRTK'], debug_print=1)
    except Exception as identifier:
        logger.exception('', exc_info=f'{identifier}')
