import csv
import pandas as pd
import xml.etree.ElementTree as ET
from pandas import DataFrame
from pathlib import Path
import os
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from tqdm import tqdm
import logging


logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO, datefmt='%y-%m-%d %H:%M',
    format='%(asctime)s %(filename)s %(lineno)d: %(message)s')


def read_pdb(protein, peptide):
    parser = PDBParser(QUIET=True)
    protein = parser.get_structure('protein', protein)
    peptide = parser.get_structure('peptide', peptide)
    chainA = protein[0]['A']
    chainB = peptide[0]['B']
    resA = {r.id[1]:r for r in chainA.get_residues()}
    resB = {r.id[1]:r for r in chainB.get_residues()}
    residuesA = {r.id[1]:r.resname for r in chainA.get_residues()}
    residuesB = {r.id[1]:r.resname for r in chainB.get_residues()}
    return resA, resB, residuesA, residuesB

def read_xml(protein_xml):
    TM_region = {}
    EC_domain = {}
    IC_domain = {}
    tree = ET.parse(protein_xml)
    ns = {'xmlns':"http://uniprot.org/uniprot"}

    root = tree.getroot()
    for feature in root.findall('xmlns:entry/xmlns:feature',namespaces=ns):
        # print(feature)
        if feature.get('type') == 'transmembrane region':
            begin = int(feature.find('xmlns:location/xmlns:begin',namespaces=ns).get('position'))
            end = int(feature.find('xmlns:location/xmlns:end',namespaces=ns).get('position'))
            for i in range(begin, end+1):
                TM_region[i] = i
        elif feature.get('type') == 'topological domain':
            desc = feature.get('description')
            begin = int(feature.find('xmlns:location/xmlns:begin',namespaces=ns).get('position'))
            end = int(feature.find('xmlns:location/xmlns:end',namespaces=ns).get('position'))
            if desc == 'Extracellular':
                for i in range(begin, end+1):
                    EC_domain[i] = i
            elif desc == 'Cytoplasmic':
                for i in range(begin, end+1):
                    IC_domain[i] = i
    return TM_region, EC_domain, IC_domain

# 比较A链残基与上述域的关系,计算两链间残基接触
def compair_contact_residue(residuesA,resA, resB, TM_region, EC_domain, IC_domain, cutoffs):          
    TMA = {k:residuesA[k] for k in TM_region if k in residuesA}
    ECA = {k:residuesA[k] for k in EC_domain if k in residuesA} 
    ICA = {k:residuesA[k] for k in IC_domain if k in residuesA}

    contacts_pair = []
    cutoff = cutoffs

    for r1 in resA.values():
        for r2 in resB.values():
            distance = float(r1['CA'] - r2['CA'])
            if distance < cutoff:
                contacts_pair.append((r1.id[1], r2.id[1] , distance))
    return contacts_pair,TMA, ECA, ICA

# 统计不同域的接触频率  
def calc_rate(contacts_pair, TMA, ECA, ICA):
    TM, EC, IC = 0, 0, 0
    for pairs in contacts_pair:
        rA, rB,distance = pairs
        if rA in TMA.keys(): 
            TM += 1
        if rA in ECA.keys():
            EC += 1
        if rA in ICA.keys():
            IC += 1
    if len(contacts_pair) == 0:
        TM_rate = 0
        EC_rate = 0
        IC_rate = 0
    else:
        TM_rate = '{:.3f}'.format(TM/len(contacts_pair))
        EC_rate = '{:.3f}'.format(EC/len(contacts_pair))
        IC_rate = '{:.3f}'.format(IC/len(contacts_pair))
    return TM_rate, EC_rate, IC_rate


# 输出B链序列及统计率  
def output_result(output_path, residuesB, TM_rate, EC_rate, IC_rate):
    with open(output_path, 'a+', newline='') as f:
        writer = csv.writer(f) 
        seqB = seq1(''.join(str(res) for res in residuesB.values()))
        writer.writerow([ seqB, TM_rate, EC_rate, IC_rate])


def filter_EC(output_path):
    output_path=Path(output_path)
    df = pd.read_csv(output_path)
    df1=DataFrame()
    if (df['EC_rate'] > 0).any():
        df1 = df[df['EC_rate'] > 0]['Sequence']
        df1.to_csv(output_path.parent/"IN_EC_domain.csv", header=["Sequence"])


def run(protein, protein_xml, output_path, cutoffs, peptide_dir):
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Sequence', 'TM_rate', 'EC_rate', 'IC_rate'])

    for peptide_filename in tqdm(os.listdir(peptide_dir)):
        if  peptide_filename.endswith('.pdb'):
            peptide = os.path.join(peptide_dir, peptide_filename)
            resA,resB,residuesA,residuesB = read_pdb(protein, peptide)
            TM_region, EC_domain, IC_domain = read_xml(protein_xml)
            contacts_pair,TMA, ECA, ICA = compair_contact_residue(residuesA,resA,resB, TM_region, EC_domain, IC_domain, cutoffs)
            TM_rate, EC_rate, IC_rate = calc_rate(contacts_pair, TMA, ECA, ICA)
            output_result(output_path, residuesB, TM_rate, EC_rate, IC_rate)
            # break
    filter_EC(output_path)


if __name__ == '__main__':
    protein = '/mnt/nas1/lanwei-125/MC5R/Sequence/8inr-MC5R.pdb'
    MC5R=Path('/mnt/nas1/lanwei-125/MC5R/Sequence/MC5R.xml')
    output_path= "/mnt/nas1/lanwei-125/test/output.csv"
    cutoffs = 5
    # peptide ="/mnt/nas1/lanwei-125/MC5R/HPEP_best_pose/MC5R-peptide-s-control-local/DFV.pdb"
    peptide_dir='/mnt/nas1/lanwei-125/MC5R/dock/HPEP_best_pose/'
    run(protein,MC5R,output_path, cutoffs, peptide_dir)
