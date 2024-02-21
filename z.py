import os
import csv
import Bio.PDB
from Bio.SeqUtils import seq1
import xml.etree.ElementTree as ET


# 读取PDB文件,获取A链和B链残基
def read_pdb(protein, peptide):
    parser = Bio.PDB.PDBParser(QUIET=True)
    protein = parser.get_structure('protein', protein)
    peptide = parser.get_structure('peptide', peptide)
    chainA = protein[0]['A']
    chainB = peptide[0]['A']
    resA = {r.id[1]:r for r in chainA.get_residues()}
    resB = {r.id[1]:r for r in chainB.get_residues()}
    residuesA = {r.id[1]:r.resname for r in chainA.get_residues()}
    residuesB = {r.id[1]:r.resname for r in chainB.get_residues()}
    return resA, resB, residuesA, residuesB


# 解析UniProt的XML,获取跨膜域和胞内外域信息
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
    TM_rate = '{:.3f}'.format(TM/len(contacts_pair))
    EC_rate = '{:.3f}'.format(EC/len(contacts_pair))
    IC_rate = '{:.3f}'.format(IC/len(contacts_pair))
    return TM_rate, EC_rate, IC_rate

# 输出B链序列及统计率  
def output_result(output_path, residuesB, TM_rate, EC_rate, IC_rate):
#    with open(output_path, 'a', newline='') as f:
#        writer = csv.writer(f)
#        writer.writerow(['Sequence', 'TM_rate', 'EC_rate', 'IC_rate'])  
#        seqB = ''.join(str(res) for res in residuesB.values())
#        writer.writerow([seqB, TM_rate, EC_rate, IC_rate])

    if not os.path.exists(output_path):
        with open(output_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Sequence', 'TM_rate', 'EC_rate', 'IC_rate'])

    with open(output_path, 'a+', newline='') as f:
        writer = csv.writer(f) 
        seqB = seq1(''.join(str(res) for res in residuesB.values()))
        writer.writerow([ seqB, TM_rate, EC_rate, IC_rate])

def run(protein, protein_xml, output_path, cutoffs):
    peptide_dir = '/mnt/sdc/lanwei/ADIPOR/ADCP-whole-top-ranked'
    for peptides in os.listdir(peptide_dir):
        if  peptides.endswith('.pdb'):
            peptide = os.path.join(peptide_dir, peptides)
            resA,resB,residuesA,residuesB = read_pdb(protein, peptide)
            TM_region, EC_domain, IC_domain = read_xml(protein_xml)
            contacts_pair,TMA, ECA, ICA = compair_contact_residue(residuesA,resA,resB, TM_region, EC_domain, IC_domain, cutoffs)
            TM_rate, EC_rate, IC_rate = calc_rate(contacts_pair, TMA, ECA, ICA)
            output_result(output_path, residuesB, TM_rate, EC_rate, IC_rate)

if __name__ == '__main__':
    protein = '/mnt/sdc/lanwei/ADIPOR/计算结合位点/ADIPOR1-A.pdb'

    protein_xml = '/mnt/sdc/lanwei/ADIPOR/计算结合位点/Q96A54.xml'
    output_path= "/mnt/sdc/lanwei/ADIPOR/计算结合位点/output.csv"
    cutoffs = 8
    

                                                                  