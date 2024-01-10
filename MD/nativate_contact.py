from __future__ import print_function
import csv
import mdtraj as md
import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from contact_map import  ContactFrequency, ContactDifference, ContactTrajectory
from sympy import residue

def read_chain_residues(file_path):
    pdb_path = file_path/ (file_path.stem+'.pdb')
    with open(pdb_path) as f:
        lines = f.readlines()
        
    chains = {}
    for line in lines:
        if line.startswith('ATOM'):
            chain_id = line[21]
            residue_id = line[22:26].strip()
            unique_key = (chain_id, residue_id)

            if unique_key not in chains:
                chains[unique_key] = True

    residue_counts = {}
    for key in chains:
        chain_id = key[0]
        if chain_id not in residue_counts:
            residue_counts[chain_id] = 0
        residue_counts[chain_id] += 1
    residue_range = []
    for chain_id, residue_count in residue_counts.items():
        residue_range.append(residue_count)

    return residue_range


def read_pdb_trajectory(pdb, topology, trajectory):
    top = md.load_prmtop(topology)
    traj = md.load(trajectory,stride=50, top=pdb)
    return top, traj

def select_residues(prmtop, residue_range):
    receptor = prmtop.select(
        f"resid 0 to {residue_range[0]-1} and element != H"
        )
    ligand = prmtop.select(
        f"resid {residue_range[0]} to {residue_range[0]+residue_range[1]-1} and element != H"
        )
    return receptor, ligand

def make_contact_frequency(traj, protein, cutoff=0.45):
    trajectory_contacts = ContactTrajectory(traj, protein,cutoff=0.45)
    freq = trajectory_contacts.contact_frequency()
    #freq.residue_contacts.most_common_idx()
    contact_lists = []
    for i in range(len(freq.residue_contacts.most_common_idx())):
        tup = freq.residue_contacts.most_common_idx()[i]
        if tup[1] >0.5:
            print(list(tup[0])[0],tup[1])
            contact_lists.append(list(tup[0])[0])
    #contact_lists.sort()
    return set(contact_lists)

def calc_contact_frequency(pdb_path, topology, trajectory):
    "returns a set of residues that are in contact with the ligand"
    residue_range =read_chain_residues(pdb_path) 
    top, traj = read_pdb_trajectory(pdb_path, topology, trajectory)
    _, ligand = select_residues(top, residue_range)
    contact_residue_list = make_contact_frequency(traj, ligand)
    return contact_residue_list

def data_path(flie_path):
    top = flie_path / 'fixed.prmtop'
    traj = flie_path / 'fixed.mdcrd'
    pdb_path = flie_path / (flie_path.stem+'.pdb')
    return pdb_path, top, traj

def calc_contact_rate(ref_path, peptdie_MD_paths):
    ref_path = Path(ref_path)
    ref_top,  ref_traj,  ref_pdb_path = data_path(ref_path)
    ref_contact_residue_list = calc_contact_frequency(ref_top,  ref_traj,  ref_pdb_path )
    rates = []
    peptdie_MD_paths = Path(peptdie_MD_paths)
    for flies_path in peptdie_MD_paths.iterdir():
        if flies_path.is_dir():
            pdb_path, top, traj = data_path(flies_path)
            if top and traj and pdb_path is not None:
                contact_residue_list = calc_contact_frequency(pdb_path, top, traj )
                rate = len(contact_residue_list.intersection(ref_contact_residue_list))/len(ref_contact_residue_list)
                a = flies_path, rate
                rates.append(a)
    return rates

if __name__ == '__main__':
    ref_path= '/mnt/nas1/lanwei-125/FGF5/FGF5-pos/FGFR1/'
    peptdie_MD_paths = '/mnt/nas1/lanwei-125/FGF5/disulfide/MD/v1/'
    rates = calc_contact_rate(ref_path, peptdie_MD_paths)

    csv_file_path = '/mnt/nas1/lanwei-125/FGF5/disulfide/MD/v1/contact_rate.csv'
    with open(csv_file_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['File', 'Rate'])
        for flies_path, rate in rates:
            writer.writerow([str(flies_path), rate])
