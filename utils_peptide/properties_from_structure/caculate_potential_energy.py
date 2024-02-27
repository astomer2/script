import os
import shutil
from dataclasses import dataclass
from typing import Sequence
from pathlib import Path

import numpy as np
from openmm import *
from openmm.app import *
from openmm.unit import *
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from openmm import unit as u

from tqdm import tqdm
import sys
sys.path.append(os.path.abspath("."))
from utils_peptide.properties_from_structure.pdb_fixer import fixed_pdb_file

@dataclass
class EnergyUnit:
    complex_energy: float
    protein_energy: float
    peptide_energy: float
    diff_energy: float


@dataclass
class Energy:
    file_name: str
    raw: EnergyUnit
    minimized: EnergyUnit 


def copy_and_create_directory(pdb_path: str or Path) -> tuple[Path, Path]:
    """
    Copy the PDB file to a new directory and create the directory if it does not exist.

    Args:
        pdb_path (str or Path): The path to the PDB file.

    Returns:
        tuple[Path, Path]: A tuple containing the path to the new directory and the path to the copied PDB file.
    """
    pdb_path = Path(pdb_path)
    if pdb_path.suffix == ".pdb":
        workdir = pdb_path.with_suffix("") 
        workdir.mkdir(parents=True, exist_ok=True)
        copy_pdb = pdb_path.with_suffix("")/pdb_path.name
        shutil.copy(pdb_path, workdir)

        return workdir,copy_pdb


def read_pdb_chain( pdb_path:Path or str) -> dict:

    chains = {}
    input_pdb =Path(pdb_path)
    with open(input_pdb) as f:
        for line in f:
            if line.startswith("ATOM"):
                chain_id = line[21]
                chains.setdefault(chain_id, []).append(line)
    return chains


def write_pdb_chain(chains: dict, workdir: Path or str, minimize = False) -> tuple[Path, Path]:
    """
    Write PDB chain files.

    Args:
        chains (dict): A dictionary of chain IDs and their corresponding lines.
        workdir (str): The directory where the PDB chain files will be written.
        minimize (bool, optional): Whether to minimize the protein and peptide files. Defaults to False.

    Returns:
        tuple: A tuple containing the paths to the protein and peptide files.
    """

    workdir = Path(workdir)
    if minimize:
        protein = workdir / "protein_minimized.pdb"
        peptide = workdir / "peptide_minimized.pdb"
    else:
        protein = workdir / "protein.pdb"
        peptide = workdir / "peptide.pdb"
    chain_ids = len(chains)

    if chain_ids >= 2:
        all_chains = list(chains.values())

        # 设置最后一个值为peptide，其余为protein
        peptide_lines = all_chains[-1]
        protein_lines = all_chains[:-1]

        with open(peptide, "w") as f:
            f.write("".join(peptide_lines))

        with open(protein, "w") as f:
            for lines in protein_lines:
                f.write("".join(lines))
                f.write("TER\n")  # 添加TER行，表示链的结束
    else:
        pass
    return protein, peptide

def openmm_minimize(pdb_file: Path) -> Path:
    """
    Opens a PDB file and performs an energy minimization using OpenMM.

    Parameters:
        pdb_file (Path): The path to the PDB file.

    Returns:
        Path: The path to the minimized PDB file.
    """
    pdb = PDBFile(str(pdb_file))
    forcefield = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=NoCutoff,
        nonbondedCutoff=10 * angstrom,
        constraints=HBonds,
    )
    integrator = LangevinMiddleIntegrator(
        310 * kelvin, 1 / u.picosecond, 0.002 * u.picoseconds 
    )    

    platform = Platform.getPlatformByName("OpenCL")
    simulation = Simulation( pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)

    simulation.minimizeEnergy(maxIterations=100, tolerance=0.01)
    positions = simulation.context.getState(getPositions=True).getPositions()

    file_name = pdb_file.parent
    minize_pdb = file_name/(file_name.stem+"_minimized.pdb")
    PDBFile.writeFile(simulation.topology, positions, open(minize_pdb, 'w'))

    return minize_pdb


def calculate_potential_energy(pdb_file: Path) -> float:

    pdb = PDBFile(str(pdb_file))
    forcefield = ForceField("amber14-all.xml", 'implicit/gbn2.xml')
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=NoCutoff,
        nonbondedCutoff=10 * angstrom,
        constraints=HBonds,
    )
    integrator = LangevinMiddleIntegrator(
        310 * kelvin, 1 / u.picoseconds, 0.002 * u.picoseconds
    )
    platform = Platform.getPlatformByName("OpenCL")
    simulation = Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy()
    state = simulation.context.getState(getEnergy=True)
    potentialEnergy = state.getPotentialEnergy()
    return potentialEnergy.value_in_unit(kilocalories_per_mole)

def process_tmp_pdb(workdir:Path,complex_pdb:Path, ) -> tuple[Path, Path, Path]:
    chains = read_pdb_chain(complex_pdb)
    protein, peptide = write_pdb_chain(chains, workdir)
    fixed_complex_pdb = fixed_pdb_file(complex_pdb)
    fixed_protein = fixed_pdb_file(protein)
    fixed_peptide = fixed_pdb_file(peptide)
    return fixed_complex_pdb, fixed_protein, fixed_peptide

def calculate_raw_energy(complex_pdb: Path, workdir: Path) -> tuple[float, float, float, float]:

    fixed_complex_pdb, fixed_protein, fixed_peptide = process_tmp_pdb(workdir,complex_pdb )

    complex_energy = round(calculate_potential_energy(fixed_complex_pdb), 3)
    protein_energy = round(calculate_potential_energy(fixed_protein), 3)
    peptide_energy = round(calculate_potential_energy(fixed_peptide), 3)
    diff_energy = round((complex_energy - protein_energy - peptide_energy), 3)

    return  complex_energy, protein_energy, peptide_energy, diff_energy

def calc_minize_energy(complex_pdb: Path, workdir: Path) -> tuple[float, float, float, float]:
    """
    Calculate the minimized energy of a complex protein-peptide system.

    This function takes a complex protein-peptide PDB file and a working directory as input.
    It performs energy minimization on the complex PDB file using the OpenMM library.
    The minimized PDB file is then processed to extract the protein and peptide chains separately.
    The protein and peptide chains are written to separate PDB files in the working directory.
    Fixed versions of the protein and peptide PDB files are created by removing any flexible parts.
    The potential energy of the complex, protein, and peptide structures is then calculated using the calculate_potential_energy() function.
    The energies are rounded to 3 decimal places.
    The difference in energy between the complex and the sum of the protein and peptide energies is also calculated.
    The complex energy, protein energy, peptide energy, and energy difference are returned as a tuple.

    Parameters:
    - complex_pdb (Path): The path to the complex protein-peptide PDB file.
    - workdir (Path): The path to the working directory.

    Returns:
    - tuple[float, float, float, float]: A tuple containing the complex energy, protein energy, peptide energy, and energy difference.

    """
    minimize_pdb = openmm_minimize(complex_pdb)
    minimize_chains = read_pdb_chain(minimize_pdb)
    minimize_protein, minimize_peptide = write_pdb_chain(minimize_chains, workdir, minimize=True)

    fixed_protein = fixed_pdb_file(minimize_protein)
    fixed_peptide = fixed_pdb_file(minimize_peptide)

    complex_energy = round(calculate_potential_energy(minimize_pdb), 3)
    protein_energy = round(calculate_potential_energy(fixed_protein), 3)
    peptide_energy = round(calculate_potential_energy(fixed_peptide), 3)
    diff_energy = round((complex_energy - protein_energy - peptide_energy), 3)

    return complex_energy, protein_energy, peptide_energy, diff_energy    

def calc_complex_energy(pdb_path: Path or str) -> Energy:
    """
    Calculate the complex energy of a protein structure.

    Parameters:
        pdb_path (Path or str): The path to the PDB file of the protein structure.

    Returns:
        Energy: An object containing the calculated energy of the protein structure.
        
    Example output:                        
        Energy(file_name='ESPLKG_top1', 
                raw=EnergyUnit(complex_energy=-3501.436, 
                    protein_energy=-3221.508, 
                    peptide_energy=-265.188, 
                    diff_energy=-14.74), 
                minimized=EnergyUnit(complex_energy=-3500.722, 
                    protein_energy=-3214.614, 
                    peptide_energy=-265.176, 
                    diff_energy=-20.932))
    """

    workdir, complex_pdb = copy_and_create_directory(pdb_path)


    raw_energies = calculate_raw_energy(complex_pdb, workdir)
    min_energies = calc_minize_energy(complex_pdb, workdir)
    energy = Energy(
        complex_pdb.stem,
        raw=EnergyUnit(*raw_energies),
        minimized=EnergyUnit(*min_energies),
    )
    #shutil.rmtree(workdir)
    return energy


def batch_calc_and_rank_by_raw_diff_energy(pdb_files: Sequence[str])-> tuple[np.ndarray, list]:
    """
    Calculate the raw difference energy for each PDB file in the given sequence of file paths.
    
    Args:
        pdb_files (Sequence[str]): A sequence of file paths to PDB files.
        
    Returns:
        Tuple[np.ndarray, List]: A tuple containing the sorted indices and sorted energies
                                 based on the raw difference energy.
    """

    energies = []
    for file in tqdm(pdb_files):
        energy = calc_complex_energy(file)
        energies.append(energy)
    raw_diff_energies = [e.raw.diff_energy for e in energies]
    sorted_index = np.argsort(raw_diff_energies)
    sorted_energies = [energies[i] for i in sorted_index]
    return sorted_index, sorted_energies


if __name__ == "__main__":
    complex_pdb_path = "/mnt/nas1/lanwei-125/MC5R/dock/hpep-complex/DEQPL.pdb"
    print(calc_complex_energy(complex_pdb_path))
