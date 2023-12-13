import os
import shutil
from dataclasses import dataclass
from typing import Sequence

import numpy as np
from openmm import *
from openmm.app import *
from openmm.unit import *
from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile
from simtk import unit as u
from tqdm import tqdm


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


def copy_and_create_directory(pdb_path):
    """
    Copy a PDB file to a new directory with the same name as the file (without the extension).

    :param pdb_path: The path to the PDB file to be copied.
    :type pdb_path: str
    :return: The path to the newly created directory where the PDB file was copied.
    :rtype: str
    """

    if pdb_path.endswith(".pdb"):
        file_path, file_name = os.path.split(pdb_path)
        file_name, _ = os.path.splitext(file_name)  # 使用os.path.splitext获取文件名和扩展名
        workdir = os.path.join(file_path, file_name)
        os.makedirs(workdir, exist_ok=True)
        shutil.copy(pdb_path, workdir)

        return workdir


def read_pdb_chain(workdir, pdb_path):
    """
    Reads a PDB file and extracts the chains.

    Args:
        workdir (str): The working directory where the PDB file is located.
        pdb_path (str): The path to the PDB file.

    Returns:
        dict: A dictionary where the keys are chain IDs and the values are lists of ATOM records belonging to each chain.
    """

    chains = {}
    file_name = pdb_path.split("/")[-1]
    input_pdb = os.path.join(workdir, file_name)
    with open(input_pdb) as f:
        for line in f:
            if line.startswith("ATOM"):
                chain_id = line[21]
                chains.setdefault(chain_id, []).append(line)
    return chains


def write_pdb_chain(chains, workdir, minimize = False):
    """
    Write PDB chain files.

    Args:
        chains (dict): A dictionary of chain IDs and their corresponding lines.
        workdir (str): The directory where the PDB chain files will be written.
        minimize (bool, optional): Whether to minimize the protein and peptide files. Defaults to False.

    Returns:
        tuple: A tuple containing the paths to the protein and peptide files.
    """

    if minimize :
        protein = os.path.join(workdir, "minimized_protein.pdb")
        peptide = os.path.join(workdir, "minimized_peptide.pdb")
    else:
        protein = os.path.join(workdir, "protein.pdb")
        peptide = os.path.join(workdir, "peptide.pdb")
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

def fixed_pdb_file(pdb_file, workdir):
    """
    Generate a fixed PDB file.

    Args:
        pdb_file (str): The path to the PDB file.
        workdir (str): The working directory where the fixed PDB file will be saved.

    Returns:
        a fixed PDB file with the same name as the input PDB file.
    """

    _, complex_name = os.path.split(pdb_file)
    complex_fixed_pdb = os.path.join(workdir, complex_name)

    fixer = PDBFixer(filename= complex_fixed_pdb)

    fixer.findMissingResidues()
    fixer.findNonstandardResidues()

    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(True)

    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    
    PDBFile.writeFile(fixer.topology, fixer.positions, open(complex_fixed_pdb, 'w'))

    return complex_fixed_pdb

def openmm_minimize(pdb_file, workdir):
    """
    Minimizes the energy of a protein structure using OpenMM.
    
    Args:
        pdb_file (str): The path to the PDB file containing the protein structure.
        workdir (str): The directory where the minimized PDB file will be saved.
    
    Returns:
        str: The path to the minimized PDB file.
    """

    pdb = PDBFile(pdb_file)
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
    file_name = pdb_file.split("/")[-1].split(".")[0]
    minize_pdb = os.path.join(workdir, file_name + "_minimized.pdb")
    PDBFile.writeFile(simulation.topology, positions, open(minize_pdb, 'w'))

    return minize_pdb


def calculate_potential_energy(pdb_file):
    """
    Calculate the potential energy of a molecule using the provided PDB file.

    Parameters:
    - pdb_file (str): The path to the PDB file containing the molecular structure.

    Returns:
    - float: The potential energy of the molecule in kilocalories per mole.
    """


    pdb = PDBFile(pdb_file)
    forcefield = ForceField("amber14-all.xml", 'implicit/gbn2.xml')
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
    simulation = Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy()
    state = simulation.context.getState(getEnergy=True)
    potentialEnergy = state.getPotentialEnergy()
    return potentialEnergy.value_in_unit(kilocalories_per_mole)


def calculate_raw_energy(complex_pdb, workdir):
    """
    Calculate the raw energy of a complex using the provided PDB file and work directory.

    Parameters:
        complex_pdb (str): The path to the PDB file of the complex.
        workdir (str): The directory where intermediate files will be stored.

    Returns:
        tuple: A tuple containing the calculated energy values for the complex, protein, peptide, and the difference in energy.
            - complex_energy (float): The energy of the complex.
            - protein_energy (float): The energy of the protein.
            - peptide_energy (float): The energy of the peptide.
            - diff_energy (float): The difference in energy between the complex, protein, and peptide.
    """

    chains = read_pdb_chain(workdir,complex_pdb)
    protein, peptide = write_pdb_chain(chains, workdir)

    fixed_complex_pdb = fixed_pdb_file(complex_pdb, workdir)
    fixed_protein = fixed_pdb_file(protein, workdir)
    fixed_peptide = fixed_pdb_file(peptide, workdir)

    complex_energy = round(calculate_potential_energy(fixed_complex_pdb), 3)
    protein_energy = round(calculate_potential_energy(fixed_protein), 3)
    peptide_energy = round(calculate_potential_energy(fixed_peptide), 3)
    diff_energy = round((complex_energy - protein_energy - peptide_energy), 3)

    return  complex_energy, protein_energy, peptide_energy, diff_energy

def calc_minize_energy(complex_pdb, workdir):
    """
    Calculate the minimized energy of a complex.

    Args:
        complex_pdb (str): The path to the complex PDB file.
        workdir (str): The directory to store intermediate and output files.

    Returns:
        tuple: A tuple containing the complex energy, protein energy, peptide energy, and difference energy.
            - complex_energy (float): The potential energy of the complex.
            - protein_energy (float): The potential energy of the protein.
            - peptide_energy (float): The potential energy of the peptide.
            - diff_energy (float): The difference in energy between the complex and the sum of the protein and peptide energies.
    """
    _, complex_name = os.path.split(complex_pdb)
    fixed_pdb = os.path.join(workdir, complex_name)

    minimize_pdb = openmm_minimize(fixed_pdb, workdir)
    minimize_chains = read_pdb_chain(workdir, minimize_pdb)
    minimize_protein, minimize_peptide = write_pdb_chain(minimize_chains, workdir)

    fixed_protein = fixed_pdb_file(minimize_protein, workdir)
    fixed_peptide = fixed_pdb_file(minimize_peptide, workdir)

    complex_energy = round(calculate_potential_energy(minimize_pdb), 3)
    protein_energy = round(calculate_potential_energy(fixed_protein), 3)
    peptide_energy = round(calculate_potential_energy(fixed_peptide), 3)
    diff_energy = round((complex_energy - protein_energy - peptide_energy), 3)

    return complex_energy, protein_energy, peptide_energy, diff_energy    

def calc_complex_energy(complex_pdb):
    """
    Calculate the complex energy of a given PDB file.

    Parameters:
        complex_pdb (str): The path to the PDB file of the complex.

    Returns:
        Energy: An Energy object containing the calculated energy values.
            - file_name (str): The name of the PDB file.
            - raw (EnergyUnit): The raw energy values.
                - complex_energy (float): The total energy of the complex.
                - protein_energy (float): The energy of the protein.
                - peptide_energy (float): The energy of the peptide.
                - diff_energy (float): The difference in energy.
            - minimized (EnergyUnit): The minimized energy values.
                - complex_energy (float): The total energy of the complex.
                - protein_energy (float): The energy of the protein.
                - peptide_energy (float): The energy of the peptide.
                - diff_energy (float): The difference in energy.

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

    #complex_pdb = str(complex_pdb)
    workdir = copy_and_create_directory(complex_pdb)

    _, complex_name = os.path.split(complex_pdb)
    file_name = complex_name.split(".")[0]

    raw_energies = calculate_raw_energy(complex_pdb, workdir)
    min_energies = calc_minize_energy(complex_pdb, workdir)
    energy = Energy(
        file_name,
        raw=EnergyUnit(*raw_energies),
        minimized=EnergyUnit(*min_energies),
    )
    shutil.rmtree(workdir)
    return energy
    """

    """

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
    complex_pdb_path = "/mnt/sdc/lanwei/TGF/top-pdb/ESPLKG_top1.pdb"
    print(calc_complex_energy(complex_pdb_path))