import os
import shutil
from dataclasses import dataclass
from typing import Sequence

import numpy as np
from openmm import *
from openmm.app import *
from openmm.unit import *
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
    minimized: EnergyUnit = None


def copy_and_create_directory(pdb_path):
    """
    将PDB文件复制到工作目录，并创建工作目录
    :param pdb_path: PDB文件路径
    """
    if pdb_path.endswith(".pdb"):
        file_path, file_name = os.path.split(pdb_path)
        file_name, _ = os.path.splitext(file_name)  # 使用os.path.splitext获取文件名和扩展名
        workdir = os.path.join(file_path, file_name)
        os.makedirs(workdir, exist_ok=True)
        shutil.copy(pdb_path, workdir)

        return workdir


def read_pdb_chain(workdir):
    """
    读取PDB文件中的链信息
    :param workdir: 工作目录
    """
    chains = {}
    for pdb_name in [f for f in os.listdir(workdir) if f.endswith(".pdb")]:
        input_pdb = os.path.join(workdir, pdb_name)

        with open(input_pdb) as f:
            for line in f:
                if line.startswith("ATOM"):
                    chain_id = line[21]
                    chains.setdefault(chain_id, []).append(line)
    return chains


def write_pdb_chain(chains, workdir):
    """
    将PDB文件中的链信息写入到PDB文件中
    :param chains: 链信息
    :param workdir: 工作目录
    """

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


def calculate_potential_energy(pdb_file):
    """
    计算PDB文件的势能
    :param pdb_file:
    """

    pdb = PDBFile(pdb_file)
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=NoCutoff,
        nonbondedCutoff=10 * angstrom,
        constraints=HBonds,
    )
    integrator = LangevinMiddleIntegrator(
        310 * kelvin, 1 / u.picosecond, 0.004 * u.picoseconds
    )

    platform = Platform.getPlatformByName("OpenCL")
    simulation = Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy()
    state = simulation.context.getState(getEnergy=True)
    potentialEnergy = state.getPotentialEnergy()
    return potentialEnergy.value_in_unit(kilocalories_per_mole)


def calculate_differnet_energy(complex_pdb_file_path, protein_path, peptide_path):
    """
    计算多肽、蛋白、复合体以及差值的势能值,单位kcal/mol

    Args:
        complex_pdb_file_path: 输入的复合体的pdb文件路径
        protein_path, peptide_path: 临时生成的蛋白、多肽的路径
    """
    file_name = complex_pdb_file_path.split("/")[-1].split(".")[0]
    complex_energy = round(calculate_potential_energy(complex_pdb_file_path), 3)
    protein_energy = round(calculate_potential_energy(protein_path), 3)
    peptide_energy = round(calculate_potential_energy(peptide_path), 3)
    diff_energy = complex_energy - protein_energy - peptide_energy
    energy = Energy(
        file_name,
        raw=EnergyUnit(complex_energy, protein_energy, peptide_energy, diff_energy),
    )
    return energy


def calc_complex_energy(complex_pdb: str):
    """
    Advice to rank the complex energy by raw diff_energy in output.

    Args:
        pdb: pdb文件的路径

    Returns:
        Energy(file_name='WPIHHVT_monomer', raw=EnergyUnit(complex_energy=-2081.293, protein_energy=-1946.171, peptide_energy=-17.132, diff_energy=-117.99000000000007), minimized=None)

    """
    complex_pdb = str(complex_pdb)
    workdir = copy_and_create_directory(complex_pdb)
    chains = read_pdb_chain(workdir)
    protein, peptide = write_pdb_chain(chains, workdir)
    energy = calculate_differnet_energy(complex_pdb, protein, peptide)
    shutil.rmtree(workdir)
    return energy


def batch_calc_and_rank_by_raw_diff_energy(pdb_files: Sequence[str]) -> list[Energy]:
    """
    Returns:
        sorted_index: the index of pdb_files sorted by raw_diff_energies, always in ascending order.
        sorted_energies: the Energy list sorted by raw_diff_energies

    Note: pdb_files should be a list of full paths, but the file_name in Energy is just the filename.
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
    complex_pdb_path = "/mnt/nas1/lanwei-125/IL8/v4/structure/WPIHHVT_monomer.pdb"
    print(calc_complex_energy(complex_pdb_path))