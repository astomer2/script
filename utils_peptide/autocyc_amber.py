import os
import subprocess
from pathlib import Path


def extract_three_letter_code(work_path):
    """
    Extracts the three-letter amino acid sequences, sequence lengths, and file names from .pdb files in a given directory.

    :param work_path: Path to the directory containing the .pdb files.
    :return: A tuple of three lists: sequences, sequence_lengths, and file_names.
    """
    sequence_lengths = []
    sequences = []
    file_names = []
    for file in Path(work_path).glob("*.pdb"):
        with open(file, "r", encoding="utf-8") as f:
            lines = f.readlines()
        sequence = []
        residue_ids = []
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                res_id = line[22:27].strip()
                if res_id not in residue_ids:
                    residue_ids.append(res_id)
                    res_name = line[17:20]
                    sequence.append(res_name)
        sequence_str = " ".join(sequence)
        sequences.append(sequence_str)
        sequence_lengths.append(len(sequence))
        file_names.append(file.stem)
    return sequences, sequence_lengths, file_names


def run_tleap(three_letter_code, sequence_length, file_name):
    """
    Runs tleap with a given sequence, sequence length, and file name.

    :param three_letter_code: A list of three-letter amino acid codes.
    :param sequence_length: The length of the sequence.
    :param file_name: The name of the file (without the .pdb extension).
    """
    formatted_sequence = "{" + " ".join(three_letter_code) + "}"

    tleap_input = f"""
    source leaprc.protein.ff19SB
    temp={formatted_sequence}
    b=loadpdbusingseq {file_name}.pdb temp
    remove b b.1.H2
    remove b b.1.H3
    remove b b.{sequence_length}.OXT
    bond b.1.N b.{sequence_length}.C
    source leaprc.water.opc
    solvateoct b OPCBOX  10.0
    saveAmberParm b {file_name}_prmtop {file_name}_inpcrd
    quit
    """

    with open(f"{file_name}_leap.in", "w") as file:
        file.write(tleap_input)

    subprocess.run(["tleap", "-f", f"{file_name}_leap.in"])


def autocyc_amber(work_path):
    sequences, sequence_lengths, file_names = extract_three_letter_code(work_path)
    os.chdir(work_path)  # Change work directory to work_path
    for sequence, sequence_length, file_name in zip(
        sequences, sequence_lengths, file_names
    ):
        run_tleap(sequence.split(), sequence_length, file_name)


if __name__ == "__main__":
    work_path = "/home/yuliu/nas/yuliu/repos/peptide-deploy/utils_peptide/linear_pep_only"
    autocyc_amber(work_path)
    print("pass")
