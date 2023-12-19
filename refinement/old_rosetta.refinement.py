import os
import shutil
import pandas as pd
from pathlib import Path
import logging
from icecream import ic

logging.basicConfig(filename='logfile.log', level=logging.INFO)

def read_pdb(dir_path):
    cluster_dir_names = []
    pdb_names = []

    for dir_entry in Path(dir_path).iterdir():
        if dir_entry.is_dir() and dir_entry.name.startswith("cluster"):
            cluster_dir_names.append(dir_entry.name)

    for cluster_dir_name in cluster_dir_names:
        pdb_path = Path(dir_path) / cluster_dir_name
        for file_entry in pdb_path.iterdir():
            if file_entry.suffix == ".pdb":
                pdb_names.append(file_entry.stem)

    return cluster_dir_names, pdb_names

def preprocess_pdb(protein_pdb, cluster_dir_names, dir_path):
    protein_lines = []
    with open(protein_pdb) as f:
        protein_liners = f.readlines()
        for protein_line in protein_liners:
            if protein_line.startswith("ATOM") or protein_line.startswith("TER"):
                protein_lines.append(protein_line)

    for cluster_dir_name in cluster_dir_names:
        pdb_path = Path(dir_path) / cluster_dir_name

        for file_entry in pdb_path.iterdir():
            if file_entry.suffix == '.pdb' and not file_entry.stem.endswith('-ref'):
                peptide_lines = []
                with open(file_entry) as f:
                    lines = f.readlines()
                    for line in lines:
                        if line.startswith("ATOM"):
                            peptide_lines.append(line)

                new_pdb_path = pdb_path / f"{file_entry.stem}-ref.pdb"
                with open(new_pdb_path, 'w') as f:
                    f.write(''.join(protein_lines))
                    if not protein_lines[-1].startswith('TER'):
                        f.write("\nTER\n")
                    f.write(''.join(peptide_lines))
                    f.write("\nTER\n")
                    f.write("END\n")

def make_flag(cluster_dir_names, dir_path, protein_chain, peptide_chain):
    for cluster_dir_name in cluster_dir_names:
        work_dir = Path(dir_path) / cluster_dir_name
        os.chdir(work_dir)
        ic(work_dir)

        for file_entry in work_dir.iterdir():
            if file_entry.suffix == '.pdb' and not file_entry.stem.endswith('-ref'):
                ref_pdb = file_entry.stem
                peptide = ref_pdb

                with open('prepack.flags', 'w') as f:
                    f.write(f"-s {peptide}.pdb \n")
                    f.write("-flexpep_prepack \n")
                    f.write("-nstruct 1 \n")
                    f.write("-ex1 \n")
                    f.write("-ex2aro \n")
                    f.write("-mute core.io.database \n")
                    f.write("-mute core.util.prof \n")
                    f.write("-mute protocols.jd2.JobDistributor \n")
                    f.write("-scorefile pack.score.sc ")

                os.system("mpirun -np 30 FlexPepDocking.mpi.linuxgccrelease @prepack.flags")

                with open('refinement.flags', 'w') as f:
                    f.write(f"-s {peptide}_0001.pdb \n")
                    f.write(f"-flexPepDocking:receptor_chain {protein_chain} \n")
                    f.write(f"-flexPepDocking:peptide_chain {peptide_chain} \n")
                    f.write("-pep_refine \n")
                    f.write("-use_input_sc \n")
                    f.write("-nstruct 250 \n")
                    f.write("-scorefile refinement.score.sc \n")
                    f.write("-ex1 \n")
                    f.write("-ex2aro \n")
                    f.write("-mute core.io.database \n")
                    f.write("-mute core.util.prof \n")
                    f.write("-mute protocols.jd2.JobDistributor ")

                os.system("mpirun -np 30 FlexPepDocking.mpi.linuxgccrelease @refinement.flags")

def take_candidate(cluster_dir_names, dir_path, result_path):
    if not Path(result_path).exists():
        Path(result_path).mkdir(parents=True)
    for cluster_dir_name in cluster_dir_names:

        score_file = Path(dir_path) / cluster_dir_name / 'refinement.score.sc'

        with open(score_file, 'r') as f:
            headers = None
            rows = []
            lines = f.readlines()
            for line in lines:
                if line.startswith('SEQUENCE:'):
                    continue
                elif line.startswith('SCORE:'):
                    if headers is None:
                        headers = line.split()[1:]
                    else:
                        rows.append(line.split()[1:])
                else:
                    rows.append(line.split())
            df = pd.DataFrame(rows, columns=headers)
            df['reweighted_sc'] = df['reweighted_sc'].astype(float)
            df = df.sort_values(by='reweighted_sc', ignore_index=True)
            min_row = df.loc[df['reweighted_sc'].argmin()]

        candidate_path = Path(dir_path) / cluster_dir_name / f"{min_row['description']}.pdb"
        shutil.copy(candidate_path, Path(result_path))

    result_score_file = Path(result_path) / 'candidate.txt'
    with open(result_score_file, 'a') as f:
        f.write("name\tscore\n")
        f.write(f"{min_row['description']}\t{min_row['reweighted_sc']}\n")

def process_refine_pdb(result_path):
    for file_entry in Path(result_path).iterdir():
        if file_entry.suffix == '.pdb':
            processed_pdb_path = Path(result_path) / file_entry.name
            candidate_pdb_path = Path(result_path) / f"{file_entry.stem.split('-')[0]}.pdb"
            candidate_pdb = []

            with open(processed_pdb_path) as f:
                lines = f.readlines()
                for line in lines:
                    if line.startswith('ATOM') or line.startswith('TER') or line.startswith('END'):
                        candidate_pdb.append(line)

            with open(candidate_pdb_path, 'w') as f:
                f.write(''.join(candidate_pdb))


if __name__ == '__main__':
    dir_path = '/mnt/sdc/lanwei/MC1R/final_best_pose_pdb/'
    result_path = f'{dir_path}/refinement/'
    protein_pdb = '/mnt/sdc/lanwei/MC1R/MC1R.pdb'
    protein_chain = "A"
    peptide_chain = "E"

    cluster_dir_names, pdb_names = read_pdb(dir_path)
    preprocess_pdb(protein_pdb, cluster_dir_names, dir_path)
    make_flag(cluster_dir_names, dir_path, protein_chain, peptide_chain)
    take_candidate(cluster_dir_names, pdb_names, dir_path, result_path)
    process_refine_pdb(result_path)

    logging.info("Script execution completed.")
