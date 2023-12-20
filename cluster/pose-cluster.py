import os
import numpy as np

import hdbscan
from pathlib import Path
from hdbscan import validity
from sklearn.manifold import TSNE
from tqdm import tqdm
from joblib import Parallel, delayed
import logging
import time
import random
import math

import multiprocessing as mp
import matplotlib.pyplot as plt
from shutil import copyfile

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
now_times = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

# Read protein PDB and extract coordinates
def read_protein(protein_pdb):
    """
    Read a protein file in PDB format and extract the coordinates of each atom.

    Args:
        protein_pdb (str): The path to the protein PDB file.

    Returns:
        list: A list of lists, where each inner list contains the x, y, and z coordinates of an atom.
    """
    protein_coords = []
    with open(protein_pdb) as f:
        for line in f:
            if line.startswith("ATOM"):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                protein_coords.append([x, y, z])
    logger.info(f"Read protein {len(protein_coords)} atom coordinates from PDB file.")
    return protein_coords


def read_peptide(peptide_pdb_dir):
    """
    Reads the coordinates and names of peptides from PDB files in a given directory.

    Parameters:
        peptide_pdb_dir (str): The directory containing the PDB files of the peptides.

    Returns:
        tuple: A tuple containing two lists. The first list contains the coordinates of each peptide,
               where each coordinate is a list of three float values representing the x, y, and z
               coordinates of an atom. The second list contains the names of the peptides, where each
               name is a string representing the name of the PDB file without the .pdb extension.
    """
    peptide_coords = []
    peptide_names = []
    for pdb_path in Path(peptide_pdb_dir).glob("*.pdb"):
        coords = []
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM"):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
        peptide_coords.append(coords)
        peptide_names.append(pdb_path.stem)
    logger.info(f'{now_times}: read {len(peptide_coords)} peptides')
    return peptide_coords, peptide_names


# 计算多肽和蛋白的接触图谱(权重图谱)，截距设为6埃
def calc_contacts(peptide_idx, peptide_coords, protein_coords, contact_cutoff):
    """
    Calculate the number of contacts between a peptide and a protein.

    Parameters:
    - peptide_idx (int): The index of the peptide in the list of peptide coordinates.
    - peptide_coords (list): A list of lists containing the coordinates of each peptide.
    - protein_coords (list): A list containing the coordinates of the protein.
    - contact_cutoff (float): The distance threshold for considering a contact.

    Returns:
    - dists (list): A list containing the number of contacts for each protein coordinate.
    """
    dists = [0] * len(protein_coords)
    for p_idx, p_coord in enumerate(protein_coords):
        contact = 0
        for pe_coord in peptide_coords[peptide_idx]:
            dist = np.linalg.norm(np.array(p_coord) - np.array(pe_coord))
            if dist < contact_cutoff:
                contact += 1
        dists[p_idx] = contact
    
    return dists


def mutiprocessing(peptide_coords, protein_coords, contact_cutoff):
    """
    Generates contact matrices for each peptide coordinate using multiprocessing.

    Parameters:
        peptide_coords (list): The list of peptide coordinates.
        protein_coords (list): The list of protein coordinates.

    Returns:
        list: A list of contact matrices for each peptide coordinate.
    """
    contact_matrices = []
    with mp.Pool() as pool:
        results = [
            pool.apply_async(
                calc_contacts, args=(i, peptide_coords, protein_coords, contact_cutoff)
            )
            for i in range(len(peptide_coords))
        ]
        contact_matrices = [r.get() for r in results]

    logger.info(f'{now_times} cutoff:{contact_cutoff} contact_matrices length: {len(contact_matrices)}')
    return contact_matrices


# 将接触图谱转为接触矩阵
def caculate_contact_matrix(contact_matrices):
    """
    Generate a contact matrix from a list of contact matrices.

    Parameters:
        contact_matrices (list): A list of contact matrices.

    Returns:
        numpy.ndarray: The contact matrix as a numpy array.
    """
    contact_matrixs = []
    for dists in contact_matrices:
        contact_matrixs.append(dists)
    contact_matrixs = np.array(contact_matrixs)

    logger.info(f"contact_matrices length: {len(contact_matrices)}")
    return contact_matrixs


# 搜索最佳参数
def cluster_contact_matrix(contact_matrixs,DBCV_score_vs_trials):
    """
    Runs HDBSCAN clustering on a contact matrix using various parameter combinations,
    and returns the best model based on the maximum score.

    Parameters:
    - contact_matrixs (list): A list of contact matrices to be clustered.

    Returns:
    - best_model: The best HDBSCAN model based on the maximum score.

    Example Usage:
    contact_matrixs = [matrix_1, matrix_2, matrix_3]
    best_model = cluster_contact_matrix(contact_matrixs)
    """

    scores = {}

    def run_hdbscan(params):
        min_cluster_size, min_samples, eps, cluster_selection_method = params
        hdb = hdbscan.HDBSCAN(
            min_cluster_size=min_cluster_size,
            min_samples=min_samples,
            cluster_selection_epsilon=eps,
            cluster_selection_method=cluster_selection_method,
            gen_min_span_tree=True,
            core_dist_n_jobs=-1,
        ).fit(contact_matrixs)
        score = hdb.relative_validity_
        return hdb, score

    min_cluster_sizes = list(range(2, 5))
    min_samples_values = list(range(2, 5))
    eps_values = [x * 0.1 for x in range(1, 30, 2)]
    cluster_selection_methods = [
        "eom", 
        #"leaf"
        ]

    param_combinations = [
        (min_cluster_size, min_samples, eps, cluster_selection_method)
        for min_cluster_size in min_cluster_sizes
        for min_samples in min_samples_values
        for eps in eps_values
        for cluster_selection_method in cluster_selection_methods
    ]

    results = Parallel(n_jobs=-1)(
        delayed(run_hdbscan)(params) for params in tqdm(param_combinations)
    )

    for hdb, score in results:
        scores[hdb] = score

    plt.plot(list(range(len(scores))), scores.values())
    plt.xlabel("Trials")
    plt.ylabel("DBCV Score")
    plt.title("DBCV Score vs Trials")
    plt.savefig(DBCV_score_vs_trials)

    best_model = max(scores, key=scores.get)

    logger.info(f"{now_times} Best Model: {best_model}, Score: {scores[best_model]}")

    return best_model


def cluster(best_model, contact_matrixs):
    """
    Clusters a set of contact matrices using the best model.

    Parameters:
        best_model (Model): The best model for clustering.
        contact_matrixs (list): A list of contact matrices to cluster.

    Returns:
        list: The labels assigned to each contact matrix after clustering.
    """
    hdb = best_model.fit(contact_matrixs)
    labels = hdb.labels_
    score = hdb.relative_validity_
    logger.info(f'{now_times} Best Model: {best_model}, Score: {score}')
    return labels


def show_cluster(labels, peptide_names, cluster_result):
    """
    Given a list of labels, peptide names, and a cluster result file, this function
    generates a dictionary where each label is mapped to a list of peptide names
    belonging to that label. The function writes the cluster information to the
    cluster result file in the format "Cluster {label}: {peptide_names} {count}".
    Finally, the function returns the dictionary of cluster keys.

    :param labels: A list of labels.
    :type labels: List
    :param peptide_names: A list of peptide names.
    :type peptide_names: List
    :param cluster_result: The path to the cluster result file.
    :type cluster_result: str
    :return: A dictionary where each label is mapped to a list of peptide names.
    :rtype: Dict

    """
    HDB_cluster_dict = {}
    for label in set(labels):
        HDB_cluster_dict[label] = []
    for i, label in enumerate(labels):
        HDB_cluster_dict[label].append(peptide_names[i])
    for label in sorted(HDB_cluster_dict):
        with open(cluster_result, "a+") as f:
            f.write(
                f"Cluster {label}:{HDB_cluster_dict[label]}\n {len(HDB_cluster_dict[label])}\n"
            )
    logger.info(f"{now_times} Cluster information saved to {cluster_result}")
    return HDB_cluster_dict


def draw_2d_tsne(contact_matrixs, labels, HDB_cluster_dict, point_plot):
    """
    Generate a 2D t-SNE plot for the given contact matrix data.

    Parameters:
    - contact_matrixs (numpy array): An array of contact matrix data.
    - labels (numpy array): An array of labels corresponding to each data point.
    - HDB_cluster_dict (list): A list of cluster keys.
    - point_plot (str): The path to save the generated plot.

    Returns:
    None
    """
    tsne = TSNE(n_components=2, perplexity=5, random_state=42)
    X_tsne = tsne.fit_transform(contact_matrixs)
    fig, ax = plt.subplots()
    for i in range(len(HDB_cluster_dict)):
        ax.scatter(
            X_tsne[labels == i, 0],
            X_tsne[labels == i, 1],
        )
    plt.savefig(point_plot)
    logger.info(f"{now_times} 2D t-SNE plot saved to {point_plot}")

# 从聚类中提取随机蛋白，不包含噪音类(label_tag : -1)
def extract_random_peptides(HDB_cluster_dict, peptide_number):
    """
    Extracts a specified number of random peptides from a collection of cluster keys.
    
    Parameters:
        HDB_cluster_dict (dict): A dictionary mapping labels to lists of peptides. The keys represent the labels of the clusters, and the values represent the peptides in each cluster.
        random_peptides (int): The number of random peptides to extract.
        
    Returns:
        dict: A dictionary mapping labels to lists of randomly selected peptides. The keys represent the labels of the clusters, and the values represent the randomly selected peptides from each cluster.
    """
    real_cluster = {}
    for i in range(0,len(HDB_cluster_dict)-1):
        a= HDB_cluster_dict[i]
        real_cluster[i] = a  

    total_peptides = sum(len(peptides) for peptides in real_cluster.values())
    weights = {label: len(peptides) / total_peptides for label, peptides in real_cluster.items()}
    peptides_to_extract = {label: math.ceil(peptide_number * weight) for label, weight in weights.items()}

    selected_peptides = {}
    for label, peptides in real_cluster.items():
        selected_peptides[label] = random.sample(peptides, peptides_to_extract[label])
        
    peptide_number = sum(len(peptides) for peptides in selected_peptides.values())
    logger.info(f"{now_times} {peptide_number} peptides are randomly selected")

    return selected_peptides


def get_score(peptide, peptide_pdb_dir, HPEP, ADCP):
    pdb_file = Path(peptide_pdb_dir) / f"{peptide}.pdb"
    with open(pdb_file) as f:
        for line in f:
            if HPEP:
                score_marker = "REMARK ITScore"
            elif ADCP:
                score_marker = "REMARK SCORE"
            else:
                break
            if line.startswith(score_marker):
                score_str = line.split()[-1]
                try:
                    return float(score_str)
                except ValueError:
                    return float('inf')
    return float('inf')

# 按照多肽的打分提取需要的数量的多肽
def score_cluster(HDB_cluster_dict, peptide_pdb_dir, HPEP, ADCP, peptide_number):
    """
    Given a dictionary of clusters (HDB_cluster_dict), a directory containing peptide PDB files (peptide_pdb_dir),
    HPEP, ADCP, and the number of peptides to select (peptide_number), this function scores and selects peptides from each 
    cluster based on their score and returns a dictionary of selected peptides.

    Parameters:
    - HDB_cluster_dict (dict): A dictionary containing clusters of peptides.
    - peptide_pdb_dir (str): The directory path where the peptide PDB files are located.
    - HPEP (str): The HPEP parameter.
    - ADCP (str): The ADCP parameter.
    - peptide_number (int): The number of peptides to select.

    Returns:
    - selected_peptides (dict): A dictionary containing the selected peptides from each cluster.
    """
    selected_peptides = {}
    total_peptides = sum(len(peptides) for peptides in HDB_cluster_dict.values())
    weights = {label: len(peptides) / total_peptides for label, peptides in HDB_cluster_dict.items()}
    peptides_to_extract = {label: math.ceil(peptide_number * weight) for label, weight in weights.items()}

    for label, peptides in HDB_cluster_dict.items():
        min_score_peptides = []
        peptides.sort(key=lambda peptide: get_score(peptide, peptide_pdb_dir, HPEP, ADCP))
        min_score_peptides = peptides[:peptides_to_extract[label]]

        selected_peptides[label] = min_score_peptides

    count = 0
    for label, peptides in selected_peptides.items():
        count += len(peptides)

    logger.info(f"{now_times} {count} peptides are selected from all cluster")
    
    return selected_peptides

def get_peptide_pdb_dir(selected_peptides, cluster_path):
    """
    Generates the directory structure and copies PDB files for each peptide in the peptide dictionary to the specified cluster path.

    Args:
        peptide_dict (dict): A dictionary containing peptide labels as keys and tuples of peptide names and other information as values.
        cluster_path (str): The path to the cluster directory where the peptide PDB files will be copied.

    Returns:
        None

    Raises:
        None
    """
    cluster_label_pdbs = []
    cluster_laber_tags = []
    for cluster_labels_keys in selected_peptides.keys():
        cluster_laber_tags.append(int(cluster_labels_keys))
        cluster_label_pdbs.append(selected_peptides[cluster_labels_keys])

    for cluster_label_pdb, cluster_laber_tag in zip(cluster_label_pdbs, cluster_laber_tags):
        cluster_dir = Path(cluster_path) / f"cluster_{cluster_laber_tag}"

        if not cluster_dir.exists():
            cluster_dir.mkdir(parents=True, exist_ok=True)
            
        for peptide in cluster_label_pdb:
            original_pdb = Path(peptide_pdb_dir) / f"{peptide}.pdb"
            clustered_pdb = cluster_dir / f"{peptide}.pdb"

            if not clustered_pdb.exists():
                copyfile(original_pdb, clustered_pdb)
        logger.info(f"Clustered PDB files for cluster {cluster_laber_tag} number :{len(cluster_label_pdb)} pdb_name:{cluster_label_pdb} have been saved to {cluster_dir}")

def merge_pdb(HDB_cluster_dict, cluster_path):
    """
    Merge the PDB files associated with the given cluster keys and save the merged PDB file.

    Parameters:
        HDB_cluster_dict (dict): A dictionary containing the cluster keys as keys and the PDB files associated with each key as values.
        cluster_path (str): The path to the directory where the merged PDB file will be saved.

    Returns:
        None
    """
    cluster_label_pdbs = []
    cluster_laber_tags = []
    for cluster_labels_keys in HDB_cluster_dict.keys():
        cluster_laber_tags.append(int(cluster_labels_keys))
        cluster_label_pdbs.append(HDB_cluster_dict[cluster_labels_keys])

    for cluster_label_pdb, cluster_laber_tag in zip(cluster_label_pdbs, cluster_laber_tags):
        out_file = Path(cluster_path) / f"{cluster_laber_tag}.pdb"
        i = 0
        for pdb in cluster_label_pdb:
            pdb_file = Path(peptide_pdb_dir) / f"{pdb}.pdb"
            with open(out_file, "a") as f_out, open(pdb_file) as f_in:
                i = i + 1
                f_out.write(f"\nMODEL {i}\n")
                # Filter out MODEL/END/ENDMDL lines
                lines = [
                    line
                    for line in f_in
                    if "MODEL" not in line
                    and "END" not in line
                    and "ENDMDL" not in line
                ]
                f_out.write("".join(lines))
                f_out.write(f_in.read())
                f_out.write("END\n")
                f_out.write("ENDMDL\n\n")    
    logger.info(f"{now_times} Clustered PDB files saved to {cluster_path}")

def main(protein_pdb, peptide_pdb_dir, cluster_path , random_peptides):
    """
    Generate a plot and result file for cluster analysis.

    Reads protein and peptide coordinates from PDB files.
    Performs multiprocessing to calculate contact matrices.
    Calculates contact matrix for each peptide.
    Clusters the contact matrices using a clustering algorithm.
    Assigns labels to the clusters.
    Writes the cluster results to a text file.
    Scores the clusters based on a given criteria.
    Generates a 2D t-SNE plot of the contact matrices.
    Merges the PDB files for the clusters.

    Parameters:
    None

    Returns:
    None
    """
    protein_pdb = Path(protein_pdb)
    peptide_pdb_dir = Path(peptide_pdb_dir)
    cluster_path = Path(cluster_path)

    os.umask(0)
    os.makedirs(cluster_path, exist_ok=True)

    protein_coords = read_protein(protein_pdb)
    peptide_coords, peptide_names = read_peptide(peptide_pdb_dir)

    contact_matrices = mutiprocessing(peptide_coords, protein_coords, contact_cutoff)
    contact_matrixs = caculate_contact_matrix(contact_matrices)

    DBCV_score_vs_trials = cluster_path / "DBCV_score_vs_trials.png"
    best_model = cluster_contact_matrix(contact_matrixs, DBCV_score_vs_trials)

    cluster_result = cluster_path / "cluster_result.txt"
    labels = cluster(best_model, contact_matrixs)
    HDB_cluster_dict = show_cluster(labels, peptide_names, cluster_result)

    # 现在提取出的多肽不包括噪音类，也就是会少一个类
    if HPEP or ADCP:
        select_peptides = score_cluster(HDB_cluster_dict, peptide_pdb_dir, HPEP, ADCP, peptide_number)
    else:
        select_peptides = extract_random_peptides(HDB_cluster_dict, peptide_number)
    get_peptide_pdb_dir(select_peptides, cluster_path)
    
    point_plot = cluster_path / "cluster_point_plot.png"
    draw_2d_tsne(contact_matrixs, labels, HDB_cluster_dict, point_plot)

    merge_pdb(HDB_cluster_dict, cluster_path)


if __name__ == "__main__":
    protein_pdb = "/mnt/nas1/lanwei-125/FGF5/dock_prepare/FGF5.pdb"
    peptide_pdb_dir = "/mnt/nas1/lanwei-125/FGF5/disulfide_peptide/"
    cluster_path = "/mnt/nas1/lanwei-125/FGF5/disulfide_peptide_cluster/"

    HPEP = False  # if your structures are from HPEP, set it to True, else set it to False
    ADCP = True   # if your structures are from ADCP, set it to True, else set it to False

    peptide_number = 70 #set the number of peptides you want，if HPEP and ADCP is False, it will be ramdomly selected 
    contact_cutoff = 7  # 6 angstroms

    main(protein_pdb, peptide_pdb_dir, cluster_path , peptide_number)
    print("done")
