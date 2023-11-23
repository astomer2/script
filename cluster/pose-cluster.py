import os
import numpy as np

import hdbscan
from hdbscan import validity
from sklearn.manifold import TSNE
from tqdm import tqdm
from joblib import Parallel, delayed

import multiprocessing as mp
import matplotlib.pyplot as plt
from shutil import copyfile


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
    print(len(protein_coords))
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
    for pdb in os.listdir(peptide_pdb_dir):
        if pdb.endswith(".pdb"):
            coords = []
            with open(os.path.join(peptide_pdb_dir, pdb)) as f:
                for line in f:
                    if line.startswith("ATOM"):
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append([x, y, z])
            peptide_coords.append(coords)
            peptide_names.append(pdb.split(".")[0])
    print(len(peptide_coords))
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
                # break

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
    print(len(contact_matrices))

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

    # Define the parameter ranges
    min_cluster_sizes = list(range(2, 5))
    min_samples_values = list(range(2, 5))
    eps_values = [x * 0.1 for x in range(1, 10, 2)]
    cluster_selection_methods = ["eom", "leaf"]

    # Generate all combinations of parameters
    param_combinations = [
        (min_cluster_size, min_samples, eps, cluster_selection_method)
        for min_cluster_size in min_cluster_sizes
        for min_samples in min_samples_values
        for eps in eps_values
        for cluster_selection_method in cluster_selection_methods
    ]

    # Run in parallel and track progress with tqdm
    results = Parallel(n_jobs=-1)(
        delayed(run_hdbscan)(params) for params in tqdm(param_combinations)
    )

    # Store results in the scores dictionary
    for hdb, score in results:
        scores[hdb] = score

    # Plotting
    plt.plot(list(range(len(scores))), scores.values())
    plt.xlabel("Trials")
    plt.ylabel("DBCV Score")
    plt.title("DBCV Score vs Trials")
    plt.savefig(DBCV_score_vs_trials)

    # Get the best model based on the maximum score
    best_model = max(scores, key=scores.get)

    print(f"Best Model: {best_model}, Score: {scores[best_model]}")

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
    print(best_model,score)

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
    cluster_keys = {}
    for label in set(labels):
        cluster_keys[label] = []
    for i, label in enumerate(labels):
        cluster_keys[label].append(peptide_names[i])
    for label in sorted(cluster_keys):
        with open(cluster_result, "a+") as f:
            f.write(
                f"Cluster {label}:{cluster_keys[label]}\n {len(cluster_keys[label])}\n"
            )

    return cluster_keys


def draw_2d_tsne(contact_matrixs, labels, cluster_keys, point_plot):
    """
    Generate a 2D t-SNE plot for the given contact matrix data.

    Parameters:
    - contact_matrixs (numpy array): An array of contact matrix data.
    - labels (numpy array): An array of labels corresponding to each data point.
    - cluster_keys (list): A list of cluster keys.
    - point_plot (str): The path to save the generated plot.

    Returns:
    None
    """
    tsne = TSNE(n_components=2, perplexity=5, random_state=42)
    X_tsne = tsne.fit_transform(contact_matrixs)
    # 2D散点图
    fig, ax = plt.subplots()

    for i in range(len(cluster_keys)):
        ax.scatter(
            X_tsne[labels == i, 0],
            X_tsne[labels == i, 1],
        )
    plt.savefig(point_plot)


# 将聚类中的最优多肽提取出来
def score_cluster(cluster_keys, peptide_pdb_dir, cluster_path):
    """
    Calculate the minimum score and corresponding peptide for each cluster label.

    Parameters:
    - cluster_keys (dict): A dictionary containing cluster labels as keys and a list of peptides as values.
    - cluster_path (str): The path to the directory where the cluster directories will be created.

    Returns:
    - min_scores (dict): A dictionary containing cluster labels as keys and a tuple of the minimum peptide and its corresponding score as values.

    Raises:
    - FileNotFoundError: If the PDB file for a peptide does not exist.
    """
    min_scores = {}
    for label, peptides in cluster_keys.items():
        min_score = float("inf")
        min_peptide = None
        for peptide in peptides:
            pdb_file = os.path.join(peptide_pdb_dir, peptide + ".pdb")
            with open(pdb_file) as f:
                for line in f:
                    if HPEP :
                        score_marker = "REMARK ITScore"
                    else :
                        score_marker = "REMARK SCORE"
                    if line.startswith(score_marker):
                        score = float(line.split()[-1])
                        if score < min_score:
                            min_score = score
                            min_peptide = peptide
        min_scores[label] = (min_peptide, min_score)

    dir_list = []

    for label, (peptide, score) in min_scores.items():

        cluster_dir = os.path.join(cluster_path, "cluster_{}".format(label))
        dir_list.append(cluster_dir)

        if not os.path.exists(cluster_dir):
            os.makedirs(cluster_dir, exist_ok=True)

        src_pdb = os.path.join(peptide_pdb_dir, peptide + ".pdb")
        dst_pdb = os.path.join(cluster_dir, peptide + ".pdb")

        if not os.path.exists(dst_pdb):
            copyfile(src_pdb, dst_pdb)


def merge_pdb(cluster_keys, cluster_path):
    """
    Merge the PDB files associated with the given cluster keys and save the merged PDB file.

    Parameters:
        cluster_keys (dict): A dictionary containing the cluster keys as keys and the PDB files associated with each key as values.
        cluster_path (str): The path to the directory where the merged PDB file will be saved.

    Returns:
        None
    """
    cluster_labels_pdb = []
    for cluster_labels_keys in cluster_keys.keys():
        cluster_labels_pdb.append(cluster_keys[cluster_labels_keys])

    a = 0
    for cluster_label_pdb in cluster_labels_pdb:
        i = 0
        a = a + 1
        # 定义输出pdb文件名
        out_file = os.path.join(cluster_path, f"{a}.pdb")

        for pdb in cluster_label_pdb:
            pdb_file = os.path.join(peptide_pdb_dir, pdb + ".pdb")
            with open(out_file, "a") as f_out:
                with open(pdb_file) as f_in:
                    i = i + 1
                    f_out.write(f"\nMODEL {i}\n")
                    # 过滤掉MODEL/END/ENDMDL行
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


def main():
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
    DBCV_score_vs_trials = f"{cluster_path}/DBCV_score_vs_trials.png"
    point_plot = f"{cluster_path}/cluster_point_plot.png"
    cluster_result = f"{cluster_path}/cluster_result.txt"

    protein_coords = read_protein(protein_pdb)
    peptide_coords, peptide_names = read_peptide(peptide_pdb_dir)

    contact_matrices = mutiprocessing(peptide_coords, protein_coords, contact_cutoff)
    contact_matrixs = caculate_contact_matrix(contact_matrices)

    best_model = cluster_contact_matrix(contact_matrixs, DBCV_score_vs_trials)
    labels = cluster(best_model, contact_matrixs)
    cluster_keys = show_cluster(labels, peptide_names, cluster_result)

    score_cluster(cluster_keys, peptide_pdb_dir, cluster_path)
    draw_2d_tsne(contact_matrixs, labels, cluster_keys, point_plot)
    merge_pdb(cluster_keys, cluster_path)


if __name__ == "__main__":
    protein_pdb = "/mnt/nas1/lanwei-125/PRLR/PRLR.pdb"
    peptide_pdb_dir = "/mnt/nas1/lanwei-125/PRLR/HPEP/hpep-result/"
    cluster_path = "/mnt/nas1/lanwei-125/PRLR/HPEP/cluster/"

    HPEP = True  # if your structures are from HPEP, set it to True, else set it to False
    contact_cutoff = 7  # 6 angstroms

    main()
    print("done")
