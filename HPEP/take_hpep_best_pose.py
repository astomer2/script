import os

def parse_hpepdock_all_out(file_path):
    """
    Parses the specified file containing HPEP-Dock output and returns a dictionary
    containing information about peptides.

    Parameters:
    - file_path (str): The path to the file to be parsed.

    Returns:
    - peptide_info (dict): A dictionary where the keys are peptides and the values
      are lists of tuples. Each tuple contains the cluster, number, and score of a 
      peptide.

    Example:
    >>> parse_hpepdock_all_out('/path/to/file.txt')
    {'PEPTIDE1': [('CLUSTER1', 'NUMBER1', 0.5), ('CLUSTER2', 'NUMBER2', 0.7)],
     'PEPTIDE2': [('CLUSTER3', 'NUMBER3', 0.3)]}
    """
    peptide_info = {}
    with open(file_path) as f:
        for line in f:
            data = line.split()
            peptide = data[4]
            cluster = data[0]
            number = data[1]
            score = float(data[3])
            if peptide not in peptide_info:
                peptide_info[peptide] = []
            peptide_info[peptide].append((cluster, number, score))
    return peptide_info

def extract_best_model(pdb_file_path, cluster, number, input_folder):
    """
    Extracts the best model from a PDB file based on the given cluster and number.
    
    Parameters:
        pdb_file_path (str): The path to the output PDB file.
        cluster (str): The cluster to extract the model from.
        number (str): The number of the model to extract.
        input_folder (str): The path to the input folder containing the PDB file.
    
    Returns:
        None
    """
    pdb_data = []
    pdb_section = False
    pdb_file_path_full = os.path.join(input_folder, "hpepdock_all.pdb")
    
    if os.path.exists(pdb_file_path_full):
        with open(pdb_file_path_full) as pdb_file:
            for line in pdb_file:
                if line.startswith("REMARK Cluster:"):
                    _, file_cluster, file_number = line.split()[-3:]
                    if cluster == file_cluster and number == file_number:
                        pdb_section = True
                    else:
                        pdb_section = False
                elif pdb_section:
                    pdb_data.append(line)

        with open(pdb_file_path, "w") as output_pdb:
            output_pdb.writelines(pdb_data)

def process_directory(input_folder, result_dir):
    """
    Process the specified input folder and generate the results in the specified result directory.

    Args:
        input_folder (str): The path to the input folder.
        result_dir (str): The path to the result directory.

    Returns:
        None
    """
    os.makedirs(result_dir, exist_ok=True)

    for root, _, files in os.walk(input_folder):
        if "hpepdock_all.out" in files:
            peptide_info = parse_hpepdock_all_out(os.path.join(root, "hpepdock_all.out"))

            for peptide, info in peptide_info.items():
                cluster, number, _ = info[0]
                pdb_file_path = os.path.join(result_dir, f"{peptide}.pdb")
                extract_best_model(pdb_file_path, cluster, number, root)

if __name__ == "__main__":
    input_folder = "/mnt/nas1/lanwei-125/PRLR/HPEP/HPEP_dock/"
    result_dir = "/mnt/nas1/lanwei-125/PRLR/HPEP/hpep-result/"

    process_directory(input_folder, result_dir)
    print("PDB文件提取完成")
