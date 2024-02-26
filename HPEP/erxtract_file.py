import csv
import os
import gzip
import shutil
import multiprocessing
import statistics
from tqdm import tqdm
import argparse

def extract_tar_gz_files(file_args):
    """
    Extracts tar.gz files from a folder to an output folder.

    Args:
        file_args (tuple): A tuple containing file path and output folder path.

    Returns:
        None
    """
    file_path, output_folder = file_args

    target_folder = os.path.join(
        output_folder, os.path.splitext(os.path.basename(file_path))[0]
    )
    os.makedirs(target_folder, exist_ok=True)

    with gzip.open(file_path, "rb") as f_in:
        shutil.unpack_archive(
            filename=file_path, extract_dir=target_folder, format="gztar"
        )

def multiple_extraction(folder_path, output_folder):
    """
    Extracts multiple tar.gz files from a folder to an output folder using multiple processes.

    Args:
        folder_path (str): The path to the folder containing the tar.gz files.
        output_folder (str): The path to the output folder.

    Returns:
        None
    """

    pool = multiprocessing.Pool()

    all_files = os.listdir(folder_path)
    args_list = [
        (os.path.join(folder_path, file_name), output_folder)
        for file_name in all_files
        if file_name.endswith(".gz")
    ]

    with multiprocessing.Pool() as pool:
        for _ in tqdm(
            pool.imap_unordered(extract_tar_gz_files, args_list), total=len(args_list)
        ):
            pass

def analyze_scores(folder_path):
    peptide_scores = {}
    result_path = os.path.join(folder_path, "results.csv")
    for root, dirs, files in os.walk(folder_path):
        if "hpepdock_all.out" in files:
            with open(os.path.join(root, "hpepdock_all.out")) as f:
                for line in f:
                    data = line.split()
                    peptide = data[4]
                    score = float(data[3])

                    if peptide not in peptide_scores:
                        peptide_scores[peptide] = []

                    peptide_scores[peptide].append(score)
  
        with open(result_path, "w+", newline='') as csvfile:
            fieldnames = ['sequence', 'mix', 'max', 'avg', 'med', 'var', 'machine_score']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()

            for peptide, scores in peptide_scores.items():
                min_score = min(scores)
                max_score = max(scores)
                avg_score = round(statistics.mean(scores), 3)
                median_score = round(statistics.median(scores), 3)
                variance = round(statistics.variance(scores), 3)
                machine_score = round(avg_score + 2.5 * (float(variance) / float(len(scores))) ** 0.5,3)

                writer.writerow({
                    'sequence': peptide,
                    'mix': min_score,
                    'max': max_score,
                    'avg': avg_score,
                    'med': median_score,
                    'var': variance,
                    'machine_score': machine_score
                })

if __name__ == "__main__":
    folder_path = '/mnt/nas1/lanwei-125/MC5R/dock/HPEP/cpep_sequence/'
    output_folder = '/mnt/nas1/lanwei-125/MC5R/dock/HPEP/cpep_sequence/'
    multiple_extraction(folder_path, output_folder)
    analyze_scores(folder_path)