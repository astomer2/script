import os
import gzip
import shutil
import multiprocessing
from tqdm import tqdm


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
        # print(f"Extracted file {file_path} to {target_folder}")


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


if __name__ == "__main__":
    folder_path = "/mnt/nas1/lanwei-125/PRLR/HPEP/HPEP_dock/"
    output_folder = "/mnt/nas1/lanwei-125/PRLR/HPEP/HPEP_dock/"
    multiple_extraction(folder_path, output_folder)
