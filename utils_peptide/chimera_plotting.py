"""
This script employs the headless version of Chimera and operates without a Graphical User Interface (GUI) by default. 
The standard version has also been tested, demonstrating the necessity of a GUI when saving figures.
"""
import subprocess
import os
from utils_comm.log_util import logger

# Location of the Chimera binary
chimeraBinaryExecPath = (
    "/mnt/nas/yuliu/softwares/UCSF-Chimera64-1.17.3-headless/bin/chimera"
)


def run_chimera_script(chimera_path, script_path):
    """
    Run a Chimera script using a specified Chimera binary.

    Args:
        chimera_path (str): The path to the Chimera binary.
        script_path (str): The path to the Chimera script.
    """
    # Construct the command
    command = [chimera_path, "--script", script_path, "--nogui"]

    # Run the command
    subprocess.run(command, check=True)


def create_chimera_script(pdb_file, output_file):
    """
    Create a Chimera script that processes a PDB file and saves an image.

    Args:
        pdb_file (str): The path to the PDB file.
        output_file (str): The path to the output image file.

    Returns:
        str: The path to the created Chimera script.
    """
    # Create the script content
    script_content = """
from chimera import runCommand as rc
from chimera import replyobj

# Emit a status message
replyobj.status("Processing {0}")

# Open the pdb file
rc("open {0}")

# Apply the commands
rc("set bgColor white")
rc("lighting mode two-point")
rc("color blue :.a")
rc("color red :.b")

# Save image to a file
rc("copy file {1} jpeg width 2000 height 2000 supersample 3")

# Close all models
rc("close all")

# Print the output file name in case it's needed
print("Output file is: {1}")
    """.format(
        pdb_file, output_file
    )

    # Write the script content to a file
    script_path = "chimera_script.py"
    with open(script_path, "w") as f:
        f.write(script_content)

    return script_path


def chimera_plotting(workdir):
    """
    Runs a Chimera script to process all PDB files in a directory and save images.

    Args:
        workdir (str): The directory containing the PDB files.
    """
    # Specify the path to the Chimera binary
    chimera_bin_path = chimeraBinaryExecPath
    os.chdir(workdir)

    # Create a subdirectory named 'chimera_plotting'
    output_dir = os.path.join(workdir, "chimera_plotting")
    os.makedirs(output_dir, exist_ok=True)

    # Get a list of all .pdb files in the directory
    pdb_files = [f for f in os.listdir(workdir) if f.endswith(".pdb")]
    total_pdb_count = len(pdb_files)
    logger.info(f"Total number of PDB files: {total_pdb_count}")

    # Process each pdb file
    for pdb_file in pdb_files:
        # Create the output file name based on the pdb_file name
        base_name = os.path.splitext(pdb_file)[0]
        output_file = os.path.join(output_dir, f"{base_name}_plotting.jpg")

        # Create the Chimera script
        chimera_script_path = create_chimera_script(pdb_file, output_file)

        # Run the Chimera script
        run_chimera_script(chimera_bin_path, chimera_script_path)

        # Delete the Chimera script
        os.remove(chimera_script_path)


# The main entry point for the script
if __name__ == "__main__":
    logger.info("start")
    path = (
        "/mnt/nas/yuliu/repos/peptide-deploy/utils_peptide/linear_pep_protein_complex"
    )
    # Run the chimera plotting function
    chimera_plotting(path)
    logger.info("finish")
