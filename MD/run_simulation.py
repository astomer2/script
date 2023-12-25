from itertools import cycle
import os
import shutil
from pathlib import Path
import subprocess
import argparse

def extract_three_letter_code(peptide_structure):
    """
    Extracts the three-letter amino acid sequences, sequence lengths, and file names from .pdb files in a given directory.

    :param work_path: Path to the directory containing the .pdb files.
    :return: A tuple of three lists: sequences, sequence_lengths, and file_names.
    """
    sequences = []

    with open(peptide_structure, "r", encoding="utf-8") as f:
        lines = f.readlines()
    sequence = []
    residue_ids = []
    for line in lines:
        if line.startswith("ATOM") :
            res_id = line[22:27].strip()
            if res_id not in residue_ids:
                residue_ids.append(res_id)
                res_name = line[17:20]
                sequence.append(res_name)
    sequence_str = " ".join(sequence)
    sequences.append(sequence_str)
    sequence_lengths = len(sequence)
    return sequences, sequence_lengths


def reference_id(protein_structure):
    """
    Generates a reference ID based on the protein structure.
    Reads the protein structure from the given file path and retrieves the lines.
    Iterates through the lines and extracts the protein ID from each line that starts
    with 'ATOM'. If the protein ID is different from the previous one, increments the
    residue count. Finally, constructs the reference ID string and returns the protein ID
    and reference ID.
    Parameters:
        protein_structure (str): The file path of the protein structure.
    Returns:
        tuple: A tuple containing the protein ID (str) and reference ID (str).
    """
    with open(protein_structure) as f:
        lines = f.readlines()
    protein_id = 0
    residue_count = 0
    prev_res_id = 0
    for line in lines:
        if line.startswith('ATOM'):
            protein_id = line[22:26].strip()
            if protein_id != prev_res_id:
                residue_count += 1
            prev_res_id = protein_id
    reference = f"1-{protein_id}"
    print(reference)
    print(protein_id)
    return protein_id, reference

def target_id(protein_id, peptide_structure):
    """
    Function to generate a target ID based on a protein ID and peptide structure.

    Parameters:
    protein_id (int): The ID of the protein.

    Returns:
    str: The generated target ID.

    This function reads the peptide structure from a file and calculates the target ID based on the protein ID and peptide ID. It iterates over the lines in the file, extracting the peptide ID from lines starting with 'ATOM'. If the peptide ID is different from the previous one, the residue count is incremented. Finally, the target ID is calculated as the sum of the protein ID and the peptide ID, plus 1.

    Example:
    >>> target_id(5)
    '6-11'
    """
    with open(peptide_structure) as f:
        lines = f.readlines()

    residue_count = 0
    prev_res_id = 0
    peptide_id = 0 
    for line in lines:
        if line.startswith('ATOM'):
            peptide_id = line[22:26].strip()
            if peptide_id != prev_res_id:
                residue_count += 1
            prev_res_id = peptide_id

    target = f"{int(protein_id)+1}-{int(peptide_id)+int(protein_id)}"
    print(peptide_id)
    print (target)
    return target

def identify_cys_id(peptide_structure):
    """
    This function takes in the path of a peptide structure file and identifies the residue IDs of all CYS or CYX residues in the file.

    Parameters:
    - peptide_structure (str): The path of the peptide structure file.

    Returns:
    - residue_ids (list): A list of integer residue IDs corresponding to the CYS or CYX residues found in the file.
    """

    cys_residue_ids = []
    new_line = []
    prev_id = 0
    with open (peptide_structure) as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('ATOM'):
            residue = line[17:20]
            if residue == 'CYS' or residue == 'CYX':
                line = line[:17] + 'CYX' + line[20:]
                new_line.append(line)
                residue_id = line[22:26].strip()
                if residue_id != prev_id:
                    prev_id = residue_id
                    cys_residue_ids.append(int(residue_id))
            else :
                new_line.append(line)
    with open(peptide_structure, 'w') as f:
        f.writelines(new_line)
    return cys_residue_ids

def change_cys_to_cyx(peptide_structure):
    """
    Change all occurrences of 'CYS' to 'CYX' in the given peptide structure file.

    Parameters:
    - peptide_structure (str): The path to the peptide structure file.

    Returns:
    - None

    Raises:
    - FileNotFoundError: If the peptide structure file does not exist.
    """
    with open (peptide_structure) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith('ATOM'):
            residue = line[17:20]
            if residue == 'CYS':
                line = line[:17] + 'CYX' + line[20:]
            lines[i] = line

    with open (peptide_structure, 'w') as f:
        f.writelines(lines)

    

def sovlate_pdb(peptide_structure:str, protein_structure:str, induced_hydrogen:bool, solvateions:bool ,pepcyc:bool, pepcys:bool):
    """
    Generate the solvated structure of a peptide-protein complex using tleap.

    Parameters:
    - peptide_structure (str): The path to the input PDB file containing the peptide structure.
    - protein_structure (str): The path to the input PDB file containing the protein structure.
    - induced_hydrogen (bool): Flag indicating whether to add induced hydrogen to the structures.
    - solvateions (bool): Flag indicating whether to solvate the complex with water molecules and ions.
    - pepcyc (bool): Flag indicating whether to perform pepcyc modifications on the peptide.
    - pepcys (bool): Flag indicating whether to perform pepcys modifications on the peptide.

    Returns:
    use amber tleap to add sovl and water in peptide-prorein complex, return complex_solvated.prmtop and complex_solvated.inpcrd 
    """
    Volume = 0
    charge = 0
    ions =  0
    if induced_hydrogen :
        rec_name = protein_structure.split(".")[0]
        pep_name = peptide_structure.split(".")[0]

        os.system("pdb4amber -y -i  " + protein_structure + "  -o " + rec_name + "_noh.pdb")
        os.system("pdb4amber -y -i  "+ peptide_structure + "  -o " + pep_name + "_noh.pdb")
        
        #os.system("pdb4amber --add-missing-atoms -i  " +  rec_name + "_noh.pdb" + "  -o " + rec_name + "_h.pdb")
        #os.system("pdb4amber --add-missing-atoms -i  "+ pep_name + "_noh.pdb" + "  -o " + pep_name + "_h.pdb")

        os.remove(protein_structure)
        os.remove(peptide_structure)

        os.rename( protein_structure.split(".")[0] + "_noh.pdb", protein_structure )
        os.rename(peptide_structure.split(".")[0]+ "_noh.pdb", peptide_structure )


    if solvateions:
        f=open(work_path+"/"+"unsolvateions"+"_tleap.in","w")
        f.write("source leaprc.protein.ff14SB \n")
        f.write("source leaprc.water.tip3p \n")
        f.write("loadamberparams frcmod.ions1lm_1264_tip3p \n")
        f.write("protein= loadpdb " + protein_structure + "\n")
        f.write("peptide= loadpdb " + peptide_structure + "\n")
        f.write("complex = combine {protein peptide} \n")
        f.write("set default PBRadii mbondi2 \n")
        f.write("charge complex \n")
        f.write("solvatebox complex TIP3PBOX 15.0 \n")
        f.write("quit \n")
        f.close()

        os.system("tleap -f "+ work_path+ "/" +"unsolvateions_tleap.in ")

        f=open(work_path+"/"+"leap.log","r")
        ln=f.readlines()
        f.close()

        for x in range(0,len(ln)):
            if ln[x].find("Volume:") !=-1:
                Volume=float(ln[x].split()[1])   
        for x in range(0,len(ln)):
            if ln[x].find("Total perturbed charge:") !=-1:
                charge=float(ln[x].split()[3])   
                print(charge)
    f=open(work_path + "/" + "solvateions"+"_tleap.in","w")
    f.write("source leaprc.protein.ff14SB \n")
    f.write("source leaprc.water.tip3p \n")
    f.write("loadamberparams frcmod.ions1lm_1264_tip3p \n")

    if pepcyc :
        sequences, sequence_lengths = extract_three_letter_code(peptide_structure)
        formatted_sequence = "{" + " ".join(sequences) + "}"
        f.write("peptide= loadpdbusingseq  " + peptide_structure + " "+formatted_sequence +"\n")
        f.write("remove peptide peptide."+ str(sequence_lengths) +".OXT \n")
        f.write("bond peptide.1.N peptide."+ str(sequence_lengths) +".C \n")

    elif  pepcys : 
        cys_residue_ids  = identify_cys_id(peptide_structure)
        f.write("peptide = loadpdb " + peptide_structure + "\n")
        f.write("bond peptide."+ str(cys_residue_ids[0]) +".SG peptide."+ str(cys_residue_ids[1]) +".SG \n")
    else:
        f.write("peptide= loadpdb " + peptide_structure + "\n")

    f.write("protein= loadpdb " + protein_structure + "\n")
    f.write("saveamberparm protein protein.prmtop protein.inpcrd \n")  
    f.write("saveamberparm peptide peptide.prmtop peptide.inpcrd \n")  
    f.write("complex = combine {protein peptide} \n")
    f.write("saveamberparm complex complex.prmtop complex.inpcrd \n")
    f.write("set default PBRadii mbondi2 \n")

     

    if solvateions:

        f.write("solvateBox complex TIP3PBOX 15 \n") 
        if charge<=0:
            f.write("addIons complex Na+ "+str(abs(charge))+"\n") 
        else:
            f.write("addIons complex Cl- "+str(charge)+"\n")

        ions=round(abs(int(conc*6.023*0.0001*Volume))) 
        f.write("addIonsRand complex Na+ "+str(ions)+" Cl- "+ str(ions) +"\n") 
        f.write("saveamberparm complex complex_solvated.prmtop complex_solvated.inpcrd \n")  
        f.write("savepdb complex complex_solvated.pdb  \n")      
    f.write("quit")
    f.close()

    # call tleap to generate the file
    os.system("tleap -f "+work_path + "/"+"solvateions"+"_tleap.in")
    #os.system("rm "+work_path + "/"+"solvateions"+"_tleap.in")
    #os.system("rm "+work_path + "/"+"unsolvateions"+"_tleap.in")

    return

def take_paramter_flie(parmter_file:str):
    """
    Copies all files with the extension '.in' from the specified parameter file directory to the working directory.

    Parameters:
        parmter_file (str): The path to the directory containing the parameter files.

    Returns:
        None
    """
    for filename in os.listdir(parmter_file):
        if filename.endswith('.in'):
            source_path = os.path.join(parmter_file, filename)
            dest_path = os.path.join(work_path, filename)
            shutil.copy(source_path, dest_path)


def add_cmap_to_prmtop(cmap_path):
    """
    Adds cmap to the given prmtop file.

    Args:
        cmap_path (str): The path to the directory containing the cmap files.

    Returns:
        None
    """
    cmap_para_file = f'{cmap_path}/ff14IDPSFF.para'
    #CAMP_scrpt = '/home/weilan/software/ADD-CMAP-master/ADD_CMAP.py'
    prmtop_files = ['peptide.prmtop', 'protein.prmtop', 'complex.prmtop', 'complex_solvated.prmtop']
    for prmtop_file in prmtop_files:
        output_file = f'{prmtop_file.split(".")[0]}_CMAP.prmtop'
        cmd = ['python', 
               f'{cmap_path}/ADD_CMAP.py',
               '-p', prmtop_file,  
               '-c', cmap_para_file,
               '-o', output_file,
               '-s']
        
        result = subprocess.run(cmd)



def run_simulation(cuda_device_id ,CMAP ,cmap_path):
    """
    Runs a simulation using the specified CUDA device ID, CMAP flag, and CMAP path.

    Parameters:
        cuda_device_id (int): The ID of the CUDA device to use for the simulation.
        CMAP (bool): A flag indicating whether to add CMAP to the PRMTOP file.
        cmap_path (str): The path to the CMAP file.

    Returns:
        None
    """

    os.environ["CUDA_VISIBLE_DEVICES"] = str(cuda_device_id)
    if CMAP :
        add_cmap_to_prmtop(cmap_path)
        complex_prmtop = 'complex_solvated_CMAP.prmtop'
    else :
        complex_prmtop = 'complex_solvated.prmtop'
    #maybe write a loop?

    os.system(f"pmemd.cuda -O -i 01.min.in -o min1.out -p {complex_prmtop} -c complex_solvated.inpcrd -r complex_min_01.rst -ref complex_solvated.inpcrd")
    os.system(f"pmemd.cuda -O -i 02.min.in -o min2.out -p {complex_prmtop} -c complex_min_01.rst -r complex_min_02.rst -x complex_min_02.mdcrd -ref complex_min_01.rst")
    os.system(f"pmemd.cuda -O -i 03.min.in -o min3.out -p {complex_prmtop} -c complex_min_02.rst -r complex_min_03.rst -x complex_min_03.mdcrd -ref complex_min_02.rst")
    os.system(f"pmemd.cuda -O -i 04.min.in -o min4.out -p {complex_prmtop} -c complex_min_03.rst -r complex_min_04.rst -x complex_min_04.mdcrd -ref complex_min_03.rst")
    os.system(f"pmemd.cuda -O -i 05.min.in -o min5.out -p {complex_prmtop} -c complex_min_04.rst -r complex_min_05.rst -x complex_min_05.mdcrd -ref complex_min_04.rst")
    os.system(f"pmemd.cuda -O -i 06.min.in -o min6.out -p {complex_prmtop} -c complex_min_05.rst -r complex_min_06.rst -x complex_min_06.mdcrd -ref complex_min_05.rst")
    os.system(f"pmemd.cuda -O -i 07.min.in -o min7.out -p {complex_prmtop} -c complex_min_06.rst -r complex_min_07.rst -x complex_min_07.mdcrd -ref complex_min_06.rst")
    os.system(f"pmemd.cuda -O -i 08.min.in -o min8.out -p {complex_prmtop} -c complex_min_07.rst -r complex_min_08.rst -x complex_min_08.mdcrd -ref complex_min_07.rst")
    os.system(f"pmemd.cuda -O -i 09.md.in -o md1.out -p {complex_prmtop} -c complex_min_08.rst -r complex_md1.rst -x complex_md1.mdcrd -ref complex_min_08.rst")




def analysis(work_path, reference,target):
    """
    Generates a function comment for the given function body in a markdown code block with the correct language syntax.

    Args:
        work_path (str): The path of the work directory.
        reference (str): The reference string.
        target (str): The target string.

    Returns:
        None
    """
    if CMAP :
        complex_prmtop = 'complex_solvated_CMAP.prmtop'
    else :
        complex_prmtop = 'complex_solvated.prmtop'

    f=open(work_path + "/" + "cpptraj.in","w")

    f.write("parm "+ complex_prmtop + " \n")
    f.write("parmstrip :WAT,Na+,Cl- \n")
    f.write("parmbox nobox \n")
    f.write("parmwrite out fixed.prmtop \n")
    f.write("go \n")

    f.write("clear all \n")
    f.write("parm "+ complex_prmtop + " \n")  
    f.write("trajin complex_md1.mdcrd \n")  
    f.write("strip :WAT,Na+,Cl- \n")
    f.write("autoimage \n")
    f.write("check .1,2 skipbadframes silent nobondcheck \n")
    f.write("trajout fixed.mdcrd nobox \n")
    f.write("go \n")

    f.write("clear all \n")
    f.write("parm fixed.prmtop \n")
    f.write("trajin fixed.mdcrd \n")
    f.write("align :"+reference+"@CA first  \n")
    f.write("rms pep :"+target+"@CA nofit out complex_rmsd.xvg \n")
    f.write("trajout out.pdb offset 500 \n")
    f.write("hbond donormask :"+reference+" acceptormask :"+target+" out nhb.dat avgout acghb.dat \n")
    f.write("go \n")
    f.write("quit \n")
    f.close()

    os.system("cpptraj -i cpptraj.in")
    os.system("MMPBSA.py -O -i mmpbsa.in -o energy.dat -do decomp-energy.dat -sp fixed.prmtop -cp complex.prmtop -rp protein.prmtop -lp peptide.prmtop -y fixed.mdcrd")
    os.system("MMPBSA.py --clean")


def simulation(work_path, parmter_file ,cmap_path , cuda_device_id, induced_hydrogen, solvateions, CMAP, pepcyc, pepcys):

    os.chdir(work_path)
    peptide_structure = f'{work_path}/peptide.pdb'
    protein_structure = f'{work_path}/protein.pdb'
    sovlate_pdb(peptide_structure, protein_structure, induced_hydrogen, solvateions, pepcyc ,pepcys)
    take_paramter_flie(parmter_file)
    run_simulation(cuda_device_id, CMAP, cmap_path)
    protein_id ,reference = reference_id(protein_structure)
    target = target_id(protein_id, peptide_structure)
    analysis(work_path, reference, target)


if __name__ == "__main__":
    # 设置默认参数
    default_cuda_device_id = 0
    default_parmter_file = "/mnt/nas1/software/MD/amber_paramter_files"
    default_cmap_path = "/mnt/nas1/software/MD/CMAP_files"
    # 定义命令行参数
    parser = argparse.ArgumentParser(description="Run amber simulation with specified parameters, if no parameters are specified, default parameters will be used")
    parser.add_argument("-i", "--work_path", required=True, help="Path to the work directory, should have protein.pdb and peptide.pdb, must be required")
    parser.add_argument("-g", "--cuda_device_id", nargs='?',type=int, default=default_cuda_device_id, help="CUDA device ID,  default is 0")
    parser.add_argument("-p", "--parmter_file", nargs='?', type=str, default=default_parmter_file, help="Path to the amber simulation parmter files, default is /mnt/nas1/software/MD/amber_paramter_files,  must be required")
    parser.add_argument("-c", "--CMAP_path", nargs='?',type=str, default=default_cmap_path, help="Path to the CMAP master path, default is /mnt/nas1/software/MD/CMAP_files")
    parser.add_argument("-hyd", "--induced_hydrogen", nargs='?',type=int, choices=[0, 1], default=1, help="Use pdb4amber to reduce oridnally hydrogen, then add amber format hydrogens, default is True")
    parser.add_argument("-solv", "--solvateions",nargs='?', type=int, choices=[0, 1], default=1, help="If enable solvations, Add water, counterions, salt ions, and periodic boundary conditions to the structure's topology, default is True")
    parser.add_argument("-CMAP", "--CMAP", nargs='?', type=int, choices=[0, 1], default=0, help="If enable CMAP, induced charmm36 parameters for Amber prmtop file, default is False")
    parser.add_argument("-cyc", "--pepcyc", nargs='?', type=int, choices=[0, 1], default=0, help="If enable pepcyc, Connect any two amino acids at the beginning and end of a polypeptide so that their free amino and carboxyl groups form a peptide bond, must close CMAP, default is False")
    parser.add_argument("-cys", "--pepcys", nargs='?', type=int, choices=[0, 1], default=0, help="If enable pepcys, Connect any two cysteines in a polypeptide whose distance is greater than 2.05 Angstroms to form a disulfide bond between their sulfhydryl groups, support CMAP, default is False")
    # 解析命令行参数
    args = parser.parse_args()

    # 设置参数
    work_path = args.work_path
    cuda_device_id = args.cuda_device_id
    parmter_file = args.parmter_file
    cmap_path = args.CMAP_path  # 修改为 CMAP_path，与命令行参数一致
    induced_hydrogen = bool(args.induced_hydrogen)
    solvateions = bool(args.solvateions)
    CMAP = bool(args.CMAP)
    pepcyc = bool(args.pepcyc)
    pepcys = bool(args.pepcys)

    # 增加额外的判断
    if pepcyc and CMAP:
        raise ValueError("Error: When pepcyc is True, CMAP must be False.")
    if CMAP and not pepcyc:
        # 可以使用默认路径参数，也可以手动传入-c的参数
        cmap_path = args.CMAP_path

    # 设置默认值

    if parmter_file is None:
        parmter_file = default_parmter_file
    if cmap_path is None:
        cmap_path = default_cmap_path
    if cuda_device_id is None:
        cuda_device_id = default_cuda_device_id

    # 盐浓度
    conc = 0.15

    # 调用 simulation 函数
    simulation(work_path, parmter_file, cmap_path, cuda_device_id, induced_hydrogen, solvateions, CMAP, pepcyc, pepcys)