import os
import shutil
import subprocess
import argparse

def sovlate_pdb():
    if induced_hydrogen :
        os.system("pdb4amber -y -i  " + protein_structure + "  -o " + protein_structure.split(".")[0] + "_h.pdb")
        os.system("pdb4amber -y -i  "+ peptide_structure + "  -o " + peptide_structure.split(".")[0] + "_h.pdb")
        os.remove(protein_structure)
        os.remove(peptide_structure)

        os.rename( protein_structure.split(".")[0] + "_h.pdb", protein_structure )
        os.rename(peptide_structure.split(".")[0]+ "_h.pdb", peptide_structure )


    if solvateions:
        f=open(work_path+"/"+"unsolvateions"+"_tleap.in","w")
        f.write("source leaprc.protein.ff14SB \n")
        f.write("source leaprc.water.tip3p \n")
        f.write("loadamberparams frcmod.ions1lm_1264_tip3p \n")
        f.write("protein= loadpdb " + protein_structure + "\n")
        f.write("peptide= loadpdb " + peptide_structure + "\n")
        f.write("complex = combine {protein peptide} \n")
        f.write("set default PBRadii mbondi2 \n")
        #f.write("saveamberparm protein protein.prmtop protein.inpcrd \n")  
        #f.write("saveamberparm peptide peptide.prmtop peptide.inpcrd \n")  
        #f.write("saveamberparm complex complex.prmtop complex.inpcrd \n")  
        f.write("charge complex \n")
        f.write("solvatebox complex TIP3PBOX 15.0 \n")
        f.write("quit \n")
        f.close()
        os.system("tleap -f "+ work_path+ "/" +"unsolvateions_tleap.in ")
        f=open(work_path+"/"+"leap.log","r")
        ln=f.readlines()
        #print(ln)  
        f.close()

        for x in range(0,len(ln)):
            if ln[x].find("Volume:") !=-1:
                Volume=float(ln[x].split()[1])   
                print(Volume)
        for x in range(0,len(ln)):
            if ln[x].find("Total perturbed charge:") !=-1:
                charge=float(ln[x].split()[3])   
                print(charge)

    f=open(work_path + "/" + "solvateions"+"_tleap.in","w")
    f.write("source leaprc.protein.ff14SB \n")
    f.write("source leaprc.water.tip3p \n")
    f.write("loadamberparams frcmod.ions1lm_1264_tip3p \n")
    f.write("protein= loadpdb " + protein_structure + "\n")
    f.write("peptide= loadpdb " + peptide_structure + "\n")
    f.write("saveamberparm protein protein.prmtop protein.inpcrd \n")  
    f.write("saveamberparm peptide peptide.prmtop peptide.inpcrd \n")  
    f.write("complex = combine {protein peptide} \n")
    f.write("saveamberparm complex complex.prmtop complex.inpcrd \n")
    f.write("set default PBRadii mbondi2 \n")

    if solvateions:
        f.write("solvateBox complex TIP3PBOX 15 \n") 
        ions=round(abs(int(conc*6.023*0.0001*Volume)))  
        if charge<=0:
            f.write("addIons complex Na+ "+str(abs(charge))+"\n") 
        else:
            f.write("addIons complex Cl- "+str(charge)+"\n")

    f.write("addIonsRand complex Na+ "+str(ions)+" Cl- "+ str(ions) +"\n")           
    f.write("saveamberparm complex complex_solvated.prmtop complex_solvated.inpcrd \n")  
    f.write("savepdb complex complex_solvated.pdb  \n")
    f.write("quit")
    f.close()
    # call tleap to generate the file
    os.system("tleap -f "+work_path + "/"+"solvateions"+"_tleap.in")
    os.system("rm "+work_path + "/"+"solvateions"+"_tleap.in")
    os.system("rm "+work_path + "/"+"unsolvateions"+"_tleap.in")
    return

def take_paramter_flie(parmter_file):
    for filename in os.listdir(parmter_file):
        if filename.endswith('.in'):
            source_path = os.path.join(parmter_file, filename)
            dest_path = os.path.join(work_path, filename)
            shutil.copy(source_path, dest_path)


def add_cmap_to_prmtop(cmap_path):
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



def run_simulation():
    os.environ["CUDA_VISIBLE_DEVICES"] = str(cuda_device_id)
    if CMAP :
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

def reference_id():
    with open(protein_structure) as f:
        lines = f.readlines()

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
    return protein_id,reference

def target_id(protein_id):
    with open(peptide_structure) as f:
        lines = f.readlines()

    residue_count = 0
    prev_res_id = 0
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


def analysis(reference,target):
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


def simulation():

    os.chdir(work_path)
    sovlate_pdb()
    take_paramter_flie(parmter_file)
    add_cmap_to_prmtop(cmap_path)
    run_simulation()
    protein_id,reference = reference_id()
    target = target_id(protein_id)
    analysis(reference,target)


if __name__ == "__main__":
    # 判断
    induced_hydrogen = True
    solvateions = True
    CMAP = True

    parser = argparse.ArgumentParser(description="Run amber simulation with specified parameters")
    parser.add_argument("-i", "--work_path", required=True, help="Path to the work directory")
    parser.add_argument("-g", "--cuda_device_id", type=int, required=True, help="CUDA device ID")
    parser.add_argument("-p", "--parmter_file", required=True, help="Path to the parmter file")
    parser.add_argument("-c", "--cmap_path", required=True, help="Path to the CMAP master path")
    # Parse the command-line arguments
    
    args = parser.parse_args()
    work_path = args.work_path
    cuda_device_id = args.cuda_device_id
    cmap_path = args.cmap_path
    parmter_file = args.parmter_file

    peptide_structure = f'{work_path}/peptide.pdb'
    protein_structure = f'{work_path}/protein.pdb'
    conc = 0.15
    simulation()


    