import os
import shutil
def sovlate_pdb():

    if solvateions:
        f=open(work_path+"/"+"unsolvateions"+"_tleap.in","w")
        f.write("source leaprc.protein.ff14SB \n")
        f.write("source leaprc.water.tip3p \n")
        f.write("loadamberparams frcmod.ions1lm_1264_tip3p \n")
        f.write("protein= loadpdb " + protein_structure + "\n")
        f.write("peptide= loadpdb " + peptide_structure + "\n")
        f.write("complex = combine {protein peptide} \n")
        f.write("set default PBRadii mbondi2 \n")
        f.write("saveamberparm protein protein.prmtop protein.inpcrd \n")  
        f.write("saveamberparm peptide peptide.prmtop peptide.inpcrd \n")  
        f.write("saveamberparm complex complex.prmtop complex.inpcrd \n")  
        f.write("charge complex \n")
        f.write("solvatebox complex TIP3PBOX 10.0 \n")
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
    f.write("complex = combine {protein peptide} \n")
    f.write("set default PBRadii mbondi2 \n")
    if solvateions:
        f.write("solvateBox complex TIP3PBOX 10 \n") 
        ions=round(abs(int(conc*6.023*0.0001*Volume)))  
        if charge<=0:
            f.write("addIons complex Na+ "+str(charge)+"\n") 
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

def take_paramter_flie():
    for filename in os.listdir(parmter_file):
        if filename.endswith('.in'):
            source_path = os.path.join(parmter_file, filename)
            dest_path = os.path.join(work_path, filename)
            shutil.copy(source_path, dest_path)
    
def run_simulation():
    os.chdir(work_path)
    #maybe write a loop?
    os.environ["CUDA_VISIBLE_DEVICES"] = str(cuda_device_id)
    os.system("pmemd.cuda -O -i 01.min.in -o min1.out -p complex_solvated.prmtop -c complex_solvated.inpcrd -r complex_min_01.rst -ref complex_solvated.inpcrd")
    os.system("pmemd.cuda -O -i 02.min.in -o min2.out -p complex_solvated.prmtop -c complex_min_01.rst -r complex_min_02.rst -x complex_min_02.mdcrd -ref complex_min_01.rst")
    os.system("pmemd.cuda -O -i 03.min.in -o min3.out -p complex_solvated.prmtop -c complex_min_02.rst -r complex_min_03.rst -x complex_min_03.mdcrd -ref complex_min_02.rst")
    os.system("pmemd.cuda -O -i 04.min.in -o min4.out -p complex_solvated.prmtop -c complex_min_03.rst -r complex_min_04.rst -x complex_min_04.mdcrd -ref complex_min_03.rst")
    os.system("pmemd.cuda -O -i 05.min.in -o min5.out -p complex_solvated.prmtop -c complex_min_04.rst -r complex_min_05.rst -x complex_min_05.mdcrd -ref complex_min_04.rst")
    os.system("pmemd.cuda -O -i 06.min.in -o min6.out -p complex_solvated.prmtop -c complex_min_05.rst -r complex_min_06.rst -x complex_min_06.mdcrd -ref complex_min_05.rst")
    os.system("pmemd.cuda -O -i 07.min.in -o min7.out -p complex_solvated.prmtop -c complex_min_06.rst -r complex_min_07.rst -x complex_min_07.mdcrd -ref complex_min_06.rst")
    os.system("pmemd.cuda -O -i 08.min.in -o min8.out -p complex_solvated.prmtop -c complex_min_07.rst -r complex_min_08.rst -x complex_min_08.mdcrd -ref complex_min_07.rst")
    os.system("pmemd.cuda -O -i 09.md.in -o md1.out -p complex_solvated.prmtop -c complex_min_08.rst -r complex_md1.rst -x complex_md1.mdcrd -ref complex_min_08.rst")
    os.system("cpptraj -i cpptraj.in")
    os.system("MMPBSA.py -O -i mmpbsa.in -o energy.dat -do decomp-energy.dat -sp fixed.prmtop -cp complex.prmtop -rp protein.prmtop -lp peptide.prmtop")
    os.system("MMPBSA.py --clean")
if __name__ == "__main__":
    # define the input parameters
    work_path = '/mnt/nas1/lanwei-125/test'
    peptide_structure = f'{work_path}/CRS1.pdb'
    protein_structure = f'{work_path}/IL8.pdb'
    parmter_file = '/mnt/sdc/lanwei/amber_input'
    cuda_device_id = 2
    solvateions = True
    conc = 0.15
    sovlate_pdb()
    take_paramter_flie()
    run_simulation()