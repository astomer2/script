import os 
import prody
import shutil
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from prody import *
from prody.measure.contacts import findNeighbors

def copy_and_create_directory(pdb_path):
    """
    将PDB文件复制到工作目录，并创建工作目录
    :param pdb_path: PDB文件路径
    """
    if pdb_path.endswith(".pdb"):
        file_path, file_name = os.path.split(pdb_path)
        file_name, _ = os.path.splitext(file_name)  # 使用os.path.splitext获取文件名和扩展名
        workdir = os.path.join(file_path, file_name)
        os.makedirs(workdir, exist_ok=True)
        shutil.copy(pdb_path, workdir)

        return workdir


def read_pdb_chain(workdir):
    """
    读取PDB文件中的链信息
    :param workdir: 工作目录
    :return : chains 字典
    """
    chains = {}
    for pdb_name in [f for f in os.listdir(workdir) if f.endswith(".pdb")]:
        input_pdb = os.path.join(workdir, pdb_name)

        with open(input_pdb) as f:
            for line in f:
                if line.startswith("ATOM"):
                    chain_id = line[21]
                    chains.setdefault(chain_id, []).append(line)
    return chains

def write_pdb_chain(chains, workdir):
    """
    将PDB文件中的链信息写入到PDB文件中
    :param chains: 链信息
    :param workdir: 工作目录
    """

    protein = os.path.join(workdir, "protein.pdb")
    peptide = os.path.join(workdir, "peptide.pdb")
    chain_ids = len(chains)
    if chain_ids >= 2:
        all_chains = list(chains.values())

        # 设置最后一个值为peptide，其余为protein
        peptide_lines = all_chains[-1]
        protein_lines = all_chains[:-1]

        with open(peptide, "w") as f:
            f.write("".join(peptide_lines))

        with open(protein, "w") as f:
            for lines in protein_lines:
                f.write("".join(lines))
                f.write("TER\n")  # 添加TER行，表示链的结束
    else:
        pass
    return protein, peptide

def identify_interface_residues(workdir, needs_b=0):
    '''This prorams expect chain ID for receptor as 'A' 
    and for peptide as 'B'
    '''
    protein = os.path.join(workdir, "protein.pdb")
    peptide = os.path.join(workdir, "peptide.pdb")
    
    rec = parsePDB(protein)
    pep = parsePDB(peptide)

    near_n = findNeighbors(rec, 5 , pep)

    interacting_chainA_res = []
    if needs_b ==1:
        interacting_chainB_res = []
    for a1,a2,d in near_n:
        
        interacting_chainA_res.append(a1.getResindex())
        if needs_b ==1:
            interacting_chainB_res.append(a2.getResindex())        
        
    interacting_chainA_res = list(set(interacting_chainA_res))
    
    interacting_chainA_res.sort()

    return interacting_chainA_res



def _openmm_minimize( pdb_str: str,env='implicit'):
    """Minimize energy via openmm.
    Adopted from AF2 minimization code """    
    
    pdb = PDBFile(pdb_str)
    
    if env == 'implicit':
        force_field = ForceField("amber99sb.xml",'implicit/gbn2.xml')
        print("Using GBSA (gbn2) environment for energy minimization")
    else:
        force_field = ForceField("amber99sb.xml")
        print("Using in-vacuo environment for energy minimization")
    # integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
    constraints = HBonds
    
    system = force_field.createSystem(  pdb.topology, nonbondedCutoff=1*nanometer, 
                                      constraints=constraints)
    
    
    ignore_list = identify_interface_residues(workdir) # No restrained residues
    restrain_(system, pdb, ignore_list=ignore_list)


    # platform = openmm.Platform.getPlatformByName("CUDA" if use_gpu else "CPU")
    simulation = Simulation( pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
      
    ret = {}
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    ret["einit"] = state.getPotentialEnergy()
    # ret["posinit"] = state.getPositions(asNumpy=True)
    simulation.minimizeEnergy(maxIterations=100, tolerance=0.01)
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    ret["efinal"] = state.getPotentialEnergy()
    # ret["pos"] = state.getPositions(asNumpy=True)
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open(minimized_pdb, 'w'))
   
    return ret,system


def restrain_(system, pdb, chain='A',ignore_list=[]):
    '''To apply strong harmonic restrains on non-interface residues'''
    restraint = CustomExternalForce("0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")    
    restraint.addGlobalParameter('k', 20000.0*kilojoules_per_mole/nanometer*nanometer)
    
    restraint.addPerParticleParameter('x0')
    restraint.addPerParticleParameter('y0')
    restraint.addPerParticleParameter('z0')
    # print(ignore_list)
    res_list_openmm=[]
    for atom in pdb.topology.atoms():
        if len(ignore_list) > 0:
            if ignore_list.count(atom.residue.index) > 0:
                if atom.name == 'CA':
                    res_list_openmm.append(atom.residue.name)
                continue
        
        # import pdb ; pdb.set_trace()
        if chain =='A':
            if atom.element.name == 'hydrogen':
                continue
            if atom.residue.chain.id == 'A':       
                restraint.addParticle(atom.index, pdb.positions[atom.index])
        else:
            restraint.addParticle(atom.index, pdb.positions[atom.index])
    if chain =='A':       
        print("omm ignore  :",res_list_openmm)            
    system.addForce(restraint)

def get_single_point_energy(pdbfile,mode='vacuum'):
    'To get a single point energy of a pdb file'
    pdb_handle = PDBFile(pdbfile)

    if mode == 'implicit':
        force_field = ForceField("amber99sb.xml",'implicit/gbn2.xml')
        print("Using GBSA (gbn2) environment for energy calculation")
    
    else:
        force_field = ForceField("amber99sb.xml")
        print("Using In-Vacuo environment for energy calculation")
        
    system = force_field.createSystem(pdb_handle.topology, nonbondedCutoff=1*nanometer, constraints=HBonds)
    restrain_(system, pdb_handle, 'All')
    integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
    simulation = Simulation(pdb_handle.topology, system, integrator)
    simulation.context.setPositions(pdb_handle.positions)
    simulation.minimizeEnergy(maxIterations=1)
    state = simulation.context.getState(getEnergy=True)
    # import pdb; pdb.set_trace()
        
    return state.getPotentialEnergy()._value

if __name__ == "__main__":
    env='implicit'
    pdb_path = '/mnt/sdc/lanwei/TLR2/IFKKITGKLKKWIK.pdb'
    workdir = copy_and_create_directory(pdb_path)
    chains =read_pdb_chain(workdir)
    protein, peptide = write_pdb_chain( chains,workdir)
    minimized_pdb = pdb_path[:-4]+"_min.pdb"
    _openmm_minimize(pdb_path, env )
    a = get_single_point_energy(pdb_path)
    a1 = get_single_point_energy(protein)
    a2 = get_single_point_energy(peptide)
    print(a,a1,a2)san
    
    workdir1 = copy_and_create_directory(minimized_pdb)
    chains1 =read_pdb_chain(workdir1)
    protein1, peptide1 = write_pdb_chain( chains1,workdir1)
    b = get_single_point_energy(minimized_pdb)
    b1 = get_single_point_energy(protein1)
    b2 = get_single_point_energy(peptide1)
    print(b,b1,b2)