import os 
import prody
from prody.measure.contacts import findNeighbors


def identify_interface_residues(workdir, needs_b=0):
    '''This prorams expect chain ID for receptor as 'A' 
    and for peptide as 'B'
    '''
    rec = []
    pep = []
    protein = os.path.join(workdir, "protein.pdb")
    with open(protein,'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("ATOM"):
                rec.append(line)
                

    peptide = os.path.join(workdir, "peptide.pdb")
    with open(peptide,'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("ATOM"):
                pep.append(line)    

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


print(identify_interface_residues('/mnt/sdc/lanwei/TLR2/IFKKITGKLKKWIK'))
    # import pdb; pdb.set_trace()

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
    
    
    ignore_list = identify_interface_residues(pdb_str) # No restrained residues
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
    PDBFile.writeFile(simulation.topology, positions, open(pdb_str[:-4]+"_min.pdb", 'w'))
   
    return ret,system