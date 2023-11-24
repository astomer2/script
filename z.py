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
'''    res_list_prody=[]
    for i in interacting_chainA_res:
        res_list_prody.append(rec.select(' name CA and resindex %d' % i ).getResnames()[0])
    print("prody ignore:",res_list_prody)
    
    if needs_b ==1:
        interacting_chainB_res = list(set(interacting_chainB_res))
        interacting_chainB_res.sort()
        pep_list_prody=[]
        for i in interacting_chainB_res:
            pep_list_prody.append(pep.select(' name CA and resindex %d' % i ).getResnames()[0])
        print("prody ignore(B):",pep_list_prody)
        
        return interacting_chainA_res,interacting_chainB_res

    # interacting_chainA_res=np.array(interacting_chainA_res)
    return interacting_chainA_res
   '''
