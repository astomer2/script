import MDAnalysis as mda
import MDAnalysis.analysis.pca as pca
import numpy as np
from pathlib import Path
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import sys
import os
import logging 
import time
sys.path.append(os.path.abspath('.'))
from utils_comm.log_util import get_logger
from utils_peptide.pdb_utils import add_chain_id_to_pdb
logger = get_logger(__name__)

def energy_landscape(top,traj):

    topology = Path(top)
    trajectory = Path(traj)
    u= mda.Universe(topology=topology,trajectory=trajectory)

    u.transfer_to_memory(2000, len(u.trajectory),5)
    
    logger.info(f"start 2000 and skip 5 frames to get {len(u.trajectory)} frames")

    pca_ = pca.PCA(u, select='backbone')
    pca_.run()
    atomgroup = u.select_atoms('backbone')
    pca_space = pca_.transform(atomgroup, n_components=2)
    pc1 = pca_space[:, 0]
    pc2 = pca_space[:, 1]

    x = pc1
    y = pc2
    k = gaussian_kde(np.vstack([x, y]))
    xi, yi = np.mgrid[x.min():x.max():x.size ** 0.5 * 1j, y.min():y.max():y.size ** 0.5 * 1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    plt.contourf(xi, yi, zi.reshape(xi.shape), cmap='Blues')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.colorbar(label='Density')
    plt.title('Free Energy Landscape')
    plt.show()

    density = zi.reshape(xi.shape)
    max_density_index = np.unravel_index(np.argmax(density), density.shape)
    max_density_point = xi[max_density_index], yi[max_density_index]

    #ToDo: find local_max_desity_point
    distances = np.sqrt(np.sum((pca_space - np.array(max_density_point))**2, axis=1))
    min_dist_idx = np.argmin(distances)

    ag = u.select_atoms("all")
    ag.write(f"{top.parent}/most_stable_conformation.pdb",frames = [min_dist_idx])

    logger.info(f"most stable conformation saved to {top.parent}/most_stable_conformation.pdb")
    return (top.parent/'most_stable_conformation.pdb')

if __name__ == "__main__":
    top = '/mnt/nas1/lanwei-125/MC5R/MD/AFD_v4_part_2/GRTRTAPF/fixed.prmtop'
    traj = '/mnt/nas1/lanwei-125/MC5R/MD/AFD_v4_part_2/GRTRTAPF/fixed.nc'
    most_stable_conformation = energy_landscape(top, traj)
    aa_positions = 275
    fix_pdb = add_chain_id_to_pdb(most_stable_conformation, most_stable_conformation.parent/'fixed_chain_id.pdb',aa_positions )
