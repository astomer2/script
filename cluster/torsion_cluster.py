from ast import parse
from pkg_resources import resource_stream
import os
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import pandas as DataFrame
from shapely.geometry import Point, Polygon
import pandas as pd
from pandas import DataFrame
from pkg_resources import resource_stream
from math import pi
from Bio.SeqUtils import seq1
from Levenshtein import distance
from sklearn.cluster import OPTICS
import pandas as pd
import warnings
import matplotlib

def degrees(rad_angle) :
    """Converts any angle in radians to degrees.

    If the input is None, then it returns None.
    For numerical input, the output is mapped to [-180,180]
    """
    if rad_angle is None :
        return None
    angle = rad_angle * 180 / pi
    while angle > 180 :
        angle = angle - 360
    while angle < -180 :
        angle = angle + 360
    return angle

def get_ignored_res(file: str):
    x, y, ignored, output = [], [], [], {}
    for model in PDBParser().get_structure(id=None, file=file):
        for chain in model:
            peptides = PPBuilder().build_peptides(chain)
            for peptide in peptides:
                for aa, angles in zip(peptide, peptide.get_phi_psi_list()):
                    residue = aa.resname + str(aa.id[1])
                    output[residue] = angles

    for key, value in output.items():
        # Only get residues with both phi and psi angles
        if value[0] and value[1]:
            x.append(value[0] * 180 / pi)
            y.append(value[1] * 180 / pi)
        else:
            ignored.append((key, value))

    return output, ignored, x, y

def define_Rama_region():
    Z = np.fromfile(resource_stream('RamachanDraw', 'data/density_estimate.dat'))
    Z = np.reshape(Z, (100, 100))
    data = np.log(np.rot90(Z))
    plt.figure()
    ax = plt.subplot(111)
    ax.imshow(data,  extent=[-180, 180, -180, 180])
    data = np.rot90(np.fliplr(Z))
    contour_set =ax.contour(data, colors='k', linewidths=0.5,
                levels=[10 ** (i / 1.73) for i in range(-10, 1)],
                antialiased=True, extent=[-180, 180, -180, 180], alpha=0.65)
    contour_paths = contour_set.collections[0].get_paths()
    phi_psi=[]
    for path in contour_paths:
        vertices = path.vertices
        x_values, y_values = vertices[:, 0], vertices[:, 1]
        for x, y in zip(x_values, y_values):
            a=x,y
            phi_psi.append(a)

    region_A=[]
    region_B=[]
    region_B1=[]
    region_X=[]
    region_Y=[]
    region_Y1=[]
    for i in phi_psi:
        phi_values=i[0]
        psi_values = i[1]
        if phi_values<0 and -79<psi_values<52:
            region_A.append(i)
            region_X.append((-phi_values, -psi_values))
        if phi_values<0 and 52<psi_values<180:
            region_B.append(i)
            region_Y.append((-phi_values, -psi_values))
        if phi_values<0 and psi_values<-100:
            region_B1.append(i)
            region_Y1.append((-phi_values, -psi_values))
            
    polar_point=(-180,180)
    mirror_polar_point=(180,-180)
    region_B.append(polar_point)
    region_Y.append(mirror_polar_point)
    region=region_A,region_B,region_B1,region_X,region_Y,region_Y1
    return region

def points2Polygon(region):
    region_A,region_B,region_B1,region_X,region_Y,region_Y1=region
    polygon_A = Polygon(region_A)
    polygon_B = Polygon(region_B)
    polygon_X = Polygon(region_X)
    polygon_Y = Polygon(region_Y)

    polygon_B1 = Polygon(region_B1)
    polygon_Y1 = Polygon(region_Y1)
    polygon=polygon_A, polygon_B, polygon_X,polygon_Y,polygon_B1,polygon_Y1
    return polygon

def assign_region(phi, psi,region):
    point = Point(phi, psi)
    Polygon=points2Polygon(region)
    polygon_A, polygon_B, polygon_X,polygon_Y,polygon_B1,polygon_Y1=Polygon
    region_A,region_B,region_B1,region_X,region_Y,region_Y1=region
    if polygon_A.contains(point):
        return 'A'
    elif polygon_B.contains(point):
        return 'B'
    elif polygon_B1.contains(point):
        return 'B'    
    elif polygon_X.contains(point):
        return 'X'
    elif polygon_Y.contains(point):
        return 'Y'
    elif polygon_Y1.contains(point):
        return 'Y'
    else:
        nearest_point_A = min(region_A, key=lambda p: point.distance(Point(p)))
        nearest_point_B = min(region_B, key=lambda p: point.distance(Point(p)))
        nearest_point_X = min(region_X, key=lambda p: point.distance(Point(p)))
        nearest_point_Y = min(region_Y, key=lambda p: point.distance(Point(p)))
        nearest_point_B1=min(region_B1, key=lambda p: point.distance(Point(p)))
        nearest_point_Y1 = min(region_Y1, key=lambda p: point.distance(Point(p)))
        
        distances = {
            'A': point.distance(Point(nearest_point_A)),
            'B': point.distance(Point(nearest_point_B)),
            'B': point.distance(Point(nearest_point_B1)),
            'X': point.distance(Point(nearest_point_X)),
            'Y': point.distance(Point(nearest_point_Y)),
            'Y': point.distance(Point(nearest_point_Y1)),
        }

        nearest_region = min(distances, key=distances.get)
        return nearest_region


def main(pdb_dir) -> DataFrame or None:
    Sequences=[]
    torsion_sequences=[]
    region=define_Rama_region()

    for pdb in Path(pdb_dir).glob('*.pdb'):
        chains=[]
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(file=pdb, id=None)

        for model in structure:
            ppb = PPBuilder()
            for pp in ppb.build_peptides(model):
                sequence = pp.get_sequence()
                Sequences.append(str(sequence))
            for chain in model:
                chains.append(chain.id)
        if len(chains)>1:
            print("please process your comlpex structure to sigle chain")
            return
        else :
            _, _, x, y=get_ignored_res(str(pdb))
            torsion_str=[]
            for phi,psi in zip(x, y):
                tosion_region_str=assign_region(phi,psi,region)
                torsion_str.append(tosion_region_str)
            torsion_sequence="".join(torsion_str)
        torsion_sequences.append(torsion_sequence)
    df = DataFrame({'Sequences': Sequences, 'Torsion_Sequences': torsion_sequences})
    
    return df

def torsion_bin_cluster(df):
    warnings.filterwarnings("ignore", category=matplotlib.MatplotlibDeprecationWarning)
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    length_dfs = {}
    for length, group_df in df.groupby(df['Torsion_Sequences'].str.len()):
        length_dfs[length] = group_df
    for length, length_df in length_dfs.items():
        distances = []
        for i, seq1 in enumerate(length_df['Torsion_Sequences']):
            row_distances = [distance(seq1, seq2) for seq2 in length_df['Torsion_Sequences']]
            distances.append(row_distances)
        clustering =OPTICS(min_samples=2)
        clusters = clustering.fit_predict(distances)
        # Add clusters to DataFrame
        length_df['Cluster'] = clusters
        # Output Sequences and Torsion_Sequences of cluster centers to a new DataFrame
        cluster_centers = length_df.groupby('Cluster').apply(lambda group: group.loc[group['Cluster'].idxmin()])
        output_df = cluster_centers[['Sequences', 'Torsion_Sequences']]
        output_df.to_csv(Path(test_path)/f'length_{length}_clusters.csv', index=False)


if __name__ == '__main__':
    test_path='/mnt/nas/alphafold/af_out/cyclic_peptide/FGF5-all-pos_disulfide_peptide_pdbs/'
    tosion_bin=main(test_path)
    torsion_bin_cluster(tosion_bin)
    
