import pathlib
import argparse

from Bio.PDB.Polypeptide import is_aa
from Bio.PDB import PDBParser, Structure, Residue, Chain, Atom

import numpy as np
from matplotlib import pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file_path", type=pathlib.Path)
    return parser.parse_args()


def read_pdb(file_path: pathlib.Path, file_name: str) -> Structure.Structure:
    parser = PDBParser()
    structure = parser.get_structure(file_name, file_path)
    return structure


def calc_dist(atom1: Atom, atom2: Atom) -> float:
    diff_vector = atom1.coord - atom2.coord
    return np.sqrt(np.sum(diff_vector * diff_vector))


def calc_dist_matrix(structure: Structure.Structure) -> np.array:
    ca_atoms = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    ca_atoms.append(residue['CA'])
    
    num_atoms = len(ca_atoms)

    dist_matrix = np.zeros((num_atoms, num_atoms))
    
    for i, atom1 in enumerate(ca_atoms):
        for j, atom2 in enumerate(ca_atoms):
            if i != j:      
                dist_matrix[i, j] = calc_dist(atom1, atom2)
            else:
                dist_matrix[i, j] = 0.0
                
    return dist_matrix

def plot_distance_matrix(dist_matrix: list[list[float]], hot_bar: bool = True) -> plt:
    plt.figure(figsize=(10, 10)) 
    if hot_bar:
        plt.imshow(dist_matrix, cmap="hot", interpolation="nearest")
        plt.colorbar(label="Distance (A)")
    else:
        plt.imshow(dist_matrix, cmap="binary", interpolation="nearest")
        plt.colorbar(label="Distance below 0.8 (A)")
    plt.title("Distance Matrix")
    plt.xlabel("Residue Index")
    plt.ylabel("Residue Index")
    return plt


def main():
    args = parse_args()
    file_name = args.file_path.name.split(".")[0].upper()
    structure = read_pdb(args.file_path, file_name)

    chains = list(structure.get_chains())
    if not chains:
        raise ValueError("Chain not found")

    dist_matrix = calc_dist_matrix(structure)
    contact_map = dist_matrix < 8.0

    contact_map_color = plot_distance_matrix(dist_matrix)
    contact_map_color.savefig(f"{file_name}_distance_matrix_color.png")

    contact_map_black = plot_distance_matrix(contact_map, False)
    contact_map_black.savefig(f"{file_name}_contact_map_black.png")


if __name__ == "__main__":
    main()
