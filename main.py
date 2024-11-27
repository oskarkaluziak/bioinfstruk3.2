import numpy as np
from Bio import PDB
import matplotlib.pyplot as plt

pdb_id = '4YWO'
parser = PDB.PDBParser(QUIET=True)
structure = parser.get_structure(pdb_id, f"{pdb_id}.pdb")

def calc_phi_psi_with_dihedral(structure):
    phi_psi_angles = []
    for model in structure:
        for chain in model:
            residues = list(chain.get_residues())
            for i in range(1, len(residues) - 1):
                try:
                    prev_C = residues[i - 1]['C']
                    this_N = residues[i]['N']
                    this_CA = residues[i]['CA']
                    this_C = residues[i]['C']
                    next_N = residues[i + 1]['N']

                    phi = PDB.calc_dihedral(
                        prev_C.get_vector(),
                        this_N.get_vector(),
                        this_CA.get_vector(),
                        this_C.get_vector()
                    )

                    psi = PDB.calc_dihedral(
                        this_N.get_vector(),
                        this_CA.get_vector(),
                        this_C.get_vector(),
                        next_N.get_vector()
                    )

                    phi_psi_angles.append((phi, psi))

                except KeyError:
                    #brak atomu = pomijamy reszte
                    pass

    return phi_psi_angles

phi_psi_angles = calc_phi_psi_with_dihedral(structure)
phi, psi = zip(*[(np.degrees(p[0]), np.degrees(p[1])) for p in phi_psi_angles if None not in p])

plt.figure(figsize=(8, 8))
plt.scatter(phi, psi, alpha=0.5, color='black')
plt.xlim(-200, 200)
plt.ylim(-200, 200)
plt.axhline(y=180, color='red', linestyle='--', linewidth=1)
plt.axvline(x=180, color='red', linestyle='--', linewidth=1)
plt.axhline(y=-180, color='red', linestyle='--', linewidth=1)
plt.axvline(x=-180, color='red', linestyle='--', linewidth=1)
plt.title('Wykres Ramachandrana')
plt.xlabel('Kąt phi (°)')
plt.ylabel('Kąt psi (°)')
plt.grid()
plt.show()