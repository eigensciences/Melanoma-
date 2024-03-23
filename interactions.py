#this script is used to compute protein-protein hydrogen bond interactions
# usage: python interactions.py [Protein Name]

from Bio.PDB import PDBParser, NeighborSearch
import numpy as np

# Define the PDB file name and the two chains to analyze
pdb_filename = 'p.pdb'
chain_id_1 = 'A'  # Change to your first chain ID
chain_id_2 = 'B'  # Change to your second chain ID

# Define criteria for hydrogen bonds
max_distance = 4.0  # Maximum distance for a hydrogen bond (in Ångstroms)
min_angle = 90  # Minimum angle for a hydrogen bond (in degrees)

# Parse the PDB file
parser = PDBParser(QUIET=True)
structure = parser.get_structure('structure', pdb_filename)
model = structure[0]  # Assuming only one model in the PDB file

# Get the two chains
chain1 = model[chain_id_1]
chain2 = model[chain_id_2]

# List of all atoms for distance calculations
all_atoms = list(model.get_atoms())

# Create a NeighborSearch object for all atoms
searcher = NeighborSearch(all_atoms)

# Function to calculate the angle (in degrees) between three points
def calc_angle(donor, h, acceptor):
    v1 = h.get_vector() - donor.get_vector()
    v2 = acceptor.get_vector() - h.get_vector()
    angle = np.degrees(v1.angle(v2))
    return angle

# Identify hydrogen bonds
hydrogen_bonds = []
for donor in chain1.get_atoms():
    if donor.element not in ['N', 'O']:
        continue  # Focus on nitrogen and oxygen donors
    hydrogens = searcher.search(donor.coord, 1.2, level='A')  # Find nearby hydrogens
    for h in hydrogens:
        if h.element != 'H':
            continue  # Ensure it's a hydrogen
        acceptors = searcher.search(h.coord, max_distance, level='A')  # Potential acceptors for hydrogen
        for acceptor in acceptors:
            if acceptor not in chain2.get_atoms():
                continue  # Ensure the acceptor is in the second chain
            angle = calc_angle(donor, h, acceptor)
            if angle >= min_angle:
                hydrogen_bonds.append((donor, h, acceptor))

# Print out the hydrogen bonds found
for bond in hydrogen_bonds:
    donor, h, acceptor = bond
    print(f"Hydrogen bond: {donor.parent.id}{donor.name} - {h.name} ... {acceptor.parent.id}{acceptor.name}")
