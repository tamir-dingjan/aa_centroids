import MDAnalysis as mda
import numpy as np

# Load the file
u = mda.Universe("AF-Q8N5B7-F1-model_v4.pdb")

# Get centroids for aromatic sidechains
aromatic_resids = np.unique(u.select_atoms("resname PHE TRP TYR").resids)


# Save out all the centroids as dummy particles
def get_centroid_of_side_chain(atomgroup):
    """
    Get the centroid of a side chain.
    """
    # Get the coordinates of the atoms in the side chain
    coords = atomgroup.positions
    # Calculate the centroid
    centroid = np.mean(coords, axis=0)
    return centroid


centroids = []

for resid in aromatic_resids:
    # Select the side chain atoms
    side_chain = u.select_atoms(
        f"resid {resid} and (name CG or name CD1 or name CD2 or name CE1 or name CE2 or name CE3 or name NE1 or name CH2 or name CZ2 or name CZ3 or name CZ)"
    )
    # Get the centroid of the side chain
    centroid = get_centroid_of_side_chain(side_chain)
    centroids.append(centroid)


# Get all-to-all distances
def get_all_to_all_distances(centroids):
    """
    Get all-to-all distances between centroids.
    """
    # Convert centroids to a numpy array
    centroids = np.array(centroids)
    # Calculate the distances
    distances = np.linalg.norm(centroids[:, np.newaxis] - centroids, axis=2)
    return distances


centroids = np.array(centroids)

# Save centroids to PDB formatted file
with open("centroids.pdb", "w") as f:
    for i, centroid in enumerate(centroids):
        f.write(
            f"ATOM  {i+1:5d}  CEN ALA A{i+1:4d}    {centroid[0]:8.3f}{centroid[1]:8.3f}{centroid[2]:8.3f}  1.00  0.00\n"
        )


distances = get_all_to_all_distances(centroids)

# Save the distances to a file
np.savetxt("distances.txt", distances, delimiter=",", fmt="%.3f")
