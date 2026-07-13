# ==============================================================================
# gmsh2porouslab
#
# This function exports a 2D mesh created in Gmsh to a MATLAB .mat file that is
# compatible with the PorousLab framework. It supports meshes composed of 
# multiple element types (e.g., TRI3 and QUAD4), exporting the node coordinates
# and the element connectivity in the expected format.
#
# The exported format includes:
#   - node: [N x 2] array with x and y coordinates of mesh nodes
#   - elem: {Nelem x 1} cell array, where each cell is a list of node indices
#
# Authors:
#   Danilo Cavalcanti (dborges@cimne.upc.edu)
# ==============================================================================

import numpy as np
import scipy.io

def gmsh2porouslab(gmsh_model, physical_group_id, output_path="mesh.mat"):
    """
    Export Gmsh mesh to PorousLab .mat format with:
      - node: [N x 2] array of coordinates
      - elem: {Nelem x 1} cell array (each entry is a list of node indices)

    Parameters:
        gmsh_model: gmsh.model
        physical_group_id: int, physical group ID for the mesh domain
        output_path (str): name of the .mat file to export
    """

    # Get the nodes
    nodeTags, nodeCoords = gmsh_model.mesh.getNodesForPhysicalGroup(dim=2, tag=physical_group_id)

    # Reshape the node coordinates into a matrix
    nodeCoords = nodeCoords.reshape((-1,3))

    # Convert nodeCoords to a format suitable for MATLAB (e.g., a numpy array)
    NODE = np.array(nodeCoords)
    NODE = NODE[:, :2]

    # Get the 2D elements
    elemTypes, elemTags, elemNodeTags = gmsh_model.mesh.getElements(dim=2, tag=-1)

    # Fill the elem list with the elements connectivity
    elem = []
    for etype, elem_connectivity in zip(elemTypes, elemNodeTags):
        # Get number of nodes per element
        _, _, _, num_nodes_per_elem, *_ = gmsh_model.mesh.getElementProperties(etype)
        
        # Reshape and map to local indices
        elem_connectivity = elem_connectivity.reshape((-1, num_nodes_per_elem))
        
        # Append each element as a list (row vector)
        for row in elem_connectivity:
            elem.append(row.tolist())

    # Convert to column cell array for MATLAB
    ELEM = np.empty((len(elem), 1), dtype=object)
    for i, connectivity in enumerate(elem):
        ELEM[i, 0] = connectivity  

    # Save as .mat file
    scipy.io.savemat(output_path, {
        "node": NODE,
        "elem": ELEM
    })

    print(f"[export_gmsh_to_matlab] Mesh exported to '{output_path}'.")