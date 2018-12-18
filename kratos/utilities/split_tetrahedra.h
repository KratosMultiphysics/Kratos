//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//                   Ruben Zorrilla
//                    
//

#if !defined(KRATOS_SPLIT_TETRAHEDRA)
#define  KRATOS_SPLIT_TETRAHEDRA
 
namespace Kratos
{

/**
 * @brief This class performs the splitting of a tetrahedra.
 * It contains three helper functions to ease the splitting:
 * TetrahedraSplitMode, Split_Tetrahedra, and TetrahedraGetNewConnectivityGID.
 * 
 * EXAMPLE: imagine that a user would like to split a tetrahedra formed by the ids 3 9 7 12 
 * by introducing a new node 15 on the edge between 9 and 7. Then he should define:
 * 
 * int aux[11];
 * int edge_ids[6];
 * int t[56];
 * 
 * and initialize aux[11] as:
 * 
 * aux[0] = 3; 
 * aux[1] = 9; 
 * aux[2] = 7; 
 * aux[3] = 12;
 * aux[4] = -1; // Edge 01 --> edge not to be refined
 * aux[5] = 15; // Edge 02 --> this is the only edge in where we want to add a node
 * aux[6] = -1; // Edge 03 --> edge not to be refined
 * aux[7] = -1; // Edge 12 --> edge not to be refined
 * aux[8] = -1; // Edge 13 --> edge not to be refined
 * aux[9] = -1; // Edge 23 --> edge not to be refined
 * aux[10] = -1; // Internal (Steiner) node (set it to -1 for the moment)
 *
 * then call:
 *
 * TetrahedraSplitMode(aux, edge_ids);
 * int nel; //number of element generated\n
 * int number_splitted_edges; //number of splitted edges\n
 * int nint; //number of internal nodes\n
 * bool split_needed = Split_Tetrahedra(edge_ids, t, &nel, &number_splitted_edges, &nint)
 *
 * if (nint == 1) { // we need to generate a new internal node
 *     aux[10] = new node id
 * }
 *
 * the new tetrahedra ids can be then inspected by calling
 * 
 * for(int i = 0; i < nel; ++i) {
 *     int i0,i1,i2,i3;
 *     TetrahedraGetNewConnectivityGID(i, t, aux, &i0, &i1, &i2, &i3);
 * }
 */

class TetrahedraSplit
{
public:

    /**
     * @brief Returns the edges vector filled with the splitting pattern. Provided the array of nodal ids, 
     * this method computes the edge by edge intersection and collapse pattern and saves it in a c array.
     * @param aux_ids Array with 11 components containing the nodal ids. 4 components for the vertices, 
     * 6 for the edge nodes (which value is negative if the edge is not split) and an extra one for those
     * cases in where the internal (Steiner) node is needed:
     * aux_ids[0] = id of the FIRST node of the original tetra
     * aux_ids[1] = id of the SECOND node of the original tetra
     * aux_ids[2] = id of the THIRD node of the original tetra
     * aux_ids[3] = id of the FOURTH node of the original tetra
     * aux_ids[4] = id of new node to be used for the edge 01 (-1 if edge not to be splitted)
     * aux_ids[5] = id of new node to be used for the edge 02 (-1 if edge not to be splitted)
     * aux_ids[6] = id of new node to be used for the edge 03 (-1 if edge not to be splitted)
     * aux_ids[7] = id of new node to be used for the edge 12 (-1 if edge not to be splitted)
     * aux_ids[8] = id of new node to be used for the edge 13 (-1 if edge not to be splitted)
     * aux_ids[9] = id of new node to be used for the edge 23 (-1 if edge not to be splitted)
     * aux_ids[10] = id of new internal node (can be setted by the user or if needed, that is, if nint=1)
     * @param edge_ids Array with 6 components in where the collapse local numbering pattern is saved. 
     */
    static void TetrahedraSplitMode(
        int aux_ids[11],
        int edge_ids[6])
    {
        // Edge 01
        if (aux_ids[4] < 0) {
            edge_ids[0] = aux_ids[0] > aux_ids[1] ? 0 : 1;
        } else {
            edge_ids[0] = 4;
        }

        // Edge 02
        if (aux_ids[5] < 0) {
            edge_ids[1] = aux_ids[0] > aux_ids[2] ? 0 : 2;
        } else {
            edge_ids[1] = 5;
        }

        // Edge 03
        if (aux_ids[6] < 0) {
            edge_ids[2] = aux_ids[0] > aux_ids[3] ? 0 : 3;
        } else {
            edge_ids[2] = 6;
        }

        // Edge 12
        if (aux_ids[7] < 0) {
            edge_ids[3] = aux_ids[1] > aux_ids[2] ? 1 : 2;
        } else {
            edge_ids[3] = 7;
        }

        // Edge 13
        if (aux_ids[8] < 0) {
            edge_ids[4] = aux_ids[1] > aux_ids[3] ? 1 : 3;
        } else {
            edge_ids[4] = 8;
        }

        // Edge 23
        if (aux_ids[9] < 0) {
            edge_ids[5] = aux_ids[2] > aux_ids[3] ? 2 : 3;
        } else {
            edge_ids[5] = 9;
        }
    }

    /**
     * @brief Returns the ids of a subtetra
     * Provided the splitting connectivities array and the array containing
     * the nodal ids, this function returns the global ids for a given subtetrahedron
     * @param tetra_index Index of the new subtetrahedron
     * @param t Integer array containing the connectivities (provided by Split_Tetrahedra)
     * @param aux_ids Array used in constructing the edge_ids
     * @param id0 Global id of node 0 of the new tetrahedron
     * @param id1 Global id of node 1 of the new tetrahedron
     * @param id2 Global id of node 2 of the new tetrahedron
     * @param id3 Global id of node 3 of the new tetrahedron
     */
    static inline void TetrahedraGetNewConnectivityGID(
        const int tetra_index,
        const int t[56],
        const int aux_ids[11],
        int* id0,
        int* id1,
        int* id2,
        int* id3)
    {
        unsigned int base = tetra_index * 4;
        *id0 = aux_ids[t[base]];
        *id1 = aux_ids[t[base + 1]];
        *id2 = aux_ids[t[base + 2]];
        *id3 = aux_ids[t[base + 3]];
    }

    /**
     * @brief Function to split a tetrahedron
     * For a given edges splitting pattern, this function computes the internal tetrahedral subdivision. 
     * @param edges Input array of size 6 containing the edges splitting pattern.
     * Edges are enumerated as 01 02 03 12 13 23
     * @param t Integer c array of size 56 (= 14*4) containing the splitting connectivities
     * @param nel Number of elements (subdivisions) generated
     * @param splitted_edges Number of splitted edges
     * @param nint If 1, the internal (Steiner) node is needed
     * @return int 1 if splitting is needed and 0 if not
     */
    static int Split_Tetrahedra(
        const int edges[6],
        int t[56],
        int* nel,
        int* splitted_edges,
        int* nint)
    {

        // Check split edges
        *(splitted_edges) = 0;
        for (unsigned int i = 0; i < 6; i++) {
            if (edges[i] > 3) {
                *splitted_edges = *splitted_edges + 1;
            }
        }

        // The internal (Steiner) node is normally not needed so by default we set to false
        (*nint) = 0;

        // No splitting needed (keep connectivities)
        if (*splitted_edges == 0) {
            *nel = 1;
            t[0] = 0;
            t[1] = 1;
            t[2] = 2;
            t[3] = 3;
            return 0; // False
        }

        // Get the split connectivities
        if (edges[0] == 0) {
            if (edges[1] == 0) {
                if (edges[2] == 0) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 1;
								t[4] = 9;
								t[5] = 0;
								t[6] = 1;
								t[7] = 3;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 1;
								t[4] = 9;
								t[5] = 0;
								t[6] = 1;
								t[7] = 3;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 8;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 8;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 9;
								t[1] = 0;
								t[2] = 8;
								t[3] = 3;
								t[4] = 0;
								t[5] = 2;
								t[6] = 9;
								t[7] = 1;
								t[8] = 0;
								t[9] = 9;
								t[10] = 8;
								t[11] = 1;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 1;
								t[4] = 9;
								t[5] = 0;
								t[6] = 1;
								t[7] = 3;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 1;
								t[4] = 9;
								t[5] = 0;
								t[6] = 1;
								t[7] = 3;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 8;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 8;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 8;
								t[1] = 2;
								t[2] = 0;
								t[3] = 9;
								t[4] = 9;
								t[5] = 0;
								t[6] = 8;
								t[7] = 3;
								t[8] = 8;
								t[9] = 2;
								t[10] = 1;
								t[11] = 0;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 0;
								t[6] = 7;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 0;
								t[6] = 7;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 7;
								t[4] = 0;
								t[5] = 9;
								t[6] = 3;
								t[7] = 1;
								t[8] = 9;
								t[9] = 0;
								t[10] = 7;
								t[11] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 0;
								t[6] = 7;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 0;
								t[6] = 7;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 7;
								t[4] = 3;
								t[5] = 7;
								t[6] = 0;
								t[7] = 9;
								t[8] = 3;
								t[9] = 7;
								t[10] = 1;
								t[11] = 0;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 7;
								t[1] = 0;
								t[2] = 1;
								t[3] = 8;
								t[4] = 2;
								t[5] = 0;
								t[6] = 7;
								t[7] = 8;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 3;
								t[1] = 7;
								t[2] = 0;
								t[3] = 2;
								t[4] = 7;
								t[5] = 0;
								t[6] = 1;
								t[7] = 8;
								t[8] = 3;
								t[9] = 7;
								t[10] = 8;
								t[11] = 0;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 3;
								t[1] = 0;
								t[2] = 9;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 9;
								t[7] = 7;
								t[8] = 9;
								t[9] = 0;
								t[10] = 7;
								t[11] = 8;
								t[12] = 0;
								t[13] = 1;
								t[14] = 7;
								t[15] = 8;
                            }
                        }
                    }
                } else if (edges[2] == 3) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 1;
								t[4] = 9;
								t[5] = 0;
								t[6] = 1;
								t[7] = 3;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 1;
								t[4] = 9;
								t[5] = 0;
								t[6] = 1;
								t[7] = 3;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 8;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 8;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 9;
								t[1] = 0;
								t[2] = 8;
								t[3] = 3;
								t[4] = 0;
								t[5] = 2;
								t[6] = 9;
								t[7] = 1;
								t[8] = 0;
								t[9] = 9;
								t[10] = 8;
								t[11] = 1;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 1;
								t[4] = 9;
								t[5] = 0;
								t[6] = 1;
								t[7] = 3;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 1;
								t[4] = 9;
								t[5] = 0;
								t[6] = 1;
								t[7] = 3;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 8;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 8;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 8;
								t[1] = 2;
								t[2] = 0;
								t[3] = 9;
								t[4] = 9;
								t[5] = 0;
								t[6] = 8;
								t[7] = 3;
								t[8] = 8;
								t[9] = 2;
								t[10] = 1;
								t[11] = 0;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 0;
								t[6] = 7;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 0;
								t[6] = 7;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 7;
								t[4] = 0;
								t[5] = 9;
								t[6] = 3;
								t[7] = 1;
								t[8] = 9;
								t[9] = 0;
								t[10] = 7;
								t[11] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 0;
								t[6] = 7;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 0;
								t[6] = 7;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 7;
								t[4] = 3;
								t[5] = 7;
								t[6] = 0;
								t[7] = 9;
								t[8] = 3;
								t[9] = 7;
								t[10] = 1;
								t[11] = 0;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 7;
								t[1] = 0;
								t[2] = 1;
								t[3] = 8;
								t[4] = 2;
								t[5] = 0;
								t[6] = 7;
								t[7] = 8;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 3;
								t[1] = 7;
								t[2] = 0;
								t[3] = 2;
								t[4] = 7;
								t[5] = 0;
								t[6] = 1;
								t[7] = 8;
								t[8] = 3;
								t[9] = 7;
								t[10] = 8;
								t[11] = 0;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 3;
								t[1] = 0;
								t[2] = 9;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 9;
								t[7] = 7;
								t[8] = 9;
								t[9] = 0;
								t[10] = 7;
								t[11] = 8;
								t[12] = 0;
								t[13] = 1;
								t[14] = 7;
								t[15] = 8;
                            }
                        }
                    }
                } else if (edges[2] == 6) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 1;
								t[4] = 2;
								t[5] = 6;
								t[6] = 1;
								t[7] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 1;
								t[4] = 2;
								t[5] = 6;
								t[6] = 1;
								t[7] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 1;
								t[4] = 9;
								t[5] = 1;
								t[6] = 3;
								t[7] = 6;
								t[8] = 9;
								t[9] = 0;
								t[10] = 1;
								t[11] = 6;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 1;
								t[4] = 2;
								t[5] = 6;
								t[6] = 1;
								t[7] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 1;
								t[4] = 2;
								t[5] = 6;
								t[6] = 1;
								t[7] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 1;
								t[4] = 9;
								t[5] = 1;
								t[6] = 3;
								t[7] = 6;
								t[8] = 9;
								t[9] = 0;
								t[10] = 1;
								t[11] = 6;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 3;
								t[2] = 6;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 6;
								t[7] = 8;
								t[8] = 2;
								t[9] = 0;
								t[10] = 1;
								t[11] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 3;
								t[2] = 6;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 6;
								t[7] = 8;
								t[8] = 2;
								t[9] = 0;
								t[10] = 1;
								t[11] = 8;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 1;
								t[2] = 9;
								t[3] = 8;
								t[4] = 3;
								t[5] = 6;
								t[6] = 9;
								t[7] = 8;
								t[8] = 2;
								t[9] = 0;
								t[10] = 1;
								t[11] = 9;
								t[12] = 6;
								t[13] = 0;
								t[14] = 9;
								t[15] = 8;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 1;
								t[4] = 2;
								t[5] = 6;
								t[6] = 1;
								t[7] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 1;
								t[4] = 2;
								t[5] = 6;
								t[6] = 1;
								t[7] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 1;
								t[4] = 9;
								t[5] = 1;
								t[6] = 3;
								t[7] = 6;
								t[8] = 9;
								t[9] = 0;
								t[10] = 1;
								t[11] = 6;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 1;
								t[4] = 2;
								t[5] = 6;
								t[6] = 1;
								t[7] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 1;
								t[4] = 2;
								t[5] = 6;
								t[6] = 1;
								t[7] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 1;
								t[4] = 9;
								t[5] = 1;
								t[6] = 3;
								t[7] = 6;
								t[8] = 9;
								t[9] = 0;
								t[10] = 1;
								t[11] = 6;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 3;
								t[2] = 6;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 6;
								t[7] = 8;
								t[8] = 2;
								t[9] = 0;
								t[10] = 1;
								t[11] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 3;
								t[2] = 6;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 6;
								t[7] = 8;
								t[8] = 2;
								t[9] = 0;
								t[10] = 1;
								t[11] = 8;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 8;
								t[1] = 2;
								t[2] = 1;
								t[3] = 0;
								t[4] = 8;
								t[5] = 2;
								t[6] = 0;
								t[7] = 9;
								t[8] = 3;
								t[9] = 6;
								t[10] = 9;
								t[11] = 8;
								t[12] = 6;
								t[13] = 0;
								t[14] = 9;
								t[15] = 8;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 1;
								t[1] = 0;
								t[2] = 6;
								t[3] = 7;
								t[4] = 2;
								t[5] = 3;
								t[6] = 6;
								t[7] = 7;
								t[8] = 3;
								t[9] = 1;
								t[10] = 6;
								t[11] = 7;
								t[12] = 0;
								t[13] = 2;
								t[14] = 6;
								t[15] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 1;
								t[1] = 0;
								t[2] = 6;
								t[3] = 7;
								t[4] = 2;
								t[5] = 3;
								t[6] = 6;
								t[7] = 7;
								t[8] = 3;
								t[9] = 1;
								t[10] = 6;
								t[11] = 7;
								t[12] = 0;
								t[13] = 2;
								t[14] = 6;
								t[15] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 1;
								t[2] = 6;
								t[3] = 9;
								t[4] = 0;
								t[5] = 2;
								t[6] = 9;
								t[7] = 7;
								t[8] = 1;
								t[9] = 0;
								t[10] = 6;
								t[11] = 7;
								t[12] = 6;
								t[13] = 0;
								t[14] = 9;
								t[15] = 7;
								t[16] = 1;
								t[17] = 6;
								t[18] = 9;
								t[19] = 7;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 1;
								t[1] = 0;
								t[2] = 6;
								t[3] = 7;
								t[4] = 2;
								t[5] = 3;
								t[6] = 6;
								t[7] = 7;
								t[8] = 3;
								t[9] = 1;
								t[10] = 6;
								t[11] = 7;
								t[12] = 0;
								t[13] = 2;
								t[14] = 6;
								t[15] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 1;
								t[1] = 0;
								t[2] = 6;
								t[3] = 7;
								t[4] = 2;
								t[5] = 3;
								t[6] = 6;
								t[7] = 7;
								t[8] = 3;
								t[9] = 1;
								t[10] = 6;
								t[11] = 7;
								t[12] = 0;
								t[13] = 2;
								t[14] = 6;
								t[15] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 7;
								t[2] = 1;
								t[3] = 6;
								t[4] = 0;
								t[5] = 2;
								t[6] = 9;
								t[7] = 7;
								t[8] = 3;
								t[9] = 7;
								t[10] = 6;
								t[11] = 9;
								t[12] = 1;
								t[13] = 0;
								t[14] = 6;
								t[15] = 7;
								t[16] = 6;
								t[17] = 0;
								t[18] = 9;
								t[19] = 7;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 8;
								t[2] = 1;
								t[3] = 7;
								t[4] = 3;
								t[5] = 2;
								t[6] = 8;
								t[7] = 6;
								t[8] = 8;
								t[9] = 0;
								t[10] = 6;
								t[11] = 7;
								t[12] = 2;
								t[13] = 8;
								t[14] = 6;
								t[15] = 7;
								t[16] = 0;
								t[17] = 2;
								t[18] = 6;
								t[19] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 8;
								t[2] = 1;
								t[3] = 7;
								t[4] = 3;
								t[5] = 7;
								t[6] = 6;
								t[7] = 2;
								t[8] = 3;
								t[9] = 7;
								t[10] = 8;
								t[11] = 6;
								t[12] = 8;
								t[13] = 0;
								t[14] = 6;
								t[15] = 7;
								t[16] = 0;
								t[17] = 2;
								t[18] = 6;
								t[19] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 8;
								t[1] = 3;
								t[2] = 9;
								t[3] = 6;
								t[4] = 0;
								t[5] = 8;
								t[6] = 9;
								t[7] = 6;
								t[8] = 0;
								t[9] = 2;
								t[10] = 9;
								t[11] = 7;
								t[12] = 0;
								t[13] = 8;
								t[14] = 1;
								t[15] = 7;
								t[16] = 8;
								t[17] = 0;
								t[18] = 9;
								t[19] = 7;
                            }
                        }
                    }
                }
            } else if (edges[1] == 2) {
                if (edges[2] == 0) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 1;
								t[4] = 9;
								t[5] = 0;
								t[6] = 1;
								t[7] = 3;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 1;
								t[4] = 9;
								t[5] = 0;
								t[6] = 1;
								t[7] = 3;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 8;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 8;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 9;
								t[1] = 0;
								t[2] = 8;
								t[3] = 3;
								t[4] = 0;
								t[5] = 2;
								t[6] = 9;
								t[7] = 1;
								t[8] = 0;
								t[9] = 9;
								t[10] = 8;
								t[11] = 1;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 1;
								t[4] = 9;
								t[5] = 0;
								t[6] = 1;
								t[7] = 3;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 1;
								t[4] = 9;
								t[5] = 0;
								t[6] = 1;
								t[7] = 3;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 8;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 8;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 8;
								t[1] = 2;
								t[2] = 0;
								t[3] = 9;
								t[4] = 9;
								t[5] = 0;
								t[6] = 8;
								t[7] = 3;
								t[8] = 8;
								t[9] = 2;
								t[10] = 1;
								t[11] = 0;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 0;
								t[6] = 7;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 0;
								t[6] = 7;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 7;
								t[4] = 0;
								t[5] = 9;
								t[6] = 3;
								t[7] = 1;
								t[8] = 9;
								t[9] = 0;
								t[10] = 7;
								t[11] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 0;
								t[6] = 7;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 0;
								t[6] = 7;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 7;
								t[4] = 7;
								t[5] = 3;
								t[6] = 0;
								t[7] = 1;
								t[8] = 7;
								t[9] = 3;
								t[10] = 9;
								t[11] = 0;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 7;
								t[1] = 0;
								t[2] = 1;
								t[3] = 8;
								t[4] = 2;
								t[5] = 0;
								t[6] = 7;
								t[7] = 8;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 3;
								t[1] = 7;
								t[2] = 0;
								t[3] = 2;
								t[4] = 7;
								t[5] = 0;
								t[6] = 1;
								t[7] = 8;
								t[8] = 3;
								t[9] = 7;
								t[10] = 8;
								t[11] = 0;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 3;
								t[1] = 0;
								t[2] = 9;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 9;
								t[7] = 7;
								t[8] = 9;
								t[9] = 0;
								t[10] = 7;
								t[11] = 8;
								t[12] = 0;
								t[13] = 1;
								t[14] = 7;
								t[15] = 8;
                            }
                        }
                    }
                } else if (edges[2] == 3) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 1;
								t[4] = 9;
								t[5] = 0;
								t[6] = 1;
								t[7] = 3;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 1;
								t[4] = 9;
								t[5] = 0;
								t[6] = 1;
								t[7] = 3;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 8;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 8;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 9;
								t[1] = 0;
								t[2] = 8;
								t[3] = 3;
								t[4] = 0;
								t[5] = 2;
								t[6] = 9;
								t[7] = 1;
								t[8] = 0;
								t[9] = 9;
								t[10] = 8;
								t[11] = 1;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 1;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 1;
								t[4] = 9;
								t[5] = 0;
								t[6] = 1;
								t[7] = 3;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                // substitute_023242
                            } else if (edges[5] == 3) {
                                // substitute_023243
                            } else if (edges[5] == 9) {
                                // substitute_023249
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 8;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 8;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 8;
								t[2] = 0;
								t[3] = 1;
								t[4] = 9;
								t[5] = 0;
								t[6] = 8;
								t[7] = 3;
								t[8] = 2;
								t[9] = 8;
								t[10] = 9;
								t[11] = 0;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 0;
								t[6] = 7;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 0;
								t[6] = 7;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 7;
								t[4] = 0;
								t[5] = 9;
								t[6] = 3;
								t[7] = 1;
								t[8] = 9;
								t[9] = 0;
								t[10] = 7;
								t[11] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 0;
								t[6] = 7;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 0;
								t[6] = 7;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 7;
								t[4] = 7;
								t[5] = 3;
								t[6] = 0;
								t[7] = 1;
								t[8] = 7;
								t[9] = 3;
								t[10] = 9;
								t[11] = 0;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 7;
								t[1] = 0;
								t[2] = 1;
								t[3] = 8;
								t[4] = 2;
								t[5] = 0;
								t[6] = 7;
								t[7] = 8;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 3;
								t[1] = 7;
								t[2] = 0;
								t[3] = 2;
								t[4] = 7;
								t[5] = 0;
								t[6] = 1;
								t[7] = 8;
								t[8] = 3;
								t[9] = 7;
								t[10] = 8;
								t[11] = 0;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 3;
								t[1] = 0;
								t[2] = 9;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 9;
								t[7] = 7;
								t[8] = 9;
								t[9] = 0;
								t[10] = 7;
								t[11] = 8;
								t[12] = 0;
								t[13] = 1;
								t[14] = 7;
								t[15] = 8;
                            }
                        }
                    }
                } else if (edges[2] == 6) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 1;
								t[4] = 2;
								t[5] = 6;
								t[6] = 1;
								t[7] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 1;
								t[4] = 2;
								t[5] = 6;
								t[6] = 1;
								t[7] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 6;
								t[2] = 0;
								t[3] = 1;
								t[4] = 9;
								t[5] = 1;
								t[6] = 3;
								t[7] = 6;
								t[8] = 2;
								t[9] = 6;
								t[10] = 1;
								t[11] = 9;
                            } 
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 1;
								t[4] = 2;
								t[5] = 6;
								t[6] = 1;
								t[7] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 1;
								t[4] = 2;
								t[5] = 6;
								t[6] = 1;
								t[7] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 6;
								t[2] = 0;
								t[3] = 1;
								t[4] = 9;
								t[5] = 1;
								t[6] = 3;
								t[7] = 6;
								t[8] = 2;
								t[9] = 6;
								t[10] = 1;
								t[11] = 9;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 3;
								t[2] = 6;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 6;
								t[7] = 8;
								t[8] = 2;
								t[9] = 0;
								t[10] = 1;
								t[11] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 3;
								t[2] = 6;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 6;
								t[7] = 8;
								t[8] = 2;
								t[9] = 0;
								t[10] = 1;
								t[11] = 8;
                            } else if (edges[5] == 9) {
                                (*nel) = 9;
								t[0] = 1;
								t[1] = 9;
								t[2] = 2;
								t[3] = 10;
								t[4] = 8;
								t[5] = 0;
								t[6] = 10;
								t[7] = 1;
								t[8] = 2;
								t[9] = 6;
								t[10] = 0;
								t[11] = 10;
								t[12] = 2;
								t[13] = 0;
								t[14] = 1;
								t[15] = 10;
								t[16] = 3;
								t[17] = 6;
								t[18] = 9;
								t[19] = 8;
								t[20] = 6;
								t[21] = 2;
								t[22] = 9;
								t[23] = 10;
								t[24] = 6;
								t[25] = 9;
								t[26] = 8;
								t[27] = 10;
								t[28] = 8;
								t[29] = 0;
								t[30] = 6;
								t[31] = 10;
								t[32] = 9;
								t[33] = 1;
								t[34] = 8;
								t[35] = 10;
								(*nint) = 1;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 1;
								t[4] = 2;
								t[5] = 6;
								t[6] = 1;
								t[7] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 1;
								t[4] = 2;
								t[5] = 6;
								t[6] = 1;
								t[7] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 6;
								t[2] = 0;
								t[3] = 1;
								t[4] = 9;
								t[5] = 1;
								t[6] = 3;
								t[7] = 6;
								t[8] = 2;
								t[9] = 6;
								t[10] = 1;
								t[11] = 9;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 1;
								t[4] = 2;
								t[5] = 6;
								t[6] = 1;
								t[7] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 1;
								t[4] = 2;
								t[5] = 6;
								t[6] = 1;
								t[7] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 6;
								t[2] = 0;
								t[3] = 1;
								t[4] = 9;
								t[5] = 1;
								t[6] = 3;
								t[7] = 6;
								t[8] = 2;
								t[9] = 6;
								t[10] = 1;
								t[11] = 9;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 3;
								t[2] = 6;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 6;
								t[7] = 8;
								t[8] = 2;
								t[9] = 0;
								t[10] = 1;
								t[11] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 3;
								t[2] = 6;
								t[3] = 8;
								t[4] = 0;
								t[5] = 2;
								t[6] = 6;
								t[7] = 8;
								t[8] = 2;
								t[9] = 0;
								t[10] = 1;
								t[11] = 8;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 2;
								t[1] = 6;
								t[2] = 8;
								t[3] = 9;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 0;
								t[8] = 2;
								t[9] = 6;
								t[10] = 0;
								t[11] = 8;
								t[12] = 3;
								t[13] = 6;
								t[14] = 9;
								t[15] = 8;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 1;
								t[1] = 0;
								t[2] = 6;
								t[3] = 7;
								t[4] = 2;
								t[5] = 3;
								t[6] = 6;
								t[7] = 7;
								t[8] = 3;
								t[9] = 1;
								t[10] = 6;
								t[11] = 7;
								t[12] = 0;
								t[13] = 2;
								t[14] = 6;
								t[15] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 1;
								t[1] = 0;
								t[2] = 6;
								t[3] = 7;
								t[4] = 2;
								t[5] = 3;
								t[6] = 6;
								t[7] = 7;
								t[8] = 3;
								t[9] = 1;
								t[10] = 6;
								t[11] = 7;
								t[12] = 0;
								t[13] = 2;
								t[14] = 6;
								t[15] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 1;
								t[2] = 6;
								t[3] = 9;
								t[4] = 2;
								t[5] = 6;
								t[6] = 0;
								t[7] = 7;
								t[8] = 2;
								t[9] = 6;
								t[10] = 7;
								t[11] = 9;
								t[12] = 1;
								t[13] = 0;
								t[14] = 6;
								t[15] = 7;
								t[16] = 1;
								t[17] = 6;
								t[18] = 9;
								t[19] = 7;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 1;
								t[1] = 0;
								t[2] = 6;
								t[3] = 7;
								t[4] = 2;
								t[5] = 3;
								t[6] = 6;
								t[7] = 7;
								t[8] = 3;
								t[9] = 1;
								t[10] = 6;
								t[11] = 7;
								t[12] = 0;
								t[13] = 2;
								t[14] = 6;
								t[15] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 1;
								t[1] = 0;
								t[2] = 6;
								t[3] = 7;
								t[4] = 2;
								t[5] = 3;
								t[6] = 6;
								t[7] = 7;
								t[8] = 3;
								t[9] = 1;
								t[10] = 6;
								t[11] = 7;
								t[12] = 0;
								t[13] = 2;
								t[14] = 6;
								t[15] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 2;
								t[1] = 6;
								t[2] = 0;
								t[3] = 7;
								t[4] = 2;
								t[5] = 6;
								t[6] = 7;
								t[7] = 9;
								t[8] = 7;
								t[9] = 3;
								t[10] = 6;
								t[11] = 1;
								t[12] = 1;
								t[13] = 0;
								t[14] = 6;
								t[15] = 7;
								t[16] = 7;
								t[17] = 3;
								t[18] = 9;
								t[19] = 6;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 8;
								t[2] = 1;
								t[3] = 7;
								t[4] = 3;
								t[5] = 2;
								t[6] = 8;
								t[7] = 6;
								t[8] = 8;
								t[9] = 0;
								t[10] = 6;
								t[11] = 7;
								t[12] = 2;
								t[13] = 8;
								t[14] = 6;
								t[15] = 7;
								t[16] = 0;
								t[17] = 2;
								t[18] = 6;
								t[19] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 8;
								t[2] = 1;
								t[3] = 7;
								t[4] = 7;
								t[5] = 3;
								t[6] = 6;
								t[7] = 8;
								t[8] = 8;
								t[9] = 0;
								t[10] = 6;
								t[11] = 7;
								t[12] = 7;
								t[13] = 3;
								t[14] = 2;
								t[15] = 6;
								t[16] = 0;
								t[17] = 2;
								t[18] = 6;
								t[19] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 3;
								t[2] = 9;
								t[3] = 6;
								t[4] = 6;
								t[5] = 7;
								t[6] = 9;
								t[7] = 8;
								t[8] = 6;
								t[9] = 7;
								t[10] = 8;
								t[11] = 0;
								t[12] = 6;
								t[13] = 2;
								t[14] = 7;
								t[15] = 0;
								t[16] = 0;
								t[17] = 8;
								t[18] = 1;
								t[19] = 7;
								t[20] = 6;
								t[21] = 2;
								t[22] = 9;
								t[23] = 7;
                            }
                        }
                    }
                }
            } else if (edges[1] == 5) {
                if (edges[2] == 0) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 1;
								t[4] = 5;
								t[5] = 3;
								t[6] = 1;
								t[7] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 1;
								t[4] = 5;
								t[5] = 3;
								t[6] = 1;
								t[7] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 9;
								t[1] = 2;
								t[2] = 1;
								t[3] = 5;
								t[4] = 9;
								t[5] = 0;
								t[6] = 1;
								t[7] = 3;
								t[8] = 0;
								t[9] = 9;
								t[10] = 1;
								t[11] = 5;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 1;
								t[4] = 5;
								t[5] = 3;
								t[6] = 1;
								t[7] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 1;
								t[4] = 5;
								t[5] = 3;
								t[6] = 1;
								t[7] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 9;
								t[1] = 2;
								t[2] = 1;
								t[3] = 5;
								t[4] = 9;
								t[5] = 0;
								t[6] = 1;
								t[7] = 3;
								t[8] = 0;
								t[9] = 9;
								t[10] = 1;
								t[11] = 5;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 5;
								t[1] = 1;
								t[2] = 2;
								t[3] = 8;
								t[4] = 3;
								t[5] = 5;
								t[6] = 2;
								t[7] = 8;
								t[8] = 5;
								t[9] = 0;
								t[10] = 1;
								t[11] = 8;
								t[12] = 0;
								t[13] = 5;
								t[14] = 3;
								t[15] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 5;
								t[1] = 1;
								t[2] = 2;
								t[3] = 8;
								t[4] = 3;
								t[5] = 5;
								t[6] = 2;
								t[7] = 8;
								t[8] = 5;
								t[9] = 0;
								t[10] = 1;
								t[11] = 8;
								t[12] = 0;
								t[13] = 5;
								t[14] = 3;
								t[15] = 8;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 9;
								t[2] = 3;
								t[3] = 8;
								t[4] = 5;
								t[5] = 9;
								t[6] = 1;
								t[7] = 2;
								t[8] = 5;
								t[9] = 0;
								t[10] = 1;
								t[11] = 8;
								t[12] = 0;
								t[13] = 5;
								t[14] = 9;
								t[15] = 8;
								t[16] = 9;
								t[17] = 5;
								t[18] = 1;
								t[19] = 8;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 1;
								t[4] = 5;
								t[5] = 3;
								t[6] = 1;
								t[7] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 1;
								t[4] = 5;
								t[5] = 3;
								t[6] = 1;
								t[7] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 9;
								t[1] = 2;
								t[2] = 1;
								t[3] = 5;
								t[4] = 9;
								t[5] = 0;
								t[6] = 1;
								t[7] = 3;
								t[8] = 0;
								t[9] = 9;
								t[10] = 1;
								t[11] = 5;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 1;
								t[4] = 5;
								t[5] = 3;
								t[6] = 1;
								t[7] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 1;
								t[4] = 5;
								t[5] = 3;
								t[6] = 1;
								t[7] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 9;
								t[1] = 2;
								t[2] = 1;
								t[3] = 5;
								t[4] = 9;
								t[5] = 0;
								t[6] = 1;
								t[7] = 3;
								t[8] = 0;
								t[9] = 9;
								t[10] = 1;
								t[11] = 5;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 5;
								t[1] = 1;
								t[2] = 2;
								t[3] = 8;
								t[4] = 3;
								t[5] = 5;
								t[6] = 2;
								t[7] = 8;
								t[8] = 5;
								t[9] = 0;
								t[10] = 1;
								t[11] = 8;
								t[12] = 0;
								t[13] = 5;
								t[14] = 3;
								t[15] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 5;
								t[1] = 1;
								t[2] = 2;
								t[3] = 8;
								t[4] = 3;
								t[5] = 5;
								t[6] = 2;
								t[7] = 8;
								t[8] = 5;
								t[9] = 0;
								t[10] = 1;
								t[11] = 8;
								t[12] = 0;
								t[13] = 5;
								t[14] = 3;
								t[15] = 8;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 8;
								t[1] = 2;
								t[2] = 5;
								t[3] = 9;
								t[4] = 0;
								t[5] = 9;
								t[6] = 3;
								t[7] = 8;
								t[8] = 5;
								t[9] = 0;
								t[10] = 1;
								t[11] = 8;
								t[12] = 0;
								t[13] = 5;
								t[14] = 9;
								t[15] = 8;
								t[16] = 8;
								t[17] = 2;
								t[18] = 1;
								t[19] = 5;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 3;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 1;
								t[7] = 7;
								t[8] = 0;
								t[9] = 5;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 3;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 1;
								t[7] = 7;
								t[8] = 0;
								t[9] = 5;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 9;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 9;
								t[5] = 0;
								t[6] = 1;
								t[7] = 3;
								t[8] = 0;
								t[9] = 9;
								t[10] = 1;
								t[11] = 7;
								t[12] = 0;
								t[13] = 5;
								t[14] = 9;
								t[15] = 7;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 3;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 1;
								t[7] = 7;
								t[8] = 0;
								t[9] = 5;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 3;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 1;
								t[7] = 7;
								t[8] = 0;
								t[9] = 5;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 9;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 7;
								t[5] = 3;
								t[6] = 0;
								t[7] = 1;
								t[8] = 7;
								t[9] = 3;
								t[10] = 9;
								t[11] = 0;
								t[12] = 0;
								t[13] = 5;
								t[14] = 9;
								t[15] = 7;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 8;
								t[4] = 8;
								t[5] = 5;
								t[6] = 2;
								t[7] = 7;
								t[8] = 5;
								t[9] = 3;
								t[10] = 8;
								t[11] = 2;
								t[12] = 0;
								t[13] = 8;
								t[14] = 1;
								t[15] = 7;
								t[16] = 0;
								t[17] = 5;
								t[18] = 8;
								t[19] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 8;
								t[4] = 3;
								t[5] = 7;
								t[6] = 8;
								t[7] = 5;
								t[8] = 3;
								t[9] = 7;
								t[10] = 5;
								t[11] = 2;
								t[12] = 0;
								t[13] = 8;
								t[14] = 1;
								t[15] = 7;
								t[16] = 0;
								t[17] = 5;
								t[18] = 8;
								t[19] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 3;
								t[2] = 8;
								t[3] = 9;
								t[4] = 9;
								t[5] = 5;
								t[6] = 2;
								t[7] = 7;
								t[8] = 0;
								t[9] = 8;
								t[10] = 1;
								t[11] = 7;
								t[12] = 0;
								t[13] = 5;
								t[14] = 9;
								t[15] = 7;
								t[16] = 8;
								t[17] = 0;
								t[18] = 9;
								t[19] = 7;
                            }
                        }
                    }
                } else if (edges[2] == 3) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 1;
								t[4] = 5;
								t[5] = 3;
								t[6] = 1;
								t[7] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 1;
								t[4] = 5;
								t[5] = 3;
								t[6] = 1;
								t[7] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 9;
								t[1] = 2;
								t[2] = 1;
								t[3] = 5;
								t[4] = 5;
								t[5] = 3;
								t[6] = 1;
								t[7] = 9;
								t[8] = 5;
								t[9] = 3;
								t[10] = 0;
								t[11] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 1;
								t[4] = 5;
								t[5] = 3;
								t[6] = 1;
								t[7] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 1;
								t[4] = 5;
								t[5] = 3;
								t[6] = 1;
								t[7] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 9;
								t[1] = 2;
								t[2] = 1;
								t[3] = 5;
								t[4] = 5;
								t[5] = 3;
								t[6] = 1;
								t[7] = 9;
								t[8] = 5;
								t[9] = 3;
								t[10] = 0;
								t[11] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 5;
								t[1] = 1;
								t[2] = 2;
								t[3] = 8;
								t[4] = 3;
								t[5] = 5;
								t[6] = 2;
								t[7] = 8;
								t[8] = 5;
								t[9] = 0;
								t[10] = 1;
								t[11] = 8;
								t[12] = 0;
								t[13] = 5;
								t[14] = 3;
								t[15] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 5;
								t[1] = 1;
								t[2] = 2;
								t[3] = 8;
								t[4] = 3;
								t[5] = 5;
								t[6] = 2;
								t[7] = 8;
								t[8] = 5;
								t[9] = 0;
								t[10] = 1;
								t[11] = 8;
								t[12] = 0;
								t[13] = 5;
								t[14] = 3;
								t[15] = 8;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 5;
								t[1] = 3;
								t[2] = 8;
								t[3] = 9;
								t[4] = 5;
								t[5] = 9;
								t[6] = 1;
								t[7] = 2;
								t[8] = 5;
								t[9] = 0;
								t[10] = 1;
								t[11] = 8;
								t[12] = 5;
								t[13] = 3;
								t[14] = 0;
								t[15] = 8;
								t[16] = 9;
								t[17] = 5;
								t[18] = 1;
								t[19] = 8;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 1;
								t[4] = 5;
								t[5] = 3;
								t[6] = 1;
								t[7] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 1;
								t[4] = 5;
								t[5] = 3;
								t[6] = 1;
								t[7] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 9;
								t[1] = 2;
								t[2] = 1;
								t[3] = 5;
								t[4] = 5;
								t[5] = 3;
								t[6] = 1;
								t[7] = 9;
								t[8] = 5;
								t[9] = 3;
								t[10] = 0;
								t[11] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 1;
								t[4] = 5;
								t[5] = 3;
								t[6] = 1;
								t[7] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 1;
								t[4] = 5;
								t[5] = 3;
								t[6] = 1;
								t[7] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 3;
								t[0] = 9;
								t[1] = 2;
								t[2] = 1;
								t[3] = 5;
								t[4] = 5;
								t[5] = 3;
								t[6] = 1;
								t[7] = 9;
								t[8] = 5;
								t[9] = 3;
								t[10] = 0;
								t[11] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 5;
								t[1] = 1;
								t[2] = 2;
								t[3] = 8;
								t[4] = 3;
								t[5] = 5;
								t[6] = 2;
								t[7] = 8;
								t[8] = 5;
								t[9] = 0;
								t[10] = 1;
								t[11] = 8;
								t[12] = 0;
								t[13] = 5;
								t[14] = 3;
								t[15] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 5;
								t[1] = 1;
								t[2] = 2;
								t[3] = 8;
								t[4] = 3;
								t[5] = 5;
								t[6] = 2;
								t[7] = 8;
								t[8] = 5;
								t[9] = 0;
								t[10] = 1;
								t[11] = 8;
								t[12] = 0;
								t[13] = 5;
								t[14] = 3;
								t[15] = 8;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 8;
								t[1] = 2;
								t[2] = 5;
								t[3] = 9;
								t[4] = 5;
								t[5] = 3;
								t[6] = 8;
								t[7] = 9;
								t[8] = 5;
								t[9] = 0;
								t[10] = 1;
								t[11] = 8;
								t[12] = 5;
								t[13] = 3;
								t[14] = 0;
								t[15] = 8;
								t[16] = 8;
								t[17] = 2;
								t[18] = 1;
								t[19] = 5;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 3;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 1;
								t[7] = 7;
								t[8] = 0;
								t[9] = 5;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 3;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 1;
								t[7] = 7;
								t[8] = 0;
								t[9] = 5;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 9;
								t[0] = 9;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 5;
								t[5] = 3;
								t[6] = 0;
								t[7] = 10;
								t[8] = 1;
								t[9] = 9;
								t[10] = 10;
								t[11] = 3;
								t[12] = 3;
								t[13] = 5;
								t[14] = 9;
								t[15] = 10;
								t[16] = 0;
								t[17] = 7;
								t[18] = 5;
								t[19] = 10;
								t[20] = 9;
								t[21] = 5;
								t[22] = 7;
								t[23] = 10;
								t[24] = 1;
								t[25] = 0;
								t[26] = 3;
								t[27] = 10;
								t[28] = 1;
								t[29] = 9;
								t[30] = 7;
								t[31] = 10;
								t[32] = 7;
								t[33] = 0;
								t[34] = 1;
								t[35] = 10;
								(*nint) = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 3;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 1;
								t[7] = 7;
								t[8] = 0;
								t[9] = 5;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 3;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 1;
								t[7] = 7;
								t[8] = 0;
								t[9] = 5;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 9;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 3;
								t[5] = 5;
								t[6] = 7;
								t[7] = 0;
								t[8] = 3;
								t[9] = 0;
								t[10] = 7;
								t[11] = 1;
								t[12] = 5;
								t[13] = 3;
								t[14] = 7;
								t[15] = 9;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 8;
								t[4] = 8;
								t[5] = 5;
								t[6] = 2;
								t[7] = 7;
								t[8] = 5;
								t[9] = 3;
								t[10] = 8;
								t[11] = 2;
								t[12] = 0;
								t[13] = 8;
								t[14] = 1;
								t[15] = 7;
								t[16] = 0;
								t[17] = 5;
								t[18] = 8;
								t[19] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 8;
								t[4] = 3;
								t[5] = 7;
								t[6] = 8;
								t[7] = 5;
								t[8] = 3;
								t[9] = 7;
								t[10] = 5;
								t[11] = 2;
								t[12] = 0;
								t[13] = 8;
								t[14] = 1;
								t[15] = 7;
								t[16] = 0;
								t[17] = 5;
								t[18] = 8;
								t[19] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 5;
								t[1] = 8;
								t[2] = 7;
								t[3] = 9;
								t[4] = 9;
								t[5] = 5;
								t[6] = 2;
								t[7] = 7;
								t[8] = 5;
								t[9] = 3;
								t[10] = 8;
								t[11] = 9;
								t[12] = 0;
								t[13] = 8;
								t[14] = 1;
								t[15] = 7;
								t[16] = 5;
								t[17] = 8;
								t[18] = 0;
								t[19] = 7;
								t[20] = 5;
								t[21] = 3;
								t[22] = 0;
								t[23] = 8;
                            }
                        }
                    }
                } else if (edges[2] == 6) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 6;
								t[1] = 2;
								t[2] = 1;
								t[3] = 5;
								t[4] = 2;
								t[5] = 6;
								t[6] = 1;
								t[7] = 3;
								t[8] = 0;
								t[9] = 6;
								t[10] = 1;
								t[11] = 5;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 3;
								t[1] = 5;
								t[2] = 1;
								t[3] = 6;
								t[4] = 3;
								t[5] = 5;
								t[6] = 2;
								t[7] = 1;
								t[8] = 0;
								t[9] = 6;
								t[10] = 1;
								t[11] = 5;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 5;
								t[1] = 6;
								t[2] = 0;
								t[3] = 1;
								t[4] = 5;
								t[5] = 6;
								t[6] = 1;
								t[7] = 9;
								t[8] = 1;
								t[9] = 9;
								t[10] = 6;
								t[11] = 3;
								t[12] = 5;
								t[13] = 9;
								t[14] = 1;
								t[15] = 2;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 6;
								t[1] = 2;
								t[2] = 1;
								t[3] = 5;
								t[4] = 2;
								t[5] = 6;
								t[6] = 1;
								t[7] = 3;
								t[8] = 0;
								t[9] = 6;
								t[10] = 1;
								t[11] = 5;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 5;
								t[1] = 3;
								t[2] = 6;
								t[3] = 1;
								t[4] = 5;
								t[5] = 3;
								t[6] = 1;
								t[7] = 2;
								t[8] = 0;
								t[9] = 6;
								t[10] = 1;
								t[11] = 5;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 6;
								t[1] = 5;
								t[2] = 1;
								t[3] = 0;
								t[4] = 6;
								t[5] = 5;
								t[6] = 9;
								t[7] = 1;
								t[8] = 1;
								t[9] = 9;
								t[10] = 6;
								t[11] = 3;
								t[12] = 5;
								t[13] = 9;
								t[14] = 1;
								t[15] = 2;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 8;
								t[1] = 5;
								t[2] = 2;
								t[3] = 1;
								t[4] = 5;
								t[5] = 0;
								t[6] = 8;
								t[7] = 6;
								t[8] = 8;
								t[9] = 5;
								t[10] = 6;
								t[11] = 2;
								t[12] = 3;
								t[13] = 8;
								t[14] = 6;
								t[15] = 2;
								t[16] = 0;
								t[17] = 5;
								t[18] = 8;
								t[19] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 8;
								t[1] = 5;
								t[2] = 2;
								t[3] = 1;
								t[4] = 3;
								t[5] = 5;
								t[6] = 8;
								t[7] = 6;
								t[8] = 5;
								t[9] = 0;
								t[10] = 8;
								t[11] = 6;
								t[12] = 3;
								t[13] = 5;
								t[14] = 2;
								t[15] = 8;
								t[16] = 0;
								t[17] = 5;
								t[18] = 8;
								t[19] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 5;
								t[1] = 6;
								t[2] = 8;
								t[3] = 9;
								t[4] = 8;
								t[5] = 6;
								t[6] = 3;
								t[7] = 9;
								t[8] = 5;
								t[9] = 6;
								t[10] = 0;
								t[11] = 8;
								t[12] = 8;
								t[13] = 0;
								t[14] = 5;
								t[15] = 1;
								t[16] = 9;
								t[17] = 8;
								t[18] = 5;
								t[19] = 1;
								t[20] = 2;
								t[21] = 9;
								t[22] = 5;
								t[23] = 1;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 6;
								t[1] = 2;
								t[2] = 1;
								t[3] = 5;
								t[4] = 2;
								t[5] = 6;
								t[6] = 1;
								t[7] = 3;
								t[8] = 0;
								t[9] = 6;
								t[10] = 1;
								t[11] = 5;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 3;
								t[1] = 5;
								t[2] = 1;
								t[3] = 6;
								t[4] = 3;
								t[5] = 5;
								t[6] = 2;
								t[7] = 1;
								t[8] = 0;
								t[9] = 6;
								t[10] = 1;
								t[11] = 5;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 5;
								t[1] = 6;
								t[2] = 0;
								t[3] = 1;
								t[4] = 5;
								t[5] = 6;
								t[6] = 1;
								t[7] = 9;
								t[8] = 1;
								t[9] = 9;
								t[10] = 6;
								t[11] = 3;
								t[12] = 5;
								t[13] = 9;
								t[14] = 1;
								t[15] = 2;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 6;
								t[1] = 2;
								t[2] = 1;
								t[3] = 5;
								t[4] = 2;
								t[5] = 6;
								t[6] = 1;
								t[7] = 3;
								t[8] = 0;
								t[9] = 6;
								t[10] = 1;
								t[11] = 5;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 5;
								t[1] = 3;
								t[2] = 6;
								t[3] = 1;
								t[4] = 5;
								t[5] = 3;
								t[6] = 1;
								t[7] = 2;
								t[8] = 0;
								t[9] = 6;
								t[10] = 1;
								t[11] = 5;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 6;
								t[1] = 5;
								t[2] = 1;
								t[3] = 0;
								t[4] = 6;
								t[5] = 5;
								t[6] = 9;
								t[7] = 1;
								t[8] = 1;
								t[9] = 9;
								t[10] = 6;
								t[11] = 3;
								t[12] = 5;
								t[13] = 9;
								t[14] = 1;
								t[15] = 2;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 8;
								t[1] = 5;
								t[2] = 2;
								t[3] = 1;
								t[4] = 5;
								t[5] = 0;
								t[6] = 8;
								t[7] = 6;
								t[8] = 8;
								t[9] = 5;
								t[10] = 6;
								t[11] = 2;
								t[12] = 3;
								t[13] = 8;
								t[14] = 6;
								t[15] = 2;
								t[16] = 0;
								t[17] = 5;
								t[18] = 8;
								t[19] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 8;
								t[1] = 5;
								t[2] = 2;
								t[3] = 1;
								t[4] = 3;
								t[5] = 5;
								t[6] = 8;
								t[7] = 6;
								t[8] = 5;
								t[9] = 0;
								t[10] = 8;
								t[11] = 6;
								t[12] = 3;
								t[13] = 5;
								t[14] = 2;
								t[15] = 8;
								t[16] = 0;
								t[17] = 5;
								t[18] = 8;
								t[19] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 2;
								t[2] = 5;
								t[3] = 9;
								t[4] = 8;
								t[5] = 6;
								t[6] = 3;
								t[7] = 9;
								t[8] = 5;
								t[9] = 6;
								t[10] = 0;
								t[11] = 8;
								t[12] = 8;
								t[13] = 0;
								t[14] = 5;
								t[15] = 1;
								t[16] = 5;
								t[17] = 6;
								t[18] = 8;
								t[19] = 9;
								t[20] = 8;
								t[21] = 2;
								t[22] = 1;
								t[23] = 5;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 7;
								t[1] = 3;
								t[2] = 6;
								t[3] = 1;
								t[4] = 5;
								t[5] = 0;
								t[6] = 7;
								t[7] = 6;
								t[8] = 7;
								t[9] = 5;
								t[10] = 6;
								t[11] = 2;
								t[12] = 3;
								t[13] = 7;
								t[14] = 6;
								t[15] = 2;
								t[16] = 0;
								t[17] = 7;
								t[18] = 6;
								t[19] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 7;
								t[1] = 3;
								t[2] = 6;
								t[3] = 1;
								t[4] = 5;
								t[5] = 0;
								t[6] = 7;
								t[7] = 6;
								t[8] = 3;
								t[9] = 5;
								t[10] = 7;
								t[11] = 6;
								t[12] = 3;
								t[13] = 5;
								t[14] = 2;
								t[15] = 7;
								t[16] = 0;
								t[17] = 7;
								t[18] = 6;
								t[19] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 2;
								t[1] = 7;
								t[2] = 9;
								t[3] = 5;
								t[4] = 6;
								t[5] = 5;
								t[6] = 7;
								t[7] = 0;
								t[8] = 6;
								t[9] = 5;
								t[10] = 9;
								t[11] = 7;
								t[12] = 6;
								t[13] = 0;
								t[14] = 7;
								t[15] = 1;
								t[16] = 3;
								t[17] = 6;
								t[18] = 9;
								t[19] = 1;
								t[20] = 6;
								t[21] = 7;
								t[22] = 9;
								t[23] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 7;
								t[1] = 3;
								t[2] = 6;
								t[3] = 1;
								t[4] = 5;
								t[5] = 0;
								t[6] = 7;
								t[7] = 6;
								t[8] = 7;
								t[9] = 5;
								t[10] = 6;
								t[11] = 2;
								t[12] = 3;
								t[13] = 7;
								t[14] = 6;
								t[15] = 2;
								t[16] = 0;
								t[17] = 7;
								t[18] = 6;
								t[19] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 7;
								t[1] = 3;
								t[2] = 6;
								t[3] = 1;
								t[4] = 5;
								t[5] = 0;
								t[6] = 7;
								t[7] = 6;
								t[8] = 5;
								t[9] = 3;
								t[10] = 7;
								t[11] = 2;
								t[12] = 5;
								t[13] = 3;
								t[14] = 6;
								t[15] = 7;
								t[16] = 0;
								t[17] = 7;
								t[18] = 6;
								t[19] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 2;
								t[1] = 7;
								t[2] = 9;
								t[3] = 5;
								t[4] = 6;
								t[5] = 5;
								t[6] = 7;
								t[7] = 0;
								t[8] = 6;
								t[9] = 5;
								t[10] = 9;
								t[11] = 7;
								t[12] = 6;
								t[13] = 0;
								t[14] = 7;
								t[15] = 1;
								t[16] = 7;
								t[17] = 3;
								t[18] = 6;
								t[19] = 1;
								t[20] = 7;
								t[21] = 3;
								t[22] = 9;
								t[23] = 6;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 0;
								t[2] = 5;
								t[3] = 7;
								t[4] = 6;
								t[5] = 2;
								t[6] = 8;
								t[7] = 5;
								t[8] = 2;
								t[9] = 6;
								t[10] = 8;
								t[11] = 3;
								t[12] = 0;
								t[13] = 6;
								t[14] = 8;
								t[15] = 5;
								t[16] = 2;
								t[17] = 8;
								t[18] = 5;
								t[19] = 7;
								t[20] = 8;
								t[21] = 0;
								t[22] = 7;
								t[23] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 6;
								t[0] = 5;
								t[1] = 8;
								t[2] = 0;
								t[3] = 7;
								t[4] = 8;
								t[5] = 5;
								t[6] = 0;
								t[7] = 6;
								t[8] = 5;
								t[9] = 8;
								t[10] = 3;
								t[11] = 6;
								t[12] = 3;
								t[13] = 5;
								t[14] = 2;
								t[15] = 7;
								t[16] = 5;
								t[17] = 8;
								t[18] = 7;
								t[19] = 3;
								t[20] = 8;
								t[21] = 0;
								t[22] = 7;
								t[23] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 7;
								t[0] = 8;
								t[1] = 3;
								t[2] = 9;
								t[3] = 6;
								t[4] = 5;
								t[5] = 8;
								t[6] = 9;
								t[7] = 6;
								t[8] = 5;
								t[9] = 8;
								t[10] = 0;
								t[11] = 7;
								t[12] = 8;
								t[13] = 5;
								t[14] = 0;
								t[15] = 6;
								t[16] = 8;
								t[17] = 0;
								t[18] = 7;
								t[19] = 1;
								t[20] = 2;
								t[21] = 9;
								t[22] = 5;
								t[23] = 7;
								t[24] = 5;
								t[25] = 8;
								t[26] = 7;
								t[27] = 9;
                            }
                        }
                    }
                }
            }
        } else if (edges[0] == 1) {
            if (edges[1] == 0) {
                if (edges[2] == 0) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 2;
								t[5] = 0;
								t[6] = 4;
								t[7] = 8;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 2;
								t[5] = 0;
								t[6] = 4;
								t[7] = 8;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 8;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 1;
								t[5] = 9;
								t[6] = 4;
								t[7] = 8;
								t[8] = 2;
								t[9] = 1;
								t[10] = 9;
								t[11] = 4;
								t[12] = 3;
								t[13] = 0;
								t[14] = 9;
								t[15] = 8;
								t[16] = 9;
								t[17] = 0;
								t[18] = 4;
								t[19] = 8;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 2;
								t[5] = 0;
								t[6] = 4;
								t[7] = 8;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 2;
								t[5] = 0;
								t[6] = 4;
								t[7] = 8;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 8;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 8;
								t[9] = 2;
								t[10] = 4;
								t[11] = 9;
								t[12] = 3;
								t[13] = 0;
								t[14] = 9;
								t[15] = 8;
								t[16] = 9;
								t[17] = 0;
								t[18] = 4;
								t[19] = 8;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 7;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 7;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 1;
								t[1] = 3;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 1;
								t[6] = 4;
								t[7] = 7;
								t[8] = 3;
								t[9] = 0;
								t[10] = 9;
								t[11] = 4;
								t[12] = 0;
								t[13] = 9;
								t[14] = 4;
								t[15] = 7;
								t[16] = 0;
								t[17] = 2;
								t[18] = 9;
								t[19] = 7;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 7;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 7;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 7;
								t[2] = 1;
								t[3] = 4;
								t[4] = 3;
								t[5] = 0;
								t[6] = 9;
								t[7] = 4;
								t[8] = 3;
								t[9] = 7;
								t[10] = 4;
								t[11] = 9;
								t[12] = 0;
								t[13] = 9;
								t[14] = 4;
								t[15] = 7;
								t[16] = 0;
								t[17] = 2;
								t[18] = 9;
								t[19] = 7;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 8;
								t[4] = 8;
								t[5] = 1;
								t[6] = 4;
								t[7] = 7;
								t[8] = 0;
								t[9] = 8;
								t[10] = 4;
								t[11] = 7;
								t[12] = 0;
								t[13] = 2;
								t[14] = 8;
								t[15] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 3;
								t[1] = 7;
								t[2] = 8;
								t[3] = 0;
								t[4] = 8;
								t[5] = 1;
								t[6] = 4;
								t[7] = 7;
								t[8] = 3;
								t[9] = 7;
								t[10] = 0;
								t[11] = 2;
								t[12] = 0;
								t[13] = 8;
								t[14] = 4;
								t[15] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 7;
								t[4] = 9;
								t[5] = 0;
								t[6] = 8;
								t[7] = 3;
								t[8] = 8;
								t[9] = 0;
								t[10] = 7;
								t[11] = 4;
								t[12] = 1;
								t[13] = 8;
								t[14] = 7;
								t[15] = 4;
								t[16] = 0;
								t[17] = 9;
								t[18] = 8;
								t[19] = 7;
                            }
                        }
                    }
                } else if (edges[2] == 3) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 8;
								t[8] = 3;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 8;
								t[8] = 3;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 1;
								t[5] = 9;
								t[6] = 4;
								t[7] = 8;
								t[8] = 2;
								t[9] = 1;
								t[10] = 9;
								t[11] = 4;
								t[12] = 4;
								t[13] = 3;
								t[14] = 9;
								t[15] = 0;
								t[16] = 4;
								t[17] = 3;
								t[18] = 8;
								t[19] = 9;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 8;
								t[8] = 3;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 8;
								t[8] = 3;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 4;
								t[9] = 3;
								t[10] = 9;
								t[11] = 0;
								t[12] = 4;
								t[13] = 3;
								t[14] = 8;
								t[15] = 9;
								t[16] = 8;
								t[17] = 2;
								t[18] = 4;
								t[19] = 9;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 7;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 7;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 1;
								t[1] = 3;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 1;
								t[6] = 4;
								t[7] = 7;
								t[8] = 3;
								t[9] = 0;
								t[10] = 9;
								t[11] = 4;
								t[12] = 0;
								t[13] = 9;
								t[14] = 4;
								t[15] = 7;
								t[16] = 0;
								t[17] = 2;
								t[18] = 9;
								t[19] = 7;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 7;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 7;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 7;
								t[2] = 1;
								t[3] = 4;
								t[4] = 3;
								t[5] = 0;
								t[6] = 9;
								t[7] = 4;
								t[8] = 3;
								t[9] = 7;
								t[10] = 4;
								t[11] = 9;
								t[12] = 0;
								t[13] = 9;
								t[14] = 4;
								t[15] = 7;
								t[16] = 0;
								t[17] = 2;
								t[18] = 9;
								t[19] = 7;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 9;
								t[0] = 4;
								t[1] = 3;
								t[2] = 8;
								t[3] = 10;
								t[4] = 3;
								t[5] = 2;
								t[6] = 8;
								t[7] = 10;
								t[8] = 8;
								t[9] = 1;
								t[10] = 4;
								t[11] = 7;
								t[12] = 0;
								t[13] = 2;
								t[14] = 3;
								t[15] = 10;
								t[16] = 3;
								t[17] = 4;
								t[18] = 0;
								t[19] = 10;
								t[20] = 7;
								t[21] = 0;
								t[22] = 10;
								t[23] = 2;
								t[24] = 7;
								t[25] = 0;
								t[26] = 4;
								t[27] = 10;
								t[28] = 8;
								t[29] = 2;
								t[30] = 7;
								t[31] = 10;
								t[32] = 4;
								t[33] = 8;
								t[34] = 7;
								t[35] = 10;
								(*nint) = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 3;
								t[1] = 4;
								t[2] = 7;
								t[3] = 8;
								t[4] = 8;
								t[5] = 1;
								t[6] = 4;
								t[7] = 7;
								t[8] = 3;
								t[9] = 4;
								t[10] = 0;
								t[11] = 7;
								t[12] = 7;
								t[13] = 3;
								t[14] = 2;
								t[15] = 0;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 7;
								t[4] = 9;
								t[5] = 4;
								t[6] = 7;
								t[7] = 8;
								t[8] = 4;
								t[9] = 3;
								t[10] = 9;
								t[11] = 0;
								t[12] = 3;
								t[13] = 4;
								t[14] = 9;
								t[15] = 8;
								t[16] = 1;
								t[17] = 8;
								t[18] = 7;
								t[19] = 4;
								t[20] = 9;
								t[21] = 4;
								t[22] = 0;
								t[23] = 7;
                            }
                        }
                    }
                } else if (edges[2] == 6) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 6;
								t[5] = 2;
								t[6] = 3;
								t[7] = 1;
								t[8] = 2;
								t[9] = 6;
								t[10] = 4;
								t[11] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 6;
								t[5] = 2;
								t[6] = 3;
								t[7] = 1;
								t[8] = 2;
								t[9] = 6;
								t[10] = 4;
								t[11] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 3;
								t[5] = 1;
								t[6] = 6;
								t[7] = 9;
								t[8] = 1;
								t[9] = 6;
								t[10] = 9;
								t[11] = 4;
								t[12] = 6;
								t[13] = 0;
								t[14] = 9;
								t[15] = 4;
								t[16] = 2;
								t[17] = 1;
								t[18] = 9;
								t[19] = 4;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 1;
								t[8] = 3;
								t[9] = 4;
								t[10] = 6;
								t[11] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 1;
								t[8] = 3;
								t[9] = 4;
								t[10] = 6;
								t[11] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 6;
								t[7] = 9;
								t[8] = 3;
								t[9] = 4;
								t[10] = 9;
								t[11] = 1;
								t[12] = 6;
								t[13] = 0;
								t[14] = 9;
								t[15] = 4;
								t[16] = 2;
								t[17] = 1;
								t[18] = 9;
								t[19] = 4;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 6;
								t[1] = 4;
								t[2] = 2;
								t[3] = 8;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 6;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
								t[12] = 3;
								t[13] = 2;
								t[14] = 8;
								t[15] = 6;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 6;
								t[1] = 4;
								t[2] = 2;
								t[3] = 8;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 6;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
								t[12] = 3;
								t[13] = 2;
								t[14] = 8;
								t[15] = 6;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 3;
								t[2] = 9;
								t[3] = 6;
								t[4] = 6;
								t[5] = 4;
								t[6] = 9;
								t[7] = 8;
								t[8] = 8;
								t[9] = 9;
								t[10] = 1;
								t[11] = 4;
								t[12] = 6;
								t[13] = 4;
								t[14] = 0;
								t[15] = 9;
								t[16] = 9;
								t[17] = 2;
								t[18] = 1;
								t[19] = 4;
								t[20] = 0;
								t[21] = 2;
								t[22] = 9;
								t[23] = 4;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 6;
								t[5] = 2;
								t[6] = 3;
								t[7] = 1;
								t[8] = 2;
								t[9] = 6;
								t[10] = 4;
								t[11] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 6;
								t[5] = 2;
								t[6] = 3;
								t[7] = 1;
								t[8] = 2;
								t[9] = 6;
								t[10] = 4;
								t[11] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 3;
								t[5] = 1;
								t[6] = 6;
								t[7] = 9;
								t[8] = 1;
								t[9] = 6;
								t[10] = 9;
								t[11] = 4;
								t[12] = 6;
								t[13] = 0;
								t[14] = 9;
								t[15] = 4;
								t[16] = 2;
								t[17] = 1;
								t[18] = 9;
								t[19] = 4;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 1;
								t[8] = 3;
								t[9] = 4;
								t[10] = 6;
								t[11] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 1;
								t[8] = 3;
								t[9] = 4;
								t[10] = 6;
								t[11] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 6;
								t[7] = 9;
								t[8] = 3;
								t[9] = 4;
								t[10] = 9;
								t[11] = 1;
								t[12] = 6;
								t[13] = 0;
								t[14] = 9;
								t[15] = 4;
								t[16] = 2;
								t[17] = 1;
								t[18] = 9;
								t[19] = 4;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 6;
								t[1] = 4;
								t[2] = 2;
								t[3] = 8;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 6;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
								t[12] = 3;
								t[13] = 2;
								t[14] = 8;
								t[15] = 6;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 6;
								t[1] = 4;
								t[2] = 2;
								t[3] = 8;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 6;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
								t[12] = 3;
								t[13] = 2;
								t[14] = 8;
								t[15] = 6;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 2;
								t[2] = 4;
								t[3] = 9;
								t[4] = 8;
								t[5] = 3;
								t[6] = 9;
								t[7] = 6;
								t[8] = 6;
								t[9] = 4;
								t[10] = 9;
								t[11] = 8;
								t[12] = 8;
								t[13] = 2;
								t[14] = 1;
								t[15] = 4;
								t[16] = 6;
								t[17] = 4;
								t[18] = 0;
								t[19] = 9;
								t[20] = 0;
								t[21] = 2;
								t[22] = 9;
								t[23] = 4;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 7;
								t[1] = 3;
								t[2] = 6;
								t[3] = 1;
								t[4] = 2;
								t[5] = 0;
								t[6] = 7;
								t[7] = 6;
								t[8] = 3;
								t[9] = 2;
								t[10] = 7;
								t[11] = 6;
								t[12] = 6;
								t[13] = 7;
								t[14] = 1;
								t[15] = 4;
								t[16] = 0;
								t[17] = 7;
								t[18] = 6;
								t[19] = 4;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 7;
								t[1] = 3;
								t[2] = 6;
								t[3] = 1;
								t[4] = 2;
								t[5] = 0;
								t[6] = 7;
								t[7] = 6;
								t[8] = 3;
								t[9] = 2;
								t[10] = 7;
								t[11] = 6;
								t[12] = 6;
								t[13] = 7;
								t[14] = 1;
								t[15] = 4;
								t[16] = 0;
								t[17] = 7;
								t[18] = 6;
								t[19] = 4;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 2;
								t[1] = 0;
								t[2] = 7;
								t[3] = 9;
								t[4] = 9;
								t[5] = 3;
								t[6] = 6;
								t[7] = 1;
								t[8] = 9;
								t[9] = 7;
								t[10] = 1;
								t[11] = 4;
								t[12] = 6;
								t[13] = 9;
								t[14] = 1;
								t[15] = 4;
								t[16] = 0;
								t[17] = 7;
								t[18] = 9;
								t[19] = 4;
								t[20] = 0;
								t[21] = 9;
								t[22] = 6;
								t[23] = 4;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 4;
								t[2] = 6;
								t[3] = 7;
								t[4] = 2;
								t[5] = 0;
								t[6] = 7;
								t[7] = 6;
								t[8] = 3;
								t[9] = 2;
								t[10] = 7;
								t[11] = 6;
								t[12] = 0;
								t[13] = 7;
								t[14] = 6;
								t[15] = 4;
								t[16] = 3;
								t[17] = 4;
								t[18] = 7;
								t[19] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 4;
								t[2] = 6;
								t[3] = 7;
								t[4] = 2;
								t[5] = 0;
								t[6] = 7;
								t[7] = 6;
								t[8] = 3;
								t[9] = 2;
								t[10] = 7;
								t[11] = 6;
								t[12] = 0;
								t[13] = 7;
								t[14] = 6;
								t[15] = 4;
								t[16] = 3;
								t[17] = 4;
								t[18] = 7;
								t[19] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 9;
								t[1] = 4;
								t[2] = 0;
								t[3] = 7;
								t[4] = 2;
								t[5] = 0;
								t[6] = 7;
								t[7] = 9;
								t[8] = 9;
								t[9] = 4;
								t[10] = 3;
								t[11] = 6;
								t[12] = 4;
								t[13] = 9;
								t[14] = 0;
								t[15] = 6;
								t[16] = 9;
								t[17] = 4;
								t[18] = 7;
								t[19] = 3;
								t[20] = 3;
								t[21] = 4;
								t[22] = 7;
								t[23] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 7;
								t[2] = 1;
								t[3] = 4;
								t[4] = 2;
								t[5] = 0;
								t[6] = 7;
								t[7] = 6;
								t[8] = 6;
								t[9] = 4;
								t[10] = 7;
								t[11] = 8;
								t[12] = 6;
								t[13] = 4;
								t[14] = 0;
								t[15] = 7;
								t[16] = 2;
								t[17] = 3;
								t[18] = 6;
								t[19] = 8;
								t[20] = 7;
								t[21] = 2;
								t[22] = 6;
								t[23] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 7;
								t[2] = 1;
								t[3] = 4;
								t[4] = 7;
								t[5] = 3;
								t[6] = 6;
								t[7] = 8;
								t[8] = 2;
								t[9] = 0;
								t[10] = 7;
								t[11] = 6;
								t[12] = 6;
								t[13] = 4;
								t[14] = 7;
								t[15] = 8;
								t[16] = 6;
								t[17] = 4;
								t[18] = 0;
								t[19] = 7;
								t[20] = 7;
								t[21] = 3;
								t[22] = 2;
								t[23] = 6;
                            } else if (edges[5] == 9) {
                                (*nel) = 7;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 7;
								t[4] = 4;
								t[5] = 9;
								t[6] = 8;
								t[7] = 7;
								t[8] = 4;
								t[9] = 9;
								t[10] = 0;
								t[11] = 6;
								t[12] = 8;
								t[13] = 7;
								t[14] = 1;
								t[15] = 4;
								t[16] = 6;
								t[17] = 9;
								t[18] = 3;
								t[19] = 8;
								t[20] = 4;
								t[21] = 9;
								t[22] = 6;
								t[23] = 8;
								t[24] = 9;
								t[25] = 4;
								t[26] = 0;
								t[27] = 7;
                            }
                        }
                    }
                }
            } else if (edges[1] == 2) {
                if (edges[2] == 0) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 2;
								t[5] = 0;
								t[6] = 4;
								t[7] = 8;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 2;
								t[5] = 0;
								t[6] = 4;
								t[7] = 8;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 8;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 1;
								t[5] = 9;
								t[6] = 4;
								t[7] = 8;
								t[8] = 2;
								t[9] = 1;
								t[10] = 9;
								t[11] = 4;
								t[12] = 3;
								t[13] = 0;
								t[14] = 9;
								t[15] = 8;
								t[16] = 9;
								t[17] = 0;
								t[18] = 4;
								t[19] = 8;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 2;
								t[5] = 0;
								t[6] = 4;
								t[7] = 8;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 2;
								t[5] = 0;
								t[6] = 4;
								t[7] = 8;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 8;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 8;
								t[9] = 2;
								t[10] = 4;
								t[11] = 9;
								t[12] = 3;
								t[13] = 0;
								t[14] = 9;
								t[15] = 8;
								t[16] = 9;
								t[17] = 0;
								t[18] = 4;
								t[19] = 8;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 4;
								t[5] = 2;
								t[6] = 3;
								t[7] = 7;
								t[8] = 4;
								t[9] = 2;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 4;
								t[5] = 2;
								t[6] = 3;
								t[7] = 7;
								t[8] = 4;
								t[9] = 2;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 1;
								t[1] = 3;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 1;
								t[6] = 4;
								t[7] = 7;
								t[8] = 3;
								t[9] = 0;
								t[10] = 9;
								t[11] = 4;
								t[12] = 4;
								t[13] = 2;
								t[14] = 9;
								t[15] = 7;
								t[16] = 4;
								t[17] = 2;
								t[18] = 0;
								t[19] = 9;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 4;
								t[5] = 2;
								t[6] = 3;
								t[7] = 7;
								t[8] = 4;
								t[9] = 2;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 4;
								t[5] = 2;
								t[6] = 3;
								t[7] = 7;
								t[8] = 4;
								t[9] = 2;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 7;
								t[1] = 3;
								t[2] = 9;
								t[3] = 4;
								t[4] = 3;
								t[5] = 0;
								t[6] = 9;
								t[7] = 4;
								t[8] = 4;
								t[9] = 2;
								t[10] = 9;
								t[11] = 7;
								t[12] = 4;
								t[13] = 2;
								t[14] = 0;
								t[15] = 9;
								t[16] = 7;
								t[17] = 3;
								t[18] = 4;
								t[19] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 8;
								t[4] = 8;
								t[5] = 1;
								t[6] = 4;
								t[7] = 7;
								t[8] = 4;
								t[9] = 2;
								t[10] = 8;
								t[11] = 7;
								t[12] = 4;
								t[13] = 2;
								t[14] = 0;
								t[15] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 9;
								t[0] = 3;
								t[1] = 7;
								t[2] = 10;
								t[3] = 2;
								t[4] = 4;
								t[5] = 2;
								t[6] = 0;
								t[7] = 10;
								t[8] = 0;
								t[9] = 8;
								t[10] = 4;
								t[11] = 10;
								t[12] = 8;
								t[13] = 1;
								t[14] = 4;
								t[15] = 7;
								t[16] = 2;
								t[17] = 4;
								t[18] = 7;
								t[19] = 10;
								t[20] = 4;
								t[21] = 8;
								t[22] = 7;
								t[23] = 10;
								t[24] = 3;
								t[25] = 7;
								t[26] = 8;
								t[27] = 10;
								t[28] = 0;
								t[29] = 3;
								t[30] = 8;
								t[31] = 10;
								t[32] = 0;
								t[33] = 2;
								t[34] = 3;
								t[35] = 10;
								(*nint) = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 4;
								t[1] = 2;
								t[2] = 9;
								t[3] = 7;
								t[4] = 4;
								t[5] = 9;
								t[6] = 8;
								t[7] = 7;
								t[8] = 9;
								t[9] = 0;
								t[10] = 8;
								t[11] = 3;
								t[12] = 4;
								t[13] = 9;
								t[14] = 0;
								t[15] = 8;
								t[16] = 1;
								t[17] = 8;
								t[18] = 7;
								t[19] = 4;
								t[20] = 4;
								t[21] = 2;
								t[22] = 0;
								t[23] = 9;
                            }
                        }
                    }
                } else if (edges[2] == 3) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 8;
								t[8] = 3;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 8;
								t[8] = 3;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 1;
								t[5] = 9;
								t[6] = 4;
								t[7] = 8;
								t[8] = 2;
								t[9] = 1;
								t[10] = 9;
								t[11] = 4;
								t[12] = 4;
								t[13] = 3;
								t[14] = 9;
								t[15] = 0;
								t[16] = 4;
								t[17] = 3;
								t[18] = 8;
								t[19] = 9;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                // substitute_423242
                            } else if (edges[5] == 3) {
                                // substitute_423243
                            } else if (edges[5] == 9) {
                                // substitute_423249
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 8;
								t[8] = 3;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 8;
								t[8] = 3;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 2;
								t[5] = 8;
								t[6] = 9;
								t[7] = 4;
								t[8] = 3;
								t[9] = 4;
								t[10] = 9;
								t[11] = 8;
								t[12] = 2;
								t[13] = 8;
								t[14] = 4;
								t[15] = 1;
								t[16] = 3;
								t[17] = 4;
								t[18] = 0;
								t[19] = 9;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 4;
								t[5] = 2;
								t[6] = 3;
								t[7] = 7;
								t[8] = 4;
								t[9] = 2;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 4;
								t[5] = 2;
								t[6] = 3;
								t[7] = 7;
								t[8] = 4;
								t[9] = 2;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 1;
								t[1] = 3;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 1;
								t[6] = 4;
								t[7] = 7;
								t[8] = 3;
								t[9] = 0;
								t[10] = 9;
								t[11] = 4;
								t[12] = 4;
								t[13] = 2;
								t[14] = 9;
								t[15] = 7;
								t[16] = 4;
								t[17] = 2;
								t[18] = 0;
								t[19] = 9;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 4;
								t[5] = 2;
								t[6] = 3;
								t[7] = 7;
								t[8] = 4;
								t[9] = 2;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 4;
								t[5] = 2;
								t[6] = 3;
								t[7] = 7;
								t[8] = 4;
								t[9] = 2;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 7;
								t[1] = 3;
								t[2] = 9;
								t[3] = 4;
								t[4] = 4;
								t[5] = 2;
								t[6] = 9;
								t[7] = 7;
								t[8] = 3;
								t[9] = 0;
								t[10] = 9;
								t[11] = 4;
								t[12] = 7;
								t[13] = 3;
								t[14] = 4;
								t[15] = 1;
								t[16] = 4;
								t[17] = 2;
								t[18] = 0;
								t[19] = 9;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 3;
								t[1] = 4;
								t[2] = 2;
								t[3] = 8;
								t[4] = 4;
								t[5] = 3;
								t[6] = 2;
								t[7] = 0;
								t[8] = 8;
								t[9] = 1;
								t[10] = 4;
								t[11] = 7;
								t[12] = 8;
								t[13] = 4;
								t[14] = 2;
								t[15] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 8;
								t[1] = 1;
								t[2] = 4;
								t[3] = 7;
								t[4] = 3;
								t[5] = 7;
								t[6] = 4;
								t[7] = 2;
								t[8] = 3;
								t[9] = 7;
								t[10] = 8;
								t[11] = 4;
								t[12] = 4;
								t[13] = 3;
								t[14] = 2;
								t[15] = 0;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 2;
								t[1] = 4;
								t[2] = 9;
								t[3] = 0;
								t[4] = 4;
								t[5] = 3;
								t[6] = 8;
								t[7] = 9;
								t[8] = 4;
								t[9] = 2;
								t[10] = 9;
								t[11] = 7;
								t[12] = 4;
								t[13] = 3;
								t[14] = 9;
								t[15] = 0;
								t[16] = 1;
								t[17] = 8;
								t[18] = 7;
								t[19] = 4;
								t[20] = 9;
								t[21] = 4;
								t[22] = 7;
								t[23] = 8;
                            }
                        }
                    }
                } else if (edges[2] == 6) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 6;
								t[5] = 2;
								t[6] = 3;
								t[7] = 1;
								t[8] = 2;
								t[9] = 6;
								t[10] = 4;
								t[11] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 6;
								t[5] = 2;
								t[6] = 3;
								t[7] = 1;
								t[8] = 2;
								t[9] = 6;
								t[10] = 4;
								t[11] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 2;
								t[1] = 6;
								t[2] = 0;
								t[3] = 4;
								t[4] = 3;
								t[5] = 1;
								t[6] = 6;
								t[7] = 9;
								t[8] = 1;
								t[9] = 6;
								t[10] = 9;
								t[11] = 4;
								t[12] = 2;
								t[13] = 6;
								t[14] = 4;
								t[15] = 9;
								t[16] = 2;
								t[17] = 1;
								t[18] = 9;
								t[19] = 4;
                            } 
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 1;
								t[8] = 3;
								t[9] = 4;
								t[10] = 6;
								t[11] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 1;
								t[8] = 3;
								t[9] = 4;
								t[10] = 6;
								t[11] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 2;
								t[1] = 6;
								t[2] = 0;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 6;
								t[7] = 9;
								t[8] = 2;
								t[9] = 6;
								t[10] = 4;
								t[11] = 9;
								t[12] = 2;
								t[13] = 1;
								t[14] = 9;
								t[15] = 4;
								t[16] = 3;
								t[17] = 4;
								t[18] = 9;
								t[19] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 6;
								t[1] = 4;
								t[2] = 2;
								t[3] = 8;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 6;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
								t[12] = 3;
								t[13] = 2;
								t[14] = 8;
								t[15] = 6;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 6;
								t[1] = 4;
								t[2] = 2;
								t[3] = 8;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 6;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
								t[12] = 3;
								t[13] = 2;
								t[14] = 8;
								t[15] = 6;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 6;
								t[1] = 2;
								t[2] = 4;
								t[3] = 0;
								t[4] = 8;
								t[5] = 3;
								t[6] = 9;
								t[7] = 6;
								t[8] = 6;
								t[9] = 4;
								t[10] = 9;
								t[11] = 8;
								t[12] = 8;
								t[13] = 9;
								t[14] = 1;
								t[15] = 4;
								t[16] = 6;
								t[17] = 2;
								t[18] = 9;
								t[19] = 4;
								t[20] = 9;
								t[21] = 2;
								t[22] = 1;
								t[23] = 4;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 6;
								t[5] = 2;
								t[6] = 3;
								t[7] = 1;
								t[8] = 2;
								t[9] = 6;
								t[10] = 4;
								t[11] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 6;
								t[5] = 2;
								t[6] = 3;
								t[7] = 1;
								t[8] = 2;
								t[9] = 6;
								t[10] = 4;
								t[11] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 2;
								t[1] = 6;
								t[2] = 0;
								t[3] = 4;
								t[4] = 3;
								t[5] = 1;
								t[6] = 6;
								t[7] = 9;
								t[8] = 1;
								t[9] = 6;
								t[10] = 9;
								t[11] = 4;
								t[12] = 2;
								t[13] = 6;
								t[14] = 4;
								t[15] = 9;
								t[16] = 2;
								t[17] = 1;
								t[18] = 9;
								t[19] = 4;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 1;
								t[8] = 3;
								t[9] = 4;
								t[10] = 6;
								t[11] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 1;
								t[8] = 3;
								t[9] = 4;
								t[10] = 6;
								t[11] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 2;
								t[1] = 6;
								t[2] = 0;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 6;
								t[7] = 9;
								t[8] = 2;
								t[9] = 6;
								t[10] = 4;
								t[11] = 9;
								t[12] = 2;
								t[13] = 1;
								t[14] = 9;
								t[15] = 4;
								t[16] = 3;
								t[17] = 4;
								t[18] = 9;
								t[19] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 6;
								t[1] = 4;
								t[2] = 2;
								t[3] = 8;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 6;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
								t[12] = 3;
								t[13] = 2;
								t[14] = 8;
								t[15] = 6;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 6;
								t[1] = 4;
								t[2] = 2;
								t[3] = 8;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 6;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
								t[12] = 3;
								t[13] = 2;
								t[14] = 8;
								t[15] = 6;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 3;
								t[2] = 9;
								t[3] = 6;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 8;
								t[9] = 2;
								t[10] = 4;
								t[11] = 9;
								t[12] = 2;
								t[13] = 6;
								t[14] = 4;
								t[15] = 9;
								t[16] = 6;
								t[17] = 2;
								t[18] = 4;
								t[19] = 0;
								t[20] = 9;
								t[21] = 6;
								t[22] = 4;
								t[23] = 8;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 7;
								t[1] = 3;
								t[2] = 6;
								t[3] = 1;
								t[4] = 3;
								t[5] = 2;
								t[6] = 7;
								t[7] = 6;
								t[8] = 6;
								t[9] = 7;
								t[10] = 1;
								t[11] = 4;
								t[12] = 4;
								t[13] = 2;
								t[14] = 0;
								t[15] = 6;
								t[16] = 4;
								t[17] = 2;
								t[18] = 6;
								t[19] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 7;
								t[1] = 3;
								t[2] = 6;
								t[3] = 1;
								t[4] = 3;
								t[5] = 2;
								t[6] = 7;
								t[7] = 6;
								t[8] = 6;
								t[9] = 7;
								t[10] = 1;
								t[11] = 4;
								t[12] = 4;
								t[13] = 2;
								t[14] = 0;
								t[15] = 6;
								t[16] = 4;
								t[17] = 2;
								t[18] = 6;
								t[19] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 9;
								t[1] = 4;
								t[2] = 1;
								t[3] = 6;
								t[4] = 4;
								t[5] = 9;
								t[6] = 2;
								t[7] = 6;
								t[8] = 9;
								t[9] = 3;
								t[10] = 6;
								t[11] = 1;
								t[12] = 4;
								t[13] = 9;
								t[14] = 1;
								t[15] = 7;
								t[16] = 4;
								t[17] = 9;
								t[18] = 7;
								t[19] = 2;
								t[20] = 4;
								t[21] = 2;
								t[22] = 0;
								t[23] = 6;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 4;
								t[2] = 6;
								t[3] = 7;
								t[4] = 4;
								t[5] = 2;
								t[6] = 6;
								t[7] = 7;
								t[8] = 3;
								t[9] = 2;
								t[10] = 7;
								t[11] = 6;
								t[12] = 4;
								t[13] = 2;
								t[14] = 0;
								t[15] = 6;
								t[16] = 3;
								t[17] = 4;
								t[18] = 7;
								t[19] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 4;
								t[2] = 6;
								t[3] = 7;
								t[4] = 4;
								t[5] = 2;
								t[6] = 6;
								t[7] = 7;
								t[8] = 3;
								t[9] = 2;
								t[10] = 7;
								t[11] = 6;
								t[12] = 4;
								t[13] = 2;
								t[14] = 0;
								t[15] = 6;
								t[16] = 3;
								t[17] = 4;
								t[18] = 7;
								t[19] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 9;
								t[1] = 4;
								t[2] = 3;
								t[3] = 6;
								t[4] = 4;
								t[5] = 9;
								t[6] = 7;
								t[7] = 2;
								t[8] = 4;
								t[9] = 9;
								t[10] = 2;
								t[11] = 6;
								t[12] = 3;
								t[13] = 4;
								t[14] = 7;
								t[15] = 1;
								t[16] = 4;
								t[17] = 9;
								t[18] = 3;
								t[19] = 7;
								t[20] = 4;
								t[21] = 2;
								t[22] = 0;
								t[23] = 6;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 7;
								t[2] = 1;
								t[3] = 4;
								t[4] = 4;
								t[5] = 2;
								t[6] = 6;
								t[7] = 7;
								t[8] = 6;
								t[9] = 4;
								t[10] = 7;
								t[11] = 8;
								t[12] = 4;
								t[13] = 2;
								t[14] = 0;
								t[15] = 6;
								t[16] = 2;
								t[17] = 3;
								t[18] = 6;
								t[19] = 8;
								t[20] = 7;
								t[21] = 2;
								t[22] = 6;
								t[23] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 7;
								t[2] = 1;
								t[3] = 4;
								t[4] = 2;
								t[5] = 4;
								t[6] = 6;
								t[7] = 0;
								t[8] = 7;
								t[9] = 3;
								t[10] = 6;
								t[11] = 8;
								t[12] = 3;
								t[13] = 7;
								t[14] = 6;
								t[15] = 2;
								t[16] = 6;
								t[17] = 7;
								t[18] = 4;
								t[19] = 2;
								t[20] = 6;
								t[21] = 7;
								t[22] = 8;
								t[23] = 4;
                            } else if (edges[5] == 9) {
                                (*nel) = 7;
								t[0] = 2;
								t[1] = 9;
								t[2] = 4;
								t[3] = 7;
								t[4] = 9;
								t[5] = 4;
								t[6] = 7;
								t[7] = 8;
								t[8] = 2;
								t[9] = 6;
								t[10] = 4;
								t[11] = 9;
								t[12] = 6;
								t[13] = 2;
								t[14] = 4;
								t[15] = 0;
								t[16] = 8;
								t[17] = 7;
								t[18] = 1;
								t[19] = 4;
								t[20] = 6;
								t[21] = 9;
								t[22] = 3;
								t[23] = 8;
								t[24] = 6;
								t[25] = 4;
								t[26] = 9;
								t[27] = 8;
                            }
                        }
                    }
                }
            } else if (edges[1] == 5) {
                if (edges[2] == 0) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 3;
								t[2] = 4;
								t[3] = 5;
								t[4] = 3;
								t[5] = 4;
								t[6] = 5;
								t[7] = 1;
								t[8] = 2;
								t[9] = 3;
								t[10] = 5;
								t[11] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 3;
								t[2] = 4;
								t[3] = 5;
								t[4] = 3;
								t[5] = 4;
								t[6] = 5;
								t[7] = 1;
								t[8] = 2;
								t[9] = 3;
								t[10] = 5;
								t[11] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 9;
								t[1] = 5;
								t[2] = 1;
								t[3] = 4;
								t[4] = 5;
								t[5] = 9;
								t[6] = 1;
								t[7] = 2;
								t[8] = 0;
								t[9] = 5;
								t[10] = 9;
								t[11] = 4;
								t[12] = 9;
								t[13] = 1;
								t[14] = 3;
								t[15] = 4;
								t[16] = 0;
								t[17] = 9;
								t[18] = 3;
								t[19] = 4;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 3;
								t[2] = 4;
								t[3] = 5;
								t[4] = 3;
								t[5] = 4;
								t[6] = 5;
								t[7] = 1;
								t[8] = 2;
								t[9] = 3;
								t[10] = 5;
								t[11] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 3;
								t[2] = 4;
								t[3] = 5;
								t[4] = 3;
								t[5] = 4;
								t[6] = 5;
								t[7] = 1;
								t[8] = 2;
								t[9] = 3;
								t[10] = 5;
								t[11] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 9;
								t[1] = 5;
								t[2] = 1;
								t[3] = 4;
								t[4] = 5;
								t[5] = 9;
								t[6] = 1;
								t[7] = 2;
								t[8] = 0;
								t[9] = 5;
								t[10] = 9;
								t[11] = 4;
								t[12] = 9;
								t[13] = 1;
								t[14] = 3;
								t[15] = 4;
								t[16] = 0;
								t[17] = 9;
								t[18] = 3;
								t[19] = 4;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 8;
								t[4] = 8;
								t[5] = 5;
								t[6] = 1;
								t[7] = 4;
								t[8] = 0;
								t[9] = 5;
								t[10] = 8;
								t[11] = 4;
								t[12] = 5;
								t[13] = 3;
								t[14] = 8;
								t[15] = 2;
								t[16] = 8;
								t[17] = 5;
								t[18] = 2;
								t[19] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 8;
								t[4] = 8;
								t[5] = 5;
								t[6] = 1;
								t[7] = 4;
								t[8] = 0;
								t[9] = 5;
								t[10] = 8;
								t[11] = 4;
								t[12] = 5;
								t[13] = 3;
								t[14] = 8;
								t[15] = 2;
								t[16] = 8;
								t[17] = 5;
								t[18] = 2;
								t[19] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 0;
								t[1] = 3;
								t[2] = 8;
								t[3] = 9;
								t[4] = 9;
								t[5] = 5;
								t[6] = 2;
								t[7] = 1;
								t[8] = 8;
								t[9] = 0;
								t[10] = 9;
								t[11] = 4;
								t[12] = 0;
								t[13] = 5;
								t[14] = 9;
								t[15] = 4;
								t[16] = 9;
								t[17] = 5;
								t[18] = 1;
								t[19] = 4;
								t[20] = 8;
								t[21] = 9;
								t[22] = 1;
								t[23] = 4;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 5;
								t[8] = 4;
								t[9] = 2;
								t[10] = 5;
								t[11] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 3;
								t[3] = 5;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 5;
								t[8] = 2;
								t[9] = 4;
								t[10] = 1;
								t[11] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 2;
								t[1] = 4;
								t[2] = 9;
								t[3] = 5;
								t[4] = 2;
								t[5] = 4;
								t[6] = 1;
								t[7] = 9;
								t[8] = 0;
								t[9] = 5;
								t[10] = 9;
								t[11] = 4;
								t[12] = 9;
								t[13] = 1;
								t[14] = 3;
								t[15] = 4;
								t[16] = 0;
								t[17] = 9;
								t[18] = 3;
								t[19] = 4;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 5;
								t[8] = 4;
								t[9] = 2;
								t[10] = 5;
								t[11] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 5;
								t[8] = 4;
								t[9] = 2;
								t[10] = 5;
								t[11] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 2;
								t[1] = 4;
								t[2] = 9;
								t[3] = 5;
								t[4] = 2;
								t[5] = 4;
								t[6] = 1;
								t[7] = 9;
								t[8] = 0;
								t[9] = 5;
								t[10] = 9;
								t[11] = 4;
								t[12] = 9;
								t[13] = 1;
								t[14] = 3;
								t[15] = 4;
								t[16] = 0;
								t[17] = 9;
								t[18] = 3;
								t[19] = 4;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 8;
								t[4] = 4;
								t[5] = 2;
								t[6] = 5;
								t[7] = 8;
								t[8] = 0;
								t[9] = 5;
								t[10] = 8;
								t[11] = 4;
								t[12] = 5;
								t[13] = 3;
								t[14] = 8;
								t[15] = 2;
								t[16] = 4;
								t[17] = 2;
								t[18] = 8;
								t[19] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 8;
								t[4] = 4;
								t[5] = 2;
								t[6] = 5;
								t[7] = 8;
								t[8] = 0;
								t[9] = 5;
								t[10] = 8;
								t[11] = 4;
								t[12] = 5;
								t[13] = 3;
								t[14] = 8;
								t[15] = 2;
								t[16] = 4;
								t[17] = 2;
								t[18] = 8;
								t[19] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 0;
								t[1] = 3;
								t[2] = 8;
								t[3] = 9;
								t[4] = 2;
								t[5] = 8;
								t[6] = 4;
								t[7] = 1;
								t[8] = 8;
								t[9] = 0;
								t[10] = 9;
								t[11] = 4;
								t[12] = 0;
								t[13] = 5;
								t[14] = 9;
								t[15] = 4;
								t[16] = 9;
								t[17] = 2;
								t[18] = 4;
								t[19] = 5;
								t[20] = 8;
								t[21] = 2;
								t[22] = 4;
								t[23] = 9;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 4;
								t[1] = 5;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 5;
								t[9] = 3;
								t[10] = 7;
								t[11] = 2;
								t[12] = 4;
								t[13] = 5;
								t[14] = 0;
								t[15] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 4;
								t[1] = 5;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 5;
								t[9] = 3;
								t[10] = 7;
								t[11] = 2;
								t[12] = 4;
								t[13] = 5;
								t[14] = 0;
								t[15] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 3;
								t[1] = 0;
								t[2] = 9;
								t[3] = 4;
								t[4] = 3;
								t[5] = 9;
								t[6] = 1;
								t[7] = 4;
								t[8] = 7;
								t[9] = 5;
								t[10] = 9;
								t[11] = 2;
								t[12] = 4;
								t[13] = 5;
								t[14] = 9;
								t[15] = 7;
								t[16] = 9;
								t[17] = 7;
								t[18] = 1;
								t[19] = 4;
								t[20] = 4;
								t[21] = 5;
								t[22] = 0;
								t[23] = 9;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 4;
								t[1] = 5;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 5;
								t[9] = 3;
								t[10] = 7;
								t[11] = 2;
								t[12] = 4;
								t[13] = 5;
								t[14] = 0;
								t[15] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 4;
								t[1] = 5;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 5;
								t[9] = 3;
								t[10] = 7;
								t[11] = 2;
								t[12] = 4;
								t[13] = 5;
								t[14] = 0;
								t[15] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 3;
								t[1] = 0;
								t[2] = 9;
								t[3] = 4;
								t[4] = 3;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 7;
								t[9] = 5;
								t[10] = 9;
								t[11] = 2;
								t[12] = 4;
								t[13] = 5;
								t[14] = 9;
								t[15] = 7;
								t[16] = 3;
								t[17] = 7;
								t[18] = 4;
								t[19] = 9;
								t[20] = 4;
								t[21] = 5;
								t[22] = 0;
								t[23] = 9;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 7;
								t[2] = 1;
								t[3] = 4;
								t[4] = 4;
								t[5] = 5;
								t[6] = 8;
								t[7] = 7;
								t[8] = 0;
								t[9] = 5;
								t[10] = 3;
								t[11] = 8;
								t[12] = 5;
								t[13] = 7;
								t[14] = 2;
								t[15] = 8;
								t[16] = 3;
								t[17] = 5;
								t[18] = 2;
								t[19] = 8;
								t[20] = 4;
								t[21] = 5;
								t[22] = 0;
								t[23] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 7;
								t[2] = 1;
								t[3] = 4;
								t[4] = 4;
								t[5] = 5;
								t[6] = 8;
								t[7] = 7;
								t[8] = 0;
								t[9] = 5;
								t[10] = 3;
								t[11] = 8;
								t[12] = 3;
								t[13] = 7;
								t[14] = 8;
								t[15] = 5;
								t[16] = 4;
								t[17] = 5;
								t[18] = 0;
								t[19] = 8;
								t[20] = 3;
								t[21] = 7;
								t[22] = 5;
								t[23] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 7;
								t[0] = 2;
								t[1] = 7;
								t[2] = 9;
								t[3] = 5;
								t[4] = 8;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 4;
								t[9] = 9;
								t[10] = 0;
								t[11] = 8;
								t[12] = 9;
								t[13] = 4;
								t[14] = 7;
								t[15] = 8;
								t[16] = 3;
								t[17] = 0;
								t[18] = 9;
								t[19] = 8;
								t[20] = 4;
								t[21] = 9;
								t[22] = 7;
								t[23] = 5;
								t[24] = 4;
								t[25] = 9;
								t[26] = 5;
								t[27] = 0;
                            }
                        }
                    }
                } else if (edges[2] == 3) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 3;
								t[2] = 4;
								t[3] = 5;
								t[4] = 3;
								t[5] = 4;
								t[6] = 5;
								t[7] = 1;
								t[8] = 2;
								t[9] = 3;
								t[10] = 5;
								t[11] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 3;
								t[2] = 4;
								t[3] = 5;
								t[4] = 3;
								t[5] = 4;
								t[6] = 5;
								t[7] = 1;
								t[8] = 2;
								t[9] = 3;
								t[10] = 5;
								t[11] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 9;
								t[1] = 5;
								t[2] = 1;
								t[3] = 4;
								t[4] = 5;
								t[5] = 3;
								t[6] = 4;
								t[7] = 9;
								t[8] = 5;
								t[9] = 9;
								t[10] = 1;
								t[11] = 2;
								t[12] = 5;
								t[13] = 3;
								t[14] = 0;
								t[15] = 4;
								t[16] = 9;
								t[17] = 1;
								t[18] = 3;
								t[19] = 4;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 3;
								t[2] = 4;
								t[3] = 5;
								t[4] = 3;
								t[5] = 4;
								t[6] = 5;
								t[7] = 1;
								t[8] = 2;
								t[9] = 3;
								t[10] = 5;
								t[11] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 3;
								t[2] = 4;
								t[3] = 5;
								t[4] = 3;
								t[5] = 4;
								t[6] = 5;
								t[7] = 1;
								t[8] = 2;
								t[9] = 3;
								t[10] = 5;
								t[11] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 9;
								t[1] = 5;
								t[2] = 1;
								t[3] = 4;
								t[4] = 5;
								t[5] = 3;
								t[6] = 4;
								t[7] = 9;
								t[8] = 5;
								t[9] = 9;
								t[10] = 1;
								t[11] = 2;
								t[12] = 5;
								t[13] = 3;
								t[14] = 0;
								t[15] = 4;
								t[16] = 9;
								t[17] = 1;
								t[18] = 3;
								t[19] = 4;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 4;
								t[2] = 0;
								t[3] = 5;
								t[4] = 8;
								t[5] = 5;
								t[6] = 1;
								t[7] = 4;
								t[8] = 5;
								t[9] = 3;
								t[10] = 8;
								t[11] = 2;
								t[12] = 8;
								t[13] = 5;
								t[14] = 2;
								t[15] = 1;
								t[16] = 3;
								t[17] = 4;
								t[18] = 5;
								t[19] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 4;
								t[2] = 0;
								t[3] = 5;
								t[4] = 8;
								t[5] = 5;
								t[6] = 1;
								t[7] = 4;
								t[8] = 5;
								t[9] = 3;
								t[10] = 8;
								t[11] = 2;
								t[12] = 8;
								t[13] = 5;
								t[14] = 2;
								t[15] = 1;
								t[16] = 3;
								t[17] = 4;
								t[18] = 5;
								t[19] = 8;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 3;
								t[1] = 5;
								t[2] = 4;
								t[3] = 0;
								t[4] = 9;
								t[5] = 5;
								t[6] = 2;
								t[7] = 1;
								t[8] = 9;
								t[9] = 3;
								t[10] = 4;
								t[11] = 8;
								t[12] = 5;
								t[13] = 3;
								t[14] = 4;
								t[15] = 9;
								t[16] = 9;
								t[17] = 5;
								t[18] = 1;
								t[19] = 4;
								t[20] = 8;
								t[21] = 9;
								t[22] = 1;
								t[23] = 4;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 3;
								t[3] = 5;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 5;
								t[8] = 2;
								t[9] = 4;
								t[10] = 1;
								t[11] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 3;
								t[3] = 5;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 5;
								t[8] = 2;
								t[9] = 4;
								t[10] = 1;
								t[11] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 2;
								t[1] = 4;
								t[2] = 9;
								t[3] = 5;
								t[4] = 3;
								t[5] = 5;
								t[6] = 4;
								t[7] = 0;
								t[8] = 2;
								t[9] = 4;
								t[10] = 1;
								t[11] = 9;
								t[12] = 9;
								t[13] = 1;
								t[14] = 3;
								t[15] = 4;
								t[16] = 3;
								t[17] = 5;
								t[18] = 9;
								t[19] = 4;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 5;
								t[8] = 4;
								t[9] = 2;
								t[10] = 5;
								t[11] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 5;
								t[8] = 4;
								t[9] = 2;
								t[10] = 5;
								t[11] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 4;
								t[1] = 2;
								t[2] = 5;
								t[3] = 9;
								t[4] = 5;
								t[5] = 3;
								t[6] = 4;
								t[7] = 9;
								t[8] = 4;
								t[9] = 2;
								t[10] = 9;
								t[11] = 1;
								t[12] = 5;
								t[13] = 3;
								t[14] = 0;
								t[15] = 4;
								t[16] = 9;
								t[17] = 1;
								t[18] = 3;
								t[19] = 4;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 4;
								t[2] = 0;
								t[3] = 5;
								t[4] = 4;
								t[5] = 2;
								t[6] = 8;
								t[7] = 1;
								t[8] = 4;
								t[9] = 2;
								t[10] = 5;
								t[11] = 8;
								t[12] = 5;
								t[13] = 3;
								t[14] = 8;
								t[15] = 2;
								t[16] = 3;
								t[17] = 4;
								t[18] = 5;
								t[19] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 4;
								t[2] = 0;
								t[3] = 5;
								t[4] = 4;
								t[5] = 2;
								t[6] = 8;
								t[7] = 1;
								t[8] = 4;
								t[9] = 2;
								t[10] = 5;
								t[11] = 8;
								t[12] = 5;
								t[13] = 3;
								t[14] = 8;
								t[15] = 2;
								t[16] = 3;
								t[17] = 4;
								t[18] = 5;
								t[19] = 8;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 2;
								t[2] = 4;
								t[3] = 9;
								t[4] = 3;
								t[5] = 5;
								t[6] = 4;
								t[7] = 0;
								t[8] = 9;
								t[9] = 3;
								t[10] = 4;
								t[11] = 8;
								t[12] = 5;
								t[13] = 3;
								t[14] = 4;
								t[15] = 9;
								t[16] = 4;
								t[17] = 2;
								t[18] = 5;
								t[19] = 9;
								t[20] = 8;
								t[21] = 2;
								t[22] = 1;
								t[23] = 4;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 4;
								t[1] = 5;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 5;
								t[9] = 3;
								t[10] = 7;
								t[11] = 2;
								t[12] = 4;
								t[13] = 5;
								t[14] = 0;
								t[15] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 4;
								t[1] = 5;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 5;
								t[9] = 3;
								t[10] = 7;
								t[11] = 2;
								t[12] = 4;
								t[13] = 5;
								t[14] = 0;
								t[15] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 3;
								t[1] = 5;
								t[2] = 4;
								t[3] = 0;
								t[4] = 5;
								t[5] = 3;
								t[6] = 4;
								t[7] = 9;
								t[8] = 3;
								t[9] = 9;
								t[10] = 1;
								t[11] = 4;
								t[12] = 5;
								t[13] = 9;
								t[14] = 4;
								t[15] = 7;
								t[16] = 7;
								t[17] = 5;
								t[18] = 9;
								t[19] = 2;
								t[20] = 9;
								t[21] = 7;
								t[22] = 1;
								t[23] = 4;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 4;
								t[1] = 5;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 5;
								t[9] = 3;
								t[10] = 7;
								t[11] = 2;
								t[12] = 4;
								t[13] = 5;
								t[14] = 0;
								t[15] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 4;
								t[1] = 5;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 5;
								t[9] = 3;
								t[10] = 7;
								t[11] = 2;
								t[12] = 4;
								t[13] = 5;
								t[14] = 0;
								t[15] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 3;
								t[1] = 7;
								t[2] = 4;
								t[3] = 9;
								t[4] = 3;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 9;
								t[9] = 4;
								t[10] = 5;
								t[11] = 7;
								t[12] = 7;
								t[13] = 5;
								t[14] = 9;
								t[15] = 2;
								t[16] = 5;
								t[17] = 3;
								t[18] = 4;
								t[19] = 9;
								t[20] = 5;
								t[21] = 3;
								t[22] = 0;
								t[23] = 4;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 7;
								t[2] = 1;
								t[3] = 4;
								t[4] = 4;
								t[5] = 3;
								t[6] = 5;
								t[7] = 0;
								t[8] = 8;
								t[9] = 4;
								t[10] = 5;
								t[11] = 7;
								t[12] = 3;
								t[13] = 4;
								t[14] = 5;
								t[15] = 8;
								t[16] = 5;
								t[17] = 7;
								t[18] = 2;
								t[19] = 8;
								t[20] = 3;
								t[21] = 5;
								t[22] = 2;
								t[23] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 7;
								t[2] = 1;
								t[3] = 4;
								t[4] = 3;
								t[5] = 7;
								t[6] = 5;
								t[7] = 2;
								t[8] = 4;
								t[9] = 3;
								t[10] = 5;
								t[11] = 0;
								t[12] = 8;
								t[13] = 4;
								t[14] = 5;
								t[15] = 7;
								t[16] = 3;
								t[17] = 4;
								t[18] = 5;
								t[19] = 8;
								t[20] = 3;
								t[21] = 7;
								t[22] = 8;
								t[23] = 5;
                            } else if (edges[5] == 9) {
                                (*nel) = 7;
								t[0] = 2;
								t[1] = 7;
								t[2] = 9;
								t[3] = 5;
								t[4] = 8;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 4;
								t[9] = 9;
								t[10] = 7;
								t[11] = 5;
								t[12] = 9;
								t[13] = 4;
								t[14] = 7;
								t[15] = 8;
								t[16] = 4;
								t[17] = 3;
								t[18] = 5;
								t[19] = 0;
								t[20] = 4;
								t[21] = 9;
								t[22] = 5;
								t[23] = 3;
								t[24] = 4;
								t[25] = 9;
								t[26] = 3;
								t[27] = 8;
                            }
                        }
                    }
                } else if (edges[2] == 6) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 5;
								t[2] = 6;
								t[3] = 4;
								t[4] = 6;
								t[5] = 5;
								t[6] = 1;
								t[7] = 4;
								t[8] = 5;
								t[9] = 6;
								t[10] = 1;
								t[11] = 2;
								t[12] = 1;
								t[13] = 6;
								t[14] = 3;
								t[15] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 3;
								t[1] = 5;
								t[2] = 1;
								t[3] = 6;
								t[4] = 0;
								t[5] = 5;
								t[6] = 6;
								t[7] = 4;
								t[8] = 6;
								t[9] = 5;
								t[10] = 1;
								t[11] = 4;
								t[12] = 3;
								t[13] = 5;
								t[14] = 2;
								t[15] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 5;
								t[1] = 6;
								t[2] = 0;
								t[3] = 4;
								t[4] = 5;
								t[5] = 1;
								t[6] = 9;
								t[7] = 4;
								t[8] = 5;
								t[9] = 6;
								t[10] = 4;
								t[11] = 9;
								t[12] = 1;
								t[13] = 5;
								t[14] = 9;
								t[15] = 2;
								t[16] = 3;
								t[17] = 1;
								t[18] = 6;
								t[19] = 9;
								t[20] = 1;
								t[21] = 6;
								t[22] = 9;
								t[23] = 4;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 10;
								t[0] = 0;
								t[1] = 10;
								t[2] = 6;
								t[3] = 4;
								t[4] = 1;
								t[5] = 5;
								t[6] = 4;
								t[7] = 10;
								t[8] = 0;
								t[9] = 10;
								t[10] = 5;
								t[11] = 6;
								t[12] = 0;
								t[13] = 10;
								t[14] = 4;
								t[15] = 5;
								t[16] = 3;
								t[17] = 4;
								t[18] = 10;
								t[19] = 1;
								t[20] = 5;
								t[21] = 1;
								t[22] = 2;
								t[23] = 10;
								t[24] = 6;
								t[25] = 5;
								t[26] = 2;
								t[27] = 10;
								t[28] = 3;
								t[29] = 6;
								t[30] = 2;
								t[31] = 10;
								t[32] = 1;
								t[33] = 3;
								t[34] = 2;
								t[35] = 10;
								t[36] = 3;
								t[37] = 4;
								t[38] = 6;
								t[39] = 10;
								(*nint) = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 5;
								t[1] = 3;
								t[2] = 1;
								t[3] = 2;
								t[4] = 0;
								t[5] = 5;
								t[6] = 6;
								t[7] = 4;
								t[8] = 3;
								t[9] = 4;
								t[10] = 5;
								t[11] = 1;
								t[12] = 3;
								t[13] = 4;
								t[14] = 6;
								t[15] = 5;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 5;
								t[1] = 1;
								t[2] = 9;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 9;
								t[7] = 1;
								t[8] = 6;
								t[9] = 5;
								t[10] = 4;
								t[11] = 0;
								t[12] = 1;
								t[13] = 5;
								t[14] = 9;
								t[15] = 2;
								t[16] = 3;
								t[17] = 4;
								t[18] = 6;
								t[19] = 9;
								t[20] = 6;
								t[21] = 5;
								t[22] = 9;
								t[23] = 4;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 5;
								t[2] = 2;
								t[3] = 1;
								t[4] = 8;
								t[5] = 5;
								t[6] = 1;
								t[7] = 4;
								t[8] = 6;
								t[9] = 4;
								t[10] = 0;
								t[11] = 5;
								t[12] = 8;
								t[13] = 5;
								t[14] = 6;
								t[15] = 2;
								t[16] = 3;
								t[17] = 8;
								t[18] = 6;
								t[19] = 2;
								t[20] = 6;
								t[21] = 4;
								t[22] = 5;
								t[23] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 5;
								t[2] = 2;
								t[3] = 1;
								t[4] = 8;
								t[5] = 5;
								t[6] = 1;
								t[7] = 4;
								t[8] = 6;
								t[9] = 4;
								t[10] = 0;
								t[11] = 5;
								t[12] = 6;
								t[13] = 4;
								t[14] = 5;
								t[15] = 8;
								t[16] = 3;
								t[17] = 5;
								t[18] = 2;
								t[19] = 8;
								t[20] = 3;
								t[21] = 5;
								t[22] = 8;
								t[23] = 6;
                            } else if (edges[5] == 9) {
                                (*nel) = 7;
								t[0] = 8;
								t[1] = 3;
								t[2] = 9;
								t[3] = 6;
								t[4] = 6;
								t[5] = 5;
								t[6] = 4;
								t[7] = 0;
								t[8] = 2;
								t[9] = 9;
								t[10] = 5;
								t[11] = 1;
								t[12] = 9;
								t[13] = 6;
								t[14] = 4;
								t[15] = 8;
								t[16] = 5;
								t[17] = 6;
								t[18] = 4;
								t[19] = 9;
								t[20] = 9;
								t[21] = 5;
								t[22] = 1;
								t[23] = 4;
								t[24] = 8;
								t[25] = 9;
								t[26] = 1;
								t[27] = 4;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 2;
								t[1] = 4;
								t[2] = 6;
								t[3] = 5;
								t[4] = 0;
								t[5] = 5;
								t[6] = 6;
								t[7] = 4;
								t[8] = 2;
								t[9] = 4;
								t[10] = 1;
								t[11] = 6;
								t[12] = 1;
								t[13] = 6;
								t[14] = 3;
								t[15] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 10;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 10;
								t[4] = 1;
								t[5] = 3;
								t[6] = 2;
								t[7] = 10;
								t[8] = 6;
								t[9] = 1;
								t[10] = 10;
								t[11] = 3;
								t[12] = 10;
								t[13] = 0;
								t[14] = 4;
								t[15] = 6;
								t[16] = 10;
								t[17] = 0;
								t[18] = 5;
								t[19] = 4;
								t[20] = 3;
								t[21] = 5;
								t[22] = 2;
								t[23] = 10;
								t[24] = 4;
								t[25] = 2;
								t[26] = 5;
								t[27] = 10;
								t[28] = 6;
								t[29] = 1;
								t[30] = 4;
								t[31] = 10;
								t[32] = 5;
								t[33] = 3;
								t[34] = 6;
								t[35] = 10;
								t[36] = 10;
								t[37] = 0;
								t[38] = 6;
								t[39] = 5;
								(*nint) = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 6;
								t[1] = 5;
								t[2] = 4;
								t[3] = 0;
								t[4] = 2;
								t[5] = 4;
								t[6] = 1;
								t[7] = 9;
								t[8] = 3;
								t[9] = 1;
								t[10] = 6;
								t[11] = 9;
								t[12] = 2;
								t[13] = 4;
								t[14] = 9;
								t[15] = 5;
								t[16] = 1;
								t[17] = 6;
								t[18] = 9;
								t[19] = 4;
								t[20] = 6;
								t[21] = 5;
								t[22] = 9;
								t[23] = 4;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 2;
								t[1] = 4;
								t[2] = 6;
								t[3] = 5;
								t[4] = 0;
								t[5] = 5;
								t[6] = 6;
								t[7] = 4;
								t[8] = 4;
								t[9] = 3;
								t[10] = 2;
								t[11] = 6;
								t[12] = 4;
								t[13] = 3;
								t[14] = 1;
								t[15] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 2;
								t[1] = 3;
								t[2] = 4;
								t[3] = 1;
								t[4] = 0;
								t[5] = 5;
								t[6] = 6;
								t[7] = 4;
								t[8] = 5;
								t[9] = 3;
								t[10] = 4;
								t[11] = 2;
								t[12] = 3;
								t[13] = 5;
								t[14] = 4;
								t[15] = 6;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 6;
								t[1] = 5;
								t[2] = 4;
								t[3] = 0;
								t[4] = 2;
								t[5] = 4;
								t[6] = 1;
								t[7] = 9;
								t[8] = 4;
								t[9] = 3;
								t[10] = 9;
								t[11] = 6;
								t[12] = 4;
								t[13] = 3;
								t[14] = 1;
								t[15] = 9;
								t[16] = 6;
								t[17] = 5;
								t[18] = 9;
								t[19] = 4;
								t[20] = 2;
								t[21] = 4;
								t[22] = 9;
								t[23] = 5;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 6;
								t[0] = 4;
								t[1] = 2;
								t[2] = 5;
								t[3] = 8;
								t[4] = 6;
								t[5] = 4;
								t[6] = 0;
								t[7] = 5;
								t[8] = 8;
								t[9] = 5;
								t[10] = 6;
								t[11] = 2;
								t[12] = 3;
								t[13] = 8;
								t[14] = 6;
								t[15] = 2;
								t[16] = 4;
								t[17] = 2;
								t[18] = 8;
								t[19] = 1;
								t[20] = 6;
								t[21] = 4;
								t[22] = 5;
								t[23] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 6;
								t[0] = 4;
								t[1] = 2;
								t[2] = 5;
								t[3] = 8;
								t[4] = 6;
								t[5] = 4;
								t[6] = 0;
								t[7] = 5;
								t[8] = 3;
								t[9] = 5;
								t[10] = 8;
								t[11] = 6;
								t[12] = 3;
								t[13] = 5;
								t[14] = 2;
								t[15] = 8;
								t[16] = 4;
								t[17] = 2;
								t[18] = 8;
								t[19] = 1;
								t[20] = 6;
								t[21] = 4;
								t[22] = 5;
								t[23] = 8;
                            } else if (edges[5] == 9) {
                                (*nel) = 7;
								t[0] = 8;
								t[1] = 3;
								t[2] = 9;
								t[3] = 6;
								t[4] = 6;
								t[5] = 5;
								t[6] = 4;
								t[7] = 0;
								t[8] = 8;
								t[9] = 2;
								t[10] = 4;
								t[11] = 9;
								t[12] = 9;
								t[13] = 6;
								t[14] = 4;
								t[15] = 8;
								t[16] = 5;
								t[17] = 6;
								t[18] = 4;
								t[19] = 9;
								t[20] = 4;
								t[21] = 2;
								t[22] = 5;
								t[23] = 9;
								t[24] = 8;
								t[25] = 2;
								t[26] = 1;
								t[27] = 4;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 6;
								t[0] = 7;
								t[1] = 3;
								t[2] = 6;
								t[3] = 1;
								t[4] = 6;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 4;
								t[9] = 5;
								t[10] = 0;
								t[11] = 6;
								t[12] = 7;
								t[13] = 5;
								t[14] = 6;
								t[15] = 2;
								t[16] = 3;
								t[17] = 7;
								t[18] = 6;
								t[19] = 2;
								t[20] = 4;
								t[21] = 5;
								t[22] = 6;
								t[23] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 6;
								t[0] = 7;
								t[1] = 3;
								t[2] = 6;
								t[3] = 1;
								t[4] = 6;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 4;
								t[9] = 5;
								t[10] = 0;
								t[11] = 6;
								t[12] = 4;
								t[13] = 5;
								t[14] = 6;
								t[15] = 7;
								t[16] = 3;
								t[17] = 5;
								t[18] = 2;
								t[19] = 7;
								t[20] = 3;
								t[21] = 5;
								t[22] = 7;
								t[23] = 6;
                            } else if (edges[5] == 9) {
                                (*nel) = 7;
								t[0] = 2;
								t[1] = 7;
								t[2] = 9;
								t[3] = 5;
								t[4] = 5;
								t[5] = 9;
								t[6] = 4;
								t[7] = 7;
								t[8] = 5;
								t[9] = 6;
								t[10] = 4;
								t[11] = 9;
								t[12] = 9;
								t[13] = 3;
								t[14] = 6;
								t[15] = 1;
								t[16] = 9;
								t[17] = 7;
								t[18] = 1;
								t[19] = 4;
								t[20] = 6;
								t[21] = 9;
								t[22] = 1;
								t[23] = 4;
								t[24] = 6;
								t[25] = 5;
								t[26] = 4;
								t[27] = 0;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 6;
								t[0] = 3;
								t[1] = 4;
								t[2] = 6;
								t[3] = 7;
								t[4] = 4;
								t[5] = 5;
								t[6] = 0;
								t[7] = 6;
								t[8] = 7;
								t[9] = 5;
								t[10] = 6;
								t[11] = 2;
								t[12] = 3;
								t[13] = 7;
								t[14] = 6;
								t[15] = 2;
								t[16] = 4;
								t[17] = 5;
								t[18] = 6;
								t[19] = 7;
								t[20] = 3;
								t[21] = 4;
								t[22] = 7;
								t[23] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 6;
								t[0] = 3;
								t[1] = 4;
								t[2] = 7;
								t[3] = 1;
								t[4] = 3;
								t[5] = 4;
								t[6] = 6;
								t[7] = 7;
								t[8] = 4;
								t[9] = 5;
								t[10] = 6;
								t[11] = 7;
								t[12] = 4;
								t[13] = 5;
								t[14] = 0;
								t[15] = 6;
								t[16] = 5;
								t[17] = 3;
								t[18] = 6;
								t[19] = 7;
								t[20] = 5;
								t[21] = 3;
								t[22] = 7;
								t[23] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 7;
								t[0] = 2;
								t[1] = 7;
								t[2] = 9;
								t[3] = 5;
								t[4] = 9;
								t[5] = 4;
								t[6] = 3;
								t[7] = 6;
								t[8] = 4;
								t[9] = 9;
								t[10] = 5;
								t[11] = 6;
								t[12] = 9;
								t[13] = 4;
								t[14] = 7;
								t[15] = 3;
								t[16] = 3;
								t[17] = 4;
								t[18] = 7;
								t[19] = 1;
								t[20] = 9;
								t[21] = 4;
								t[22] = 5;
								t[23] = 7;
								t[24] = 6;
								t[25] = 5;
								t[26] = 4;
								t[27] = 0;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 7;
								t[0] = 2;
								t[1] = 3;
								t[2] = 6;
								t[3] = 8;
								t[4] = 7;
								t[5] = 1;
								t[6] = 8;
								t[7] = 4;
								t[8] = 4;
								t[9] = 6;
								t[10] = 5;
								t[11] = 0;
								t[12] = 6;
								t[13] = 4;
								t[14] = 5;
								t[15] = 8;
								t[16] = 8;
								t[17] = 4;
								t[18] = 5;
								t[19] = 7;
								t[20] = 2;
								t[21] = 6;
								t[22] = 5;
								t[23] = 8;
								t[24] = 2;
								t[25] = 5;
								t[26] = 7;
								t[27] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 7;
								t[0] = 3;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 5;
								t[5] = 8;
								t[6] = 4;
								t[7] = 7;
								t[8] = 8;
								t[9] = 5;
								t[10] = 3;
								t[11] = 7;
								t[12] = 7;
								t[13] = 1;
								t[14] = 8;
								t[15] = 4;
								t[16] = 5;
								t[17] = 8;
								t[18] = 6;
								t[19] = 4;
								t[20] = 4;
								t[21] = 5;
								t[22] = 0;
								t[23] = 6;
								t[24] = 5;
								t[25] = 8;
								t[26] = 3;
								t[27] = 6;
                            } else if (edges[5] == 9) {
                                (*nel) = 8;
								t[0] = 9;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 9;
								t[5] = 3;
								t[6] = 6;
								t[7] = 8;
								t[8] = 6;
								t[9] = 4;
								t[10] = 0;
								t[11] = 5;
								t[12] = 1;
								t[13] = 4;
								t[14] = 8;
								t[15] = 7;
								t[16] = 9;
								t[17] = 4;
								t[18] = 7;
								t[19] = 8;
								t[20] = 4;
								t[21] = 9;
								t[22] = 6;
								t[23] = 8;
								t[24] = 4;
								t[25] = 9;
								t[26] = 5;
								t[27] = 6;
								t[28] = 4;
								t[29] = 9;
								t[30] = 7;
								t[31] = 5;
                            }
                        }
                    }
                }
            }
        } if (edges[0] == 4) {
            if (edges[1] == 0) {
                if (edges[2] == 0) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 2;
								t[5] = 0;
								t[6] = 4;
								t[7] = 8;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 2;
								t[5] = 0;
								t[6] = 4;
								t[7] = 8;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 8;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 1;
								t[5] = 9;
								t[6] = 4;
								t[7] = 8;
								t[8] = 2;
								t[9] = 1;
								t[10] = 9;
								t[11] = 4;
								t[12] = 3;
								t[13] = 0;
								t[14] = 9;
								t[15] = 8;
								t[16] = 9;
								t[17] = 0;
								t[18] = 4;
								t[19] = 8;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 2;
								t[5] = 0;
								t[6] = 4;
								t[7] = 8;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 2;
								t[5] = 0;
								t[6] = 4;
								t[7] = 8;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 8;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 8;
								t[9] = 2;
								t[10] = 4;
								t[11] = 9;
								t[12] = 3;
								t[13] = 0;
								t[14] = 9;
								t[15] = 8;
								t[16] = 9;
								t[17] = 0;
								t[18] = 4;
								t[19] = 8;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 7;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 7;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 1;
								t[1] = 3;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 1;
								t[6] = 4;
								t[7] = 7;
								t[8] = 3;
								t[9] = 0;
								t[10] = 9;
								t[11] = 4;
								t[12] = 0;
								t[13] = 9;
								t[14] = 4;
								t[15] = 7;
								t[16] = 0;
								t[17] = 2;
								t[18] = 9;
								t[19] = 7;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 7;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 7;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 7;
								t[2] = 1;
								t[3] = 4;
								t[4] = 3;
								t[5] = 0;
								t[6] = 9;
								t[7] = 4;
								t[8] = 3;
								t[9] = 7;
								t[10] = 4;
								t[11] = 9;
								t[12] = 0;
								t[13] = 9;
								t[14] = 4;
								t[15] = 7;
								t[16] = 0;
								t[17] = 2;
								t[18] = 9;
								t[19] = 7;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 8;
								t[4] = 8;
								t[5] = 1;
								t[6] = 4;
								t[7] = 7;
								t[8] = 0;
								t[9] = 8;
								t[10] = 4;
								t[11] = 7;
								t[12] = 0;
								t[13] = 2;
								t[14] = 8;
								t[15] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 3;
								t[1] = 7;
								t[2] = 8;
								t[3] = 0;
								t[4] = 8;
								t[5] = 1;
								t[6] = 4;
								t[7] = 7;
								t[8] = 3;
								t[9] = 7;
								t[10] = 0;
								t[11] = 2;
								t[12] = 0;
								t[13] = 8;
								t[14] = 4;
								t[15] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 7;
								t[4] = 9;
								t[5] = 0;
								t[6] = 8;
								t[7] = 3;
								t[8] = 8;
								t[9] = 0;
								t[10] = 7;
								t[11] = 4;
								t[12] = 1;
								t[13] = 8;
								t[14] = 7;
								t[15] = 4;
								t[16] = 0;
								t[17] = 9;
								t[18] = 8;
								t[19] = 7;
                            }
                        }
                    }
                } else if (edges[2] == 3) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 8;
								t[8] = 3;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 8;
								t[8] = 3;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 1;
								t[5] = 9;
								t[6] = 4;
								t[7] = 8;
								t[8] = 2;
								t[9] = 1;
								t[10] = 9;
								t[11] = 4;
								t[12] = 4;
								t[13] = 3;
								t[14] = 9;
								t[15] = 0;
								t[16] = 4;
								t[17] = 3;
								t[18] = 8;
								t[19] = 9;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 8;
								t[8] = 3;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 8;
								t[8] = 3;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 4;
								t[9] = 3;
								t[10] = 9;
								t[11] = 0;
								t[12] = 4;
								t[13] = 3;
								t[14] = 8;
								t[15] = 9;
								t[16] = 8;
								t[17] = 2;
								t[18] = 4;
								t[19] = 9;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 7;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 7;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 1;
								t[1] = 3;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 1;
								t[6] = 4;
								t[7] = 7;
								t[8] = 3;
								t[9] = 0;
								t[10] = 9;
								t[11] = 4;
								t[12] = 0;
								t[13] = 9;
								t[14] = 4;
								t[15] = 7;
								t[16] = 0;
								t[17] = 2;
								t[18] = 9;
								t[19] = 7;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 7;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 7;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 7;
								t[2] = 1;
								t[3] = 4;
								t[4] = 3;
								t[5] = 0;
								t[6] = 9;
								t[7] = 4;
								t[8] = 3;
								t[9] = 7;
								t[10] = 4;
								t[11] = 9;
								t[12] = 0;
								t[13] = 9;
								t[14] = 4;
								t[15] = 7;
								t[16] = 0;
								t[17] = 2;
								t[18] = 9;
								t[19] = 7;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 9;
								t[0] = 4;
								t[1] = 3;
								t[2] = 8;
								t[3] = 10;
								t[4] = 3;
								t[5] = 2;
								t[6] = 8;
								t[7] = 10;
								t[8] = 8;
								t[9] = 1;
								t[10] = 4;
								t[11] = 7;
								t[12] = 0;
								t[13] = 2;
								t[14] = 3;
								t[15] = 10;
								t[16] = 3;
								t[17] = 4;
								t[18] = 0;
								t[19] = 10;
								t[20] = 7;
								t[21] = 0;
								t[22] = 10;
								t[23] = 2;
								t[24] = 7;
								t[25] = 0;
								t[26] = 4;
								t[27] = 10;
								t[28] = 8;
								t[29] = 2;
								t[30] = 7;
								t[31] = 10;
								t[32] = 4;
								t[33] = 8;
								t[34] = 7;
								t[35] = 10;
								(*nint) = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 3;
								t[1] = 4;
								t[2] = 7;
								t[3] = 8;
								t[4] = 8;
								t[5] = 1;
								t[6] = 4;
								t[7] = 7;
								t[8] = 3;
								t[9] = 4;
								t[10] = 0;
								t[11] = 7;
								t[12] = 7;
								t[13] = 3;
								t[14] = 2;
								t[15] = 0;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 7;
								t[4] = 9;
								t[5] = 4;
								t[6] = 7;
								t[7] = 8;
								t[8] = 4;
								t[9] = 3;
								t[10] = 9;
								t[11] = 0;
								t[12] = 3;
								t[13] = 4;
								t[14] = 9;
								t[15] = 8;
								t[16] = 1;
								t[17] = 8;
								t[18] = 7;
								t[19] = 4;
								t[20] = 9;
								t[21] = 4;
								t[22] = 0;
								t[23] = 7;
                            }
                        }
                    }
                } else if (edges[2] == 6) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 6;
								t[5] = 2;
								t[6] = 3;
								t[7] = 1;
								t[8] = 2;
								t[9] = 6;
								t[10] = 4;
								t[11] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 6;
								t[5] = 2;
								t[6] = 3;
								t[7] = 1;
								t[8] = 2;
								t[9] = 6;
								t[10] = 4;
								t[11] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 3;
								t[5] = 1;
								t[6] = 6;
								t[7] = 9;
								t[8] = 1;
								t[9] = 6;
								t[10] = 9;
								t[11] = 4;
								t[12] = 6;
								t[13] = 0;
								t[14] = 9;
								t[15] = 4;
								t[16] = 2;
								t[17] = 1;
								t[18] = 9;
								t[19] = 4;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 1;
								t[8] = 3;
								t[9] = 4;
								t[10] = 6;
								t[11] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 1;
								t[8] = 3;
								t[9] = 4;
								t[10] = 6;
								t[11] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 6;
								t[7] = 9;
								t[8] = 3;
								t[9] = 4;
								t[10] = 9;
								t[11] = 1;
								t[12] = 6;
								t[13] = 0;
								t[14] = 9;
								t[15] = 4;
								t[16] = 2;
								t[17] = 1;
								t[18] = 9;
								t[19] = 4;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 6;
								t[1] = 4;
								t[2] = 2;
								t[3] = 8;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 6;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
								t[12] = 3;
								t[13] = 2;
								t[14] = 8;
								t[15] = 6;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 6;
								t[1] = 4;
								t[2] = 2;
								t[3] = 8;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 6;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
								t[12] = 3;
								t[13] = 2;
								t[14] = 8;
								t[15] = 6;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 3;
								t[2] = 9;
								t[3] = 6;
								t[4] = 6;
								t[5] = 4;
								t[6] = 9;
								t[7] = 8;
								t[8] = 8;
								t[9] = 9;
								t[10] = 1;
								t[11] = 4;
								t[12] = 6;
								t[13] = 4;
								t[14] = 0;
								t[15] = 9;
								t[16] = 9;
								t[17] = 2;
								t[18] = 1;
								t[19] = 4;
								t[20] = 0;
								t[21] = 2;
								t[22] = 9;
								t[23] = 4;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 6;
								t[5] = 2;
								t[6] = 3;
								t[7] = 1;
								t[8] = 2;
								t[9] = 6;
								t[10] = 4;
								t[11] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 6;
								t[5] = 2;
								t[6] = 3;
								t[7] = 1;
								t[8] = 2;
								t[9] = 6;
								t[10] = 4;
								t[11] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 3;
								t[5] = 1;
								t[6] = 6;
								t[7] = 9;
								t[8] = 1;
								t[9] = 6;
								t[10] = 9;
								t[11] = 4;
								t[12] = 6;
								t[13] = 0;
								t[14] = 9;
								t[15] = 4;
								t[16] = 2;
								t[17] = 1;
								t[18] = 9;
								t[19] = 4;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 1;
								t[8] = 3;
								t[9] = 4;
								t[10] = 6;
								t[11] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 1;
								t[8] = 3;
								t[9] = 4;
								t[10] = 6;
								t[11] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 6;
								t[7] = 9;
								t[8] = 3;
								t[9] = 4;
								t[10] = 9;
								t[11] = 1;
								t[12] = 6;
								t[13] = 0;
								t[14] = 9;
								t[15] = 4;
								t[16] = 2;
								t[17] = 1;
								t[18] = 9;
								t[19] = 4;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 6;
								t[1] = 4;
								t[2] = 2;
								t[3] = 8;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 6;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
								t[12] = 3;
								t[13] = 2;
								t[14] = 8;
								t[15] = 6;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 6;
								t[1] = 4;
								t[2] = 2;
								t[3] = 8;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 6;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
								t[12] = 3;
								t[13] = 2;
								t[14] = 8;
								t[15] = 6;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 2;
								t[2] = 4;
								t[3] = 9;
								t[4] = 8;
								t[5] = 3;
								t[6] = 9;
								t[7] = 6;
								t[8] = 6;
								t[9] = 4;
								t[10] = 9;
								t[11] = 8;
								t[12] = 8;
								t[13] = 2;
								t[14] = 1;
								t[15] = 4;
								t[16] = 6;
								t[17] = 4;
								t[18] = 0;
								t[19] = 9;
								t[20] = 0;
								t[21] = 2;
								t[22] = 9;
								t[23] = 4;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 7;
								t[1] = 3;
								t[2] = 6;
								t[3] = 1;
								t[4] = 2;
								t[5] = 0;
								t[6] = 7;
								t[7] = 6;
								t[8] = 3;
								t[9] = 2;
								t[10] = 7;
								t[11] = 6;
								t[12] = 6;
								t[13] = 7;
								t[14] = 1;
								t[15] = 4;
								t[16] = 0;
								t[17] = 7;
								t[18] = 6;
								t[19] = 4;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 7;
								t[1] = 3;
								t[2] = 6;
								t[3] = 1;
								t[4] = 2;
								t[5] = 0;
								t[6] = 7;
								t[7] = 6;
								t[8] = 3;
								t[9] = 2;
								t[10] = 7;
								t[11] = 6;
								t[12] = 6;
								t[13] = 7;
								t[14] = 1;
								t[15] = 4;
								t[16] = 0;
								t[17] = 7;
								t[18] = 6;
								t[19] = 4;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 2;
								t[1] = 0;
								t[2] = 7;
								t[3] = 9;
								t[4] = 9;
								t[5] = 3;
								t[6] = 6;
								t[7] = 1;
								t[8] = 9;
								t[9] = 7;
								t[10] = 1;
								t[11] = 4;
								t[12] = 6;
								t[13] = 9;
								t[14] = 1;
								t[15] = 4;
								t[16] = 0;
								t[17] = 7;
								t[18] = 9;
								t[19] = 4;
								t[20] = 0;
								t[21] = 9;
								t[22] = 6;
								t[23] = 4;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 4;
								t[2] = 6;
								t[3] = 7;
								t[4] = 2;
								t[5] = 0;
								t[6] = 7;
								t[7] = 6;
								t[8] = 3;
								t[9] = 2;
								t[10] = 7;
								t[11] = 6;
								t[12] = 0;
								t[13] = 7;
								t[14] = 6;
								t[15] = 4;
								t[16] = 3;
								t[17] = 4;
								t[18] = 7;
								t[19] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 4;
								t[2] = 6;
								t[3] = 7;
								t[4] = 2;
								t[5] = 0;
								t[6] = 7;
								t[7] = 6;
								t[8] = 3;
								t[9] = 2;
								t[10] = 7;
								t[11] = 6;
								t[12] = 0;
								t[13] = 7;
								t[14] = 6;
								t[15] = 4;
								t[16] = 3;
								t[17] = 4;
								t[18] = 7;
								t[19] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 9;
								t[1] = 4;
								t[2] = 0;
								t[3] = 7;
								t[4] = 2;
								t[5] = 0;
								t[6] = 7;
								t[7] = 9;
								t[8] = 9;
								t[9] = 4;
								t[10] = 3;
								t[11] = 6;
								t[12] = 4;
								t[13] = 9;
								t[14] = 0;
								t[15] = 6;
								t[16] = 9;
								t[17] = 4;
								t[18] = 7;
								t[19] = 3;
								t[20] = 3;
								t[21] = 4;
								t[22] = 7;
								t[23] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 7;
								t[2] = 1;
								t[3] = 4;
								t[4] = 2;
								t[5] = 0;
								t[6] = 7;
								t[7] = 6;
								t[8] = 6;
								t[9] = 4;
								t[10] = 7;
								t[11] = 8;
								t[12] = 6;
								t[13] = 4;
								t[14] = 0;
								t[15] = 7;
								t[16] = 2;
								t[17] = 3;
								t[18] = 6;
								t[19] = 8;
								t[20] = 7;
								t[21] = 2;
								t[22] = 6;
								t[23] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 7;
								t[2] = 1;
								t[3] = 4;
								t[4] = 7;
								t[5] = 3;
								t[6] = 6;
								t[7] = 8;
								t[8] = 2;
								t[9] = 0;
								t[10] = 7;
								t[11] = 6;
								t[12] = 6;
								t[13] = 4;
								t[14] = 7;
								t[15] = 8;
								t[16] = 6;
								t[17] = 4;
								t[18] = 0;
								t[19] = 7;
								t[20] = 7;
								t[21] = 3;
								t[22] = 2;
								t[23] = 6;
                            } else if (edges[5] == 9) {
                                (*nel) = 7;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 7;
								t[4] = 4;
								t[5] = 9;
								t[6] = 8;
								t[7] = 7;
								t[8] = 4;
								t[9] = 9;
								t[10] = 0;
								t[11] = 6;
								t[12] = 8;
								t[13] = 7;
								t[14] = 1;
								t[15] = 4;
								t[16] = 6;
								t[17] = 9;
								t[18] = 3;
								t[19] = 8;
								t[20] = 4;
								t[21] = 9;
								t[22] = 6;
								t[23] = 8;
								t[24] = 9;
								t[25] = 4;
								t[26] = 0;
								t[27] = 7;
                            }
                        }
                    }
                }
            } else if (edges[1] == 2) {
                if (edges[2] == 0) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 2;
								t[5] = 0;
								t[6] = 4;
								t[7] = 8;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 2;
								t[5] = 0;
								t[6] = 4;
								t[7] = 8;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 8;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 1;
								t[5] = 9;
								t[6] = 4;
								t[7] = 8;
								t[8] = 2;
								t[9] = 1;
								t[10] = 9;
								t[11] = 4;
								t[12] = 3;
								t[13] = 0;
								t[14] = 9;
								t[15] = 8;
								t[16] = 9;
								t[17] = 0;
								t[18] = 4;
								t[19] = 8;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 2;
								t[5] = 0;
								t[6] = 4;
								t[7] = 8;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 2;
								t[5] = 0;
								t[6] = 4;
								t[7] = 8;
								t[8] = 0;
								t[9] = 2;
								t[10] = 3;
								t[11] = 8;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 8;
								t[9] = 2;
								t[10] = 4;
								t[11] = 9;
								t[12] = 3;
								t[13] = 0;
								t[14] = 9;
								t[15] = 8;
								t[16] = 9;
								t[17] = 0;
								t[18] = 4;
								t[19] = 8;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 4;
								t[5] = 2;
								t[6] = 3;
								t[7] = 7;
								t[8] = 4;
								t[9] = 2;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 4;
								t[5] = 2;
								t[6] = 3;
								t[7] = 7;
								t[8] = 4;
								t[9] = 2;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 1;
								t[1] = 3;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 1;
								t[6] = 4;
								t[7] = 7;
								t[8] = 3;
								t[9] = 0;
								t[10] = 9;
								t[11] = 4;
								t[12] = 4;
								t[13] = 2;
								t[14] = 9;
								t[15] = 7;
								t[16] = 4;
								t[17] = 2;
								t[18] = 0;
								t[19] = 9;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 4;
								t[5] = 2;
								t[6] = 3;
								t[7] = 7;
								t[8] = 4;
								t[9] = 2;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 4;
								t[5] = 2;
								t[6] = 3;
								t[7] = 7;
								t[8] = 4;
								t[9] = 2;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 7;
								t[1] = 3;
								t[2] = 9;
								t[3] = 4;
								t[4] = 3;
								t[5] = 0;
								t[6] = 9;
								t[7] = 4;
								t[8] = 4;
								t[9] = 2;
								t[10] = 9;
								t[11] = 7;
								t[12] = 4;
								t[13] = 2;
								t[14] = 0;
								t[15] = 9;
								t[16] = 7;
								t[17] = 3;
								t[18] = 4;
								t[19] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 8;
								t[4] = 8;
								t[5] = 1;
								t[6] = 4;
								t[7] = 7;
								t[8] = 4;
								t[9] = 2;
								t[10] = 8;
								t[11] = 7;
								t[12] = 4;
								t[13] = 2;
								t[14] = 0;
								t[15] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 9;
								t[0] = 3;
								t[1] = 7;
								t[2] = 10;
								t[3] = 2;
								t[4] = 4;
								t[5] = 2;
								t[6] = 0;
								t[7] = 10;
								t[8] = 0;
								t[9] = 8;
								t[10] = 4;
								t[11] = 10;
								t[12] = 8;
								t[13] = 1;
								t[14] = 4;
								t[15] = 7;
								t[16] = 2;
								t[17] = 4;
								t[18] = 7;
								t[19] = 10;
								t[20] = 4;
								t[21] = 8;
								t[22] = 7;
								t[23] = 10;
								t[24] = 3;
								t[25] = 7;
								t[26] = 8;
								t[27] = 10;
								t[28] = 0;
								t[29] = 3;
								t[30] = 8;
								t[31] = 10;
								t[32] = 0;
								t[33] = 2;
								t[34] = 3;
								t[35] = 10;
								(*nint) = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 4;
								t[1] = 2;
								t[2] = 9;
								t[3] = 7;
								t[4] = 4;
								t[5] = 9;
								t[6] = 8;
								t[7] = 7;
								t[8] = 9;
								t[9] = 0;
								t[10] = 8;
								t[11] = 3;
								t[12] = 4;
								t[13] = 9;
								t[14] = 0;
								t[15] = 8;
								t[16] = 1;
								t[17] = 8;
								t[18] = 7;
								t[19] = 4;
								t[20] = 4;
								t[21] = 2;
								t[22] = 0;
								t[23] = 9;
                            }
                        }
                    }
                } else if (edges[2] == 3) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 8;
								t[8] = 3;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 8;
								t[8] = 3;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 1;
								t[5] = 9;
								t[6] = 4;
								t[7] = 8;
								t[8] = 2;
								t[9] = 1;
								t[10] = 9;
								t[11] = 4;
								t[12] = 4;
								t[13] = 3;
								t[14] = 9;
								t[15] = 0;
								t[16] = 4;
								t[17] = 3;
								t[18] = 8;
								t[19] = 9;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 2;
								t[0] = 0;
								t[1] = 2;
								t[2] = 3;
								t[3] = 4;
								t[4] = 2;
								t[5] = 3;
								t[6] = 4;
								t[7] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 0;
								t[6] = 4;
								t[7] = 3;
								t[8] = 4;
								t[9] = 9;
								t[10] = 3;
								t[11] = 1;
								t[12] = 2;
								t[13] = 9;
								t[14] = 4;
								t[15] = 1;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                // substitute_423242
                            } else if (edges[5] == 3) {
                                // substitute_423243
                            } else if (edges[5] == 9) {
                                // substitute_423249
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 8;
								t[8] = 3;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 8;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 8;
								t[8] = 3;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 2;
								t[2] = 9;
								t[3] = 4;
								t[4] = 2;
								t[5] = 8;
								t[6] = 9;
								t[7] = 4;
								t[8] = 3;
								t[9] = 4;
								t[10] = 9;
								t[11] = 8;
								t[12] = 2;
								t[13] = 8;
								t[14] = 4;
								t[15] = 1;
								t[16] = 3;
								t[17] = 4;
								t[18] = 0;
								t[19] = 9;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 4;
								t[5] = 2;
								t[6] = 3;
								t[7] = 7;
								t[8] = 4;
								t[9] = 2;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 4;
								t[5] = 2;
								t[6] = 3;
								t[7] = 7;
								t[8] = 4;
								t[9] = 2;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 1;
								t[1] = 3;
								t[2] = 9;
								t[3] = 4;
								t[4] = 9;
								t[5] = 1;
								t[6] = 4;
								t[7] = 7;
								t[8] = 3;
								t[9] = 0;
								t[10] = 9;
								t[11] = 4;
								t[12] = 4;
								t[13] = 2;
								t[14] = 9;
								t[15] = 7;
								t[16] = 4;
								t[17] = 2;
								t[18] = 0;
								t[19] = 9;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 4;
								t[5] = 2;
								t[6] = 3;
								t[7] = 7;
								t[8] = 4;
								t[9] = 2;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 3;
								t[2] = 1;
								t[3] = 7;
								t[4] = 4;
								t[5] = 2;
								t[6] = 3;
								t[7] = 7;
								t[8] = 4;
								t[9] = 2;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 7;
								t[1] = 3;
								t[2] = 9;
								t[3] = 4;
								t[4] = 4;
								t[5] = 2;
								t[6] = 9;
								t[7] = 7;
								t[8] = 3;
								t[9] = 0;
								t[10] = 9;
								t[11] = 4;
								t[12] = 7;
								t[13] = 3;
								t[14] = 4;
								t[15] = 1;
								t[16] = 4;
								t[17] = 2;
								t[18] = 0;
								t[19] = 9;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 3;
								t[1] = 4;
								t[2] = 2;
								t[3] = 8;
								t[4] = 4;
								t[5] = 3;
								t[6] = 2;
								t[7] = 0;
								t[8] = 8;
								t[9] = 1;
								t[10] = 4;
								t[11] = 7;
								t[12] = 8;
								t[13] = 4;
								t[14] = 2;
								t[15] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 8;
								t[1] = 1;
								t[2] = 4;
								t[3] = 7;
								t[4] = 3;
								t[5] = 7;
								t[6] = 4;
								t[7] = 2;
								t[8] = 3;
								t[9] = 7;
								t[10] = 8;
								t[11] = 4;
								t[12] = 4;
								t[13] = 3;
								t[14] = 2;
								t[15] = 0;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 2;
								t[1] = 4;
								t[2] = 9;
								t[3] = 0;
								t[4] = 4;
								t[5] = 3;
								t[6] = 8;
								t[7] = 9;
								t[8] = 4;
								t[9] = 2;
								t[10] = 9;
								t[11] = 7;
								t[12] = 4;
								t[13] = 3;
								t[14] = 9;
								t[15] = 0;
								t[16] = 1;
								t[17] = 8;
								t[18] = 7;
								t[19] = 4;
								t[20] = 9;
								t[21] = 4;
								t[22] = 7;
								t[23] = 8;
                            }
                        }
                    }
                } else if (edges[2] == 6) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 6;
								t[5] = 2;
								t[6] = 3;
								t[7] = 1;
								t[8] = 2;
								t[9] = 6;
								t[10] = 4;
								t[11] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 6;
								t[5] = 2;
								t[6] = 3;
								t[7] = 1;
								t[8] = 2;
								t[9] = 6;
								t[10] = 4;
								t[11] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 2;
								t[1] = 6;
								t[2] = 0;
								t[3] = 4;
								t[4] = 3;
								t[5] = 1;
								t[6] = 6;
								t[7] = 9;
								t[8] = 1;
								t[9] = 6;
								t[10] = 9;
								t[11] = 4;
								t[12] = 2;
								t[13] = 6;
								t[14] = 4;
								t[15] = 9;
								t[16] = 2;
								t[17] = 1;
								t[18] = 9;
								t[19] = 4;
                            } 
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 1;
								t[8] = 3;
								t[9] = 4;
								t[10] = 6;
								t[11] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 1;
								t[8] = 3;
								t[9] = 4;
								t[10] = 6;
								t[11] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 2;
								t[1] = 6;
								t[2] = 0;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 6;
								t[7] = 9;
								t[8] = 2;
								t[9] = 6;
								t[10] = 4;
								t[11] = 9;
								t[12] = 2;
								t[13] = 1;
								t[14] = 9;
								t[15] = 4;
								t[16] = 3;
								t[17] = 4;
								t[18] = 9;
								t[19] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 6;
								t[1] = 4;
								t[2] = 2;
								t[3] = 8;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 6;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
								t[12] = 3;
								t[13] = 2;
								t[14] = 8;
								t[15] = 6;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 6;
								t[1] = 4;
								t[2] = 2;
								t[3] = 8;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 6;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
								t[12] = 3;
								t[13] = 2;
								t[14] = 8;
								t[15] = 6;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 6;
								t[1] = 2;
								t[2] = 4;
								t[3] = 0;
								t[4] = 8;
								t[5] = 3;
								t[6] = 9;
								t[7] = 6;
								t[8] = 6;
								t[9] = 4;
								t[10] = 9;
								t[11] = 8;
								t[12] = 8;
								t[13] = 9;
								t[14] = 1;
								t[15] = 4;
								t[16] = 6;
								t[17] = 2;
								t[18] = 9;
								t[19] = 4;
								t[20] = 9;
								t[21] = 2;
								t[22] = 1;
								t[23] = 4;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 6;
								t[5] = 2;
								t[6] = 3;
								t[7] = 1;
								t[8] = 2;
								t[9] = 6;
								t[10] = 4;
								t[11] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 6;
								t[5] = 2;
								t[6] = 3;
								t[7] = 1;
								t[8] = 2;
								t[9] = 6;
								t[10] = 4;
								t[11] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 2;
								t[1] = 6;
								t[2] = 0;
								t[3] = 4;
								t[4] = 3;
								t[5] = 1;
								t[6] = 6;
								t[7] = 9;
								t[8] = 1;
								t[9] = 6;
								t[10] = 9;
								t[11] = 4;
								t[12] = 2;
								t[13] = 6;
								t[14] = 4;
								t[15] = 9;
								t[16] = 2;
								t[17] = 1;
								t[18] = 9;
								t[19] = 4;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 1;
								t[8] = 3;
								t[9] = 4;
								t[10] = 6;
								t[11] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 2;
								t[2] = 6;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 2;
								t[7] = 1;
								t[8] = 3;
								t[9] = 4;
								t[10] = 6;
								t[11] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 2;
								t[1] = 6;
								t[2] = 0;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 6;
								t[7] = 9;
								t[8] = 2;
								t[9] = 6;
								t[10] = 4;
								t[11] = 9;
								t[12] = 2;
								t[13] = 1;
								t[14] = 9;
								t[15] = 4;
								t[16] = 3;
								t[17] = 4;
								t[18] = 9;
								t[19] = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 6;
								t[1] = 4;
								t[2] = 2;
								t[3] = 8;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 6;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
								t[12] = 3;
								t[13] = 2;
								t[14] = 8;
								t[15] = 6;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 6;
								t[1] = 4;
								t[2] = 2;
								t[3] = 8;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 6;
								t[9] = 4;
								t[10] = 0;
								t[11] = 2;
								t[12] = 3;
								t[13] = 2;
								t[14] = 8;
								t[15] = 6;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 3;
								t[2] = 9;
								t[3] = 6;
								t[4] = 8;
								t[5] = 2;
								t[6] = 1;
								t[7] = 4;
								t[8] = 8;
								t[9] = 2;
								t[10] = 4;
								t[11] = 9;
								t[12] = 2;
								t[13] = 6;
								t[14] = 4;
								t[15] = 9;
								t[16] = 6;
								t[17] = 2;
								t[18] = 4;
								t[19] = 0;
								t[20] = 9;
								t[21] = 6;
								t[22] = 4;
								t[23] = 8;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 7;
								t[1] = 3;
								t[2] = 6;
								t[3] = 1;
								t[4] = 3;
								t[5] = 2;
								t[6] = 7;
								t[7] = 6;
								t[8] = 6;
								t[9] = 7;
								t[10] = 1;
								t[11] = 4;
								t[12] = 4;
								t[13] = 2;
								t[14] = 0;
								t[15] = 6;
								t[16] = 4;
								t[17] = 2;
								t[18] = 6;
								t[19] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 7;
								t[1] = 3;
								t[2] = 6;
								t[3] = 1;
								t[4] = 3;
								t[5] = 2;
								t[6] = 7;
								t[7] = 6;
								t[8] = 6;
								t[9] = 7;
								t[10] = 1;
								t[11] = 4;
								t[12] = 4;
								t[13] = 2;
								t[14] = 0;
								t[15] = 6;
								t[16] = 4;
								t[17] = 2;
								t[18] = 6;
								t[19] = 7;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 9;
								t[1] = 4;
								t[2] = 1;
								t[3] = 6;
								t[4] = 4;
								t[5] = 9;
								t[6] = 2;
								t[7] = 6;
								t[8] = 9;
								t[9] = 3;
								t[10] = 6;
								t[11] = 1;
								t[12] = 4;
								t[13] = 9;
								t[14] = 1;
								t[15] = 7;
								t[16] = 4;
								t[17] = 9;
								t[18] = 7;
								t[19] = 2;
								t[20] = 4;
								t[21] = 2;
								t[22] = 0;
								t[23] = 6;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 4;
								t[2] = 6;
								t[3] = 7;
								t[4] = 4;
								t[5] = 2;
								t[6] = 6;
								t[7] = 7;
								t[8] = 3;
								t[9] = 2;
								t[10] = 7;
								t[11] = 6;
								t[12] = 4;
								t[13] = 2;
								t[14] = 0;
								t[15] = 6;
								t[16] = 3;
								t[17] = 4;
								t[18] = 7;
								t[19] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 4;
								t[2] = 6;
								t[3] = 7;
								t[4] = 4;
								t[5] = 2;
								t[6] = 6;
								t[7] = 7;
								t[8] = 3;
								t[9] = 2;
								t[10] = 7;
								t[11] = 6;
								t[12] = 4;
								t[13] = 2;
								t[14] = 0;
								t[15] = 6;
								t[16] = 3;
								t[17] = 4;
								t[18] = 7;
								t[19] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 9;
								t[1] = 4;
								t[2] = 3;
								t[3] = 6;
								t[4] = 4;
								t[5] = 9;
								t[6] = 7;
								t[7] = 2;
								t[8] = 4;
								t[9] = 9;
								t[10] = 2;
								t[11] = 6;
								t[12] = 3;
								t[13] = 4;
								t[14] = 7;
								t[15] = 1;
								t[16] = 4;
								t[17] = 9;
								t[18] = 3;
								t[19] = 7;
								t[20] = 4;
								t[21] = 2;
								t[22] = 0;
								t[23] = 6;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 7;
								t[2] = 1;
								t[3] = 4;
								t[4] = 4;
								t[5] = 2;
								t[6] = 6;
								t[7] = 7;
								t[8] = 6;
								t[9] = 4;
								t[10] = 7;
								t[11] = 8;
								t[12] = 4;
								t[13] = 2;
								t[14] = 0;
								t[15] = 6;
								t[16] = 2;
								t[17] = 3;
								t[18] = 6;
								t[19] = 8;
								t[20] = 7;
								t[21] = 2;
								t[22] = 6;
								t[23] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 7;
								t[2] = 1;
								t[3] = 4;
								t[4] = 2;
								t[5] = 4;
								t[6] = 6;
								t[7] = 0;
								t[8] = 7;
								t[9] = 3;
								t[10] = 6;
								t[11] = 8;
								t[12] = 3;
								t[13] = 7;
								t[14] = 6;
								t[15] = 2;
								t[16] = 6;
								t[17] = 7;
								t[18] = 4;
								t[19] = 2;
								t[20] = 6;
								t[21] = 7;
								t[22] = 8;
								t[23] = 4;
                            } else if (edges[5] == 9) {
                                (*nel) = 7;
								t[0] = 2;
								t[1] = 9;
								t[2] = 4;
								t[3] = 7;
								t[4] = 9;
								t[5] = 4;
								t[6] = 7;
								t[7] = 8;
								t[8] = 2;
								t[9] = 6;
								t[10] = 4;
								t[11] = 9;
								t[12] = 6;
								t[13] = 2;
								t[14] = 4;
								t[15] = 0;
								t[16] = 8;
								t[17] = 7;
								t[18] = 1;
								t[19] = 4;
								t[20] = 6;
								t[21] = 9;
								t[22] = 3;
								t[23] = 8;
								t[24] = 6;
								t[25] = 4;
								t[26] = 9;
								t[27] = 8;
                            }
                        }
                    }
                }
            } else if (edges[1] == 5) {
                if (edges[2] == 0) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 3;
								t[2] = 4;
								t[3] = 5;
								t[4] = 3;
								t[5] = 4;
								t[6] = 5;
								t[7] = 1;
								t[8] = 2;
								t[9] = 3;
								t[10] = 5;
								t[11] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 3;
								t[2] = 4;
								t[3] = 5;
								t[4] = 3;
								t[5] = 4;
								t[6] = 5;
								t[7] = 1;
								t[8] = 2;
								t[9] = 3;
								t[10] = 5;
								t[11] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 9;
								t[1] = 5;
								t[2] = 1;
								t[3] = 4;
								t[4] = 5;
								t[5] = 9;
								t[6] = 1;
								t[7] = 2;
								t[8] = 0;
								t[9] = 5;
								t[10] = 9;
								t[11] = 4;
								t[12] = 9;
								t[13] = 1;
								t[14] = 3;
								t[15] = 4;
								t[16] = 0;
								t[17] = 9;
								t[18] = 3;
								t[19] = 4;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 3;
								t[2] = 4;
								t[3] = 5;
								t[4] = 3;
								t[5] = 4;
								t[6] = 5;
								t[7] = 1;
								t[8] = 2;
								t[9] = 3;
								t[10] = 5;
								t[11] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 3;
								t[2] = 4;
								t[3] = 5;
								t[4] = 3;
								t[5] = 4;
								t[6] = 5;
								t[7] = 1;
								t[8] = 2;
								t[9] = 3;
								t[10] = 5;
								t[11] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 9;
								t[1] = 5;
								t[2] = 1;
								t[3] = 4;
								t[4] = 5;
								t[5] = 9;
								t[6] = 1;
								t[7] = 2;
								t[8] = 0;
								t[9] = 5;
								t[10] = 9;
								t[11] = 4;
								t[12] = 9;
								t[13] = 1;
								t[14] = 3;
								t[15] = 4;
								t[16] = 0;
								t[17] = 9;
								t[18] = 3;
								t[19] = 4;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 8;
								t[4] = 8;
								t[5] = 5;
								t[6] = 1;
								t[7] = 4;
								t[8] = 0;
								t[9] = 5;
								t[10] = 8;
								t[11] = 4;
								t[12] = 5;
								t[13] = 3;
								t[14] = 8;
								t[15] = 2;
								t[16] = 8;
								t[17] = 5;
								t[18] = 2;
								t[19] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 8;
								t[4] = 8;
								t[5] = 5;
								t[6] = 1;
								t[7] = 4;
								t[8] = 0;
								t[9] = 5;
								t[10] = 8;
								t[11] = 4;
								t[12] = 5;
								t[13] = 3;
								t[14] = 8;
								t[15] = 2;
								t[16] = 8;
								t[17] = 5;
								t[18] = 2;
								t[19] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 0;
								t[1] = 3;
								t[2] = 8;
								t[3] = 9;
								t[4] = 9;
								t[5] = 5;
								t[6] = 2;
								t[7] = 1;
								t[8] = 8;
								t[9] = 0;
								t[10] = 9;
								t[11] = 4;
								t[12] = 0;
								t[13] = 5;
								t[14] = 9;
								t[15] = 4;
								t[16] = 9;
								t[17] = 5;
								t[18] = 1;
								t[19] = 4;
								t[20] = 8;
								t[21] = 9;
								t[22] = 1;
								t[23] = 4;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 5;
								t[8] = 4;
								t[9] = 2;
								t[10] = 5;
								t[11] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 3;
								t[3] = 5;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 5;
								t[8] = 2;
								t[9] = 4;
								t[10] = 1;
								t[11] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 2;
								t[1] = 4;
								t[2] = 9;
								t[3] = 5;
								t[4] = 2;
								t[5] = 4;
								t[6] = 1;
								t[7] = 9;
								t[8] = 0;
								t[9] = 5;
								t[10] = 9;
								t[11] = 4;
								t[12] = 9;
								t[13] = 1;
								t[14] = 3;
								t[15] = 4;
								t[16] = 0;
								t[17] = 9;
								t[18] = 3;
								t[19] = 4;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 5;
								t[8] = 4;
								t[9] = 2;
								t[10] = 5;
								t[11] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 5;
								t[8] = 4;
								t[9] = 2;
								t[10] = 5;
								t[11] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 2;
								t[1] = 4;
								t[2] = 9;
								t[3] = 5;
								t[4] = 2;
								t[5] = 4;
								t[6] = 1;
								t[7] = 9;
								t[8] = 0;
								t[9] = 5;
								t[10] = 9;
								t[11] = 4;
								t[12] = 9;
								t[13] = 1;
								t[14] = 3;
								t[15] = 4;
								t[16] = 0;
								t[17] = 9;
								t[18] = 3;
								t[19] = 4;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 8;
								t[4] = 4;
								t[5] = 2;
								t[6] = 5;
								t[7] = 8;
								t[8] = 0;
								t[9] = 5;
								t[10] = 8;
								t[11] = 4;
								t[12] = 5;
								t[13] = 3;
								t[14] = 8;
								t[15] = 2;
								t[16] = 4;
								t[17] = 2;
								t[18] = 8;
								t[19] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 8;
								t[4] = 4;
								t[5] = 2;
								t[6] = 5;
								t[7] = 8;
								t[8] = 0;
								t[9] = 5;
								t[10] = 8;
								t[11] = 4;
								t[12] = 5;
								t[13] = 3;
								t[14] = 8;
								t[15] = 2;
								t[16] = 4;
								t[17] = 2;
								t[18] = 8;
								t[19] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 0;
								t[1] = 3;
								t[2] = 8;
								t[3] = 9;
								t[4] = 2;
								t[5] = 8;
								t[6] = 4;
								t[7] = 1;
								t[8] = 8;
								t[9] = 0;
								t[10] = 9;
								t[11] = 4;
								t[12] = 0;
								t[13] = 5;
								t[14] = 9;
								t[15] = 4;
								t[16] = 9;
								t[17] = 2;
								t[18] = 4;
								t[19] = 5;
								t[20] = 8;
								t[21] = 2;
								t[22] = 4;
								t[23] = 9;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 4;
								t[1] = 5;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 5;
								t[9] = 3;
								t[10] = 7;
								t[11] = 2;
								t[12] = 4;
								t[13] = 5;
								t[14] = 0;
								t[15] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 4;
								t[1] = 5;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 5;
								t[9] = 3;
								t[10] = 7;
								t[11] = 2;
								t[12] = 4;
								t[13] = 5;
								t[14] = 0;
								t[15] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 3;
								t[1] = 0;
								t[2] = 9;
								t[3] = 4;
								t[4] = 3;
								t[5] = 9;
								t[6] = 1;
								t[7] = 4;
								t[8] = 7;
								t[9] = 5;
								t[10] = 9;
								t[11] = 2;
								t[12] = 4;
								t[13] = 5;
								t[14] = 9;
								t[15] = 7;
								t[16] = 9;
								t[17] = 7;
								t[18] = 1;
								t[19] = 4;
								t[20] = 4;
								t[21] = 5;
								t[22] = 0;
								t[23] = 9;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 4;
								t[1] = 5;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 5;
								t[9] = 3;
								t[10] = 7;
								t[11] = 2;
								t[12] = 4;
								t[13] = 5;
								t[14] = 0;
								t[15] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 4;
								t[1] = 5;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 5;
								t[9] = 3;
								t[10] = 7;
								t[11] = 2;
								t[12] = 4;
								t[13] = 5;
								t[14] = 0;
								t[15] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 3;
								t[1] = 0;
								t[2] = 9;
								t[3] = 4;
								t[4] = 3;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 7;
								t[9] = 5;
								t[10] = 9;
								t[11] = 2;
								t[12] = 4;
								t[13] = 5;
								t[14] = 9;
								t[15] = 7;
								t[16] = 3;
								t[17] = 7;
								t[18] = 4;
								t[19] = 9;
								t[20] = 4;
								t[21] = 5;
								t[22] = 0;
								t[23] = 9;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 7;
								t[2] = 1;
								t[3] = 4;
								t[4] = 4;
								t[5] = 5;
								t[6] = 8;
								t[7] = 7;
								t[8] = 0;
								t[9] = 5;
								t[10] = 3;
								t[11] = 8;
								t[12] = 5;
								t[13] = 7;
								t[14] = 2;
								t[15] = 8;
								t[16] = 3;
								t[17] = 5;
								t[18] = 2;
								t[19] = 8;
								t[20] = 4;
								t[21] = 5;
								t[22] = 0;
								t[23] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 7;
								t[2] = 1;
								t[3] = 4;
								t[4] = 4;
								t[5] = 5;
								t[6] = 8;
								t[7] = 7;
								t[8] = 0;
								t[9] = 5;
								t[10] = 3;
								t[11] = 8;
								t[12] = 3;
								t[13] = 7;
								t[14] = 8;
								t[15] = 5;
								t[16] = 4;
								t[17] = 5;
								t[18] = 0;
								t[19] = 8;
								t[20] = 3;
								t[21] = 7;
								t[22] = 5;
								t[23] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 7;
								t[0] = 2;
								t[1] = 7;
								t[2] = 9;
								t[3] = 5;
								t[4] = 8;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 4;
								t[9] = 9;
								t[10] = 0;
								t[11] = 8;
								t[12] = 9;
								t[13] = 4;
								t[14] = 7;
								t[15] = 8;
								t[16] = 3;
								t[17] = 0;
								t[18] = 9;
								t[19] = 8;
								t[20] = 4;
								t[21] = 9;
								t[22] = 7;
								t[23] = 5;
								t[24] = 4;
								t[25] = 9;
								t[26] = 5;
								t[27] = 0;
                            }
                        }
                    }
                } else if (edges[2] == 3) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 3;
								t[2] = 4;
								t[3] = 5;
								t[4] = 3;
								t[5] = 4;
								t[6] = 5;
								t[7] = 1;
								t[8] = 2;
								t[9] = 3;
								t[10] = 5;
								t[11] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 3;
								t[2] = 4;
								t[3] = 5;
								t[4] = 3;
								t[5] = 4;
								t[6] = 5;
								t[7] = 1;
								t[8] = 2;
								t[9] = 3;
								t[10] = 5;
								t[11] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 9;
								t[1] = 5;
								t[2] = 1;
								t[3] = 4;
								t[4] = 5;
								t[5] = 3;
								t[6] = 4;
								t[7] = 9;
								t[8] = 5;
								t[9] = 9;
								t[10] = 1;
								t[11] = 2;
								t[12] = 5;
								t[13] = 3;
								t[14] = 0;
								t[15] = 4;
								t[16] = 9;
								t[17] = 1;
								t[18] = 3;
								t[19] = 4;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 3;
								t[2] = 4;
								t[3] = 5;
								t[4] = 3;
								t[5] = 4;
								t[6] = 5;
								t[7] = 1;
								t[8] = 2;
								t[9] = 3;
								t[10] = 5;
								t[11] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 0;
								t[1] = 3;
								t[2] = 4;
								t[3] = 5;
								t[4] = 3;
								t[5] = 4;
								t[6] = 5;
								t[7] = 1;
								t[8] = 2;
								t[9] = 3;
								t[10] = 5;
								t[11] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 9;
								t[1] = 5;
								t[2] = 1;
								t[3] = 4;
								t[4] = 5;
								t[5] = 3;
								t[6] = 4;
								t[7] = 9;
								t[8] = 5;
								t[9] = 9;
								t[10] = 1;
								t[11] = 2;
								t[12] = 5;
								t[13] = 3;
								t[14] = 0;
								t[15] = 4;
								t[16] = 9;
								t[17] = 1;
								t[18] = 3;
								t[19] = 4;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 4;
								t[2] = 0;
								t[3] = 5;
								t[4] = 8;
								t[5] = 5;
								t[6] = 1;
								t[7] = 4;
								t[8] = 5;
								t[9] = 3;
								t[10] = 8;
								t[11] = 2;
								t[12] = 8;
								t[13] = 5;
								t[14] = 2;
								t[15] = 1;
								t[16] = 3;
								t[17] = 4;
								t[18] = 5;
								t[19] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 4;
								t[2] = 0;
								t[3] = 5;
								t[4] = 8;
								t[5] = 5;
								t[6] = 1;
								t[7] = 4;
								t[8] = 5;
								t[9] = 3;
								t[10] = 8;
								t[11] = 2;
								t[12] = 8;
								t[13] = 5;
								t[14] = 2;
								t[15] = 1;
								t[16] = 3;
								t[17] = 4;
								t[18] = 5;
								t[19] = 8;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 3;
								t[1] = 5;
								t[2] = 4;
								t[3] = 0;
								t[4] = 9;
								t[5] = 5;
								t[6] = 2;
								t[7] = 1;
								t[8] = 9;
								t[9] = 3;
								t[10] = 4;
								t[11] = 8;
								t[12] = 5;
								t[13] = 3;
								t[14] = 4;
								t[15] = 9;
								t[16] = 9;
								t[17] = 5;
								t[18] = 1;
								t[19] = 4;
								t[20] = 8;
								t[21] = 9;
								t[22] = 1;
								t[23] = 4;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 3;
								t[3] = 5;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 5;
								t[8] = 2;
								t[9] = 4;
								t[10] = 1;
								t[11] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 2;
								t[1] = 4;
								t[2] = 3;
								t[3] = 5;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 5;
								t[8] = 2;
								t[9] = 4;
								t[10] = 1;
								t[11] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 2;
								t[1] = 4;
								t[2] = 9;
								t[3] = 5;
								t[4] = 3;
								t[5] = 5;
								t[6] = 4;
								t[7] = 0;
								t[8] = 2;
								t[9] = 4;
								t[10] = 1;
								t[11] = 9;
								t[12] = 9;
								t[13] = 1;
								t[14] = 3;
								t[15] = 4;
								t[16] = 3;
								t[17] = 5;
								t[18] = 9;
								t[19] = 4;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 5;
								t[8] = 4;
								t[9] = 2;
								t[10] = 5;
								t[11] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 3;
								t[0] = 4;
								t[1] = 2;
								t[2] = 3;
								t[3] = 1;
								t[4] = 0;
								t[5] = 3;
								t[6] = 4;
								t[7] = 5;
								t[8] = 4;
								t[9] = 2;
								t[10] = 5;
								t[11] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 5;
								t[0] = 4;
								t[1] = 2;
								t[2] = 5;
								t[3] = 9;
								t[4] = 5;
								t[5] = 3;
								t[6] = 4;
								t[7] = 9;
								t[8] = 4;
								t[9] = 2;
								t[10] = 9;
								t[11] = 1;
								t[12] = 5;
								t[13] = 3;
								t[14] = 0;
								t[15] = 4;
								t[16] = 9;
								t[17] = 1;
								t[18] = 3;
								t[19] = 4;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 4;
								t[2] = 0;
								t[3] = 5;
								t[4] = 4;
								t[5] = 2;
								t[6] = 8;
								t[7] = 1;
								t[8] = 4;
								t[9] = 2;
								t[10] = 5;
								t[11] = 8;
								t[12] = 5;
								t[13] = 3;
								t[14] = 8;
								t[15] = 2;
								t[16] = 3;
								t[17] = 4;
								t[18] = 5;
								t[19] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 5;
								t[0] = 3;
								t[1] = 4;
								t[2] = 0;
								t[3] = 5;
								t[4] = 4;
								t[5] = 2;
								t[6] = 8;
								t[7] = 1;
								t[8] = 4;
								t[9] = 2;
								t[10] = 5;
								t[11] = 8;
								t[12] = 5;
								t[13] = 3;
								t[14] = 8;
								t[15] = 2;
								t[16] = 3;
								t[17] = 4;
								t[18] = 5;
								t[19] = 8;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 2;
								t[2] = 4;
								t[3] = 9;
								t[4] = 3;
								t[5] = 5;
								t[6] = 4;
								t[7] = 0;
								t[8] = 9;
								t[9] = 3;
								t[10] = 4;
								t[11] = 8;
								t[12] = 5;
								t[13] = 3;
								t[14] = 4;
								t[15] = 9;
								t[16] = 4;
								t[17] = 2;
								t[18] = 5;
								t[19] = 9;
								t[20] = 8;
								t[21] = 2;
								t[22] = 1;
								t[23] = 4;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 4;
								t[1] = 5;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 5;
								t[9] = 3;
								t[10] = 7;
								t[11] = 2;
								t[12] = 4;
								t[13] = 5;
								t[14] = 0;
								t[15] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 4;
								t[1] = 5;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 5;
								t[9] = 3;
								t[10] = 7;
								t[11] = 2;
								t[12] = 4;
								t[13] = 5;
								t[14] = 0;
								t[15] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 3;
								t[1] = 5;
								t[2] = 4;
								t[3] = 0;
								t[4] = 5;
								t[5] = 3;
								t[6] = 4;
								t[7] = 9;
								t[8] = 3;
								t[9] = 9;
								t[10] = 1;
								t[11] = 4;
								t[12] = 5;
								t[13] = 9;
								t[14] = 4;
								t[15] = 7;
								t[16] = 7;
								t[17] = 5;
								t[18] = 9;
								t[19] = 2;
								t[20] = 9;
								t[21] = 7;
								t[22] = 1;
								t[23] = 4;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 4;
								t[1] = 5;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 5;
								t[9] = 3;
								t[10] = 7;
								t[11] = 2;
								t[12] = 4;
								t[13] = 5;
								t[14] = 0;
								t[15] = 3;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 4;
								t[1] = 5;
								t[2] = 3;
								t[3] = 7;
								t[4] = 3;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 5;
								t[9] = 3;
								t[10] = 7;
								t[11] = 2;
								t[12] = 4;
								t[13] = 5;
								t[14] = 0;
								t[15] = 3;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 3;
								t[1] = 7;
								t[2] = 4;
								t[3] = 9;
								t[4] = 3;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 9;
								t[9] = 4;
								t[10] = 5;
								t[11] = 7;
								t[12] = 7;
								t[13] = 5;
								t[14] = 9;
								t[15] = 2;
								t[16] = 5;
								t[17] = 3;
								t[18] = 4;
								t[19] = 9;
								t[20] = 5;
								t[21] = 3;
								t[22] = 0;
								t[23] = 4;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 7;
								t[2] = 1;
								t[3] = 4;
								t[4] = 4;
								t[5] = 3;
								t[6] = 5;
								t[7] = 0;
								t[8] = 8;
								t[9] = 4;
								t[10] = 5;
								t[11] = 7;
								t[12] = 3;
								t[13] = 4;
								t[14] = 5;
								t[15] = 8;
								t[16] = 5;
								t[17] = 7;
								t[18] = 2;
								t[19] = 8;
								t[20] = 3;
								t[21] = 5;
								t[22] = 2;
								t[23] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 7;
								t[2] = 1;
								t[3] = 4;
								t[4] = 3;
								t[5] = 7;
								t[6] = 5;
								t[7] = 2;
								t[8] = 4;
								t[9] = 3;
								t[10] = 5;
								t[11] = 0;
								t[12] = 8;
								t[13] = 4;
								t[14] = 5;
								t[15] = 7;
								t[16] = 3;
								t[17] = 4;
								t[18] = 5;
								t[19] = 8;
								t[20] = 3;
								t[21] = 7;
								t[22] = 8;
								t[23] = 5;
                            } else if (edges[5] == 9) {
                                (*nel) = 7;
								t[0] = 2;
								t[1] = 7;
								t[2] = 9;
								t[3] = 5;
								t[4] = 8;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 4;
								t[9] = 9;
								t[10] = 7;
								t[11] = 5;
								t[12] = 9;
								t[13] = 4;
								t[14] = 7;
								t[15] = 8;
								t[16] = 4;
								t[17] = 3;
								t[18] = 5;
								t[19] = 0;
								t[20] = 4;
								t[21] = 9;
								t[22] = 5;
								t[23] = 3;
								t[24] = 4;
								t[25] = 9;
								t[26] = 3;
								t[27] = 8;
                            }
                        }
                    }
                } else if (edges[2] == 6) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 0;
								t[1] = 5;
								t[2] = 6;
								t[3] = 4;
								t[4] = 6;
								t[5] = 5;
								t[6] = 1;
								t[7] = 4;
								t[8] = 5;
								t[9] = 6;
								t[10] = 1;
								t[11] = 2;
								t[12] = 1;
								t[13] = 6;
								t[14] = 3;
								t[15] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 3;
								t[1] = 5;
								t[2] = 1;
								t[3] = 6;
								t[4] = 0;
								t[5] = 5;
								t[6] = 6;
								t[7] = 4;
								t[8] = 6;
								t[9] = 5;
								t[10] = 1;
								t[11] = 4;
								t[12] = 3;
								t[13] = 5;
								t[14] = 2;
								t[15] = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 5;
								t[1] = 6;
								t[2] = 0;
								t[3] = 4;
								t[4] = 5;
								t[5] = 1;
								t[6] = 9;
								t[7] = 4;
								t[8] = 5;
								t[9] = 6;
								t[10] = 4;
								t[11] = 9;
								t[12] = 1;
								t[13] = 5;
								t[14] = 9;
								t[15] = 2;
								t[16] = 3;
								t[17] = 1;
								t[18] = 6;
								t[19] = 9;
								t[20] = 1;
								t[21] = 6;
								t[22] = 9;
								t[23] = 4;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 10;
								t[0] = 0;
								t[1] = 10;
								t[2] = 6;
								t[3] = 4;
								t[4] = 1;
								t[5] = 5;
								t[6] = 4;
								t[7] = 10;
								t[8] = 0;
								t[9] = 10;
								t[10] = 5;
								t[11] = 6;
								t[12] = 0;
								t[13] = 10;
								t[14] = 4;
								t[15] = 5;
								t[16] = 3;
								t[17] = 4;
								t[18] = 10;
								t[19] = 1;
								t[20] = 5;
								t[21] = 1;
								t[22] = 2;
								t[23] = 10;
								t[24] = 6;
								t[25] = 5;
								t[26] = 2;
								t[27] = 10;
								t[28] = 3;
								t[29] = 6;
								t[30] = 2;
								t[31] = 10;
								t[32] = 1;
								t[33] = 3;
								t[34] = 2;
								t[35] = 10;
								t[36] = 3;
								t[37] = 4;
								t[38] = 6;
								t[39] = 10;
								(*nint) = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 5;
								t[1] = 3;
								t[2] = 1;
								t[3] = 2;
								t[4] = 0;
								t[5] = 5;
								t[6] = 6;
								t[7] = 4;
								t[8] = 3;
								t[9] = 4;
								t[10] = 5;
								t[11] = 1;
								t[12] = 3;
								t[13] = 4;
								t[14] = 6;
								t[15] = 5;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 5;
								t[1] = 1;
								t[2] = 9;
								t[3] = 4;
								t[4] = 3;
								t[5] = 4;
								t[6] = 9;
								t[7] = 1;
								t[8] = 6;
								t[9] = 5;
								t[10] = 4;
								t[11] = 0;
								t[12] = 1;
								t[13] = 5;
								t[14] = 9;
								t[15] = 2;
								t[16] = 3;
								t[17] = 4;
								t[18] = 6;
								t[19] = 9;
								t[20] = 6;
								t[21] = 5;
								t[22] = 9;
								t[23] = 4;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 5;
								t[2] = 2;
								t[3] = 1;
								t[4] = 8;
								t[5] = 5;
								t[6] = 1;
								t[7] = 4;
								t[8] = 6;
								t[9] = 4;
								t[10] = 0;
								t[11] = 5;
								t[12] = 8;
								t[13] = 5;
								t[14] = 6;
								t[15] = 2;
								t[16] = 3;
								t[17] = 8;
								t[18] = 6;
								t[19] = 2;
								t[20] = 6;
								t[21] = 4;
								t[22] = 5;
								t[23] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 6;
								t[0] = 8;
								t[1] = 5;
								t[2] = 2;
								t[3] = 1;
								t[4] = 8;
								t[5] = 5;
								t[6] = 1;
								t[7] = 4;
								t[8] = 6;
								t[9] = 4;
								t[10] = 0;
								t[11] = 5;
								t[12] = 6;
								t[13] = 4;
								t[14] = 5;
								t[15] = 8;
								t[16] = 3;
								t[17] = 5;
								t[18] = 2;
								t[19] = 8;
								t[20] = 3;
								t[21] = 5;
								t[22] = 8;
								t[23] = 6;
                            } else if (edges[5] == 9) {
                                (*nel) = 7;
								t[0] = 8;
								t[1] = 3;
								t[2] = 9;
								t[3] = 6;
								t[4] = 6;
								t[5] = 5;
								t[6] = 4;
								t[7] = 0;
								t[8] = 2;
								t[9] = 9;
								t[10] = 5;
								t[11] = 1;
								t[12] = 9;
								t[13] = 6;
								t[14] = 4;
								t[15] = 8;
								t[16] = 5;
								t[17] = 6;
								t[18] = 4;
								t[19] = 9;
								t[20] = 9;
								t[21] = 5;
								t[22] = 1;
								t[23] = 4;
								t[24] = 8;
								t[25] = 9;
								t[26] = 1;
								t[27] = 4;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 2;
								t[1] = 4;
								t[2] = 6;
								t[3] = 5;
								t[4] = 0;
								t[5] = 5;
								t[6] = 6;
								t[7] = 4;
								t[8] = 2;
								t[9] = 4;
								t[10] = 1;
								t[11] = 6;
								t[12] = 1;
								t[13] = 6;
								t[14] = 3;
								t[15] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 10;
								t[0] = 2;
								t[1] = 4;
								t[2] = 1;
								t[3] = 10;
								t[4] = 1;
								t[5] = 3;
								t[6] = 2;
								t[7] = 10;
								t[8] = 6;
								t[9] = 1;
								t[10] = 10;
								t[11] = 3;
								t[12] = 10;
								t[13] = 0;
								t[14] = 4;
								t[15] = 6;
								t[16] = 10;
								t[17] = 0;
								t[18] = 5;
								t[19] = 4;
								t[20] = 3;
								t[21] = 5;
								t[22] = 2;
								t[23] = 10;
								t[24] = 4;
								t[25] = 2;
								t[26] = 5;
								t[27] = 10;
								t[28] = 6;
								t[29] = 1;
								t[30] = 4;
								t[31] = 10;
								t[32] = 5;
								t[33] = 3;
								t[34] = 6;
								t[35] = 10;
								t[36] = 10;
								t[37] = 0;
								t[38] = 6;
								t[39] = 5;
								(*nint) = 1;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 6;
								t[1] = 5;
								t[2] = 4;
								t[3] = 0;
								t[4] = 2;
								t[5] = 4;
								t[6] = 1;
								t[7] = 9;
								t[8] = 3;
								t[9] = 1;
								t[10] = 6;
								t[11] = 9;
								t[12] = 2;
								t[13] = 4;
								t[14] = 9;
								t[15] = 5;
								t[16] = 1;
								t[17] = 6;
								t[18] = 9;
								t[19] = 4;
								t[20] = 6;
								t[21] = 5;
								t[22] = 9;
								t[23] = 4;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 4;
								t[0] = 2;
								t[1] = 4;
								t[2] = 6;
								t[3] = 5;
								t[4] = 0;
								t[5] = 5;
								t[6] = 6;
								t[7] = 4;
								t[8] = 4;
								t[9] = 3;
								t[10] = 2;
								t[11] = 6;
								t[12] = 4;
								t[13] = 3;
								t[14] = 1;
								t[15] = 2;
                            } else if (edges[5] == 3) {
                                (*nel) = 4;
								t[0] = 2;
								t[1] = 3;
								t[2] = 4;
								t[3] = 1;
								t[4] = 0;
								t[5] = 5;
								t[6] = 6;
								t[7] = 4;
								t[8] = 5;
								t[9] = 3;
								t[10] = 4;
								t[11] = 2;
								t[12] = 3;
								t[13] = 5;
								t[14] = 4;
								t[15] = 6;
                            } else if (edges[5] == 9) {
                                (*nel) = 6;
								t[0] = 6;
								t[1] = 5;
								t[2] = 4;
								t[3] = 0;
								t[4] = 2;
								t[5] = 4;
								t[6] = 1;
								t[7] = 9;
								t[8] = 4;
								t[9] = 3;
								t[10] = 9;
								t[11] = 6;
								t[12] = 4;
								t[13] = 3;
								t[14] = 1;
								t[15] = 9;
								t[16] = 6;
								t[17] = 5;
								t[18] = 9;
								t[19] = 4;
								t[20] = 2;
								t[21] = 4;
								t[22] = 9;
								t[23] = 5;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 6;
								t[0] = 4;
								t[1] = 2;
								t[2] = 5;
								t[3] = 8;
								t[4] = 6;
								t[5] = 4;
								t[6] = 0;
								t[7] = 5;
								t[8] = 8;
								t[9] = 5;
								t[10] = 6;
								t[11] = 2;
								t[12] = 3;
								t[13] = 8;
								t[14] = 6;
								t[15] = 2;
								t[16] = 4;
								t[17] = 2;
								t[18] = 8;
								t[19] = 1;
								t[20] = 6;
								t[21] = 4;
								t[22] = 5;
								t[23] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 6;
								t[0] = 4;
								t[1] = 2;
								t[2] = 5;
								t[3] = 8;
								t[4] = 6;
								t[5] = 4;
								t[6] = 0;
								t[7] = 5;
								t[8] = 3;
								t[9] = 5;
								t[10] = 8;
								t[11] = 6;
								t[12] = 3;
								t[13] = 5;
								t[14] = 2;
								t[15] = 8;
								t[16] = 4;
								t[17] = 2;
								t[18] = 8;
								t[19] = 1;
								t[20] = 6;
								t[21] = 4;
								t[22] = 5;
								t[23] = 8;
                            } else if (edges[5] == 9) {
                                (*nel) = 7;
								t[0] = 8;
								t[1] = 3;
								t[2] = 9;
								t[3] = 6;
								t[4] = 6;
								t[5] = 5;
								t[6] = 4;
								t[7] = 0;
								t[8] = 8;
								t[9] = 2;
								t[10] = 4;
								t[11] = 9;
								t[12] = 9;
								t[13] = 6;
								t[14] = 4;
								t[15] = 8;
								t[16] = 5;
								t[17] = 6;
								t[18] = 4;
								t[19] = 9;
								t[20] = 4;
								t[21] = 2;
								t[22] = 5;
								t[23] = 9;
								t[24] = 8;
								t[25] = 2;
								t[26] = 1;
								t[27] = 4;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                (*nel) = 6;
								t[0] = 7;
								t[1] = 3;
								t[2] = 6;
								t[3] = 1;
								t[4] = 6;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 4;
								t[9] = 5;
								t[10] = 0;
								t[11] = 6;
								t[12] = 7;
								t[13] = 5;
								t[14] = 6;
								t[15] = 2;
								t[16] = 3;
								t[17] = 7;
								t[18] = 6;
								t[19] = 2;
								t[20] = 4;
								t[21] = 5;
								t[22] = 6;
								t[23] = 7;
                            } else if (edges[5] == 3) {
                                (*nel) = 6;
								t[0] = 7;
								t[1] = 3;
								t[2] = 6;
								t[3] = 1;
								t[4] = 6;
								t[5] = 7;
								t[6] = 1;
								t[7] = 4;
								t[8] = 4;
								t[9] = 5;
								t[10] = 0;
								t[11] = 6;
								t[12] = 4;
								t[13] = 5;
								t[14] = 6;
								t[15] = 7;
								t[16] = 3;
								t[17] = 5;
								t[18] = 2;
								t[19] = 7;
								t[20] = 3;
								t[21] = 5;
								t[22] = 7;
								t[23] = 6;
                            } else if (edges[5] == 9) {
                                (*nel) = 7;
								t[0] = 2;
								t[1] = 7;
								t[2] = 9;
								t[3] = 5;
								t[4] = 5;
								t[5] = 9;
								t[6] = 4;
								t[7] = 7;
								t[8] = 5;
								t[9] = 6;
								t[10] = 4;
								t[11] = 9;
								t[12] = 9;
								t[13] = 3;
								t[14] = 6;
								t[15] = 1;
								t[16] = 9;
								t[17] = 7;
								t[18] = 1;
								t[19] = 4;
								t[20] = 6;
								t[21] = 9;
								t[22] = 1;
								t[23] = 4;
								t[24] = 6;
								t[25] = 5;
								t[26] = 4;
								t[27] = 0;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                (*nel) = 6;
								t[0] = 3;
								t[1] = 4;
								t[2] = 6;
								t[3] = 7;
								t[4] = 4;
								t[5] = 5;
								t[6] = 0;
								t[7] = 6;
								t[8] = 7;
								t[9] = 5;
								t[10] = 6;
								t[11] = 2;
								t[12] = 3;
								t[13] = 7;
								t[14] = 6;
								t[15] = 2;
								t[16] = 4;
								t[17] = 5;
								t[18] = 6;
								t[19] = 7;
								t[20] = 3;
								t[21] = 4;
								t[22] = 7;
								t[23] = 1;
                            } else if (edges[5] == 3) {
                                (*nel) = 6;
								t[0] = 3;
								t[1] = 4;
								t[2] = 7;
								t[3] = 1;
								t[4] = 3;
								t[5] = 4;
								t[6] = 6;
								t[7] = 7;
								t[8] = 4;
								t[9] = 5;
								t[10] = 6;
								t[11] = 7;
								t[12] = 4;
								t[13] = 5;
								t[14] = 0;
								t[15] = 6;
								t[16] = 5;
								t[17] = 3;
								t[18] = 6;
								t[19] = 7;
								t[20] = 5;
								t[21] = 3;
								t[22] = 7;
								t[23] = 2;
                            } else if (edges[5] == 9) {
                                (*nel) = 7;
								t[0] = 2;
								t[1] = 7;
								t[2] = 9;
								t[3] = 5;
								t[4] = 9;
								t[5] = 4;
								t[6] = 3;
								t[7] = 6;
								t[8] = 4;
								t[9] = 9;
								t[10] = 5;
								t[11] = 6;
								t[12] = 9;
								t[13] = 4;
								t[14] = 7;
								t[15] = 3;
								t[16] = 3;
								t[17] = 4;
								t[18] = 7;
								t[19] = 1;
								t[20] = 9;
								t[21] = 4;
								t[22] = 5;
								t[23] = 7;
								t[24] = 6;
								t[25] = 5;
								t[26] = 4;
								t[27] = 0;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                (*nel) = 7;
								t[0] = 2;
								t[1] = 3;
								t[2] = 6;
								t[3] = 8;
								t[4] = 7;
								t[5] = 1;
								t[6] = 8;
								t[7] = 4;
								t[8] = 4;
								t[9] = 6;
								t[10] = 5;
								t[11] = 0;
								t[12] = 6;
								t[13] = 4;
								t[14] = 5;
								t[15] = 8;
								t[16] = 8;
								t[17] = 4;
								t[18] = 5;
								t[19] = 7;
								t[20] = 2;
								t[21] = 6;
								t[22] = 5;
								t[23] = 8;
								t[24] = 2;
								t[25] = 5;
								t[26] = 7;
								t[27] = 8;
                            } else if (edges[5] == 3) {
                                (*nel) = 7;
								t[0] = 3;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 5;
								t[5] = 8;
								t[6] = 4;
								t[7] = 7;
								t[8] = 8;
								t[9] = 5;
								t[10] = 3;
								t[11] = 7;
								t[12] = 7;
								t[13] = 1;
								t[14] = 8;
								t[15] = 4;
								t[16] = 5;
								t[17] = 8;
								t[18] = 6;
								t[19] = 4;
								t[20] = 4;
								t[21] = 5;
								t[22] = 0;
								t[23] = 6;
								t[24] = 5;
								t[25] = 8;
								t[26] = 3;
								t[27] = 6;
                            } else if (edges[5] == 9) {
                                (*nel) = 8;
								t[0] = 9;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 9;
								t[5] = 3;
								t[6] = 6;
								t[7] = 8;
								t[8] = 6;
								t[9] = 4;
								t[10] = 0;
								t[11] = 5;
								t[12] = 1;
								t[13] = 4;
								t[14] = 8;
								t[15] = 7;
								t[16] = 9;
								t[17] = 4;
								t[18] = 7;
								t[19] = 8;
								t[20] = 4;
								t[21] = 9;
								t[22] = 6;
								t[23] = 8;
								t[24] = 4;
								t[25] = 9;
								t[26] = 5;
								t[27] = 6;
								t[28] = 4;
								t[29] = 9;
								t[30] = 7;
								t[31] = 5;
                            }
                        }
                    }
                }
            }
        }
        return 1; // True
    }
};

}

#endif /* KRATOS_SPLIT_TETRAHEDRA  defined */
