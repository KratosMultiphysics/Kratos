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
                                // 000112
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 000113
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 000119
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
                                // 000132
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 000133
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 000139
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
                                // 000182
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
                                // 000183
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
                                // 000189
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
                                // 000212
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 000213
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 000219
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
                                // 000232
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 000233
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 000239
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
                                // 000282
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
                                // 000283
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
                                // 000289
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
                                // 000712
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
                                // 000713
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
                                // 000719
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
                                // 000732
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
                                // 000733
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
                                // 000739
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
                                // 000782
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
                                // 000783
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
                                // 000789
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
                                // 003112
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 003113
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 003119
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
                                // 003132
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 003133
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 003139
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
                                // 003182
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
                                // 003183
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
                                // 003189
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
                                // 003212
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 003213
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 003219
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
                                // 003232
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 003233
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 003239
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
                                // 003282
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
                                // 003283
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
                                // 003289
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
                                // 003712
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
                                // 003713
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
                                // 003719
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
                                // 003732
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
                                // 003733
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
                                // 003739
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
                                // 003782
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
                                // 003783
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
                                // 003789
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
                                // 006112
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
                                // 006113
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
                                // 006119
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
                                // 006132
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
                                // 006133
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
                                // 006139
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
                                // 006182
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
                                // 006183
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
                                // 006189
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
                                // 006212
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
                                // 006213
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
                                // 006219
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
                                // 006232
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
                                // 006233
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
                                // 006239
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
                                // 006282
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
                                // 006283
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
                                // 006289
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
                                // 006712
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
                                // 006713
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
                                // 006719
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
                                // 006732
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
                                // 006733
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
                                // 006739
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
                                // 006782
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
                                // 006783
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
                                // 006789
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
                                // 020112
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 020113
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 020119
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
                                // 020132
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 020133
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 020139
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
                                // 020182
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
                                // 020183
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
                                // 020189
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
                                // 020212
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 020213
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 020219
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
                                // 020232
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 020233
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 020239
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
                                // 020282
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
                                // 020283
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
                                // 020289
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
                                // 020712
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
                                // 020713
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
                                // 020719
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
                                // 020732
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
                                // 020733
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
                                // 020739
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
                                // 020782
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
                                // 020783
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
                                // 020789
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
                                // 023112
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 023113
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 023119
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
                                // 023132
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 023133
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 023139
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
                                // 023182
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
                                // 023183
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
                                // 023189
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
                                // 023212
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 023213
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 023219
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
                                // 023232
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 023233
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 023239
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
                                // 023282
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
                                // 023283
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
                                // 023289
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
                                // 023712
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
                                // 023713
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
                                // 023719
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
                                // 023732
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
                                // 023733
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
                                // 023739
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
                                // 023782
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
                                // 023783
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
                                // 023789
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
                                // 026112
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
                                // 026113
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
                                // 026119
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
                                // 026132
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
                                // 026133
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
                                // 026139
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
                                // 026182
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
                                // 026183
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
                                // 026189
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
                                // 026212
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
                                // 026213
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
                                // 026219
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
                                // 026232
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
                                // 026233
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
                                // 026239
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
                                // 026282
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
                                // 026283
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
                                // 026289
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
                                // 026712
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
                                // 026713
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
                                // 026719
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
                                // 026732
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
                                // 026733
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
                                // 026739
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
                                // 026782
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
                                // 026783
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
                                // 026789
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
                                // 050112
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
                                // 050113
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
                                // 050119
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
                                // 050132
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
                                // 050133
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
                                // 050139
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
                                // 050182
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
                                // 050183
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
                                // 050189
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
                                // 050212
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
                                // 050213
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
                                // 050219
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
                                // 050232
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
                                // 050233
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
                                // 050239
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
                                // 050282
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
                                // 050283
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
                                // 050289
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
                                // 050712
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
                                // 050713
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
                                // 050719
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
                                // 050732
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
                                // 050733
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
                                // 050739
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
                                // 050782
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
                                // 050783
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
                                // 050789
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
                                // 053112
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
                                // 053113
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
                                // 053119
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
                                // 053132
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
                                // 053133
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
                                // 053139
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
                                // 053182
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
                                // 053183
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
                                // 053189
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
                                // 053212
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
                                // 053213
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
                                // 053219
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
                                // 053232
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
                                // 053233
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
                                // 053239
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
                                // 053282
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
                                // 053283
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
                                // 053289
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
                                // 053712
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
                                // 053713
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
                                // 053719
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
                                // 053732
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
                                // 053733
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
                                // 053739
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
                                // 053782
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
                                // 053783
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
                                // 053789
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
                                // 056112
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
                                // 056113
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
                                // 056119
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
                                // 056132
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
                                // 056133
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
                                // 056139
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
                                // 056182
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
                                // 056183
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
                                // 056189
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
                                // 056212
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
                                // 056213
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
                                // 056219
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
                                // 056232
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
                                // 056233
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
                                // 056239
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
                                // 056282
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
                                // 056283
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
                                // 056289
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
                                // 056712
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
                                // 056713
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
                                // 056719
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
                                // 056732
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
                                // 056733
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
                                // 056739
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
                                // 056782
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
                                // 056783
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
                                // 056789
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
                                // 100112
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 100113
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 100119
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
                                // 100132
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 100133
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 100139
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
                                // 100182
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
                                // 100183
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
                                // 100189
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
                                // 100212
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 100213
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 100219
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
                                // 100232
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 100233
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 100239
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
                                // 100282
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
                                // 100283
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
                                // 100289
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
                                // 100712
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
                                // 100713
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
                                // 100719
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
                                // 100732
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
                                // 100733
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
                                // 100739
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
                                // 100782
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
                                // 100783
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
                                // 100789
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
                                // 103112
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 103113
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 103119
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
                                // 103132
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 103133
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 103139
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
                                // 103182
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
                                // 103183
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
                                // 103189
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
                                // 103212
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 103213
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 103219
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
                                // 103232
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 103233
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 103239
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
                                // 103282
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
                                // 103283
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
                                // 103289
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
                                // 103712
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
                                // 103713
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
                                // 103719
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
                                // 103732
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
                                // 103733
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
                                // 103739
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
                                // 103782
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
                                // 103783
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
                                // 103789
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
                                // 106112
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
                                // 106113
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
                                // 106119
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
                                // 106132
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
                                // 106133
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
                                // 106139
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
                                // 106182
								(*nel) = 3;
								t[0] = 2;
								t[1] = 3;
								t[2] = 6;
								t[3] = 8;
								t[4] = 6;
								t[5] = 1;
								t[6] = 2;
								t[7] = 8;
								t[8] = 6;
								t[9] = 1;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 3) {
                                // 106183
								(*nel) = 3;
								t[0] = 2;
								t[1] = 3;
								t[2] = 6;
								t[3] = 8;
								t[4] = 6;
								t[5] = 1;
								t[6] = 2;
								t[7] = 8;
								t[8] = 6;
								t[9] = 1;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 9) {
                                // 106189
								(*nel) = 4;
								t[0] = 6;
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
								t[13] = 1;
								t[14] = 0;
								t[15] = 9;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                // 106212
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
                                // 106213
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
                                // 106219
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
                                // 106232
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
                                // 106233
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
                                // 106239
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
                                // 106282
								(*nel) = 3;
								t[0] = 2;
								t[1] = 3;
								t[2] = 6;
								t[3] = 8;
								t[4] = 6;
								t[5] = 1;
								t[6] = 2;
								t[7] = 8;
								t[8] = 6;
								t[9] = 1;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 3) {
                                // 106283
								(*nel) = 3;
								t[0] = 2;
								t[1] = 3;
								t[2] = 6;
								t[3] = 8;
								t[4] = 6;
								t[5] = 1;
								t[6] = 2;
								t[7] = 8;
								t[8] = 6;
								t[9] = 1;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 9) {
                                // 106289
								(*nel) = 9;
								t[0] = 6;
								t[1] = 1;
								t[2] = 0;
								t[3] = 10;
								t[4] = 8;
								t[5] = 2;
								t[6] = 10;
								t[7] = 9;
								t[8] = 2;
								t[9] = 0;
								t[10] = 1;
								t[11] = 10;
								t[12] = 3;
								t[13] = 6;
								t[14] = 9;
								t[15] = 8;
								t[16] = 6;
								t[17] = 0;
								t[18] = 9;
								t[19] = 10;
								t[20] = 0;
								t[21] = 2;
								t[22] = 9;
								t[23] = 10;
								t[24] = 8;
								t[25] = 2;
								t[26] = 1;
								t[27] = 10;
								t[28] = 6;
								t[29] = 9;
								t[30] = 8;
								t[31] = 10;
								t[32] = 1;
								t[33] = 6;
								t[34] = 8;
								t[35] = 10;
								(*nint) = 1;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                // 106712
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
                                // 106713
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
                                // 106719
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
                                // 106732
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
                                // 106733
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
                                // 106739
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
                                // 106782
								(*nel) = 5;
								t[0] = 1;
								t[1] = 6;
								t[2] = 8;
								t[3] = 7;
								t[4] = 1;
								t[5] = 6;
								t[6] = 7;
								t[7] = 0;
								t[8] = 3;
								t[9] = 2;
								t[10] = 8;
								t[11] = 6;
								t[12] = 2;
								t[13] = 8;
								t[14] = 6;
								t[15] = 7;
								t[16] = 0;
								t[17] = 2;
								t[18] = 6;
								t[19] = 7;
                            } else if (edges[5] == 3) {
                                // 106783
								(*nel) = 5;
								t[0] = 3;
								t[1] = 7;
								t[2] = 6;
								t[3] = 2;
								t[4] = 1;
								t[5] = 6;
								t[6] = 8;
								t[7] = 7;
								t[8] = 1;
								t[9] = 6;
								t[10] = 7;
								t[11] = 0;
								t[12] = 3;
								t[13] = 7;
								t[14] = 8;
								t[15] = 6;
								t[16] = 0;
								t[17] = 2;
								t[18] = 6;
								t[19] = 7;
                            } else if (edges[5] == 9) {
                                // 106789
								(*nel) = 6;
								t[0] = 8;
								t[1] = 3;
								t[2] = 9;
								t[3] = 6;
								t[4] = 6;
								t[5] = 7;
								t[6] = 0;
								t[7] = 9;
								t[8] = 0;
								t[9] = 2;
								t[10] = 9;
								t[11] = 7;
								t[12] = 6;
								t[13] = 7;
								t[14] = 9;
								t[15] = 8;
								t[16] = 6;
								t[17] = 1;
								t[18] = 7;
								t[19] = 8;
								t[20] = 6;
								t[21] = 1;
								t[22] = 0;
								t[23] = 7;
                            }
                        }
                    }
                }
            } else if (edges[1] == 2) {
                if (edges[2] == 0) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                // 120112
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 120113
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 120119
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
                                // 120132
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 120133
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 120139
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
                                // 120182
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
                                // 120183
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
                                // 120189
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
                                // 120212
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 120213
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 120219
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
                                // 120232
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 120233
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 120239
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
                                // 120282
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
                                // 120283
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
                                // 120289
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
                                // 120712
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
                                // 120713
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
                                // 120719
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
                                // 120732
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
                                // 120733
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
                                // 120739
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
                                // 120782
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
                                // 120783
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
                                // 120789
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
                                // 123112
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 123113
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 123119
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
                                // 123132
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 123133
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 123139
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
                                // 123182
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
                                // 123183
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
                                // 123189
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
                                // 123212
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 123213
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 123219
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
                                // 123232
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 3) {
                                // 123233
								(*nel) = 1;
								t[0] = 0;
								t[1] = 1;
								t[2] = 2;
								t[3] = 3;
                            } else if (edges[5] == 9) {
                                // 123239
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
                                // 123282
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
                                // 123283
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
                                // 123289
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
                                // 123712
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
                                // 123713
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
                                // 123719
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
                                // 123732
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
                                // 123733
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
                                // 123739
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
                                // 123782
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
                                // 123783
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
                                // 123789
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
                                // 126112
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
                                // 126113
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
                                // 126119
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
                                // 126132
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
                                // 126133
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
                                // 126139
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
                                // 126182
								(*nel) = 3;
								t[0] = 2;
								t[1] = 3;
								t[2] = 6;
								t[3] = 8;
								t[4] = 6;
								t[5] = 1;
								t[6] = 2;
								t[7] = 8;
								t[8] = 6;
								t[9] = 1;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 3) {
                                // 126183
								(*nel) = 3;
								t[0] = 2;
								t[1] = 3;
								t[2] = 6;
								t[3] = 8;
								t[4] = 6;
								t[5] = 1;
								t[6] = 2;
								t[7] = 8;
								t[8] = 6;
								t[9] = 1;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 9) {
                                // 126189
								(*nel) = 4;
								t[0] = 1;
								t[1] = 6;
								t[2] = 8;
								t[3] = 9;
								t[4] = 6;
								t[5] = 2;
								t[6] = 9;
								t[7] = 1;
								t[8] = 3;
								t[9] = 6;
								t[10] = 9;
								t[11] = 8;
								t[12] = 6;
								t[13] = 2;
								t[14] = 1;
								t[15] = 0;
                            }
                        }
                    } else if (edges[3] == 2) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                // 126212
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
                                // 126213
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
                                // 126219
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
                                // 126232
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
                                // 126233
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
                                // 126239
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
                                // 126282
								(*nel) = 3;
								t[0] = 2;
								t[1] = 3;
								t[2] = 6;
								t[3] = 8;
								t[4] = 6;
								t[5] = 1;
								t[6] = 2;
								t[7] = 8;
								t[8] = 6;
								t[9] = 1;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 3) {
                                // 126283
								(*nel) = 3;
								t[0] = 2;
								t[1] = 3;
								t[2] = 6;
								t[3] = 8;
								t[4] = 6;
								t[5] = 1;
								t[6] = 2;
								t[7] = 8;
								t[8] = 6;
								t[9] = 1;
								t[10] = 0;
								t[11] = 2;
                            } else if (edges[5] == 9) {
                                // 126289
								(*nel) = 4;
								t[0] = 8;
								t[1] = 2;
								t[2] = 1;
								t[3] = 6;
								t[4] = 3;
								t[5] = 6;
								t[6] = 9;
								t[7] = 8;
								t[8] = 6;
								t[9] = 2;
								t[10] = 1;
								t[11] = 0;
								t[12] = 8;
								t[13] = 2;
								t[14] = 6;
								t[15] = 9;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                // 126712
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
                                // 126713
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
                                // 126719
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
                                // 126732
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
                                // 126733
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
                                // 126739
								(*nel) = 5;
								t[0] = 3;
								t[1] = 7;
								t[2] = 1;
								t[3] = 6;
								t[4] = 2;
								t[5] = 6;
								t[6] = 0;
								t[7] = 7;
								t[8] = 2;
								t[9] = 6;
								t[10] = 7;
								t[11] = 9;
								t[12] = 3;
								t[13] = 7;
								t[14] = 6;
								t[15] = 9;
								t[16] = 1;
								t[17] = 0;
								t[18] = 6;
								t[19] = 7;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                // 126782
								(*nel) = 5;
								t[0] = 6;
								t[1] = 1;
								t[2] = 7;
								t[3] = 8;
								t[4] = 3;
								t[5] = 2;
								t[6] = 8;
								t[7] = 6;
								t[8] = 6;
								t[9] = 1;
								t[10] = 0;
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
                                // 126783
								(*nel) = 5;
								t[0] = 3;
								t[1] = 7;
								t[2] = 6;
								t[3] = 2;
								t[4] = 1;
								t[5] = 6;
								t[6] = 8;
								t[7] = 7;
								t[8] = 1;
								t[9] = 6;
								t[10] = 7;
								t[11] = 0;
								t[12] = 3;
								t[13] = 7;
								t[14] = 8;
								t[15] = 6;
								t[16] = 0;
								t[17] = 2;
								t[18] = 6;
								t[19] = 7;
                            } else if (edges[5] == 9) {
                                // 126789
								(*nel) = 6;
								t[0] = 8;
								t[1] = 3;
								t[2] = 9;
								t[3] = 6;
								t[4] = 6;
								t[5] = 2;
								t[6] = 7;
								t[7] = 0;
								t[8] = 6;
								t[9] = 1;
								t[10] = 7;
								t[11] = 8;
								t[12] = 6;
								t[13] = 2;
								t[14] = 9;
								t[15] = 7;
								t[16] = 1;
								t[17] = 6;
								t[18] = 7;
								t[19] = 0;
								t[20] = 7;
								t[21] = 6;
								t[22] = 8;
								t[23] = 9;
                            }
                        }
                    }
                }
            } else if (edges[1] == 5) {
                if (edges[2] == 0) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                // 150112
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
                                // 150113
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
                                // 150119
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
                                // 150132
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
                                // 150133
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
                                // 150139
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
                                // 150182
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
                                // 150183
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
                                // 150189
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
                                // 150212
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
                                // 150213
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
                                // 150219
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
                                // 150232
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
                                // 150233
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
                                // 150239
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
                                // 150282
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
                                // 150283
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
                                // 150289
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
                                // 150712
								(*nel) = 3;
								t[0] = 3;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 1;
								t[5] = 5;
								t[6] = 3;
								t[7] = 7;
								t[8] = 1;
								t[9] = 5;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 3) {
                                // 150713
								(*nel) = 3;
								t[0] = 3;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 5;
								t[5] = 1;
								t[6] = 3;
								t[7] = 0;
								t[8] = 5;
								t[9] = 1;
								t[10] = 7;
								t[11] = 3;
                            } else if (edges[5] == 9) {
                                // 150719
								(*nel) = 4;
								t[0] = 9;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 5;
								t[5] = 1;
								t[6] = 9;
								t[7] = 0;
								t[8] = 9;
								t[9] = 0;
								t[10] = 1;
								t[11] = 3;
								t[12] = 5;
								t[13] = 1;
								t[14] = 7;
								t[15] = 9;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                // 150732
								(*nel) = 3;
								t[0] = 3;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 1;
								t[5] = 5;
								t[6] = 3;
								t[7] = 7;
								t[8] = 1;
								t[9] = 5;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 3) {
                                // 150733
								(*nel) = 3;
								t[0] = 3;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 1;
								t[5] = 5;
								t[6] = 3;
								t[7] = 7;
								t[8] = 1;
								t[9] = 5;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 9) {
                                // 150739
								(*nel) = 9;
								t[0] = 9;
								t[1] = 5;
								t[2] = 7;
								t[3] = 10;
								t[4] = 9;
								t[5] = 5;
								t[6] = 2;
								t[7] = 7;
								t[8] = 0;
								t[9] = 5;
								t[10] = 9;
								t[11] = 10;
								t[12] = 3;
								t[13] = 7;
								t[14] = 10;
								t[15] = 9;
								t[16] = 5;
								t[17] = 1;
								t[18] = 7;
								t[19] = 10;
								t[20] = 1;
								t[21] = 5;
								t[22] = 0;
								t[23] = 10;
								t[24] = 1;
								t[25] = 0;
								t[26] = 3;
								t[27] = 10;
								t[28] = 3;
								t[29] = 7;
								t[30] = 1;
								t[31] = 10;
								t[32] = 0;
								t[33] = 9;
								t[34] = 3;
								t[35] = 10;
								(*nint) = 1;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                // 150782
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
								t[12] = 1;
								t[13] = 5;
								t[14] = 8;
								t[15] = 7;
								t[16] = 1;
								t[17] = 5;
								t[18] = 0;
								t[19] = 8;
                            } else if (edges[5] == 3) {
                                // 150783
								(*nel) = 5;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 8;
								t[4] = 3;
								t[5] = 7;
								t[6] = 8;
								t[7] = 5;
								t[8] = 1;
								t[9] = 5;
								t[10] = 8;
								t[11] = 7;
								t[12] = 1;
								t[13] = 5;
								t[14] = 0;
								t[15] = 8;
								t[16] = 3;
								t[17] = 7;
								t[18] = 5;
								t[19] = 2;
                            } else if (edges[5] == 9) {
                                // 150789
								(*nel) = 6;
								t[0] = 5;
								t[1] = 1;
								t[2] = 8;
								t[3] = 0;
								t[4] = 0;
								t[5] = 3;
								t[6] = 8;
								t[7] = 9;
								t[8] = 8;
								t[9] = 5;
								t[10] = 9;
								t[11] = 7;
								t[12] = 9;
								t[13] = 5;
								t[14] = 2;
								t[15] = 7;
								t[16] = 1;
								t[17] = 5;
								t[18] = 8;
								t[19] = 7;
								t[20] = 8;
								t[21] = 5;
								t[22] = 0;
								t[23] = 9;
                            }
                        }
                    }
                } else if (edges[2] == 3) {
                    if (edges[3] == 1) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                // 153112
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
                                // 153113
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
                                // 153119
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
                                // 153132
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
                                // 153133
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
                                // 153139
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
                                // 153182
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
                                // 153183
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
                                // 153189
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
                                // 153212
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
                                // 153213
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
                                // 153219
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
                                // 153232
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
                                // 153233
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
                                // 153239
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
                                // 153282
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
                                // 153283
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
                                // 153289
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
                                // 153712
								(*nel) = 3;
								t[0] = 3;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 1;
								t[5] = 5;
								t[6] = 3;
								t[7] = 7;
								t[8] = 1;
								t[9] = 5;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 3) {
                                // 153713
								(*nel) = 3;
								t[0] = 3;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 1;
								t[5] = 5;
								t[6] = 3;
								t[7] = 7;
								t[8] = 1;
								t[9] = 5;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 9) {
                                // 153719
								(*nel) = 4;
								t[0] = 9;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 1;
								t[5] = 9;
								t[6] = 7;
								t[7] = 5;
								t[8] = 3;
								t[9] = 5;
								t[10] = 1;
								t[11] = 0;
								t[12] = 1;
								t[13] = 9;
								t[14] = 5;
								t[15] = 3;
                            }
                        } else if (edges[4] == 3) {
                            if (edges[5] == 2) {
                                // 153732
								(*nel) = 3;
								t[0] = 3;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 1;
								t[5] = 5;
								t[6] = 3;
								t[7] = 7;
								t[8] = 1;
								t[9] = 5;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 3) {
                                // 153733
								(*nel) = 3;
								t[0] = 3;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 1;
								t[5] = 5;
								t[6] = 3;
								t[7] = 7;
								t[8] = 1;
								t[9] = 5;
								t[10] = 0;
								t[11] = 3;
                            } else if (edges[5] == 9) {
                                // 153739
								(*nel) = 4;
								t[0] = 9;
								t[1] = 5;
								t[2] = 2;
								t[3] = 7;
								t[4] = 1;
								t[5] = 5;
								t[6] = 0;
								t[7] = 3;
								t[8] = 1;
								t[9] = 5;
								t[10] = 3;
								t[11] = 7;
								t[12] = 5;
								t[13] = 3;
								t[14] = 7;
								t[15] = 9;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                // 153782
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
								t[12] = 1;
								t[13] = 5;
								t[14] = 8;
								t[15] = 7;
								t[16] = 1;
								t[17] = 5;
								t[18] = 0;
								t[19] = 8;
                            } else if (edges[5] == 3) {
                                // 153783
								(*nel) = 5;
								t[0] = 0;
								t[1] = 5;
								t[2] = 3;
								t[3] = 8;
								t[4] = 1;
								t[5] = 5;
								t[6] = 8;
								t[7] = 7;
								t[8] = 3;
								t[9] = 7;
								t[10] = 8;
								t[11] = 5;
								t[12] = 3;
								t[13] = 7;
								t[14] = 5;
								t[15] = 2;
								t[16] = 1;
								t[17] = 5;
								t[18] = 0;
								t[19] = 8;
                            } else if (edges[5] == 9) {
                                // 153789
								(*nel) = 6;
								t[0] = 1;
								t[1] = 5;
								t[2] = 8;
								t[3] = 7;
								t[4] = 5;
								t[5] = 8;
								t[6] = 7;
								t[7] = 9;
								t[8] = 9;
								t[9] = 5;
								t[10] = 2;
								t[11] = 7;
								t[12] = 5;
								t[13] = 3;
								t[14] = 8;
								t[15] = 9;
								t[16] = 1;
								t[17] = 5;
								t[18] = 0;
								t[19] = 8;
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
                                // 156112
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
                                // 156113
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
                                // 156119
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
                                // 156132
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
                                // 156133
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
                                // 156139
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
                                // 156182
								(*nel) = 5;
								t[0] = 8;
								t[1] = 5;
								t[2] = 2;
								t[3] = 1;
								t[4] = 6;
								t[5] = 1;
								t[6] = 5;
								t[7] = 8;
								t[8] = 6;
								t[9] = 1;
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
                            } else if (edges[5] == 3) {
                                // 156183
								(*nel) = 5;
								t[0] = 8;
								t[1] = 5;
								t[2] = 2;
								t[3] = 1;
								t[4] = 6;
								t[5] = 1;
								t[6] = 5;
								t[7] = 8;
								t[8] = 6;
								t[9] = 1;
								t[10] = 0;
								t[11] = 5;
								t[12] = 3;
								t[13] = 5;
								t[14] = 8;
								t[15] = 6;
								t[16] = 3;
								t[17] = 5;
								t[18] = 2;
								t[19] = 8;
                            } else if (edges[5] == 9) {
                                // 156189
								(*nel) = 6;
								t[0] = 5;
								t[1] = 6;
								t[2] = 8;
								t[3] = 9;
								t[4] = 8;
								t[5] = 6;
								t[6] = 3;
								t[7] = 9;
								t[8] = 6;
								t[9] = 1;
								t[10] = 0;
								t[11] = 5;
								t[12] = 6;
								t[13] = 1;
								t[14] = 5;
								t[15] = 8;
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
                                // 156212
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
                                // 156213
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
                                // 156219
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
                                // 156232
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
                                // 156233
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
                                // 156239
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
                                // 156282
								(*nel) = 5;
								t[0] = 8;
								t[1] = 5;
								t[2] = 2;
								t[3] = 1;
								t[4] = 6;
								t[5] = 1;
								t[6] = 5;
								t[7] = 8;
								t[8] = 6;
								t[9] = 1;
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
                            } else if (edges[5] == 3) {
                                // 156283
								(*nel) = 5;
								t[0] = 8;
								t[1] = 5;
								t[2] = 2;
								t[3] = 1;
								t[4] = 6;
								t[5] = 1;
								t[6] = 5;
								t[7] = 8;
								t[8] = 6;
								t[9] = 1;
								t[10] = 0;
								t[11] = 5;
								t[12] = 3;
								t[13] = 5;
								t[14] = 8;
								t[15] = 6;
								t[16] = 3;
								t[17] = 5;
								t[18] = 2;
								t[19] = 8;
                            } else if (edges[5] == 9) {
                                // 156289
								(*nel) = 6;
								t[0] = 5;
								t[1] = 6;
								t[2] = 8;
								t[3] = 9;
								t[4] = 8;
								t[5] = 6;
								t[6] = 3;
								t[7] = 9;
								t[8] = 6;
								t[9] = 1;
								t[10] = 0;
								t[11] = 5;
								t[12] = 8;
								t[13] = 2;
								t[14] = 5;
								t[15] = 9;
								t[16] = 6;
								t[17] = 1;
								t[18] = 5;
								t[19] = 8;
								t[20] = 8;
								t[21] = 2;
								t[22] = 1;
								t[23] = 5;
                            }
                        }
                    } else if (edges[3] == 7) {
                        if (edges[4] == 1) {
                            if (edges[5] == 2) {
                                // 156712
								(*nel) = 5;
								t[0] = 7;
								t[1] = 3;
								t[2] = 6;
								t[3] = 1;
								t[4] = 1;
								t[5] = 5;
								t[6] = 6;
								t[7] = 7;
								t[8] = 7;
								t[9] = 5;
								t[10] = 6;
								t[11] = 2;
								t[12] = 3;
								t[13] = 7;
								t[14] = 6;
								t[15] = 2;
								t[16] = 1;
								t[17] = 5;
								t[18] = 0;
								t[19] = 6;
                            } else if (edges[5] == 3) {
                                // 156713
								(*nel) = 5;
								t[0] = 7;
								t[1] = 3;
								t[2] = 6;
								t[3] = 1;
								t[4] = 3;
								t[5] = 5;
								t[6] = 7;
								t[7] = 6;
								t[8] = 1;
								t[9] = 5;
								t[10] = 6;
								t[11] = 7;
								t[12] = 3;
								t[13] = 5;
								t[14] = 2;
								t[15] = 7;
								t[16] = 1;
								t[17] = 5;
								t[18] = 0;
								t[19] = 6;
                            } else if (edges[5] == 9) {
                                // 156719
								(*nel) = 6;
								t[0] = 2;
								t[1] = 7;
								t[2] = 9;
								t[3] = 5;
								t[4] = 1;
								t[5] = 5;
								t[6] = 6;
								t[7] = 7;
								t[8] = 5;
								t[9] = 6;
								t[10] = 7;
								t[11] = 9;
								t[12] = 1;
								t[13] = 5;
								t[14] = 0;
								t[15] = 6;
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
                                // 156732
								(*nel) = 5;
								t[0] = 7;
								t[1] = 3;
								t[2] = 6;
								t[3] = 1;
								t[4] = 1;
								t[5] = 5;
								t[6] = 6;
								t[7] = 7;
								t[8] = 7;
								t[9] = 5;
								t[10] = 6;
								t[11] = 2;
								t[12] = 3;
								t[13] = 7;
								t[14] = 6;
								t[15] = 2;
								t[16] = 1;
								t[17] = 5;
								t[18] = 0;
								t[19] = 6;
                            } else if (edges[5] == 3) {
                                // 156733
								(*nel) = 5;
								t[0] = 1;
								t[1] = 5;
								t[2] = 6;
								t[3] = 7;
								t[4] = 7;
								t[5] = 3;
								t[6] = 6;
								t[7] = 1;
								t[8] = 5;
								t[9] = 3;
								t[10] = 7;
								t[11] = 2;
								t[12] = 5;
								t[13] = 3;
								t[14] = 6;
								t[15] = 7;
								t[16] = 1;
								t[17] = 5;
								t[18] = 0;
								t[19] = 6;
                            } else if (edges[5] == 9) {
                                // 156739
								(*nel) = 6;
								t[0] = 2;
								t[1] = 7;
								t[2] = 9;
								t[3] = 5;
								t[4] = 6;
								t[5] = 5;
								t[6] = 9;
								t[7] = 7;
								t[8] = 1;
								t[9] = 5;
								t[10] = 0;
								t[11] = 6;
								t[12] = 3;
								t[13] = 7;
								t[14] = 6;
								t[15] = 9;
								t[16] = 3;
								t[17] = 7;
								t[18] = 1;
								t[19] = 6;
								t[20] = 1;
								t[21] = 5;
								t[22] = 6;
								t[23] = 7;
                            }
                        } else if (edges[4] == 8) {
                            if (edges[5] == 2) {
                                // 156782
								(*nel) = 6;
								t[0] = 8;
								t[1] = 1;
								t[2] = 5;
								t[3] = 7;
								t[4] = 6;
								t[5] = 2;
								t[6] = 8;
								t[7] = 5;
								t[8] = 1;
								t[9] = 6;
								t[10] = 5;
								t[11] = 0;
								t[12] = 2;
								t[13] = 6;
								t[14] = 8;
								t[15] = 3;
								t[16] = 6;
								t[17] = 1;
								t[18] = 5;
								t[19] = 8;
								t[20] = 2;
								t[21] = 8;
								t[22] = 5;
								t[23] = 7;
                            } else if (edges[5] == 3) {
                                // 156783
								(*nel) = 6;
								t[0] = 8;
								t[1] = 5;
								t[2] = 3;
								t[3] = 7;
								t[4] = 5;
								t[5] = 8;
								t[6] = 3;
								t[7] = 6;
								t[8] = 3;
								t[9] = 5;
								t[10] = 2;
								t[11] = 7;
								t[12] = 1;
								t[13] = 5;
								t[14] = 0;
								t[15] = 6;
								t[16] = 8;
								t[17] = 5;
								t[18] = 1;
								t[19] = 6;
								t[20] = 8;
								t[21] = 5;
								t[22] = 7;
								t[23] = 1;
                            } else if (edges[5] == 9) {
                                // 156789
								(*nel) = 7;
								t[0] = 5;
								t[1] = 8;
								t[2] = 1;
								t[3] = 7;
								t[4] = 8;
								t[5] = 3;
								t[6] = 9;
								t[7] = 6;
								t[8] = 5;
								t[9] = 8;
								t[10] = 9;
								t[11] = 6;
								t[12] = 8;
								t[13] = 5;
								t[14] = 1;
								t[15] = 6;
								t[16] = 2;
								t[17] = 9;
								t[18] = 5;
								t[19] = 7;
								t[20] = 1;
								t[21] = 5;
								t[22] = 0;
								t[23] = 6;
								t[24] = 5;
								t[25] = 8;
								t[26] = 7;
								t[27] = 9;
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
                                // 400112
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
                                // 400113
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
                                // 400119
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
                                // 400132
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
                                // 400133
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
                                // 400139
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
                                // 400182
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
                                // 400183
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
                                // 400189
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
                                // 400212
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
                                // 400213
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
                                // 400219
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
                                // 400232
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
                                // 400233
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
                                // 400239
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
                                // 400282
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
                                // 400283
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
                                // 400289
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
                                // 400712
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
                                // 400713
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
                                // 400719
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
                                // 400732
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
                                // 400733
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
                                // 400739
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
                                // 400782
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
                                // 400783
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
                                // 400789
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
                                // 403112
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
                                // 403113
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
                                // 403119
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
                                // 403132
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
                                // 403133
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
                                // 403139
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
                                // 403182
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
                                // 403183
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
                                // 403189
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
                                // 403212
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
                                // 403213
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
                                // 403219
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
                                // 403232
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
                                // 403233
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
                                // 403239
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
                                // 403282
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
                                // 403283
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
                                // 403289
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
                                // 403712
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
                                // 403713
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
                                // 403719
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
                                // 403732
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
                                // 403733
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
                                // 403739
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
                                // 403782
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
                                // 403783
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
                                // 403789
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
                                // 406112
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
                                // 406113
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
                                // 406119
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
                                // 406132
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
                                // 406133
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
                                // 406139
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
                                // 406182
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
                                // 406183
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
                                // 406189
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
                                // 406212
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
                                // 406213
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
                                // 406219
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
                                // 406232
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
                                // 406233
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
                                // 406239
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
                                // 406282
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
                                // 406283
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
                                // 406289
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
                                // 406712
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
                                // 406713
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
                                // 406719
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
                                // 406732
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
                                // 406733
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
                                // 406739
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
                                // 406782
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
                                // 406783
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
                                // 406789
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
                                // 420112
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
                                // 420113
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
                                // 420119
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
                                // 420132
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
                                // 420133
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
                                // 420139
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
                                // 420182
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
                                // 420183
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
                                // 420189
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
                                // 420212
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
                                // 420213
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
                                // 420219
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
                                // 420232
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
                                // 420233
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
                                // 420239
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
                                // 420282
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
                                // 420283
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
                                // 420289
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
                                // 420712
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
                                // 420713
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
                                // 420719
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
                                // 420732
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
                                // 420733
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
                                // 420739
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
                                // 420782
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
                                // 420783
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
                                // 420789
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
                                // 423112
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
                                // 423113
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
                                // 423119
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
                                // 423132
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
                                // 423133
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
                                // 423139
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
                                // 423182
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
                                // 423183
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
                                // 423189
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
                                // 423212
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
                                // 423213
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
                                // 423219
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
                                // 423232
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
                                // 423233
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
                                // 423239
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
                                // 423282
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
                                // 423283
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
                                // 423289
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
                                // 423712
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
                                // 423713
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
                                // 423719
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
                                // 423732
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
                                // 423733
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
                                // 423739
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
                                // 423782
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
                                // 423783
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
                                // 423789
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
                                // 426112
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
                                // 426113
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
                                // 426119
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
                                // 426132
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
                                // 426133
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
                                // 426139
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
                                // 426182
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
                                // 426183
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
                                // 426189
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
                                // 426212
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
                                // 426213
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
                                // 426219
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
                                // 426232
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
                                // 426233
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
                                // 426239
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
                                // 426282
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
                                // 426283
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
                                // 426289
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
                                // 426712
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
                                // 426713
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
                                // 426719
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
                                // 426732
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
                                // 426733
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
                                // 426739
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
                                // 426782
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
                                // 426783
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
                                // 426789
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
                                // 450112
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
                                // 450113
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
                                // 450119
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
                                // 450132
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
                                // 450133
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
                                // 450139
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
                                // 450182
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
                                // 450183
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
                                // 450189
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
                                // 450212
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
                                // 450213
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
                                // 450219
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
                                // 450232
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
                                // 450233
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
                                // 450239
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
                                // 450282
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
                                // 450283
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
                                // 450289
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
                                // 450712
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
                                // 450713
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
                                // 450719
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
                                // 450732
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
                                // 450733
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
                                // 450739
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
                                // 450782
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
                                // 450783
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
                                // 450789
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
                                // 453112
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
                                // 453113
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
                                // 453119
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
                                // 453132
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
                                // 453133
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
                                // 453139
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
                                // 453182
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
                                // 453183
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
                                // 453189
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
                                // 453212
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
                                // 453213
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
                                // 453219
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
                                // 453232
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
                                // 453233
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
                                // 453239
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
                                // 453282
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
                                // 453283
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
                                // 453289
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
                                // 453712
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
                                // 453713
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
                                // 453719
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
                                // 453732
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
                                // 453733
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
                                // 453739
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
                                // 453782
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
                                // 453783
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
                                // 453789
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
                                // 456112
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
                                // 456113
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
                                // 456119
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
                                // 456132
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
                                // 456133
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
                                // 456139
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
                                // 456182
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
                                // 456183
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
                                // 456189
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
                                // 456212
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
                                // 456213
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
                                // 456219
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
                                // 456232
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
                                // 456233
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
                                // 456239
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
                                // 456282
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
                                // 456283
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
                                // 456289
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
                                // 456712
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
                                // 456713
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
                                // 456719
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
                                // 456732
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
                                // 456733
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
                                // 456739
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
                                // 456782
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
                                // 456783
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
                                // 456789
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
