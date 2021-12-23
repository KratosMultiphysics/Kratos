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
//
//

#if !defined(KRATOS_SPLIT_PRISM)
#define  KRATOS_SPLIT_PRISM

/** @file split_prism.hpp
 * @brief The class contains three helper functions to ease the splitting:
 * PrismSplitMode, Split_Prism, and PrismGetNewConnectivityGID
 * EXAMPLE: imagine that an user would like to split a prism wit a face formed
 * by the ids 3 9 7 by introducing a new node 15 on the edge between 3 and 9
 * he should define
 * int aux[12]
 * int edge_ids[6]
 * int t[24]
 * then initialize
 * aux[0] = 3; aux[1] = 9; aux[2] = 7;
 * aux[3] = 15; //node on edge 01: edge to be refined
 * aux[4] = -1; //node on edge 12:edge not to be refined
 * aux[5] = -1; //node on edge 20:edge not to be refined
 * then call
 *
 * PrismSplitMode(edge_ids,aux)
 * int number_elem; //number of nodes generated
 * int number_splitted_edges; //number of splitted edges
 * int nint; //number of internal nodes
 * int split_needed = Split_Prism(edge_ids,t, &number_elem, &number_splitted_edges, &nint)
 *
 * The split of the other face is straightforward, in fact this header is an adaptation of the split_triangle.h
 *
 * the new prisms ids can be then inspected by
 * for(int i=0; i<number_elem; i++)
 * {
 *     int i0,i1,i2, i3,i4,i5;
 *     PrismGetNewConnectivityGID(i, t, aux, &i0,&i1,&i2,&i3,&i4,&i5);
 * }
 */

/***********************************************************************************/
/***********************************************************************************/

/** This function computes the splitting mode for the prism
 * @param aux_ids Contains a vector with the input Ids, organized as follows
 *   aux_ids[0] = id of FIRST node of the original prism
 *   aux_ids[1] = id of SECOND node of the original prism
 *   aux_ids[2] = id of THIRD node of the original prism
 *   aux_ids[3] = id of FOURTH node of the original prism
 *   aux_ids[4] = id of FIFTH node of the original prism
 *   aux_ids[5] = id of SIXTH node of the original prism
 *   aux_ids[6] = id of new node to be used for the lower edge 01 (-1 if edge not to be splitted)
 *   aux_ids[7] = id of new node to be used for the lower edge 12 (-1 if edge not to be splitted)
 *   aux_ids[8] = id of new node to be used for the lower edge 20 (-1 if edge not to be splitted)
 *   aux_ids[9] = id of new node to be used for the upper edge 01 (-1 if edge not to be splitted)
 *   aux_ids[10] = id of new node to be used for the upper edge 12 (-1 if edge not to be splitted)
 *   aux_ids[11] = id of new node to be used for the upper edge 20 (-1 if edge not to be splitted)
 *   Given this data it fills an auxiliary vector of size 6 that will be used in the splitting
 * @param edge_ids This is an auxiliary array with the local numbering. It is necessary for the split_prism function
 */

void PrismSplitMode(
  const int aux_ids[12],
  int edge_ids[6]
)
{
    /* Edge 01*/
    if (aux_ids[6] < 0)
    {
        if (aux_ids[0] > aux_ids[1])
	{
	  edge_ids[0] = 0;
	  edge_ids[3] = 3;
	}
        else
	{
	  edge_ids[0] = 1;
	  edge_ids[3] = 4;
	}
    }
    else
    {
        edge_ids[0] = 6;
        edge_ids[3] = 9;
    }

    /* Edge 12*/
    if (aux_ids[7] < 0)
        if (aux_ids[1] > aux_ids[2])
	{
	  edge_ids[1] = 1;
	  edge_ids[4] = 4;
	}
        else
	{
	  edge_ids[1] = 2;
	  edge_ids[4] = 5;
	}
    else
    {
        edge_ids[1] = 7;
        edge_ids[4] = 10;
    }

    /* Edge 20*/
    if (aux_ids[8] < 0)
    {
        if (aux_ids[2] > aux_ids[0])
	{
	  edge_ids[2] = 2;
	  edge_ids[5] = 5;
	}
        else
	{
	  edge_ids[2] = 0;
	  edge_ids[5] = 3;
	}
    }
    else
    {
        edge_ids[2] = 8;
        edge_ids[5] = 11;
    }
}

/***********************************************************************************/
/***********************************************************************************/

/**
* Utility function to get the global ids for the new prisms to be generated
* @param prism_index the index of the new prism to be generated (Should be less than the number number_elem provided by Split_Prism)
 *@param t integer array provided by Split_Prism
* @param aux_ids array used in constructing the edge_ids (contains the Global Ids of the new nodes)
* @param id0 Global ID of node0 of the new prism
* @param id1 Global ID of node1 of the new prism
* @param id2 Global ID of node2 of the new prism
* @param id3 Global ID of node3 of the new prism
* @param id4 Global ID of node4 of the new prism
* @param id5 Global ID of node5 of the new prism
*/

inline void PrismGetNewConnectivityGID(
  const int prism_index,
  const int t[24],
  const int aux_ids[12],
  int* id0,
  int* id1,
  int* id2,
  int* id3,
  int* id4,
  int* id5
)
{
    unsigned int base = prism_index * 6;
    *id0 = aux_ids[t[base]];
    *id1 = aux_ids[t[base + 1]];
    *id2 = aux_ids[t[base + 2]];
    *id3 = aux_ids[t[base + 3]];
    *id4 = aux_ids[t[base + 4]];
    *id5 = aux_ids[t[base + 5]];
}

/**
* Utility to split prisms
* @param edges (input) int c array of size 6
* @param t (output) int c array of size 24 (6*4)
* @param number_elem (output) number of elements in the subdivision
* @param splitted_edges (output) provides the number of splitted edges
* @param nint (output)  internal node (not needed for prisms)
* @return 1-->splitting needed    0-->no splitting needed
*/

int Split_Prism(
  const int  edges[6],
  int  t[24],
  int* number_elem,
  int* splitted_edges,
  int* nint
)
{
    *splitted_edges = 0;
    int topology[6];
    topology[0] = 0;
    topology[1] = 0;
    topology[2] = 0;
    topology[3] = 0;
    topology[4] = 0;
    topology[5] = 0;

    for (unsigned int i = 0; i < 6; i++)
    {
        if (edges[i] > 2)
        {
            topology[i] = 1;
            *splitted_edges = *splitted_edges + 1;
        }
    }

    if (*splitted_edges == 0 && *nint == 0)
    {
        /* No splitting needed */
        *number_elem = 1;
        // Lower face
        t[0] = 0;
        t[1] = 1;
        t[2] = 2;

        // Upper face
        t[3] = 3;
        t[4] = 4;
        t[5] = 5;

        return 0;
    }
    /*WARNING = case new central node needed*/
    else if (*splitted_edges == 0 && *nint == 1)
    {
        *number_elem = 3;
        // Lower face
        t[0] = 6;
        t[1] = 0;
        t[2] = 1;

        t[6] = 6;
        t[7] = 1;
        t[8] = 2;

        t[8] = 6;
        t[8] = 2;
        t[8] = 0;

        // Upper face
        t[3] = 9;
        t[4] = 3;
        t[5] = 4;

        t[9]  = 9;
        t[10] = 4;
        t[11] = 5;

        t[11] = 9;
        t[11] = 5;
        t[11] = 3;

        return 1;
    }
    else if (*splitted_edges == 2)
    {
        *number_elem = 2;
        /* Case 1*/
        if (topology[0] == 1)
        {
            // Lower face
            t[0] = 6;
            t[1] = 2;
            t[2] = 0;

            t[6] = 6;
            t[7] = 1;
            t[8] = 2;

            // Upper face
            t[3] = 9;
            t[4] = 5;
            t[5] = 3;

            t[9]  = 9;
            t[10] = 4;
            t[11] = 5;
        }
        /* Case 2*/
        else if (topology[1] == 1)
        {
            // Lower face
            t[0] = 7;
            t[1] = 0;
            t[2] = 1;

            t[6] = 7;
            t[7] = 2;
            t[8] = 0;

            // Upper face
            t[3] = 10;
            t[4] = 3;
            t[5] = 4;

            t[9]  = 10;
            t[10] = 5;
            t[11] = 3;
        }
        /* Case 3*/
        else if (topology[2] == 1)
        {
            // Lower face
            t[0] = 8;
            t[1] = 1;
            t[2] = 2;

            t[6] = 8;
            t[7] = 0;
            t[8] = 1;

            // Upper face
            t[3] = 11;
            t[4] = 4;
            t[5] = 5;

            t[9]  = 11;
            t[10] = 3;
            t[11] = 4;
        }

        return 1;
    }
    else if (*splitted_edges == 4)
    {
        *number_elem = 3;
        /* Case 4*/
        if (topology[0] == 1 && topology[1] == 1)
        {
            if (edges[2] == 0) // If I colapse to the node 0 local
            {
                // Lower face
                t[0] = 7;
                t[1] = 6;
                t[2] = 1;

                t[6] = 7;
                t[7] = 0;
                t[8] = 6;

                t[12] = 7;
                t[13] = 2;
                t[14] = 0;

                // Upper face
                t[3] = 10;
                t[4] = 9;
                t[5] = 4;

                t[9]  = 10;
                t[10] = 3;
                t[11] = 9;

                t[15] = 10;
                t[16] = 5;
                t[17] = 3;

            }
            else if (edges[2] == 2) // If I colapse to the node 2 local
            {
                // Lower face
                t[0] = 7;
                t[1] = 6;
                t[2] = 1;

                t[6] = 7;
                t[7] = 2;
                t[8] = 6;

                t[12] = 6;
                t[13] = 2;
                t[14] = 0;

                // Upper face
                t[3] = 10;
                t[4] = 9;
                t[5] = 4;

                t[9]  = 10;
                t[10] = 5;
                t[11] = 9;

                t[15] = 9;
                t[16] = 5;
                t[17] = 3;
            }
        }
        /* Case 5*/
        else if (topology[1] == 1 && topology[2] == 1)
        {
            if (edges[0] == 0) // If I colapse to the node 0 local
            {
                // Lower face
                t[0] = 5;
                t[1] = 4;
                t[2] = 2;

                t[6] = 5;
                t[7] = 0;
                t[8] = 4;

                t[12] = 4;
                t[13] = 0;
                t[14] = 1;
            }
            else if (edges[0] == 1)
            {
                // Lower face
                t[0] = 5;
                t[1] = 4;
                t[2] = 2;

                t[6] = 5;
                t[7] = 1;
                t[8] = 4;

                t[12] = 5;
                t[13] = 0;
                t[14] = 1;
            }
        }
        /* Case 6 */
        else if (topology[0] == 1 && topology[2] == 1)
        {
            if (edges[1] == 1)
            {
                // Lower face
                t[0] = 8;
                t[1] = 0;
                t[2] = 6;

                t[6] = 8;
                t[7] = 6;
                t[8] = 1;

                t[12] = 8;
                t[13] = 1;
                t[14] = 2;

                // Upper face
                t[3] = 11;
                t[4] = 3;
                t[5] = 9;

                t[9] = 11;
                t[10] = 9;
                t[11] = 4;

                t[15] = 11;
                t[14] = 4;
                t[15] = 5;
            }
            else if (edges[1] == 2)
            {
                // Lower face
                t[0] = 8;
                t[1] = 0;
                t[2] = 6;

                t[6] = 8;
                t[7] = 6;
                t[8] = 2;

                t[12] = 6;
                t[13] = 1;
                t[14] = 2;

                // Upper face
                t[3] = 11;
                t[4] = 3;
                t[5] = 9;

                t[9]  = 11;
                t[10] = 9;
                t[11] = 5;

                t[15] = 9;
                t[16] = 4;
                t[17] = 5;
            }
        }
        return 1;
    }
    else if (*splitted_edges == 6)
    {
        *number_elem = 4;
        // Lower face
        t[0] = 8;
        t[1] = 0;
        t[2] = 6;

        t[6] = 8;
        t[7] = 6;
        t[8] = 7;

        t[12] = 7;
        t[13] = 6;
        t[14] = 1;

        t[18] = 8;
        t[19] = 7;
        t[20] = 2;

        // Upper face
        t[3] = 11;
        t[4] = 3;
        t[5] = 9;

        t[9]  = 11;
        t[10] = 9;
        t[11] = 10;

        t[15] = 10;
        t[16] = 9;
        t[17] = 4;

        t[21] = 11;
        t[22] = 10;
        t[23] = 5;

        return 1;
    }
    else
    {
        return 0;
    }

}

#endif /* KRATOS_SPLIT_PRISM  defined */

