#if !defined(KRATOS_SPLIT_TRIANGLE)
#define  KRATOS_SPLIT_TRIANGLE

/*VERSION 1.0 17 Feb 2011*/

/* Copyright (C) 2010 Riccardo Rossi, Pooyan Dadvand, Nelson Maireni
 Email contact: rrossi@cimne.upc.edu
 The current triangle splitting library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/
/** @file split_triangle.c
 * @brief The class contains three helper functions to ease the splitting: \n
 * TriangleSplitMode, Split_Triangle, and TriangleGetNewConnectivityGID\n
 * EXAMPLE: imagine that an user would like to split a triangle formed\n
 * by the ids 3 9 7 by introducing a new node 15 on the edge between 3 and 9\n
 * he should define\n
 * int aux[6]\n
 * int edge_ids[3]\n
 * int t[12]\n
 * then initialize\n
 * aux[0] = 3; aux[1] = 9; aux[2] = 7;\n
 * aux[3] = 15; //node on edge 01 --> edge to be refined\n
 * aux[4] = -1; //node on edge 12 -->edge not to be refined\n
 * aux[5] = -1; //node on edge 20 -->edge not to be refined\n
 * then call\n
 *
 * TriangleSplitMode(edge_ids,aux)\n
 * int nel; //number of nodes generated\n
 * int number_splitted_edges; //number of splitted edges\n
 * int nint; //number of internal nodes\n
 * bools split_needed = Split_Triangle(edge_ids,t, &nel, &number_splitted_edges, &nint)\n
 *
 * the new triangles ids can be then inspected by\n
 * for(int i=0; i<nel; i++)\n
 * {\n
 *     int i0,i1,i2;\n
 *     TriangleGetNewConnectivityGID(i, t, aux, &i0,&i1,&i2);\n
 * }\n
 */

/**
 * this function computes the splitting mode for the triangle\n
 *@param aux_ids contains a vector with the input Ids, organized as follows:\n
 *		aux_ids[0] = id of FIRST node of the original triangle\n
 *		aux_ids[1] = id of SECOND node of the original triangle\n
 *		aux_ids[2] = id of THIRD node of the original triangle\n
 *		aux_ids[4] = id of new node to be used for the edge 01 (-1 if edge not to be splitted)\n
 *		aux_ids[5] = id of new node to be used for the edge 12 (-1 if edge not to be splitted)\n
 *		aux_ids[6] = id of new node to be used for the edge 20 (-1 if edge not to be splitted)\n
 *		given this data it fills an auxiliary vector of size 3 that will be used in the splitting\n
 *@param edge_ids this is an auxiliary array with the local numbering. It is necessary for\n
 * 		the split_triangle function\n
 */
void TriangleSplitMode(const int aux_ids[6], int edge_ids[3]) {
    //edge 01
    if (aux_ids[3] < 0)
        if (aux_ids[0] > aux_ids[1]) edge_ids[0] = 0;
        else edge_ids[0] = 1;
    else
        edge_ids[0] = 3;

    //edge 12
    if (aux_ids[4] < 0)
        if (aux_ids[1] > aux_ids[2]) edge_ids[1] = 1;
        else edge_ids[1] = 2;
    else
        edge_ids[1] = 4;

    //edge 20
    if (aux_ids[5] < 0)
        if (aux_ids[2] > aux_ids[0]) edge_ids[2] = 2;
        else edge_ids[2] = 0;
    else
        edge_ids[2] = 5;
}

/**utility function to get the global ids for the new triangles to be generated\n
 *@param triangle_index --> the index of the new triangle to be generated\n
 *		(Should be less than the number nel provided by Split_Triangle)\n
 *@param t --> integer array provided by Split_Triangle\n
 *@param aux_ids --> array used in constructing the edge_ids (contains the Global Ids of the new nodes)\n
 *@param id0 --> Global ID of node0 of the new triangle\n
 *@param id1 --> Global ID of node1 of the new triangle\n
 *@param id2 --> Global ID of node2 of the new triangle\n
 */
inline void TriangleGetNewConnectivityGID(const int triangle_index,
        const int t[12],
        const int aux_ids[6],
        int* id0, int* id1, int* id2) {
    unsigned int base = triangle_index * 3;
    *id0 = aux_ids[t[base]];
    *id1 = aux_ids[t[base + 1]];
    *id2 = aux_ids[t[base + 2]];
}


///Utility to split triangles

/**
 * @param edges --> (input) int c array of size 3\n
 * @param t  --> (output) int c array of size 12 (3*4)\n
 * @param nel --> (output) number of elements in the subdivision\n
 * @param splitted_edges --> (output) provides the number of splitted edges\n
 * @param nint --> (output)  internal node (not needed for triangles)\n
 * @return true->splitting needed    false-->no splitting needed\n
 */
bool Split_Triangle(const int  edges[3], int t[12], int* nel, int* splitted_edges, int* nint) {
    *splitted_edges = 0;
    bool topology[3];
    topology[0] = false;
    topology[1] = false;
    topology[2] = false;
    for (unsigned int i = 0; i < 3; i++) {
        if (edges[i] > 2) {
            topology[i] = true;
            *splitted_edges = *splitted_edges + 1;
        }
    }

    if (*splitted_edges == 0 && *nint == 0) {
        //no splitting needed
        *nel = 1;
        t[0] = 0;
        t[1] = 1;
        t[2] = 2;
        return false;
    }
   //WARNING = case new central node needed
    else if (*splitted_edges == 0 && *nint == 1) {
        *nel = 3;

        t[0] = 3;
        t[1] = 0;
        t[2] = 1;

        t[3] = 3;
        t[4] = 1;
        t[5] = 2;

        t[5] = 3;
        t[5] = 2;
        t[5] = 0;
        return true;
    }

    else if (*splitted_edges == 1) {
        *nel = 2;
        // caso 1
        if (topology[0] == true) {
            t[0] = 3;
            t[1] = 2;
            t[2] = 0;

            t[3] = 3;
            t[4] = 1;
            t[5] = 2;
        }

            // caso 2
        else if (topology[1] == true) {
            t[0] = 4;
            t[1] = 0;
            t[2] = 1;

            t[3] = 4;
            t[4] = 2;
            t[5] = 0;
        }
            // caso 3
        else if (topology[2] == true) {
            t[0] = 5;
            t[1] = 1;
            t[2] = 2;

            t[3] = 5;
            t[4] = 0;
            t[5] = 1;
        }

        return true;

    }

    else if (*splitted_edges == 2) {
        *nel = 3;
        // caso 4
        if (topology[0] == true && topology[1] == true) {
            if (edges[2] == 0) // si colapso al nodo 0 local
            {
                t[0] = 4;
                t[1] = 3;
                t[2] = 1;

                t[3] = 4;
                t[4] = 0;
                t[5] = 3;

                t[6] = 4;
                t[7] = 2;
                t[8] = 0;
            }
            else if (edges[2] == 2) // si colapso al nodo 2 local
            {
                t[0] = 4;
                t[1] = 3;
                t[2] = 1;

                t[3] = 4;
                t[4] = 2;
                t[5] = 3;

                t[6] = 3;
                t[7] = 2;
                t[8] = 0;
            }

        }
            // caso 5
        else if (topology[1] == true && topology[2] == true) {
            if (edges[0] == 0) // si colapso al nodo 0 local
            {
                t[0] = 5;
                t[1] = 4;
                t[2] = 2;

                t[3] = 5;
                t[4] = 0;
                t[5] = 4;

                t[6] = 4;
                t[7] = 0;
                t[8] = 1;
            }
            else if (edges[0] == 1) // si colapso al nodo 2 local
            {
                t[0] = 5;
                t[1] = 4;
                t[2] = 2;

                t[3] = 5;
                t[4] = 1;
                t[5] = 4;

                t[6] = 5;
                t[7] = 0;
                t[8] = 1;
            }
        }


            /// caso 3
        else if (topology[0] == true && topology[2] == true) {
            if (edges[1] == 1) // si colapso al nodo 0 local
            {
                t[0] = 5;
                t[1] = 0;
                t[2] = 3;

                t[3] = 5;
                t[4] = 3;
                t[5] = 1;

                t[6] = 5;
                t[7] = 1;
                t[8] = 2;
            }
            else if (edges[1] == 2) // si colapso al nodo 2 local
            {
                t[0] = 5;
                t[1] = 0;
                t[2] = 3;

                t[3] = 5;
                t[4] = 3;
                t[5] = 2;

                t[6] = 3;
                t[7] = 1;
                t[8] = 2;
            }

        }

        return true;
    }
    else if (*splitted_edges == 3) {
        *nel = 4;
        t[0] = 5;
        t[1] = 0;
        t[2] = 3;

        t[3] = 5;
        t[4] = 3;
        t[5] = 4;

        t[6] = 4;
        t[7] = 3;
        t[8] = 1;

        t[9] = 5;
        t[10] = 4;
        t[11] = 2;

        return true;
    } else {
        return false;
    }

}

#endif // KRATOS_SPLIT_TRIANGLE  defined 

