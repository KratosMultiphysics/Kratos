/*
==============================================================================
KratosIncompressibleFluidApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: antonia $
//   Date:                $Date: 2009-01-13 16:40:58 $
//   Revision:            $Revision: 1.24 $
//
//


#if !defined(KRATOS_ENRICHMENT_UTILITIES_INCLUDED )
#define  KRATOS_ENRICHMENT_UTILITIES_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <algorithm>
#include <limits>

// External includes 


// Project includes
#include "includes/define.h"
//#include "utilities/split_tetrahedra.c"


namespace Kratos
{

    /** This utility can be used to calculate the enriched shape function for tetrahedra element.
     *  The metodology consists in partitioning the tetrahedra in a set of sub-tetrahedra and
     *  cacluate the enrichment information using these partitions.
     */
    class EnrichmentUtilities
    {
    public:

        /**
         * The method to calculate the ernriched shape functions for given tetrahedra.
         * @param rPoints A 4x3 matrix where row i has the coordinates of node i.
         * @param DN_DX The gradient of the shape functions Ni respect to the reference coordinates
         * @param rDistances is an input  vector of 4 size which holds relative distance (not need to be exact) for each node.
         *        it is used internally to mark the position of the zero level
         * @param rVolumes Result vector with size 6 (maximumn number of partitions) holding the volume of each partition
         * @param rShapeFunctionValues Result 6x4 matrix where each row represents a partition and holds the shape functions N1 to N4
         *        of the original tetrahedra evaluated in the gauss point (center) of the partition.
         * @param rPartitionsSign A result vector of 6 holding the sign of the distance for the partition.
         *        The value -1 represents the negative distance sign, 1 represents positive distance and 0 stands for not used partition
         * @param rGradientsValue Restult vector of size 6 holding the gradient of the enriched shape funciton for each volume.
         *        Each element of vector is a 1x3 matrix representing the gradient of enriched shape function. The use of
         *        matrix is for possible future improvement.
         * @param Nenriched is a Matrix that contains for every gauss point the values of the enriched shape functions at the position of the gauss point
         *        so that Nenriched(1,0) contains the value of the enriched shape function "0" at the gauss point "1"
         * @return number of partitions created which can be from 1 to 6.
         *         1 holds for only 1 partition which is the original element. (No partitioning needed)
         */
        template<class TMatrixType, class TVectorType, class TGradientType>
        static int CalculateTetrahedraEnrichedShapeFuncions(TMatrixType const& rPoints, TGradientType const& DN_DX,
        TVectorType const& rDistances, TVectorType& rVolumes, TMatrixType& rShapeFunctionValues,
        TVectorType& rPartitionsSign, std::vector<TMatrixType>& rGradientsValue, TMatrixType& Nenriched)
        {
            KRATOS_TRY

                    const int n_nodes = 4; // it works only for tetrahedra
            const int n_edges = 6; // it works only for tetrahedra

            const int edge_i[] = {0, 0, 0, 1, 1, 2};
            const int edge_j[] = {1, 2, 3, 2, 3, 3};

            const int edges[4][4] = {
                {-1, 0, 1, 2},
                { 0, -1, 3, 4},
                { 1, 3, -1, 5},
                { 2, 4, 5, -1}};

            const double epsilon = 1.00e-9;
            const double near_factor = 1.00e-6;

            int number_of_partitions = 1;

            array_1d<double, n_edges> edges_dx; // It will be initialize later
            array_1d<double, n_edges> edges_dy; // It will be initialize later
            array_1d<double, n_edges> edges_dz; // It will be initialize later
            array_1d<double, n_edges> edges_length; // It will be initialize later
            // The divided part length from first node of edge respect to the edge length
            array_1d<double, n_edges> edge_division_i = ZeroVector(n_edges); // The 0 is for no split
            // The divided part length from second node of edge respect to the edge length
            array_1d<double, n_edges> edge_division_j = ZeroVector(n_edges); // The 0 is for no split

            int split_edge[] = {0, 1, 2, 3, -1, -1, -1, -1, -1, -1, -1};
            int new_node_id = 4;
            bounded_matrix<double, 4, 4 > length = ZeroMatrix(4, 4);

            int n_zero_distance_nodes = 0;
            int n_negative_distance_nodes = 0;
            int n_positive_distance_nodes = 0;
            int signs[] = {-2, -2, -2, -2};
            int zero_distance_nodes[] = {-1, -1, -1, -1};
            int negative_distance_nodes[] = {-1, -1, -1, -1};
            int positive_distance_nodes[] = {-1, -1, -1, -1};

            for (int i = 0; i < 6; i++)
                for (int j = 0; j < n_nodes; j++)
                    rShapeFunctionValues(i, j) = 0.25;



            //compute the gradient of the distance and normalize it
            array_1d<double, 3 > grad_d;
            noalias(grad_d) = prod(trans(DN_DX), rDistances);
            double norm = norm_2(grad_d);
            if (norm > epsilon)
                grad_d /= norm;

            array_1d<double, n_nodes> exact_distance = rDistances;
            array_1d<double, n_nodes> abs_distance = ZeroVector(n_nodes);




            for (int edge = 0; edge < n_edges; edge++)
            {
                const int i = edge_i[edge];
                const int j = edge_j[edge];

                double dx = rPoints(j, 0) - rPoints(i, 0);
                double dy = rPoints(j, 1) - rPoints(i, 1);
                double dz = rPoints(j, 2) - rPoints(i, 2);

                double l = sqrt(dx * dx + dy * dy + dz * dz);

                edges_dx[edge] = dx;
                edges_dy[edge] = dy;
                edges_dz[edge] = dz;
                edges_length[edge] = l;

                double relatively_near = l * near_factor;

                //               KRATOS_WATCH(i);
                //               KRATOS_WATCH(j);
                //               KRATOS_WATCH(l);
                //               KRATOS_WATCH(relatively_near);
                //               KRATOS_WATCH(rDistances[i]);
                //               KRATOS_WATCH(rDistances[j]);


                if (fabs(rDistances[i]) < relatively_near)
                {
                    signs[i] = 0;
                } else if (fabs(rDistances[j]) < relatively_near)
                {
                    signs[j] = 0;
                } else if (rDistances[i] * rDistances[j] < -epsilon)
                {
                    split_edge[edge + 4] = new_node_id++;
                    edge_division_i[edge] = fabs(rDistances[i]) / (fabs(rDistances[i]) + fabs(rDistances[j]));
                    edge_division_j[edge] = 1.00 - edge_division_i[edge];
                    const double d = fabs(edges_dx[edge] * grad_d[0] + edges_dy[edge] * grad_d[1] + edges_dz[edge] * grad_d[2]);

                    abs_distance[i] = d * edge_division_i[edge];
                    abs_distance[j] = d * edge_division_j[edge];
                }
            }

            for (int i_node = 0; i_node < n_nodes; i_node++)
            {
                if (signs[i_node] == 0)
                    zero_distance_nodes[n_zero_distance_nodes++] = i_node;
                else if (rDistances[i_node] < 0.00)
                {
                    signs[i_node] = -1;
                    negative_distance_nodes[n_negative_distance_nodes++] = i_node;
                } else
                {
                    signs[i_node] = 1;
                    positive_distance_nodes[n_positive_distance_nodes++] = i_node;
                }

            }

            //            exact_distance[0] =  (edges_dx[0] * grad_d[0] + edges_dy[0] * grad_d[1] + edges_dz[0] * grad_d[2]) * edge_division_i[0];
            //            abs_distance[0] = fabs(exact_distance[0]);

            for (int i = 0; i < n_nodes; i++)
            {

                if (rDistances[i] < 0)
                    exact_distance[i] = -abs_distance[i];
                else
                    exact_distance[i] = abs_distance[i];
                //exact_distance[i] =  (edges_dx[i] * grad_d[0] + edges_dy[i] * grad_d[1] + edges_dz[i] * grad_d[2]) * edge_division_j[i];
                //                abs_distance[i] = fabs(exact_distance[i]);
            }

            array_1d<double, 3 > exact_distance_gradient;
            noalias(exact_distance_gradient) = prod(trans(DN_DX), exact_distance);

            array_1d<double, 3 > abs_distance_gradient;
            noalias(abs_distance_gradient) = prod(trans(DN_DX), abs_distance);

            KRATOS_WATCH(rDistances);
            KRATOS_WATCH(exact_distance);
            KRATOS_WATCH(exact_distance_gradient);
            //            KRATOS_WATCH(grad_d);
            //            KRATOS_WATCH(edge_division_i)
            //            KRATOS_WATCH(edge_division_j);
            //            KRATOS_WATCH(n_zero_distance_nodes);
            //            KRATOS_WATCH(n_negative_distance_nodes);
            //            KRATOS_WATCH(n_positive_distance_nodes);
            //
            //            std::cout << "zero_distance_nodes : [" << zero_distance_nodes[0] << ", " <<  zero_distance_nodes[1] << ", "  << zero_distance_nodes[2] << ", " << zero_distance_nodes[3] << "]"<< std::endl;
            //            std::cout << "negative_distance_nodes : [" << negative_distance_nodes[0] << ", " <<  negative_distance_nodes[1] << ", "  << negative_distance_nodes[2] << ", " << negative_distance_nodes[3] << "]"<< std::endl;
            //            std::cout << "positive_distance_nodes : [" << positive_distance_nodes[0] << ", " <<  positive_distance_nodes[1] << ", "  << positive_distance_nodes[2] << ", " << positive_distance_nodes[3] << "]"<< std::endl;

            int number_of_splitted_edges = new_node_id - 4; //number of splitted edges

            double volume = edges_dx[0] * edges_dy[1] * edges_dz[2] -
                    edges_dx[0] * edges_dz[1] * edges_dy[2] +
                    edges_dy[0] * edges_dz[1] * edges_dx[2] -
                    edges_dy[0] * edges_dx[1] * edges_dz[2] +
                    edges_dz[0] * edges_dx[1] * edges_dy[2] -
                    edges_dz[0] * edges_dy[1] * edges_dx[2];

            const double one_sixth = 1.00 / 6.00;
            volume *= one_sixth;


            KRATOS_WATCH(volume)
            if (number_of_splitted_edges == 0) // no spliting
            {
                rVolumes[0] = volume;
                // Looking for the first node with sign not zero to get the sign of the element.
                for (int i_node = 0; i_node < n_nodes; i_node++)
                    if (signs[i_node] != 0)
                    {
                        rPartitionsSign[0] = signs[i_node];
                        break;
                    }
                number_of_partitions = 1; // There is only one partition, the element itself.
            } else if (number_of_splitted_edges == 1) // The surface passes between two nodes with zero distance
            {
                const int edge = edges[negative_distance_nodes[0]][positive_distance_nodes[0]];
                const int volume1_id = int(edge_i[edge] == positive_distance_nodes[0]);
                const int volume2_id = int(edge_j[edge] == positive_distance_nodes[0]);
                Divide1To2(edge_division_i, edge_division_j, edge, volume1_id, volume2_id, volume, rShapeFunctionValues, rVolumes);
                rPartitionsSign[0] = -1;
                rPartitionsSign[1] = 1;
                number_of_partitions = 2; // There are two partitions
            } else if (number_of_splitted_edges == 2) // The surface passes through a node with zero distance
            {
                // The node 0 is the node with zero distance, 
                // node1 is the one with sign different than the other two.
                int node1;
                int node2;
                int node3;
                if (n_negative_distance_nodes > 1)
                {
                    node1 = positive_distance_nodes[0];
                    node2 = negative_distance_nodes[0];
                    node3 = negative_distance_nodes[1];
                    rPartitionsSign[0] = -1;
                    rPartitionsSign[1] = -1;
                    rPartitionsSign[2] = 1;
                } else
                {
                    node1 = negative_distance_nodes[0];
                    node2 = positive_distance_nodes[0];
                    node3 = positive_distance_nodes[1];
                    rPartitionsSign[0] = 1;
                    rPartitionsSign[1] = 1;
                    rPartitionsSign[2] = -1;
                }
                KRATOS_WATCH(node1);
                KRATOS_WATCH(node2);
                // dividing volume to volume 0 and 1
                int edge = edges[node1][node2];
                int volume1_id = int(edge_i[edge] == node1);
                int volume2_id = int(edge_j[edge] == node1);
                Divide1To2(edge_division_i, edge_division_j, edge, volume1_id, volume2_id, volume, rShapeFunctionValues, rVolumes);

                // dividing volume 1 to volume 1 and 2
                edge = edges[node1][node3];
                volume1_id = 1 + int(edge_i[edge] == node1);
                volume2_id = 1 + int(edge_j[edge] == node1);
                volume = rVolumes[1];
                CopyShapeFunctionsValues(rShapeFunctionValues, 1, 2);
                Divide1To2(edge_division_i, edge_division_j, edge, volume1_id, volume2_id, volume, rShapeFunctionValues, rVolumes);

                number_of_partitions = 3; // There are three partitions
            } else if (number_of_splitted_edges == 3) // One node on one side the rest on the other side
            {

                // The node 0 is the node with different sign than the rest,
                // node1 is the one with sign different than the other two.
                int node0;
                int node1;
                int node2;
                int node3;
                if (n_negative_distance_nodes > 1)
                {
                    node0 = positive_distance_nodes[0];
                    node1 = negative_distance_nodes[0];
                    node2 = negative_distance_nodes[1];
                    node3 = negative_distance_nodes[2];
                    rPartitionsSign[0] = -1;
                    rPartitionsSign[1] = -1;
                    rPartitionsSign[2] = -1;
                    rPartitionsSign[3] = 1;
                } else
                {
                    node0 = negative_distance_nodes[0];
                    node1 = positive_distance_nodes[0];
                    node2 = positive_distance_nodes[1];
                    node3 = positive_distance_nodes[2];
                    rPartitionsSign[0] = 1;
                    rPartitionsSign[1] = 1;
                    rPartitionsSign[2] = 1;
                    rPartitionsSign[3] = -1;
                }

                // dividing volume to volume 0 and 1
                int edge = edges[node0][node1];
                int volume1_id = int(edge_i[edge] == node0);
                int volume2_id = int(edge_j[edge] == node0);

                Divide1To2(edge_division_i, edge_division_j, edge, volume1_id, volume2_id, volume, rShapeFunctionValues, rVolumes);

                // dividing volume 1 to volume 1 and 2
                edge = edges[node0][node2];
                volume1_id = 1 + int(edge_i[edge] == node0);
                volume2_id = 1 + int(edge_j[edge] == node0);
                volume = rVolumes[1];
                CopyShapeFunctionsValues(rShapeFunctionValues, 1, 2);
                Divide1To2(edge_division_i, edge_division_j, edge, volume1_id, volume2_id, volume, rShapeFunctionValues, rVolumes);

                // dividing volume 2 to volume 2 and 3
                edge = edges[node0][node3];
                volume1_id = 2 + int(edge_i[edge] == node0);
                volume2_id = 2 + int(edge_j[edge] == node0);
                volume = rVolumes[2];
                CopyShapeFunctionsValues(rShapeFunctionValues, 2, 3);
                Divide1To2(edge_division_i, edge_division_j, edge, volume1_id, volume2_id, volume, rShapeFunctionValues, rVolumes);
                number_of_partitions = 4; // There are four partitions
            } else if (number_of_splitted_edges == 4)
            {

                // The node 0 is the first negative node,
                // node1 is the one with sign different than the other two.
                int node0 = positive_distance_nodes[0];
                int node1 = positive_distance_nodes[1];
                int node2 = negative_distance_nodes[0];
                int node3 = negative_distance_nodes[1];

                rPartitionsSign[0] = 1;
                rPartitionsSign[1] = 1;
                rPartitionsSign[2] = 1;
                rPartitionsSign[3] = -1;
                rPartitionsSign[4] = -1;
                rPartitionsSign[5] = -1;



                // dividing volume  to volume 1 and 2
                int edge = edges[node0][node2];
                int volume1_id = int(edge_i[edge] == node2);
                int volume2_id = int(edge_j[edge] == node2);
                Divide1To2(edge_division_i, edge_division_j, edge, volume1_id, volume2_id, volume, rShapeFunctionValues, rVolumes);

                // dividing volume 0 to volume 0 and 2
                edge = edges[node0][node3];
                volume1_id = 2 * int(edge_i[edge] == node3);
                volume2_id = 2 * int(edge_j[edge] == node3);
                volume = rVolumes[0];
                CopyShapeFunctionsValues(rShapeFunctionValues, 0, 2);
                Divide1To2(edge_division_i, edge_division_j, edge, volume1_id, volume2_id, volume, rShapeFunctionValues, rVolumes);

                // dividing volume 1 to volume 1 and 3
                edge = edges[node1][node2];
                volume1_id = 1 + 2 * int(edge_i[edge] == node2);
                volume2_id = 1 + 2 * int(edge_j[edge] == node2);
                volume = rVolumes[1];
                CopyShapeFunctionsValues(rShapeFunctionValues, 1, 3);
                Divide1To2(edge_division_i, edge_division_j, edge, volume1_id, volume2_id, volume, rShapeFunctionValues, rVolumes);

                // dividing volume 1 to volume 1 and 4
                edge = edges[node1][node3];
                volume1_id = 1 + 3 * int(edge_i[edge] == node3);
                volume2_id = 1 + 3 * int(edge_j[edge] == node3);
                volume = rVolumes[1];
                CopyShapeFunctionsValues(rShapeFunctionValues, 1, 4);
                Divide1To2(edge_division_i, edge_division_j, edge, volume1_id, volume2_id, volume, rShapeFunctionValues, rVolumes);

                // dividing volume 2 to volume 2 and 5
                edge = edges[node1][node3];
                volume1_id = 2 + 3 * int(edge_i[edge] == node3);
                volume2_id = 2 + 3 * int(edge_j[edge] == node3);
                volume = rVolumes[2];
                CopyShapeFunctionsValues(rShapeFunctionValues, 2, 5);
                Divide1To2(edge_division_i, edge_division_j, edge, volume1_id, volume2_id, volume, rShapeFunctionValues, rVolumes);

                number_of_partitions = 6; // There are six partitions
            }


            for (int i = 0; i < number_of_partitions; i++)
            {
                //compute enriched shape function values
                double dist = 0.0;
                double abs_dist = 0.0
                for (int j = 0; j < 3; j++)
                {
                    dist += rShapeFunctionValues(i, j) * exact_distance;
                    abs_dist += rShapeFunctionValues(i, j) * abs_distance;
                }

                Nenriched(i, 0) = 0.5 * (rPartitionsSign[i] * abs_dist - dist);

                //compute shape function gradients
                for (int j = 0; j < 3; j++)
                    rGradientsValue[i](0, j) = 0.5 * (rPartitionsSign[i] * abs_distance_gradient[j] - exact_distance_gradient[j]);
            }

            /*
                        int edge_ids[6];
                        TetrahedraSplitMode(split_edge, edge_ids);
                        int n_new_elements; //number of elements generated
                        int n_splitted_edges; //number of splitted edges
                        int center_node_needed; //number of internal nodes
                        int t[56];
                        bool split_needed = Split_Tetrahedra(edge_ids, t, &n_new_elements, &n_splitted_edges, &center_node_needed);

                        if(center_node_needed == 1) // we need to generate a new internal node
                        {
                            //generate new node for the center
                            split_edge[10] = new_node_id++;
                        }



                        if (number_of_splitted_edges == 2) // The surface passes through a node with zero distance
                        {
                            int splited_edge = -1;
                            for (int i = 0; i < 6; i++)
                                for (int j = 0; j < n_nodes; j++)
                                    rShapeFunctionValues(i, j) = 0.25;
                
                            for (int i = 0; i < n_new_elements; i++) {
                                int nodes[]={0,0,0,0};
                                TetrahedraGetNewConnectivityGID(i, t, split_edge, nodes, nodes + 1, nodes + 2, nodes + 3);

                                std::cout << "element #" << i << " : " << nodes[0] << ", " << nodes[1] << ", " << nodes[2] << ", " << nodes[3] << std::endl;
                            }
                
                            for (int edge = 0; edge < n_edges; edge++)
                                {
                                if(edge_division_i[edge] != 0.00)
                                {
                                    splited_edge = edge;
                                    rVolumes[0] = edge_division_i[edge] * volume;
                                    rVolumes[1] = edge_division_j[edge] * volume;

                                    const int i = edge_i[edge];
                                    const int j = edge_j[edge];

                                    rShapeFunctionValues(0,i) *= edge_division_i[edge] * 0.25;
                                    rShapeFunctionValues(0,j) *= 0.5 - 0.25 * edge_division_i[edge];
                                    rShapeFunctionValues(1,i) *= 0.5 - 0.25 * edge_division_j[edge];
                                    rShapeFunctionValues(1,j) *= edge_division_j[edge] * 0.25;

                                    break; // There is only one splited edge
                                }
                            }

                            //return 3; // There are three partitions
                        }

                        std::cout << "split edge = ";

                        for(int i = 0 ; i < 11 ; i++)
                            std::cout << split_edge[i] << ",";

                        std::cout << std::endl;

                        std::cout << "edge ids = ";

                        for(int i = 0 ; i < 6 ; i++)
                            std::cout << edge_ids[i] << ",";

                        std::cout << std::endl;

                        KRATOS_WATCH(edge_division_i);
                        KRATOS_WATCH(n_new_elements);
                        KRATOS_WATCH(n_splitted_edges);
                        KRATOS_WATCH(center_node_needed);
                        KRATOS_WATCH(edges_length);

                        for (int i = 0; i < n_new_elements; i++) {
                            int i0, i1, i2, i3;
                            TetrahedraGetNewConnectivityGID(i, t, split_edge, &i0, &i1, &i2, &i3);

                            std::cout << "element #" << i << " : "  << i0 << ", "   << i1 << ", "   << i2 << ", "   << i3 << std::endl;
                        }

                        for(int i = 0 ; i < n_nodes ; i++)
                        {
                            for(int j = i ; j < n_nodes ; j++)
                            {
                                double dx=rPoints(i,0) - rPoints(j,0);
                                double dy=rPoints(i,1) - rPoints(j,1);
                                double dz=rPoints(i,2) - rPoints(j,2);
                                double l = sqrt(dx*dx+dy*dy+dz*dz);

                                length(i,j) = l;
                                length(j,i) = l;
                            }
                        }
             */
            //            KRATOS_WATCH(length)
            //
            //            for(int i = 0 ; i < n_nodes ; i++)
            //            {
            //                for(int j = i ; j < n_nodes ; j++)
            //                    if(fabs(rDistances[i]) < epsilon)
            //                    {
            //
            //                    }
            //            }

            //            for(int i = 0 ; i < n_nodes ; i++)
            //            {
            //                for(int j = i ; j < n_nodes ; j++)
            //                n_zero_distance_nodes += (fabs(rDistances[i]) < epsilon);
            //            }
            //
            //
            //            signs[0] = int(rDistances[0] < -epsilon);
            //
            //            for(int i = 1 ; i < n_nodes ; i++)
            //            {
            //                int different_sign = (rDistances[0]*rDistances[i] < -epsilon );
            //                signs[i]= different_sign;
            //                n_different_sign += different_sign;
            //            }



            //            KRATOS_WATCH(signs);
            //            KRATOS_WATCH(n_different_sign);
            //

            return number_of_partitions;
            KRATOS_CATCH("");
        }


    private:

        template<class TMatrixType>
        static void CopyShapeFunctionsValues(TMatrixType& rShapeFunctionValues, int OriginId, int DestinationId)
        {
            const int n_nodes = 4;
            for (int i = 0; i < n_nodes; i++)
                rShapeFunctionValues(DestinationId, i) = rShapeFunctionValues(OriginId, i);
        }

        template<class TMatrixType, class TVectorType>
        static void Divide1To2(array_1d<double, 6 > const& EdgeDivisionI, array_1d<double, 6 > const& EdgeDivisionJ, int Edge,
        int Volume1Id, int Volume2Id, double Volume, TMatrixType& rShapeFunctionValues, TVectorType& rVolumes)
        {
            const int edge_i[] = {0, 0, 0, 1, 1, 2};
            const int edge_j[] = {1, 2, 3, 2, 3, 3};

            const double division_i = EdgeDivisionI[Edge];
            const double division_j = EdgeDivisionJ[Edge];

            rVolumes[Volume1Id] = division_i * Volume;
            rVolumes[Volume2Id] = division_j * Volume;

            const int i = edge_i[Edge];
            const int j = edge_j[Edge];

            std::cout << "splitting edge" << i << " " << j << std::endl;
            KRATOS_WATCH(Volume1Id);
            KRATOS_WATCH(Volume2Id);
            KRATOS_WATCH(rShapeFunctionValues)

                    double delta1 = rShapeFunctionValues(Volume1Id, i) * (1.00 - division_i);
            rShapeFunctionValues(Volume1Id, i) += delta1;
            rShapeFunctionValues(Volume1Id, j) -= delta1;

            double delta2 = rShapeFunctionValues(Volume2Id, i) * (1.00 - division_j);
            rShapeFunctionValues(Volume2Id, j) += delta2;
            rShapeFunctionValues(Volume2Id, i) -= delta2;

            KRATOS_WATCH(delta1)
            KRATOS_WATCH(delta2)

            KRATOS_WATCH(rShapeFunctionValues)



                    //            rShapeFunctionValues(Volume1Id, i) += rShapeFunctionValues(Volume1Id, j) * (1.00 - division_i);
                    //            rShapeFunctionValues(Volume1Id, j) *= division_i;
                    //            rShapeFunctionValues(Volume2Id, j) += rShapeFunctionValues(Volume2Id, i) * (1.00 - division_j);
                    //            rShapeFunctionValues(Volume2Id, i) *= division_j;
                    //
                    //            rShapeFunctionValues(Volume1Id, i) = division_i * 0.25;
                    //            rShapeFunctionValues(Volume1Id, j) = 0.5 - 0.25 * division_i;
                    //            rShapeFunctionValues(Volume2Id, i) = 0.5 - 0.25 * division_j;
                    //            rShapeFunctionValues(Volume2Id, j) = division_j * 0.25;
        }


    };

} // namespace Kratos.

#endif // KRATOS_ENRICHMENT_UTILITIES_INCLUDED  defined


