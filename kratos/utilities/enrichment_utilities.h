//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Pablo Becker
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
#include "utilities/split_tetrahedra.h"


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
     *        so that it is  N(gauss_index, node_index)
     * @param rPartitionsSign A result vector of 6 holding the sign of the distance for the partition.
     *        The value -1 represents the negative distance sign, 1 represents positive distance and 0 stands for not used partition
     * @param rGradientsValue Restult vector of size 6 holding the gradient of the enriched shape funciton for each volume.
     *        Each element of vector is a 1x3 matrix representing the gradient of enriched shape function. The use of
     *        matrix is for possible future improvement.
     * @param Nenriched is a Matrix that contains for every gauss point the values of the enriched shape functions at the position of the gauss point
     *        so that Nenriched(1,0) contains the value of the enriched shape function "0" at the gauss point "1"
     * @param edge_areas an array of size 6 containing a positive double value for each edge. This value is zero if the edge is not cut, it corresponds to the
     *        area of the cut surface associtated to the edge otherwise
     * @return number of partitions created which can be from 1 to 6.
     *         1 holds for only 1 partition which is the original element. (No partitioning needed)
     */
    template<class TMatrixType, class TVectorType, class TGradientType>
    static int CalculateTetrahedraEnrichedShapeFuncions(TMatrixType const& rPoints,
            TGradientType const& DN_DX,
            TVectorType rDistances, TVectorType& rVolumes,
            TMatrixType& rShapeFunctionValues,
            TVectorType& rPartitionsSign,
            std::vector<TMatrixType>& rGradientsValue,
            TMatrixType& NEnriched,
            array_1d<double,6>& edge_areas
                                                       )
    {
        KRATOS_TRY

        const int n_nodes = 4; // it works only for tetrahedra
        const int n_edges = 6; // it works only for tetrahedra

        const int edge_i[] = {0, 0, 0, 1, 1, 2};
        const int edge_j[] = {1, 2, 3, 2, 3, 3};

// KRATOS_WATCH("                                                ");

        /*const int edges[4][4] = {
            {-1, 0, 1, 2},
            { 0, -1, 3, 4},
            { 1, 3, -1, 5},
            { 2, 4, 5, -1}
        };*/

        const double epsilon = 1e-15; //1.00e-9;
        const double near_factor = 1.00e-12;

        int number_of_partitions = 1;

        array_1d<double, n_edges> edges_dx; // It will be initialize later
        array_1d<double, n_edges> edges_dy; // It will be initialize later
        array_1d<double, n_edges> edges_dz; // It will be initialize later
        array_1d<double, n_edges> edges_length; // It will be initialize later
        // The divided part length from first node of edge respect to the edge length
        array_1d<double, n_edges> edge_division_i = ZeroVector(n_edges); // The 0 is for no split
        // The divided part length from second node of edge respect to the edge length
        array_1d<double, n_edges> edge_division_j = ZeroVector(n_edges); // The 0 is for no split

        BoundedMatrix<double, 8, 3 > aux_coordinates; //8 is the max number of nodes and aux_nodes
        for (unsigned int i = 0; i < 4; i++)
            for (unsigned int j = 0; j < 3; j++)
                aux_coordinates(i, j) = rPoints(i, j);
        for (unsigned int i = 4; i < 8; i++)
            for (unsigned int j = 0; j < 3; j++)
                aux_coordinates(i, j) = -10000.0; //set to a large number so that errors will be evident

        int split_edge[] = {0, 1, 2, 3, -1, -1, -1, -1, -1, -1, -1, -1};
        int new_node_id = 4;
        BoundedMatrix<double, 4, 4 > length = ZeroMatrix(4, 4);

        //int n_zero_distance_nodes = 0;
        int n_negative_distance_nodes = 0;
        int n_positive_distance_nodes = 0;
        array_1d<int,4> signs(4,-2);//[] = {-2, -2, -2, -2};
        //int zero_distance_nodes[] = {-1, -1, -1, -1};
        array_1d<int,4> negative_distance_nodes(4,-1);//[] = {-1, -1, -1, -1};
        array_1d<int,4> positive_distance_nodes(4,-1);//[] = {-1, -1, -1, -1};

        for (int i = 0; i < 6; i++)
            for (int j = 0; j < n_nodes; j++)
                rShapeFunctionValues(i, j) = 0.25;



        //compute the gradient of the distance and normalize it
        array_1d<double, 3 > grad_d;
        noalias(grad_d) = prod(trans(DN_DX), rDistances);
        double norm = norm_2(grad_d);
        if (norm > epsilon)
            grad_d /= (norm);

        array_1d<double, n_nodes> exact_distance = rDistances;
        array_1d<double, n_nodes> abs_distance = ZeroVector(n_nodes);
        double sub_volumes_sum = 0.00;

        //compute edge lenghts and max_lenght
        double max_lenght = 0.0;
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

            if(l > max_lenght)
                max_lenght = l;
        }

        double relatively_close = near_factor*max_lenght;
        array_1d<bool,4> collapsed_node;
        //identify collapsed nodes
        for (unsigned int i=0; i<4; i++)
        {
            if(fabs(rDistances[i]) < relatively_close)
            {
                collapsed_node[i] = true;
                rDistances[i] = 0.0;
// 		 KRATOS_WATCH("********************************!!!!!!!!!!!!!!!!!!!!!!!!!!!! collapsed node")
            }
            else
                collapsed_node[i] = false;

            abs_distance[i] = fabs(rDistances[i]);
        }

        //now decide splitting pattern
        for (int edge = 0; edge < n_edges; edge++)
        {
            const int i = edge_i[edge];
            const int j = edge_j[edge];
            if (rDistances[i] * rDistances[j] < 0.0)
            {
                const double tmp = fabs(rDistances[i]) / (fabs(rDistances[i]) + fabs(rDistances[j]));
// 		    const double d = fabs(edges_dx[edge] * grad_d[0] + edges_dy[edge] * grad_d[1] + edges_dz[edge] * grad_d[2]);
// 		    abs_distance[i] = d * tmp;
// 		    abs_distance[j] = d * (1.0 - tmp);

                if (collapsed_node[i] == false && collapsed_node[j] == false)
                {
                    split_edge[edge + 4] = new_node_id;
                    edge_division_i[edge] = tmp;
                    edge_division_j[edge] = 1.00 - tmp;

                    //compute the position of the edge node
                    for (unsigned int k = 0; k < 3; k++)
                        aux_coordinates(new_node_id, k) = rPoints(i, k) * edge_division_j[edge] + rPoints(j, k) * edge_division_i[edge];

                    new_node_id++;
                }
// 		    else
// 		    {
// 		       if(collapsed_node[i] == true) split_edge[edge + 4] = i;
// 		       else split_edge[edge + 4] = j;
// 		    }
            }
        }

        //compute the abs exact distance for all of the nodes
        if(new_node_id > 4) //at least one edge is cut
        {
            array_1d<double,3> base_point;
            base_point[0] = aux_coordinates(4,0);
            base_point[1] = aux_coordinates(4,1);
            base_point[2] = aux_coordinates(4,2);


            for (int i_node = 0; i_node < n_nodes; i_node++)
            {
                double d =    (rPoints(i_node,0) - base_point[0]) * grad_d[0] +
                              (rPoints(i_node,1) - base_point[1]) * grad_d[1] +
                              (rPoints(i_node,2) - base_point[2]) * grad_d[2] ;
                abs_distance[i_node] = fabs(d);
            }

        }


        for (int i_node = 0; i_node < n_nodes; i_node++)
        {
//                 if (collapsed_node[i_node] == true)
// 		{
// 		    abs_distance[i_node] = 0.0;
// 		    signs[i_node] = 1;
// 		    positive_distance_nodes[n_negative_distance_nodes++] = i_node;
// //                     zero_distance_nodes[n_zero_distance_nodes++] = i_node;
// 		}
//                 else
            if (rDistances[i_node] < 0.00)
            {
                signs[i_node] = -1;
                negative_distance_nodes[n_negative_distance_nodes++] = i_node;
            }
            else
            {
                signs[i_node] = 1;
                positive_distance_nodes[n_positive_distance_nodes++] = i_node;
            }
        }

        //assign correct sign to exact distance
        for (int i = 0; i < n_nodes; i++)
        {
            if (rDistances[i] < 0.0)
                exact_distance[i] = -abs_distance[i];
            else
                exact_distance[i] = abs_distance[i];
        }

        //compute exact distance gradients
        array_1d<double, 3 > exact_distance_gradient = prod(trans(DN_DX), exact_distance);

        array_1d<double, 3 > abs_distance_gradient = prod(trans(DN_DX), abs_distance);

        int number_of_splitted_edges = new_node_id - 4; //number of splitted edges

        double volume = edges_dx[0] * edges_dy[1] * edges_dz[2] -
                        edges_dx[0] * edges_dz[1] * edges_dy[2] +
                        edges_dy[0] * edges_dz[1] * edges_dx[2] -
                        edges_dy[0] * edges_dx[1] * edges_dz[2] +
                        edges_dz[0] * edges_dx[1] * edges_dy[2] -
                        edges_dz[0] * edges_dy[1] * edges_dx[2];

        const double one_sixth = 1.00 / 6.00;
        volume *= one_sixth;


        //            KRATOS_WATCH(volume)
        if (number_of_splitted_edges == 0) // no splitting
        {

            rVolumes[0] = volume;
            sub_volumes_sum = volume;
//                 // Looking for the first node with sign not zero to get the sign of the element.
//                 for (int i_node = 0; i_node < n_nodes; i_node++)
//                     if (signs[i_node] != 0) {
//                         rPartitionsSign[0] = signs[i_node];
//                         break;
//                     }
            //take the sign from the node with min distance
            double min_distance = 1e9;
            for (int j = 0; j < 4; j++)
                if(exact_distance[j] < min_distance) min_distance = exact_distance[j];

            if(min_distance < 0.0)
                rPartitionsSign[0] = -1.0;
            else
                rPartitionsSign[0] = 1.0;

            number_of_partitions = 1;

            for (int j = 0; j < 4; j++)
                rShapeFunctionValues(0, j) = 0.25;

            for (int j = 0; j < number_of_partitions; j++)
                NEnriched(j, 0) = 0.0;

            rGradientsValue[0] = ZeroMatrix(1,3);
        }
        else //if (number_of_splitted_edges == 4)
        {
            //define the splitting mode for the tetrahedra
            array_1d<double,4> local_subtet_indices;

            int edge_ids[6];
            TetrahedraSplit::TetrahedraSplitMode(split_edge, edge_ids);
            int nel = 0; //number of elements generated
            int n_splitted_edges; //number of splitted edges
            int nint; //number of internal nodes
            int t[56];
            TetrahedraSplit::Split_Tetrahedra(edge_ids, t, &nel, &n_splitted_edges, &nint);



            if (nint != 0)
                KRATOS_THROW_ERROR(std::logic_error, "requiring an internal node for splitting ... can not accept this", "");


            //now obtain the tetras and compute their center coordinates and volume
            noalias(edge_areas) = ZeroVector(6);
            array_1d<double, 3 > center_position;
            for (int i = 0; i < nel; i++)
            {
                int i0, i1, i2, i3; //indices of the subtetrahedra
                TetrahedraSplit::TetrahedraGetNewConnectivityGID(i, t, split_edge, &i0, &i1, &i2, &i3);


                double sub_volume = ComputeSubTetraVolumeAndCenter(aux_coordinates, center_position, i0, i1, i2, i3);


                //here we add to the areas of the cut adges
                local_subtet_indices[0] = t[i*4];
                local_subtet_indices[1] = t[i*4+1];
                local_subtet_indices[2] = t[i*4+2];
                local_subtet_indices[3] = t[i*4+3];
                AddToEdgeAreas<3>(edge_areas,exact_distance,local_subtet_indices,sub_volume);


                rVolumes[i] = sub_volume;

                sub_volumes_sum += sub_volume;

                array_1d<double, 4 > N;
                ComputeElementCoordinates(N, center_position, rPoints, volume);
                for (int j = 0; j < 4; j++)
                    rShapeFunctionValues(i, j) = N[j];
            }

            number_of_partitions = nel;

        }

//             if(fabs(sub_volumes_sum/volume - 1.0) > 1e-9)
// 	    {
// 	      KRATOS_WATCH(volume);
// 	      KRATOS_WATCH(rVolumes);
// 	      KRATOS_WATCH(sub_volumes_sum);
// 	      KRATOS_THROW_ERROR(std::logic_error,"the elemental volume does not match the sum of the sub volumes","")
// 	    }
// KRATOS_WATCH(exact_distance);
// KRATOS_WATCH(abs_distance);

        //double verify_volume = 0.0;
        if (number_of_partitions > 1)   // we won't calculate the N and its gradients for element without partitions
        {

            //compute the maximum absolute distance on the cut so to normalize the shape functions
            //now decide splitting pattern
            double max_aux_dist_on_cut = -1;
            for (int edge = 0; edge < n_edges; edge++)
            {
                const int i = edge_i[edge];
                const int j = edge_j[edge];
                if (rDistances[i] * rDistances[j] < 0.0)
                {
                    const double tmp = fabs(rDistances[i]) / (fabs(rDistances[i]) + fabs(rDistances[j]));

                    //compute the position of the edge node
                    double abs_dist_on_cut = abs_distance[i] * tmp + abs_distance[j] * (1.00 - tmp);

                    if(abs_dist_on_cut > max_aux_dist_on_cut) max_aux_dist_on_cut = abs_dist_on_cut;

                }
            }

            /*		if(max_aux_dist_on_cut < 1e-10)
            		  max_aux_dist_on_cut = 1e-10;*/
            if(max_aux_dist_on_cut < 1e-9*max_lenght)
                max_aux_dist_on_cut =  1e-9*max_lenght;

            for (int i = 0; i < number_of_partitions; i++)
            {
                //compute enriched shape function values
                double dist = 0.0;
                double abs_dist = 0.0;
                for (int j = 0; j < 4; j++)
                {
                    dist += rShapeFunctionValues(i, j) * exact_distance[j];
                    abs_dist += rShapeFunctionValues(i, j) * abs_distance[j];
                }

                if (dist < 0.0)
                    rPartitionsSign[i] = -1.0;
                else
                    rPartitionsSign[i] = 1.0;

                NEnriched(i, 0) = 0.5 * (abs_dist - rPartitionsSign[i] * dist);

                //normalizing
//                 NEnriched(i, 0) /= max_aux_dist_on_cut;
                /*KRATOS_WATCH(abs_dist);
                KRATOS_WATCH(dist);
                KRATOS_WATCH(rPartitionsSign);
                KRATOS_WATCH(abs_distance_gradient);
                KRATOS_WATCH(exact_distance_gradient);*/
                //compute shape function gradients
                for (int j = 0; j < 3; j++)
                {
                    rGradientsValue[i](0, j) = (0.5/*/max_aux_dist_on_cut*/) * (abs_distance_gradient[j] - rPartitionsSign[i] * exact_distance_gradient[j]);
                }


            }
        }
        else
        {
            NEnriched(0,0) = 0.0;

            for (int j = 0; j < 3; j++)
                rGradientsValue[0](0, j) = 0.0;
        }

// 	    for (unsigned int i=0; i<4; i++)
// 	    {
// 	       if( collapsed_node[i] == true)
// 	       {
// 		 KRATOS_WATCH(number_of_partitions);
// 		 for (int k = 0; k < number_of_partitions; k++)
// 		 {
// 		 KRATOS_WATCH(rVolumes[k]);
//
// 		 KRATOS_WATCH(NEnriched(k,0));
// 		 KRATOS_WATCH(rGradientsValue[k]);
// 		 }
// 	       }
// 	    }





        return number_of_partitions;
        KRATOS_CATCH("");
    }

    /**
     * The method to calculate 2 ernrichment shape functions for given tetrahedra.
     * the first one is the gradient, just like the previous function, and the second is a jump.
     * another difference with the previous function is that fixed size arrays and matrices are used if possible.
     * @param rPoints A 4x3 matrix where row i has the coordinates of node i.
     * @param DN_DX The gradient of the shape functions Ni respect to the reference coordinates
     * @param rDistances is an input  vector of 4 size which holds relative distance (not need to be exact) for each node.
     *        it is used internally to mark the position of the zero level
     * @param rVolumes Result vector with size 6 (maximumn number of partitions) holding the volume of each partition
     * @param rShapeFunctionValues Result 6x4 matrix where each row represents a partition and holds the shape functions N1 to N4
     *        of the original tetrahedra evaluated in the gauss point (center) of the partition.
     *        so that it is  N(gauss_index, node_index)
     * @param rPartitionsSign A result vector of 6 holding the sign of the distance for the partition.
     *        The value -1 represents the negative distance sign, 1 represents positive distance and 0 stands for not used partition
     * @param rGradientsValue Restult vector of size 6 holding the gradient of the enriched shape funcitons for each volume.
     *        Each element of vector is a 2x3 matrix representing the gradients of enriched shape functions. The row indicates enrichemnt function
     *        and the colum the direction of the gradient
     * @param Nenriched is a Matrix of size 6x2 that contains for every gauss point the values of the enriched shape functions at the position of the gauss point
     *        so that Nenriched(1,0) contains the value of the enriched shape function "0" at the gauss point "1"
     * @return number of partitions created which can be from 1 to 6.
     *         1 holds for only 1 partition which is the original element. (No partitioning needed)
     */


    //3d: (tetrahedrons)
    static int CalculateEnrichedShapeFuncions(BoundedMatrix<double,4, 3 >& rPoints,
            BoundedMatrix<double, 4, 3 >& DN_DX,
            array_1d<double,4>& rDistances,
            array_1d<double,6>& rVolumes,
            BoundedMatrix<double, 6, 4 >& rShapeFunctionValues,
            array_1d<double,6>& rPartitionsSign,
            std::vector<Matrix>& rGradientsValue,
            BoundedMatrix<double,6, 2>& NEnriched)
    {
        KRATOS_TRY

		//the first thing we do is checking if the rGradientsValue has the right lenght
		if (rGradientsValue.size()!=6)
			rGradientsValue.resize(6);
		// now we check if matrices contained in rGradientsValue have the right size:
		for (unsigned int i = 0; i < 6; i++)
		{
			if (rGradientsValue[i].size1() != 2 || rGradientsValue[i].size2() != 3)
			{
				std::cout << "resizing rGradientsValue[i] to 2x3, please modify your element" << std::endl;
				rGradientsValue[i].resize(2, 3, false);  //2 values of the 2 shape functions, and derivates in (xyz) direction).
			}
		}

        const int n_nodes = 4; // it works only for tetrahedra
        const int n_edges = 6; // it works only for tetrahedra

        const int edge_i[] = {0, 0, 0, 1, 1, 2};
        const int edge_j[] = {1, 2, 3, 2, 3, 3};

// KRATOS_WATCH("                                                ");

        /*const int edges[4][4] = {
            {-1, 0, 1, 2},
            { 0, -1, 3, 4},
            { 1, 3, -1, 5},
            { 2, 4, 5, -1}
        };*/

        const double epsilon = 1e-15; //1.00e-9;
        const double near_factor = 1.00e-12;

        int number_of_partitions = 1;

        array_1d<double, n_edges> edges_dx; // It will be initialize later
        array_1d<double, n_edges> edges_dy; // It will be initialize later
        array_1d<double, n_edges> edges_dz; // It will be initialize later
        array_1d<double, n_edges> edges_length; // It will be initialize later
        // The divided part length from first node of edge respect to the edge length
        array_1d<double, n_edges> edge_division_i = ZeroVector(n_edges); // The 0 is for no split
        // The divided part length from second node of edge respect to the edge length
        array_1d<double, n_edges> edge_division_j = ZeroVector(n_edges); // The 0 is for no split

        BoundedMatrix<double, 8, 3 > aux_coordinates; //8 is the max number of nodes and aux_nodes
        for (unsigned int i = 0; i < 4; i++)
            for (unsigned int j = 0; j < 3; j++)
                aux_coordinates(i, j) = rPoints(i, j);
        for (unsigned int i = 4; i < 8; i++)
            for (unsigned int j = 0; j < 3; j++)
                aux_coordinates(i, j) = -10000.0; //set to a large number so that errors will be evident

        int split_edge[] = {0, 1, 2, 3, -1, -1, -1, -1, -1, -1, -1, -1};
        int new_node_id = 4;
        BoundedMatrix<double, 4, 4 > length = ZeroMatrix(4, 4);

        //int n_zero_distance_nodes = 0;
        int n_negative_distance_nodes = 0;
        int n_positive_distance_nodes = 0;
        array_1d<int,4> signs(4,-2);//[] = {-2, -2, -2, -2};
        //int zero_distance_nodes[] = {-1, -1, -1, -1};
        array_1d<int,4> negative_distance_nodes(4,-1);//[] = {-1, -1, -1, -1};
        array_1d<int,4> positive_distance_nodes(4,-1);//[] = {-1, -1, -1, -1};

        for (int i = 0; i < 6; i++)
            for (int j = 0; j < n_nodes; j++)
                rShapeFunctionValues(i, j) = 0.25;



        //compute the gradient of the distance and normalize it
        array_1d<double, 3 > grad_d;
        noalias(grad_d) = prod(trans(DN_DX), rDistances);
        double norm = norm_2(grad_d);
        if (norm > epsilon)
            grad_d /= (norm);

        array_1d<double, n_nodes> exact_distance = rDistances;
        array_1d<double, n_nodes> abs_distance = ZeroVector(n_nodes);
        double sub_volumes_sum = 0.00;

        //compute edge lenghts and max_lenght
        double max_lenght = 0.0;
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

            if(l > max_lenght)
                max_lenght = l;
        }

        double relatively_close = near_factor*max_lenght;
        array_1d<bool,4> collapsed_node;
        //identify collapsed nodes
        for (unsigned int i=0; i<4; i++)
        {
            if(fabs(rDistances[i]) < relatively_close)
            {
                collapsed_node[i] = true;
                rDistances[i] = 0.0;
// 		 KRATOS_WATCH("********************************!!!!!!!!!!!!!!!!!!!!!!!!!!!! collapsed node")
            }
            else
                collapsed_node[i] = false;

            abs_distance[i] = fabs(rDistances[i]);
        }

        //now decide splitting pattern
        for (int edge = 0; edge < n_edges; edge++)
        {
            const int i = edge_i[edge];
            const int j = edge_j[edge];
            if (rDistances[i] * rDistances[j] < 0.0)
            {
                const double tmp = fabs(rDistances[i]) / (fabs(rDistances[i]) + fabs(rDistances[j]));
// 		    const double d = fabs(edges_dx[edge] * grad_d[0] + edges_dy[edge] * grad_d[1] + edges_dz[edge] * grad_d[2]);
// 		    abs_distance[i] = d * tmp;
// 		    abs_distance[j] = d * (1.0 - tmp);

                if (collapsed_node[i] == false && collapsed_node[j] == false)
                {
                    split_edge[edge + 4] = new_node_id;
                    edge_division_i[edge] = tmp;
                    edge_division_j[edge] = 1.00 - tmp;

                    //compute the position of the edge node
                    for (unsigned int k = 0; k < 3; k++)
                        aux_coordinates(new_node_id, k) = rPoints(i, k) * edge_division_j[edge] + rPoints(j, k) * edge_division_i[edge];

                    new_node_id++;
                }
// 		    else
// 		    {
// 		       if(collapsed_node[i] == true) split_edge[edge + 4] = i;
// 		       else split_edge[edge + 4] = j;
// 		    }
            }
        }

        //compute the abs exact distance for all of the nodes
        if(new_node_id > 4) //at least one edge is cut
        {
            array_1d<double,3> base_point;
            base_point[0] = aux_coordinates(4,0);
            base_point[1] = aux_coordinates(4,1);
            base_point[2] = aux_coordinates(4,2);


            for (int i_node = 0; i_node < n_nodes; i_node++)
            {
                double d =    (rPoints(i_node,0) - base_point[0]) * grad_d[0] +
                              (rPoints(i_node,1) - base_point[1]) * grad_d[1] +
                              (rPoints(i_node,2) - base_point[2]) * grad_d[2] ;
                abs_distance[i_node] = fabs(d);
            }

        }


        for (int i_node = 0; i_node < n_nodes; i_node++)
        {
//                 if (collapsed_node[i_node] == true)
// 		{
// 		    abs_distance[i_node] = 0.0;
// 		    signs[i_node] = 1;
// 		    positive_distance_nodes[n_negative_distance_nodes++] = i_node;
// //                     zero_distance_nodes[n_zero_distance_nodes++] = i_node;
// 		}
//                 else
            if (rDistances[i_node] < 0.00)
            {
                signs[i_node] = -1;
                negative_distance_nodes[n_negative_distance_nodes++] = i_node;
            }
            else
            {
                signs[i_node] = 1;
                positive_distance_nodes[n_positive_distance_nodes++] = i_node;
            }
        }

        //assign correct sign to exact distance
        for (int i = 0; i < n_nodes; i++)
        {
            if (rDistances[i] < 0.0)
                exact_distance[i] = -abs_distance[i];
            else
                exact_distance[i] = abs_distance[i];
        }

        //compute exact distance gradients
        array_1d<double, 3 > exact_distance_gradient;
        noalias(exact_distance_gradient) = prod(trans(DN_DX), exact_distance);

        array_1d<double, 3 > abs_distance_gradient;
        noalias(abs_distance_gradient) = prod(trans(DN_DX), abs_distance);

        int number_of_splitted_edges = new_node_id - 4; //number of splitted edges

        double volume = edges_dx[0] * edges_dy[1] * edges_dz[2] -
                        edges_dx[0] * edges_dz[1] * edges_dy[2] +
                        edges_dy[0] * edges_dz[1] * edges_dx[2] -
                        edges_dy[0] * edges_dx[1] * edges_dz[2] +
                        edges_dz[0] * edges_dx[1] * edges_dy[2] -
                        edges_dz[0] * edges_dy[1] * edges_dx[2];

        const double one_sixth = 1.00 / 6.00;
        volume *= one_sixth;


        //            KRATOS_WATCH(volume)
        if (number_of_splitted_edges == 0) // no splitting
        {

            rVolumes[0] = volume;
            sub_volumes_sum = volume;
//                 // Looking for the first node with sign not zero to get the sign of the element.
//                 for (int i_node = 0; i_node < n_nodes; i_node++)
//                     if (signs[i_node] != 0) {
//                         rPartitionsSign[0] = signs[i_node];
//                         break;
//                     }
            //take the sign from the node with min distance
            double min_distance = 1e9;
            for (int j = 0; j < 4; j++)
                if(exact_distance[j] < min_distance) min_distance = exact_distance[j];

            if(min_distance < 0.0)
                rPartitionsSign[0] = -1.0;
            else
                rPartitionsSign[0] = 1.0;

            number_of_partitions = 1;

            for (int j = 0; j < 4; j++)
                rShapeFunctionValues(0, j) = 0.25;

            for (int j = 0; j < number_of_partitions; j++)
                NEnriched(j, 0) = 0.0;

            rGradientsValue[0] = ZeroMatrix(1,3);
        }
        else //if (number_of_splitted_edges == 4)
        {
            //define the splitting mode for the tetrahedra
            int edge_ids[6];
            TetrahedraSplit::TetrahedraSplitMode(split_edge, edge_ids);
            int nel=0; //number of elements generated
            int n_splitted_edges; //number of splitted edges
            int nint; //number of internal nodes
            int t[56];
            TetrahedraSplit::Split_Tetrahedra(edge_ids, t, &nel, &n_splitted_edges, &nint);



            if (nint != 0)
                KRATOS_THROW_ERROR(std::logic_error, "requiring an internal node for splitting ... can not accept this", "");


            //now obtain the tetras and compute their center coordinates and volume
            array_1d<double, 3 > center_position;
            for (int i = 0; i < nel; i++)
            {
                int i0, i1, i2, i3; //indices of the subtetrahedra
                TetrahedraSplit::TetrahedraGetNewConnectivityGID(i, t, split_edge, &i0, &i1, &i2, &i3);


                double sub_volume = ComputeSubTetraVolumeAndCenter(aux_coordinates, center_position, i0, i1, i2, i3);

                rVolumes[i] = sub_volume;

                sub_volumes_sum += sub_volume;

                array_1d<double, 4 > N;
                ComputeElementCoordinates(N, center_position, rPoints, volume);
                for (int j = 0; j < 4; j++)
                    rShapeFunctionValues(i, j) = N[j];
            }

            number_of_partitions = nel;

        }

//             if(fabs(sub_volumes_sum/volume - 1.0) > 1e-9)
// 	    {
// 	      KRATOS_WATCH(volume);
// 	      KRATOS_WATCH(rVolumes);
// 	      KRATOS_WATCH(sub_volumes_sum);
// 	      KRATOS_THROW_ERROR(std::logic_error,"the elemental volume does not match the sum of the sub volumes","")
// 	    }
// KRATOS_WATCH(exact_distance);
// KRATOS_WATCH(abs_distance);

        //double verify_volume = 0.0;
        if (number_of_partitions > 1)   // we won't calculate the N and its gradients for element without partitions
        {

            //compute the maximum absolute distance on the cut so to normalize the shape functions
            //now decide splitting pattern
            double max_aux_dist_on_cut = -1;
            for (int edge = 0; edge < n_edges; edge++)
            {
                const int i = edge_i[edge];
                const int j = edge_j[edge];
                if (rDistances[i] * rDistances[j] < 0.0)
                {
                    const double tmp = fabs(rDistances[i]) / (fabs(rDistances[i]) + fabs(rDistances[j]));

                    //compute the position of the edge node
                    double abs_dist_on_cut = abs_distance[i] * tmp + abs_distance[j] * (1.00 - tmp);

                    if(abs_dist_on_cut > max_aux_dist_on_cut) max_aux_dist_on_cut = abs_dist_on_cut;

                }
            }

            /*		if(max_aux_dist_on_cut < 1e-10)
            		  max_aux_dist_on_cut = 1e-10;*/
            if(max_aux_dist_on_cut < 1e-9*max_lenght)
                max_aux_dist_on_cut =  1e-9*max_lenght;

            for (int i = 0; i < number_of_partitions; i++)
            {
                //compute enriched shape function values
                double dist = 0.0;
                double abs_dist = 0.0;
                for (int j = 0; j < 4; j++)
                {
                    dist += rShapeFunctionValues(i, j) * exact_distance[j];
                    abs_dist += rShapeFunctionValues(i, j) * abs_distance[j];
                }

                if (dist < 0.0)
                    rPartitionsSign[i] = -1.0;
                else
                    rPartitionsSign[i] = 1.0;

                NEnriched(i, 0) = 0.5 * (abs_dist - rPartitionsSign[i] * dist);

                //normalizing
                NEnriched(i, 0) /= max_aux_dist_on_cut;
                NEnriched(i, 1) = (rPartitionsSign[i])*NEnriched(i, 0);
                /*KRATOS_WATCH(abs_dist);
                KRATOS_WATCH(dist);
                KRATOS_WATCH(rPartitionsSign);
                KRATOS_WATCH(abs_distance_gradient);
                KRATOS_WATCH(exact_distance_gradient);*/
                //compute shape function gradients
                for (int j = 0; j < 3; j++)
                {
                    rGradientsValue[i](0, j) = (0.5/max_aux_dist_on_cut) * (abs_distance_gradient[j] - rPartitionsSign[i] * exact_distance_gradient[j]);
                    rGradientsValue[i](1, j) = (rPartitionsSign[i])* rGradientsValue[i](0, j);
                }


            }
        }
        else
        {
            NEnriched(0,0) = 0.0;

            for (int j = 0; j < 3; j++)
                rGradientsValue[0](0, j) = 0.0;
        }

// 	    for (unsigned int i=0; i<4; i++)
// 	    {
// 	       if( collapsed_node[i] == true)
// 	       {
// 		 KRATOS_WATCH(number_of_partitions);
// 		 for (int k = 0; k < number_of_partitions; k++)
// 		 {
// 		 KRATOS_WATCH(rVolumes[k]);
//
// 		 KRATOS_WATCH(NEnriched(k,0));
// 		 KRATOS_WATCH(rGradientsValue[k]);
// 		 }
// 	       }
// 	    }





        return number_of_partitions;
        KRATOS_CATCH("");
    }

    //2D
    static int CalculateEnrichedShapeFuncions(
        BoundedMatrix<double,3, 2 >& rPoints,
                BoundedMatrix<double, 3, 2 >& DN_DX,
                 array_1d<double,3>& rDistances,
                 array_1d<double,3>& rVolumes,
                 BoundedMatrix<double, 3, 3 >& rGPShapeFunctionValues,
                 array_1d<double,3>& rPartitionsSign,
                 std::vector<Matrix>& rGradientsValue,
                 BoundedMatrix<double,3, 2>& NEnriched)
    {
        KRATOS_TRY

        const double one_third=1.0/3.0;
        BoundedMatrix<double,3,2> aux_points; //for auxiliary nodes 4(between 1 and 2) ,5(between 2 and 3) ,6 (between 3 and 1)
        BoundedMatrix<double, 3, 2 > coord_subdomain; //used to pass arguments when we must calculate areas, shape functions, etc
        BoundedMatrix<double,3,2> DN_DX_subdomain; //used to retrieve derivatives

        double most_common_sign=0; //the side of the cut in which two nodes are found (same sign) will be the ones that remains unchanged when builing the discontinuity
        double Area;//area of the complete element
        rGPShapeFunctionValues(0,0)=one_third;
        rGPShapeFunctionValues(0,1)=one_third;
        rGPShapeFunctionValues(0,2)=one_third; //default, when no interfase has been found
        Area = CalculateVolume2D( rPoints );
        array_1d<bool,3> cut_edges;
        array_1d<double,3> aux_nodes_relative_locations;
        BoundedMatrix<int,3,2> aux_nodes_father_nodes;

        //to begin with we must check whether our element is cut or not by the interfase.
        if( (rDistances(0)*rDistances(1))>0.0 && (rDistances(0)*rDistances(2))>0.0 ) //it means that this element IS NOT cut by the interfase. we must return data of a normal, non-enriched element
        {
            rVolumes(0)=Area;
            rGPShapeFunctionValues(0,0)=one_third;
            rGPShapeFunctionValues(0,1)=one_third;
            rGPShapeFunctionValues(0,2)=one_third;
            NEnriched(0,0) = 0.0;
            //type_of_cut=1;
            for (int j = 0; j < 3; j++)
                rGradientsValue[0](0, j) = 0.0;
            if (rDistances(0) < 0.0) rPartitionsSign[0] = -1.0;
            else rPartitionsSign[0] = 1.0;
            //KRATOS_WATCH("one element not in the intefase")
            return 1;
        }

        //else //we must create the enrichement, it can be in 2 or 3 parts. we'll start with 3 always.


        //const double epsilon = 1e-15; //1.00e-9;
        //compute the gradient of the distance and normalize it
        array_1d<double, 2 > grad_d;
        noalias(grad_d) = prod(trans(DN_DX), rDistances);
        /*
        double norm = norm_2(grad_d);
        if (norm > epsilon)
        	grad_d /= (norm);
        */
        array_1d<double, 3> exact_distance = rDistances;
        array_1d<double, 3> abs_distance = ZeroVector(3);


        //KRATOS_WATCH("one element IS in the intefase")
        if ((rDistances(0)*rDistances(1))<0.0) //edge 12 is cut
            cut_edges[0]=true;
        else
            cut_edges[0]=false;
        if ((rDistances(1)*rDistances(2))<0.0) //edge 23 is cut.
            cut_edges[1]=true;
        else
            cut_edges[1]=false;
        if ((rDistances(2)*rDistances(0))<0.0) //edge 23 is cut.
            cut_edges[2]=true;
        else
            cut_edges[2]=false;

        //'TRICK' TO AVOID HAVING THE INTERFASE TOO CLOSE TO THE NODES:
        //since we cannot collapse node because we have to contemplate the possibility of discontinuities, we will move a little the intefase so that it is not that close.
        const double unsigned_distance0=fabs(rDistances(0));
        const double unsigned_distance1=fabs(rDistances(1));
        const double unsigned_distance2=fabs(rDistances(2));
        //we begin by finding the largest distance:
        double longest_distance=fabs(unsigned_distance0);
        if (unsigned_distance1>longest_distance)
            longest_distance=unsigned_distance1;
        if (unsigned_distance2>longest_distance)
            longest_distance=unsigned_distance2;
        //Now we set a maximum relative distance
        const double tolerable_distance =longest_distance*0.001;	// (1/1,000,000 seems to have good results)
        //and now we check if a distance is too small:
        if (unsigned_distance0<tolerable_distance)
            rDistances[0]=tolerable_distance*(rDistances[0]/fabs(rDistances[0]));
        if (unsigned_distance1<tolerable_distance)
            rDistances[1]=tolerable_distance*(rDistances[1]/fabs(rDistances[1]));
        if (unsigned_distance2<tolerable_distance)
            rDistances[2]=tolerable_distance*(rDistances[2]/fabs(rDistances[2]));
        //END OF TRICK. REMEMBER TO OVERWRITE THE DISTANCE VARIABLE IN THE ELEMENT IN CASE THESE LINES HAVE MODIFIED THEM (distances)


        //for (int jj = 0; jj < 3; jj++)
        //	KRATOS_WATCH(rDistances(jj));
        for (unsigned int i=0; i<3; i++) //we go over the 3 edges:
        {
            int edge_begin_node=i;
            int edge_end_node=i+1;
            if (edge_end_node==3) edge_end_node=0; //it's a triangle, so node 3 is actually node 0

            if(cut_edges(i)==true)
            {
                aux_nodes_relative_locations(i)=fabs(rDistances(edge_end_node)/(rDistances(edge_end_node)-rDistances(edge_begin_node) ) ) ; //position in 'natural' coordinates of edge 12, 1 when it passes over node 1. (it is over the edge 01)
                aux_nodes_father_nodes(i,0)=edge_begin_node;
                aux_nodes_father_nodes(i,1)=edge_end_node;
            }
            else
            {
                if(fabs(rDistances(edge_end_node))>fabs(rDistances(edge_begin_node))) //if edge is not cut, we collapse the aux node into the node which has the highest absolute value to have "nicer" (less "slivery") subelements
                {
                    aux_nodes_relative_locations(i)=0.0;
                    aux_nodes_father_nodes(i,0)=edge_end_node;
                    aux_nodes_father_nodes(i,1)=edge_end_node;
                }
                else
                {
                    aux_nodes_relative_locations(i)=1.0;
                    aux_nodes_father_nodes(i,0)=edge_begin_node;
                    aux_nodes_father_nodes(i,1)=edge_begin_node;
                }
            }

            //and we save the coordinate of the new aux nodes:
            for (unsigned int j=0; j<2; j++)	//x,y coordinates
                aux_points(i,j)= rPoints(edge_begin_node,j) * aux_nodes_relative_locations(i) + rPoints(edge_end_node,j) * (1.0- aux_nodes_relative_locations(i));
        }


        array_1d<double,2> base_point;
        if (cut_edges(0)==true) // it means it is a cut edge, if it was 0.0 or 1.0 then it would be an uncut edge
        {
            base_point[0] = aux_points(0,0);
            base_point[1] = aux_points(0,1);
        }
        else //it means aux_point 0 is a clone of other point, so we go to the second edge.
        {
            base_point[0] = aux_points(1,0);
            base_point[1] = aux_points(1,1);
        }

        for (int i_node = 0; i_node < 3; i_node++)
        {
            double d =    (rPoints(i_node,0) - base_point[0]) * grad_d[0] +
                          (rPoints(i_node,1) - base_point[1]) * grad_d[1] ;
            abs_distance[i_node] = fabs(d);
        }

        //assign correct sign to exact distance
        for (int i = 0; i < 3; i++)
        {
            if (rDistances[i] < 0.0)
            {
                exact_distance[i] = -abs_distance[i];
                --most_common_sign;
            }
            else
            {
                exact_distance[i] = abs_distance[i];
                ++most_common_sign;
            }
        }





        //compute exact distance gradients
        array_1d<double, 2 > exact_distance_gradient;
        noalias(exact_distance_gradient) = prod(trans(DN_DX), exact_distance);

        array_1d<double, 2 > abs_distance_gradient;
        noalias(abs_distance_gradient) = prod(trans(DN_DX), abs_distance);


        double max_aux_dist_on_cut = -1;
        for (int edge = 0; edge < 3; edge++)
        {
            const int i = edge;
            int j = edge+1;
            if (j==3) j=0;
            if (rDistances[i] * rDistances[j] < 0.0)
            {
                const double tmp = fabs(rDistances[i]) / (fabs(rDistances[i]) + fabs(rDistances[j]));
                //compute the position of the edge node
                double abs_dist_on_cut = abs_distance[i] * tmp + abs_distance[j] * (1.00 - tmp);
                if(abs_dist_on_cut > max_aux_dist_on_cut) max_aux_dist_on_cut = abs_dist_on_cut;
            }
        }


        //we reset all data:
        rGradientsValue[0]=ZeroMatrix(2,2);
        rGradientsValue[1]=ZeroMatrix(2,2);
        rGradientsValue[2]=ZeroMatrix(2,2);
        NEnriched=ZeroMatrix(3,2);
        rGPShapeFunctionValues=ZeroMatrix(3,3);


        //now we must check the 4 created partitions of the domain.
        //one has been collapsed, so we discard it and therefore save only one.
        unsigned int partition_number=0;		//
        //the 3 first partitions are  created using 2 auxiliary nodes and a normal node. at least one of these will be discarded due to zero area
        //the last one is composed by the 3 auxiliary nodes. it 'looks' wrong, but since at least one has been collapsed, it actually has a normal node.
        bool found_empty_partition=false;
        for (unsigned int i=0; i<4; i++) //i partition
        {
            unsigned int j_aux = i + 2;
            if (j_aux>2) j_aux -= 3;
            BoundedMatrix<int,3,2> partition_father_nodes;
            array_1d<double,3> N;
            if (i<3)
            {
                partition_father_nodes(0,0)=i;
                partition_father_nodes(0,1)=i;
                partition_father_nodes(1,0)=aux_nodes_father_nodes(i,0); //we are using i aux node
                partition_father_nodes(1,1)=aux_nodes_father_nodes(i,1); //we are using i aux node
                partition_father_nodes(2,0)=aux_nodes_father_nodes(j_aux,0); //we are using j_aux node
                partition_father_nodes(2,1)=aux_nodes_father_nodes(j_aux,1); //we are using j_aux node

                coord_subdomain(0,0)=rPoints(i,0);
                coord_subdomain(0,1)=rPoints(i,1);
                coord_subdomain(1,0)=aux_points(i,0);
                coord_subdomain(1,1)=aux_points(i,1);
                coord_subdomain(2,0)=aux_points(j_aux,0);
                coord_subdomain(2,1)=aux_points(j_aux,1);
            }
            else
            {
                //the last partition, made by the 3 aux nodes.
                partition_father_nodes=aux_nodes_father_nodes;
                coord_subdomain=aux_points;
                //found_last_partition=true;
            }
            //calculate data of this partition
            double temp_area;
            CalculateGeometryData(coord_subdomain, DN_DX_subdomain, temp_area);
            if (temp_area > 1.0e-20) //ok, it does not have zero area
            {
                rVolumes(partition_number)=temp_area;
                //we look for the gauss point of the partition:
                double x_GP_partition =  one_third * ( coord_subdomain(0,0) + coord_subdomain(1,0) + coord_subdomain(2,0) );
                double y_GP_partition =  one_third * ( coord_subdomain(0,1) + coord_subdomain(1,1) + coord_subdomain(2,1) );
                double z_GP_partition  = 0.0;
                //we reset the coord_subdomain matrix so that we have the whole element again:
                coord_subdomain = rPoints;
                //and we calculate its shape function values
                CalculatePosition ( coord_subdomain , x_GP_partition ,y_GP_partition ,z_GP_partition , N);
                //we check the partition sign.
                const double partition_sign = (N(0)*rDistances(0) + N(1)*rDistances(1) + N(2)*rDistances(2))/fabs(N(0)*rDistances(0) + N(1)*rDistances(1) + N(2)*rDistances(2));
                //rPartitionsSign(partition_number)=partition_sign;

                rGPShapeFunctionValues(partition_number,0)=N(0);
                rGPShapeFunctionValues(partition_number,1)=N(1);
                rGPShapeFunctionValues(partition_number,2)=N(2);

                //compute enriched shape function values
                double dist = 0.0;
                double abs_dist = 0.0;
                for (int j = 0; j < 3; j++)
                {
                    dist += rGPShapeFunctionValues(partition_number, j) * exact_distance[j];
                    abs_dist += rGPShapeFunctionValues(partition_number, j) * abs_distance[j];
                }

                if (partition_sign < 0.0)
                    rPartitionsSign[partition_number] = -1.0;
                else
                    rPartitionsSign[partition_number] = 1.0;

                //enrichment for the gradient (continuous pressure, discontinuous gradient. Coppola-Codina enrichment
                NEnriched(partition_number, 0) = 0.5/max_aux_dist_on_cut * (abs_dist - rPartitionsSign[partition_number] * dist);
                for (int j = 0; j < 2; j++)
                {
                    rGradientsValue[partition_number](0, j) = (0.5/max_aux_dist_on_cut) * (abs_distance_gradient[j] - rPartitionsSign[partition_number] * exact_distance_gradient[j]);
                }

                //enrichment to add a constant discontinuity along the interfase.
                if (rPartitionsSign[partition_number]*most_common_sign > 0.0) //on this side we will copy the standard shape functions. (maybe we change the sign, but keeping the sign
                {
                    NEnriched(partition_number, 1) = -1.0*rPartitionsSign[partition_number]*NEnriched(partition_number, 0) ;
                    for (int j = 0; j < 2; j++)
                    {
                        rGradientsValue[partition_number](1, j) = -1.0*rPartitionsSign[partition_number]*rGradientsValue[partition_number](0, j);
                    }
                }
                else //we have to construct the shape function to guarantee a constant jump to 2:
                {
                    //notice that the partition to be changed must one the subdomains created with a real node. so the i index is also the real node. (used 4 lines below to recover the distance.

                    NEnriched(partition_number, 1) = rPartitionsSign[partition_number]*( 1.0*NEnriched(partition_number, 0) - 0.6666666666666666666666666666 ) ;
                    for (int j = 0; j < 2; j++)
                    {
                        rGradientsValue[partition_number](1, j) =  (rPartitionsSign[partition_number]*1.0*rGradientsValue[partition_number](0, j) - exact_distance_gradient[j]*1.0/(abs_distance[i]));
                    }
                }


                partition_number++;

            }
            else
                found_empty_partition=true;

        }
        if (found_empty_partition==false)
            KRATOS_WATCH("WROOOONGGGGGGGGGGG");

        //KRATOS_WATCH(NEnriched)
        return 3;
        KRATOS_CATCH("");

    }

    //2D but using 2 enrichments for the gradient + 1 for the jump + a single one for gradient discontinuity
    static int CalculateEnrichedShapeFuncionsExtended(BoundedMatrix<double,(2+1), 2 >& rPoints, BoundedMatrix<double, (2+1), 2 >& DN_DX,
            array_1d<double,(2+1)>& rDistances, array_1d<double,(3*(2-1))>& rVolumes, BoundedMatrix<double, 3*(2-1), (2+1) >& rGPShapeFunctionValues,
            array_1d<double,(3*(2-1))>& rPartitionsSign, std::vector<Matrix>& rGradientsValue, BoundedMatrix<double,3*(2-1), (4)>& NEnriched)
    {
        KRATOS_TRY

        const double one_third=1.0/3.0;
        BoundedMatrix<double,3,2> aux_points; //for auxiliary nodes 4(between 1 and 2) ,5(between 2 and 3) ,6 (between 3 and 1)
        BoundedMatrix<double, 3, 2 > coord_subdomain; //used to pass arguments when we must calculate areas, shape functions, etc
        BoundedMatrix<double,3,2> DN_DX_subdomain; //used to retrieve derivatives

        double most_common_sign=0; //the side of the cut in which two nodes are found (same sign) will be the ones that remains unchanged when builing the discontinuity
        double Area;//area of the complete element
        rGPShapeFunctionValues(0,0)=one_third;
        rGPShapeFunctionValues(0,1)=one_third;
        rGPShapeFunctionValues(0,2)=one_third; //default, when no interfase has been found
        Area = CalculateVolume2D( rPoints );
        array_1d<bool,3> cut_edges;
        array_1d<double,3> aux_nodes_relative_locations;
        BoundedMatrix<int,3,2> aux_nodes_father_nodes;

        //to begin with we must check whether our element is cut or not by the interfase.
        if( (rDistances(0)*rDistances(1))>0.0 && (rDistances(0)*rDistances(2))>0.0 ) //it means that this element IS NOT cut by the interfase. we must return data of a normal, non-enriched element
        {
            rVolumes(0)=Area;
            rGPShapeFunctionValues(0,0)=one_third;
            rGPShapeFunctionValues(0,1)=one_third;
            rGPShapeFunctionValues(0,2)=one_third;
            NEnriched(0,0) = 0.0;
            //type_of_cut=1;
            for (int j = 0; j < 2; j++)
                rGradientsValue[0](0, j) = 0.0;
            if (rDistances(0) < 0.0) rPartitionsSign[0] = -1.0;
            else rPartitionsSign[0] = 1.0;
            //KRATOS_WATCH("one element not in the intefase")
            return 1;
        }

        //else //we must create the enrichement, it can be in 2 or 3 parts. we'll start with 3 always.


        //const double epsilon = 1e-15; //1.00e-9;
        //compute the gradient of the distance and normalize it
        array_1d<double, 2 > grad_d;
        noalias(grad_d) = prod(trans(DN_DX), rDistances);
        /*
        double norm = norm_2(grad_d);
        if (norm > epsilon)
        	grad_d /= (norm);
        */
        array_1d<double, 3> exact_distance = rDistances;
        array_1d<double, 3> abs_distance = ZeroVector(3);


        //KRATOS_WATCH("one element IS in the intefase")
        if ((rDistances(0)*rDistances(1))<0.0) //edge 12 is cut
            cut_edges[0]=true;
        else
            cut_edges[0]=false;
        if ((rDistances(1)*rDistances(2))<0.0) //edge 23 is cut.
            cut_edges[1]=true;
        else
            cut_edges[1]=false;
        if ((rDistances(2)*rDistances(0))<0.0) //edge 23 is cut.
            cut_edges[2]=true;
        else
            cut_edges[2]=false;

        //'TRICK' TO AVOID HAVING THE INTERFASE TOO CLOSE TO THE NODES:
        //since we cannot collapse node because we have to contemplate the possibility of discontinuities, we will move a little the intefase so that it is not that close.
        const double unsigned_distance0=fabs(rDistances(0));
        const double unsigned_distance1=fabs(rDistances(1));
        const double unsigned_distance2=fabs(rDistances(2));
        //we begin by finding the largest distance:
        double longest_distance=fabs(unsigned_distance0);
        if (unsigned_distance1>longest_distance)
            longest_distance=unsigned_distance1;
        if (unsigned_distance2>longest_distance)
            longest_distance=unsigned_distance2;
        //Now we set a maximum relative distance
        const double tolerable_distance =longest_distance*0.001;	// (1/1,000,000 seems to have good results)
        //and now we check if a distance is too small:
        if (unsigned_distance0<tolerable_distance)
            rDistances[0]=tolerable_distance*(rDistances[0]/fabs(rDistances[0]));
        if (unsigned_distance1<tolerable_distance)
            rDistances[1]=tolerable_distance*(rDistances[1]/fabs(rDistances[1]));
        if (unsigned_distance2<tolerable_distance)
            rDistances[2]=tolerable_distance*(rDistances[2]/fabs(rDistances[2]));
        //END OF TRICK. REMEMBER TO OVERWRITE THE DISTANCE VARIABLE IN THE ELEMENT IN CASE THESE LINES HAVE MODIFIED THEM (distances)


        //for (int jj = 0; jj < 3; jj++)
        //	KRATOS_WATCH(rDistances(jj));


        //We have 3 edges, meaning we created 3 aux nodes. But one of them actually matches the position of a real node (the one that is not on an interface edge is displaced to one of the ends (a node)
        //the new shape functions are built by setting the values in all (real and aux) nodes to zero except in one of the interphase nodes. in the array aux_node_shape_function_index we assign this value
        array_1d<int, 3 > aux_node_enrichment_shape_function_index; //when not used, it must be -1;

        int shape_function_id=0;

        for (unsigned int i=0; i<3; i++) //we go over the 3 edges:
        {
            int edge_begin_node=i;
            int edge_end_node=i+1;
            if (edge_end_node==3) edge_end_node=0; //it's a triangle, so node 3 is actually node 0

            if(cut_edges(i)==true)
            {
                aux_nodes_relative_locations(i)=fabs(rDistances(edge_end_node)/(rDistances(edge_end_node)-rDistances(edge_begin_node) ) ) ; //position in 'natural' coordinates of edge 12, 1 when it passes over node 1. (it is over the edge 01)
                aux_nodes_father_nodes(i,0)=edge_begin_node;
                aux_nodes_father_nodes(i,1)=edge_end_node;

                aux_node_enrichment_shape_function_index(i)=shape_function_id;
                shape_function_id++;
            }
            else
            {
                if(fabs(rDistances(edge_end_node))>fabs(rDistances(edge_begin_node))) //if edge is not cut, we collapse the aux node into the node which has the highest absolute value to have "nicer" (less "slivery") subelements
                {
                    aux_nodes_relative_locations(i)=0.0;
                    aux_nodes_father_nodes(i,0)=edge_end_node;
                    aux_nodes_father_nodes(i,1)=edge_end_node;
                }
                else
                {
                    aux_nodes_relative_locations(i)=1.0;
                    aux_nodes_father_nodes(i,0)=edge_begin_node;
                    aux_nodes_father_nodes(i,1)=edge_begin_node;
                }

                aux_node_enrichment_shape_function_index(i)=-1;
            }

            //and we save the coordinate of the new aux nodes:
            for (unsigned int j=0; j<2; j++)	//x,y coordinates
                aux_points(i,j)= rPoints(edge_begin_node,j) * aux_nodes_relative_locations(i) + rPoints(edge_end_node,j) * (1.0- aux_nodes_relative_locations(i));
        }


        array_1d<double,2> base_point;
        if (cut_edges(0)==true) // it means it is a cut edge, if it was 0.0 or 1.0 then it would be an uncut edge
        {
            base_point[0] = aux_points(0,0);
            base_point[1] = aux_points(0,1);
        }
        else //it means aux_point 0 is a clone of other point, so we go to the second edge.
        {
            base_point[0] = aux_points(1,0);
            base_point[1] = aux_points(1,1);
        }

        for (int i_node = 0; i_node < 3; i_node++)
        {
            double d =    (rPoints(i_node,0) - base_point[0]) * grad_d[0] +
                          (rPoints(i_node,1) - base_point[1]) * grad_d[1] ;
            abs_distance[i_node] = fabs(d);
        }

        //assign correct sign to exact distance
        for (int i = 0; i < 3; i++)
        {
            if (rDistances[i] < 0.0)
            {
                exact_distance[i] = -abs_distance[i];
                --most_common_sign;
            }
            else
            {
                exact_distance[i] = abs_distance[i];
                ++most_common_sign;
            }
        }

        //compute exact distance gradients
        array_1d<double, 2 > exact_distance_gradient;
        noalias(exact_distance_gradient) = prod(trans(DN_DX), exact_distance);

        array_1d<double, 2 > abs_distance_gradient;
        noalias(abs_distance_gradient) = prod(trans(DN_DX), abs_distance);


        double max_aux_dist_on_cut = -1;
        for (int edge = 0; edge < 3; edge++)
        {
            const int i = edge;
            int j = edge+1;
            if (j==3) j=0;
            if (rDistances[i] * rDistances[j] < 0.0)
            {
                const double tmp = fabs(rDistances[i]) / (fabs(rDistances[i]) + fabs(rDistances[j]));
                //compute the position of the edge node
                double abs_dist_on_cut = abs_distance[i] * tmp + abs_distance[j] * (1.00 - tmp);
                if(abs_dist_on_cut > max_aux_dist_on_cut) max_aux_dist_on_cut = abs_dist_on_cut;
            }
        }


        //we reset all data:
        rGradientsValue[0]=ZeroMatrix(4,2);
        rGradientsValue[1]=ZeroMatrix(4,2);
        rGradientsValue[2]=ZeroMatrix(4,2);
        NEnriched=ZeroMatrix(3,4);
        rGPShapeFunctionValues=ZeroMatrix(3,3);


        //now we must check the 4 created partitions of the domain.
        //one has been collapsed, so we discard it and therefore save only one.
        unsigned int partition_number=0;		//
        //the 3 first partitions are  created using 2 auxiliary nodes and a normal node. at least one of these will be discarded due to zero area
        //the last one is composed by the 3 auxiliary nodes. it 'looks' wrong, but since at least one has been collapsed, it actually has a normal node.
        bool found_empty_partition=false;

        //the enrichment is directly the shape functions created by using the partition.
        //we have to save for the enrichment 0 and 1 which is the node that will be active, that is, whose shape function value is zero (for some partitions all 3 shape functions are inactive)



        for (unsigned int i=0; i<4; i++) //i partition
        {

            array_1d<int, 2 > active_node_in_enrichment_shape_function;
            active_node_in_enrichment_shape_function(0)=-1;
            active_node_in_enrichment_shape_function(1)=-1; //initialized as if all are inactive -> gradient=0;
            //the same but for the replacement shape functions
            array_1d<int, 3 > active_node_in_replacement_shape_function;
            active_node_in_replacement_shape_function(0)=-1;
            active_node_in_replacement_shape_function(1)=-1;
            active_node_in_replacement_shape_function(2)=-1; //initialized as if all are inactive -> gradient=0;

            unsigned int j_aux = i + 2;
            if (j_aux>2) j_aux -= 3;
            BoundedMatrix<int,3,2> partition_father_nodes;
            array_1d<double,3> N;
            if (i<3)
            {
                partition_father_nodes(0,0)=i;
                partition_father_nodes(0,1)=i;
                partition_father_nodes(1,0)=aux_nodes_father_nodes(i,0); //we are using i aux node
                partition_father_nodes(1,1)=aux_nodes_father_nodes(i,1); //we are using i aux node
                partition_father_nodes(2,0)=aux_nodes_father_nodes(j_aux,0); //we are using j_aux node
                partition_father_nodes(2,1)=aux_nodes_father_nodes(j_aux,1); //we are using j_aux node

                coord_subdomain(0,0)=rPoints(i,0);
                coord_subdomain(0,1)=rPoints(i,1);
                coord_subdomain(1,0)=aux_points(i,0);
                coord_subdomain(1,1)=aux_points(i,1);
                coord_subdomain(2,0)=aux_points(j_aux,0);
                coord_subdomain(2,1)=aux_points(j_aux,1);

                //notice that local nodes 2 and 3 and the possible candidates, with indexes i and j_aux:
                if (aux_node_enrichment_shape_function_index(i)> -1) //that is, local node 2 it is a useful node:
                    active_node_in_enrichment_shape_function(  aux_node_enrichment_shape_function_index(i)  )=1;	//saving that local node 2 will be active for either enrichment 1 or 2.
                // else we do nothing, we are not saving this node and the -1 stays

                //now the same for the local node 3 (j_aux)
                if (aux_node_enrichment_shape_function_index(j_aux)> -1) //that is, local node 3 it is a useful node:
                    active_node_in_enrichment_shape_function( aux_node_enrichment_shape_function_index(j_aux) )=2;	//saving that local node 3 will be active for either enrichment 1 or 2.
                // else we do nothing, we are not saving this node and the -1 stays

                active_node_in_replacement_shape_function(i)=0; //standard shape function i will be replaced by the one of local subelement node 1
                //now local nodes 2 and 3
                if (aux_nodes_father_nodes(i,0)==aux_nodes_father_nodes(i,1))
                    active_node_in_replacement_shape_function(aux_nodes_father_nodes(i,0))=1;
                if (aux_nodes_father_nodes(j_aux,0)==aux_nodes_father_nodes(j_aux,1))
                    active_node_in_replacement_shape_function(aux_nodes_father_nodes(j_aux,0))=2;

            }
            else
            {
                //the last partition, made by the 3 aux nodes.
                partition_father_nodes=aux_nodes_father_nodes;
                coord_subdomain=aux_points;

                //we have to check the 3 of them:
                for (int j = 0; j < 3; j++)
                {
                    if (aux_node_enrichment_shape_function_index(j)> -1) //that is, local node j it is a useful node:
                        active_node_in_enrichment_shape_function(  aux_node_enrichment_shape_function_index(j)  ) = j;

                    //to replace the standard shape functions:
                    if (aux_nodes_father_nodes(j,0)==aux_nodes_father_nodes(j,1))
                        active_node_in_replacement_shape_function(aux_nodes_father_nodes(j,0))=j;
                }

                //found_last_partition=true;
            }
            //calculate data of this partition
            double temp_area;
            CalculateGeometryData(coord_subdomain, DN_DX_subdomain, temp_area);
            if (temp_area > 1.0e-20) //ok, it does not have zero area
            {
                rVolumes(partition_number)=temp_area;
                //we look for the gauss point of the partition:
                double x_GP_partition =  one_third * ( coord_subdomain(0,0) + coord_subdomain(1,0) + coord_subdomain(2,0) );
                double y_GP_partition =  one_third * ( coord_subdomain(0,1) + coord_subdomain(1,1) + coord_subdomain(2,1) );
                double z_GP_partition  = 0.0;
                //we reset the coord_subdomain matrix so that we have the whole element again:
                coord_subdomain = rPoints;
                //and we calculate its shape function values
                CalculatePosition ( coord_subdomain , x_GP_partition ,y_GP_partition ,z_GP_partition , N);
                //we check the partition sign.
                const double partition_sign = (N(0)*rDistances(0) + N(1)*rDistances(1) + N(2)*rDistances(2))/fabs(N(0)*rDistances(0) + N(1)*rDistances(1) + N(2)*rDistances(2));
                //rPartitionsSign(partition_number)=partition_sign;

                rGPShapeFunctionValues(partition_number,0)=N(0);
                rGPShapeFunctionValues(partition_number,1)=N(1);
                rGPShapeFunctionValues(partition_number,2)=N(2);

                //compute enriched shape function values
                double dist = 0.0;
                double abs_dist = 0.0;
                for (int j = 0; j < 3; j++)
                {
                    dist += rGPShapeFunctionValues(partition_number, j) * exact_distance[j];
                    abs_dist += rGPShapeFunctionValues(partition_number, j) * abs_distance[j];
                }

                if (partition_sign < 0.0)
                    rPartitionsSign[partition_number] = -1.0;
                else
                    rPartitionsSign[partition_number] = 1.0;

                //We use the sublement shape functions and derivatives:
                //we loop the 2 enrichment shape functions:
                for (int index_shape_function = 0; index_shape_function < 2; index_shape_function++) //enrichment shape function
                {
                    if (active_node_in_enrichment_shape_function(index_shape_function) > -1) //recall some of them are inactive:
                    {
                        NEnriched(partition_number, index_shape_function ) = one_third ; //only one gauss point. if more were to be used, it would still be simple (1/2,1/2,0);(0,1/2,1/2);(1/2,0,1/2);
                        //NEnriched(partition_number, index_shape_function+2 ) = one_third*rPartitionsSign[partition_number] ; //only one gauss point. if more were to be used, it would still be simple (1/2,1/2,0);(0,1/2,1/2);(1/2,0,1/2);
                        for (int j = 0; j < 2; j++) //x,y,(z)
                        {
                            rGradientsValue[partition_number](index_shape_function, j) = DN_DX_subdomain(active_node_in_enrichment_shape_function(index_shape_function),j);
                            //rGradientsValue[partition_number](index_shape_function+2, j) = DN_DX_subdomain(active_node_in_enrichment_shape_function(index_shape_function),j)*rPartitionsSign[partition_number];
                        }

                    }
                    else
                    {
                        NEnriched(partition_number, index_shape_function ) = 0.0;
                        //NEnriched(partition_number, index_shape_function +2) = 0.0;
                        for (int j = 0; j < 2; j++) //x,y,(z)
                        {
                            rGradientsValue[partition_number](index_shape_function, j) = 0.0;
                            //rGradientsValue[partition_number](index_shape_function+2, j) = 0.0;
                        }
                    }

                }

                //We use the sublement shape functions and derivatives:
                //we loop the 3 replacement shape functions:
                unsigned int number_of_real_nodes=0; //(each partition can have either 1 or 2);
                for (int index_shape_function = 0; index_shape_function < 3; index_shape_function++) //enrichment shape function
                {
                    if (active_node_in_replacement_shape_function(index_shape_function) > -1) //recall some of them are inactive:
                        number_of_real_nodes++;
                }

                if (number_of_real_nodes==2) //the jump enrichment has to be constructed from the auxiliar node.
                {
                    for (int index_shape_function = 0; index_shape_function < 2; index_shape_function++) //enrichment shape function
                    {
                        if (active_node_in_enrichment_shape_function(index_shape_function) > -1)
                        {
                            NEnriched(partition_number, 2 ) = one_third*rPartitionsSign[partition_number] ; //at interface, values are 1 for positive side and -1 for the negative . Constant jump of 2.
                            NEnriched(partition_number, 3 ) = one_third;									//at interface, both sides are constant 1. 		No jump in the interface, discontinuous gradient.
                            for (int j = 0; j < 2; j++) //x,y,(z)
                            {
                                rGradientsValue[partition_number](2, j) = DN_DX_subdomain(active_node_in_enrichment_shape_function(index_shape_function),j)*rPartitionsSign[partition_number]; //at interface, values are 1 for positive side and -1 for the negative
                                rGradientsValue[partition_number](3, j) = DN_DX_subdomain(active_node_in_enrichment_shape_function(index_shape_function),j); //at interface, both sides are constant 1.
                            }
                            break;
                        }
                    }
                }
                else //there is one real node, so we construct the jump enrichment from this function but changing the sign
                {
                    for (int index_shape_function = 0; index_shape_function < 3; index_shape_function++) //enrichment shape function
                    {
                        if (active_node_in_replacement_shape_function(index_shape_function) > -1)
                        {
                            NEnriched(partition_number, 2 ) = one_third*rPartitionsSign[partition_number]*2.0 ;
                            NEnriched(partition_number, 3 ) = one_third*2.0 ;
                            for (int j = 0; j < 2; j++) //x,y,(z)
                            {
                                rGradientsValue[partition_number](2, j) =  - DN_DX_subdomain(active_node_in_replacement_shape_function(index_shape_function),j)*rPartitionsSign[partition_number]; //at interface, values are 1 for positive side and -1 for the negative
                                rGradientsValue[partition_number](3, j) =  - DN_DX_subdomain(active_node_in_replacement_shape_function(index_shape_function),j); //at interface, both sides are constant 1.
                            }
                            break;
                        }
                    }
                }
                partition_number++;

            }
            else
                found_empty_partition=true;
        }
        if (found_empty_partition==false)
            KRATOS_WATCH("WROOOONGGGGGGGGGGG");

        //KRATOS_WATCH(NEnriched)
        return 3;
        KRATOS_CATCH("");

    }


    //2D with information on the interfase
    static int CalculateEnrichedShapeFuncions(BoundedMatrix<double,(2+1), 2 >& rPoints, BoundedMatrix<double, (2+1), 2 >& DN_DX,
            array_1d<double,(2+1)>& rDistances, array_1d<double,(3*(2-1))>& rVolumes, BoundedMatrix<double, 3*(2-1), (2+1) >& rGPShapeFunctionValues,
            array_1d<double,(3*(2-1))>& rPartitionsSign, std::vector<Matrix>& rGradientsValue, BoundedMatrix<double,3*(2-1), (2)>& NEnriched,
            //information about the interfase: Standard shape functions values, Enrichment Functions (first(gradient), second (positive_side,negative side) ) and Area
            array_1d<double,(3)>&  rGPShapeFunctionValues_in_interfase, array_1d<double,(3)>&  NEnriched_in_interfase, double & InterfaseArea)
    {
        KRATOS_TRY

        const double one_third=1.0/3.0;
        BoundedMatrix<double,3,2> aux_points; //for auxiliary nodes 4(between 1 and 2) ,5(between 2 and 3) ,6 (between 3 and 1)
        BoundedMatrix<double, 3, 2 > coord_subdomain; //used to pass arguments when we must calculate areas, shape functions, etc
        BoundedMatrix<double,3,2> DN_DX_subdomain; //used to retrieve derivatives

        double most_common_sign=0; //the side of the cut in which two nodes are found (same sign) will be the ones that remains unchanged when builing the discontinuity
        double Area;//area of the complete element
        rGPShapeFunctionValues(0,0)=one_third;
        rGPShapeFunctionValues(0,1)=one_third;
        rGPShapeFunctionValues(0,2)=one_third; //default, when no interfase has been found
        Area = CalculateVolume2D( rPoints );
        array_1d<bool,3> cut_edges;
        array_1d<double,3> aux_nodes_relative_locations;
        BoundedMatrix<int,3,2> aux_nodes_father_nodes;

        //to begin with we must check whether our element is cut or not by the interfase.
        if( (rDistances(0)*rDistances(1))>0.0 && (rDistances(0)*rDistances(2))>0.0 ) //it means that this element IS NOT cut by the interfase. we must return data of a normal, non-enriched element
        {
            rVolumes(0)=Area;
            rGPShapeFunctionValues(0,0)=one_third;
            rGPShapeFunctionValues(0,1)=one_third;
            rGPShapeFunctionValues(0,2)=one_third;
            NEnriched(0,0) = 0.0;
            //type_of_cut=1;
            for (int j = 0; j < 3; j++)
                rGradientsValue[0](0, j) = 0.0;
            if (rDistances(0) < 0.0) rPartitionsSign[0] = -1.0;
            else rPartitionsSign[0] = 1.0;
            //KRATOS_WATCH("one element not in the intefase")
            return 1;
        }

        //else //we must create the enrichement, it can be in 2 or 3 parts. we'll start with 3 always.


        //const double epsilon = 1e-15; //1.00e-9;
        //compute the gradient of the distance and normalize it
        array_1d<double, 2 > grad_d;
        noalias(grad_d) = prod(trans(DN_DX), rDistances);
        /*
        double norm = norm_2(grad_d);
        if (norm > epsilon)
        	grad_d /= (norm);
        */
        array_1d<double, 3> exact_distance = rDistances;
        array_1d<double, 3> abs_distance = ZeroVector(3);


        //KRATOS_WATCH("one element IS in the intefase")
        if ((rDistances(0)*rDistances(1))<0.0) //edge 12 is cut
            cut_edges[0]=true;
        else
            cut_edges[0]=false;
        if ((rDistances(1)*rDistances(2))<0.0) //edge 23 is cut.
            cut_edges[1]=true;
        else
            cut_edges[1]=false;
        if ((rDistances(2)*rDistances(0))<0.0) //edge 23 is cut.
            cut_edges[2]=true;
        else
            cut_edges[2]=false;

        //'TRICK' TO AVOID HAVING THE INTERFASE TOO CLOSE TO THE NODES:
        //since we cannot collapse node because we have to contemplate the possibility of discontinuities, we will move a little the intefase so that it is not that close.
        const double unsigned_distance0=fabs(rDistances(0));
        const double unsigned_distance1=fabs(rDistances(1));
        const double unsigned_distance2=fabs(rDistances(2));
        //we begin by finding the largest distance:
        double longest_distance=fabs(unsigned_distance0);
        if (unsigned_distance1>longest_distance)
            longest_distance=unsigned_distance1;
        if (unsigned_distance2>longest_distance)
            longest_distance=unsigned_distance2;
        //Now we set a maximum relative distance
        const double tolerable_distance =longest_distance*0.001;	// (1/1,000,000 seems to have good results)
        //and now we check if a distance is too small:
        if (unsigned_distance0<tolerable_distance)
            rDistances[0]=tolerable_distance*(rDistances[0]/fabs(rDistances[0]));
        if (unsigned_distance1<tolerable_distance)
            rDistances[1]=tolerable_distance*(rDistances[1]/fabs(rDistances[1]));
        if (unsigned_distance2<tolerable_distance)
            rDistances[2]=tolerable_distance*(rDistances[2]/fabs(rDistances[2]));
        //END OF TRICK. REMEMBER TO OVERWRITE THE DISTANCE VARIABLE IN THE ELEMENT IN CASE THESE LINES HAVE MODIFIED THEM (distances)


        //for (int jj = 0; jj < 3; jj++)
        //	KRATOS_WATCH(rDistances(jj));
        for (unsigned int i=0; i<3; i++) //we go over the 3 edges:
        {
            int edge_begin_node=i;
            int edge_end_node=i+1;
            if (edge_end_node==3) edge_end_node=0; //it's a triangle, so node 3 is actually node 0

            if(cut_edges(i)==true)
            {
                aux_nodes_relative_locations(i)=fabs(rDistances(edge_end_node)/(rDistances(edge_end_node)-rDistances(edge_begin_node) ) ) ; //position in 'natural' coordinates of edge 12, 1 when it passes over node 1. (it is over the edge 01)
                aux_nodes_father_nodes(i,0)=edge_begin_node;
                aux_nodes_father_nodes(i,1)=edge_end_node;
            }
            else
            {
                if(fabs(rDistances(edge_end_node))>fabs(rDistances(edge_begin_node))) //if edge is not cut, we collapse the aux node into the node which has the highest absolute value to have "nicer" (less "slivery") subelements
                {
                    aux_nodes_relative_locations(i)=0.0;
                    aux_nodes_father_nodes(i,0)=edge_end_node;
                    aux_nodes_father_nodes(i,1)=edge_end_node;
                }
                else
                {
                    aux_nodes_relative_locations(i)=1.0;
                    aux_nodes_father_nodes(i,0)=edge_begin_node;
                    aux_nodes_father_nodes(i,1)=edge_begin_node;
                }
            }

            //and we save the coordinate of the new aux nodes:
            for (unsigned int j=0; j<2; j++)	//x,y coordinates
                aux_points(i,j)= rPoints(edge_begin_node,j) * aux_nodes_relative_locations(i) + rPoints(edge_end_node,j) * (1.0- aux_nodes_relative_locations(i));
        }


        array_1d<double,2> base_point;
        if (cut_edges(0)==true) // it means it is a cut edge, if it was 0.0 or 1.0 then it would be an uncut edge
        {
            base_point[0] = aux_points(0,0);
            base_point[1] = aux_points(0,1);

        }
        else //it means aux_point 0 is a clone of other point, so we go to the second edge.
        {
            base_point[0] = aux_points(1,0);
            base_point[1] = aux_points(1,1);
        }

        for (int i_node = 0; i_node < 3; i_node++)
        {
            double d =    (rPoints(i_node,0) - base_point[0]) * grad_d[0] +
                          (rPoints(i_node,1) - base_point[1]) * grad_d[1] ;
            abs_distance[i_node] = fabs(d);
        }

        //assign correct sign to exact distance
        for (int i = 0; i < 3; i++)
        {
            if (rDistances[i] < 0.0)
            {
                exact_distance[i] = -abs_distance[i];
                --most_common_sign;
            }
            else
            {
                exact_distance[i] = abs_distance[i];
                ++most_common_sign;
            }
        }





        //compute exact distance gradients
        array_1d<double, 2 > exact_distance_gradient;
        noalias(exact_distance_gradient) = prod(trans(DN_DX), exact_distance);

        array_1d<double, 2 > abs_distance_gradient;
        noalias(abs_distance_gradient) = prod(trans(DN_DX), abs_distance);


        double max_aux_dist_on_cut = -1;
        for (int edge = 0; edge < 3; edge++)
        {
            const int i = edge;
            int j = edge+1;
            if (j==3) j=0;
            if (rDistances[i] * rDistances[j] < 0.0)
            {
                const double tmp = fabs(rDistances[i]) / (fabs(rDistances[i]) + fabs(rDistances[j]));
                //compute the position of the edge node
                double abs_dist_on_cut = abs_distance[i] * tmp + abs_distance[j] * (1.00 - tmp);
                if(abs_dist_on_cut > max_aux_dist_on_cut) max_aux_dist_on_cut = abs_dist_on_cut;
            }
        }


        //we reset all data:
        rGradientsValue[0]=ZeroMatrix(2,2);
        rGradientsValue[1]=ZeroMatrix(2,2);
        rGradientsValue[2]=ZeroMatrix(2,2);
        NEnriched=ZeroMatrix(3,2);
        rGPShapeFunctionValues=ZeroMatrix(3,3);


        //now we must check the 4 created partitions of the domain.
        //one has been collapsed, so we discard it and therefore save only one.
        unsigned int partition_number=0;		//
        //the 3 first partitions are  created using 2 auxiliary nodes and a normal node. at least one of these will be discarded due to zero area
        //the last one is composed by the 3 auxiliary nodes. it 'looks' wrong, but since at least one has been collapsed, it actually has a normal node.
        bool found_empty_partition=false;
        for (unsigned int i=0; i<4; i++) //i partition
        {
            unsigned int j_aux = i + 2;
            if (j_aux>2) j_aux -= 3;
            BoundedMatrix<int,3,2> partition_father_nodes;
            array_1d<double,3> N;
            if (i<3)
            {
                partition_father_nodes(0,0)=i;
                partition_father_nodes(0,1)=i;
                partition_father_nodes(1,0)=aux_nodes_father_nodes(i,0); //we are using i aux node
                partition_father_nodes(1,1)=aux_nodes_father_nodes(i,1); //we are using i aux node
                partition_father_nodes(2,0)=aux_nodes_father_nodes(j_aux,0); //we are using j_aux node
                partition_father_nodes(2,1)=aux_nodes_father_nodes(j_aux,1); //we are using j_aux node

                coord_subdomain(0,0)=rPoints(i,0);
                coord_subdomain(0,1)=rPoints(i,1);
                coord_subdomain(1,0)=aux_points(i,0);
                coord_subdomain(1,1)=aux_points(i,1);
                coord_subdomain(2,0)=aux_points(j_aux,0);
                coord_subdomain(2,1)=aux_points(j_aux,1);
            }
            else
            {
                //the last partition, made by the 3 aux nodes.
                partition_father_nodes=aux_nodes_father_nodes;
                coord_subdomain=aux_points;
                //found_last_partition=true;
            }
            //calculate data of this partition
            double temp_area;
            CalculateGeometryData(coord_subdomain, DN_DX_subdomain, temp_area);
            if (temp_area > 1.0e-20) //ok, it does not have zero area
            {
                rVolumes(partition_number)=temp_area;
                //we look for the gauss point of the partition:
                double x_GP_partition =  one_third * ( coord_subdomain(0,0) + coord_subdomain(1,0) + coord_subdomain(2,0) );
                double y_GP_partition =  one_third * ( coord_subdomain(0,1) + coord_subdomain(1,1) + coord_subdomain(2,1) );
                double z_GP_partition  = 0.0;
                //we reset the coord_subdomain matrix so that we have the whole element again:
                coord_subdomain = rPoints;
                //and we calculate its shape function values
                CalculatePosition ( coord_subdomain , x_GP_partition ,y_GP_partition ,z_GP_partition , N);
                //we check the partition sign.
                const double partition_sign = (N(0)*rDistances(0) + N(1)*rDistances(1) + N(2)*rDistances(2))/fabs(N(0)*rDistances(0) + N(1)*rDistances(1) + N(2)*rDistances(2));
                //rPartitionsSign(partition_number)=partition_sign;

                rGPShapeFunctionValues(partition_number,0)=N(0);
                rGPShapeFunctionValues(partition_number,1)=N(1);
                rGPShapeFunctionValues(partition_number,2)=N(2);

                //compute enriched shape function values
                double dist = 0.0;
                double abs_dist = 0.0;
                for (int j = 0; j < 3; j++)
                {
                    dist += rGPShapeFunctionValues(partition_number, j) * exact_distance[j];
                    abs_dist += rGPShapeFunctionValues(partition_number, j) * abs_distance[j];
                }

                if (partition_sign < 0.0)
                    rPartitionsSign[partition_number] = -1.0;
                else
                    rPartitionsSign[partition_number] = 1.0;

                //enrichment for the gradient (continuous pressure, discontinuous gradient. Coppola-Codina enrichment
                NEnriched(partition_number, 0) = 0.5/max_aux_dist_on_cut * (abs_dist - rPartitionsSign[partition_number] * dist);
                for (int j = 0; j < 2; j++)
                {
                    rGradientsValue[partition_number](0, j) = (0.5/max_aux_dist_on_cut) * (abs_distance_gradient[j] - rPartitionsSign[partition_number] * exact_distance_gradient[j]);
                }

                //enrichment to add a constant discontinuity along the interfase.
                if (rPartitionsSign[partition_number]*most_common_sign > 0.0) //on this side we will copy the standard shape functions. (maybe we change the sign, but keeping the sign
                {
                    NEnriched(partition_number, 1) = -1.0*rPartitionsSign[partition_number]*NEnriched(partition_number, 0) ;
                    for (int j = 0; j < 2; j++)
                    {
                        rGradientsValue[partition_number](1, j) = -1.0*rPartitionsSign[partition_number]*rGradientsValue[partition_number](0, j);
                    }
                }
                else //we have to construct the shape function to guarantee a constant jump to 2:
                {
                    //notice that the partition to be changed must one the subdomains created with a real node. so the i index is also the real node. (used 4 lines below to recover the distance.

                    NEnriched(partition_number, 1) = rPartitionsSign[partition_number]*( 1.0*NEnriched(partition_number, 0) - 0.6666666666666666666666666666 ) ;
                    for (int j = 0; j < 2; j++)
                    {
                        rGradientsValue[partition_number](1, j) =  (rPartitionsSign[partition_number]*1.0*rGradientsValue[partition_number](0, j) - exact_distance_gradient[j]*1.0/(abs_distance[i]));
                    }

                }


                partition_number++;

            }
            else
                found_empty_partition=true;

        }
        if (found_empty_partition==false)
            KRATOS_WATCH("WROOOONGGGGGGGGGGG");

        // finally the interfase:
        if (cut_edges[0])
        {
            if (cut_edges[1])
            {
                InterfaseArea = sqrt(pow(aux_points(0,0)-aux_points(1,0),2)+pow(aux_points(0,1)-aux_points(1,1),2));
                base_point[0] = (aux_points(0,0)+aux_points(1,0))*0.5;
                base_point[1] = (aux_points(0,1)+aux_points(1,1))*0.5;
                CalculatePosition ( rPoints , base_point[0] ,base_point[1] , 0.0 , rGPShapeFunctionValues_in_interfase);
                double abs_dist_on_interfase = 0.0;
                for (int j = 0; j < 3; j++)
                    abs_dist_on_interfase += rGPShapeFunctionValues_in_interfase ( j) * abs_distance[j];
                NEnriched_in_interfase(0) = 0.5/max_aux_dist_on_cut * abs_dist_on_interfase;
            }
            else
            {
                InterfaseArea= sqrt(pow(aux_points(0,0)-aux_points(2,0),2)+pow(aux_points(0,1)-aux_points(2,1),2));
                base_point[0] = (aux_points(0,0)+aux_points(2,0))*0.5;
                base_point[1] = (aux_points(0,1)+aux_points(2,1))*0.5;
                CalculatePosition ( rPoints , base_point[0] ,base_point[1] , 0.0 , rGPShapeFunctionValues_in_interfase);
                double abs_dist_on_interfase = 0.0;
                for (int j = 0; j < 3; j++)
                    abs_dist_on_interfase += rGPShapeFunctionValues_in_interfase ( j) * abs_distance[j];
                NEnriched_in_interfase(0) = 0.5/max_aux_dist_on_cut * abs_dist_on_interfase;
            }
        }
        else
        {
            InterfaseArea= sqrt(pow(aux_points(2,0)-aux_points(1,0),2)+pow(aux_points(2,1)-aux_points(1,1),2));
            base_point[0] = (aux_points(2,0)+aux_points(1,0))*0.5;
            base_point[1] = (aux_points(2,1)+aux_points(1,1))*0.5;
            CalculatePosition ( rPoints , base_point[0] ,base_point[1] , 0.0 , rGPShapeFunctionValues_in_interfase);
            double abs_dist_on_interfase = 0.0;
            for (int j = 0; j < 3; j++)
                abs_dist_on_interfase += rGPShapeFunctionValues_in_interfase ( j) * abs_distance[j];
            NEnriched_in_interfase(0) = 0.5/max_aux_dist_on_cut * abs_dist_on_interfase;
        }


        //KRATOS_WATCH(NEnriched)
        return 3;
        KRATOS_CATCH("");

    }



    //2D
    static int CalculateEnrichedShapeFuncionsInLocalAxis(BoundedMatrix<double,(2+1), 2 >& rOriginalPoints, BoundedMatrix<double, (2+1), 2 >& DN_DX_original,
            array_1d<double,(2+1)>& rDistances, array_1d<double,(3*(2-1))>& rVolumes, BoundedMatrix<double, 3*(2-1), (2+1) >& rGPShapeFunctionValues,
            array_1d<double,(3*(2-1))>& rPartitionsSign, std::vector<Matrix>& rGradientsValue, BoundedMatrix<double,3*(2-1), (2)>& NEnriched, BoundedMatrix<double,(2), 2 >& rRotationMatrix, BoundedMatrix<double, (2+1), 2 > & DN_DX_in_local_axis)

    {
        KRATOS_TRY

        const double one_third=1.0/3.0;
        BoundedMatrix<double,3,2> aux_points; //for auxiliary nodes 4(between 1 and 2) ,5(between 2 and 3) ,6 (between 3 and 1)
        BoundedMatrix<double, 3, 2 > coord_subdomain; //used to pass arguments when we must calculate areas, shape functions, etc
        BoundedMatrix<double,3,2> DN_DX_subdomain; //used to retrieve derivatives

        BoundedMatrix<double,(2+1), 2 > rRotatedPoints;

        double most_common_sign=0; //the side of the cut in which two nodes are found (same sign) will be the ones that remains unchanged when builing the discontinuity
        double Area;//area of the complete element
        rGPShapeFunctionValues(0,0)=one_third;
        rGPShapeFunctionValues(0,1)=one_third;
        rGPShapeFunctionValues(0,2)=one_third; //default, when no interfase has been found
        Area = CalculateVolume2D( rOriginalPoints );
        array_1d<bool,3> cut_edges;
        array_1d<double,3> aux_nodes_relative_locations;
        BoundedMatrix<int,3,2> aux_nodes_father_nodes;

        //to begin with we must check whether our element is cut or not by the interfase.
        if( (rDistances(0)*rDistances(1))>0.0 && (rDistances(0)*rDistances(2))>0.0 ) //it means that this element IS NOT cut by the interfase. we must return data of a normal, non-enriched element
        {
            rVolumes(0)=Area;
            rGPShapeFunctionValues(0,0)=one_third;
            rGPShapeFunctionValues(0,1)=one_third;
            rGPShapeFunctionValues(0,2)=one_third;
            NEnriched(0,0) = 0.0;
            //type_of_cut=1;
            for (int j = 0; j < 3; j++)
                rGradientsValue[0](0, j) = 0.0;
            if (rDistances(0) < 0.0) rPartitionsSign[0] = -1.0;
            else rPartitionsSign[0] = 1.0;
            //KRATOS_WATCH("one element not in the intefase")
            return 1;
        }

        //else //we must create the enrichement, it can be in 2 or 3 parts. we'll start with 3 always.


        //const double epsilon = 1e-15; //1.00e-9;
        //compute the gradient of the distance and normalize it

        //double norm = norm_2(grad_d);
        //if (norm > epsilon)
        //	grad_d /= (norm);

        array_1d<double, 3> exact_distance = rDistances;
        array_1d<double, 3> abs_distance = ZeroVector(3);


        //KRATOS_WATCH("one element IS in the intefase")
        if ((rDistances(0)*rDistances(1))<0.0) //edge 12 is cut
            cut_edges[0]=true;
        else
            cut_edges[0]=false;
        if ((rDistances(1)*rDistances(2))<0.0) //edge 23 is cut.
            cut_edges[1]=true;
        else
            cut_edges[1]=false;
        if ((rDistances(2)*rDistances(0))<0.0) //edge 23 is cut.
            cut_edges[2]=true;
        else
            cut_edges[2]=false;



        //'TRICK' TO AVOID HAVING THE INTERFASE TOO CLOSE TO THE NODES:
        //since we cannot collapse node because we have to contemplate the possibility of discontinuities, we will move a little the intefase so that it is not that close.
        const double unsigned_distance0=fabs(rDistances(0));
        const double unsigned_distance1=fabs(rDistances(1));
        const double unsigned_distance2=fabs(rDistances(2));
        //we begin by finding the largest distance:
        double longest_distance=fabs(unsigned_distance0);
        if (unsigned_distance1>longest_distance)
            longest_distance=unsigned_distance1;
        if (unsigned_distance2>longest_distance)
            longest_distance=unsigned_distance2;
        //Now we set a maximum relative distance
        const double tolerable_distance =longest_distance*0.001;	// (1/1,000,000 seems to have good results)
        //and now we check if a distance is too small:
        if (unsigned_distance0<tolerable_distance)
            rDistances[0]=tolerable_distance*(rDistances[0]/fabs(rDistances[0]));
        if (unsigned_distance1<tolerable_distance)
            rDistances[1]=tolerable_distance*(rDistances[1]/fabs(rDistances[1]));
        if (unsigned_distance2<tolerable_distance)
            rDistances[2]=tolerable_distance*(rDistances[2]/fabs(rDistances[2]));
        //END OF TRICK. REMEMBER TO OVERWRITE THE DISTANCE VARIABLE IN THE ELEMENT IN CASE THESE LINES HAVE MODIFIED THEM (distances)


        //for (int jj = 0; jj < 3; jj++)
        //	KRATOS_WATCH(rDistances(jj));
        for (unsigned int i=0; i<3; i++) //we go over the 3 edges:
        {
            int edge_begin_node=i;
            int edge_end_node=i+1;
            if (edge_end_node==3) edge_end_node=0; //it's a triangle, so node 3 is actually node 0

            if(cut_edges(i)==true)
            {
                aux_nodes_relative_locations(i)=fabs(rDistances(edge_end_node)/(rDistances(edge_end_node)-rDistances(edge_begin_node) ) ) ; //position in 'natural' coordinates of edge 12, 1 when it passes over node 1. (it is over the edge 01)
                aux_nodes_father_nodes(i,0)=edge_begin_node;
                aux_nodes_father_nodes(i,1)=edge_end_node;
            }
            else
            {
                if(fabs(rDistances(edge_end_node))>fabs(rDistances(edge_begin_node))) //if edge is not cut, we collapse the aux node into the node which has the highest absolute value to have "nicer" (less "slivery") subelements
                {
                    aux_nodes_relative_locations(i)=0.0;
                    aux_nodes_father_nodes(i,0)=edge_end_node;
                    aux_nodes_father_nodes(i,1)=edge_end_node;
                }
                else
                {
                    aux_nodes_relative_locations(i)=1.0;
                    aux_nodes_father_nodes(i,0)=edge_begin_node;
                    aux_nodes_father_nodes(i,1)=edge_begin_node;
                }
            }

            //and we save the coordinate of the new aux nodes:
            for (unsigned int j=0; j<2; j++)	//x,y coordinates
                aux_points(i,j)= rOriginalPoints(edge_begin_node,j) * aux_nodes_relative_locations(i) + rOriginalPoints(edge_end_node,j) * (1.0- aux_nodes_relative_locations(i));
        }

        //having gathered all the points location in local coordinates, we will proceed to move everything to a local system. the reference point will be the first point of aux_points.
        double x_reference=aux_points(0,0);
        double y_reference=aux_points(0,1);
        double cosinus=0.0;
        double sinus=0.0;
        if (cut_edges[0]==false) //then the segment is defined by aux_points 1 and 2. so the values used to initialize were wrong. changing them:
        {
            x_reference=aux_points(1,0);
            y_reference=aux_points(1,1);
            const double one_over_interfase_lenght = 1.0/sqrt( pow((aux_points(2,0) - x_reference),2) + pow((aux_points(2,1) - y_reference),2));
            cosinus = (aux_points(2,0) - x_reference)*one_over_interfase_lenght;
            sinus = - (aux_points(2,1) - y_reference)*one_over_interfase_lenght; //WARNING; THERE IS A MINUS IN FRONT TO GET THE INVERSE ROTATION (FROM REFERENCE TO LOCAL)
        }
        else if(cut_edges[1]==true)
        {
            const double one_over_interfase_lenght = 1.0/sqrt( pow((aux_points(1,0) - x_reference),2) + pow((aux_points(1,1) - y_reference),2));
            cosinus = (aux_points(1,0) - x_reference)*one_over_interfase_lenght;
            sinus = - (aux_points(1,1) - y_reference)*one_over_interfase_lenght; //WARNING; THERE IS A MINUS IN FRONT TO GET THE INVERSE ROTATION (FROM REFERENCE TO LOCAL)
        }
        else
        {
            const double one_over_interfase_lenght = 1.0/sqrt( pow((aux_points(2,0) - x_reference),2) + pow((aux_points(2,1) - y_reference),2));
            cosinus = (aux_points(2,0) - x_reference)*one_over_interfase_lenght;
            sinus = - (aux_points(2,1) - y_reference)*one_over_interfase_lenght; //WARNING; THERE IS A MINUS IN FRONT TO GET THE INVERSE ROTATION (FROM REFERENCE TO LOCAL)
        }

        for (unsigned int i=0; i<3; i++) //we go over the 3 nodes and 3 aux nodes to move them to the new coordinates:
        {
            rRotatedPoints(i,0)= cosinus * (rOriginalPoints(i,0)-x_reference) - sinus * (rOriginalPoints(i,1)-y_reference);
            rRotatedPoints(i,1)= cosinus * (rOriginalPoints(i,1)-y_reference) + sinus * (rOriginalPoints(i,0)-x_reference);

            //and we directly change the aux points, anyway they are used only locally so it's fine:
            double aux_x_coord=aux_points(i,0);
            aux_points(i,0)= cosinus * (aux_x_coord-x_reference) - sinus * (aux_points(i,1)-y_reference);
            aux_points(i,1)= cosinus * (aux_points(i,1)-y_reference) + sinus * (aux_x_coord-x_reference);
        }

        rRotationMatrix(0,0)=cosinus;
        rRotationMatrix(0,1)= sinus;
        rRotationMatrix(1,0)= -sinus;
        rRotationMatrix(1,1)=cosinus;

        double temp_area;
        //to calculate the new rigidity matrix in local coordinates, the element will need the derivated in the rotated axis and the rotation matrix:
        CalculateGeometryData(rRotatedPoints, DN_DX_in_local_axis, temp_area);



        //KRATOS_WATCH(rRotatedPoints)
        //KRATOS_WATCH(aux_points)

        array_1d<double, 2 > grad_d;
        noalias(grad_d) = prod(trans(DN_DX_in_local_axis), rDistances);


        array_1d<double,2> base_point;
        if (cut_edges(0)==true) // it means it is a cut edge, if it was 0.0 or 1.0 then it would be an uncut edge
        {
            base_point[0] = aux_points(0,0);
            base_point[1] = aux_points(0,1);
        }
        else //it means aux_point 0 is a clone of other point, so we go to the second edge.
        {
            base_point[0] = aux_points(1,0);
            base_point[1] = aux_points(1,1);
        }

        for (int i_node = 0; i_node < 3; i_node++)
        {
            double d =    (rRotatedPoints(i_node,0) - base_point[0]) * grad_d[0] +
                          (rRotatedPoints(i_node,1) - base_point[1]) * grad_d[1] ;
            abs_distance[i_node] = fabs(d);
        }

        //assign correct sign to exact distance
        for (int i = 0; i < 3; i++)
        {
            if (rDistances[i] < 0.0)
            {
                exact_distance[i] = -abs_distance[i];
                --most_common_sign;
            }
            else
            {
                exact_distance[i] = abs_distance[i];
                ++most_common_sign;
            }
        }


        //compute exact distance gradients
        array_1d<double, 2 > exact_distance_gradient;
        noalias(exact_distance_gradient) = prod(trans(DN_DX_in_local_axis), exact_distance);

        array_1d<double, 2 > abs_distance_gradient;
        noalias(abs_distance_gradient) = prod(trans(DN_DX_in_local_axis), abs_distance);


        double max_aux_dist_on_cut = -1;
        for (int edge = 0; edge < 3; edge++)
        {
            const int i = edge;
            int j = edge+1;
            if (j==3) j=0;
            if (rDistances[i] * rDistances[j] < 0.0)
            {
                const double tmp = fabs(rDistances[i]) / (fabs(rDistances[i]) + fabs(rDistances[j]));
                //compute the position of the edge node
                double abs_dist_on_cut = abs_distance[i] * tmp + abs_distance[j] * (1.00 - tmp);
                if(abs_dist_on_cut > max_aux_dist_on_cut) max_aux_dist_on_cut = abs_dist_on_cut;
            }
        }


        //we reset all data:
        rGradientsValue[0]=ZeroMatrix(2,2);
        rGradientsValue[1]=ZeroMatrix(2,2);
        rGradientsValue[2]=ZeroMatrix(2,2);
        NEnriched=ZeroMatrix(3,2);
        rGPShapeFunctionValues=ZeroMatrix(3,3);


        //now we must check the 4 created partitions of the domain.
        //one has been collapsed, so we discard it and therefore save only one.
        unsigned int partition_number=0;		//
        //the 3 first partitions are  created using 2 auxiliary nodes and a normal node. at least one of these will be discarded due to zero area
        //the last one is composed by the 3 auxiliary nodes. it 'looks' wrong, but since at least one has been collapsed, it actually has a normal node.
        bool found_empty_partition=false;
        for (unsigned int i=0; i<4; i++) //i partition
        {
            unsigned int j_aux = i + 2;
            if (j_aux>2) j_aux -= 3;
            BoundedMatrix<int,3,2> partition_father_nodes;
            array_1d<double,3> N;
            if (i<3)
            {
                partition_father_nodes(0,0)=i;
                partition_father_nodes(0,1)=i;
                partition_father_nodes(1,0)=aux_nodes_father_nodes(i,0); //we are using i aux node
                partition_father_nodes(1,1)=aux_nodes_father_nodes(i,1); //we are using i aux node
                partition_father_nodes(2,0)=aux_nodes_father_nodes(j_aux,0); //we are using j_aux node
                partition_father_nodes(2,1)=aux_nodes_father_nodes(j_aux,1); //we are using j_aux node

                coord_subdomain(0,0)=rRotatedPoints(i,0);
                coord_subdomain(0,1)=rRotatedPoints(i,1);
                coord_subdomain(1,0)=aux_points(i,0);
                coord_subdomain(1,1)=aux_points(i,1);
                coord_subdomain(2,0)=aux_points(j_aux,0);
                coord_subdomain(2,1)=aux_points(j_aux,1);
            }
            else
            {
                //the last partition, made by the 3 aux nodes.
                partition_father_nodes=aux_nodes_father_nodes;
                coord_subdomain=aux_points;
                //found_last_partition=true;
            }
            //calculate data of this partition

            CalculateGeometryData(coord_subdomain, DN_DX_subdomain, temp_area);
            if (temp_area > 1.0e-20) //ok, it does not have zero area
            {
                rVolumes(partition_number)=temp_area;
                //we look for the gauss point of the partition:
                double x_GP_partition =  one_third * ( coord_subdomain(0,0) + coord_subdomain(1,0) + coord_subdomain(2,0) );
                double y_GP_partition =  one_third * ( coord_subdomain(0,1) + coord_subdomain(1,1) + coord_subdomain(2,1) );
                double z_GP_partition  = 0.0;
                //we reset the coord_subdomain matrix so that we have the whole element again:
                coord_subdomain = rRotatedPoints;
                //and we calculate its shape function values
                CalculatePosition ( coord_subdomain , x_GP_partition ,y_GP_partition ,z_GP_partition , N);
                //we check the partition sign.
                const double partition_sign = (N(0)*rDistances(0) + N(1)*rDistances(1) + N(2)*rDistances(2))/fabs(N(0)*rDistances(0) + N(1)*rDistances(1) + N(2)*rDistances(2));
                //rPartitionsSign(partition_number)=partition_sign;

                rGPShapeFunctionValues(partition_number,0)=N(0);
                rGPShapeFunctionValues(partition_number,1)=N(1);
                rGPShapeFunctionValues(partition_number,2)=N(2);

                //compute enriched shape function values
                double dist = 0.0;
                double abs_dist = 0.0;
                for (int j = 0; j < 3; j++)
                {
                    dist += rGPShapeFunctionValues(partition_number, j) * exact_distance[j];
                    abs_dist += rGPShapeFunctionValues(partition_number, j) * abs_distance[j];
                }

                if (partition_sign < 0.0)
                    rPartitionsSign[partition_number] = -1.0;
                else
                    rPartitionsSign[partition_number] = 1.0;

                //enrichment for the gradient (continuous pressure, discontinuous gradient. Coppola-Codina enrichment
                NEnriched(partition_number, 0) = 0.5/max_aux_dist_on_cut * (abs_dist - rPartitionsSign[partition_number] * dist);
                for (int j = 0; j < 2; j++)
                {
                    rGradientsValue[partition_number](0, j) = (0.5/max_aux_dist_on_cut) * (abs_distance_gradient[j] - rPartitionsSign[partition_number] * exact_distance_gradient[j]);
                }

                //enrichment to add a constant discontinuity along the interfase.
                if (rPartitionsSign[partition_number]*most_common_sign > 0.0) //on this side we will copy the standard shape functions. (maybe we change the sign, but keeping the sign
                {
                    NEnriched(partition_number, 1) = -1.0*rPartitionsSign[partition_number]*NEnriched(partition_number, 0) ;
                    for (int j = 0; j < 2; j++)
                    {
                        rGradientsValue[partition_number](1, j) = -1.0*rPartitionsSign[partition_number]*rGradientsValue[partition_number](0, j);
                    }
                }
                else //we have to construct the shape function to guarantee a constant jump to 2:
                {
                    //notice that the partition to be changed must one the subdomains created with a real node. so the i index is also the real node. (used 4 lines below to recover the distance.

                    NEnriched(partition_number, 1) = rPartitionsSign[partition_number]*( 1.0*NEnriched(partition_number, 0) - 0.6666666666666666666666666666 ) ;
                    for (int j = 0; j < 2; j++)
                    {
                        rGradientsValue[partition_number](1, j) =  (rPartitionsSign[partition_number]*1.0*rGradientsValue[partition_number](0, j) - exact_distance_gradient[j]*1.0/(abs_distance[i]));
                    }
                }

                partition_number++;

            }
            else
                found_empty_partition=true;

        }
        if (found_empty_partition==false)
            KRATOS_WATCH("WROOOONGGGGGGGGGGG");
        return 3;
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

        //            std::cout << "splitting edge" << i << " " << j << std::endl;
        //            KRATOS_WATCH(Volume1Id);
        //            KRATOS_WATCH(Volume2Id);
        //            KRATOS_WATCH(rShapeFunctionValues)

        double delta1 = rShapeFunctionValues(Volume1Id, i) * (1.00 - division_i);
        rShapeFunctionValues(Volume1Id, i) += delta1;
        rShapeFunctionValues(Volume1Id, j) -= delta1;

        double delta2 = rShapeFunctionValues(Volume2Id, i) * (1.00 - division_j);
        rShapeFunctionValues(Volume2Id, j) += delta2;
        rShapeFunctionValues(Volume2Id, i) -= delta2;

        //            KRATOS_WATCH(delta1)
        //            KRATOS_WATCH(delta2)
        //
        //            KRATOS_WATCH(rShapeFunctionValues)



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

    static double ComputeSubTetraVolumeAndCenter(const BoundedMatrix<double, 3, 8 > & aux_coordinates,
            array_1d<double, 3 > & center_position,
            const int i0, const int i1, const int i2, const int i3)
    {
        double x10 = aux_coordinates(i1, 0) - aux_coordinates(i0, 0); //geom[1].X() - geom[0].X();
        double y10 = aux_coordinates(i1, 1) - aux_coordinates(i0, 1); // geom[1].Y() - geom[0].Y();
        double z10 = aux_coordinates(i1, 2) - aux_coordinates(i0, 2); // geom[1].Z() - geom[0].Z();

        double x20 = aux_coordinates(i2, 0) - aux_coordinates(i0, 0); // geom[2].X() - geom[0].X();
        double y20 = aux_coordinates(i2, 1) - aux_coordinates(i0, 1); // geom[2].Y() - geom[0].Y();
        double z20 = aux_coordinates(i2, 2) - aux_coordinates(i0, 2); // geom[2].Z() - geom[0].Z();

        double x30 = aux_coordinates(i3, 0) - aux_coordinates(i0, 0); // geom[3].X() - geom[0].X();
        double y30 = aux_coordinates(i3, 1) - aux_coordinates(i0, 1); // geom[3].Y() - geom[0].Y();
        double z30 = aux_coordinates(i3, 2) - aux_coordinates(i0, 2); // geom[3].Z() - geom[0].Z();

        double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
        double vol = detJ * 0.1666666666666666666667;

        for (unsigned int i = 0; i < 3; i++)
        {
            center_position[i] = aux_coordinates(i0, i) + aux_coordinates(i1, i) + aux_coordinates(i2, i) + aux_coordinates(i3, i);
        }
        center_position *= 0.25;

        return vol;
    }

    template<class TMatrixType>
    static void ComputeElementCoordinates(array_1d<double, 4 > & N, const array_1d<double, 3 > & center_position, const TMatrixType& rPoints, const double vol)
    {
        double x0 = rPoints(0, 0); //geom[0].X();
        double y0 = rPoints(0, 1); //geom[0].Y();
        double z0 = rPoints(0, 2); //geom[0].Z();
        double x1 = rPoints(1, 0); //geom[1].X();
        double y1 = rPoints(1, 1); //geom[1].Y();
        double z1 = rPoints(1, 2); //geom[1].Z();
        double x2 = rPoints(2, 0); //geom[2].X();
        double y2 = rPoints(2, 1); //geom[2].Y();
        double z2 = rPoints(2, 2); //geom[2].Z();
        double x3 = rPoints(3, 0); //geom[3].X();
        double y3 = rPoints(3, 1); //geom[3].Y();
        double z3 = rPoints(3, 2); //geom[3].Z();

        double xc = center_position[0];
        double yc = center_position[1];
        double zc = center_position[2];

        double inv_vol = 1.0 / vol;
        //            N[0] = CalculateVol(x1, y1, z1, x3, y3, z3, x2, y2, z2, xc, yc, zc) * inv_vol;
        //            N[1] = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, xc, yc, zc) * inv_vol;
        //            N[2] = CalculateVol(x3, y3, z3, x1, y1, z1, x0, y0, z0, xc, yc, zc) * inv_vol;
        //            N[3] = CalculateVol(x3, y3, z3, x0, y0, z0, x2, y2, z2, xc, yc, zc) * inv_vol;
        N[0] = CalculateVol(x1, y1, z1, x3, y3, z3, x2, y2, z2, xc, yc, zc) * inv_vol;
        N[1] = CalculateVol(x0, y0, z0, x2, y2, z2, x3, y3, z3, xc, yc, zc) * inv_vol;
        N[2] = CalculateVol(x3, y3, z3, x1, y1, z1, x0, y0, z0, xc, yc, zc) * inv_vol;
        N[3] = CalculateVol(x1, y1, z1, x2, y2, z2, x0, y0, z0, xc, yc, zc) * inv_vol;

    }

    static inline double CalculateVol(const double x0, const double y0, const double z0,
                                      const double x1, const double y1, const double z1,
                                      const double x2, const double y2, const double z2,
                                      const double x3, const double y3, const double z3
                                     )
    {
        double x10 = x1 - x0;
        double y10 = y1 - y0;
        double z10 = z1 - z0;

        double x20 = x2 - x0;
        double y20 = y2 - y0;
        double z20 = z2 - z0;

        double x30 = x3 - x0;
        double y30 = y3 - y0;
        double z30 = z3 - z0;

        double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
        return detJ * 0.1666666666666666666667;
    }

    //2d
    static inline void CalculateGeometryData(
        const BoundedMatrix<double, 3, 3 > & coordinates,
        BoundedMatrix<double,3,2>& DN_DX,
        array_1d<double,3>& N,
        double& Area)
    {
        double x10 = coordinates(1,0) - coordinates(0,0);
        double y10 = coordinates(1,1) - coordinates(0,1);

        double x20 = coordinates(2,0) - coordinates(0,0);
        double y20 = coordinates(2,1) - coordinates (0,1);

        //Jacobian is calculated:
        //  |dx/dxi  dx/deta|	|x1-x0   x2-x0|
        //J=|				|=	|			  |
        //  |dy/dxi  dy/deta|	|y1-y0   y2-y0|


        double detJ = x10 * y20-y10 * x20;

        DN_DX(0,0) = -y20 + y10;
        DN_DX(0,1) = x20 - x10;
        DN_DX(1,0) =  y20	   ;
        DN_DX(1,1) = -x20     ;
        DN_DX(2,0) = -y10	   ;
        DN_DX(2,1) = x10	   ;

        DN_DX /= detJ;
        N[0] = 0.333333333333333;
        N[1] = 0.333333333333333;
        N[2] = 0.333333333333333;

        Area = 0.5*detJ;
    }

    //template<class TMatrixType, class TVectorType, class TGradientType>
    static inline double CalculateVolume2D(
        const BoundedMatrix<double, 3, 3 > & coordinates)
    {
        double x10 = coordinates(1,0) - coordinates(0,0);
        double y10 = coordinates(1,1) - coordinates(0,1);

        double x20 = coordinates(2,0) - coordinates(0,0);
        double y20 = coordinates(2,1) - coordinates (0,1);
        double detJ = x10 * y20-y10 * x20;
        return 0.5*detJ;
    }

    static inline bool CalculatePosition(const BoundedMatrix<double, 3, 3 > & coordinates,
                                         const double xc, const double yc, const double zc,
                                         array_1d<double, 3 > & N
                                        )
    {
        double x0 = coordinates(0,0);
        double y0 = coordinates(0,1);
        double x1 = coordinates(1,0);
        double y1 = coordinates(1,1);
        double x2 = coordinates(2,0);
        double y2 = coordinates(2,1);

        double area = CalculateVol(x0, y0, x1, y1, x2, y2);
        double inv_area = 0.0;
        if (area == 0.0)
        {
            KRATOS_THROW_ERROR(std::logic_error, "element with zero area found", "");
        }
        else
        {
            inv_area = 1.0 / area;
        }


        N[0] = CalculateVol(x1, y1, x2, y2, xc, yc) * inv_area;
        N[1] = CalculateVol(x2, y2, x0, y0, xc, yc) * inv_area;
        N[2] = CalculateVol(x0, y0, x1, y1, xc, yc) * inv_area;
        //KRATOS_WATCH(N);

        if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0) //if the xc yc is inside the triangle return true
            return true;

        return false;
    }

    static inline double CalculateVol(const double x0, const double y0,
                                      const double x1, const double y1,
                                      const double x2, const double y2
                                     )
    {
        return 0.5 * ((x1 - x0)*(y2 - y0)- (y1 - y0)*(x2 - x0));
    }

    static inline void CalculateGeometryData(
        const BoundedMatrix<double, 3, 3 > & coordinates,
        BoundedMatrix<double,3,2>& DN_DX,
        double& Area)
    {
        double x10 = coordinates(1,0) - coordinates(0,0);
        double y10 = coordinates(1,1) - coordinates(0,1);

        double x20 = coordinates(2,0) - coordinates(0,0);
        double y20 = coordinates(2,1) - coordinates (0,1);

        //Jacobian is calculated:
        //  |dx/dxi  dx/deta|	|x1-x0   x2-x0|
        //J=|				|=	|			  |
        //  |dy/dxi  dy/deta|	|y1-y0   y2-y0|


        double detJ = x10 * y20-y10 * x20;

        DN_DX(0,0) = -y20 + y10;
        DN_DX(0,1) = x20 - x10;
        DN_DX(1,0) =  y20	   ;
        DN_DX(1,1) = -x20     ;
        DN_DX(2,0) = -y10	   ;
        DN_DX(2,1) = x10	   ;

        DN_DX /= detJ;

        Area = 0.5*detJ;
    }


    template <int TDim>
    static void AddToEdgeAreas(array_1d<double, (TDim-1)*3 >& edge_areas,
                               const array_1d<double, TDim+1 >& exact_distance,
                               const array_1d<double, TDim+1 >& indices,
                               const double sub_volume)
    {
        //check if the element has 3 "nodes" on the cut surface and if the remaining one is positive
        //to do so, remember that edge nodes are marked with an id greater than 3
        unsigned int ncut=0, pos=0, positive_pos=0;
        for(unsigned int i=0; i<TDim+1; i++)
        {
            if(indices[i] > TDim) ncut++;
            else if(exact_distance[indices[i]] > 0)
            {
                positive_pos = indices[i];
                pos++;
            }
        }

        if(ncut == TDim && pos==1) //cut face with a positive node!!
        {
            double edge_area = sub_volume*3.0/fabs(exact_distance[positive_pos]);
            edge_area /= static_cast<double>(TDim);
            for(unsigned int i=0; i<TDim+1; i++)
            {
                if( indices[i] > TDim)
                {
                    int edge_index = indices[i] - TDim - 1;
                    edge_areas[edge_index] += edge_area;
                }
            }
        }
    }


};

} // namespace Kratos.

#endif // KRATOS_ENRICHMENT_UTILITIES_INCLUDED  defined


