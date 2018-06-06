//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pablo Becker
//                    
//


#if !defined(KRATOS_ENRICHMENT_UTILITIES_DUPLICATE_DOFS_INCLUDED )
#define  KRATOS_ENRICHMENT_UTILITIES_DUPLICATE_DOFS_INCLUDED



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

/** This function implements the calculation of enrichment functions, where the position of the discontinuity
* is prescribed by giving a level set function prescribing the distance the enrichmenet surface.
* WARNING: the function works only for simplicial elements (Triangles and tetrahedras)
* TODO: only implemented for Tetrahedras so far
*
* The essential idea is that the element is divided by the level set functions into 2 parts one + and the other -
* let's suppose for a moment that nodes 1 and 2 (in local numeration) lay on the negative side and nodes 3,4 on the positive side
* if N1, N2, N3, N4 are the shape functions "centered" in the node with corresponding id, the enrichment shape functions are defined as
* follows:
*
* POSITIVE SIDE
* N1enriched = N1, N2enriched = N2, N3enriched=0, N4enriched=0
*
* NEGATIVE SIDE
* N1enriched = 0, N2enriched = 0, N3enriched=N3, N4enriched=N4
*
* the combination of the original element shape functions and of thos enrichements can model a discontinuity and a discontinuity of the gradient
* of the function to be enriched.
*
* WARNING: difficulties can be expected if the division prescribed passes exactly in one node or if one of the two sides has a volume
* approximately 0. In this case the definition above works correctly but the use of the enrichment shape functions requires some care

*
*
*/
class EnrichmentUtilitiesDuplicateDofs
{
public:

    /**
     * The method to calculate the ernriched shape functions for given tetrahedra.
     * 4 ENRICHED SHAPE FUNCTIONS ARE GENERATED!!
     *
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
     *        Each element of vector is a 4*3 matrix representing the gradient of enriched shape functions. The use of
     *        matrix is for possible future improvement.
     * @param Nenriched is a Matrix that contains for every gauss point the values of the enriched shape functions at the position of the gauss point
     *        so that Nenriched(1,0) contains the value of the enriched shape function "0" at the gauss point "1"
     * @return number of partitions created which can be from 1 to 6.
     *         1 holds for only 1 partition which is the original element. (No partitioning needed)
     */
    template<class TMatrixType, class TVectorType, class TGradientType>
    static int CalculateTetrahedraEnrichedShapeFuncions(TMatrixType const& rPoints, TGradientType const& DN_DX,
            TVectorType& rDistances, TVectorType& rVolumes, TMatrixType& rShapeFunctionValues,
            TVectorType& rPartitionsSign, std::vector<TMatrixType>& rGradientsValue, TMatrixType& NEnriched)
    {
        KRATOS_TRY

        const int n_nodes = 4; // it works only for tetrahedra
        const int n_edges = 6; // it works only for tetrahedra

        const int edge_i[] = {0, 0, 0, 1, 1, 2};
        const int edge_j[] = {1, 2, 3, 2, 3, 3};


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

//         for (int i = 0; i < 6; i++)
//             for (int j = 0; j < n_nodes; j++)
//                 rShapeFunctionValues(i, j) = 0.25;
        
                




        array_1d<double, n_nodes> exact_distance = rDistances;
        array_1d<double, n_nodes> abs_distance = ZeroVector(n_nodes);
        //double sub_volumes_sum = 0.00;

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

                //modify the distances to avoid zeros
        for(unsigned int i=0; i<4;i++)
            if(fabs(rDistances[i]) < 1e-4*max_lenght) rDistances[i]=1e-4*max_lenght;

        
        
        for (unsigned int i = 0; i < 4; i++)
            abs_distance[i] = fabs(rDistances[i]);

        //compute the gradient of the distance and normalize it
        array_1d<double, 3 > grad_d;
        noalias(grad_d) = prod(trans(DN_DX), rDistances);
        double norm = norm_2(grad_d);
        
        if(norm < 1e-10) norm=1e-10;
        
        grad_d /= (norm);

    //now decide splitting pattern
    for (int edge = 0; edge < n_edges; edge++)
    {
        const int i = edge_i[edge];
        const int j = edge_j[edge];
        if (rDistances[i] * rDistances[j] < 0.0)
        {
            const double tmp = fabs(rDistances[i]) / (fabs(rDistances[i]) + fabs(rDistances[j]));

            split_edge[edge + 4] = new_node_id;
            edge_division_i[edge] = tmp;
            edge_division_j[edge] = 1.00 - tmp;

            //compute the position of the edge node
            for (unsigned int k = 0; k < 3; k++)
                aux_coordinates(new_node_id, k) = rPoints(i, k) * edge_division_j[edge] + rPoints(j, k) * edge_division_i[edge];

            new_node_id++;

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
//         rVolumes.resize(1,false)
//         rVolumes[0] = volume;
//        sub_volumes_sum = volume;
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

//         KRATOS_/*WATCH*/("line 358")


        if(min_distance < 0.0)
            rPartitionsSign[0] = -1.0;
        else
            rPartitionsSign[0] = 1.0;

        number_of_partitions = 1;

        if(rShapeFunctionValues.size1() != 4 || rShapeFunctionValues.size2() != 4)
            rShapeFunctionValues.resize(4,4,false);

//         for (int j = 0; j < 4; j++)
//             rShapeFunctionValues(0, j) = 0.25;
//         KRATOS_WATCH("line 374")

        if(NEnriched.size1() != 0 || NEnriched.size2() != 0)
            NEnriched.resize(0,0,false);
//         KRATOS_WATCH("line 377")

        if(rGradientsValue.size() !=0)
            rGradientsValue.resize(0);
//         KRATOS_WATCH("line 382")

//             for (int j = 0; j < number_of_partitions; j++)
//                 NEnriched(j, 0) = 0.0;

//             rGradientsValue[0] = ZeroMatrix(1,3);
    }
    else //if (number_of_splitted_edges == 4)
    {
        //define the splitting mode for the tetrahedra
        int edge_ids[6];
        TetrahedraSplit::TetrahedraSplitMode(split_edge, edge_ids);
        int nel; //number of elements generated
        int n_splitted_edges; //number of splitted edges
        int nint; //number of internal nodes
        int t[56];
        TetrahedraSplit::Split_Tetrahedra(edge_ids, t, &nel, &n_splitted_edges, &nint);

//         KRATOS_WATCH("line 395")


        if (nint != 0)
            KRATOS_THROW_ERROR(std::logic_error, "requiring an internal node for splitting ... can not accept this", "");


        //now obtain the tetras and compute their center coordinates and volume
        if(rShapeFunctionValues.size1() != (unsigned int)(nel)*4 || rShapeFunctionValues.size2() != 4)
            rShapeFunctionValues.resize(nel*4,4,false);

        std::vector< array_1d<double, 3 > >center_position(4);
        Vector gauss_volumes(4);
        if(rVolumes.size() != (unsigned int)(nel)*4)
            rVolumes.resize(nel*4,false);

        for (int i = 0; i < nel; i++)
        {
            int i0, i1, i2, i3; //indices of the subtetrahedra
            TetrahedraSplit::TetrahedraGetNewConnectivityGID(i, t, split_edge, &i0, &i1, &i2, &i3);


//                 double sub_volume = ComputeSubTetraVolumeAndCenter(aux_coordinates, center_position, i0, i1, i2, i3);
            ComputeSubTetraVolumeAndGaussPoints(aux_coordinates, gauss_volumes, center_position, i0, i1, i2, i3);
//             KRATOS_WATCH("line 414")
            for(unsigned int k=0; k<4; k++)
            {
                rVolumes[i*4+k] = gauss_volumes[k];

                array_1d<double, 4 > N;
                ComputeElementCoordinates(N, center_position[k], rPoints, volume);

                for (int j = 0; j < 4; j++)
                    rShapeFunctionValues(i*4+k, j) = N[j];
            }

//             KRATOS_WATCH("line 426")



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
      if(NEnriched.size1() != (unsigned int)(number_of_partitions)*4 || NEnriched.size2() != 4)
            NEnriched.resize(number_of_partitions*4,4,false);

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

        if(rGradientsValue.size() != (unsigned int)(number_of_partitions)*4)
            rGradientsValue.resize(number_of_partitions*4);

	if(rPartitionsSign.size() != (unsigned int)(number_of_partitions)*4)
            rPartitionsSign.resize(number_of_partitions*4,false);

//         KRATOS_WATCH(rShapeFunctionValues)

        for (int i = 0; i < number_of_partitions*4; i++) //i is the gauss point number
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

            rGradientsValue[i] = DN_DX;

            double disc_factor;
            for( unsigned int j = 0 ; j<rShapeFunctionValues.size2(); j++)
            {
                if(exact_distance[j]*dist > 0) //shape function corresponding to a node on the same side as the gauss point
                    disc_factor = 0.0;
                else
                    disc_factor = 1.0;

                //if( fabs(dist) < relatively_close)
                //    disc_factor = 0.0;

                NEnriched(i,j) = disc_factor * rShapeFunctionValues(i,j);


                for(unsigned int k=0; k<3; k++)
                {
                    rGradientsValue[i](j, k) *= disc_factor;
                }
            }




        }
    }
//         else
//         {
//            if(NEnriched.size1() != 6 || NEnriched.size2() != 1)
//                 NEnriched.resize(1,1,false);
//
//             //compute multiplier
//             NEnriched(0,0) = 0.0;
//
//             for (int j = 0; j < 3; j++)
//                 rGradientsValue[0](0, j) = 0.0;
//         }

    return number_of_partitions;
    KRATOS_CATCH("");
}



private:


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

//compute 4 gauss point per subtetra
static void ComputeSubTetraVolumeAndGaussPoints(const BoundedMatrix<double, 3, 8 > & aux_coordinates,
        Vector& volumes,
        std::vector< array_1d<double, 3 > >& center_position,
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

    for (unsigned int i = 0; i < 4; i++)
    {
        volumes[i] = 0.25*vol;
    }

    //gauss point 0
    for (unsigned int i = 0; i < 3; i++)
    {
        center_position[0][i] = 0.58541020*aux_coordinates(i0, i) + 0.13819660*aux_coordinates(i1, i) + 0.13819660*aux_coordinates(i2, i) + 0.13819660*aux_coordinates(i3, i);
    }
    //gauss point 1
    for (unsigned int i = 0; i < 3; i++)
    {
        center_position[1][i] = 0.13819660*aux_coordinates(i0, i) + 0.58541020*aux_coordinates(i1, i) + 0.13819660*aux_coordinates(i2, i) + 0.13819660*aux_coordinates(i3, i);
    }

    //gauss point 2
    for (unsigned int i = 0; i < 3; i++)
    {
        center_position[2][i] = 0.13819660*aux_coordinates(i0, i) + 0.13819660*aux_coordinates(i1, i) + 0.58541020*aux_coordinates(i2, i) + 0.13819660*aux_coordinates(i3, i);
    }

    //gauss point 3
    for (unsigned int i = 0; i < 3; i++)
    {
        center_position[3][i] = 0.13819660*aux_coordinates(i0, i) + 0.13819660*aux_coordinates(i1, i) + 0.13819660*aux_coordinates(i2, i) + 0.58541020*aux_coordinates(i3, i);
    }

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


};

} // namespace Kratos.

#endif // KRATOS_ENRICHMENT_UTILITIES_DUPLICATE_DOFS_INCLUDED  defined


