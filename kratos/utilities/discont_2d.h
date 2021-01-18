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


#if !defined(KRATOS_DISCONTINUOUS_2D_UTILITIES_INCLUDED )
#define  KRATOS_DISCONTINUOUS_2D_UTILITIES_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>
#include <limits>

// External includes


// Project includes
#include "includes/define.h"
//#include "utilities/split_tetrahedra.h"

namespace Kratos
{

/** @brief This utility can be used to calculate the enriched shape function for tetrahedra element.
 *  @details The metodology consists in partitioning the tetrahedra in a set of sub-tetrahedra and
 *  cacluate the enrichment information using these partitions.
 */
class DiscontinuousShapeFunctionsUtilities_2D
{
public:

    /**
     * @brief The method to calculate the enriched shape functions for given triangle
     * @details Basically, two shape functions are provided, 
     * -# One is the enrichment to capturate discontinuities in the gradients of the pressure 
     *       ("Improving Eulerian two-phase flow finite element approximation with discontinuous gradient (i.e pressure) shape functions" Coppola-Owen and Codina) 
     * -# The second one is to capturate discontinuities in the varialbe(ie. pressure). it is a shape function that is zero on the nodes and has a constant discontinuity along the found interfase)
     * 
     * @param rPoints A 3x3 matrix where row i has the coordinates of node i.
     * @param DN_DX The gradient of the shape functions Ni respect to the reference coordinates
     * @param rDistances is an input  vector of 3 size which holds relative distance for each node.
     *        it is used internally to mark the position of the zero level
     * @param rVolumes Result vector with size (3+1) (maximum number of partitions) holding the volume of each partition. 
     * 		  Partitions are given in the following order: first the one that is "alone". meaning that it is on the side of the interfase where there's only a triangle
     * 		   and the other 2 are the ones in the other side. So the paritition signs will be allways -1,1,1 or 1,-1,-1
     * 		  Note:volume 4 is not a partition itself, but is added to be used in the element. it would be the result of creating the 2 triangles of the cuadrilateral paritition with the opposite edge as the current configuration (joining node k_aux with node_4 instead of the current node_j-node5)
     * @param rShapeFunctionValues Result 3x3 matrix where each row represents a partition and holds the shape functions N1 to N3 ( the cut)
     *        of the original triangle evaluated in the gauss point (center) of the partition.
     *        so that it is  N(gauss_index, node_index)
     * @param rPartitionsSign A result vector of 3 holding the sign of the distance for the partition. Read rVolumes 5 lines above
     *        The value -1 represents the negative distance sign, 1 represents positive distance and 0 stands for not used partition
     * @param rGradientsValue Restult vector of size 3 holding the gradient of the enriched shape funciton for each volume.
     *        Each element of vector is a 1x3 matrix representing the gradient of enriched shape function. The use of
     *        matrix is for possible future improvement.
     * @param Nenriched is a Matrix that contains for every gauss point the values of the enriched shape functions at the position of the gauss point
     *        so that Nenriched(1,0) contains the value of the enriched shape function "0" at the gauss point "1"
     * @param face_gauss_N is the location of the (single) integration point of the interfase: its midpoint.
     * @param face_gauss_N_enrich is the value of the enrichment shape functions in the integration point
     * 		  actually it's value is always the same so no need to use it: the shape functions were defined to make it 1 in the first shape function,
     * 		  And 1 and -1 the second shape function (it's discontinous, so it has these two values in the interfase)
     * 		  @warning Therefore the discontinuity in the shape function is equal to 2.
     * @param type_of_cut The partition that is 'alone': the one that is on one side of the shape function
     * 		  the other two are the ones in the other side, meaning they have the same derivatives and , for example, densities.
     * 
     * @return number of partitions created which can be from 1 to 3.
     *         1 holds for only 1 partition which is the original element. (No partitioning needed)
     */
   
    //with some added vectors, to be used when we need information about the interfase between the 2 elements
    //template<class TMatrixType, class TVectorType, class TGradientType>
    /*
    static int CalculateTriangleDiscontinuousShapeFunctions(BoundedMatrix<double,(TDim+1), TDim >& rPoints, BoundedMatrix<double, (TDim+1), TDim >& DN_DX,
            array_1d<double,(TDim+1)>& rDistances, array_1d<double,(3*(TDim-1))>& rVolumes, BoundedMatrix<double, 3*(TDim-1), (TDim+1) >& rGPShapeFunctionValues,
            array_1d<double,(3*(TDim-1))>& rPartitionsSign, std::vector<TMatrixType>& rGradientsValue, BoundedMatrix<double,3*(TDim-1), (TDim+1)>& NEnriched, //and information about the interfase:
            array_1d<double,(3)>& face_gauss_N, array_1d<double,(3)>& face_gauss_Nenriched, double& face_Area, array_1d<double,(3)>& face_n ,unsigned int& type_of_cut)    
            *
/*
    static int CalculateTriangleDiscontinuousShapeFunctions(BoundedMatrix<double,(2+1), 2 >& rPoints, BoundedMatrix<double, (2+1), 2 >& DN_DX,
            array_1d<double,(2+1)>& rDistances, array_1d<double,(3*(2-1))>& rVolumes, BoundedMatrix<double, 3*(2-1), (2+1) >& rGPShapeFunctionValues,
            array_1d<double,(3*(2-1))>& rPartitionsSign, std::vector<Matrix>& rGradientsValue, BoundedMatrix<double,3*(2-1), (2+1)>& NEnriched, //and information about the interfase:
            array_1d<double,(3)>& face_gauss_N, array_1d<double,(3)>& face_gauss_Nenriched, double& face_Area, array_1d<double,(3)>& face_n ,unsigned int& type_of_cut)   
 
   */ 
	template<class TMatrixType, class TVectorType, class TGradientType>
	static int CalculateTriangleDiscontinuousShapeFunctions(TMatrixType const& rPoints, TGradientType const& DN_DX,
            TVectorType rDistances, TVectorType& rVolumes, TMatrixType& rGPShapeFunctionValues,
            TVectorType& rPartitionsSign, std::vector<TMatrixType>& rGradientsValue, TMatrixType& Nenriched, //and information about the interfase:
            TVectorType& face_gauss_N, TVectorType& face_gauss_Nenriched, double& face_Area, TVectorType& face_n ,unsigned int& type_of_cut)    
    
	{
        KRATOS_TRY
	
		//unsigned int i,j,k;
		//unsigned int i_aux,j_aux,k_aux; //
		type_of_cut = 0;   // 0 means no cuts, 1 means element is cut through edges ij,ik;    2 ij,jk ;    3 ik , kj ;   INTERFASES ON nodes are not contemplated   
		const double one_third=1.0/3.0;
		BoundedMatrix<double,3,2> aux_points; //for auxiliary nodes 4(between 1 and 2) ,5(between 2 and 3) ,6 (between 3 and 1)
		BoundedMatrix<double, 3, 2 > coord_subdomain; //used to pass arguments when we must calculate areas, shape functions, etc
		BoundedMatrix<double,3,2> DN_DX_subdomain; //used to retrieve derivatives
		
		
		double Area;//area of the complete element
		rGPShapeFunctionValues(0,0)=one_third; rGPShapeFunctionValues(0,1)=one_third; rGPShapeFunctionValues(0,2)=one_third; //default, when no interfase has been found
		Area = CalculateVolume2D( rPoints );
		array_1d<bool,3> cut_edges;
		array_1d<double,3> aux_nodes_relative_locations;
		BoundedMatrix<int,3,2> aux_nodes_father_nodes;

        //to begin with we must check whether our element is cut or not by the interfase.
        if( (rDistances(0)*rDistances(1))>0.0 && (rDistances(0)*rDistances(2))>0.0 ) //it means that this element IS NOT cut by the interfase. we must return data of a normal, non-enriched element
		{
			rVolumes(0)=Area;
			rGPShapeFunctionValues(0,0)=one_third; rGPShapeFunctionValues(0,1)=one_third; rGPShapeFunctionValues(0,2)=one_third;
			Nenriched(0,0) = 0.0;
			type_of_cut=1;
            for (int j = 0; j < 3; j++)
                rGradientsValue[0](0, j) = 0.0;
            if (rDistances(0) < 0.0) rPartitionsSign[0] = -1.0;
            else rPartitionsSign[0] = 1.0;
			//KRATOS_WATCH("one element not in the intefase")
			return 1;
		}
		
		else //we must create the enrichement, it can be in 2 or 3 parts. we'll start with 3 always.
		{
			//to begin with we reset the NEnriched:
			Nenriched=ZeroMatrix(3,3);
			
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
			
		}
		
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
			for (unsigned int j=0;j<2;j++)	//x,y coordinates
				aux_points(i,j)= rPoints(edge_begin_node,j) * aux_nodes_relative_locations(i) + rPoints(edge_end_node,j) * (1.0- aux_nodes_relative_locations(i));
		}

		
		//we reset all data:
		rGradientsValue[0]=ZeroMatrix(3,2);
		rGradientsValue[1]=ZeroMatrix(3,2);
		rGradientsValue[2]=ZeroMatrix(3,2);
		Nenriched=ZeroMatrix(3,3);
		rGPShapeFunctionValues=ZeroMatrix(3,3);
		

		 //now we must check the 4 created partitions of the domain.	
		 //one has been collapsed, so we discard it and therefore save only one.
		 unsigned int partition_number=0;		//	
		 //the 3 first partitions are  created using 2 auxiliary nodes and a normal node. at least one of these will be discarded due to zero area		
		 //the last one is composed by the 3 auxiliary nodes. it 'looks' wrong, but since at least one has been collapsed, it actually has a normal node.      
		 for (unsigned int i=0; i<4; i++) //i partition	
		 {	 
			 unsigned int j_aux = i + 2;
			 if (j_aux>2) j_aux -= 3; 
			 BoundedMatrix<int,3,2> partition_father_nodes;
			 array_1d<double,3> N;
			 if (i<4)
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
				 bool is_found = CalculatePosition ( coord_subdomain , x_GP_partition ,y_GP_partition ,z_GP_partition , N);
				 //we check the partition sign.
				 const double partition_sign = (N(0)*rDistances(0) + N(1)*rDistances(1) + N(2)*rDistances(2))/fabs(N(0)*rDistances(0) + N(1)*rDistances(1) + N(2)*rDistances(2));
				 rPartitionsSign(partition_number)=partition_sign;
				 //now we must add the contribution to the normal nodes only if they have the same sign:
				 for (unsigned int j=0;j<3;j++) //j (real) node
				 {
					 if((partition_sign*rDistances(j))>0) //ok. add contribution! 
					 {
						 //now we loop in the nodes that define the partition.
						 for (unsigned int k=0;k<3;k++) //partition	
							 if (partition_father_nodes(k,0)==j || partition_father_nodes(k,1)==j )// (first node)
							 {
								 Nenriched(partition_number,j)+=one_third; //partition, shape function
								 rGradientsValue[partition_number](j,0)+=DN_DX_subdomain(k,0); //[i_partition], (shape function gradient,direction(x,y))
								 rGradientsValue[partition_number](j,1)+=DN_DX_subdomain(k,1); //[i_partition], (shape function gradient,direction(x,y))
							 }
 
					 }
					 //else //do nothing. it simply can't add to a  node that is not in the same side, since we are creating discontinous shape functions
				 }
				 
				 rGPShapeFunctionValues(partition_number,0)=N(0);
				 rGPShapeFunctionValues(partition_number,1)=N(1);
				 rGPShapeFunctionValues(partition_number,2)=N(2);
				 
				 partition_number++;
				 
			 }
			 
		 }

		 
		 

		 return 3;
		 KRATOS_CATCH("");
        
    }
    
    //with some added vectors, to be used when we need information about the interfase between the 2 elements
    template<class TMatrixType, class TVectorType, class TGradientType>
    static int CalculateTriangleDiscontinuousShapeFunctions_ZeroInBoundary(TMatrixType const& rPoints, TGradientType const& DN_DX,
            TVectorType rDistances, TVectorType& rVolumes, TMatrixType& rGPShapeFunctionValues,
            TVectorType& rPartitionsSign, std::vector<TMatrixType>& rGradientsValue, TMatrixType& NEnriched, //and information about the interfase:
            TVectorType& face_gauss_N, TVectorType& face_gauss_Nenriched, double& face_Area, TVectorType& face_n ,unsigned int& type_of_cut)    
    {
        KRATOS_TRY
	
		//unsigned int i,j,k;
		unsigned int i_aux,j_aux,k_aux; //
		type_of_cut = 0;   // 0 means no cuts, 1 means element is cut through edges ij,ik;    2 ij,jk ;    3 ik , kj ;   INTERFASES ON nodes are not contemplated   
		const double one_third=1.0/3.0;
		BoundedMatrix<double, 3, 2 > coord_subdomain; //used to pass arguments when we must calculate areas, shape functions, etc
		BoundedMatrix<double,3,2> DN_DX_subdomain; //used to retrieve derivatives
		double Area;//area of the complete element
		rGPShapeFunctionValues(0,0)=one_third; rGPShapeFunctionValues(0,1)=one_third; rGPShapeFunctionValues(0,2)=one_third; //default, when no interfase has been found
		Area = CalculateVolume2D( rPoints );


        //to begin with we must check whether our element is cut or not by the interfase.
        if( (rDistances(0)*rDistances(1))>0.0 && (rDistances(0)*rDistances(2))>0.0 ) //it means that this element IS NOT cut by the interfase. we must return data of a normal, non-enriched element
		{
			rVolumes(0)=Area;
			rGPShapeFunctionValues(0,0)=one_third; rGPShapeFunctionValues(0,1)=one_third; rGPShapeFunctionValues(0,2)=one_third;
			NEnriched(0,0) = 0.0;
			type_of_cut=1;
            for (int j = 0; j < 3; j++)
                rGradientsValue[0](0, j) = 0.0;
            if (rDistances(0) < 0.0) rPartitionsSign[0] = -1.0;
            else rPartitionsSign[0] = 1.0;
			//KRATOS_WATCH("one element not in the intefase")
			return 1;
		}
		
		else //we must create the enrichement, it can be in 2 or 3 parts. we'll start with 3 always.
		{
			//to begin with we reset the NEnriched:
			NEnriched=ZeroMatrix(3,2);
			
			//KRATOS_WATCH("one element IS in the intefase")
			if ((rDistances(0)*rDistances(1))<0.0) //edge 12 is cut
			{
				if ((rDistances(0)*rDistances(2))<0.0) //edge 13 is cut. 
				{
					//check distances later to see if they're close to nodes and belong to case 4,5or6. for the moment we suppose all edges are cut (cases 1-3)
					type_of_cut=1;
				}
				else
				{
					type_of_cut=2;
				}
			}
			else //edge ij is not cut
			{
				type_of_cut=3;
			}
		}
		
		//KRATOS_WATCH(type_of_cut) 
		switch (type_of_cut)
		{
			case 1: //  edge 12 and 13 is cut
				i_aux=0;
				j_aux=1;
				k_aux=2;	
				
				break;
			case 2:  //  edge 12 and 23 is cut
				i_aux=1;
				j_aux=2;
				k_aux=0;	
				break;
			case 3:  //  edge 23 and 13 is cut
				i_aux=2;
				j_aux=0;
				k_aux=1;					
				break;
		}
		/*
		KRATOS_WATCH(i_aux)
		KRATOS_WATCH(j_aux)
		KRATOS_WATCH(k_aux)
		*/
		 //const double dist12=abs(rDistances(0)-rDistances(1) );
		 if (rDistances(i_aux) < 0.0) 
		 {
			 rPartitionsSign[0] = -1.0;
			 rPartitionsSign[1] =  1.0;
			 rPartitionsSign[2] =  1.0;
		 }
		 else
		 {
			 rPartitionsSign[0] =  1.0;
			 rPartitionsSign[1] = -1.0;
			 rPartitionsSign[2] = -1.0;
		 }
		 
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
		const double tolerable_distance =longest_distance*0.01;	// (1/1,000,000 seems to have good results)
			
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
		 
		 const double node4_relative_position=fabs(rDistances(j_aux)/(rDistances(i_aux)-rDistances(j_aux) ) ) ; //position in 'natural' coordinates of edge 12, 0 when it passes over node 2. (it is over the edge 12)
		 const double node5_relative_position=fabs(rDistances(k_aux)/(rDistances(i_aux)-rDistances(k_aux) ) ) ; //position in 'natural' coordinates of edge 12, 0 when it passes over node 2. (it is over the edge 23)
		 //KRATOS_WATCH(node4_relative_position);
		 //KRATOS_WATCH(node5_relative_position);
		 
		 //Standard Shape function values in the 'new nodes' created by the interfase:
		 const double Ni_aux_node4 = node4_relative_position;
		 const double Nj_aux_node4 = (1.0-node4_relative_position);
		 //Nk_aux_node_4 is zero (we are on the ij edge)
		 const double Ni_aux_node5 = node5_relative_position;
		 //Nj_aux_node_5 is zero (we are on the ik edge)
		 const double Nk_aux_node5 = (1.0-node5_relative_position);
		 
		 //location of the midpoint of the interface : halfway between the two points
		 const double x_midpoint=((rPoints(i_aux,0)*Ni_aux_node4+rPoints(j_aux,0)*Nj_aux_node4)+
							(rPoints(i_aux,0)*Ni_aux_node5+rPoints(k_aux,0)*Nk_aux_node5))*0.5;
		 const double y_midpoint=((rPoints(i_aux,1)*Ni_aux_node4+rPoints(j_aux,1)*Nj_aux_node4)+
							(rPoints(i_aux,1)*Ni_aux_node5+rPoints(k_aux,1)*Nk_aux_node5))*0.5;
		 const double z_midpoint=0.0;
		
		//now we look for the shape function values in the midpoint:
		 coord_subdomain = rPoints;			
		 array_1d<double, 3 >  N;		
		 bool is_found = CalculatePosition ( coord_subdomain , x_midpoint,y_midpoint,z_midpoint, N);
		 face_gauss_N = N;
		 
		  //we will adimensionalize the new shape functions so that the value is one on the midpoint on the interfase:
		  const double adim_Nenriched_i_aux = 1.0 /(face_gauss_N(i_aux)); //for partitions 2 and 3
		  //also the other 2 shape fuctions must be resized so that the new enrichment is constant along the interfase:
		  const double adim_Nenriched_j_aux = adim_Nenriched_i_aux * Ni_aux_node4 / Nj_aux_node4 ; //to have a constant value in the interfase at node 4: =N_j_aux_node4 / N_j_aux when adimensionalized
		  const double adim_Nenriched_k_aux =  adim_Nenriched_i_aux * Ni_aux_node5 / Nk_aux_node5 ; //to have a constant value in the interfase at node 5: = N_k_aux_node4 / N_k_aux when adimensionalized
		  
		  //for the jump, we will create a shape function that holds a constant difference of 2 along the interfase: 
		  const double adim_Nenriched_j_aux_b = (2.0 - node4_relative_position * adim_Nenriched_i_aux) / Nj_aux_node4;
		  const double adim_Nenriched_k_aux_b = (2.0 - node5_relative_position * adim_Nenriched_i_aux) / Nk_aux_node5;
		  
		 
		  //value of the shape functions in the interfase
		  face_gauss_Nenriched(i_aux)= 1.0;
		  face_gauss_Nenriched(j_aux)= 0.0; //actualy it's not defined, it's a discontinous function
		 
		 //first partition
		 //now we must calculate the position of the new nodes to get the area.
		 coord_subdomain(0,0)=rPoints(i_aux,0);
		 coord_subdomain(0,1)=rPoints(i_aux,1);
		 coord_subdomain(1,0)=rPoints(i_aux,0)*Ni_aux_node4+rPoints(j_aux,0)*Nj_aux_node4;
		 coord_subdomain(1,1)=rPoints(i_aux,1)*Ni_aux_node4+rPoints(j_aux,1)*Nj_aux_node4;
		 coord_subdomain(2,0)=rPoints(i_aux,0)*Ni_aux_node5+rPoints(k_aux,0)*Nk_aux_node5;
		 coord_subdomain(2,1)=rPoints(i_aux,1)*Ni_aux_node5+rPoints(k_aux,1)*Nk_aux_node5;
		 //rVolumes(0)=CalculateVolume2D(coord_subdomain);
		 
		 CalculateGeometryData(coord_subdomain, DN_DX_subdomain, rVolumes(0));
		 
		 rGradientsValue[0]=ZeroMatrix(3,2);
		 rGradientsValue[0](i_aux,0)=DN_DX_subdomain(0,0);
		 rGradientsValue[0](i_aux,1)=DN_DX_subdomain(0,1);		 
		 //all the others are zero!!
		 
		 //we look for the gauss point of the partition:
		 double x_GP_partition =  one_third * ( coord_subdomain(0,0) + coord_subdomain(1,0) + coord_subdomain(2,0) );
		 double y_GP_partition =  one_third * ( coord_subdomain(0,1) + coord_subdomain(1,1) + coord_subdomain(2,1) );
		 double z_GP_partition  = 0.0;
		 //we reset the coord_subdomain matrix so that we have the whole element again:
		 coord_subdomain = rPoints;	
		 //and we calculate its shape function values
		 is_found = CalculatePosition ( coord_subdomain , x_GP_partition ,y_GP_partition ,z_GP_partition , N);
		 rGPShapeFunctionValues(0,0)=N(0);
		 rGPShapeFunctionValues(0,1)=N(1);
		 rGPShapeFunctionValues(0,2)=N(2);
		 
		 NEnriched(0,i_aux)=one_third;
		 //the others are zero
		 //NEnriched(0,j_aux)=0.0;
		 //NEnriched(0,k_aux)=0.0;
		 
		 
		 //KRATOS_WATCH(rVolumes(0));
		 
		 //now the face area(actually it's just the distance from point 4 to 5.
		 face_Area=sqrt(pow((coord_subdomain(2,0)-coord_subdomain(1,0)),2)+pow((coord_subdomain(2,1)-coord_subdomain(1,1)),2));
		 //and the normal vector. face_Area already has the modulus of the vector, so:
		 face_n(0)=(coord_subdomain(2,1)-coord_subdomain(1,1))/face_Area;
		 face_n(1)=-(coord_subdomain(2,0)-coord_subdomain(1,0))/face_Area;
		 
		 /*
		 KRATOS_WATCH(coord_subdomain(0,0));
		 KRATOS_WATCH(coord_subdomain(0,1));
		 KRATOS_WATCH(coord_subdomain(1,0));
		 KRATOS_WATCH(coord_subdomain(1,1));
		 KRATOS_WATCH(coord_subdomain(2,0));
		 KRATOS_WATCH(coord_subdomain(2,1));
		 KRATOS_WATCH(rGPShapeFunctionValues(0,i_aux));
		 KRATOS_WATCH(rGPShapeFunctionValues(0,j_aux));
		 KRATOS_WATCH(rGPShapeFunctionValues(0,k_aux));
		 KRATOS_WATCH( rGradientsValue[0](0,0))
		 KRATOS_WATCH( rGradientsValue[0](0,1))
		*/
		
		
		
		//second partition and second GP 
		 
 
		 //now we must calculate the position of the new nodes to get the area.
		 //coord_subdomain = rPoints; //easier to start this way. node 2 is already ok.
		 coord_subdomain(0,0) = rPoints(j_aux,0);
		 coord_subdomain(0,1) = rPoints(j_aux,1);
		 coord_subdomain(1,0) = rPoints(k_aux,0);
		 coord_subdomain(1,1) = rPoints(k_aux,1);
		 coord_subdomain(2,0) = rPoints(i_aux,0)*Ni_aux_node5+rPoints(k_aux,0)*Nk_aux_node5;
		 coord_subdomain(2,1) = rPoints(i_aux,1)*Ni_aux_node5+rPoints(k_aux,1)*Nk_aux_node5;
		 
		 //rVolumes(1)=CalculateVolume2D(coord_subdomain);
		 CalculateGeometryData(coord_subdomain, DN_DX_subdomain, rVolumes(1));
		 
		 //for the first Shape Funct(i_aux) they art zero
		 rGradientsValue[1]=ZeroMatrix(3,2);
		 rGradientsValue[1](j_aux,0)=DN_DX_subdomain(0,0);
		 rGradientsValue[1](j_aux,1)=DN_DX_subdomain(0,1);
		 rGradientsValue[1](k_aux,0)= DN_DX_subdomain(1,0);
		 rGradientsValue[1](k_aux,1)= DN_DX_subdomain(1,1);
		 
		 
		 
		 //we look for the gauss point of the partition:
		 x_GP_partition =  one_third * ( coord_subdomain(0,0) + coord_subdomain(1,0) + coord_subdomain(2,0) );
		 y_GP_partition =  one_third * ( coord_subdomain(0,1) + coord_subdomain(1,1) + coord_subdomain(2,1) );
		 z_GP_partition  = 0.0;
		 //we reset the coord_subdomain matrix so that we have the whole element again:
		 coord_subdomain = rPoints;
		 //and we calculate the shape functions at it	
		 is_found = CalculatePosition ( coord_subdomain , x_GP_partition ,y_GP_partition ,z_GP_partition , N);
		 rGPShapeFunctionValues(1,0)=N(0);
		 rGPShapeFunctionValues(1,1)=N(1);
		 rGPShapeFunctionValues(1,2)=N(2);
		 //now the values of the shape functions. both are zero except the first one
		 //NEnriched(1,i_aux) = 0.0;
		 NEnriched(1,j_aux) = one_third;   
		 NEnriched(1,k_aux) = one_third;   
		 /*
		 KRATOS_WATCH(rVolumes(1)); 
		 		 KRATOS_WATCH(coord_subdomain(0,0));
		 KRATOS_WATCH(coord_subdomain(0,1));
		 KRATOS_WATCH(coord_subdomain(1,0));
		 KRATOS_WATCH(coord_subdomain(1,1));
		 KRATOS_WATCH(coord_subdomain(2,0));
		 KRATOS_WATCH(coord_subdomain(2,1));
		 KRATOS_WATCH(rGPShapeFunctionValues(1,i_aux));
		 KRATOS_WATCH(rGPShapeFunctionValues(1,j_aux));
		 KRATOS_WATCH(rGPShapeFunctionValues(1,k_aux));
		 KRATOS_WATCH( rGradientsValue[1](0,0))
		 KRATOS_WATCH( rGradientsValue[1](0,1))
		 */
		  
		 //NOT A PARTITION, just the volume of a different conectivity:
		 //just replacing the third node. Useful only to recover information, for example N(i,j,k) of node5
		 coord_subdomain(2,0) = rPoints(i_aux,0)*Ni_aux_node4+rPoints(j_aux,0)*Nj_aux_node4;
		 coord_subdomain(2,1) = rPoints(i_aux,1)*Ni_aux_node4+rPoints(j_aux,1)*Nj_aux_node4;
		 rVolumes(3)=CalculateVolume2D(coord_subdomain);

		  //third partition:

		 coord_subdomain(0,0) = rPoints(j_aux,0);
		 coord_subdomain(0,1) = rPoints(j_aux,1);
		 coord_subdomain(1,0)=rPoints(i_aux,0)*Ni_aux_node5+rPoints(k_aux,0)*Nk_aux_node5;
		 coord_subdomain(1,1)=rPoints(i_aux,1)*Ni_aux_node5+rPoints(k_aux,1)*Nk_aux_node5;
		 coord_subdomain(2,0)=rPoints(i_aux,0)*Ni_aux_node4+rPoints(j_aux,0)*Nj_aux_node4;
		 coord_subdomain(2,1)=rPoints(i_aux,1)*Ni_aux_node4+rPoints(j_aux,1)*Nj_aux_node4;
		 //rVolumes(2)=Area-rVolumes(0)-rVolumes(1);
		 CalculateGeometryData(coord_subdomain, DN_DX_subdomain, rVolumes(2));
		  //for the first Shape Funct(i_aux) and second(j_aux) they art zero
		 rGradientsValue[2]=ZeroMatrix(3,2);
		 rGradientsValue[2](j_aux,0)= DN_DX_subdomain(0,0);
		 rGradientsValue[2](j_aux,1)= DN_DX_subdomain(0,1);
		 //now the values of the shape functions. both are zero except the second one
		 //NEnriched(1,i_aux) = 0.0;
		 //NEnriched(1,j_aux) = one_third;   
		 NEnriched(2,j_aux) = one_third;   
		 
		 
		 x_GP_partition =  one_third * ( coord_subdomain(0,0) + coord_subdomain(1,0) + coord_subdomain(2,0) );
		 y_GP_partition =  one_third * ( coord_subdomain(0,1) + coord_subdomain(1,1) + coord_subdomain(2,1) );
		 z_GP_partition  = 0.0;
		 coord_subdomain = rPoints;	
		 is_found = CalculatePosition ( coord_subdomain , x_GP_partition ,y_GP_partition ,z_GP_partition , N);
		 rGPShapeFunctionValues(2,0)=N(0);
		 rGPShapeFunctionValues(2,1)=N(1);
		 rGPShapeFunctionValues(2,2)=N(2);
		 
		 /*
			KRATOS_WATCH(rVolumes(2)); 
		 		 KRATOS_WATCH(coord_subdomain(0,0));
		 KRATOS_WATCH(coord_subdomain(0,1));
		 KRATOS_WATCH(coord_subdomain(1,0));
		 KRATOS_WATCH(coord_subdomain(1,1));
		 KRATOS_WATCH(coord_subdomain(2,0));
		 KRATOS_WATCH(coord_subdomain(2,1));
		 KRATOS_WATCH(rGPShapeFunctionValues(2,i_aux));
		 KRATOS_WATCH(rGPShapeFunctionValues(2,j_aux));
		 KRATOS_WATCH(rGPShapeFunctionValues(2,k_aux));
		 KRATOS_WATCH( rGradientsValue[2](0,0))
		 KRATOS_WATCH( rGradientsValue[2](0,1))
		 */ 
		 /*
		 std::cout <<"GAUSS POINTS" << '\n';
		 for (unsigned int ji=0; ji<3; ji++)
			std::cout <<rGPShapeFunctionValues(ji,0)<< "  " << rGPShapeFunctionValues(ji,1) << "  " << rGPShapeFunctionValues(ji,2) <<'\n';
			
		 std::cout <<"GRADIENTS" << '\n'; 
		 for (unsigned int ji=0; ji<3; ji++)
		 {
			std::cout << "first Shape function" << '\n';
			std::cout << rGradientsValue[ji](0,0) << "  " <<rGradientsValue[ji](0,1) <<'\n';
			std::cout << "second Shape function" << '\n';
			std::cout << rGradientsValue[ji](1,0) << "  "  <<rGradientsValue[ji](1,1) <<'\n';
		 }
		 */
		 return 3;
		 KRATOS_CATCH("");
        
    }
    
    
    
    
    
    /**this function is designed to compute the shape function derivatives, shape functions and volume in 3D
     * @param geom it is the array of nodes. It is expected to be a triangle
     * @param a stack matrix of size 3*2 to store the shape function's derivatives
     * @param an array_1d to store the shape functions at the barycenter
     * @param the volume of the element
     */
    private:
    
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
            } else
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
        

    
    };

} // namespace Kratos.

#endif // KRATOS_DISCONTINUOUS_2D_UTILITIES_INCLUDED  defined


