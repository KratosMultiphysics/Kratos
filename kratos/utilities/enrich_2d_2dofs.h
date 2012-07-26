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
//   Last Modified by:    $Author: pablo $
//   Date:                $Date: 2009-01-13 16:40:58 $
//   Revision:            $Revision: 1.24 $
//
//


#if !defined(KRATOS_ENRICHMENT_2D_UTILITIES_INCLUDED )
#define  KRATOS_ENRICHMENT_2D_UTILITIES_INCLUDED



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

/** This utility can be used to calculate the enriched shape function for tetrahedra element.
 *  The metodology consists in partitioning the tetrahedra in a set of sub-tetrahedra and
 *  cacluate the enrichment information using these partitions.
 */
class EnrichmentUtilities_2D
{
public:

    /**
     * The method to calculate the ernriched shape functions for given tetrahedra.
     * @param rPoints A 3x3 matrix where row i has the coordinates of node i.
     * @param DN_DX The gradient of the shape functions Ni respect to the reference coordinates
     * @param rDistances is an input  vector of 3 size which holds relative distance for each node.
     *        it is used internally to mark the position of the zero level
     * @param rVolumes Result vector with size 3 (maximumn number of partitions) holding the volume of each partition ?                  O DEBE SER 2????????????????????????????????
     * @param rShapeFunctionValues Result 3x3 matrix where each row represents a partition and holds the shape functions N1 to N3 ( the cut)
     *        of the original triangle evaluated in the gauss point (center) of the partition.
     *        so that it is  N(gauss_index, node_index)
     * @param rPartitionsSign A result vector of 6 holding the sign of the distance for the partition.
     *        The value -1 represents the negative distance sign, 1 represents positive distance and 0 stands for not used partition
     * @param rGradientsValue Restult vector of size 3 holding the gradient of the enriched shape funciton for each volume.
     *        Each element of vector is a 1x3 matrix representing the gradient of enriched shape function. The use of
     *        matrix is for possible future improvement.
     * @param Nenriched is a Matrix that contains for every gauss point the values of the enriched shape functions at the position of the gauss point
     *        so that Nenriched(1,0) contains the value of the enriched shape function "0" at the gauss point "1"
     * @return number of partitions created which can be from 1 to 3.
     *         1 holds for only 1 partition which is the original element. (No partitioning needed)
     */
   
    
    //with some added vectors, to be used when we need information about the interfase between the 2 elements
    template<class TMatrixType, class TVectorType, class TGradientType>
    static int CalculateTriangleEnrichedShapeFuncions(TMatrixType const& rPoints, TGradientType const& DN_DX,
            TVectorType rDistances, TVectorType& rVolumes, TMatrixType& rGPShapeFunctionValues,
            TVectorType& rPartitionsSign, std::vector<TMatrixType>& rGradientsValue, TMatrixType& NEnriched, //and information about the interfase:
            TVectorType& face_gauss_N, TVectorType& face_gauss_Nenriched, double& face_Area, TVectorType& face_n ,unsigned int& type_of_cut)    
    {
        KRATOS_TRY

		unsigned int i,j,k;
		unsigned int i_aux,j_aux,k_aux;
 //       const int n_nodes = 3; // it works only for triangles
 //       const int n_edges = 3; // it works only for triangles
		type_of_cut = 0;   // 0 means no cuts, 1 means element is cut through edges ij,ik;    2 ij,jk ;    3 ik , kj ;     4 only ij ; 5 only jk ;  6 only ik   
							   //                  (somehow they look like the quadratic shape functions for a triangle
		const double one_third=1.0/3.0;
		bounded_matrix<double, 3, 3 > coord_subdomain; //used to pass arguments when we must calculate areas, shape functions, etc
		//array_1d<double,3>& GP1, GP2, GP3;
		double Area;//, Area1, Area2, Area3;
        //boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX_element;
        //array_1d<double,3>& N_element;
		rGPShapeFunctionValues(0,0)=one_third; rGPShapeFunctionValues(0,1)=one_third; rGPShapeFunctionValues(0,2)=one_third;
		Area = CalculateVolume2D( rPoints );

        //const int edge_i[] = {0, 0, 0, 1, 1, 2};
        //const int edge_j[] = {1, 2, 3, 2, 3, 3};    
        
        //to begin with we must check whether our element is cut or not by the interfase.
        if( (rDistances(0)*rDistances(1))>0.0 && (rDistances(0)*rDistances(2))>0.0 ) //it means that this element IS NOT cut by the interfase. we must return data of a normal, non-enriched element
		{
			//CalculateGeometryData(rPoints, DN_DX, rShapeFunctionValues,rVolumes);   //mmmm revisar!
			rVolumes(0)=Area;
			rGPShapeFunctionValues(0,0)=one_third; rGPShapeFunctionValues(0,1)=one_third; rGPShapeFunctionValues(0,2)=one_third;
			NEnriched(0,0) = 0.0;
            for (int j = 0; j < 3; j++)
                rGradientsValue[0](0, j) = 0.0;
            if (rDistances(0) < 0.0) rPartitionsSign[0] = -1;
            else rPartitionsSign[0] = 1.0;
			//KRATOS_WATCH("one element not in the intefase")
			return 1;
		}
		
		else //we must create the enrichement, it can be in 2 or 3 parts. we'll start with 3 always.
		{
			KRATOS_WATCH("one element IS in the intefase")
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
		
		
		 //const double dist12=abs(rDistances(0)-rDistances(1) );
		 if (rDistances(i_aux) < 0.0) 
		 {
			 rPartitionsSign[0] = -1;
			 rPartitionsSign[1] =  1;
			 rPartitionsSign[2] =  1;
		 }
		 else
		 {
			 rPartitionsSign[0] =  1;
			 rPartitionsSign[1] = -1;
			 rPartitionsSign[2] = -1;
		 }
		 
		 //for (int jj = 0; jj < 3; jj++)
		 //	KRATOS_WATCH(rDistances(jj));
		 
		 const double node4_relative_position=fabs(rDistances(j_aux)/(rDistances(i_aux)-rDistances(j_aux) ) ) ; //position in 'natural' coordinates of edge 12, 0 when it passes over node 2. (it is over the edge 12)
		 const double node5_relative_position=fabs(rDistances(k_aux)/(rDistances(i_aux)-rDistances(k_aux) ) ) ; //position in 'natural' coordinates of edge 12, 0 when it passes over node 2. (it is over the edge 23)
		 KRATOS_WATCH(node4_relative_position);
		 KRATOS_WATCH(node5_relative_position);
		 
		 
		 		 //location of the midpoint of the interface
		 face_gauss_N(i_aux)= 0.5*(node4_relative_position+node5_relative_position);
		 face_gauss_N(j_aux)= 0.5*( 1.0 - node4_relative_position);
		 face_gauss_N(k_aux)= 0.5*( 1.0 - node5_relative_position);
		 
		 
		 //we should define them like:
		 
		 ///face_gauss_Nenriched(i_aux)= face_gauss_N(j_aux)+face_gauss_N(k_aux);
		 ///face_gauss_Nenriched(j_aux)= face_gauss_N(i_aux);
		 //but to be able to impose strongly the jump the must be one. instead we'll use the previous parameters to scale the shape functions so that they're one in the midpoint of the interfase
		 //this way the jump can be imposed (strongly) directly.
		  const double adim_Nenriched_i_aux = 1.0 /(face_gauss_N(j_aux)+face_gauss_N(k_aux));
		  const double adim_Nenriched_j_aux = 1.0 /(face_gauss_N(i_aux));
		  const double adim_Nenriched_node4 = node4_relative_position * adim_Nenriched_j_aux / (1.0 - node4_relative_position ) ;
		  const double adim_Nenriched_node5 = node5_relative_position * adim_Nenriched_j_aux / (1.0 - node5_relative_position ) ;
		  face_gauss_Nenriched(i_aux)= 1.0;
		  face_gauss_Nenriched(j_aux)= 0.0; //actualy it's not defined, it's a discontinous function
		 
		 //first partition
		 rGPShapeFunctionValues(0,i_aux)=one_third*(1.0 + node4_relative_position + node5_relative_position);  //we create 3 gauss points, 
		 rGPShapeFunctionValues(0,j_aux)=one_third*(1.0 - node4_relative_position);  //the triangle from gauss point 1 is indepentent (has its own plane), the other two have the same shape function. 
		 rGPShapeFunctionValues(0,k_aux)=one_third*(1.0 - node5_relative_position);  //the triangle from gauss point 1 is indepentent (has its own plane), the other two have the same shape function. 
		 //rGradientsValue[0](0,0)=(DN_DX(j_aux,0)+DN_DX(k_aux,0))*adim_Nenriched_i_aux;
		 //rGradientsValue[0](0,1)=(DN_DX(j_aux,1)+DN_DX(k_aux,1))*adim_Nenriched_i_aux;		  //	      i   j,k				 i    j,k
		 rGradientsValue[0](0,0)=DN_DX(j_aux,0)*adim_Nenriched_node4+DN_DX(k_aux,0)*adim_Nenriched_node5;
		 rGradientsValue[0](0,1)=DN_DX(j_aux,1)*adim_Nenriched_node4+DN_DX(k_aux,1)*adim_Nenriched_node5;		  //	      i   j,k				 i    j,k
		 
		 rGradientsValue[0](1,0)=(DN_DX(j_aux,0)+DN_DX(k_aux,0))*adim_Nenriched_i_aux;                 //the shape function are: 1:   ___/\____       2:	    ___/ ___ 
		 rGradientsValue[0](1,1)=(DN_DX(j_aux,1)+DN_DX(k_aux,1))*adim_Nenriched_i_aux;			    	//															/
		 NEnriched(0,0)=rGPShapeFunctionValues(0,j_aux)*adim_Nenriched_node4+rGPShapeFunctionValues(0,k_aux)*adim_Nenriched_node5;
		 NEnriched(0,1)=(rGPShapeFunctionValues(0,j_aux)+rGPShapeFunctionValues(0,k_aux))*adim_Nenriched_i_aux;;  
		 //now we must calculate the position of the new nodes to get the area.
		 //coord_subdomain=rPoints; //easier to start this way. node 1 is already ok.
		 coord_subdomain(0,0)=rPoints(i_aux,0);
		 coord_subdomain(0,1)=rPoints(i_aux,1);
		 coord_subdomain(1,0)=rPoints(i_aux,0)*(node4_relative_position)+rPoints(j_aux,0)*(1.0-node4_relative_position);
		 coord_subdomain(1,1)=rPoints(i_aux,1)*(node4_relative_position)+rPoints(j_aux,1)*(1.0-node4_relative_position);
		 coord_subdomain(2,0)=rPoints(i_aux,0)*(node5_relative_position)+rPoints(k_aux,0)*(1.0-node5_relative_position);
		 coord_subdomain(2,1)=rPoints(i_aux,1)*(node5_relative_position)+rPoints(k_aux,1)*(1.0-node5_relative_position);
		 rVolumes(0)=CalculateVolume2D(coord_subdomain);
		 KRATOS_WATCH(rVolumes(0));
		 
		 //now the face area(actually it's just the distance from point 4 to 5.
		 face_Area=sqrt(pow((coord_subdomain(2,0)-coord_subdomain(1,0)),2)+pow((coord_subdomain(2,1)-coord_subdomain(1,1)),2));
		 //and the normal vector. face_Area already has the modulus of the vector, so:
		 face_n(0)=-(coord_subdomain(2,0)-coord_subdomain(1,0))/face_Area;
		 face_n(1)=(coord_subdomain(2,1)-coord_subdomain(1,1))/face_Area;
		 
		 /*
		 KRATOS_WATCH(coord_subdomain(0,0));
		 KRATOS_WATCH(coord_subdomain(0,1));
		 KRATOS_WATCH(coord_subdomain(1,0));
		 KRATOS_WATCH(coord_subdomain(1,1));
		 KRATOS_WATCH(coord_subdomain(2,0));
		 KRATOS_WATCH(coord_subdomain(2,1));
		*/
		
		//second partition and second GP
		 rGPShapeFunctionValues(1,i_aux)=one_third*(node5_relative_position);  
		 rGPShapeFunctionValues(1,j_aux)=one_third;  
		 rGPShapeFunctionValues(1,k_aux)=one_third*(2.0 - node5_relative_position);
		 rGradientsValue[1](0,0)=DN_DX(i_aux,0)*adim_Nenriched_j_aux;
		 rGradientsValue[1](0,1)=DN_DX(i_aux,1)*adim_Nenriched_j_aux;
		 
		 rGradientsValue[1](1,0)= -rGradientsValue[1](0,0);
		 rGradientsValue[1](1,1)= -rGradientsValue[1](0,1);
		 
		 NEnriched(1,0) = rGPShapeFunctionValues(1,0)*adim_Nenriched_j_aux;
		 NEnriched(1,1) = -NEnriched(1,0);    
		 //now we must calculate the position of the new nodes to get the area.
		 //coord_subdomain = rPoints; //easier to start this way. node 2 is already ok.
		 coord_subdomain(0,0) = rPoints(j_aux,0);
		 coord_subdomain(0,1) = rPoints(j_aux,1);
		 coord_subdomain(1,0) = rPoints(k_aux,0);
		 coord_subdomain(1,1) = rPoints(k_aux,1);
		 coord_subdomain(2,0) = rPoints(i_aux,0)*(node5_relative_position)+rPoints(j_aux,0)*(1.0-node5_relative_position);
		 coord_subdomain(2,1) = rPoints(i_aux,1)*(node5_relative_position)+rPoints(j_aux,1)*(1.0-node5_relative_position);
		 rVolumes(1)=CalculateVolume2D(coord_subdomain);
		 /*
		 KRATOS_WATCH(coord_subdomain(0,0));
		 KRATOS_WATCH(coord_subdomain(0,1));
		 KRATOS_WATCH(coord_subdomain(1,0));
		 KRATOS_WATCH(coord_subdomain(1,1));
		 KRATOS_WATCH(coord_subdomain(2,0));
		 KRATOS_WATCH(coord_subdomain(2,1));
		 */
		 KRATOS_WATCH(rVolumes(1));
		 //third partition and third GP
		 rGPShapeFunctionValues(2,i_aux)=one_third*(node4_relative_position + node5_relative_position);  
		 rGPShapeFunctionValues(2,j_aux)=one_third*(2.0 - node4_relative_position);  
		 rGPShapeFunctionValues(2,k_aux)=one_third*(1.0 - node5_relative_position);  
		 rGradientsValue[2](0,0)=DN_DX(i_aux,0)*adim_Nenriched_j_aux;
		 rGradientsValue[2](0,1)=DN_DX(i_aux,1)*adim_Nenriched_j_aux; 	
		 
		 rGradientsValue[2](1,0)= -rGradientsValue[1](0,0);
		 rGradientsValue[2](1,1)= -rGradientsValue[1](0,1);
		 
		 NEnriched(2,0)= rGPShapeFunctionValues(1,0)*adim_Nenriched_j_aux;
		 NEnriched(2,1)= -NEnriched(2,0);
		 //now we must calculate the position of the new nodes to get the area.
		 /*
		 coord_subdomain=rPoints; //easier to start this way. node 2 and 3 is already ok.
		 coord_subdomain(0,0)=rPoints(0,0)*(node5_relative_postion)+rPoints(2,0)*(1.0-node5_relative_postion);
		 coord_subdomain(0,1)=rPoints(0,0)*(node5_relative_postion)+rPoints(2,2)*(1.0-node5_relative_postion);
		 Area1=CalculateVolume2D(coord_subdomain);
		 */
		 

		 
		 
		 
		 
		 
		 rVolumes(2)=Area-rVolumes(0)-rVolumes(1);
		 KRATOS_WATCH(rVolumes(2));
		 
		 
		 
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
			const bounded_matrix<double, 3, 3 > & coordinates,
			boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,
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
			const bounded_matrix<double, 3, 3 > & coordinates)
		{
			double x10 = coordinates(1,0) - coordinates(0,0);
			double y10 = coordinates(1,1) - coordinates(0,1);

			double x20 = coordinates(2,0) - coordinates(0,0);
			double y20 = coordinates(2,1) - coordinates (0,1);
			double detJ = x10 * y20-y10 * x20;
			return 0.5*detJ;
		}
    
    };

} // namespace Kratos.

#endif // KRATOS_ENRICHMENT_UTILITIES_INCLUDED  defined


