//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Julio Marti.
//


#if !defined(KRATOS_ENRICHMENTUTILITIES_INCLUDED )
#define  KRATOS_ENRICHMENTUTILITIES_INCLUDED



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
  class EnrichmentUtilitiesforPFEM2
  {
  public:

    //2D with information on the interfase
    static int CalculateEnrichedShapeFuncions(BoundedMatrix<double,(2+1), 2 >& rPoints, BoundedMatrix<double, (2+1), 2 >& DN_DX,array_1d<double,(2+1)>& rDistances, array_1d<double,(3*(2-1))>& rVolumes, BoundedMatrix<double, 3*(2-1), (2+1) >& rGPShapeFunctionValues,array_1d<double,(3*(2-1))>& rPartitionsSign, std::vector<Matrix>& rGradientsValue, BoundedMatrix<double,3*(2-1), (2)>& NEnriched,array_1d<double,(3)>&  rGPShapeFunctionValues_in_interfase, array_1d<double,(3)>&  NEnriched_in_interfase, double & InterfaseArea)
    {
      KRATOS_TRY

	const double one_third=1.0/3.0;
      BoundedMatrix<double,3,2> aux_points; //for auxiliary nodes 4(between 1 and 2) ,5(between 2 and 3) ,6 (between 3 and 1)
      BoundedMatrix<double, 3, 2 > coord_subdomain; //used to pass arguments when we must calculate areas, shape functions, etc
      BoundedMatrix<double,3,2> DN_DX_subdomain; //used to retrieve derivatives

      double most_common_sign=0; //the side of the cut in which two nodes are found (same sign) will be the ones that remains unchanged when builing the discontinuity
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


    static int CalculateEnrichedShapeFuncionsExtendedmodified(BoundedMatrix<double,(2+1), 2 >& rPoints, BoundedMatrix<double, (2+1), 2 >& DN_DX, array_1d<double,(2+1)>& rDistances, array_1d<double,(3*(2-1))>& rVolumes, BoundedMatrix<double, 3*(2-1), (2+1) >& rGPShapeFunctionValues, array_1d<double,(3*(2-1))>& rPartitionsSign, std::vector<Matrix>& rGradientsValue, BoundedMatrix<double,3*(2-1), (5)>& NEnriched,BoundedMatrix<double,10, 2>& rGradientpositive,BoundedMatrix<double,10, 2>& rGradientnegative ,BoundedMatrix<int,3,3>& father_nodes)
    {
      KRATOS_TRY

	const double one_third=1.0/3.0;
      BoundedMatrix<double,3,2> aux_points; //for auxiliary nodes 4(between 1 and 2) ,5(between 2 and 3) ,6 (between 3 and 1)
      BoundedMatrix<double, 3, 2 > coord_subdomain; //used to pass arguments when we must calculate areas, shape functions, etc
      BoundedMatrix<double,3,2> DN_DX_subdomain; //used to retrieve derivatives

      double most_common_sign=0; //the side of the cut in which two nodes are found (same sign) will be the ones that remains unchanged when builing the discontinuity
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
	  NEnriched(0,0) = 0.0;
	  //type_of_cut=1;
	  for (int j = 0; j < 2; j++)
	    rGradientsValue[0](0, j) = 0.0;
	  if (rDistances(0) < 0.0) rPartitionsSign[0] = -1.0;
	  else rPartitionsSign[0] = 1.0;
	  return 1;
	}

      //else //we must create the enrichement, it can be in 2 or 3 parts. we'll start with 3 always.

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
      //const double tolerable_distance =longest_distance*0.001;// 0.001 	// (1/1,000,000 seems to have good results)
      //and now we check if a distance is too small:

      //END OF TRICK. REMEMBER TO OVERWRITE THE DISTANCE VARIABLE IN THE ELEMENT IN CASE THESE LINES HAVE MODIFIED THEM (distances)


      //const double epsilon = 1e-15; //1.00e-9;
      //compute the gradient of the distance and normalize it
      array_1d<double, 2 > grad_d;
      noalias(grad_d) = prod(trans(DN_DX), rDistances);


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

      //We have 3 edges, meaning we created 3 aux nodes. But one of them actually matches the position of a real node (the one that is not on an interface edge is displaced to one of the ends (a node)
      //the new shape functions are built by setting the values in all (real and aux) nodes to zero except in one of the interphase nodes. in the array aux_node_shape_function_index we assign this value
      array_1d<int, 3 > aux_node_enrichment_shape_function_index; //when not used, it must be -1;

      int shape_function_id=0;
      father_nodes(0,0)=-1;
      father_nodes(0,1)=-1;
      father_nodes(0,2)=-1;
      father_nodes(1,0)=-1;
      father_nodes(1,1)=-1;
      father_nodes(1,2)=-1;
      father_nodes(2,0)=-1;
      father_nodes(2,1)=-1;
      father_nodes(2,2)=-1;

      //KRATOS_WATCH(father_nodes);
      for (unsigned int i=0; i<3; i++) //we go over the 3 edges:
	{
	  int edge_begin_node=i;
	  int edge_end_node=i+1;
	  if (edge_end_node==3) edge_end_node=0; //it's a triangle, so node 3 is actually node 0

	  if(cut_edges(i)==true)
	    {
	      aux_nodes_relative_locations(i)=fabs(rDistances(edge_end_node)/(rDistances(edge_end_node)-rDistances(edge_begin_node) ) ) ; //position in 'natural' coordinates of edge 12, 1 when it passes over node 1. (it is over the edge 01)
	      //KRATOS_WATCH(aux_nodes_relative_locations(i));

	      aux_nodes_father_nodes(i,0)=edge_begin_node;
	      aux_nodes_father_nodes(i,1)=edge_end_node;

	      aux_node_enrichment_shape_function_index(i)=shape_function_id;
	      father_nodes(i,0)=edge_begin_node;
	      father_nodes(i,1)=edge_end_node;
	      father_nodes(i,2)=shape_function_id;
	      shape_function_id++;

	    }
	  else
	    {//<
	      if(fabs(rDistances(edge_end_node))<fabs(rDistances(edge_begin_node))) //if edge is not cut, we collapse the aux node into the node which has the highest absolute value to have "nicer" (less "slivery") subelements
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
	      //KRATOS_WATCH("compileddd");
	      aux_node_enrichment_shape_function_index(i)=-1;
	    }

	  //and we save the coordinate of the new aux nodes:
	  for (unsigned int j=0;j<2;j++){	//x,y coordinates
	    aux_points(i,j)= rPoints(edge_begin_node,j) * aux_nodes_relative_locations(i) + rPoints(edge_end_node,j) * (1.0- aux_nodes_relative_locations(i));
	  }
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
	  double d =    (rPoints(i_node,0) - base_point[0]) * grad_d[0] + (rPoints(i_node,1) - base_point[1]) * grad_d[1] ;
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
      rGradientsValue[0]=ZeroMatrix(5,2);
      rGradientsValue[1]=ZeroMatrix(5,2);
      rGradientsValue[2]=ZeroMatrix(5,2);

      NEnriched=ZeroMatrix(3,5);
      rGPShapeFunctionValues=ZeroMatrix(3,3);


      //now we must check the 4 created partitions of the domain.
      //one has been collapsed, so we discard it and therefore save only one.
      unsigned int partition_number=0;		//
      //the 3 first partitions are  created using 2 auxiliary nodes and a normal node. at least one of these will be discarded due to zero area
      //the last one is composed by the 3 auxiliary nodes. it 'looks' wrong, but since at least one has been collapsed, it actually has a normal node.
      bool found_empty_partition=false;

      //the enrichment is directly the shape functions created by using the partition.
      //we have to save for the enrichment 0 and 1 which is the node that will be active, that is, whose shape function value is zero (for some partitions all 3 shape functions are inactive)



      for (int i=0; i<4; i++) //i partition
	{

	  array_1d<int, 2 > active_node_in_enrichment_shape_function;
	  active_node_in_enrichment_shape_function(0)=-1;  active_node_in_enrichment_shape_function(1)=-1; //initialized as if all are inactive -> gradient=0;
	  //the same but for the replacement shape functions
	  array_1d<int, 3 > active_node_in_replacement_shape_function;
	  active_node_in_replacement_shape_function(0)=-1;  active_node_in_replacement_shape_function(1)=-1; active_node_in_replacement_shape_function(2)=-1; //initialized as if all are inactive -> gradient=0;

	  int j_aux = i + 2;
	  if (j_aux>2) j_aux -= 3;
	  BoundedMatrix<int,3,2> partition_father_nodes;
	  array_1d<double,3> N; //bool useful=false;

	  int useful_node_for_N0star=-1;
	  int useful_node_for_N1star=-1;

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

	      //KRATOS_WATCH(coord_subdomain)
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
	    //  if( (aux_nodes_father_nodes(i,0)!=aux_nodes_father_nodes(i,1)) && (aux_nodes_father_nodes(j_aux,0)!=aux_nodes_father_nodes(j_aux,1)))
		//{
		//  useful=true;
		//}

	      coord_subdomain(0,0)=rPoints(i,0);
	      coord_subdomain(0,1)=rPoints(i,1);
	      coord_subdomain(1,0)=aux_points(i,0);
	      coord_subdomain(1,1)=aux_points(i,1);
	      coord_subdomain(2,0)=aux_points(j_aux,0);
	      coord_subdomain(2,1)=aux_points(j_aux,1);

	      //an aux_node is useful when one of its father nodes is the real node of the subelement. that means that the edge is part of the subelement.
	      if(partition_father_nodes(1,0)==i || partition_father_nodes(1,1)==i) //if one of the father nodes of node_aux_i is equal to the real node i
		{
		  if(aux_node_enrichment_shape_function_index(i)==0)
		    useful_node_for_N0star=1;
		  if(aux_node_enrichment_shape_function_index(i)==1)
		    useful_node_for_N1star=1;
		}
	      if(partition_father_nodes(2,0)==j_aux || partition_father_nodes(2,1)==j_aux) //if one of the father nodes of node_aux_i is equal to the real node i
		{
		  if(aux_node_enrichment_shape_function_index(j_aux)==0)
		    useful_node_for_N0star=2;
		  if(aux_node_enrichment_shape_function_index(j_aux)==1)
		    useful_node_for_N1star=2;
		}
	    }
	  else
	    {
	      partition_father_nodes=aux_nodes_father_nodes;
	      coord_subdomain=aux_points;
	      //useful=true;
	      //int non_aux_node=-1; //one of the nodes of this partition is actually a real node, which has repeated father node
	      int non_aux_node_father_node=-1; //the real node id
	      //we have to check the 3 of them:
	      for (int j = 0; j < 3; j++)
		{
		  if(partition_father_nodes(j,0)==partition_father_nodes(j,1))
		    {
		      //non_aux_node=j;
		      non_aux_node_father_node=partition_father_nodes(j,0);
		    }
		}
	      for (int j = 0; j < 3; j++)
		{
		  if (aux_node_enrichment_shape_function_index(j)> -1) //that is, local node j it is a useful node:
		    {
		      active_node_in_enrichment_shape_function(  aux_node_enrichment_shape_function_index(j)  ) = j;

		      if(partition_father_nodes(j,0)==non_aux_node_father_node || partition_father_nodes(j,1)==non_aux_node_father_node)
			{
			  if (aux_node_enrichment_shape_function_index(j)==0)
			    useful_node_for_N0star=j;
			  if (aux_node_enrichment_shape_function_index(j)==1)
			    useful_node_for_N1star=j;
			}
		    }
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
	      //KRATOS_ERROR(std::logic_error, "method not implemented", "");
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
		      NEnriched(partition_number, index_shape_function+3 ) = one_third ; //only one gauss point. if more were to be used, it would still be simple (1/2,1/2,0);(0,1/2,1/2);(1/2,0,1/2);
		      for (int j = 0; j < 2; j++) //x,y,(z)
			{
			  rGradientsValue[partition_number](index_shape_function+3, j) = DN_DX_subdomain(active_node_in_enrichment_shape_function(index_shape_function),j);
			}

		    }
		  else
		    {
		      NEnriched(partition_number, index_shape_function+3 ) = 0.0;
		      for (int j = 0; j < 2; j++) //x,y,(z)
			{
			  rGradientsValue[partition_number](index_shape_function+3, j) = 0.0;
			}
		    }

		}

	      array_1d<int,(2+1)> replacement_shape_function_nodes = ZeroVector(3);
	      for (int index_shape_function = 0; index_shape_function < 3; index_shape_function++) //replacement shape function
		{
		  int active_node=-1;
		  for (int j = 0; j < 3; j++) //number_of_nodes in the subelement
		    if (partition_father_nodes(j,0)==index_shape_function && partition_father_nodes(j,1)==index_shape_function)
		      {
			active_node=j;
			break;
		      }
		  if(active_node> -1)
		    {
		      for (int j = 0; j < 2; j++) //x,y,(z)
			rGradientsValue[partition_number](index_shape_function, j) = DN_DX_subdomain(active_node,j);
		      NEnriched(partition_number, index_shape_function ) = one_third;
		      replacement_shape_function_nodes(index_shape_function) = active_node;
		    }
		  else
		    {
		      for (int j = 0; j < 2; j++) //x,y,(z)
			rGradientsValue[partition_number](index_shape_function, j) = 0.0;
		      replacement_shape_function_nodes(index_shape_function) = -1;
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

	      if(useful_node_for_N0star > -1)
		{
		  for (int j = 0; j < 2; j++) //x,y,(z)
		    {
		      if(partition_sign>0)
			{
			  //first two rows are for the side where N*1 = 1
			  //row 1 is for gradN*1
			  rGradientpositive(3, j)=DN_DX_subdomain(useful_node_for_N0star, j);
			  //row 2 is for gradN*2
			  if (active_node_in_enrichment_shape_function(1) > -1)
			    rGradientpositive(4, j)=DN_DX_subdomain(active_node_in_enrichment_shape_function(1), j);

			  for (int index_shape_function = 0; index_shape_function < 3; index_shape_function++) //replacement shape function
			    {
			      if(replacement_shape_function_nodes(index_shape_function)>-1)
				{
				  rGradientpositive(index_shape_function, j) = DN_DX_subdomain(replacement_shape_function_nodes(index_shape_function), j);
				}
			    }
			}
		      else
			{
			  rGradientnegative(3, j) =DN_DX_subdomain(useful_node_for_N0star, j);
			  if (active_node_in_enrichment_shape_function(1) > -1)
			    rGradientnegative(4, j) =DN_DX_subdomain(active_node_in_enrichment_shape_function(1), j);

			  for (int index_shape_function = 0; index_shape_function < 3; index_shape_function++) //replacement shape function
			    {
			      if(replacement_shape_function_nodes(index_shape_function)>-1)
				{
				  rGradientnegative(index_shape_function, j) = DN_DX_subdomain(replacement_shape_function_nodes(index_shape_function), j);
				}
			    }
			}
		    }
		}
	      if(useful_node_for_N1star > -1)
		{
		  for (int j = 0; j < 2; j++) //x,y,(z)
		    {
		      if(partition_sign>0)
			{
			  //rows 3 and 4 are for the side where N*2 = 1
			  //row 4 is for gradN*2
			  rGradientpositive(9, j)=DN_DX_subdomain(useful_node_for_N1star, j);
			  //row 3 is for gradN*1
			  if (active_node_in_enrichment_shape_function(0) > -1)
			    rGradientpositive(8, j)=DN_DX_subdomain(active_node_in_enrichment_shape_function(0), j);

			  for (int index_shape_function = 0; index_shape_function < 3; index_shape_function++) //replacement shape function
			    {
			      if(replacement_shape_function_nodes(index_shape_function)>-1)
				{
				  rGradientpositive(index_shape_function+5, j) = DN_DX_subdomain(replacement_shape_function_nodes(index_shape_function), j);
				}
			    }
			}
		      else
			{
			  rGradientnegative(9, j)=DN_DX_subdomain(useful_node_for_N1star, j);
			  if(active_node_in_enrichment_shape_function(0) > -1)
			    rGradientnegative(8, j)=DN_DX_subdomain(active_node_in_enrichment_shape_function(0), j);

			  for (int index_shape_function = 0; index_shape_function < 3; index_shape_function++) //replacement shape function
			    {
			      if(replacement_shape_function_nodes(index_shape_function)>-1)
				{
				  rGradientnegative(index_shape_function+5, j) = DN_DX_subdomain(replacement_shape_function_nodes(index_shape_function), j);
				}
			    }
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


    //2D: 2 enrichment functions for capturing weak discontinities. All the shape functions follow the criteria of Partition of Unity.
    static int CalculateEnrichedShapeFuncionsExtendedmodified_gausspoints(Geometry< Node<3> >& trianglegeom,BoundedMatrix<double,(2+1), 2 >& rPoints, BoundedMatrix<double, (2+1), 2 >& DN_DX,array_1d<double,(2+1)>& rDistances, array_1d<double,(3*(2-1))>& rVolumes, BoundedMatrix<double, 3*(2-1), (2+1) >& rGPShapeFunctionValues, array_1d<double,(3*(2-1))>& rPartitionsSign, std::vector<Matrix>& rGradientsValue, BoundedMatrix<double,3*(2-1), (5)>& NEnriched,BoundedMatrix<double,10, 2>& rGradientpositive,BoundedMatrix<double,10, 2>& rGradientnegative ,BoundedMatrix<int,3,3>& father_nodes,std::vector<Matrix>& PRUEBA, array_1d<double,6>& weight)
    {
      KRATOS_TRY

	const double one_third=1.0/3.0;
      BoundedMatrix<double,3,2> aux_points; //for auxiliary nodes 4(between 1 and 2) ,5(between 2 and 3) ,6 (between 3 and 1)
      BoundedMatrix<double, 3, 2 > coord_subdomain; //used to pass arguments when we must calculate areas, shape functions, etc
      BoundedMatrix<double, 3, 2 > coord_subdomain_aux;
      BoundedMatrix<double,3,2> DN_DX_subdomain; //used to retrieve derivatives

      double most_common_sign=0; //the side of the cut in which two nodes are found (same sign) will be the ones that remains unchanged when builing the discontinuity
      double Area;//area of the complete element
      rGPShapeFunctionValues(0,0)=one_third; rGPShapeFunctionValues(0,1)=one_third; rGPShapeFunctionValues(0,2)=one_third; //default, when no interfase has been found
      Area = CalculateVolume2D( rPoints );
      array_1d<bool,3> cut_edges;
      array_1d<double,3> aux_nodes_relative_locations;
      BoundedMatrix<int,3,2> aux_nodes_father_nodes;
      BoundedMatrix<double,3,2> DN_DX_subdomainaux_1aux; //used to retrieve derivatives

      //to begin with we must check whether our element is cut or not by the interfase.
      if( (rDistances(0)*rDistances(1))>0.0 && (rDistances(0)*rDistances(2))>0.0 ) //it means that this element IS NOT cut by the interfase. we must return data of a normal, non-enriched element
	{
	  rVolumes(0)=Area;
	  rGPShapeFunctionValues(0,0)=one_third; rGPShapeFunctionValues(0,1)=one_third; rGPShapeFunctionValues(0,2)=one_third;
	  NEnriched(0,0) = 0.0;
	  //type_of_cut=1;
	  for (int j = 0; j < 2; j++)
	    rGradientsValue[0](0, j) = 0.0;
	  if (rDistances(0) < 0.0) rPartitionsSign[0] = -1.0;
	  else rPartitionsSign[0] = 1.0;
	  //KRATOS_WATCH("one element not in the intefase")
	  return 1;
	}

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
      //const double tolerable_distance =longest_distance*0.001;// 0.001 	// (1/1,000,000 seems to have good results)

      //const double epsilon = 1e-15; //1.00e-9;
      //compute the gradient of the distance and normalize it
      array_1d<double, 2 > grad_d;
      noalias(grad_d) = prod(trans(DN_DX), rDistances);

      array_1d<double, 3> exact_distance = rDistances;
      array_1d<double, 3> abs_distance = ZeroVector(3);


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

      //We have 3 edges, meaning we created 3 aux nodes. But one of them actually matches the position of a real node (the one that is not on an interface edge is displaced to one of the ends (a node)
      //the new shape functions are built by setting the values in all (real and aux) nodes to zero except in one of the interphase nodes. in the array aux_node_shape_function_index we assign this value
      array_1d<int, 3 > aux_node_enrichment_shape_function_index; //when not used, it must be -1;

      int shape_function_id=0;
      father_nodes(0,0)=-1;
      father_nodes(0,1)=-1;
      father_nodes(0,2)=-1;
      father_nodes(1,0)=-1;
      father_nodes(1,1)=-1;
      father_nodes(1,2)=-1;
      father_nodes(2,0)=-1;
      father_nodes(2,1)=-1;
      father_nodes(2,2)=-1;

      //KRATOS_WATCH(father_nodes);
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
	      father_nodes(i,0)=edge_begin_node;
	      father_nodes(i,1)=edge_end_node;
	      father_nodes(i,2)=shape_function_id;
	      shape_function_id++;


	    }
	  else
	    {//<
	      if(fabs(rDistances(edge_end_node))<fabs(rDistances(edge_begin_node))) //if edge is not cut, we collapse the aux node into the node which has the highest absolute value to have "nicer" (less "slivery") subelements
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
	  for (unsigned int j=0;j<2;j++){	//x,y coordinates
	    aux_points(i,j)= rPoints(edge_begin_node,j) * aux_nodes_relative_locations(i) + rPoints(edge_end_node,j) * (1.0- aux_nodes_relative_locations(i));
	  }
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
	  double d =    (rPoints(i_node,0) - base_point[0]) * grad_d[0] + (rPoints(i_node,1) - base_point[1]) * grad_d[1] ;
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
      rGradientsValue[0]=ZeroMatrix(5,2);
      rGradientsValue[1]=ZeroMatrix(5,2);
      rGradientsValue[2]=ZeroMatrix(5,2);

      NEnriched=ZeroMatrix(3,5);
      rGPShapeFunctionValues=ZeroMatrix(3,3);

      //now we must check the 4 created partitions of the domain.
      //one has been collapsed, so we discard it and therefore save only one.
      unsigned int partition_number=0;		//
      //the 3 first partitions are  created using 2 auxiliary nodes and a normal node. at least one of these will be discarded due to zero area
      //the last one is composed by the 3 auxiliary nodes. it 'looks' wrong, but since at least one has been collapsed, it actually has a normal node.
      bool found_empty_partition=false;

      //the enrichment is directly the shape functions created by using the partition.
      //we have to save for the enrichment 0 and 1 which is the node that will be active, that is, whose shape function value is zero (for some partitions all 3 shape functions are inactive)

      for ( int i=0; i<4; i++) //i partition
	{

	  array_1d<int, 2 > active_node_in_enrichment_shape_function;
	  active_node_in_enrichment_shape_function(0)=-1;  active_node_in_enrichment_shape_function(1)=-1; //initialized as if all are inactive -> gradient=0;
	  //the same but for the replacement shape functions
	  array_1d<int, 3 > active_node_in_replacement_shape_function;
	  active_node_in_replacement_shape_function(0)=-1;  active_node_in_replacement_shape_function(1)=-1; active_node_in_replacement_shape_function(2)=-1; //initialized as if all are inactive -> gradient=0;

	  int j_aux = i + 2;
	  if (j_aux>2) j_aux -= 3;
	  BoundedMatrix<int,3,2> partition_father_nodes;
	  array_1d<double,3> N; //bool useful=false;

	  int useful_node_for_N0star=-1;
	  int useful_node_for_N1star=-1;

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

	      //KRATOS_WATCH(coord_subdomain)
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
	    //  if( (aux_nodes_father_nodes(i,0)!=aux_nodes_father_nodes(i,1)) && (aux_nodes_father_nodes(j_aux,0)!=aux_nodes_father_nodes(j_aux,1)))
		//{
		//  useful=true;
		//}

	      coord_subdomain(0,0)=rPoints(i,0);
	      coord_subdomain(0,1)=rPoints(i,1);
	      coord_subdomain(1,0)=aux_points(i,0);
	      coord_subdomain(1,1)=aux_points(i,1);
	      coord_subdomain(2,0)=aux_points(j_aux,0);
	      coord_subdomain(2,1)=aux_points(j_aux,1);


	      //an aux_node is useful when one of its father nodes is the real node of the subelement. that means that the edge is part of the subelement.
	      if(partition_father_nodes(1,0)==i || partition_father_nodes(1,1)==i) //if one of the father nodes of node_aux_i is equal to the real node i
		{
		  if(aux_node_enrichment_shape_function_index(i)==0)
		    useful_node_for_N0star=1;
		  if(aux_node_enrichment_shape_function_index(i)==1)
		    useful_node_for_N1star=1;
		}
	      if(partition_father_nodes(2,0)==j_aux || partition_father_nodes(2,1)==j_aux) //if one of the father nodes of node_aux_i is equal to the real node i
		{
		  if(aux_node_enrichment_shape_function_index(j_aux)==0)
		    useful_node_for_N0star=2;
		  if(aux_node_enrichment_shape_function_index(j_aux)==1)
		    useful_node_for_N1star=2;
		}
	    }
	  else
	    {
	      //the last partition, made by the 3 aux nodes.
	      partition_father_nodes=aux_nodes_father_nodes;
	      coord_subdomain=aux_points;
	      //useful=true;
	      //int non_aux_node=-1; //one of the nodes of this partition is actually a real node, which has repeated father node
	      int non_aux_node_father_node=-1; //the real node id
	      //we have to check the 3 of them:
	      for (int j = 0; j < 3; j++)
		{
		  if(partition_father_nodes(j,0)==partition_father_nodes(j,1))
		    {
		      //non_aux_node=j;
		      non_aux_node_father_node=partition_father_nodes(j,0);
		    }
		}
	      for (int j = 0; j < 3; j++)
		{
		  if (aux_node_enrichment_shape_function_index(j)> -1) //that is, local node j it is a useful node:
		    {
		      active_node_in_enrichment_shape_function(  aux_node_enrichment_shape_function_index(j)  ) = j;

		      if(partition_father_nodes(j,0)==non_aux_node_father_node || partition_father_nodes(j,1)==non_aux_node_father_node)
			{
			  if (aux_node_enrichment_shape_function_index(j)==0)
			    useful_node_for_N0star=j;
			  if (aux_node_enrichment_shape_function_index(j)==1)
			    useful_node_for_N1star=j;
			}
		    }

		  //to replace the standard shape functions:
		  if (aux_nodes_father_nodes(j,0)==aux_nodes_father_nodes(j,1))
		    active_node_in_replacement_shape_function(aux_nodes_father_nodes(j,0))=j;
		}

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
	      coord_subdomain_aux = coord_subdomain;
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

	      typedef GeometryData::IntegrationMethod IntegrationMethod;
	      IntegrationMethod mThisIntegrationMethod;
	      std::vector< Matrix > mInvJ0;
	      Vector mDetJ0;

	      Geometry< Node<3> >::PointsArrayType NewPoints;
	      Node<3>::Pointer p_new_node;
	      int id_local=0;
	      for (unsigned int j=0; j!=3; j++)
		{
		  p_new_node = Node<3>::Pointer(new Node<3>(id_local, coord_subdomain_aux(j,0), coord_subdomain_aux(j,1), coord_subdomain_aux(j,2)));
		  NewPoints.push_back(p_new_node);
		  id_local++;
		}

	      Geometry< Node<3> >::Pointer pGeom = trianglegeom.Create(NewPoints);
	      //const unsigned int number_of_points = pGeom->size();
	      //KRATOS_WATCH(number_of_points);
	      mThisIntegrationMethod= GeometryData::GI_GAUSS_2;

	      typedef Geometry<Node<3> >::IntegrationPointsArrayType IntegrationPointsArrayType;
	      //const GeomtryType::IntegrationPointsArrayType& integration_points = pGeom->IntegrationPoints(mThisIntegrationMethod);
	      const IntegrationPointsArrayType& integration_points = pGeom->IntegrationPoints(mThisIntegrationMethod);

	      mInvJ0.resize(integration_points.size());
	      mDetJ0.resize(integration_points.size(),false);

	      //KRATOS_WATCH(integration_points.size());
	      //KRATOS_ERROR(std::logic_error, "element with zero area found", "");
	      Element::GeometryType::JacobiansType J0;
	      J0 = pGeom->Jacobian(J0, mThisIntegrationMethod);

	      const Matrix& Ncontainer = pGeom->ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

	      double xp=0.0;
	      double yp=0.0;
	      double temp_areaaux_2=0.0;
	      bounded_matrix<double, 4, 8 > N_new=ZeroMatrix(4,8);

	      for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
	  	{
		  //unsigned int number_of_nodes = pGeom->PointsNumber();
		  //unsigned int dimension = pGeom->WorkingSpaceDimension();
		  const Vector& N=row(Ncontainer,PointNumber);

		  MathUtils<double>::InvertMatrix(J0[PointNumber],mInvJ0[PointNumber],mDetJ0[PointNumber]);
		  double Weight = integration_points[PointNumber].Weight()* mDetJ0[PointNumber];

		  xp=0.0;
		  yp=0.0;

		  for (unsigned int tt=0; tt<3; tt++)
		    {
		      xp += N[tt] * coord_subdomain_aux(tt,0);
		      yp += N[tt] * coord_subdomain_aux(tt,1);
		    }

		  for (unsigned int h=0; h<3; h++)  //loop por los nodos del subelemento
		    {

		      bounded_matrix<double, 3, 2 > original_coordinates; //8 is the max number of nodes and aux_nodes
		      for (unsigned int i = 0; i < 3; i++)
			for (unsigned int j = 0; j < 2; j++)
			  original_coordinates(i, j) = rPoints(i, j);


		      original_coordinates(h,0) = xp;
		      original_coordinates(h,1) = yp;


		      CalculateGeometryData(original_coordinates, DN_DX_subdomainaux_1aux, temp_areaaux_2);//222

		      PRUEBA[partition_number](PointNumber,h) = temp_areaaux_2/Area;// * Weight;
		      weight(partition_number)=Weight;
		    }

		}

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
		      NEnriched(partition_number, index_shape_function+3 ) = one_third ; //only one gauss point. if more were to be used, it would still be simple (1/2,1/2,0);(0,1/2,1/2);(1/2,0,1/2);
		      for (int j = 0; j < 2; j++) //x,y,(z)
			{
			  rGradientsValue[partition_number](index_shape_function+3, j) = DN_DX_subdomain(active_node_in_enrichment_shape_function(index_shape_function),j);
			}

		    }
		  else
		    {
		      NEnriched(partition_number, index_shape_function+3 ) = 0.0;
		      //NEnriched(partition_number, index_shape_function +2) = 0.0;
		      for (int j = 0; j < 2; j++) //x,y,(z)
			{
			  rGradientsValue[partition_number](index_shape_function+3, j) = 0.0;
			  //KRATOS_WATCH(rGradientsValue);
			  //rGradientsValue[partition_number](index_shape_function+2, j) = 0.0;
			}
		    }

		}

	      array_1d<int,(2+1)> replacement_shape_function_nodes = ZeroVector(3);
	      for (int index_shape_function = 0; index_shape_function < 3; index_shape_function++) //replacement shape function
		{
		  int active_node=-1;
		  for (int j = 0; j < 3; j++) //number_of_nodes in the subelement
		    if (partition_father_nodes(j,0)==index_shape_function && partition_father_nodes(j,1)==index_shape_function)
		      {
			active_node=j;
			break;
		      }
		  if(active_node> -1)
		    {
		      for (int j = 0; j < 2; j++) //x,y,(z)
			rGradientsValue[partition_number](index_shape_function, j) = DN_DX_subdomain(active_node,j);
		      NEnriched(partition_number, index_shape_function ) = one_third;
		      replacement_shape_function_nodes(index_shape_function) = active_node;
		    }
		  else
		    {
		      for (int j = 0; j < 2; j++) //x,y,(z)
			rGradientsValue[partition_number](index_shape_function, j) = 0.0;
		      replacement_shape_function_nodes(index_shape_function) = -1;
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

	      if(useful_node_for_N0star > -1)
		{
		  for (int j = 0; j < 2; j++) //x,y,(z)
		    {
		      if(partition_sign>0)
			{
			  //first two rows are for the side where N*1 = 1
			  //row 1 is for gradN*1
			  rGradientpositive(3, j)=DN_DX_subdomain(useful_node_for_N0star, j);
			  //row 2 is for gradN*2
			  if (active_node_in_enrichment_shape_function(1) > -1)
			    rGradientpositive(4, j)=DN_DX_subdomain(active_node_in_enrichment_shape_function(1), j);

			  for (int index_shape_function = 0; index_shape_function < 3; index_shape_function++) //replacement shape function
			    {
			      if(replacement_shape_function_nodes(index_shape_function)>-1)
				{
				  rGradientpositive(index_shape_function, j) = DN_DX_subdomain(replacement_shape_function_nodes(index_shape_function), j);
				}
			    }
			}
		      else
			{
			  rGradientnegative(3, j) =DN_DX_subdomain(useful_node_for_N0star, j);
			  //KRATOS_WATCH(rGradientnegative);
			  if (active_node_in_enrichment_shape_function(1) > -1)
			    rGradientnegative(4, j) =DN_DX_subdomain(active_node_in_enrichment_shape_function(1), j);

			  for (int index_shape_function = 0; index_shape_function < 3; index_shape_function++) //replacement shape function
			    {
			      if(replacement_shape_function_nodes(index_shape_function)>-1)
				{
				  rGradientnegative(index_shape_function, j) = DN_DX_subdomain(replacement_shape_function_nodes(index_shape_function), j);
				}
			    }
			}
		    }
		}
	      if(useful_node_for_N1star > -1)
		{
		  for (int j = 0; j < 2; j++) //x,y,(z)
		    {
		      if(partition_sign>0)
			{
			  //rows 3 and 4 are for the side where N*2 = 1
			  //row 4 is for gradN*2
			  rGradientpositive(9, j)=DN_DX_subdomain(useful_node_for_N1star, j);
			  //row 3 is for gradN*1
			  if (active_node_in_enrichment_shape_function(0) > -1)
			    rGradientpositive(8, j)=DN_DX_subdomain(active_node_in_enrichment_shape_function(0), j);
			  //KRATOS_WATCH(rGradientpositive);

			  for (int index_shape_function = 0; index_shape_function < 3; index_shape_function++) //replacement shape function
			    {
			      if(replacement_shape_function_nodes(index_shape_function)>-1)
				{
				  rGradientpositive(index_shape_function+5, j) = DN_DX_subdomain(replacement_shape_function_nodes(index_shape_function), j);
				}
			    }
			}
		      else
			{
			  rGradientnegative(9, j)=DN_DX_subdomain(useful_node_for_N1star, j);
			  if(active_node_in_enrichment_shape_function(0) > -1)
			    rGradientnegative(8, j)=DN_DX_subdomain(active_node_in_enrichment_shape_function(0), j);
			  //KRATOS_WATCH(rGradientnegative);

			  for (int index_shape_function = 0; index_shape_function < 3; index_shape_function++) //replacement shape function
			    {
			      if(replacement_shape_function_nodes(index_shape_function)>-1)
				{
				  rGradientnegative(index_shape_function+5, j) = DN_DX_subdomain(replacement_shape_function_nodes(index_shape_function), j);
				}
			    }
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

      return 3;
      KRATOS_CATCH("");

    }

    //for 3D
    static int CalculateEnrichedShapeFuncions_Simplified(Geometry< Node<3> >& rGeom,BoundedMatrix<double,(3+1), 3 >& rPoints, BoundedMatrix<double, (3+1), 3 >& DN_DX, array_1d<double,(3+1)>& rDistances, array_1d<double,(3*(3-1))>& rVolumes, BoundedMatrix<double, 3*(3-1), (3+1) >& rShapeFunctionValues,array_1d<double,(3*(3-1))>& rPartitionsSign, std::vector<Matrix>& rGradientsValue, std::vector<Matrix>& rGradientsValueaux,     BoundedMatrix<double,3*(3-1), (2)>& NEnriched,int& number_interface_elements,BoundedMatrix<double, 2, 3 >& coord_interface_nodes, array_1d<double,6>& area_interface, array_1d<double,6>& area_inter, array_1d<double,6>& N_Star, bool& switch_off_e, std::vector<Matrix>&  edges_t,std::vector<Matrix>&  nodes,std::vector<Matrix>&  original_edges, std::vector<Matrix>& rGradientaux1, int& totalnodes, std::vector<Matrix>& interface_nodes, BoundedMatrix<double, 3*(3-1), 8 >& Ngauss_new, std::vector<Matrix>& Tres, std::vector<Matrix>& PRUEBA, array_1d<double,6>& weight)
    {
      KRATOS_TRY

        const int n_nodes = 4; // it works only for tetrahedra
      const int n_edges = 6; // it works only for tetrahedra
      bool switch_off_ee=false;
      const int edge_i[] = {0, 0, 0, 1, 1, 2};
      const int edge_j[] = {1, 2, 3, 2, 3, 3};
      //int z=0;
      //int n=0;
      std::vector< Matrix > edges_o(8);

      for (unsigned int i = 0; i < 8; i++)
	{
	  edges_o[i].resize(5, 3, false);  //2 values of the 2 shape functions, and derivates in (xyz) direction).
	}

      const double epsilon = 1e-15; //1.00e-9;
      const double near_factor = 1.00e-12;

      bounded_matrix<double, 4, 3 > coord_node=ZeroMatrix(4,3);

      int number_of_partitions = 1;
      for (unsigned int i = 0; i < 8; i++) //// 5
	for (unsigned int j = 0; j < 5; j++)
	  for (unsigned int k = 0; k < 3; k++)
	    edges_t[i](j,k) =-1.0;

      for (unsigned int i = 0; i < 8; i++) //// 5
	for (unsigned int j = 0; j < 5; j++)
	  for (unsigned int k = 0; k < 3; k++)
	    original_edges[i](j,k) =-1.0;

      for (unsigned int i = 0; i < 6; i++)    //4
      	for (unsigned int j = 0; j < 8; j++)
	  for (unsigned int k = 0; k < 3; k++)
	    rGradientaux1[i](j,k) =0.0;

      array_1d<double, n_edges> edges_dx; // It will be initialize later
      array_1d<double, n_edges> edges_dy; // It will be initialize later
      array_1d<double, n_edges> edges_dz; // It will be initialize later
      array_1d<double, n_edges> edges_length; // It will be initialize later
      // The divided part length from first node of edge respect to the edge length
      array_1d<double, n_edges> edge_division_i = ZeroVector(n_edges); // The 0 is for no split
      // The divided part length from second node of edge respect to the edge length
      array_1d<double, n_edges> edge_division_j = ZeroVector(n_edges); // The 0 is for no split

      bounded_matrix<double, 8, 3 > aux_coordinates; //8 is the max number of nodes and aux_nodes
      for (unsigned int i = 0; i < 4; i++)
	for (unsigned int j = 0; j < 3; j++)
	  aux_coordinates(i, j) = rPoints(i, j);
      for (unsigned int i = 4; i < 8; i++)
	for (unsigned int j = 0; j < 3; j++)
	  aux_coordinates(i, j) = -10000.0; //set to a large number so that errors will be evident

      int split_edge[] = {0, 1, 2, 3, -1, -1, -1, -1, -1, -1, -1, -1};
      int new_node_id = 4;
      bounded_matrix<double, 4, 4 > length = ZeroMatrix(4, 4);

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

      for (int i = 0; i < 6; i++)
	for (int j = 0; j < 8; j++)
	  Ngauss_new(i, j) = 0.0;
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
            }
        }
      totalnodes=new_node_id;
      for(int aux = 0; aux < 4; aux++)
        {
	  nodes[aux](0,0)=rPoints(aux, 0);
	  nodes[aux](0,1)=rPoints(aux, 1);
	  nodes[aux](0,2)=rPoints(aux, 2);

	}

      for(int aux = 4; aux < new_node_id; aux++)
        {
	  nodes[aux](0,0)=aux_coordinates(aux, 0);
	  nodes[aux](0,1)=aux_coordinates(aux, 1);
	  nodes[aux](0,2)=aux_coordinates(aux, 2);

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

      if (number_of_splitted_edges == 0) // no splitting
        {

	  rVolumes[0] = volume;
	  sub_volumes_sum = volume;
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
	  int nel; //number of elements generated
	  int n_splitted_edges; //number of splitted edges
	  int nint; //number of internal nodes
	  int t[56];
	  TetrahedraSplit::Split_Tetrahedra(edge_ids, t, &nel, &n_splitted_edges, &nint);

	  if (nint != 0)
	    KRATOS_ERROR<<"requiring an internal node for splitting ... can not accept this";

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

      BoundedMatrix<double,4,3> DN_DX_subdomainaux1; //used to retrieve derivatives
      BoundedMatrix<double, 4, 3 > coord_subdomainaux1;
      double temp_areaaux1=0.0;


      CalculateGeometryData(rPoints, DN_DX_subdomainaux1, temp_areaaux1);//222

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

	  if(max_aux_dist_on_cut < 1e-9*max_lenght)
	    max_aux_dist_on_cut =  1e-9*max_lenght;

	  int edge_ids[6];
	  TetrahedraSplit::TetrahedraSplitMode(split_edge, edge_ids);
	  int nel; //number of elements generated
	  int n_splitted_edges; //number of splitted edges
	  int nint; //number of internal nodes
	  int t[56];
	  TetrahedraSplit::Split_Tetrahedra(edge_ids, t, &nel, &n_splitted_edges, &nint);

	  if (nint != 0)
	    KRATOS_ERROR<<"requiring an internal node for splitting ... can not accept this";
	  coord_node=ZeroMatrix(2,3);
	  coord_interface_nodes=ZeroMatrix(2,3);

	  array_1d<double, 3 > center_position;
	  N_Star= ZeroVector(6);

	  area_interface= ZeroVector(6);
	  area_inter=ZeroVector(6);

	  BoundedMatrix<double,4,3> DN_DX_subdomainaux; //used to retrieve derivatives
	  BoundedMatrix<double,4,3> DN_DX_subdomainaux_1; //used to retrieve derivatives
	  BoundedMatrix<double,4,3> DN_DX_subdomainaux_1aux; //used to retrieve derivatives
	  BoundedMatrix<double, 4, 3 > coord_subdomainaux;
	  BoundedMatrix<double, 4, 3 > coord_subdomainaux_aux;
	  double temp_areaaux=0.0;
	  int local=0;

	  for (unsigned int i = 0; i < 6; i++)
	    {
	      PRUEBA[i].resize(4, 4, false);
	      PRUEBA[i] *=0.0;
	    }

	  for (unsigned int i = 0; i < 6; i++)
	    {
	      Tres[i] *=0.0;
	    }

	  array_1d<double, 4 > N;
	  bounded_matrix<double, 4, 3 > original_coordinates; //8 is the max number of nodes and aux_nodes
	  for (unsigned int i = 0; i < 4; i++)
	    for (unsigned int j = 0; j < 3; j++)
	      original_coordinates(i, j) = rPoints(i, j);

	  ///para el elemento global
	  double temp_areaaux_1=0.0;
	  CalculateGeometryData(rPoints, DN_DX_subdomainaux_1aux, temp_areaaux_1);//222


	  for (int i = 0; i < number_of_partitions; i++)
	    {
	      //compute enriched shape function values

	      int i0, i1, i2, i3; //indices of the subtetrahedra
	      TetrahedraSplit::TetrahedraGetNewConnectivityGID(i, t, split_edge, &i0, &i1, &i2, &i3);

	      BoundedMatrix<double, 4, 3 > coord_subdomain; //used to pass arguments when we must calculate areas, shape functions, etc                					BoundedMatrix<double,4,3> DN_DX_subdomain; //used to retrieve derivatives
	      BoundedMatrix<double, 4, 3 > coord_xg; //use

	      coord_xg=ZeroMatrix(4,3);

	      array_1d<int,4> partition_nodes;
	      //double temp_area;
	      partition_nodes[0]=i0;
	      partition_nodes[1]=i1;
	      partition_nodes[2]=i2;
	      partition_nodes[3]=i3;
	      //double interface_area=0.0;

	      int nummm=0;
	      for (unsigned int j=0; j!=4; j++)
		{//4 nodes of the partition
		  const int node_id=partition_nodes[j];
		  if(node_id>3){
		    nummm++;
		  }
		}

	      std::vector<array_1d<double,3> > PointsOfFSTriangle;
	      PointsOfFSTriangle.reserve(3);

	      //int position=0;

	      for (unsigned int j=0; j!=4; j++)
		{
		  for (unsigned int k=0; k!=3; k++) {//x,y,z
		    const int node_id=partition_nodes[j];
		    coord_subdomainaux(j,k)= aux_coordinates( node_id ,k );
		  }
		}


	      for (unsigned int i = 0; i < 3; i++)
		{
		  center_position[i] = aux_coordinates(0, i) + aux_coordinates(1, i) + aux_coordinates(2, i) + aux_coordinates(3, i);
		}
	      center_position *= 0.25;

	      CalculateGeometryData(coord_subdomainaux, DN_DX_subdomainaux, temp_areaaux);//222

	      typedef GeometryData::IntegrationMethod IntegrationMethod;
	      IntegrationMethod mThisIntegrationMethod;
              std::vector< Matrix > mInvJ0;
              Vector mDetJ0;


	      array_1d<double, 4 > msN;
              BoundedMatrix<double, 4, 3 > msDN_DX;
	      //double Area=0.0;

	      Geometry< Node<3> >::PointsArrayType NewPoints;
	      Node<3>::Pointer p_new_node;
	      int id_local=0;
	      for (unsigned int j=0; j!=4; j++)
		{
		  id_local=partition_nodes[j];
		  p_new_node = Node<3>::Pointer(new Node<3>(id_local, coord_subdomainaux(j,0), coord_subdomainaux(j,1), coord_subdomainaux(j,2)));
		  NewPoints.push_back(p_new_node);
		  id_local++;
		}
              Geometry< Node<3> >::Pointer pGeom = rGeom.Create(NewPoints);
      	      //const unsigned int number_of_points = pGeom->size();
	      //number of gauss points
	      mThisIntegrationMethod= GeometryData::GI_GAUSS_2;

	      typedef Geometry<Node<3> >::IntegrationPointsArrayType IntegrationPointsArrayType;

	      const IntegrationPointsArrayType& integration_points = pGeom->IntegrationPoints(mThisIntegrationMethod);

	      mInvJ0.resize(integration_points.size());
	      mDetJ0.resize(integration_points.size(),false);

	      Element::GeometryType::JacobiansType J0;
	      J0 = pGeom->Jacobian(J0, mThisIntegrationMethod);

	      const Matrix& Ncontainer = pGeom->ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
	      double xp=0.0;
	      double yp=0.0;
	      double zp=0.0;

	      double temp_areaaux_2=0.0;
	      bounded_matrix<double, 4, 8 > N_new=ZeroMatrix(4,8);

	      for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
	  	{
		  //unsigned int number_of_nodes = pGeom->PointsNumber();
		  //unsigned int dimension = pGeom->WorkingSpaceDimension();
		  const Vector& N=row(Ncontainer,PointNumber);

		  MathUtils<double>::InvertMatrix(J0[PointNumber],mInvJ0[PointNumber],mDetJ0[PointNumber]);
		  double Weight = integration_points[PointNumber].Weight()* mDetJ0[PointNumber];

		  xp=0.0;
		  yp=0.0;
		  zp=0.0;

		  for (unsigned int tt=0; tt!=4; tt++)
		    {
		      xp += N[tt] * coord_subdomainaux(tt,0);
		      yp += N[tt] * coord_subdomainaux(tt,1);
		      zp += N[tt] * coord_subdomainaux(tt,2);
		    }

		  for (unsigned int j=0; j!=4; j++)  //loop por los nodos del subelemento
		    {
		      bounded_matrix<double, 8, 3 > aux_coordinates_aux;
		      aux_coordinates_aux=aux_coordinates;

		      bounded_matrix<double, 4, 3 > original_coordinates; //8 is the max number of nodes and aux_nodes
		      for (unsigned int i = 0; i < 4; i++)
	  		for (unsigned int j = 0; j < 3; j++)
			  original_coordinates(i, j) = rPoints(i, j);


		      const int node_idaux=partition_nodes[j];

		      aux_coordinates_aux(node_idaux,0)=xp;
		      aux_coordinates_aux(node_idaux,1)=yp;
		      aux_coordinates_aux(node_idaux,2)=zp;

		      original_coordinates(j,0) = xp;
		      original_coordinates(j,1) = yp;
		      original_coordinates(j,2) = zp;

		      for (unsigned int jj=0; jj!=4; jj++)
			{
			  for (unsigned int k=0; k!=3; k++)
			    {//x,y,z
			      const int node_id=partition_nodes[jj];
			      coord_subdomainaux_aux(jj,k)= aux_coordinates_aux( node_id ,k );
			    }
			}

		      ///para las coordenadas globales
		      CalculateGeometryData(original_coordinates, DN_DX_subdomainaux_1aux, temp_areaaux_2);//222

		      PRUEBA[i](PointNumber,j) = temp_areaaux_2/temp_areaaux_1;
		      weight(i)=Weight;

		    }

		  //int id_local_aux1=0;
		  //int id_local_aux2=0;

		}



	      bounded_matrix<double, 4, 3 > edges=ZeroMatrix(4,3);
	      bounded_matrix<double, 4, 3 > edged_enriched= ZeroMatrix(4,3);

	      array_1d<int, 2 > edgeone;
	      array_1d<int, 2 > edgetwo;
	      array_1d<int,3> edgetriangle;
	      array_1d<int, 2 > active_node_in_enrichment_shape_function;
	      active_node_in_enrichment_shape_function(0)=-1;  active_node_in_enrichment_shape_function(1)=-1;
	      array_1d<int, 2 > active_node_in_replacement_shape_function;
	      active_node_in_replacement_shape_function(0)=-1;  active_node_in_replacement_shape_function(1)=-1;


	      edges(0,0)=partition_nodes[0];
	      edges(0,1)=partition_nodes[2];
	      edges(0,2)=partition_nodes[1];

	      edges(1,0)=partition_nodes[0];
	      edges(1,1)=partition_nodes[3];
	      edges(1,2)=partition_nodes[2];

	      edges(2,0)=partition_nodes[0];
	      edges(2,1)=partition_nodes[1];
	      edges(2,2)=partition_nodes[3];

	      edges(3,0)=partition_nodes[2];
	      edges(3,1)=partition_nodes[3];
	      edges(3,2)=partition_nodes[1];


	      int shape_function_id=0;
	      int shape_function_aux_id=0;
	      for(int i_aux=0; i_aux<4;i_aux++)
		{
		  shape_function_id=0;
		  shape_function_aux_id=0;
		  active_node_in_enrichment_shape_function(0)=-1;  active_node_in_enrichment_shape_function(1)=-1;

		  active_node_in_replacement_shape_function(0)=-1;  active_node_in_replacement_shape_function(1)=-1;

		  edgeone(0)=-1;  edgeone(1)=-1;

		  edgetwo(0)=-1;  edgetwo(1)=-1;

		  edgetriangle(0)=-1;edgetriangle(1)=-1;edgetriangle(2)=-1;

		  for (int j=0;j<3;j++)
		    {
		      if(edges(i_aux,j)<4)
			{
			  active_node_in_replacement_shape_function(shape_function_id)=edges(i_aux,j);
			  shape_function_id++;
			}
		      else
			{
			  active_node_in_enrichment_shape_function(shape_function_aux_id)=edges(i_aux,j);
			  shape_function_aux_id++;
			}
		    }


		  //int t_max = shape_function_id;

		  unsigned int index=shape_function_aux_id;

		  for (unsigned int pos=0; pos<index; pos++)
		    {
		      for (unsigned int edge_aux=0; edge_aux<6; edge_aux++)
			{
			  if (split_edge[4+edge_aux]> -1) //that is, loca
			    {
			      if(active_node_in_enrichment_shape_function(pos)==split_edge[4+edge_aux])
				{
				  if(pos==0)
				    {
				      edgeone(0)=edge_i[edge_aux];
				      edgeone(1)=edge_j[edge_aux];


				    }
				  else
				    {
				      edgetwo(0)=edge_i[edge_aux];
				      edgetwo(1)=edge_j[edge_aux];
				    }
				}
			    }
			}
		    }
		  bool well=false;
		  for (unsigned int pos=0; pos<index; pos++)
		    {
		      for (unsigned int k=0; k<2; k++)
			{
			  if(edgeone(pos)==edgetwo(k)) well=true;
			}
		    }

		  edgetriangle(0)=edgeone(0);
		  edgetriangle(1)=edgeone(1);
		  for (unsigned int pos=0; pos<2; pos++)
		    {

		      if(edgetriangle(pos)==edgetwo(0))
			{
			  edgetriangle(2)=edgetwo(1);
			}
		      if(edgetriangle(pos)==edgetwo(1)) {
			edgetriangle(2)=edgetwo(0);
		      }
		    }

		  if(index==1) well=true;
		  if(index==1)
		    {

		      for (unsigned int pos=0; pos<2; pos++)
			{

			  if(edgetriangle(pos)==active_node_in_replacement_shape_function(0))
			    {
			      edgetriangle(2)= active_node_in_replacement_shape_function(1);
			    }
			  if(edgetriangle(pos)==active_node_in_replacement_shape_function(1))
			    {
			      edgetriangle(2)=active_node_in_replacement_shape_function(0);
			    }
			}
		    }
		  if(well==true)
		    {
		      for (unsigned int jj=0; jj<3; jj++){
			original_edges[i](i_aux,jj)=edgetriangle(jj);
		      }
		    }

		  for(int aux=0; aux<2;aux++)
		    {
		      for (unsigned int edge_aux=0; edge_aux<6; edge_aux++)
			{
			  if(active_node_in_replacement_shape_function(aux)==edge_i[edge_aux] or active_node_in_replacement_shape_function(aux)==edge_j[edge_aux] )
			    {
			      if (split_edge[4+edge_aux]> -1)
				{

				  if(active_node_in_enrichment_shape_function(0)==split_edge[4+edge_aux])
				    {
				      if(well==true)
					{
					  for (int j=0;j<3;j++)
					    {
					      edges_t[i](i_aux,j)=edges(i_aux,j);

					    }
					}
				      edges_o[i](i_aux,0)=edge_i[edge_aux];
				      edges_o[i](i_aux,1)=edge_j[edge_aux];
				      edges_o[i](i_aux,0)=active_node_in_enrichment_shape_function(1);
				    }
				  else if(active_node_in_enrichment_shape_function(1)==split_edge[4+edge_aux])
				    {
				      if(well==true)
					{
					  for (int j=0;j<3;j++)
					    {
					      edges_t[i](i_aux,j)=edges(i_aux,j);
					    }
					}
				      edges_o[i](i_aux,0)=edge_i[edge_aux];
				      edges_o[i](i_aux,1)=edge_j[edge_aux];
				      edges_o[i](i_aux,0)=active_node_in_enrichment_shape_function(0);
				    }

				}//
			    }
			}
		    }
		}
	      for (unsigned int j=0; j!=4; j++)
		{
		  for (unsigned int k=0; k!=3; k++) {
		    const int node_id=partition_nodes[j];
		    rGradientaux1[i](node_id,k)= DN_DX_subdomainaux(j,k);
		    Ngauss_new(i,node_id)=0.25;

		  }
		}

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

	      if(nummm==3)
	 	{
		  if(rPartitionsSign[i] == 1.0 )
		    {
		      int index=0;
		      for (unsigned int j=0; j!=4; j++)
			{
			  const int node_id=partition_nodes[j];
			  if(node_id>3) {
			    interface_nodes[local](0,index)=node_id;

			    index++;
			  }
			}
		      local++;
		    }

		}

	      NEnriched(i, 0) = 0.5 * (abs_dist - rPartitionsSign[i] * dist);
	      //normalizing
	      NEnriched(i, 0) /= max_aux_dist_on_cut;

	      for (int j = 0; j < 3; j++)
		{
		  rGradientsValue[i](0, j) = (0.5/max_aux_dist_on_cut) * (abs_distance_gradient[j] - rPartitionsSign[i] * exact_distance_gradient[j]);
		}
	    }
	}
      else
	{
	  NEnriched(0,0) = 0.0;
	  switch_off_ee=true;
	  switch_off_e=switch_off_ee;
	  for (int j = 0; j < 3; j++)
	    rGradientsValue[0](0, j) = 0.0;
	}

      return number_of_partitions;
      KRATOS_CATCH("");
    }



  private:


    static void ComputeElementCoordinates(array_1d<double, 4 > & N, const array_1d<double, 3 > & center_position, BoundedMatrix<double, 4, 3 > &  rPoints, const double vol)
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

    static inline double CalculateVol(const double x0, const double y0, const double z0,const double x1, const double y1, const double z1,const double x2, const double y2, const double z2,const double x3, const double y3, const double z3)
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
    static inline void CalculateGeometryData(const bounded_matrix<double, 3, 3 > & coordinates,BoundedMatrix<double,3,2>& DN_DX,array_1d<double,3>& N,double& Area)
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
    static inline double CalculateVolume2D( const bounded_matrix<double, 3, 3 > & coordinates)
    {
      double x10 = coordinates(1,0) - coordinates(0,0);
      double y10 = coordinates(1,1) - coordinates(0,1);

      double x20 = coordinates(2,0) - coordinates(0,0);
      double y20 = coordinates(2,1) - coordinates (0,1);
      double detJ = x10 * y20-y10 * x20;
      return 0.5*detJ;
    }

    static inline bool CalculatePosition(const bounded_matrix<double, 3, 3 > & coordinates,const double xc, const double yc, const double zc, array_1d<double, 3 > & N )
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
	  KRATOS_ERROR<<"element with zero area found";
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

    static inline double CalculateVol(const double x0, const double y0,const double x1, const double y1,const double x2, const double y2)
    {
      return 0.5 * ((x1 - x0)*(y2 - y0)- (y1 - y0)*(x2 - x0));
    }

    static inline void CalculateGeometryData(const bounded_matrix<double, 3, 3 > & coordinates,BoundedMatrix<double,3,2>& DN_DX, double& Area)
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


    static inline void CalculateGeometryData( BoundedMatrix<double, 4, 3 > & coordinates, BoundedMatrix<double,4,3>& DN_DX,double& Volume)
    {
      double x10 = coordinates(1,0) - coordinates(0,0);
      double y10 = coordinates(1,1) - coordinates(0,1);
      double z10 = coordinates(1,2) - coordinates(0,2);

      double x20 = coordinates(2,0) - coordinates(0,0);
      double y20 = coordinates(2,1) - coordinates (0,1);
      double z20 = coordinates(2,2) - coordinates (0,2);

      double x30 = coordinates(3,0) - coordinates(0,0);
      double y30 = coordinates(3,1) - coordinates (0,1);
      double z30 = coordinates(3,2) - coordinates (0,2);

      double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;

      DN_DX(0,0) = -y20 * z30 + y30 * z20 + y10 * z30 - z10 * y30 - y10 * z20 + z10 * y20;
      DN_DX(0,1) = -z20 * x30 + x20 * z30 - x10 * z30 + z10 * x30 + x10 * z20 - z10 * x20;
      DN_DX(0,2) = -x20 * y30 + y20 * x30 + x10 * y30 - y10 * x30 - x10 * y20 + y10 * x20;
      DN_DX(1,0) = y20 * z30 - y30 * z20;
      DN_DX(1,1) = z20 * x30 - x20 * z30;
      DN_DX(1,2) = x20 * y30 - y20 * x30;
      DN_DX(2,0) = -y10 * z30 + z10 * y30;
      DN_DX(2,1) = x10 * z30 - z10 * x30;
      DN_DX(2,2) = -x10 * y30 + y10 * x30;
      DN_DX(3,0) = y10 * z20 - z10 * y20;
      DN_DX(3,1) = -x10 * z20 + z10 * x20;
      DN_DX(3,2) = x10 * y20 - y10 * x20;

      DN_DX /= detJ;

      Volume = detJ*0.1666666666666666666667;
    }

    static double ComputeSubTetraVolumeAndCenter(const bounded_matrix<double, 3, 8 > & aux_coordinates, array_1d<double, 3 > & center_position, const int i0, const int i1, const int i2, const int i3)
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




  };

} // namespace Kratos.

#endif // KRATOS_ENRICHMENT_UTILITIES_INCLUDED  defined
