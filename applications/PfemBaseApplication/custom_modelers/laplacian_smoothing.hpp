//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_LAPLACIAN_SMOOTHING_H_INCLUDED)
#define  KRATOS_LAPLACIAN_SMOOTHING_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <stdlib.h>

#include <boost/timer.hpp>

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_flags.h"
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"
#include "custom_utilities/mesh_data_transfer_utilities.hpp"

#include "pfem_base_application_variables.h"

#ifdef   SINGLE
#define  REAL float
#else    // not SINGLE
#define  REAL double
#endif   // not SINGLE


#if !defined(KRATOS_TRIANGLE_EXTERNAL_H_INCLUDED)
#define  KRATOS_TRIANGLE_EXTERNAL_H_INCLUDED
#include "triangle.h"
#endif

//VARIABLES used:
//Data:     NEIGHBOUR_NODES
//StepData: DISPLACEMENT, CONTACT_FORCE, NORMAL, OFFSET
//Flags:    (checked) BOUNDARY, TO_ERASE, VISITED 
//          (set) 
//          (modified) VISITED 
//          (reset)

namespace Kratos
{
  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{

  ///@}
  ///@name  Enum's
  ///@{

  ///@}
  ///@name  Functions
  ///@{

  ///@}
  ///@name Kratos Classes
  ///@{

  /// Short class definition.
  /** Applies a recollocation of the nodes improving the mesh shape
   *  variables are interpolated to the new positions
   */
  class LaplacianSmoothing
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of data transfer
    KRATOS_CLASS_POINTER_DEFINITION( LaplacianSmoothing );
 
    typedef ModelPart::PropertiesType                                PropertiesType;
    typedef ModelPart::MeshType                                            MeshType;
    typedef ModelPart::ElementsContainerType                  ElementsContainerType;
    typedef ModelPart::NodesContainerType                        NodesContainerType;
    typedef ModelPart::MeshType::GeometryType::PointsArrayType      PointsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LaplacianSmoothing(ModelPart& rModelPart)
    {
      for (unsigned int i = 0; i < rModelPart.NumberOfMeshes(); i++)
	mPreviousMeshes.push_back(boost::make_shared<MeshType>((rModelPart.GetMesh(i)).Clone())); //in fact do not clones

    } //

    /// Destructor.
    virtual ~LaplacianSmoothing() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    //*******************************************************************************************
    //*******************************************************************************************

    void ApplyMeshSmoothing(ModelPart& rModelPart,
			    Flags Options,
			    ModelPart::IndexType MeshId=0)
    {

      //if(Meshes[MeshId].Is( MeshType::LAPLACIAN_SMOOTHING )){
      IncrementalSmoothing(rModelPart,MeshId);
      //}
	
    }

    //*******************************************************************************************
    //*******************************************************************************************


    void IncrementalSmoothing(ModelPart& rModelPart,
			      ModelPart::IndexType MeshId=0)
    {

      KRATOS_TRY
		    
      KRATOS_WATCH( "Start Laplacian Smoothing INCREMENTAL" )

      //defintions for spatial search
      typedef Node<3>                                  PointType;
      typedef Node<3>::Pointer                  PointPointerType;
      typedef std::vector<PointPointerType>          PointVector;
      typedef PointVector::iterator                PointIterator;
      typedef std::vector<double>                 DistanceVector;
      typedef std::vector<double>::iterator     DistanceIterator;

      typedef Bucket<3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > BucketType;
      typedef Tree< KDTreePartition<BucketType> >     KdtreeType; //Kdtree
      //defintions for spatial search
  
      double convergence_tol =0.001;
      double smoothing_factor=0.1;
      double smoothing_iters =3;
      double iters=0;

      bool simple = true; //weight = 1;

      NodesContainerType& rNodes = rModelPart.Nodes(MeshId);
      std::vector<array_1d<double,3> > initial_nodes_position;
      initial_nodes_position.resize(rNodes.size()+1);
      initial_nodes_position[0].clear();

      std::vector<int> nodes_ids;
      nodes_ids.resize(rModelPart.NumberOfNodes()+1); //mesh 0
      std::fill( nodes_ids.begin(), nodes_ids.end(), 0 );

      //initial loop to set the position to MESH_DISPLACEMENT
      rNodes.Sort();
      rNodes.Unique();
	  
      //std::cout<<" Nodes Size "<<rNodes.size()+1<<std::endl;

      unsigned int id=1;
      for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
	{
	  //point position
	  //std::cout<<" Initial Position ["<<in->Id()<<"]: ("<<in->X()<<", "<<in->Y()<<", "<<in->Z()<<") "<<std::endl;
		
	  nodes_ids[in->Id()] = id;
	  initial_nodes_position[id].clear();
	  initial_nodes_position[id][0] = in->X();	
	  initial_nodes_position[id][1] = in->Y();
	  initial_nodes_position[id][2] = in->Z();		
		
	  id++;
	}

	  
      bool converged=false;
      double MaxLength=0;
      double NewMaxLength=0;


      while ( iters<smoothing_iters && converged==false ){ 

	//std::cout<<" Iter "<< iters <<std::endl;

	array_1d<double,3> P;
	array_1d<double,3> Q;//neighbour position
	array_1d<double,3> D;
	    	    
	    
	double TotalWeight = 0;
	double Weight = 0;
	array_1d<double,3> TotalDistance;


	//convergence variables
	double Length = 0;
	MaxLength     = NewMaxLength;
	NewMaxLength  = 0;

	for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
	  {
	    if(in->IsNot(BOUNDARY))// && in->IsNot(STRUCTURE) )
	      {
		WeakPointerVector<Node<3> >& rNeighbours = in->GetValue(NEIGHBOUR_NODES);

		unsigned int NumberOfNeighbours = rNeighbours.size();
	      
		TotalDistance.clear();
		TotalWeight = 0;
		Weight = 0;

		//point position
		P[0] = in->X();
		P[1] = in->Y();
		P[2] = in->Z();
		     
		    
		//std::cout<<" Initial Position: "<<P<<std::endl;
		Length = 0;

		for(unsigned int i = 0; i < NumberOfNeighbours; i++)
		  {
		    //neighbour position
		    Q[0] = rNeighbours[i].X();
		    Q[1] = rNeighbours[i].Y();
		    Q[2] = rNeighbours[i].Z();
		    	
		    D = P-Q;
			
		    Length =sqrt(D[0]*D[0]+D[1]*D[1]+D[2]*D[2]);
			  

		    if( simple ){

		      Weight = 1;

		    }
		    else{
			  		
		      if(Length !=0)
			Weight = ( 1.0/Length );
		      else
			Weight = 0;
		    }

		    if(NewMaxLength<Length)
		      NewMaxLength = Length;
			
		    TotalDistance += (Weight*(Q-P)) ;  		  
		    TotalWeight   += Weight ;
		  
		  }
		    
	     
		if(TotalWeight!=0)
		  D = ( smoothing_factor / TotalWeight ) * TotalDistance;
		else
		  D.clear();
		    

		P += D;

		in->X() = P[0];
		in->Y() = P[1];
		in->Z() = P[2];


		//std::cout<<" Performed Position ["<<in->Id()<<":"<<iters<<"]: ("<<P[0]-D[0]<<", "<<P[1]-D[1]<<", "<<P[2]-D[2]<<") to ("<<in->X()<<", "<<in->Y()<<", "<<in->Z()<<") -> Total Weight: "<<TotalWeight<<std::endl;

	      }
	  }
	  

	if( (NewMaxLength-MaxLength)/NewMaxLength < convergence_tol ){
	  converged = true;
	  if( GetEchoLevel() > 0 )
	    std::cout<<" Laplacian smoothing convergence achieved "<<std::endl;
	}


	iters++;

	      
	  
      }

      if(iters==smoothing_iters && !converged)
	std::cout<<"   WARNING: Laplacian smoothing convergence NOT achieved "<<std::endl;

      bool transfer=true; //transfer active or inactive

      if(transfer==true){

	//create the list of the nodes to be check during the search
	PointVector list_of_nodes;	  
	list_of_nodes.reserve(rModelPart.NumberOfNodes(MeshId));

	for(NodesContainerType::iterator i_node = rModelPart.NodesBegin(MeshId) ; i_node != rModelPart.NodesEnd(MeshId) ; i_node++)
	  {

	    (list_of_nodes).push_back(*(i_node.base()));
	  }


	//Find out where the new nodes belong to:
	array_1d<double,3> N;

	std::vector<VariablesListDataValueContainer> VariablesListVector (list_of_nodes.size()+1);
	std::vector<int>                 UniquePosition (list_of_nodes.size()+1);
	std::fill( UniquePosition.begin(), UniquePosition.end(), 0 );

	unsigned int   MaximumNumberOfResults = list_of_nodes.size();
	PointVector    Results              (MaximumNumberOfResults);
	DistanceVector ResultsDistances     (MaximumNumberOfResults);

	unsigned int   bucket_size = 20;  
	KdtreeType  nodes_tree(list_of_nodes.begin(),list_of_nodes.end(),bucket_size);
	Node<3> work_point(0,0.0,0.0,0.0);    

	VariablesList& variables_list = rModelPart.GetNodalSolutionStepVariablesList();
	//unsigned int  step_data_size = rModelPart.GetNodalSolutionStepDataSize();

	//find the center and "radius" of the element
	double xc = 0;
	double yc = 0;
	double zc = 0;
	double radius = 0;

	for(ElementsContainerType::iterator ie = rModelPart.ElementsBegin(MeshId); ie != rModelPart.ElementsEnd(MeshId); ie++)
	  {

	    PointsArrayType& vertices = ie->GetGeometry().Points();
		    
	    //Idea, guardar coordenadas iniciales y desplazamiento total *** y corregir las coordenadas de los vertices a sus originales !! ---- eso: La corrección mediante variable auxiliar MESH_DISPLACEMENT.

	    array_1d<double,3>& V1 = initial_nodes_position[nodes_ids[vertices[0].Id()]];
	    array_1d<double,3>& V2 = initial_nodes_position[nodes_ids[vertices[1].Id()]];
	    array_1d<double,3>& V3 = initial_nodes_position[nodes_ids[vertices[2].Id()]];

	    // std::cout<<" V1 "<<V1<<" =  ("<<vertices[0].X()<<", "<<vertices[0].Y()<<",0 )"<<std::endl;
	    // std::cout<<" V2 "<<V2<<" =  ("<<vertices[1].X()<<", "<<vertices[1].Y()<<",0 )"<<std::endl;
	    // std::cout<<" V3 "<<V3<<" =  ("<<vertices[2].X()<<", "<<vertices[2].Y()<<",0 )"<<std::endl;

	    CalculateCenterAndSearchRadius( V1[0], V1[1],
					    V2[0], V2[1],
					    V3[0], V3[1],
					    xc,yc,radius);


	    //find all of the new nodes within the radius
	    work_point.X() = xc;
	    work_point.Y() = yc;
	    work_point.Z() = zc;
	      
	    // std::cout<<" Work Point  : "<<work_point<<std::endl;

	    int number_of_points_in_radius = nodes_tree.SearchInRadius (work_point, radius*1.01, Results.begin(), ResultsDistances.begin(),  MaximumNumberOfResults);


	    //check if inside and eventually interpolate
	    for( PointIterator it_found = Results.begin(); it_found != Results.begin() + number_of_points_in_radius; it_found++)
	      {
		//if((*it_found)->IsNot(STRUCTURE)){
		bool is_inside = false;
		is_inside = CalculatePosition( V1[0], V1[1],
					       V2[0], V2[1],
					       V3[0], V3[1],
					       (*it_found)->X(),(*it_found)->Y(),N);


		if(is_inside == true)
		  {
		    if(UniquePosition [nodes_ids[(*it_found)->Id()]] == 0){
		      UniquePosition [nodes_ids[(*it_found)->Id()]] = 1;

		      double alpha = 0.25; //[0,1] //smoothing level of the interpolation
		      MeshDataTransferUtilities DataTransferUtilities;
		      VariablesListVector  [nodes_ids[(*it_found)->Id()]] = DataTransferUtilities.InterpolateVariables( ie->GetGeometry(), N, variables_list, (*it_found), alpha );


		    }
		    else{
		      std::cout<<" WARNING LS: The node is relocated in a new element:: Something is Wrong "<<std::endl;
		      //std::cout<<" OldN "<<ShapeFunctions[nodes_ids[(*it_found)->Id()]]<<" NewN "<<N<<std::endl;
			      
		    }

		  }
		//}
	      }
	  }
	
	rModelPart.Nodes(MeshId).Sort();

	//*******************************************************************
	//CREATE NEW NODES:

	if( GetEchoLevel() > 0 )
	  std::cout<<" Create New Nodes "<<std::endl;

	int id=0;
	for(NodesContainerType::iterator i_node = rModelPart.NodesBegin(MeshId) ; i_node != rModelPart.NodesEnd(MeshId) ; i_node++)
	  {
	    if(i_node->IsNot(BOUNDARY)){
	      //recover the original position of the node
	      id = i_node->Id();
		  
	      i_node->SolutionStepData() = VariablesListVector[nodes_ids[id]];
		  
	      const array_1d<double,3>& disp = i_node->FastGetSolutionStepValue(DISPLACEMENT);
	      i_node->X0() = i_node->X() - disp[0];
	      i_node->Y0() = i_node->Y() - disp[1];
	      i_node->Z0() = i_node->Z() - disp[2];
	    }
		
	  }
	           
      }
      else{

	for(NodesContainerType::iterator i_node = rModelPart.NodesBegin(MeshId) ; i_node != rModelPart.NodesEnd(MeshId) ; i_node++)
	  {
		
	    //recover the original position of the node
	    const array_1d<double,3>& disp = i_node->FastGetSolutionStepValue(DISPLACEMENT);
	    i_node->X0() = i_node->X() - disp[0];
	    i_node->Y0() = i_node->Y() - disp[1];
	    i_node->Z0() = i_node->Z() - disp[2];

	  }
	    
      }
  
      KRATOS_WATCH( "Finished Laplacian Smoothing" )
    
      KRATOS_CATCH( "" )

   }
	
 
	
    void ApplyMeshSmoothing(ModelPart& rModelPart,
			    std::vector<int> & PreservedElements,
			    struct triangulateio &out,
			    std::vector<Geometry<Node<3> > >& list_of_element_vertices,
			    ModelPart::IndexType MeshId=0)
    {
      
      //defintions for spatial search
      typedef Node<3>                                  PointType;
      typedef Node<3>::Pointer                  PointPointerType;
      typedef std::vector<PointPointerType>          PointVector;
      typedef PointVector::iterator                PointIterator;
      typedef std::vector<double>                 DistanceVector;
      typedef std::vector<double>::iterator     DistanceIterator;
	
      typedef Bucket<3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > BucketType;
      typedef Tree< KDTreePartition<BucketType> >     KdtreeType; //Kdtree
      //defintions for spatial search

	
      //*******************************************************************
      //NEIGHBOUR NODES:
      std::vector<std::vector<int> >  list_of_neighbor_nodes= SetNeighborNodes(PreservedElements,out);
		
      //*******************************************************************
      //MOVE NODES: LAPLACIAN SMOOTHING:
	 
      double convergence_tol =0.001;
      double smoothing_factor=0.1;
      double smoothing_iters =3;
      double iters=0;

      bool simple = true; //weight = 1;

      NodesContainerType& rNodes = rModelPart.Nodes(MeshId);
      std::vector<array_1d<double,3> > initial_nodes_position;
      initial_nodes_position.resize(rNodes.size()+1);
      initial_nodes_position[0].clear();

      std::vector<int> nodes_ids;
      nodes_ids.resize(rModelPart.NumberOfNodes()+1); //mesh 0
      std::fill( nodes_ids.begin(), nodes_ids.end(), 0 );

		  
      //initial loop to set the position to MESH_DISPLACEMENT
      rNodes.Sort();
	    
      //create the list of original nodes position
      int dimension = 2;
      for(int in = 0; in<out.numberofpoints; in++)
	{
	  //point position
	  initial_nodes_position[in+1].clear();
	  for(int pn=0; pn<dimension; pn++)
	    {
	      initial_nodes_position[in+1][pn] = out.pointlist[in*dimension+pn];	
	    }

	}

	  
      bool converged=false;
      double MaxLength=0;
      double NewMaxLength=0;


      while ( iters<smoothing_iters && converged==false ){ 

	//std::cout<<" Iter "<< iters <<std::endl;

	array_1d<double,3> P;
	array_1d<double,3> Q;//neighbour position
	array_1d<double,3> D;
	    	    
	    
	double TotalWeight = 0;
	double Weight = 0;
	array_1d<double,3> TotalDistance;


	//convergence variables
	double Length = 0;
	MaxLength     = NewMaxLength;
	NewMaxLength  = 0;
	   
	for(int in = 0; in<out.numberofpoints; in++)
	  {
	    if(rNodes[in+1].IsNot(BOUNDARY) && rNodes[in+1].IsNot(TO_ERASE))// && rNodes[in+1].IsNot(STRUCTURE) )
	      {
		unsigned int NumberOfNeighbours = list_of_neighbor_nodes[in+1].size();
	    	      
		TotalDistance.clear();
		TotalWeight = 0;
		Weight = 0;

		//point position
		P[0] = out.pointlist[in*dimension];
		P[1] = out.pointlist[in*dimension+1];
		P[2] = 0;
		     
		    
		//std::cout<<" Initial Position: "<<P<<std::endl;
		Length = 0;

		for(unsigned int i = 0; i < NumberOfNeighbours; i++)
		  {
		    //neighbour position
		    Q[0] = out.pointlist[(list_of_neighbor_nodes[in+1][i]-1)*dimension];
		    Q[1] = out.pointlist[(list_of_neighbor_nodes[in+1][i]-1)*dimension+1];
		    Q[2] = 0;
		    	
		    D = P-Q;
			
		    Length =sqrt(D[0]*D[0]+D[1]*D[1]+D[2]*D[2]);
			  

		    if( simple ){

		      Weight = 1;

		    }
		    else{
			  		
		      if(Length !=0)
			Weight = ( 1.0/Length );
		      else
			Weight = 0;
		    }

		    if(NewMaxLength<Length)
		      NewMaxLength = Length;
			
		    TotalDistance += (Weight*(Q-P)) ;  		  
		    TotalWeight   += Weight ;
		  
		  }
		    
	     
		if(TotalWeight!=0)
		  D = ( smoothing_factor / TotalWeight ) * TotalDistance;
		else
		  D.clear();
		    

		P += D;

		out.pointlist[in*dimension]   = P[0];
		out.pointlist[in*dimension+1] = P[1];
		    

	      }
	  }
	  

	if( (NewMaxLength-MaxLength)/NewMaxLength < convergence_tol ){
	  converged = true;
	  if( GetEchoLevel() > 0 )
	    std::cout<<"   Laplacian smoothing convergence achieved "<<std::endl;
	}


	iters++;

      }

      if(iters==smoothing_iters && !converged)
	std::cout<<"   WARNING: Laplacian smoothing convergence NOT achieved "<<std::endl;



      //*******************************************************************
      //MOVE NODES: BOUNDARY SMOOTHING
	 
      SetBoundarySmoothing (rModelPart, rNodes, PreservedElements, out, MeshId);


      //*******************************************************************
      //MOVE NODES: BOUNDARY PROJECTION
	 
      //SetInsideProjection (rModelPart, out, list_of_neighbor_nodes, MeshId);


      //*******************************************************************
      //TRANSFER VARIABLES TO NEW NODES POSITION:


      bool transfer=true; //transfer active or inactive

      if(transfer==true){

	//create the list of the nodes to be check during the search (new positions after the laplacian smoothing)
	PointVector list_of_nodes;	  
	list_of_nodes.reserve(rModelPart.NumberOfNodes(MeshId));
	    
	  
	for(NodesContainerType::iterator i_node = rModelPart.NodesBegin(MeshId) ; i_node != rModelPart.NodesEnd(MeshId) ; i_node++)
	  {
	    if( i_node->IsNot(TO_ERASE) ){
	      //set new positions:
	      int id =i_node->Id()-1;
	      i_node->X() = out.pointlist[id*dimension];
	      i_node->Y() = out.pointlist[id*dimension+1];
	      i_node->Z() = 0;
	      
	      (list_of_nodes).push_back(*(i_node.base()));
	    }
	    else {
	      std::cout <<" LLM:PSEUDOERROR node to erase : " << i_node->Id() << std::endl;
	    }
	  }
	    
	//Find out where the new nodes belong to:
	array_1d<double,3> N;
	array_1d<int,3>    ID;
	    
	std::vector<VariablesListDataValueContainer> VariablesListVector(rModelPart.NumberOfNodes(MeshId)+1);
	std::vector<int>                 UniquePosition (rModelPart.NumberOfNodes(MeshId)+1);
	std::fill( UniquePosition.begin(), UniquePosition.end(), 0 );
	    
	unsigned int   MaximumNumberOfResults = list_of_nodes.size();
	PointVector    Results              (MaximumNumberOfResults);
	DistanceVector ResultsDistances     (MaximumNumberOfResults);

	unsigned int   bucket_size = 20;
	KdtreeType     nodes_tree(list_of_nodes.begin(),list_of_nodes.end(),bucket_size);
	ModelPart::NodesContainerType::iterator nodes_begin = rModelPart.NodesBegin(MeshId);

	//unsigned int  step_data_size = rModelPart.GetNodalSolutionStepDataSize();
	VariablesList&  variables_list = rModelPart.GetNodalSolutionStepVariablesList();

	//find the center and "radius" of the element
	double xc = 0;
	double yc = 0;
	double zc = 0;
	double radius = 0;
	Node<3> work_point (0,0.0,0.0,0.0);

	//geometry
	double x1,x2,x3,y1,y2,y3;

	for(int el = 0; el<out.numberoftriangles; el++)
	  {
	    if(PreservedElements[el])
	      {
		    
		ID[0] = out.trianglelist[el*3+0];
		ID[1] = out.trianglelist[el*3+1];
		ID[2] = out.trianglelist[el*3+2];


		x1 = initial_nodes_position[ID[0]][0];
		y1 = initial_nodes_position[ID[0]][1];

		x2 = initial_nodes_position[ID[1]][0];
		y2 = initial_nodes_position[ID[1]][1];

		x3 = initial_nodes_position[ID[2]][0];
		y3 = initial_nodes_position[ID[2]][1];
		    

		CalculateCenterAndSearchRadius( x1, y1,
						x2, y2,
						x3, y3,
						xc, yc, radius);


		//find all of the new nodes within the radius
		work_point.X() = xc;
		work_point.Y() = yc;
		work_point.Z() = zc;
		    
		//std::cout<<" Work Point  : "<<work_point<<std::endl;

		int number_of_points_in_radius = nodes_tree.SearchInRadius (work_point, radius*1.01, Results.begin(), ResultsDistances.begin(),  MaximumNumberOfResults);

		//std::cout<<"[ID:"<<el<<"]: NumberOfPointsInRadius "<<number_of_points_in_radius<<std::endl;

		//check if inside and eventually interpolate
		for( PointIterator it_found = Results.begin(); it_found != Results.begin() + number_of_points_in_radius; it_found++)
		  {

		    if((*it_found)->IsNot(TO_ERASE)){

		      //std::cout<<" Found ID "<<(*it_found)->Id()<<std::endl;
		      bool is_inside = false;
		      is_inside = CalculatePosition( x1, y1,
						     x2, y2,
						     x3, y3,
						     (*it_found)->X(),(*it_found)->Y(),N);


		      PointsArrayType  PointsArray;
		      PointsArray.push_back( *( (nodes_begin +  ID[0]-1).base() ) ); 
		      PointsArray.push_back( *( (nodes_begin +  ID[1]-1).base() ) ); 
		      PointsArray.push_back( *( (nodes_begin +  ID[2]-1).base() ) ); 

		      Geometry<Node<3> > geom( PointsArray );

		      if(is_inside == true)
			{
			  //std::cout<<"  Node interpolation: "<<(*it_found)->Id()<<" VariablesList size "<<VariablesListVector.size()<<std::endl;
			  if(UniquePosition [(*it_found)->Id()] == 0){
			    
			    UniquePosition [(*it_found)->Id()] = 1;

			    double alpha = 0.25; //[0,1] //smoothing level of the interpolation
			    MeshDataTransferUtilities DataTransferUtilities;
			    VariablesListVector  [(*it_found)->Id()] = DataTransferUtilities.InterpolateVariables( geom, N, variables_list, (*it_found), alpha );
			    
			  }
			  // else{
			  //   std::cout<<" The node "<<(*it_found)->Id()<<" is relocated in a new element:: Something is Wrong "<<std::endl;
			  //   std::cout<<" OldN "<<ShapeFunctions [(*it_found)->Id()]<<" NewN "<<N<<std::endl;
			  // }

			}
		    }
		  }
	      }
	  }
	    
	//the search moves the nodes order using its PreIds
	rModelPart.Nodes(MeshId).Sort();

	//*******************************************************************
	//CREATE NEW NODE INFORMATION:

	int id=0;
	for(NodesContainerType::iterator i_node = rModelPart.NodesBegin(MeshId) ; i_node != rModelPart.NodesEnd(MeshId) ; i_node++)
	  {
	    //std::cout<<" ID "<<i_node->Id()<<std::endl;
	    if(i_node->IsNot(BOUNDARY) && i_node->IsNot(TO_ERASE)){
	      //recover the original position of the node
	      id = i_node->Id();

	      // double PressurePrev = (i_node)->FastGetSolutionStepValue(PRESSURE); 

	      //i_node->SolutionStepData() = VariablesListVector[id];

	      if ( VariablesListVector[id].Has(DISPLACEMENT) == false)
		{
		  std::cout << " PSEUDOERROR: IN this line, there is a node that does not have displacement" << std::endl;
		  std::cout << " Laplacian. ThisNode new information does not have displacement " << i_node->Id() << std::endl;
		  std::cout << " THIS IS BECAUSE THE NODE is out of the DOMAIN and the interpolation is wrong" << std::endl;
		  std::cout << "    X: " << i_node->X() << " Y: " << i_node->Y() << std::endl;
		}
	      else {
		i_node->SolutionStepData() = VariablesListVector[id];
	      }
		  
	      // double PressurePost = (i_node)->FastGetSolutionStepValue(PRESSURE);
	      // std::cout<<" PRESSURE PREV "<<PressurePrev<<" PRESSURE POST "<<PressurePost<<std::endl;

	      const array_1d<double,3>& disp = i_node->FastGetSolutionStepValue(DISPLACEMENT);

	      i_node->X0() = i_node->X() - disp[0];
	      i_node->Y0() = i_node->Y() - disp[1];
	      i_node->Z0() = i_node->Z() - disp[2];

	    }
	    // Set the position of boundary laplacian 
	    else if ( i_node->Is(BOUNDARY) && i_node->IsNot(TO_ERASE) && i_node->Is(VISITED) )
	      {
		i_node->Set(VISITED, false);

		//recover the original position of the node
		id = i_node->Id();

		// double PressurePrev = (i_node)->FastGetSolutionStepValue(PRESSURE); 

		//i_node->SolutionStepData() = VariablesListVector[id];

		if ( VariablesListVector[id].Has(DISPLACEMENT) == false)
		  {
		    std::cout << " OUT::PSEUDOERROR: IN this line, there is a node that does not have displacement" << std::endl;
		    std::cout << " Laplacian. ThisNode new information does not have displacement " << i_node->Id() << std::endl;
		    std::cout << " THIS IS BECAUSE THE NODE is out of the DOMAIN and the interpolation is wrong" << std::endl;
		    std::cout << "    X: " << i_node->X() << " Y: " << i_node->Y() << std::endl;
		  }
		else {
		  i_node->SolutionStepData() = VariablesListVector[id];
		}

		// double PressurePost = (i_node)->FastGetSolutionStepValue(PRESSURE);
		// std::cout<<" PRESSURE PREV "<<PressurePrev<<" PRESSURE POST "<<PressurePost<<std::endl;

		if ( i_node->SolutionStepsDataHas(DISPLACEMENT) == false)
		  {
		    std::cout << " AFTER WIERD " << std::endl;
		    std::cout << " Laplacian. ThisNode Does not have displacemenet " << i_node->Id() << std::endl;
		    std::cout << "    X: " << i_node->X() << " Y: " << i_node->Y() << std::endl;
		  }

		const array_1d<double,3>& disp = i_node->FastGetSolutionStepValue(DISPLACEMENT);

		bool MoveFixedNodes = false; 
		if (MoveFixedNodes)
		  {
		    i_node->X0() = i_node->X() - disp[0];
		    i_node->Y0() = i_node->Y() - disp[1];
		    i_node->Z0() = i_node->Z() - disp[2];
		  }
		else {
		  if ( i_node->pGetDof(DISPLACEMENT_X)->IsFixed() == false) {
		    i_node->X0() = i_node->X() - disp[0];
		  }
		  if ( i_node->pGetDof(DISPLACEMENT_Y)->IsFixed() == false) {
		    i_node->Y0() = i_node->Y() - disp[1];
		  }

		  if ( i_node->pGetDof(DISPLACEMENT_Z)->IsFixed() == false) {
		    i_node->Z0() = i_node->Z() - disp[2];
		  }
		}

	      }
	  }

	//New nodes created, element vertices must be set again: NOT NOW
	// std::vector<Geometry<Node<3> > > EmptyList
	// list_of_element_vertices.swap(EmptyList); //destroy all elements
	// list_of_element_vertices.clear(); 

	// for(int el = 0; el<out.numberoftriangles; el++)
	//   {
	// 	if(PreservedElements[el])
	// 	  {
	// 	    Geometry<Node<3> > vertices;
	// 	    std::vector<int >  neighbours (3);
	
	// 	    for(int pn=0; pn<3; pn++)
	// 	      {
	// 		vertices.push_back(*(nodes_begin + out.trianglelist[el*3+pn]-1).base());
	// 		//vertices.push_back(rModelPart.pGetNode(out.trianglelist[el*3+pn],MeshId));
	// 	      }
		  
	// 	    list_of_element_vertices.push_back( vertices );
	// 	  }
	//   }

	    	
      }
      else{

	for(NodesContainerType::iterator i_node = rModelPart.NodesBegin(MeshId) ; i_node != rModelPart.NodesEnd(MeshId) ; i_node++)
	  {
		
	    //set the new position to the node
	    int id = i_node->Id()-1;
	    i_node->X() = out.pointlist[id*dimension];
	    i_node->Y() = out.pointlist[id*dimension+1];
	    i_node->Z() = 0;

	    //recover the original position of the node
	    array_1d<double,3>& disp = i_node->FastGetSolutionStepValue(DISPLACEMENT);
	    if(norm_2(disp)>0){
	      i_node->X0() = i_node->X() - disp[0];
	      i_node->Y0() = i_node->Y() - disp[1];
	      i_node->Z0() = i_node->Z() - disp[2];
	    }
       //Set the position of boundary laplacian (Reset the flag)
       if ( i_node->Is(BOUNDARY) && i_node->IsNot(TO_ERASE) && i_node->Is(VISITED) )
       {
            i_node->Set(VISITED, false); //LMV.
       }

	  }
	    
      }

    }


    /**
     * level of echo for the mesh smoothing
     */
    virtual void SetEchoLevel(int Level)
    {
      mEchoLevel = Level;
    }
      
    int GetEchoLevel()
    {
      return mEchoLevel;
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
      return "";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


    ///@}
    ///@name Friends
    ///@{


    ///@}

  protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

  private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart::MeshesContainerType mPreviousMeshes;

    int mEchoLevel;

    ///@}
    ///@name Private Operators
    ///@{

    /// Assignment operator.
    LaplacianSmoothing& operator=(LaplacianSmoothing const& rOther);

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Unaccessible methods
    ///@{


    std::vector<std::vector<int> > SetNeighborNodes (std::vector<int> & PreservedElements,struct triangulateio &out)
    {
		 
      std::vector<int> empty_vector(0);
      std::vector<std::vector<int> >  list_of_neighbor_nodes(out.numberofpoints+1);
      std::fill( list_of_neighbor_nodes.begin(), list_of_neighbor_nodes.end(), empty_vector );

      bool neighb_set  = false;
      int  neighb_size = 0;
      for(int el = 0; el<out.numberoftriangles; el++)
	{
	  if(PreservedElements[el])
	    {
	      //a) Create list of node neighbors (list_of_neighbor_nodes)
	      for(int ipn=0; ipn<3; ipn++)
		{
		  
		  for(int jpn=0; jpn<3; jpn++)
		    {
		      if(ipn!=jpn){
			//add unique node neighbor
			neighb_size = list_of_neighbor_nodes[out.trianglelist[el*3+ipn]].size();
			neighb_set = false;
			for(int npn=0; npn<neighb_size; npn++)
			  {
			    if( list_of_neighbor_nodes[out.trianglelist[el*3+ipn]][npn]==(out.trianglelist[el*3+jpn]) ){
			      neighb_set=true;
			    }
			  }
			if(neighb_set==false){
			  list_of_neighbor_nodes[out.trianglelist[el*3+ipn]].push_back(out.trianglelist[el*3+jpn]);
			}
		      }					  
		    }
		}
	    }
	}

      return list_of_neighbor_nodes;
    }


    /**
     * boundary smoothing
     */
 
    std::vector<std::vector<int> > SetBoundaryNeighborNodes (ModelPart& rModelPart, 
							     std::vector<int> & PreservedElements, 
							     struct triangulateio &out, 
							     ModelPart::IndexType MeshId=0)
    {
		 
      std::vector<int> empty_vector(0);
      std::vector<std::vector<int> >  list_of_neighbor_nodes(out.numberofpoints+1);
      std::fill( list_of_neighbor_nodes.begin(), list_of_neighbor_nodes.end(), empty_vector );

      bool neighb_set  = false;
      int  neighb_size = 0;
      for(int el = 0; el<out.numberoftriangles; el++)
	{
	  if(PreservedElements[el])
	    {
	      //a) Create list of node neighbors (list_of_neighbor_nodes)
	      for(int ipn=0; ipn<3; ipn++)
		{
		  int nodei_id = out.trianglelist[el*3+ipn];
		  if(rModelPart.Nodes(MeshId)[nodei_id].Is(BOUNDARY)){

		    for(int jpn=0; jpn<3; jpn++)
		      {
			int nodej_id = out.trianglelist[el*3+jpn];
			if(ipn!=jpn && rModelPart.Nodes(MeshId)[nodej_id].Is(BOUNDARY) ){

			  //add unique node neighbor
			  neighb_size = list_of_neighbor_nodes[out.trianglelist[el*3+ipn]].size();
			  neighb_set = false;
			  for(int npn=0; npn<neighb_size; npn++)
			    {
			      if( list_of_neighbor_nodes[out.trianglelist[el*3+ipn]][npn]==(out.trianglelist[el*3+jpn]) ){
				neighb_set=true;
			      }
			    }
			  if(neighb_set==false){
			    list_of_neighbor_nodes[out.trianglelist[el*3+ipn]].push_back(out.trianglelist[el*3+jpn]);
			  }
			 			  
			}
		      }
		  }
		}
	    }
	}

      return list_of_neighbor_nodes;
    }


    std::vector<double>  SetRanks (ModelPart& rModelPart,
				   struct triangulateio &out,
				   std::vector<std::vector<int> > & list_of_neighbor_nodes,
				   ModelPart::IndexType MeshId=0)
    {
	
      //set ranks
	
      std::vector<double> nodes_ranks;
      nodes_ranks.resize(rModelPart.NumberOfNodes(MeshId)+1); //mesh 0
      std::fill( nodes_ranks.begin(), nodes_ranks.end(), 0 );
	
	
      //initial assignation
      for(int in = 0; in<out.numberofpoints; in++)
	{
	  array_1d<double, 3 > & ContactForceNormal  = rModelPart.Nodes(MeshId)[in+1].FastGetSolutionStepValue(CONTACT_FORCE);
	  if(norm_2(ContactForceNormal)==0 && rModelPart.Nodes(MeshId)[in+1].IsNot(BOUNDARY)){
	    nodes_ranks[in+1]=5;
	  }	    
	  // else if( norm_2(ContactForceNormal)==0 && rModelPart.Nodes(MeshId)[in+1].Is(BOUNDARY)){
	  //   nodes_ranks[in+1]=1;
	  // }
	    
	}
	
	
      //RANK assignation: 
      double rang_assign = 0;
      double rang_top = 5; //3;
	
      while (rang_assign<rang_top){
	  
	for(int in = 0; in<out.numberofpoints; in++)
	  {
	    if(nodes_ranks[in+1]==rang_assign){
		
	      //Rank 0
	      unsigned int shared_node=1;
	      unsigned int NumberOfNeighbours = list_of_neighbor_nodes[in+1].size();	
	      for(unsigned int i = 0; i < NumberOfNeighbours; i++)
		{
		  shared_node = list_of_neighbor_nodes[in+1][i];
		    
		  if(nodes_ranks[shared_node]>rang_assign)
		    nodes_ranks[shared_node]=rang_assign+1;
		}
	    }
	      
	  }
	  
	rang_assign++;
      }
	
      return nodes_ranks;
	  
    }
      

    std::vector<double>  SetFirstLayer (ModelPart& rModelPart,
					struct triangulateio &out,
					std::vector<std::vector<int> > & list_of_neighbor_nodes,
					ModelPart::IndexType MeshId=0)
    {
	
      //set ranks
	
      std::vector<double> nodes_layer;
      nodes_layer.resize(rModelPart.NumberOfNodes(MeshId)+1); //mesh 0
      std::fill( nodes_layer.begin(), nodes_layer.end(), 0 );
	
	
      //initial assignation
      for(int in = 0; in<out.numberofpoints; in++)
	{
	  if(rModelPart.Nodes(MeshId)[in+1].IsNot(BOUNDARY)){
	    nodes_layer[in+1]=2;
	  }	    
 	    
	}
	
	
      //LAYER assignation: 
      double layer_assign = 0;
      double layer_top    = 1; 
	
      while (layer_assign<layer_top){
	  
	for(int in = 0; in<out.numberofpoints; in++)
	  {
	    if(nodes_layer[in+1]==layer_assign){
		
	      //Rank 0
	      unsigned int shared_node=1;
	      unsigned int NumberOfNeighbours = list_of_neighbor_nodes[in+1].size();	
	      for(unsigned int i = 0; i < NumberOfNeighbours; i++)
		{
		  shared_node = list_of_neighbor_nodes[in+1][i];
		    
		  if(nodes_layer[shared_node]>layer_assign)
		    nodes_layer[shared_node]=layer_assign+1;
		}
	    }
	      
	  }
	  
	layer_assign++;
      }
	
      return nodes_layer;
	  
    }
      



    void SetBoundarySmoothing(ModelPart& rModelPart,
			      NodesContainerType& rNodes,
			      std::vector<int> & PreservedElements,
			      struct triangulateio &out,
			      ModelPart::IndexType MeshId=0)
    {
      
      //defintions for spatial search
      typedef Node<3>                                  PointType;
      typedef Node<3>::Pointer                  PointPointerType;
      typedef std::vector<PointPointerType>          PointVector;
      typedef PointVector::iterator                PointIterator;
      typedef std::vector<double>                 DistanceVector;
      typedef std::vector<double>::iterator     DistanceIterator;
	
      typedef Bucket<3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > BucketType;
      typedef Tree< KDTreePartition<BucketType> >     KdtreeType; //Kdtree
      //defintions for spatial search

	
      //*******************************************************************
      //NEIGHBOUR NODES:
      std::vector<std::vector<int> >  list_of_neighbor_nodes= SetBoundaryNeighborNodes(rModelPart,PreservedElements,out,MeshId);
		
      //*******************************************************************
      //MOVE BOUNDARY NODES: LAPLACIAN SMOOTHING:
	 
      double convergence_tol =0.001;
      double smoothing_factor=0.1; //0.1
      double smoothing_iters =4; //3
      double iters=0;

      bool simple = true; //weight = 1;  
      bool converged=false;

      double MaxLength=0;
      double NewMaxLength=0;


      while ( iters<smoothing_iters && converged==false ){ 

	//std::cout<<" Iter "<< iters <<std::endl;

	array_1d<double,3> P;
	array_1d<double,3> Q;//neighbour position
	array_1d<double,3> D;
	    	    
	    
	double TotalWeight = 0;
	double Weight = 0;
	array_1d<double,3> TotalDistance;


	//convergence variables
	double Length = 0;
	MaxLength     = NewMaxLength;
	NewMaxLength  = 0;
	   
	int dimension = 2;
	for(int in = 0; in<out.numberofpoints; in++)
	  {
	    if(rNodes[in+1].Is(BOUNDARY) && rNodes[in+1].IsNot(TO_ERASE) && rNodes[in+1].Is(VISITED) )
	      {
		unsigned int NumberOfNeighbours = list_of_neighbor_nodes[in+1].size();
	    	      
		TotalDistance.clear();
		TotalWeight = 0;
		Weight = 0;

		//point position
		P[0] = out.pointlist[in*dimension];
		P[1] = out.pointlist[in*dimension+1];
		P[2] = 0;
		     
		    
		//std::cout<<" Initial Position: "<<P<<std::endl;
		Length = 0;

		for(unsigned int i = 0; i < NumberOfNeighbours; i++)
		  {
		    //neighbour position
		    Q[0] = out.pointlist[(list_of_neighbor_nodes[in+1][i]-1)*dimension];
		    Q[1] = out.pointlist[(list_of_neighbor_nodes[in+1][i]-1)*dimension+1];
		    Q[2] = 0;
		    	
		    D = P-Q;
			
		    Length =sqrt(D[0]*D[0]+D[1]*D[1]+D[2]*D[2]);
			  

		    if( simple ){

		      Weight = 1;

		    }
		    else{
			  		
		      if(Length !=0)
			Weight = ( 1.0/Length );
		      else
			Weight = 0;
		    }

		    if(NewMaxLength<Length)
		      NewMaxLength = Length;
			
		    TotalDistance += (Weight*(Q-P)) ;  		  
		    TotalWeight   += Weight ;
		  
		  }
		    
	     
		if(TotalWeight!=0)
		  D = ( smoothing_factor / TotalWeight ) * TotalDistance;
		else
		  D.clear();
		    

		P += D;


		//std::cout<<" SET BOUNDARY SMOOTHING "<<P<<" from ["<<out.pointlist[in*dimension]<<","<<out.pointlist[in*dimension+1]<<"]"<<std::endl;

		out.pointlist[in*dimension]   = P[0];
		out.pointlist[in*dimension+1] = P[1];
		    

	      }

	    //rNodes[in+1].Set(VISITED,false); //LMV: Reset the flag after interpolation. Indeed, if the flag is set, only one iteration takes place
	  }
	  

	if( (NewMaxLength-MaxLength)/NewMaxLength < convergence_tol ){
	  converged = true;
	  if( GetEchoLevel() > 0 )
	    std::cout<<"   Laplacian smoothing convergence achieved "<<std::endl;
	}


	iters++;

      }

      if(iters==smoothing_iters && !converged)
	std::cout<<"   WARNING: Boundary Laplacian smoothing convergence NOT achieved (iters:"<<iters<<")"<<std::endl;


    }



    //Note : to increase de robustness I propose to detect the layer nodes, via setting ranks to nodes 
    // then ensure that the layer nodes have a certain distance to the boundary. Bigger than the offset applied


    void SetInsideProjection (ModelPart& rModelPart,
			      struct triangulateio &out,
			      std::vector<std::vector<int> > & list_of_neighbor_nodes,
			      ModelPart::IndexType MeshId=0)
    {
	
	
      std::vector<double> nodes_ranks = SetRanks(rModelPart,out,list_of_neighbor_nodes,MeshId);

      std::vector<double> nodes_layer = SetFirstLayer(rModelPart,out,list_of_neighbor_nodes,MeshId);

      double movement_factor = 1.2;
      double contact_factor  = 2.0;
      const array_1d<double,3> ZeroPoint(3,0.0);
      std::vector<array_1d<double,3> > initial_nodes_distances (rModelPart.NumberOfNodes(MeshId)+1);
      std::fill( initial_nodes_distances.begin(), initial_nodes_distances.end(), ZeroPoint );

      for(int in = 0; in<out.numberofpoints; in++)
	{
	    
	  //if(nodes_ranks[in+1]<=1)
	  if(nodes_ranks[in+1]<1)
	    {
	      array_1d<double, 3 >  DeltaDisplacement  = rModelPart.Nodes(MeshId)[in+1].FastGetSolutionStepValue(DISPLACEMENT) - rModelPart.Nodes(MeshId)[in+1].FastGetSolutionStepValue(DISPLACEMENT,1);

	      array_1d<double, 3>&  Normal= rModelPart.Nodes(MeshId)[in+1].FastGetSolutionStepValue(NORMAL); //BOUNDARY_NORMAL must be set as nodal variable

	      double projection=inner_prod(DeltaDisplacement,Normal);

	      array_1d<double, 3 > & ContactForce = rModelPart.Nodes(MeshId)[in+1].FastGetSolutionStepValue(CONTACT_FORCE);
	      if(norm_2(ContactForce)!=0){
		initial_nodes_distances[in+1] = (-1)*(movement_factor*contact_factor)*fabs(projection)*Normal;
	      }
	      else{
		initial_nodes_distances[in+1] = (-1)*(movement_factor)*fabs(projection)*Normal;
	      }
		
	    }

	  //layer modification

	  array_1d<double,3> P;
	  array_1d<double,3> Q;//neighbour position
	  array_1d<double,3> D;
	  int dimension = 2;


	  if(nodes_layer[in+1]==1)
	    {

	      unsigned int NumberOfNeighbours = list_of_neighbor_nodes[in+1].size();
	    
	      //point position
	      P[0] = out.pointlist[in*dimension];
	      P[1] = out.pointlist[in*dimension+1];
	      P[2] = 0;
		  
	      unsigned int shared_node = 1;
	      for(unsigned int i = 0; i < NumberOfNeighbours; i++)
		{
		  shared_node = list_of_neighbor_nodes[in+1][i];

		  if(rModelPart.Nodes(MeshId)[shared_node].Is(BOUNDARY)){

		    //neighbour position
		    Q[0] = out.pointlist[(list_of_neighbor_nodes[in+1][i]-1)*dimension];
		    Q[1] = out.pointlist[(list_of_neighbor_nodes[in+1][i]-1)*dimension+1];
		    Q[2] = 0;
		    array_1d<double, 3>&  Normal= rModelPart.Nodes(MeshId)[shared_node].FastGetSolutionStepValue(NORMAL); //BOUNDARY_NORMAL must be set as nodal variable
			
		    D = Q-P;

		    double projection=inner_prod(D,Normal);

		    array_1d<double, 3>& Offset = rModelPart.Nodes(MeshId)[shared_node].FastGetSolutionStepValue(OFFSET);
		    double offset = norm_2(Offset);
			
		    double secure_offset_factor = 1.1;

		    if(projection<offset && offset!=0)
		      initial_nodes_distances[in+1] = (-1)*(secure_offset_factor)*(offset-projection)*Normal;
		  }
		      
		}
		  
	    }
	}


      double smoothing_factor=0.15;
      double smoothing_iters =10;

      double iters=0;  

      while ( iters<smoothing_iters ){ 

	//std::cout<<" Iter "<< iters <<std::endl;
	    
	double TotalWeight = 0;
	double Weight      = 0;
	array_1d<double,3> TotalDistance;
	array_1d<double,3> Distance;
	array_1d<double,3> OffsetDistance;

	//convergence variables	   
	for(int in = 0; in<out.numberofpoints; in++)
	  {
	      
	    TotalDistance = initial_nodes_distances[in+1];
	    OffsetDistance.clear();


	    if(nodes_ranks[in+1]>0)
	      {
		  
		if(nodes_layer[in+1]==1)
		  OffsetDistance = TotalDistance;
		    
		unsigned int NumberOfNeighbours = list_of_neighbor_nodes[in+1].size();
	    	      
		unsigned int shared_node=1;

		TotalWeight = 0;
		Weight      = 0;
		Distance.clear();
		    

		for(unsigned int i = 0; i < NumberOfNeighbours; i++)
		  {

		    //neighbour
		    shared_node    = list_of_neighbor_nodes[in+1][i];

		    Weight         = 1.0 / (nodes_ranks[shared_node]+1.0);
		    TotalWeight   += Weight ;	
		    Distance      += initial_nodes_distances[shared_node] * Weight;  	      	  
		  }
		    
	     
		if(TotalWeight!=0)
		  Distance *= ( 1.0 / TotalWeight );	    
		else
		  Distance = initial_nodes_distances[in+1];

		TotalDistance += smoothing_factor*(Distance-initial_nodes_distances[in+1]);
	      }


	    if( nodes_layer[in+1]==1 && norm_2(OffsetDistance)>norm_2(TotalDistance)+1e-10 ){
	      TotalDistance = OffsetDistance;
	      if( GetEchoLevel() > 0 )
		std::cout<<" Layer Correction "<<norm_2(OffsetDistance)<<" > "<<norm_2(TotalDistance)<<std::endl; 
	    }

	    initial_nodes_distances[in+1] =  TotalDistance;


	  }
	  

	iters++;

      }

      int dimension = 2;
      for(int in = 0; in<out.numberofpoints; in++)
	{
	    
	  // array_1d<double, 3>&  Projection= rModelPart.Nodes(MeshId)[in+1].FastGetSolutionStepValue(FORCE_EXTERNAL); //BOUNDARY_NORMAL must be set as nodal variable
	  // Projection = initial_nodes_distances[in+1];

	  if(nodes_ranks[in+1]>0)	      
	    {
	      //std::cout<<" Projection Set "<<initial_nodes_distances[in+1]<<std::endl;
	      out.pointlist[in*dimension]   += initial_nodes_distances[in+1][0];
	      out.pointlist[in*dimension+1] += initial_nodes_distances[in+1][1];
	    }

	    
	}



    }



    inline void CalculateCenterAndSearchRadius(const double x0, const double y0,
					       const double x1, const double y1,
					       const double x2, const double y2,
					       double& xc, double& yc, double& R)

    {
      xc = 0.3333333333333333333*(x0+x1+x2);
      yc = 0.3333333333333333333*(y0+y1+y2);

      double R1 = (xc-x0)*(xc-x0) + (yc-y0)*(yc-y0);
      double R2 = (xc-x1)*(xc-x1) + (yc-y1)*(yc-y1);
      double R3 = (xc-x2)*(xc-x2) + (yc-y2)*(yc-y2);

      R = R1;
      if(R2 > R) R = R2;
      if(R3 > R) R = R3;

      R = sqrt(R);
    }

    inline double CalculateVol(const double x0, const double y0,
			       const double x1, const double y1,
			       const double x2, const double y2)
    {
      return 0.5*( (x1-x0)*(y2-y0)- (y1-y0)*(x2-x0) );
    }

    inline bool CalculatePosition(const double x0, const double y0,
				  const double x1, const double y1,
				  const double x2, const double y2,
				  const double xc, const double yc,
				  array_1d<double,3>& N)
    {
      double area = CalculateVol(x0,y0,x1,y1,x2,y2);

      //std::cout<<" Area "<<area<<std::endl;
	    
      if(area < 1e-15)
	{
	  //KRATOS_THROW_ERROR( std::logic_error,"element with zero area found", "" );
	  std::cout<<" ERROR LS: element with zero area found: "<<area<<" position ("<<x0<<", "<<y0<<") ("<<x1<<", "<<y1<<") ("<<x2<<", "<<y2<<") "<<std::endl;
	}

      N[0] = CalculateVol(x1,y1,x2,y2,xc,yc)  / area;
      N[1] = CalculateVol(x2,y2,x0,y0,xc,yc)  / area;
      N[2] = CalculateVol(x0,y0,x1,y1,xc,yc)  / area;

      double tol = 1e-5;
      double upper_limit = 1.0+tol;
      double lower_limit = -tol;

      if(N[0] >= lower_limit && N[1] >= lower_limit && N[2] >= lower_limit && N[0] <= upper_limit && N[1] <= upper_limit && N[2] <= upper_limit) //if the xc yc is inside the triangle
	return true;

      return false;
    }


    inline void Clear(ModelPart::NodesContainerType::iterator node_it,  int step_data_size )
    {
      unsigned int buffer_size = node_it->GetBufferSize();
	
      for(unsigned int step = 0; step<buffer_size; step++)
	{
	  //getting the data of the solution step
	  double* step_data = (node_it)->SolutionStepData().Data(step);
	    
	  //copying this data in the position of the vector we are interested in
	  for(int j= 0; j< step_data_size; j++)
	    {
	      step_data[j] = 0.0;
	    }
	}
	
    }
      
    inline void ClearVariables(ModelPart::NodesContainerType::iterator node_it , Variable<array_1d<double,3> >& rVariable)
    {
      array_1d<double, 3>& Aux_var = node_it->FastGetSolutionStepValue(rVariable, 0);
	
      noalias(Aux_var) = ZeroVector(3);
	
    }
      
      
    inline void ClearVariables(ModelPart::NodesContainerType::iterator node_it,  Variable<double>& rVariable)
    {
      double& Aux_var = node_it->FastGetSolutionStepValue(rVariable, 0);
	
      Aux_var = 0.0;
	
    }

	
    ///@}

  }; // Class LaplacianSmoothing

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    LaplacianSmoothing& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const LaplacianSmoothing& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif // KRATOS_LAPLACIAN_SMOOTHING_H_INCLUDED  defined 


