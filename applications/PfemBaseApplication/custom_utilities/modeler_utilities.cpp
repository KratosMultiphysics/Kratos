//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

// System includes

// External includes

// Project includes
#include "includes/kratos_flags.h"
#include "custom_utilities/modeler_utilities.hpp"


namespace Kratos
{

  /**
   * Flags related to the meshing parameters
   */

  //meshing options

  //(configuration)
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, REMESH,               0 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, REFINE,               1 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, RECONNECT,            2 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, CONSTRAINED,          3 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, CONTACT_SEARCH,       4 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, MESH_SMOOTHING,       5 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, VARIABLES_SMOOTHING,  6 );

  //execution options (tessellation)
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, NEIGHBOURS_SEARCH,    7 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, BOUNDARIES_SEARCH,    8 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, SET_DOF,              9 );

  //removing options

  //(configuration)
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, REMOVE_NODES,                       0 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, REMOVE_NODES_ON_DISTANCE,           1 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, REMOVE_NODES_ON_ERROR,              2 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, REMOVE_NODES_ON_THRESHOLD,          3 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, REMOVE_BOUNDARY_NODES,              4 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, REMOVE_BOUNDARY_NODES_ON_DISTANCE,  5 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, REMOVE_BOUNDARY_NODES_ON_ERROR,     6 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, REMOVE_BOUNDARY_NODES_ON_THRESHOLD, 7 );

  //refining options

  //(configuration)
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, REFINE_ADD_NODES,             0 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, REFINE_INSERT_NODES,          1 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, REFINE_ELEMENTS,              2 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, REFINE_ELEMENTS_ON_DISTANCE,  3 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, REFINE_ELEMENTS_ON_ERROR,     4 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, REFINE_ELEMENTS_ON_THRESHOLD, 5 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, REFINE_BOUNDARY,              6 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, REFINE_BOUNDARY_ON_DISTANCE,  7 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, REFINE_BOUNDARY_ON_ERROR,     8 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, REFINE_BOUNDARY_ON_THRESHOLD, 9 );

  //execution options

  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, SELECT_ELEMENTS,       0 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, PASS_ALPHA_SHAPE,      1 );
  KRATOS_CREATE_LOCAL_FLAG ( ModelerUtilities, ENGAGED_NODES,         2 );


  //*******************************************************************************************
  //*******************************************************************************************

  void ModelerUtilities::SetDomainLabels (ModelPart& rModelPart)
  {

    unsigned int start=0;
    unsigned int NumberOfMeshes=rModelPart.NumberOfMeshes();
    if(NumberOfMeshes>1) 
      start=1;
      

    for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
      {
	for(ModelPart::NodesContainerType::iterator i_node = rModelPart.NodesBegin(MeshId) ; i_node != rModelPart.NodesEnd(MeshId) ; i_node++)
	  {

	    i_node->SetValue(DOMAIN_LABEL,MeshId);
	  }
      }

  }


  //*******************************************************************************************
  //*******************************************************************************************

  void ModelerUtilities::BuildTotalMesh (ModelPart& rModelPart, int EchoLevel)
  {
    //Mesh Id=0

    if( EchoLevel > 0 )
      std::cout<<"   [ START MESH [Id=0] [Elems=:"<<rModelPart.NumberOfElements()<<"|Nodes="<<rModelPart.NumberOfNodes()<<"|Conds="<<rModelPart.NumberOfConditions()<<"] ] "<<std::endl;      

    rModelPart.Nodes().clear();
    rModelPart.Elements().clear();

    //contact conditions are located on Mesh_0
    ModelPart::ConditionsContainerType KeepConditions;


    //std::cout<<" [ Number of Meshes "<<rModelPart.GetMeshes().size()-1<<" ]"<<std::endl;

    unsigned int nodeId=1;
    unsigned int elemId=1;
    unsigned int condId=1;

    unsigned int start=0;
    unsigned int NumberOfMeshes=rModelPart.NumberOfMeshes();
    if(NumberOfMeshes>1) 
      start=1;
      

    for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
      {
	if( EchoLevel > 0 )
	  std::cout<<"    [ CHILD MESH [Id:"<<MeshId<<"] [Elems="<<rModelPart.NumberOfElements(MeshId)<<"|Nodes="<<rModelPart.NumberOfNodes(MeshId)<<"|Conds="<<rModelPart.NumberOfConditions(MeshId)<<"] ] "<<std::endl;


	for(ModelPart::ElementsContainerType::iterator i_elem = rModelPart.ElementsBegin(MeshId) ; i_elem != rModelPart.ElementsEnd(MeshId) ; i_elem++)
	  {
	    PointsArrayType& vertices=i_elem->GetGeometry().Points();
	    for(unsigned int i=0; i<vertices.size(); i++)
	      {
		vertices[i].Set(BLOCKED);
	      }

	    (rModelPart.Elements()).push_back(*(i_elem.base()));	
	    rModelPart.Elements().back().SetId(elemId);
	    elemId+=1;
	  }



	//Clean Nodes when redefining the total mesh:
	const array_1d<double,3> ZeroNormal(3,0.0);
	ModelPart::NodesContainerType temporal_nodes;
	temporal_nodes.reserve(rModelPart.Nodes(MeshId).size());
	temporal_nodes.swap(rModelPart.Nodes(MeshId));
	  
	for(ModelPart::NodesContainerType::iterator i_node = temporal_nodes.begin() ; i_node != temporal_nodes.end() ; i_node++)
	  {
	    //i_node->PrintInfo(std::cout);
	    //std::cout<<std::endl;
	    if(i_node->Is(BLOCKED)){
		
	      i_node->Reset(NEW_ENTITY); //reset here if the node is labeled as insert 
	      i_node->Reset(TO_REFINE); //reset here if the node is labeled as refine (to not duplicate boundary conditions)
	      i_node->Reset(BLOCKED); 

	      (rModelPart.Nodes(MeshId)).push_back(*(i_node.base()));
	      (rModelPart.Nodes()).push_back(*(i_node.base()));	
	      rModelPart.Nodes().back().SetId(nodeId);
	      nodeId+=1;

	    }
	    else if (!rModelPart.NumberOfElements(MeshId)){
		
	      (rModelPart.Nodes(MeshId)).push_back(*(i_node.base()));
	      (rModelPart.Nodes()).push_back(*(i_node.base()));	
	      rModelPart.Nodes().back().SetId(nodeId);
	      nodeId+=1;
	    }
	    else{
	      std::cout<<" NOT ENGAGED NODE "<<i_node->Id()<<std::endl;
	    }
		
	    if(i_node->IsNot(BOUNDARY))
	      {
		noalias(i_node->GetSolutionStepValue(CONTACT_FORCE)) = ZeroNormal;
		noalias(i_node->GetSolutionStepValue(CONTACT_FORCE)) = ZeroNormal;
	      }


	  }
	  
	//rModelPart.Nodes(MeshId).Sort();  

	for(ModelPart::ConditionsContainerType::iterator i_cond = rModelPart.ConditionsBegin(MeshId) ; i_cond != rModelPart.ConditionsEnd(MeshId) ; i_cond++)
	  {
	    i_cond->Reset(NEW_ENTITY); //reset here if the node is inserted
	    KeepConditions.push_back(*(i_cond.base()));
	    KeepConditions.back().SetId(condId);
	    condId+=1;	
	  }
      }


    for(ModelPart::ConditionsContainerType::iterator i_cond = rModelPart.ConditionsBegin(); i_cond!= rModelPart.ConditionsEnd(); i_cond++)
      {
	if(i_cond->Is(CONTACT)){
	  KeepConditions.push_back(*(i_cond.base()));
	  KeepConditions.back().SetId(condId);
	  condId+=1;
	}
      }
      
    rModelPart.Conditions().swap(KeepConditions);


    //Sort
    rModelPart.Nodes().Sort();
    rModelPart.Elements().Sort();
    rModelPart.Conditions().Sort();     

    //Unique
    rModelPart.Nodes().Unique();
    rModelPart.Elements().Unique();
    rModelPart.Conditions().Unique();

    //Sort Again to have coherent numeration for nodes (mesh with shared nodes)
    unsigned int consecutive_index = 1;
    for(ModelPart::NodesContainerType::iterator in = rModelPart.NodesBegin(0) ; in != rModelPart.NodesEnd(0) ; in++)
      in->SetId(consecutive_index++);
      
    if( EchoLevel > 0 )
      std::cout<<"   [ END MESH [Id=0] [Elems=:"<<rModelPart.NumberOfElements()<<"|Nodes="<<rModelPart.NumberOfNodes()<<"|Conds="<<rModelPart.NumberOfConditions()<<"] ] "<<std::endl;      
 
  }
  

  //*******************************************************************************************
  //*******************************************************************************************

  void ModelerUtilities::CleanMeshFlags(ModelPart& rModelPart,ModelPart::IndexType MeshId)
  {
    
    for(ModelPart::NodesContainerType::const_iterator i_node = rModelPart.NodesBegin(MeshId); i_node != rModelPart.NodesEnd(MeshId); i_node++)
      {

	i_node->Reset(NEW_ENTITY); //reset here if the node is labeled as insert 
	i_node->Reset(TO_REFINE);  //reset here if the node is labeled as refine (to not duplicate bo

      }

    for(ModelPart::ConditionsContainerType::iterator i_cond = rModelPart.ConditionsBegin(MeshId) ; i_cond != rModelPart.ConditionsEnd(MeshId) ; i_cond++)
      {
	i_cond->Reset(NEW_ENTITY); //reset here if the node is inserted
      }
  }


  //*******************************************************************************************
  //*******************************************************************************************


  void ModelerUtilities::CleanRemovedNodes(ModelPart& rModelPart,ModelPart::IndexType MeshId)
  {
    KRATOS_TRY

    //MESH 0 total domain mesh
    ModelPart::NodesContainerType temporal_nodes;
    temporal_nodes.reserve(rModelPart.Nodes(MeshId).size());
	
    temporal_nodes.swap(rModelPart.Nodes(MeshId));
	    
    for(ModelPart::NodesContainerType::iterator i_node = temporal_nodes.begin() ; i_node != temporal_nodes.end() ; i_node++)
      {
	if( i_node->IsNot(TO_ERASE) ){
	  (rModelPart.Nodes(MeshId)).push_back(*(i_node.base()));	
	}
	else{
	  if( i_node->Is(BOUNDARY) )
	    std::cout<<"   BOUNDARY NODE RELEASED "<<i_node->Id()<<std::endl;
	}
      }
	
	
    rModelPart.Nodes(MeshId).Sort();
	

    KRATOS_CATCH( "" )
  }


  //*******************************************************************************************
  //*******************************************************************************************

 
  bool ModelerUtilities::CheckSubdomain(Geometry<Node<3> >& rGeometry)
  {

    KRATOS_TRY

    unsigned int DomainLabel = rGeometry[0].GetValue(DOMAIN_LABEL); //DOMAIN_LABEL must be set as nodal variable
      
    int samesbd=0;
      
    for(int pn=1; pn<3; pn++)
      {
	if(DomainLabel!=rGeometry[pn].GetValue(DOMAIN_LABEL))
	  {
	    samesbd++;
	  }
      }
      
      
    if(samesbd>0)
      return false;

    return true;

    KRATOS_CATCH( "" )
  }
  

  //*******************************************************************************************
  //*******************************************************************************************

  bool ModelerUtilities::CheckInnerCentre(Geometry<Node<3> >& rGeometry)
  {

    KRATOS_TRY

    bool inner = true;

    unsigned int BoundaryNodes = 0;
    const unsigned int size = rGeometry.size();

    for(unsigned int i = 0; i < size; i++)
      {
	if(rGeometry[i].Is(BOUNDARY))
	  {
	    BoundaryNodes += 1;
	  }
      }


    if(BoundaryNodes == size)
      {

	//Baricenter
	array_1d<double, 3>  Center;
	Center.clear();
	array_1d<double, 3>  Normal;

	std::vector<array_1d<double, 3> > Vertices;
	array_1d<double, 3>  Vertex;
	
	
	for(unsigned int i = 0; i < size; i++)
	  {    
	    Vertex  = rGeometry[i].Coordinates();

	    Vertices.push_back(Vertex);

	    Center += Vertex;
	  }
	

	Center /= (double)size;
	    
	array_1d<double, 3> Corner;

	double tolerance = 0.05;
	int numouter     = 0;


	for(unsigned int i = 0; i < size; i++)
	  {

	    Normal = rGeometry[i].FastGetSolutionStepValue(NORMAL); 

	    double NormNormal = norm_2(Normal);
	    if( NormNormal != 0)
	      Normal /= NormNormal;

	    //change position to be the vector from the vertex to the geometry center
	    Corner = Center-Vertices[i];

	    double NormCorner = norm_2(Corner);
	    if( NormCorner != 0 )
	      Corner/= NormCorner;

	    double projection = inner_prod(Corner,Normal);

	    if( projection > tolerance )
	      {
		numouter++;
	      }
	  }


	if( numouter > 0 )
	  inner = false;
	
      }

    return inner; //if is inside the body domain returns true

    KRATOS_CATCH( "" )
  }


  //*******************************************************************************************
  //*******************************************************************************************

  bool ModelerUtilities::CheckOuterCentre(Geometry<Node<3> >& rGeometry, double& rOffsetFactor, bool& rSelfContact)
  {
    KRATOS_TRY

    bool outer=false;

    unsigned int BoundaryNodes = 0;
    const unsigned int size = rGeometry.size();

    for(unsigned int i = 0; i < size; i++)
      {
	if(rGeometry[i].Is(BOUNDARY))
	  {
	    BoundaryNodes += 1;
	  }
      }


    if(BoundaryNodes == size)
      {

	//Baricenter
	array_1d<double, 3>  Center;
	Center.clear();
	std::vector<array_1d<double, 3> > Vertices;
	array_1d<double, 3>  Vertex;
	array_1d<double, 3>  Normal;
	double Shrink = 0;
	
	for(unsigned int i = 0; i < size; i++)
	  {
	    Normal  = rGeometry[i].FastGetSolutionStepValue(NORMAL);

	    double NormNormal = norm_2(Normal);
	    if( NormNormal != 0)
	      Normal /= NormNormal;

	    Shrink  = rGeometry[i].FastGetSolutionStepValue(SHRINK_FACTOR);
	    
	    Normal *= Shrink * rOffsetFactor;
	    
	    Vertex  = rGeometry[i].Coordinates() - Normal;

	    Vertices.push_back(Vertex);

	    Center  += Vertex;
	  }
	
	Center /= (double)size;
	    
	double ortho  = 0.15;
	double slope  = 0.25;  //error assumed for some elements in the corners < 45 degrees
	double extra  = 0.95;

	int numouter         =0;
	int numextra         =0;
	int numcoplanar      =0;
	int numsamedirection =0;
	int numorthogonal    =0;

	array_1d<double, 3> Coplanar = rGeometry[0].FastGetSolutionStepValue(NORMAL); 
	double NormCoplanar = norm_2(Coplanar);
	if( NormCoplanar != 0)
	  Coplanar /= NormCoplanar;

	array_1d<double, 3> Corner;

	for(unsigned int i = 0; i < size; i++)
	  {
	   
	    //std::cout<<" V: ("<<rGeometry[i].Id()<<"): ["<<rGeometry[i].X()<<", "<<rGeometry[i].Y()<<"]  Normal:"<<rGeometry[i].FastGetSolutionStepValue(NORMAL)<<std::endl;

	    Normal = rGeometry[i].FastGetSolutionStepValue(NORMAL); 

	    double NormNormal = norm_2(Normal);
	    if( NormNormal != 0)
	      Normal /= NormNormal;
				

	    //change position to be the vector from the vertex to the geometry center
	    Corner = Vertices[i]-Center;

	    if(norm_2(Corner))
	      Corner/= norm_2(Corner);

	    double projection = inner_prod(Corner,Normal);

	    if( projection < slope)
	      {
		numouter++;
	      }
	    else
	      {
		if( projection < extra)
		  numextra++;
	      }
		
	    double coplanar = inner_prod(Coplanar,Normal);

	    //std::cout<<" V["<<i<<"]: "<<rGeometry[i]<<" Normal:"<<Normal<<" coplanar "<<fabs(coplanar)<<" < "<<ortho<<std::endl;

	    if(coplanar>0){
	      numsamedirection++;
	    }
		
	    if(coplanar>extra){
	      numcoplanar++;
	    }
		
	    if(fabs(coplanar)<=ortho){
	      numorthogonal++;
	    }

	  }

	//std::cout<<std::endl;
	//std::cout<<"  [ no:"<<numouter<<";ne:"<<numextra<<";nc:"<<numcoplanar<<";ns: "<<numsamedirection<<";nor:"<<numorthogonal<<"]"<<std::endl;
	
	int num = (int)size;

	if(numouter==num)
	  outer=true;

	if(numouter==(num-1) && numextra==1)
	  outer=true;
	
	if(numouter>0 && (numextra>0 && numorthogonal>0) && !rSelfContact){
	  outer=true;
	  std::cout<<"   Element with "<<num<<" corners accepted:case1 "<<std::endl;
	}
	
	if(numouter==0 && (numextra>(num-2) && numorthogonal>0) && !rSelfContact){
	  outer=true;
	  std::cout<<"   Element with "<<num<<" corners accepted:case2 "<<std::endl;
	}

	if(numcoplanar==num)
	  outer=false;
	    
	if(numsamedirection==num && numorthogonal==0)
	  outer=false;
	    
	// if(numorthogonal>=1)
	//   outer=false;
	    

      }

    return outer; //if is outside the body domain returns true

    KRATOS_CATCH( "" )
  }

  
  //*******************************************************************************************
  //*******************************************************************************************

  bool ModelerUtilities::CheckGeometryShape(Geometry<Node<3> >& rGeometry, double& rShape)
  {
    KRATOS_TRY
      
    bool sliver    = false;
    bool distorted = false;

    const unsigned int  size = rGeometry.size();

    //check if is a sliver and the distorsion of the sliver

    double Volume = 0;
    double MaximumSideLength = 0;
    //double MinimumSideLength = 0;
    
    double CriticalRelativeSideLength = (double)size; //edge relative length
    
    Volume = rGeometry.Volume();
    
    //compare side lengths
    double RelativeSideLength = 0;
    //RelativeSideLength = CompareSideLengths(rGeometry, MaximumSideLength, MinimumSideLength);
 
    if( RelativeSideLength > CriticalRelativeSideLength ){
      distorted = true;
    }
      
    double GeometrySideVolume =  1e-6 * pow( MaximumSideLength, size );

    if( Volume < GeometrySideVolume ){
      sliver = true;
    }
    
    //check if is an SELF CONTACT element (contact domain definition)
    std::vector<int> SlaveVertices;
    //CheckContactElement(rGeometry, SlaveVertices);

    //check if is an EDGE or a FACE element (contact domain definition)
    unsigned int slaves = SlaveVertices.size();
    if( slaves == 1 ){ //FACE
      
      //check the projection of the slave vertex on the geometry face
      double AreaTolerance = 2; //if is outside of the face (only a 2*FaceArea deviation is allowed)
      double FaceArea = 0;
      double ProjectedArea = 0;
      //FaceArea = ComputeFaceArea(rGeometry, SlaveVertex);
      //ProjectedArea = ComputePointToFaceProjection(rGeometry, SlaveVertex);

      if( ProjectedArea < 0 ){ // projection outside of the face
	
	if( FaceArea < AreaTolerance * fabs(ProjectedArea) ) 
	  distorted = true;
      }

    }
    else if( slaves > 1 ){ //EDGE

      //compare vertex normals (detect coplanar faces and orthogonal faces)
      //if( !CheckVertexNormals(rGeometry) )
      distorted = true;
    }
    else if( slaves == 0 ){ //SELF CONTACT  :: element with crossed diagonals (strange situation)
      
      //compare vertex normals (detect coplanar faces and orthogonal faces)
      //if( !CheckVertexNormals(rGeometry) )
      distorted = true;
    }

    if( sliver )
      rShape = 1;
    else
      rShape = 0;
   
    if( distorted )
      return false; //returns false if the geometry has not the correct shape and has to be removed
    else
      return true;  //returns true if the geometry has the correct shape and is kept
    
    
    KRATOS_CATCH( "" )
  }


  //*******************************************************************************************
  //*******************************************************************************************

  double ModelerUtilities::FindBoundaryH (Node<3>& BoundaryPoint)
  {
    KRATOS_TRY

    double havg = 0.00;
      
    if((BoundaryPoint.GetValue(NEIGHBOUR_NODES)).size() != 0)
      {
	double xc = BoundaryPoint.X();
	double yc = BoundaryPoint.Y();
	double zc = BoundaryPoint.Z();

	double h_nodes = 0;
	double h = 1000.0;
	for( WeakPointerVector< Node<3> >::iterator i = BoundaryPoint.GetValue(NEIGHBOUR_NODES).begin();
	     i !=  BoundaryPoint.GetValue(NEIGHBOUR_NODES).end(); i++)
	  {
	    if( i->Is(BOUNDARY) ){
	      double x = i->X();
	      double y = i->Y();
	      double z = i->Z();
	      double l = (x-xc)*(x-xc);
	      l += (y-yc)*(y-yc);
	      l += (z-zc)*(z-zc);
				    
	      if(l<h) h = l;
		
	      h = sqrt(h);
	      havg += h;
	      h_nodes += 1;
	    }
	  }
			      
	//calculate average h
	if(h_nodes == 0)
	  KRATOS_THROW_ERROR( std::logic_error,"no node has neighbours!!!!", "" )
			      
	havg /= h_nodes;
      }
			  
    return havg;

    KRATOS_CATCH( "" )
  }




  //*******************************************************************************************
  //*******************************************************************************************

  //returns false if it should be removed
  double& ModelerUtilities::ComputeRadius(double& rRadius, double& rVolume, std::vector<Vector >& rVertices,const unsigned int& dimension)
  {
    KRATOS_TRY
        
    if( dimension == 2 ){

      bounded_matrix<double,2,2> mJ;    //local jacobian
      bounded_matrix<double,2,2> mJinv; //inverse jacobian

      //calculation of the jacobian  //coordinate center point 0
      for(unsigned int i = 0; i < dimension; i++)
	{
	  for(unsigned int j = 0; j < dimension; j++)
	    {
	      mJ(i,j)=rVertices[i+1][j]-rVertices[0][j];
	    }
	}
      
      mJ *= 2.0;

      //calculation of the determinant (volume/2)
      rVolume = mJ(0,0)*mJ(1,1)-mJ(0,1)*mJ(1,0); //detJ

      //calculation of the inverse of the jacobian
      mJinv(0,0) =  mJ(1,1);
      mJinv(0,1) = -mJ(0,1);
      mJinv(1,0) = -mJ(1,0);
      mJinv(1,1) =  mJ(0,0);
        
      mJinv /= rVolume;
  
      //calculation of the center
      Vector Center = ZeroVector(3);    //center pos

      //center point 0
      for(unsigned int i = 0; i < dimension; i++)
	{
	  for(unsigned int j = 0; j < dimension; j++)
	    {
	      Center[i]  += (rVertices[i+1][j] * rVertices[i+1][j]);
	      Center[i]	 -= (rVertices[0][j]   * rVertices[0][j]  );
	    }
	}

      Center = prod(mJinv,Center);

      //calculate the element radius
      Center -= rVertices[0];

      rRadius = norm_2(Center);
      
    }
    else if( dimension == 3 ){

      bounded_vector<double,3>   mRHS;  //center pos
      bounded_matrix<double,3,3> mJ;    //local jacobian
      bounded_matrix<double,3,3> mJinv; //inverse jacobian


      //calculation of the jacobian  //coordinate center point 0
      for(unsigned int i = 0; i < dimension; i++)
	{
	  for(unsigned int j = 0; j < dimension; j++)
	    {
	      mJ(i,j)=rVertices[i+1][j]-rVertices[0][j];
	    }
	}

      //calculation of the inverse of the jacobian
      //first column
      mJinv(0,0) =  mJ(1,1)*mJ(2,2) - mJ(1,2)*mJ(2,1);
      mJinv(1,0) = -mJ(1,0)*mJ(2,2) + mJ(1,2)*mJ(2,0);
      mJinv(2,0) =  mJ(1,0)*mJ(2,1) - mJ(1,1)*mJ(2,0);		
      //second column
      mJinv(0,1) = -mJ(0,1)*mJ(2,2) + mJ(0,2)*mJ(2,1);
      mJinv(1,1) =  mJ(0,0)*mJ(2,2) - mJ(0,2)*mJ(2,0);
      mJinv(2,1) = -mJ(0,0)*mJ(2,1) + mJ(0,1)*mJ(2,0);
      //third column
      mJinv(0,2) =  mJ(0,1)*mJ(1,2) - mJ(0,2)*mJ(1,1);
      mJinv(1,2) = -mJ(0,0)*mJ(1,2) + mJ(0,2)*mJ(1,0);
      mJinv(2,2) =  mJ(0,0)*mJ(1,1) - mJ(0,1)*mJ(1,0);
      

      //calculation of the determinant (volume/6) 
      rVolume = mJ(0,0)*mJinv(0,0)+ mJ(0,1)*mJinv(1,0)+ mJ(0,2)*mJinv(2,0); //detJ


      //calculation of the center
      Vector Center = ZeroVector(3);    //center pos
  
      mRHS[0]= (mJ(0,0)*mJ(0,0)+ mJ(0,1)*mJ(0,1)+ mJ(0,2)*mJ(0,2));
      mRHS[1]= (mJ(1,0)*mJ(1,0)+ mJ(1,1)*mJ(1,1)+ mJ(1,2)*mJ(1,2));
      mRHS[2]= (mJ(2,0)*mJ(2,0)+ mJ(2,1)*mJ(2,1)+ mJ(2,2)*mJ(2,2));


      //bxc
      Center[0] = (+mRHS[0])*(mJ(1,1)*mJ(2,2)-mJ(1,2)*mJ(2,1));
      Center[1] = (-mRHS[0])*(mJ(1,0)*mJ(2,2)-mJ(1,2)*mJ(2,0));
      Center[2] = (+mRHS[0])*(mJ(1,0)*mJ(2,1)-mJ(1,1)*mJ(2,0));
  
      //cxa
      Center[0]+= (+mRHS[1])*(mJ(2,1)*mJ(0,2)-mJ(2,2)*mJ(0,1));
      Center[1]+= (-mRHS[1])*(mJ(2,0)*mJ(0,2)-mJ(2,2)*mJ(0,0));
      Center[2]+= (+mRHS[1])*(mJ(2,0)*mJ(0,1)-mJ(2,1)*mJ(0,0));
  
      //axb
      Center[0]+= (+mRHS[2])*(mJ(0,1)*mJ(1,2)-mJ(0,2)*mJ(1,1));
      Center[1]+= (-mRHS[2])*(mJ(0,0)*mJ(1,2)-mJ(0,2)*mJ(1,0));
      Center[2]+= (+mRHS[2])*(mJ(0,0)*mJ(1,1)-mJ(0,1)*mJ(1,0));
  
      Center /= (2.0*rVolume);
    
      //calculate the element radius
      rRadius = norm_2(Center);

     }
    else{

      rRadius = 0;
      KRATOS_THROW_ERROR(std::logic_error, "WorkingSpaceDimension not Correct in Radius AlphaShape calculation", "" )
      
    }

    return rRadius;


    KRATOS_CATCH( "" )      
  }


  //*******************************************************************************************
  //*******************************************************************************************

  
  //returns false if it should be removed
  bool ModelerUtilities::AlphaShape(double AlphaParameter, Geometry<Node<3> >& rGeometry, const unsigned int dimension)
  {
    KRATOS_TRY
   
    const unsigned int size      = rGeometry.size();

    //calculate geometry radius and volume
    double Radius = 0;
    double Volume = 0;
    std::vector<Vector > Vertices;
    Vector Vertex = ZeroVector(3);
    for(unsigned int i = 0; i < size; i++)
      {
	Vertex[0] = rGeometry[i].X();
	Vertex[1] = rGeometry[i].Y();
	Vertex[2] = rGeometry[i].Z();

	Vertices.push_back(Vertex);
      }
    
    Radius = ComputeRadius(Radius, Volume, Vertices, dimension);

    //calculate average h
    double h = 0;
    
    for( unsigned int i = 0; i<size; i++ )
      {
	h += rGeometry[i].FastGetSolutionStepValue(NODAL_H);
      }

    h /= (double)size;

    double CriticalVolume = 1e-6 * pow(h, size-1);
    double AlphaRadius    = AlphaParameter * h;


    if( Volume < CriticalVolume ) //sliver
      {
	return false;
      }
    else
      {

	if( Radius < AlphaRadius )
	  {
	    return true;
	  }
	else
	  {
	    return false;
	  }
      }


    KRATOS_CATCH( "" )
  }


  //*******************************************************************************************
  //*******************************************************************************************

  //returns false if it should be removed
  bool ModelerUtilities::ShrankAlphaShape(double AlphaParameter, Geometry<Node<3> >& rGeometry,double& rOffsetFactor, const unsigned int dimension)
  {
    KRATOS_TRY

    const unsigned int size      = rGeometry.size();

    //calculate geometry radius and volume
    double Radius = 0;
    double Volume = 0;
    std::vector<Vector > Vertices;

    Vector Vertex = ZeroVector(3);
    array_1d<double, 3>  Normal;
    double Shrink = 0;

    for(unsigned int i = 0; i < size; i++)
      {
	Normal    = rGeometry[i].FastGetSolutionStepValue(NORMAL);

	double NormNormal = norm_2(Normal);
	if( NormNormal != 0)
	  Normal /= NormNormal;

	Shrink    = rGeometry[i].FastGetSolutionStepValue(SHRINK_FACTOR);
	
	Normal *= Shrink * rOffsetFactor;

	Vertex[0] = rGeometry[i].X() - Normal[0];
	Vertex[1] = rGeometry[i].Y() - Normal[1];
	Vertex[2] = rGeometry[i].Z() - Normal[2];

	Vertices.push_back(Vertex);
      }
    
    Radius = ComputeRadius(Radius, Volume, Vertices, dimension);


    //calculate average h and  average h of boundary face
    double h = 0;
    double h_face = 0;

    for( unsigned int i = 0; i<size; i++ )
      {
	h += rGeometry[i].FastGetSolutionStepValue(NODAL_H);
	h_face += FindBoundaryH(rGeometry[i]);
      }

    h /= (double)size;
    h_face /= (double)size;

   
    if(h_face>h)
      h = h_face;

    double ExtraAlpha = 1.4;
    
    double CriticalVolume = 1e-6 * pow(h, size-1);
    double AlphaRadius    = AlphaParameter * h * ExtraAlpha;


    if( Volume < CriticalVolume ) //sliver
      {
	return false;
      }
    else
      {

	if( Radius < AlphaRadius )
	  {
	    return true;
	  }
	else
	  {
	    return false;
	  }
      }


    KRATOS_CATCH( "" )
  }



  //*******************************************************************************************
  //*******************************************************************************************

  void ModelerUtilities::CheckParticles (ModelPart& rModelPart,ModelPart::IndexType MeshId)
  {
    KRATOS_TRY

    int NumberOfNodes = rModelPart.NumberOfNodes();
    std::cout<<" Number of Nodes "<<NumberOfNodes<<std::endl;
    for(int id=1; id<=NumberOfNodes; id++)
      {
	std::cout<<" Check node: "<<id<<std::endl;
	if(rModelPart.Nodes()[id].Is(BOUNDARY)){
	  std::cout<<" Node : "<<id<<" is boundary "<<std::endl;
	}
	else{
	  std::cout<<" Node : "<<id<<" is not boundary "<<std::endl;
	}  

      }

    KRATOS_CATCH( "" )
  }


  //*******************************************************************************************
  //*******************************************************************************************

  bool ModelerUtilities::CheckConditionInBox(Condition::Pointer& pCondition, SpatialBoundingBox& rRefiningBox, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    bool inside = true;
    Vector Point(3);
    Geometry< Node<3> >& rGeometry = pCondition->GetGeometry();

    for(unsigned int i=0; i<rGeometry.size(); i++)
      {
	Point[0] = rGeometry[i].X();
	Point[1] = rGeometry[i].Y();
	Point[2] = rGeometry[i].Z();
	if( !rRefiningBox.IsInside(Point,rCurrentProcessInfo[TIME]) ){
	  inside = false;
	  break;
	}
      }

    return inside;

    KRATOS_CATCH( "" )
  }

  //*******************************************************************************************
  //*******************************************************************************************

  bool ModelerUtilities::CheckElementInBox(Element::Pointer& pElement, SpatialBoundingBox& rRefiningBox, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY
    
    bool inside = true;
    Vector Point(3);
    Geometry< Node<3> >& rGeometry = pElement->GetGeometry();

    for(unsigned int i=0; i<rGeometry.size(); i++)
      {
	Point[0] = rGeometry[i].X();
	Point[1] = rGeometry[i].Y();
	Point[2] = rGeometry[i].Z();
	if( !rRefiningBox.IsInside(Point,rCurrentProcessInfo[TIME]) ){
	  inside = false;
	  break;
	}
      }

    return inside;

    KRATOS_CATCH( "" )
  }

  //*******************************************************************************************
  //*******************************************************************************************

  bool ModelerUtilities::CheckVerticesInBox(Geometry<Node<3> >& rGeometry, SpatialBoundingBox& rRefiningBox, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    bool inside = true;
    Vector Point(3);

    for(unsigned int i=0; i<rGeometry.size(); i++)
      {
	Point[0] = rGeometry[i].X();
	Point[1] = rGeometry[i].Y();
	Point[2] = rGeometry[i].Z();
	if( !rRefiningBox.IsInside(Point,rCurrentProcessInfo[TIME]) ){
	  inside = false;
	  break;
	}
      }

    return inside;

    KRATOS_CATCH( "" )
  }


  //*******************************************************************************************
  //*******************************************************************************************

  Condition::Pointer ModelerUtilities::FindMasterCondition(Condition::Pointer& pCondition, ModelPart::ConditionsContainerType & rModelConditions,bool & condition_found)
  {
    KRATOS_TRY
		
    Condition::Pointer pMasterCondition;

    Geometry< Node<3> >& rGeometry = pCondition->GetGeometry();
    boost::numeric::ublas::matrix<unsigned int> lpofa; //points that define the faces
    rGeometry.NodesInFaces(lpofa);   
		
    //std::cout<<" lpofa "<<lpofa<<std::endl;
    //std::cout<<" rGeometry "<<rGeometry<<std::endl;

    condition_found=false;
    for(ModelPart::ConditionsContainerType::iterator ic=rModelConditions.begin(); ic!=rModelConditions.end(); ic++)
      {
	//2D edges:
	if(ic->IsNot(CONTACT)){

	  Geometry< Node<3> >& rConditionGeom = ic->GetGeometry();
			
	  for(unsigned int i=0; i<lpofa.size2();i++)
	    {
	      // std::cout<<" General Conditions IDs ["<<rConditionGeom[0].Id()<<"] ["<<rConditionGeom[1].Id()<<"] "<<std::endl;
	      // std::cout<<" Local Conditions IDs ("<<i<<"):["<<rGeometry[lpofa(1,i)].Id()<<"] ["<<rGeometry[lpofa(2,i)].Id()<<"] "<<std::endl;

	      if( (   rConditionGeom[0].Id() == rGeometry[lpofa(1,i)].Id() 
		      && rConditionGeom[1].Id() == rGeometry[lpofa(2,i)].Id() ) || 
		  (   rConditionGeom[0].Id() == rGeometry[lpofa(2,i)].Id() 
		      && rConditionGeom[1].Id() == rGeometry[lpofa(1,i)].Id() ) )
		{
		  pMasterCondition= *(ic.base());
		  condition_found=true;
		  break;
		}
			       
	    }

	}
	if(condition_found)
	  {
	    break;
	  }
			
      }

    if(!condition_found) {
      //   //KRATOS_THROW_ERROR(std::logic_error, "Boundary Condition NOT FOUND after CONTACT MESHING SEARCH", "" )
      std::cout<<" Boundary Condition NOT FOUND after CONTACT MESHING SEARCH "<<std::endl;
      //   std::cout<<" rGeometry "<<rGeometry<<std::endl;

      //   for(ModelPart::ConditionsContainerType::iterator ic=rModelConditions.begin(); ic!=rModelConditions.end(); ic++)
      //     {
      //       //2D edges:
			
      //       Geometry< Node<3> >& rConditionGeom = ic->GetGeometry();
			
      for(unsigned int i=0; i<lpofa.size2();i++)
	{
	  // 	  std::cout<<" General Conditions IDs ["<<rConditionGeom[0].Id()<<"] ["<<rConditionGeom[1].Id()<<"] "<<std::endl;
	  std::cout<<" Local Conditions IDs ("<<i<<"):["<<rGeometry[lpofa(1,i)].Id()<<"] ["<<rGeometry[lpofa(2,i)].Id()<<"] "<<std::endl;
		
	}

      //     }
    }

    return pMasterCondition;
	
    KRATOS_CATCH( "" )
  }
  
  //*******************************************************************************************
  //*******************************************************************************************

  Condition::Pointer ModelerUtilities::FindMasterCondition(Condition::Pointer& pCondition, PointType& pSlaveNode, ModelPart::ConditionsContainerType & rModelConditions,bool & condition_found)
  {
    KRATOS_TRY    
    
    Condition::Pointer pMasterCondition;

    Geometry< Node<3> >& rGeometry = pCondition->GetGeometry();
    boost::numeric::ublas::matrix<unsigned int> lpofa; //points that define the faces
    rGeometry.NodesInFaces(lpofa);   
		
    //std::cout<<" lpofa "<<lpofa<<std::endl;
    //std::cout<<" rGeometry "<<rGeometry<<std::endl;

    condition_found=false;
    for(ModelPart::ConditionsContainerType::iterator ic=rModelConditions.begin(); ic!=rModelConditions.end(); ic++)
      {
	//2D edges:
	if(ic->IsNot(CONTACT)){
				      
	  Geometry< Node<3> >& rConditionGeom = ic->GetGeometry();
			
	  for(unsigned int i=0; i<lpofa.size2();i++)
	    {
	      // std::cout<<" General Conditions IDs ["<<rConditionGeom[0].Id()<<"] ["<<rConditionGeom[1].Id()<<"] "<<std::endl;
	      // std::cout<<" Local Conditions IDs ("<<i<<"):["<<rGeometry[lpofa(1,i)].Id()<<"] ["<<rGeometry[lpofa(2,i)].Id()<<"] "<<std::endl;

	      if( (   rConditionGeom[0].Id() == rGeometry[lpofa(1,i)].Id() 
		      && rConditionGeom[1].Id() == rGeometry[lpofa(2,i)].Id() ) || 
		  (   rConditionGeom[0].Id() == rGeometry[lpofa(2,i)].Id() 
		      && rConditionGeom[1].Id() == rGeometry[lpofa(1,i)].Id() ) )
		{
		  pMasterCondition = *(ic.base());
		  pSlaveNode = rGeometry[lpofa(0,i)];
		  //std::cout<<"   Slave_Node: found: "<<rGeometry[lpofa(0,i)].Id()<<std::endl;
		  condition_found=true;
		  break;
		}
			       
	    }
	}
	if(condition_found)
	  {
	    break;
	  }
			
      }

    // if(!found) 
    //     KRATOS_THROW_ERROR( std::logic_error, "Boundary Condition NOT FOUND after CONTACT MESHING SEARCH", "" )

    return pMasterCondition;

    KRATOS_CATCH( "" )

  }




  //*******************************************************************************************
  //*******************************************************************************************

 
  bool ModelerUtilities::CheckContactActive(GeometryType& rConditionGeometry, bool& rSemiActiveContact, std::vector<bool>& rSemiActiveNodes)
  {
    KRATOS_TRY

    unsigned int size = rConditionGeometry.size();
    unsigned int counter = 0;

    rSemiActiveContact = false;
    rSemiActiveNodes.resize(size);
    std::fill( rSemiActiveNodes.begin(), rSemiActiveNodes.end(), false );
     
    for(unsigned int i=0; i<rConditionGeometry.size(); i++){
            
      array_1d<double, 3 > & ContactForceNormal  = rConditionGeometry[i].FastGetSolutionStepValue(CONTACT_FORCE);
      
      if(norm_2(ContactForceNormal)>0){
	rSemiActiveContact  = true;
	rSemiActiveNodes[i] = true;
	counter++;
      }
    }

    if(counter == size)
      return true;
    else
      return false;
      
    KRATOS_CATCH( "" )
  }



  //*******************************************************************************************
  //*******************************************************************************************

  bool ModelerUtilities::CheckNodeCloseWallTip(std::vector<SpatialBoundingBox::Pointer>& rRigidWalls, PointType& rNode, ProcessInfo& rCurrentProcessInfo, double& rFactor)
  {
    KRATOS_TRY
    
    // double tip_radius = 0;
    // Vector tip_center = ZeroVector(3);
    
    Vector Point(3);
    Point[0] = rNode.X();
    Point[1] = rNode.Y();
    Point[2] = rNode.Z();

    int ContactFace = 0; //free surface
    int AuxContactFace = 0;
    for( unsigned int i = 0; i < rRigidWalls.size(); i++ )
      {
	ContactFace = 0;
        AuxContactFace = 0;
	if( rRigidWalls[i]->IsInside( Point, rCurrentProcessInfo[TIME], AuxContactFace ) ){
      ContactFace = AuxContactFace;  // OBS: in the case of circles, ContactFace = 2 and IsInside = false/true. Then, it does not break but it continues.
	  // tip_radius = rRigidWalls[i]->GetRadius();
	  // tip_center = rRigidWalls[i]->GetCenter();
	  break;
	}
      }

    // Point[0]-=tip_center[0];
    // Point[1]-=tip_center[1];
    // Point[2]-=tip_center[2];
    
    // double distance=norm_2(Point);
    
    //criterion A:
    if( ContactFace == 2 )
      return true;

    //criterion B:
    // if( ((1-rFactor)*tip_radius < distance &&  (1+rFactor)*tip_radius > distance) )
    //   return true;
    // else
    //   return false;
      
    return false;

    KRATOS_CATCH( "" )
  }


  //*******************************************************************************************
  //*******************************************************************************************

  double ModelerUtilities::CheckCriticalRadius (ModelPart& rModelPart, double& rCriticalRadius, unsigned int MeshId)
  {
    KRATOS_TRY

    double minimum_h = rCriticalRadius;
    
    for(ModelPart::NodesContainerType::iterator i_node = rModelPart.NodesBegin(MeshId) ; i_node != rModelPart.NodesEnd(MeshId) ; i_node++)
      {
      
	double nodal_h = i_node->FastGetSolutionStepValue(NODAL_H);
	if( nodal_h < rCriticalRadius )
	  minimum_h = nodal_h; 
  
      }

    if( minimum_h < rCriticalRadius )
      std::cout<<"  [ FAULT :: CRITICAL MESH SIZE :: supplied size "<<rCriticalRadius<<" is bigger than initial mesh size "<<minimum_h<<" ] "<<std::endl;

    return minimum_h;
    
    KRATOS_CATCH( "" )

  }

} // Namespace Kratos

