//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_utilities/mesher_utilities.hpp"

namespace Kratos
{

  /**
   * Flags related to the meshing parameters
   */

  //meshing options

  //(configuration)
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, REMESH,               0 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, REFINE,               1 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, RECONNECT,            2 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, TRANSFER,             3 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, CONSTRAINED,          4 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, CONTACT_SEARCH,       5 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, MESH_SMOOTHING,       6 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, VARIABLES_SMOOTHING,  7 );

  //removing options

  //(configuration)
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, REMOVE_NODES,                       0 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, REMOVE_NODES_ON_DISTANCE,           1 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, REMOVE_NODES_ON_ERROR,              2 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, REMOVE_NODES_ON_THRESHOLD,          3 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, REMOVE_BOUNDARY_NODES,              4 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, REMOVE_BOUNDARY_NODES_ON_DISTANCE,  5 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, REMOVE_BOUNDARY_NODES_ON_ERROR,     6 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, REMOVE_BOUNDARY_NODES_ON_THRESHOLD, 7 );

  //refining options

  //(configuration)
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, REFINE_ADD_NODES,             0 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, REFINE_INSERT_NODES,          1 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, REFINE_ELEMENTS,              2 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, REFINE_ELEMENTS_ON_DISTANCE,  3 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, REFINE_ELEMENTS_ON_ERROR,     4 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, REFINE_ELEMENTS_ON_THRESHOLD, 5 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, REFINE_BOUNDARY,              6 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, REFINE_BOUNDARY_ON_DISTANCE,  7 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, REFINE_BOUNDARY_ON_ERROR,     8 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, REFINE_BOUNDARY_ON_THRESHOLD, 9 );

  //execution options

  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, INITIALIZE_MESHER_INPUT,              0 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, FINALIZE_MESHER_INPUT,                1 );

  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, TRANSFER_KRATOS_NODES_TO_MESHER,      2 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, TRANSFER_KRATOS_ELEMENTS_TO_MESHER,   3 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER, 4 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, TRANSFER_KRATOS_FACES_TO_MESHER,      5 );

  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, SELECT_TESSELLATION_ELEMENTS,         6 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, KEEP_ISOLATED_NODES,                  7 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, REFINE_WALL_CORNER,                   8 );

  //execution options (tessellation)
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, NEIGHBOURS_SEARCH,                    8 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, BOUNDARIES_SEARCH,                    9 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, SET_DOF,                             10 );
  KRATOS_CREATE_LOCAL_FLAG ( MesherUtilities, PASS_ALPHA_SHAPE,                    11 );

  //*******************************************************************************************
  //*******************************************************************************************

  void MesherUtilities::SetModelPartNameToElements(ModelPart& rModelPart)
  {

    unsigned int start=0;
    unsigned int NumberOfSubModelParts=rModelPart.NumberOfSubModelParts();

    if(NumberOfSubModelParts>0)
    {
      for(auto& i_mp : rModelPart.SubModelParts())
      {
        if( i_mp.NumberOfElements() != 0 ){
          if( i_mp.Is(BOUNDARY) || i_mp.IsNot(ACTIVE) ){ //wall elements or domain elements (unique model part)
            for(auto& i_elem : i_mp.Elements())
            {
              i_elem.SetValue(MODEL_PART_NAME,i_mp.Name());
            }
          }

        }
      }
    }

  }

  //*******************************************************************************************
  //*******************************************************************************************

  void MesherUtilities::SetModelPartNameToConditions(ModelPart& rModelPart)
  {

    unsigned int start=0;
    unsigned int NumberOfSubModelParts=rModelPart.NumberOfSubModelParts();

    if(NumberOfSubModelParts>0){
      for(auto& i_mp : rModelPart.SubModelParts())
      {
        if( i_mp.NumberOfConditions() != 0 ){
          if( i_mp.Is(BOUNDARY) && i_mp.NumberOfElements() == 0 ){ // only model parts with conditions (unique model part)
            for(auto& i_cond : i_mp.Conditions())
            {
              i_cond.SetValue(MODEL_PART_NAME,i_mp.Name());
            }

          }
        }
      }
    }

  }


  //*******************************************************************************************
  //*******************************************************************************************

  void MesherUtilities::SetModelPartNameToNodes (ModelPart& rModelPart)
  {

    unsigned int start=0;
    unsigned int NumberOfSubModelParts=rModelPart.NumberOfSubModelParts();


    if(NumberOfSubModelParts>0){
      for(auto& i_mp : rModelPart.SubModelParts())
      {

        if( i_mp.NumberOfNodes() != 0 ){
          if( i_mp.Is(BOUNDARY) ){ // shared model parts for nodes in boundary conditions
            for(auto& i_node : i_mp.Nodes())
            {
              i_node.GetValue(MODEL_PART_NAMES).push_back(i_mp.Name());
            }

          }
          else if( i_mp.IsNot(ACTIVE) && i_mp.IsNot(BOUNDARY) ){ //unique domain model part
            for(auto& i_node : i_mp.Nodes())
            {
              i_node.SetValue(MODEL_PART_NAME,i_mp.Name());
            }
          }
        }
      }
    }


  }


  //*******************************************************************************************
  //*******************************************************************************************

  void MesherUtilities::SetFlagsToNodes(ModelPart& rModelPart, const std::vector<Flags> rControlFlags, const std::vector<Flags> rAssignFlags)
  {
    const int nnodes = rModelPart.Nodes().size();
    ModelPart::NodesContainerType::iterator it_begin = rModelPart.NodesBegin();

    #pragma omp parallel for
    for (int i = 0; i < nnodes; i++)
    {
      ModelPart::NodesContainerType::iterator it = it_begin + i;

      for(unsigned int i = 0; i<rControlFlags.size(); i++)
      {
        if( it->Is(rControlFlags[i]) ){
          for(unsigned int i = 0; i<rAssignFlags.size(); i++)
            it->Set(rAssignFlags[i]);
        }
      }
    }
  }

  //*******************************************************************************************
  //*******************************************************************************************

  bool MesherUtilities::CheckSubdomain(Geometry<Node<3> >& rGeometry)
  {

    KRATOS_TRY

    std::string DomainName = rGeometry[0].GetValue(MODEL_PART_NAME); //MODEL_PART_NAME must be set as nodal variable

    int samesbd=0;

    const unsigned int size = rGeometry.size();

    for(unsigned int i=0; i<size; ++i)
      {
	if(DomainName!=rGeometry[i].GetValue(MODEL_PART_NAME))
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

  double MesherUtilities::ComputeModelPartVolume(ModelPart& rModelPart)
  {
    KRATOS_TRY

    const unsigned int dimension = rModelPart.GetProcessInfo()[SPACE_DIMENSION];
    double ModelPartVolume = 0;
    if( dimension == 2 ){

      for(auto& i_elem : rModelPart.Elements())
      {
        if( i_elem.GetGeometry().Dimension() == 2 )
          ModelPartVolume += i_elem.GetGeometry().Area();
      }
    }
    else{ //dimension == 3

      for(auto& i_elem : rModelPart.Elements())
	{
	  if( i_elem.GetGeometry().Dimension() == 3 )
	    ModelPartVolume += i_elem.GetGeometry().Volume();
	}
     }

    return ModelPartVolume;

    KRATOS_CATCH(" ")

  }


  //*******************************************************************************************
  //*******************************************************************************************

  bool MesherUtilities::CheckRigidOuterCentre(Geometry<Node<3> >& rGeometry)
  {

    KRATOS_TRY

    bool outer = false;

    unsigned int RigidNodes = 0;
    const unsigned int size = rGeometry.size();

    for(unsigned int i = 0; i < size; ++i)
      {
	if(rGeometry[i].Is(RIGID))
	  {
	    RigidNodes += 1;
	  }
      }


    if(RigidNodes >= size-1)
    {

      //Baricenter
      array_1d<double, 3>  Center;
      Center.clear();
      array_1d<double, 3>  Normal;

      std::vector<array_1d<double, 3> > Vertices;
      array_1d<double, 3>  Vertex;


      for(unsigned int i = 0; i < size; ++i)
      {
        Vertex  = rGeometry[i].Coordinates();

        Vertices.push_back(Vertex);

        Center += Vertex;
      }

      Center /= (double)size;

      array_1d<double, 3> Corner;

      double tolerance = 0.05;
      int numouter     = 0;

      int numnodes = 0;
      for(unsigned int i = 0; i < size; ++i)
      {
        if(rGeometry[i].Is(RIGID)){

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
            ++numouter;
          }
          ++numnodes;
        }
      }

      if(RigidNodes == size){
        if(numouter > 0)
          outer = true;
      }
      else if(RigidNodes == size-1){
        if(numouter = numouter )
          outer = true;
      }

    }

    return outer; //if is outside the body

    KRATOS_CATCH( "" )
  }

  //*******************************************************************************************
  //*******************************************************************************************

  bool MesherUtilities::CheckInnerCentre(Geometry<Node<3> >& rGeometry)
  {

    KRATOS_TRY

    bool inner = true;

    unsigned int BoundaryNodes = 0;
    const unsigned int size = rGeometry.size();

    for(unsigned int i = 0; i < size; ++i)
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
      noalias(Center) = ZeroVector(3);
      array_1d<double, 3>  Normal;

      std::vector<array_1d<double, 3> > Vertices;
      array_1d<double, 3>  Vertex;


      for(unsigned int i = 0; i < size; ++i)
      {
        Vertex  = rGeometry[i].Coordinates();

        Vertices.push_back(Vertex);

        Center += Vertex;
      }


      Center /= (double)size;

	array_1d<double, 3> Corner;

	double tolerance = 0.05;
	int numouter     = 0;


	for(unsigned int i = 0; i < size; ++i)
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

  bool MesherUtilities::CheckOuterCentre(Geometry<Node<3> >& rGeometry, double& rOffsetFactor, bool& rSelfContact)
  {
    KRATOS_TRY

    bool outer=false;

    unsigned int BoundaryNodes = 0;
    const unsigned int size = rGeometry.size();

    for(unsigned int i = 0; i < size; ++i)
      {
	if(rGeometry[i].Is(BOUNDARY))
	  {
	    BoundaryNodes += 1;
	  }
      }


    if(BoundaryNodes == size){

	//Baricenter
	array_1d<double, 3>  Center;
	Center.clear();
	std::vector<array_1d<double, 3> > Vertices;
	array_1d<double, 3>  Vertex;
	array_1d<double, 3>  Normal;
	double Shrink = 0;

	for(unsigned int i = 0; i < size; ++i)
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

	for(unsigned int i = 0; i < size; ++i)
	{

	    //std::cout<<" V: ("<<rGeometry[i].Id()<<"): ["<<rGeometry[i].X()<<", "<<rGeometry[i].Y()<<", "<<rGeometry[i].Z()<<"]  Normal:"<<rGeometry[i].FastGetSolutionStepValue(NORMAL)<<std::endl;

	    Normal = rGeometry[i].FastGetSolutionStepValue(NORMAL);

	    double NormNormal = norm_2(Normal);
	    if( NormNormal != 0)
		Normal /= NormNormal;


	    //change position to be the vector from the vertex to the geometry center
	    Corner = Center-Vertices[i];

	    if(norm_2(Corner))
		Corner/= norm_2(Corner);

	    double projection = inner_prod(Corner,Normal);

	    if( projection > 0 ){

		if( projection < slope)
		{
		    numouter++;
		}
		else
		{
		    if( projection < extra)
			numextra++;
		}
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

	//std::cout<<std::endl;
	//std::cout<<"  [ no:"<<numouter<<";ne:"<<numextra<<";nc:"<<numcoplanar<<";ns: "<<numsamedirection<<";nor:"<<numorthogonal<<"] ACCEPTED: "<<outer<<std::endl;
    }
    else{
	//std::cout<<" No boundary Element "<<BoundaryNodes<<std::endl;
    }

    return outer; //if is outside the body domain returns true

    KRATOS_CATCH( "" )
  }

    //*******************************************************************************************
  //*******************************************************************************************

  MesherUtilities::ContactElementType MesherUtilities::CheckContactElement(Geometry<Node<3> >& rGeometry, std::vector<int>& rSlaveVertices)
  {

    KRATOS_TRY

    const unsigned int  size = rGeometry.size();

    //Identify subdomains: (non selfcontact elements)
    for(unsigned int i=0; i<size; ++i)
      if(rGeometry[i].IsNot(BOUNDARY))
	return MesherUtilities::NonContact;

    unsigned int NumberOfSlaves = 0;

    if(rSlaveVertices.size() != size)
      rSlaveVertices.resize(size);

    std::fill(rSlaveVertices.begin(),rSlaveVertices.end(),0);


    //Identify subdomains: (non selfcontact elements)
    for(unsigned int i=0; i<size; ++i)
      {
	for(unsigned int j=i+1; j<size; ++j)
	  {
	    if( rGeometry[i].GetValue(MODEL_PART_NAME) == rGeometry[j].GetValue(MODEL_PART_NAME) )
	      {
		//std::cout<<" MP name "<<rGeometry[i].GetValue(MODEL_PART_NAME)<<" "<<rGeometry[j].GetValue(MODEL_PART_NAME)<<std::endl;
		rSlaveVertices[i]+=1;
		rSlaveVertices[j]+=1;
	      }
	  }

	NumberOfSlaves+=rSlaveVertices[i];
      }

    //NonContact Elements or Selfcontact elements (2D/3D): Number of Slaves = size * (size-1);
    if(NumberOfSlaves == size*(size-1)){

      std::vector<int> NeighbourVertices(size);
      std::fill(NeighbourVertices.begin(),NeighbourVertices.end(),0);
      unsigned int NumberOfNeighbours = 0;

      for(unsigned int i=0; i<size; ++i)
	{

	  if( rGeometry[i].Is(NEW_ENTITY) )
	    return MesherUtilities::Undefined;

	  NodeWeakPtrVectorType& nNodes = rGeometry[i].GetValue(NEIGHBOUR_NODES);

	  for(auto& i_nnode : nNodes)
          {
            for(unsigned int j=i+1; j<size; ++j)
            {
              if( i_nnode.Id() == rGeometry[j].Id() )
              {
                NeighbourVertices[i] +=1;
                NeighbourVertices[j] +=1;
              }
            }
          }

	  NumberOfNeighbours += NeighbourVertices[i];
	}

      //Node to face elements (2D/3D): Number of Slaves = (size-1) * (size-2);
      if(NumberOfNeighbours == (size-1)*(size-2))
	return MesherUtilities::PointToFace;
      //Edge to edge elements (3D): Number of Slaves = size;
      if(NumberOfNeighbours == size)
	return MesherUtilities::EdgeToEdge;

      return MesherUtilities::NonContact;
    }

    //Node to face elements (2D/3D): Number of Slaves = (size-1) * (size-2);
    if(NumberOfSlaves == (size-1)*(size-2))
      return MesherUtilities::PointToFace;
    //Edge to edge elements (3D): Number of Slaves = size;
    if(NumberOfSlaves == size)
      return MesherUtilities::EdgeToEdge;
    //Elements with one vertex to each domain(2D/3D): Number of Slaves = 0
    if(NumberOfSlaves == 0)
      return MesherUtilities::PointToPoint;


    return MesherUtilities::NonContact;

    KRATOS_CATCH( "" )

  }

  //*******************************************************************************************
  //*******************************************************************************************

  bool MesherUtilities::CheckSliver(Geometry<Node<3> >& rGeometry)
  {

    KRATOS_TRY

    Vector VectorZero(3);
    noalias(VectorZero) = ZeroVector(3);

    std::vector<Vector> FaceNormals(rGeometry.FacesNumber());
    std::fill(FaceNormals.begin(),FaceNormals.end(),VectorZero);

    std::vector<double> FaceAreas(rGeometry.FacesNumber());
    std::fill(FaceAreas.begin(),FaceAreas.end(), 0.0 );

    double MaximumFaceArea = std::numeric_limits<double>::min();
    double MinimumFaceArea = std::numeric_limits<double>::max();

    DenseMatrix<unsigned int> lpofa;  //points that define the faces
    rGeometry.NodesInFaces(lpofa);
    DenseVector<unsigned int> lnofa;  //number of nodes per face (3)
    rGeometry.NumberNodesInFaces(lnofa);

    //calculate face normals
    Vector FirstVectorPlane(3);
    noalias(FirstVectorPlane) = VectorZero;
    Vector SecondVectorPlane(3);
    noalias(SecondVectorPlane)= VectorZero;

    double FaceArea = 0;
    for(unsigned int i=0; i<rGeometry.FacesNumber(); ++i)
      {
	std::vector<Vector> FaceCoordinates(lnofa[i]); // (3)
	std::fill(FaceCoordinates.begin(),FaceCoordinates.end(),ZeroVector(3));
	for(unsigned int j=0; j<lnofa[i]; ++j)
	  {
	    for(unsigned int k=0; k<3; ++k)
	      {
		FaceCoordinates[j][k] = rGeometry[lpofa(j+1,i)].Coordinates()[k];
	      }
	  }

	if( lnofa[i] > 2 ){
	  FirstVectorPlane  =  FaceCoordinates[1]-FaceCoordinates.front();
	  SecondVectorPlane =  FaceCoordinates.back()-FaceCoordinates.front();
	}
	else{
	  KRATOS_THROW_ERROR( std::logic_error,"2D check sliver not implemented", "" )
	}

	FaceNormals[i] = MathUtils<double>::CrossProduct(FirstVectorPlane,SecondVectorPlane);

	FaceArea = norm_2(FaceNormals[i]);

	if(FaceArea<MinimumFaceArea)
	  MinimumFaceArea = FaceArea;
	if(FaceArea>MaximumFaceArea)
	  MaximumFaceArea = FaceArea;
      }

    //check areas
    if( MaximumFaceArea >= MinimumFaceArea * 1.0e2 ){
      return true;
    }

    //normalize normals
    for(unsigned int i=0; i<FaceNormals.size(); ++i)
      if(norm_2(FaceNormals[i])!=0)
	FaceNormals[i]/=norm_2(FaceNormals[i]);

    //check coincident normals
    std::vector<int> FaceCoincidentNormals(rGeometry.FacesNumber());
    std::fill(FaceCoincidentNormals.begin(),FaceCoincidentNormals.end(), 0 );

    unsigned int CoincidentNormals = 0;
    for(unsigned int i=0; i<lpofa.size2(); ++i)
      {
	for(unsigned int j=i+1; j<lpofa.size2(); ++j)
	  {
	    double projection = inner_prod(FaceNormals[i],FaceNormals[j]);
	    if( fabs(projection) >= 0.99 ){
	      FaceCoincidentNormals[i] +=1;
	      FaceCoincidentNormals[j] +=1;
	    }
	  }

	CoincidentNormals +=FaceCoincidentNormals[i];
      }

    //CoincidentNormals: 12 => all-faces coincident
    //CoincidentNormals:  6 => 3-faces coincident
    //CoincidentNormals:  4 => 2-faces/2-faces coincident
    //CoincidentNormals:  2 => 2-faces coincident

    const unsigned int  size = rGeometry.size();
    unsigned int NumberOfBoundaryNodes = 0;
    for(unsigned int i=0; i<size; ++i)
      if(rGeometry[i].Is(BOUNDARY))
	NumberOfBoundaryNodes+=1;

    // for(unsigned int i=0; i<FaceNormals.size(); ++i)
    //   {
    //     std::cout<<"FaceNormal ["<<i<<"] "<<FaceNormals[i]<<std::endl;
    //   }

    if( NumberOfBoundaryNodes == size ){ //boundary elements
      if( CoincidentNormals >= 4 ){
	return true;
      }
    }
    else if(NumberOfBoundaryNodes >=2){
      if( CoincidentNormals >= 6 ){
	return true;
      }
    }
    else{ //inside elements
      if( CoincidentNormals >= 12 ){
     	return true;
      }
    }

    return false;

    KRATOS_CATCH( "" )

  }


  //*******************************************************************************************
  //*******************************************************************************************

  double MesherUtilities::GetAndCompareSideLenghts(Geometry<Node<3> >& rGeometry, double& rMaximumSideLength, double& rMinimumSideLength)
  {

    KRATOS_TRY

    rMaximumSideLength = std::numeric_limits<double>::min();
    rMinimumSideLength = std::numeric_limits<double>::max();

    DenseMatrix<unsigned int> lpofa;
    rGeometry.NodesInFaces(lpofa);

    double SideLength = 0;
    for(unsigned int i=0; i<lpofa.size2(); ++i)
      {

	for(unsigned int j=1; j<lpofa.size1(); ++j)
	  {
	    SideLength = norm_2(rGeometry[lpofa(0,i)].Coordinates() - rGeometry[lpofa(j,i)].Coordinates());

	    if( SideLength < rMinimumSideLength )
	      rMinimumSideLength = SideLength;

	    if( SideLength > rMaximumSideLength )
	      rMaximumSideLength = SideLength;
	  }
      }


    return (rMaximumSideLength/rMinimumSideLength);

    KRATOS_CATCH( "" )

  }


  //*******************************************************************************************
  //*******************************************************************************************

  bool MesherUtilities::CheckGeometryShape(Geometry<Node<3> >& rGeometry, int& rShape)
  {
    KRATOS_TRY

    bool sliver    = false;
    bool distorted = false;

    const unsigned int  size = rGeometry.size();

    //check if is a sliver and the distorsion of the sliver

    double Volume = 0;
    double MaximumSideLength = 0;
    double MinimumSideLength = 0;

    double CriticalRelativeSideLength = (double)size*5; //edge relative length (3,4)

    Volume = rGeometry.Volume();

    //compare side lengths
    double RelativeSideLength = 0;
    RelativeSideLength = GetAndCompareSideLenghts(rGeometry, MaximumSideLength, MinimumSideLength);

    if( RelativeSideLength > CriticalRelativeSideLength ){
      //std::cout<<" RelativeSideLength "<<RelativeSideLength<<std::endl;
      distorted = true;
    }

    double CriticalVolume =  1e-12 * pow( MinimumSideLength, size-1 );

    //check sliver (volume)
    if( Volume < CriticalVolume ){
      //std::cout<<" Sliver (volume) "<<Volume<<" "<<CriticalVolume<<std::endl;
      sliver = true;
    }

    //check sliver (volume + normals and faces)
    if( !sliver ){
      sliver = CheckSliver(rGeometry);
      // if(sliver)
      // 	std::cout<<" Sliver (faces) "<<sliver<<std::endl;
    }



    //check if it is a contact element (contact domain definition)
    std::vector<int> SlaveVertices;
    ContactElementType ContactType= CheckContactElement(rGeometry, SlaveVertices);

    if( ContactType != NonContact ){

      //std::cout<<" contact type "<<std::endl;

      if( ContactType == PointToFace ){ //POINT_FACE

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
      else if( ContactType == EdgeToEdge ){ //EDGE_EDGE

	//compare vertex normals (detect coplanar faces and orthogonal faces)
	//if( !CheckVertexNormals(rGeometry) )
	distorted = true;
      }
      else if( ContactType == PointToPoint ){ //POINT_POINT

	//compare vertex normals (detect coplanar faces and orthogonal faces)
	//if( !CheckVertexNormals(rGeometry) )
	distorted = true;
      }

    }

    //std::cout<<" DISTORTED "<<distorted<<" SLIVER "<<sliver<<std::endl;

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

  double MesherUtilities::FindBoundaryH (Node<3>& BoundaryPoint)
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

        NodeWeakPtrVectorType& nNodes = BoundaryPoint.GetValue(NEIGHBOUR_NODES);
        for(auto& i_nnode : nNodes)
	  {
	    if( i_nnode.Is(BOUNDARY) ){
	      double x = i_nnode.X();
	      double y = i_nnode.Y();
	      double z = i_nnode.Z();
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
  double& MesherUtilities::ComputeRadius(double& rRadius, double& rVolume, std::vector<Vector >& rVertices,const unsigned int& dimension)
  {
    KRATOS_TRY

    if( dimension == 2 ){

      BoundedMatrix<double,2,2> mJ;    //local jacobian
      BoundedMatrix<double,2,2> mJinv; //inverse jacobian

      //calculation of the jacobian  //coordinate center point 0
      for(unsigned int i = 0; i < dimension; ++i)
	{
	  for(unsigned int j = 0; j < dimension; ++j)
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
      Vector Center = ZeroVector(2);    //center pos

      //center point 0
      for(unsigned int i = 0; i < dimension; ++i)
	{
	  for(unsigned int j = 0; j < dimension; ++j)
	    {
	      Center[i]  += (rVertices[i+1][j] * rVertices[i+1][j]);
	      Center[i]	 -= (rVertices[0][j]   * rVertices[0][j]  );
	    }
	}

      Center = prod(mJinv,Center);

      //calculate the element radius
      Center[0] -= rVertices[0][0];
      Center[1] -= rVertices[0][1];

      rRadius = norm_2(Center);

    }
    else if( dimension == 3 ){

      BoundedVector<double,3>   mRHS;  //center pos
      BoundedMatrix<double,3,3> mJ;    //local jacobian
      BoundedMatrix<double,3,3> mJinv; //inverse jacobian


      //calculation of the jacobian  //coordinate center point 0
      for(unsigned int i = 0; i < dimension; ++i)
	{
	  for(unsigned int j = 0; j < dimension; ++j)
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
  bool MesherUtilities::AlphaShape(double AlphaParameter, Geometry<Node<3> >& rGeometry, const unsigned int dimension, const double MeanMeshSize)
  {
    KRATOS_TRY

    const unsigned int size      = rGeometry.size();

    //calculate geometry radius and volume
    double Radius = 0;
    double Volume = 0;
    std::vector<Vector > Vertices;
    Vector Vertex = ZeroVector(3);
    for(unsigned int i = 0; i < size; ++i)
      {
	Vertex[0] = rGeometry[i].X();
	Vertex[1] = rGeometry[i].Y();
	Vertex[2] = rGeometry[i].Z();

	Vertices.push_back(Vertex);
      }


    Radius = ComputeRadius(Radius, Volume, Vertices, dimension);

    // double CriticalVolume = 1e-12 * pow(h, size-1);
    double AlphaRadius    = AlphaParameter * MeanMeshSize;

    if( Radius<0 ) //degenerated element
      {
	std::cout<<" Sliver (radius) "<<Radius<<" (alpha_volume) "<<Volume<<std::endl;
	return false;
      }
    else
      {
	if( Radius < AlphaRadius )
	  {
	    // std::cout<<"  ACCEPTED!"<<std::endl;
	    return true;
	  }
	else
	  {
	    //std::cout<<"  ERASED! "<<Radius<<" < "<<AlphaRadius<<" MeanMeshSize "<<MeanMeshSize<<" Alpha "<<AlphaParameter<<std::endl;
	    return false;
	  }
      }


    KRATOS_CATCH( "" )
  }


  //returns false if it should be removed
  bool MesherUtilities::AlphaShape(double AlphaParameter, Geometry<Node<3> >& rGeometry, const unsigned int dimension)
  {
    KRATOS_TRY

    const unsigned int size      = rGeometry.size();

    //calculate geometry radius and volume
    double Radius = 0;
    double Volume = 0;
    std::vector<Vector > Vertices;
    Vector Vertex = ZeroVector(3);
    for(unsigned int i = 0; i < size; ++i)
      {
	Vertex[0] = rGeometry[i].X();
	Vertex[1] = rGeometry[i].Y();
	Vertex[2] = rGeometry[i].Z();

	Vertices.push_back(Vertex);
      }

    Radius = ComputeRadius(Radius, Volume, Vertices, dimension);

    //calculate average h
    double h = 0;

    for( unsigned int i = 0; i<size; ++i )
      {
	h += rGeometry[i].FastGetSolutionStepValue(NODAL_H);
      }

    h /= (double)size;

    double CriticalVolume = 1e-12 * pow(h, size-1);
    double AlphaRadius    = AlphaParameter * h;

    if( Volume < CriticalVolume ) //sliver
      {
	//std::cout<<" Sliver (alpha_volume) "<<Volume<<" "<<CriticalVolume<<std::endl;
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
  bool MesherUtilities::ShrankAlphaShape(double AlphaParameter, Geometry<Node<3> >& rGeometry,double& rOffsetFactor, const unsigned int dimension)
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

    for(unsigned int i = 0; i < size; ++i)
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

    for( unsigned int i = 0; i<size; ++i )
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

  void MesherUtilities::CheckParticles (ModelPart& rModelPart)
  {
    KRATOS_TRY

    int NumberOfNodes = rModelPart.NumberOfNodes();
    std::cout<<" Number of Nodes "<<NumberOfNodes<<std::endl;
    for(int id=1; id<=NumberOfNodes; ++id)
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

  bool MesherUtilities::CheckRelativeVelocities(Geometry<Node<3> >& rVertices, const double& rRelativeFactor)
  {
    KRATOS_TRY

    const unsigned int NumberOfVertices = rVertices.size();

    std::vector<double> VelocityModulus(NumberOfVertices);

    for(unsigned int i = 0; i<NumberOfVertices; ++i)
    {
      VelocityModulus[i]=norm_2(rVertices[i].FastGetSolutionStepValue(VELOCITY));
    }

    for(unsigned int i = 0; i<NumberOfVertices-1; ++i)
    {
      for(unsigned int j = i+1; j<NumberOfVertices; ++j)
      {
        if( VelocityModulus[i]/VelocityModulus[j]>rRelativeFactor || VelocityModulus[j]/VelocityModulus[i]>rRelativeFactor ){
          return true;
        }
      }
    }

    return false;

    KRATOS_CATCH( "" )
  }

  //*******************************************************************************************
  //*******************************************************************************************

  bool MesherUtilities::CheckVolumeDecrease(GeometryType& rVertices, const unsigned int& rDimension,const double& rTolerance, double& VolumeChange)
  {
    bool accepted = false;
    if(rDimension==2){
      Triangle2D3<Node<3> > CurrentTriangle(rVertices);
      double CurrentArea = CurrentTriangle.Area();

      //new volume with a 1.0 * DeltaDisplacement
      double MovedArea = GetMovedVolume(rVertices,rDimension,1.0);

      //std::cout<<" control fluid  "<<MovedArea<<" "<<CurrentArea<<std::endl;
      VolumeChange = CurrentArea-MovedArea;

      if(MovedArea+rTolerance<CurrentArea){
        accepted = true;
      }

    }
    else if(rDimension==3){

      Tetrahedra3D4<Node<3> > CurrentTetrahedron(rVertices);
      double CurrentVolume = CurrentTetrahedron.Volume();

      //new volume with a 1.0 * DeltaDisplacement
      double MovedVolume = GetMovedVolume(rVertices,rDimension,1.0);

      //std::cout<<" control fluid  "<<MovedVolume<<" "<<CurrentVolume<<std::endl;
      VolumeChange = CurrentVolume-MovedVolume;

      if(MovedVolume+rTolerance<CurrentVolume){
        accepted = true;
      }
    }

    return accepted;
  }

  //*******************************************************************************************
  //*******************************************************************************************

  double MesherUtilities::GetMovedVolume(GeometryType& rVertices, const unsigned int& rDimension, double MovementFactor)
  {
    double MovedVolume = 0.0;
    if(rDimension==2){

      //Triangle geometry
      array_1d<double,3> P0;
      noalias(P0) = rVertices[0].Coordinates() + MovementFactor * (rVertices[0].FastGetSolutionStepValue(DISPLACEMENT) - rVertices[0].FastGetSolutionStepValue(DISPLACEMENT,1));
      array_1d<double,3> P1;
      noalias(P1) = rVertices[1].Coordinates() + MovementFactor * (rVertices[1].FastGetSolutionStepValue(DISPLACEMENT) - rVertices[1].FastGetSolutionStepValue(DISPLACEMENT,1));

      double x10 = P1[0] - P0[0];
      double y10 = P1[1] - P0[1];

      noalias(P1) = rVertices[2].Coordinates() + MovementFactor * (rVertices[2].FastGetSolutionStepValue(DISPLACEMENT) - rVertices[2].FastGetSolutionStepValue(DISPLACEMENT,1));

      double x20 = P1[0] - P0[0];
      double y20 = P1[1] - P0[1];

      MovedVolume = 0.5 * ( x10 * y20 - y10 * x20 );

    }
    else if(rDimension==3){

      //Tetrahedron geometry
      const double onesixth = 1.0/6.0;

      array_1d<double,3> P0;
      noalias(P0) = rVertices[0].Coordinates() + MovementFactor * (rVertices[0].FastGetSolutionStepValue(DISPLACEMENT) - rVertices[0].FastGetSolutionStepValue(DISPLACEMENT,1));

      array_1d<double,3> P1;
      noalias(P1) = rVertices[1].Coordinates() + MovementFactor * (rVertices[1].FastGetSolutionStepValue(DISPLACEMENT) - rVertices[1].FastGetSolutionStepValue(DISPLACEMENT,1));

      double x10 = P1[0] - P0[0];
      double y10 = P1[1] - P0[1];
      double z10 = P1[2] - P0[2];

      noalias(P1) = rVertices[2].Coordinates() + MovementFactor * (rVertices[2].FastGetSolutionStepValue(DISPLACEMENT) - rVertices[2].FastGetSolutionStepValue(DISPLACEMENT,1));

      double x20 = P1[0] - P0[0];
      double y20 = P1[1] - P0[1];
      double z20 = P1[2] - P0[2];

     noalias(P1)  = rVertices[3].Coordinates() + MovementFactor * (rVertices[3].FastGetSolutionStepValue(DISPLACEMENT) - rVertices[3].FastGetSolutionStepValue(DISPLACEMENT,1));

      double x30 = P1[0] - P0[0];
      double y30 = P1[1] - P0[1];
      double z30 = P1[2] - P0[2];

      MovedVolume = onesixth * (x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30);

      // if(MovedVolume<0)
      //   std::cout<<" VOLUME negative "<<std::endl;
    }

    return MovedVolume;
  }


  //*******************************************************************************************
  //*******************************************************************************************

  double MesherUtilities::GetDeformationGradientDeterminant(GeometryType& rVertices, const unsigned int& rDimension)
  {
    //Deformation Gradient determinant
    unsigned int number_of_nodes = rVertices.size();

    //Configuration increment
    Matrix DeltaPosition(number_of_nodes,rDimension);
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        const array_1d<double, 3 > & CurrentDisplacement  = rVertices[i].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3 > & PreviousDisplacement = rVertices[i].FastGetSolutionStepValue(DISPLACEMENT,1);

        for ( unsigned int j = 0; j < rDimension; j++ )
        {
          DeltaPosition(i,j) = CurrentDisplacement[j]-PreviousDisplacement[j];
        }
    }

    //Compute cartesian derivatives [dN/dx_n]
    Matrix DN_DX;

    if(rDimension==2){
      Triangle2D3<Node<3> > Triangle(rVertices);

      DenseVector<Matrix> J;
      J = Triangle.Jacobian( J, GeometryData::GI_GAUSS_1, DeltaPosition );

      //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
      Matrix InvJ;
      double detJ;

      MathUtils<double>::InvertMatrix2( J[0], InvJ, detJ);

      const Matrix& DN_De = Triangle.ShapeFunctionLocalGradient(0,GeometryData::GI_GAUSS_1);
      DN_DX = prod( DN_De, InvJ );

    }
    else if(rDimension==3){

      Tetrahedra3D4<Node<3> > Tetrahedron(rVertices);

      DenseVector<Matrix> J;
      J = Tetrahedron.Jacobian( J, GeometryData::GI_GAUSS_1, DeltaPosition );

      //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
      Matrix InvJ;
      double detJ;
      MathUtils<double>::InvertMatrix3( J[0], InvJ, detJ);

      const Matrix& DN_De = Tetrahedron.ShapeFunctionLocalGradient(0,GeometryData::GI_GAUSS_1);
      DN_DX = prod( DN_De, InvJ );
    }


    Matrix F(rDimension,rDimension);
    noalias(F) = ZeroMatrix(rDimension,rDimension);
    for (unsigned int i = 0; i < rDimension; i++)
    {
      for (unsigned int j = 0; j < rDimension; j++)
      {
        for (unsigned int k = 0; k < number_of_nodes; k++)
        {
          F(i,j)+= rVertices[k].Coordinates()[i]*DN_DX(k,j);
        }
      }
    }
    double detF  = MathUtils<double>::Det(F);

    if(detF<0)
      std::cout<<" NEGATIVE ELEMENT (DET_F: "<<detF<<")"<<std::endl;

    return detF;
  }


  //*******************************************************************************************
  //*******************************************************************************************

  bool MesherUtilities::CheckConditionInBox(Condition::Pointer& pCondition, SpatialBoundingBox& rRefiningBox, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    bool inside = true;
    Vector Point(3);
    Geometry< Node<3> >& rGeometry = pCondition->GetGeometry();

    for(unsigned int i=0; i<rGeometry.size(); ++i)
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

  bool MesherUtilities::CheckElementInBox(Element::Pointer& pElement, SpatialBoundingBox& rRefiningBox, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    bool inside = true;
    Vector Point(3);
    Geometry< Node<3> >& rGeometry = pElement->GetGeometry();

    for(unsigned int i=0; i<rGeometry.size(); ++i)
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

  bool MesherUtilities::CheckVerticesInBox(Geometry<Node<3> >& rGeometry, SpatialBoundingBox& rRefiningBox, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    bool inside = true;
    Vector Point(3);

    for(unsigned int i=0; i<rGeometry.size(); ++i)
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

  Condition::Pointer MesherUtilities::FindMasterCondition(Condition::Pointer& pCondition, ModelPart::ConditionsContainerType & rConditions,bool & condition_found)
  {
    KRATOS_TRY

    Condition::Pointer pMasterCondition;

    Geometry< Node<3> >& rGeometry = pCondition->GetGeometry();
    DenseMatrix<unsigned int> lpofa; //points that define the faces
    rGeometry.NodesInFaces(lpofa);

    //std::cout<<" lpofa "<<lpofa<<std::endl;
    //std::cout<<" rGeometry "<<rGeometry<<std::endl;

    unsigned int face_elements = 0;
    unsigned int edge_elements = 0;


    if( rGeometry.size() == 3 )
    { //triangles of 3 nodes

      condition_found=false;
      for(auto i_cond(rConditions.begin()); i_cond != rConditions.end(); ++i_cond)
      {

        //2D edges:
        if(i_cond->IsNot(CONTACT)){

          Geometry<Node<3> >& rConditionGeometry = i_cond->GetGeometry();

          for(unsigned int iface=0; iface<lpofa.size2(); ++iface)
          {
            if( (   rConditionGeometry[0].Id() == rGeometry[lpofa(1,iface)].Id()
                    && rConditionGeometry[1].Id() == rGeometry[lpofa(2,iface)].Id() ) ||
                (   rConditionGeometry[0].Id() == rGeometry[lpofa(2,iface)].Id()
                    && rConditionGeometry[1].Id() == rGeometry[lpofa(1,iface)].Id() ) )
            {
              pMasterCondition=*i_cond.base();
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
    }
    else if( rGeometry.size() == 4 )
    { //tetraheda of 4 nodes

      face_elements = 0;
      edge_elements = 0;

      condition_found=false;
      for(auto i_cond(rConditions.begin()); i_cond != rConditions.end(); ++i_cond)
      {

        //3D faces:
        if(i_cond->IsNot(CONTACT)){

          Geometry< Node<3> >& rConditionGeometry = i_cond->GetGeometry();

          for(unsigned int iface=0; iface<lpofa.size2(); ++iface)
          {
            //detection for contact elements clockwise numeration of the contact geometry.
            if( (   rConditionGeometry[2].Id() == rGeometry[lpofa(1,iface)].Id()
                    && rConditionGeometry[1].Id() == rGeometry[lpofa(2,iface)].Id()
                    && rConditionGeometry[0].Id() == rGeometry[lpofa(3,iface)].Id() ) ||
                (   rConditionGeometry[2].Id() == rGeometry[lpofa(3,iface)].Id()
                    && rConditionGeometry[1].Id() == rGeometry[lpofa(1,iface)].Id()
                    && rConditionGeometry[0].Id() == rGeometry[lpofa(2,iface)].Id() ) ||
                (   rConditionGeometry[2].Id() == rGeometry[lpofa(2,iface)].Id()
                    && rConditionGeometry[1].Id() == rGeometry[lpofa(3,iface)].Id()
                    && rConditionGeometry[0].Id() == rGeometry[lpofa(1,iface)].Id() ) )
            {
              pMasterCondition = *i_cond.base();
              condition_found = true;
              break;
            }

          }

        }

        if(condition_found)
        {
          pCondition->Set(SELECTED.AsFalse()); //meaning that is a element that shares faces
          face_elements++;
          break;
        }

      }

      if(!condition_found) {

	//check if it is EDGE_TO_EDGE element sharing only edges with the conditions

	condition_found=false;
        for(auto i_cond(rConditions.begin()); i_cond != rConditions.end(); ++i_cond)
        {

          //3D edges: there are 4 possibilities, it takes the first one that matches
          if(i_cond->IsNot(CONTACT)){

            Geometry< Node<3> >& rConditionGeometry = i_cond->GetGeometry();

            for(unsigned int iface=0; iface<lpofa.size2()-1; ++iface)
            {
              if( (   rConditionGeometry[0].Id() == rGeometry[lpofa(1,iface)].Id()
                      && rConditionGeometry[1].Id() == rGeometry[lpofa(2,iface)].Id() )||
                  (   rConditionGeometry[1].Id() == rGeometry[lpofa(1,iface)].Id()
                      && rConditionGeometry[2].Id() == rGeometry[lpofa(2,iface)].Id() )||
                  (   rConditionGeometry[2].Id() == rGeometry[lpofa(1,iface)].Id()
                      && rConditionGeometry[0].Id() == rGeometry[lpofa(2,iface)].Id() )||

                  (   rConditionGeometry[0].Id() == rGeometry[lpofa(2,iface)].Id()
                      && rConditionGeometry[1].Id() == rGeometry[lpofa(3,iface)].Id() )||
                  (   rConditionGeometry[1].Id() == rGeometry[lpofa(2,iface)].Id()
                      && rConditionGeometry[2].Id() == rGeometry[lpofa(3,iface)].Id() )||
                  (   rConditionGeometry[2].Id() == rGeometry[lpofa(2,iface)].Id()
                      && rConditionGeometry[0].Id() == rGeometry[lpofa(3,iface)].Id() ) )
              {
                pMasterCondition= *i_cond.base();
                condition_found=true;
                break;
              }

            }

          }

          if(condition_found)
          {
            pCondition->Set(SELECTED); //meaning that is a element that shares edges instead of faces
            edge_elements++;
            break;
          }

        }

      }

    }


    if(!condition_found)
    {
      std::cout<<" WARNING:: Boundary Condition NOT FOUND after CONTACT MESHING SEARCH "<<std::endl;

      std::cout<<" Condition Nodes[ ";
      for(unsigned int i=0; i<rGeometry.size();++i)
	{
	  std::cout<<" "<<rGeometry[i].Id();
	}
      std::cout<<" ]"<<std::endl;
    }
    // else{
    //   std::cout<<"    [Face Elements: "<<face_elements<<" Edge Elements: "<<edge_elements<<"]"<<std::endl;
    // }

    return pMasterCondition;

    KRATOS_CATCH( "" )
  }

  //*******************************************************************************************
  //*******************************************************************************************

  Condition::Pointer MesherUtilities::FindMasterCondition(Condition::Pointer& pCondition, PointType& pSlaveNode, ModelPart::ConditionsContainerType & rConditions,bool & condition_found)
  {
    KRATOS_TRY

    Condition::Pointer pMasterCondition;

    Geometry< Node<3> >& rGeometry = pCondition->GetGeometry();
    DenseMatrix<unsigned int> lpofa; //points that define the faces
    rGeometry.NodesInFaces(lpofa);

    //std::cout<<" lpofa "<<lpofa<<std::endl;
    //std::cout<<" rGeometry "<<rGeometry<<std::endl;

    condition_found=false;
    for(auto i_cond(rConditions.begin()); i_cond != rConditions.end(); ++i_cond)
    {
      //2D edges:
      if(i_cond->IsNot(CONTACT)){

        Geometry< Node<3> >& rConditionGeom = i_cond->GetGeometry();

        for(unsigned int i=0; i<lpofa.size2();++i)
        {
          // std::cout<<" General Conditions IDs ["<<rConditionGeom[0].Id()<<"] ["<<rConditionGeom[1].Id()<<"] "<<std::endl;
          // std::cout<<" Local Conditions IDs ("<<i<<"):["<<rGeometry[lpofa(1,i)].Id()<<"] ["<<rGeometry[lpofa(2,i)].Id()<<"] "<<std::endl;

          if( (   rConditionGeom[0].Id() == rGeometry[lpofa(1,i)].Id()
                  && rConditionGeom[1].Id() == rGeometry[lpofa(2,i)].Id() ) ||
              (   rConditionGeom[0].Id() == rGeometry[lpofa(2,i)].Id()
                  && rConditionGeom[1].Id() == rGeometry[lpofa(1,i)].Id() ) )
          {
            pMasterCondition = *i_cond.base();
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

  bool MesherUtilities::CheckContactActive(GeometryType& rConditionGeometry, bool& rSemiActiveContact, std::vector<bool>& rSemiActiveNodes)
  {
    KRATOS_TRY

    unsigned int size = rConditionGeometry.size();
    unsigned int counter = 0;

    rSemiActiveContact = false;
    rSemiActiveNodes.resize(size);
    std::fill( rSemiActiveNodes.begin(), rSemiActiveNodes.end(), false );

    for(unsigned int i=0; i<size; ++i){

      bool contact_active = false;
      if( rConditionGeometry[i].SolutionStepsDataHas(CONTACT_FORCE) ){
	array_1d<double, 3 > & ContactForceNormal  = rConditionGeometry[i].FastGetSolutionStepValue(CONTACT_FORCE);

	if(norm_2(ContactForceNormal)>0)
	  contact_active = true;
      }

      if( contact_active ){
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


  bool MesherUtilities::CheckContactCurvature(GeometryType& rConditionGeometry, std::vector<array_1d<double, 3> >& rContactNormals)
  {
    KRATOS_TRY

    unsigned int size = rConditionGeometry.size();
    unsigned int counter = 0;

    array_1d<double, 3> Normal;
    Normal.clear();
    rContactNormals.resize(size);
    std::fill( rContactNormals.begin(), rContactNormals.end(), Normal );

    double modulus = 1.0;
    for(unsigned int i=0; i<size; ++i){

      bool contact_active = false;
      if( rConditionGeometry[i].SolutionStepsDataHas(CONTACT_NORMAL) ){
	array_1d<double, 3 > & ContactNormal  = rConditionGeometry[i].FastGetSolutionStepValue(CONTACT_NORMAL);

	modulus = norm_2(ContactNormal);
	if( modulus )
	  rContactNormals[i] = (1.0/modulus) * ContactNormal;
	else
	  counter++;
      }

    }

    // if no CONTACT_NORMAL is assigned in some condition nodes
    // the curvature can not be evaluated
    if( counter )
      return false;

    counter = 0;
    double tolerance = 0.05;

    for(unsigned int i=1; i<size; ++i){
      modulus = inner_prod(rContactNormals[i-1], rContactNormals[i]);
      if( modulus < 0.95 )
	counter++;
    }
    modulus = inner_prod(rContactNormals[0], rContactNormals[size-1]);
    if( modulus < 0.95 )
      counter++;

    if(counter == size)
      return true;
    else
      return false;

    KRATOS_CATCH( "" )
  }


  //*******************************************************************************************
  //*******************************************************************************************

  double MesherUtilities::CheckCriticalRadius (ModelPart& rModelPart, double rCriticalRadius)
  {
    KRATOS_TRY

    double minimum_h = rCriticalRadius;

    for(auto& i_node : rModelPart.Nodes())
      {
	double nodal_h = i_node.FastGetSolutionStepValue(NODAL_H);
	if( nodal_h < rCriticalRadius )
	  minimum_h = nodal_h;
      }

    //if( minimum_h < rCriticalRadius )
    //  KRATOS_INFO(" CRITICAL MESH SIZE ")<<" supplied size "<<rCriticalRadius<<" is bigger than initial mesh size "<<minimum_h<<" ] "<<std::endl;

    return minimum_h;

    KRATOS_CATCH( "" )

  }

  //**************************************************************************
  //**************************************************************************

  bool MesherUtilities::FindCondition(Geometry< Node<3> >& rConditionGeometry ,Geometry< Node<3> >& rGeometry, DenseMatrix<unsigned int>& lpofa, DenseVector<unsigned int>& lnofa, unsigned int& iface)
  {

    KRATOS_TRY
    // not equivalent geometry sizes for boundary conditions:
    if( rConditionGeometry.size() != lnofa[iface] )
      return false;

    // line boundary condition:
    if( lnofa[iface] == 2 )
      {
	if( (   rConditionGeometry[0].Id() == rGeometry[lpofa(1,iface)].Id()
		&& rConditionGeometry[1].Id() == rGeometry[lpofa(2,iface)].Id() ) ||
	    (   rConditionGeometry[0].Id() == rGeometry[lpofa(2,iface)].Id()
		&& rConditionGeometry[1].Id() == rGeometry[lpofa(1,iface)].Id() ) )
	  {
	    return true;
	  }
	else
	  {
	    return false;
	  }

      }

    //3D faces:
    if(  lnofa[iface] == 3 )
      {
	if( (   rConditionGeometry[0].Id() == rGeometry[lpofa(1,iface)].Id()
		&& rConditionGeometry[1].Id() == rGeometry[lpofa(2,iface)].Id()
		&& rConditionGeometry[2].Id() == rGeometry[lpofa(3,iface)].Id() ) ||
	    (   rConditionGeometry[0].Id() == rGeometry[lpofa(3,iface)].Id()
		&& rConditionGeometry[1].Id() == rGeometry[lpofa(1,iface)].Id()
		&& rConditionGeometry[2].Id() == rGeometry[lpofa(2,iface)].Id() ) ||
	    (   rConditionGeometry[0].Id() == rGeometry[lpofa(2,iface)].Id()
		&& rConditionGeometry[1].Id() == rGeometry[lpofa(3,iface)].Id()
		&& rConditionGeometry[2].Id() == rGeometry[lpofa(1,iface)].Id() ) )
	  {
	    return true;
	  }
	else
	  {
	    return false;
	  }

      }

    if(  lnofa[iface] > 3 )
      {
	KRATOS_THROW_ERROR( std::logic_error, "Wrong Condition Number of Face Nodes",*this );
      }

    return false;

    KRATOS_CATCH( "" )

  }

  //*******************************************************************************************
  //*******************************************************************************************

  void MesherUtilities::SetNodes(ModelPart& rModelPart,
				  MeshingParameters& rMeshingVariables)
  {
    KRATOS_TRY

    const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

    //*********************************************************************
    //input mesh: NODES
    MesherUtilities::MeshContainer& InMesh = rMeshingVariables.InMesh;

    InMesh.CreatePointList(rModelPart.Nodes().size(), dimension);

    double* PointList     = InMesh.GetPointList();
    int& NumberOfPoints   = InMesh.GetNumberOfPoints();

    if(!rMeshingVariables.InputInitializedFlag){

      rMeshingVariables.NodeMaxId = 0;
      if((int)rMeshingVariables.NodalPreIds.size() != NumberOfPoints+1)
	rMeshingVariables.NodalPreIds.resize(NumberOfPoints+1);

      std::fill( rMeshingVariables.NodalPreIds.begin(), rMeshingVariables.NodalPreIds.end(), 0 );
    }

    //writing the points coordinates in a vector and reordening the Id's
    ModelPart::NodesContainerType::iterator nodes_begin = rModelPart.NodesBegin();

    int base   = 0;
    int direct = 1;

    for(int i = 0; i<NumberOfPoints; ++i)
    {
      //from now on it is consecutive
      if(!rMeshingVariables.InputInitializedFlag){
        rMeshingVariables.NodalPreIds[direct]=(nodes_begin + i)->Id();
        (nodes_begin + i)->SetId(direct);
        if( rMeshingVariables.NodalPreIds[direct] > rMeshingVariables.NodeMaxId)
          rMeshingVariables.NodeMaxId = rMeshingVariables.NodalPreIds[direct];
      }

      array_1d<double, 3>& Coordinates = (nodes_begin + i)->Coordinates();

      if(rMeshingVariables.Options.Is(MesherUtilities::CONSTRAINED)){

        if( (nodes_begin + i)->Is(BOUNDARY) ){

          array_1d<double, 3>&  Normal=(nodes_begin + i)->FastGetSolutionStepValue(NORMAL); //BOUNDARY_NORMAL must be set as nodal variable
          double Shrink = (nodes_begin + i)->FastGetSolutionStepValue(SHRINK_FACTOR);   //SHRINK_FACTOR   must be set as nodal variable

          array_1d<double, 3> Offset;

          Normal /= norm_2(Normal);
          for(unsigned int j=0; j<dimension; ++j){
            Offset[j] = ( (-1) * Normal[j] * Shrink * rMeshingVariables.OffsetFactor * 0.25 );
          }

          for(unsigned int j=0; j<dimension; ++j){
            PointList[base+j]   = Coordinates[j] + Offset[j];
          }

          //std::cout<<" Node ["<<(nodes_begin + i)->Id()<<"] "<<Coordinates + Offset<<std::endl;
        }
        else{
          for(unsigned int j=0; j<dimension; ++j){
            PointList[base+j]   = Coordinates[j];
          }

          //std::cout<<" Node ["<<(nodes_begin + i)->Id()<<"] "<<Coordinates<<std::endl;
        }

      }
      else{
        for(unsigned int j=0; j<dimension; ++j){
          PointList[base+j]   = Coordinates[j];
        }
      }

      base+=dimension;
      direct+=1;
    }

    //InMesh.SetPointList(PointList);

    KRATOS_CATCH( "" )

  }


  //*******************************************************************************************
  //*******************************************************************************************

  void MesherUtilities::SetElements(ModelPart& rModelPart,
				     MeshingParameters& rMeshingVariables)
  {
    KRATOS_TRY

    //*********************************************************************
    //input mesh: ELEMENTS
    ModelPart::ElementsContainerType::iterator element_begin = rModelPart.ElementsBegin();

    const unsigned int nds       = element_begin->GetGeometry().size();

    MesherUtilities::MeshContainer& InMesh = rMeshingVariables.InMesh;

    InMesh.CreateElementList(rModelPart.Elements().size(), nds);

    int* ElementList      = InMesh.GetElementList();
    int& NumberOfElements = InMesh.GetNumberOfElements();


    int base=0;
    for(unsigned int el = 0; el<(unsigned int)NumberOfElements; ++el)
    {
      Geometry<Node<3> >& geom = (element_begin+el)->GetGeometry();

      for(unsigned int i=0; i<nds; ++i)
      {
        ElementList[base+i] = geom[i].Id();
      }
      base+=nds;
    }

    KRATOS_CATCH( "" )

  }



} // Namespace Kratos
