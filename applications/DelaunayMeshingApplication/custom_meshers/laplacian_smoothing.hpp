//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined(KRATOS_LAPLACIAN_SMOOTHING_H_INCLUDED)
#define  KRATOS_LAPLACIAN_SMOOTHING_H_INCLUDED

// System includes

// Project includes
#include "includes/variables.h"
#include "spatial_containers/spatial_containers.h"
#include "custom_utilities/mesher_utilities.hpp"

#include "delaunay_meshing_application_variables.h"

#ifdef   SINGLE
#define  REAL float
#else    // not SINGLE
#define  REAL double
#endif   // not SINGLE


#if !defined(KRATOS_TRIANGLE_EXTERNAL_H_INCLUDED)
#define  KRATOS_TRIANGLE_EXTERNAL_H_INCLUDED
#include "triangle.h"
#endif

// External includes
#if !defined(KRATOS_TETGEN_EXTERNAL_H_INCLUDED)
#define  KRATOS_TETGEN_EXTERNAL_H_INCLUDED
#include "tetgen.h"
#endif


//VARIABLES used:
//Data:     NEIGHBOUR_NODES
//StepData: DISPLACEMENT, CONTACT_FORCE, NORMAL, OFFSET
//Flags:    (checked) BOUNDARY, TO_ERASE, INSIDE
//          (set)
//          (modified) INSIDE
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
      for (unsigned int i = 0; i < rModelPart.NumberOfMeshes(); ++i)
	mPreviousMeshes.push_back(Kratos::make_shared<MeshType>((rModelPart.GetMesh(i)).Clone())); //in fact do not clones

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
			    std::vector<int> & PreservedElements,
			    const int* pElementsList,
			    const int& NumberOfPoints)
    {

      KRATOS_TRY

      NodesContainerType& rNodes = rModelPart.Nodes();

      ModelPart::ElementsContainerType::iterator element_begin = rModelPart.ElementsBegin();

      const unsigned int nds = element_begin->GetGeometry().size();

      //*******************************************************************
      //NEIGHBOUR NODES:

      std::vector<int> EmptyVector(0);
      std::vector<std::vector<int> >  NeighborNodesList(rNodes.size());
      std::fill( NeighborNodesList.begin(), NeighborNodesList.end(), EmptyVector );

      this->GetNeighborNodes(NeighborNodesList, PreservedElements, pElementsList, NumberOfPoints, nds);


      //*******************************************************************
      //MOVE NODES: LAPLACIAN SMOOTHING:


      double convergence_tol  = 0.001;
      double smoothing_factor = 0.1;
      double smoothing_iters  = 3;//4
      double iters = 0;

      bool simple = true; //weight = 1;

      bool   converged = false;
      double MaxLength = 0;
      double NewMaxLength = 0;

      bool contact_active = false;
      double boundary_weight = 0.9; //(0,1]
      double contact_weight  = 0.8; //(0,1]

      unsigned int number_of_nodes = 0;

      ModelPart::NodesContainerType::iterator nodes_begin = rModelPart.NodesBegin();

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


       	number_of_nodes = 0;
	for(unsigned int in = 0; in<rNodes.size(); in++)
	  {

	    unsigned int NumberOfNeighbours = NeighborNodesList[in+1].size();

	    if(rNodes[in+1].IsNot(BOUNDARY) && rNodes[in+1].IsNot(RIGID) && rNodes[in+1].IsNot(TO_ERASE) && NumberOfNeighbours>1)
	      {

		TotalDistance.clear();
		TotalWeight = 0;
		Weight = 0;

		//point position
		P[0] = (nodes_begin+in)->X();
		P[1] = (nodes_begin+in)->Y();
		P[2] = (nodes_begin+in)->Z();

		//std::cout<<" Initial Position: "<<P<<std::endl;
		Length = 0;

		for(unsigned int i = 0; i < NumberOfNeighbours; ++i)
		  {
		    //neighbour position
		    Q[0] = (nodes_begin+(NeighborNodesList[in+1][i]-1))->X();
		    Q[1] = (nodes_begin+(NeighborNodesList[in+1][i]-1))->Y();
		    Q[2] = (nodes_begin+(NeighborNodesList[in+1][i]-1))->Z();

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

		    if( (nodes_begin+(NeighborNodesList[in+1][i]-1))->Is(BOUNDARY) ){

		      contact_active = false;
		      if( (nodes_begin+(NeighborNodesList[in+1][i]-1))->SolutionStepsDataHas(CONTACT_FORCE) ){
			array_1d<double, 3 > & ContactForce = (nodes_begin+(NeighborNodesList[in+1][i]-1))->FastGetSolutionStepValue(CONTACT_FORCE);
			if( norm_2(ContactForce) !=0 )
			  contact_active = true;
		      }

		      if( contact_active )
			Weight *= contact_weight;
		      else
			Weight *= boundary_weight;

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

		(nodes_begin+in)->X() = P[0];
		(nodes_begin+in)->Y() = P[1];
		(nodes_begin+in)->Z() = P[2];

		number_of_nodes +=1;
	      }

	  }


	if( (NewMaxLength-MaxLength)/NewMaxLength < convergence_tol ){
	  converged = true;
	  if( GetEchoLevel() > 0 )
	    std::cout<<"   Laplacian smoothing convergence achieved "<<std::endl;
	}


	iters++;

      }

      if(iters==smoothing_iters && !converged){
	if( GetEchoLevel() > 0 )
	  std::cout<<"     WARNING: Laplacian smoothing convergence NOT achieved (iters:"<<iters<<")"<<std::endl;
      }

      //*******************************************************************
      //MOVE NODES: BOUNDARY SMOOTHING
      SetBoundarySmoothing (rModelPart, PreservedElements, pElementsList, NumberOfPoints);


      //*******************************************************************
      //MOVE NODES: BOUNDARY PROJECTION
      //SetInsideProjection (rModelPart, out, NeighborNodesList);

      //*******************************************************************
      //TRANSFER VARIABLES TO NEW NODES POSITION:
      ProjectVariablesToNewPositions(rModelPart);

      KRATOS_CATCH( "" )

    }


    //*******************************************************************************************
    //*******************************************************************************************

    void ProjectVariablesToNewPositions(ModelPart& rModelPart)
    {

      KRATOS_TRY

      bool transfer=true; //transfer active or inactive

      if(transfer==true){

	//create the list of the nodes to be check during the search (new positions after the laplacian smoothing)

	std::vector<Node<3>::Pointer> list_of_nodes;

	for(NodesContainerType::iterator i_node = rModelPart.NodesBegin() ; i_node != rModelPart.NodesEnd() ; i_node++)
	  {

	    //std::cout<<" ID: "<<i_node->Id()<<" NODAL_H "<<i_node->FastGetSolutionStepValue(NODAL_H)<<std::endl;

	    //if(rNodes[in+1].IsNot(BOUNDARY) && rNodes[in+1].IsNot(TO_ERASE) && NumberOfNeighbours>1)
	    if( i_node->IsNot(TO_ERASE) ){
	      (list_of_nodes).push_back(*(i_node.base()));
	    }
	    // else {
	    //   std::cout <<" LLM:PSEUDOERROR node to erase : " << i_node->Id() << std::endl;
	    // }
	  }


	//defintions for spatial search
	typedef Node<3>                                  PointType;
	typedef Node<3>::Pointer                  PointPointerType;
	typedef std::vector<PointPointerType>          PointVector;
	typedef PointVector::iterator                PointIterator;
	//typedef std::vector<double>                 DistanceVector;
	typedef std::vector<double>::iterator     DistanceIterator;

	typedef Bucket<3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > BucketType;
	typedef Tree< KDTreePartition<BucketType> >     KdtreeType; //Kdtree
	//defintions for spatial search

	unsigned int  bucket_size = 20;
	KdtreeType    NodesTree(list_of_nodes.begin(),list_of_nodes.end(),bucket_size);

	//Find out where the new nodes belong to:
	unsigned int number_of_nodes = list_of_nodes.size()+1;

	ModelPart::ElementsContainerType::iterator element_begin = rModelPart.ElementsBegin();
	const unsigned int nds = element_begin->GetGeometry().size();

	std::vector<double> ShapeFunctionsN(nds);


	if( number_of_nodes < rModelPart.Nodes().size()+1 ){
	  number_of_nodes = rModelPart.Nodes().size()+1;
	  std::cout<<" WARNING: ISOLATED NODES DELETED "<<std::endl;
	}


	std::vector<VariablesListDataValueContainer> VariablesListVector(number_of_nodes);
	std::vector<int> UniquePosition (number_of_nodes);
	std::fill( UniquePosition.begin(), UniquePosition.end(), 0 );

	//unsigned int  step_data_size = rModelPart.GetNodalSolutionStepDataSize();
	VariablesList&  rVariablesList = rModelPart.GetNodalSolutionStepVariablesList();

	//find the center and "radius" of the element
	double  radius = 0;
	Node<3> center(0,0.0,0.0,0.0);

	unsigned int MaximumNumberOfPointsInRadius = list_of_nodes.size();
	std::vector<Node<3>::Pointer> PointsInRadius (MaximumNumberOfPointsInRadius);
	std::vector<double>  PointsInRadiusDistances (MaximumNumberOfPointsInRadius);

	//geometry
	std::vector<std::vector<double> > ElementPointCoordinates(nds);
	std::vector<double> PointCoordinates(3); //dimension=3
	std::fill( PointCoordinates.begin(), PointCoordinates.end(), 0.0 );
	std::fill( ElementPointCoordinates.begin(), ElementPointCoordinates.end(), PointCoordinates );


	for(ModelPart::ElementsContainerType::const_iterator ie = rModelPart.ElementsBegin();
	    ie != rModelPart.ElementsEnd(); ie++)
	  {

	    for(unsigned int i=0; i<nds; ++i)
	      {
		//ID[cn] = rElementsList[el][cn].Id();
		const array_1d<double,3>& Displacement = ie->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
		PointCoordinates[0] = ie->GetGeometry()[i].X0() + Displacement[0];
		PointCoordinates[1] = ie->GetGeometry()[i].Y0() + Displacement[1];
		PointCoordinates[2] = ie->GetGeometry()[i].Z0() + Displacement[2];

		ElementPointCoordinates[i] = PointCoordinates;
	      }

	    std::fill( PointCoordinates.begin(), PointCoordinates.end(), 0.0 );
	    MeshDataTransferUtilities DataTransferUtilities;
	    DataTransferUtilities.CalculateCenterAndSearchRadius( ElementPointCoordinates, PointCoordinates, radius );

	    //find all of the new nodes within the radius
	    center.X() = PointCoordinates[0];
	    center.Y() = PointCoordinates[1];
	    center.Z() = PointCoordinates[2];

	    double Radius = radius * 1.15;
	    int NumberOfPointsInRadius = NodesTree.SearchInRadius (center, Radius, PointsInRadius.begin(), PointsInRadiusDistances.begin(),  MaximumNumberOfPointsInRadius);

	    //std::cout<<"[ID:"<<ie->Id()<<"]: NumberOfPointsInRadius "<<NumberOfPointsInRadius<<" Radius "<<Radius<<std::endl;

	    //check if inside and eventually interpolate
	    for(std::vector<Node<3>::Pointer>::iterator it_found = PointsInRadius.begin(); it_found != (PointsInRadius.begin() + NumberOfPointsInRadius) ; ++it_found)
	      {

		if((*it_found)->IsNot(TO_ERASE)){

		  //std::cout<<" Found ID "<<(*it_found)->Id()<<std::endl;

		  PointCoordinates[0] = (*it_found)->X();
		  PointCoordinates[1] = (*it_found)->Y();
		  PointCoordinates[2] = (*it_found)->Z();

		  bool is_inside = MesherUtilities::CalculatePosition( ElementPointCoordinates, PointCoordinates, ShapeFunctionsN );

		  if(is_inside == true)
		    {
		      //std::cout<<"  Node interpolation: "<<(*it_found)->Id()<<" VariablesList size "<<VariablesListVector.size()<<std::endl;
		      //std::cout<<"  Node interpolation: "<<(*it_found)->Id()<<" N ["<<ie->GetGeometry()[0].Id()<<"] "<<ShapeFunctionsN[0]<<" ["<<ie->GetGeometry()[1].Id()<<"] "<<ShapeFunctionsN[1]<<" ["<<ie->GetGeometry()[2].Id()<<"] "<<ShapeFunctionsN[2]<<std::endl;
		      if(UniquePosition [(*it_found)->Id()] == 0){

			UniquePosition [(*it_found)->Id()] = 1;

			double alpha = 0.25; //[0,1] //smoothing level of the interpolation

			MeshDataTransferUtilities DataTransferUtilities;
			VariablesListVector[(*it_found)->Id()] = DataTransferUtilities.InterpolateVariables( ie->GetGeometry(), ShapeFunctionsN, rVariablesList, (*it_found), alpha );

		      }
		      else{
			UniquePosition [(*it_found)->Id()] += 1;
			//std::cout<<" Node "<<(*it_found)->Id()<<" is relocated again in a element:: "<<ie->Id()<<" num locations "<<UniquePosition [(*it_found)->Id()]<<std::endl;
			//std::cout<<" ShapeFunctionsN "<<ShapeFunctionsN.size()<<std::endl;
		      }

		    }
		}
	      }

	  }



	//the search moves the nodes order using its PreIds
	rModelPart.Nodes().Sort();


	//*******************************************************************
	//CREATE NEW NODE INFORMATION:

	int id=0;
	for(NodesContainerType::iterator i_node = rModelPart.NodesBegin() ; i_node != rModelPart.NodesEnd() ; ++i_node)
	  {

	    if( UniquePosition[i_node->Id()] && i_node->IsNot(TO_ERASE) ){

	      if ( i_node->SolutionStepsDataHas(DISPLACEMENT) == false)
		{
		  std::cout << " WIERD " << std::endl;
		  std::cout << " Laplacian. ThisNode Does not have displacemenet " << i_node->Id() << std::endl;
		  std::cout << "    X: " << i_node->X() << " Y: " << i_node->Y() <<  " Z: " << i_node->Z() << std::endl;
		}

	      if( i_node->IsNot(BOUNDARY) ){

		//recover the original position of the node
		id = i_node->Id();

		i_node->SolutionStepData() = VariablesListVector[id];

		const array_1d<double,3>& disp = i_node->FastGetSolutionStepValue(DISPLACEMENT);

		i_node->X0() = i_node->X() - disp[0];
		i_node->Y0() = i_node->Y() - disp[1];
		i_node->Z0() = i_node->Z() - disp[2];

	      }
	      else if ( i_node->Is(BOUNDARY) && i_node->Is(INSIDE) ){ // Set the position of boundary laplacian

		i_node->Set(INSIDE, false);

		//recover the original position of the node
		id = i_node->Id();

		i_node->SolutionStepData() = VariablesListVector[id];

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
	    else{


	      if ( i_node->SolutionStepsDataHas(DISPLACEMENT) == false)
		{
		  std::cout << " WIERD " << std::endl;
		  std::cout << " Laplacian. ThisNode Does not have displacemenet " << i_node->Id() << std::endl;
		  std::cout << "    X: " << i_node->X() << " Y: " << i_node->Y() <<  " Z: " << i_node->Z() << std::endl;
		}

	      //recover the original position of the node
	      id = i_node->Id();

	      const array_1d<double,3>& disp = i_node->FastGetSolutionStepValue(DISPLACEMENT);

	      i_node->X0() = i_node->X() - disp[0];
	      i_node->Y0() = i_node->Y() - disp[1];
	      i_node->Z0() = i_node->Z() - disp[2];

	      if( GetEchoLevel() > 1 ){
		std::cout << " OUT::PSEUDOERROR: IN this line, there is a node that does not have displacement" << std::endl;
		std::cout << " Laplacian. ThisNode new information does not have displacement " << i_node->Id() << std::endl;
		std::cout << " THIS IS BECAUSE THE NODE is out of the DOMAIN and the interpolation is wrong" << std::endl;
		std::cout << "    X: " << i_node->X() << " Y: " << i_node->Y() <<  " Z: " << i_node->Z() << std::endl;
	      }

	    }

	    //std::cout<<"(B) ID: "<<i_node->Id()<<" NODAL_H "<<i_node->FastGetSolutionStepValue(NODAL_H)<<std::endl;
	  }


      }
      else{

	for(NodesContainerType::iterator i_node = rModelPart.NodesBegin() ; i_node != rModelPart.NodesEnd() ; ++i_node)
	  {

	    //recover the original position of the node
	    array_1d<double,3>& disp = i_node->FastGetSolutionStepValue(DISPLACEMENT);
	    if(norm_2(disp)>0){
	      i_node->X0() = i_node->X() - disp[0];
	      i_node->Y0() = i_node->Y() - disp[1];
	      i_node->Z0() = i_node->Z() - disp[2];
	    }

	    //Set the position of boundary laplacian (Reset the flag)
	    if ( i_node->Is(BOUNDARY) && i_node->IsNot(TO_ERASE) && i_node->Is(INSIDE) )
	      {
		i_node->Set(INSIDE, false); //LMV.
	      }

	  }

      }


      KRATOS_CATCH( "" )
    }



    //*******************************************************************************************
    //*******************************************************************************************

    void SetBoundarySmoothing(ModelPart& rModelPart,
			      std::vector<int> & PreservedElements,
			      const int* pElementsList,
			      const int& NumberOfPoints)
    {
      KRATOS_TRY

      NodesContainerType& rNodes = rModelPart.Nodes();

      ModelPart::ElementsContainerType::iterator element_begin = rModelPart.ElementsBegin();

      const unsigned int nds = element_begin->GetGeometry().size();


      //*******************************************************************
      //NEIGHBOUR NODES:

      std::vector<int> EmptyVector(0);
      std::vector<std::vector<int> >  NeighborNodesList(rNodes.size());
      std::fill( NeighborNodesList.begin(), NeighborNodesList.end(), EmptyVector );

      this->GetBoundaryNeighborNodes(rNodes, NeighborNodesList, PreservedElements, pElementsList, NumberOfPoints, nds);

      //*******************************************************************
      //MOVE BOUNDARY NODES: LAPLACIAN SMOOTHING:

      double convergence_tol =0.001;
      double smoothing_factor=0.1; //0.1
      double smoothing_iters =4; //3,4
      double iters=0;

      bool   simple = true; //weight = 1;
      bool   converged = false;

      double MaxLength=0;
      double NewMaxLength=0;

      unsigned int number_of_nodes = 0;

      ModelPart::NodesContainerType::iterator nodes_begin = rModelPart.NodesBegin();

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

	for(unsigned int in = 0; in<rNodes.size(); ++in)
	  {
	    unsigned int NumberOfNeighbours = NeighborNodesList[in+1].size();

	    if(rNodes[in+1].Is(BOUNDARY) && rNodes[in+1].IsNot(TO_ERASE) && rNodes[in+1].IsNot(BOUNDARY) &&
	       rNodes[in+1].Is(INSIDE) && NumberOfNeighbours>1 )
	      {
		TotalDistance.clear();
		TotalWeight = 0;
		Weight = 0;

		//point position
		P[0] = (nodes_begin+in)->X();
		P[1] = (nodes_begin+in)->Y();
		P[2] = (nodes_begin+in)->Z();


		//std::cout<<" Initial Position: "<<P<<std::endl;
		Length = 0;

		for(unsigned int i = 0; i < NumberOfNeighbours; ++i)
		  {
		    //neighbour position
		    Q[0] = (nodes_begin+(NeighborNodesList[in+1][i]-1))->X();
		    Q[1] = (nodes_begin+(NeighborNodesList[in+1][i]-1))->Y();
		    Q[2] = (nodes_begin+(NeighborNodesList[in+1][i]-1))->Z();


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


		(nodes_begin+in)->X() = P[0];
		(nodes_begin+in)->Y() = P[1];
		(nodes_begin+in)->Z() = P[2];


		number_of_nodes +=1;

	      }

	    //rNodes[in+1].Set(INSIDE,false); //LMV: Reset the flag after interpolation. Indeed, if the flag is set, only one iteration takes place
	  }


	if( (NewMaxLength-MaxLength)/NewMaxLength < convergence_tol ){
	  converged = true;
	  if( GetEchoLevel() > 0 )
	    std::cout<<"   Laplacian smoothing convergence achieved "<<std::endl;
	}


	iters++;

      }

      if(iters==smoothing_iters && !converged){
	if( GetEchoLevel() > 0 )
	  std::cout<<"     WARNING: Boundary Laplacian smoothing convergence NOT achieved (iters:"<<iters<<")"<<std::endl;
      }

      KRATOS_CATCH( "" )
    }

    //*******************************************************************************************
    //*******************************************************************************************

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

    //*******************************************************************************************
    //*******************************************************************************************

    void GetNeighborNodes (std::vector<std::vector<int> >& list_of_neighbor_nodes, std::vector<int> & PreservedElements,const int* ElementList, const int& NumberOfPoints, const unsigned int& nds)
    {

      KRATOS_TRY

      if( (int)list_of_neighbor_nodes.size() != NumberOfPoints+1 ){
	list_of_neighbor_nodes.resize(NumberOfPoints+1);
      }
      std::vector<int> empty_vector(0);
      std::fill( list_of_neighbor_nodes.begin(), list_of_neighbor_nodes.end(), empty_vector );

      bool neighb_set  = false;
      int  neighb_size = 0;

      for(unsigned int el = 0; el<PreservedElements.size(); ++el)
	{
	  if(PreservedElements[el])
	    {
	      //a) Create list of node neighbors (list_of_neighbor_nodes)
	      for(unsigned int ipn=0; ipn<nds; ++ipn)
		{

		  for(unsigned int jpn=0; jpn<nds; ++jpn)
		    {
		      if(ipn!=jpn){
			//add unique node neighbor
			neighb_size = list_of_neighbor_nodes[ElementList[el*nds+ipn]].size();
			neighb_set = false;
			for(int npn=0; npn<neighb_size; ++npn)
			  {
			    if( list_of_neighbor_nodes[ElementList[el*nds+ipn]][npn]==(ElementList[el*nds+jpn]) ){
			      neighb_set=true;
			    }
			  }
			if(neighb_set==false){
			  list_of_neighbor_nodes[ElementList[el*nds+ipn]].push_back(ElementList[el*nds+jpn]);
			}
		      }
		    }
		}
	    }
	}


      KRATOS_CATCH( "" )

    }

    //*******************************************************************************************
    //*******************************************************************************************

    void GetBoundaryNeighborNodes (NodesContainerType& rNodes, std::vector<std::vector<int> >& list_of_neighbor_nodes, std::vector<int> & PreservedElements,const int* ElementList, const int& NumberOfPoints, const unsigned int& nds)
    {

      KRATOS_TRY

      if( (int)list_of_neighbor_nodes.size() != NumberOfPoints+1 ){
	list_of_neighbor_nodes.resize(NumberOfPoints+1);
	std::vector<int> empty_vector(0);
	std::fill( list_of_neighbor_nodes.begin(), list_of_neighbor_nodes.end(), empty_vector );
      }

      bool neighb_set  = false;
      int  neighb_size = 0;

      for(unsigned int el = 0; el<PreservedElements.size(); ++el)
	{
	  if(PreservedElements[el])
	    {
	      //a) Create list of node neighbors (list_of_neighbor_nodes)
	      for(unsigned int ipn=0; ipn<nds; ++ipn)
		{
		  if(rNodes[ElementList[el*nds+ipn]].Is(BOUNDARY)){
		      for(unsigned int jpn=0; jpn<nds; ++jpn)
			{
			  if(ipn!=jpn && rNodes[ElementList[el*nds+jpn]].Is(BOUNDARY)){
			    //add unique node neighbor
			    neighb_size = list_of_neighbor_nodes[ElementList[el*nds+ipn]].size();
			    neighb_set = false;
			    for(int npn=0; npn<neighb_size; ++npn)
			      {
				if( list_of_neighbor_nodes[ElementList[el*nds+ipn]][npn]==(ElementList[el*nds+jpn]) ){
				  neighb_set=true;
				}
			      }
			    if(neighb_set==false){
			      list_of_neighbor_nodes[ElementList[el*nds+ipn]].push_back(ElementList[el*nds+jpn]);
			    }
			  }
			}
		    }
		}
	    }
	}


      KRATOS_CATCH( "" )

    }


    //*******************************************************************************************
    //*******************************************************************************************

    std::vector<double>  SetRanks (ModelPart& rModelPart,
				   struct triangulateio &out,
				   std::vector<std::vector<int> > & list_of_neighbor_nodes)
    {

      KRATOS_TRY

      //set ranks

      std::vector<double> nodes_ranks;
      nodes_ranks.resize(MesherUtilities::GetMaxNodeId(rModelPart)+1);
      //nodes_ranks.resize(rModelPart.NumberOfNodes()+1); //mesh 0
      std::fill( nodes_ranks.begin(), nodes_ranks.end(), 0 );


      //initial assignation
      for(int in = 0; in<out.numberofpoints; ++in)
	{
	  bool contact_active = false;

	  if( rModelPart.Nodes()[in+1].SolutionStepsDataHas(CONTACT_FORCE) ){
	    array_1d<double, 3 > & ContactForceNormal  = rModelPart.Nodes()[in+1].FastGetSolutionStepValue(CONTACT_FORCE);
	    if(norm_2(ContactForceNormal) !=0 )
	      contact_active = true;
	  }

	  if( contact_active && rModelPart.Nodes()[in+1].IsNot(BOUNDARY)){
	    nodes_ranks[in+1]=5;
	  }
	  // else if( norm_2(ContactForceNormal)==0 && rModelPart.Nodes()[in+1].Is(BOUNDARY)){
	  //   nodes_ranks[in+1]=1;
	  // }

	}


      //RANK assignation:
      double rang_assign = 0;
      double rang_top = 5; //3;

      while (rang_assign<rang_top){

	for(int in = 0; in<out.numberofpoints; ++in)
	  {
	    if(nodes_ranks[in+1]==rang_assign){

	      //Rank 0
	      unsigned int shared_node=1;
	      unsigned int NumberOfNeighbours = list_of_neighbor_nodes[in+1].size();
	      for(unsigned int i = 0; i < NumberOfNeighbours; ++i)
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

      KRATOS_CATCH( "" )

    }

    //*******************************************************************************************
    //*******************************************************************************************

    std::vector<double>  SetFirstLayer (ModelPart& rModelPart,
					struct triangulateio &out,
					std::vector<std::vector<int> > & list_of_neighbor_nodes)
    {
      KRATOS_TRY

      //set ranks

      std::vector<double> nodes_layer;
      nodes_layer.resize(MesherUtilities::GetMaxNodeId(rModelPart)+1);
      //nodes_layer.resize(rModelPart.NumberOfNodes()+1); //mesh 0
      std::fill( nodes_layer.begin(), nodes_layer.end(), 0 );


      //initial assignation
      for(int in = 0; in<out.numberofpoints; ++in)
	{
	  if(rModelPart.Nodes()[in+1].IsNot(BOUNDARY)){
	    nodes_layer[in+1]=2;
	  }

	}


      //LAYER assignation:
      double layer_assign = 0;
      double layer_top    = 1;

      while (layer_assign<layer_top){

	for(int in = 0; in<out.numberofpoints; ++in)
	  {
	    if(nodes_layer[in+1]==layer_assign){

	      //Rank 0
	      unsigned int shared_node=1;
	      unsigned int NumberOfNeighbours = list_of_neighbor_nodes[in+1].size();
	      for(unsigned int i = 0; i < NumberOfNeighbours; ++i)
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

      KRATOS_CATCH( "" )

    }


    //*******************************************************************************************
    //*******************************************************************************************

    //GENERAL :: TODO
    //Note : to increase de robustness I propose to detect the layer nodes, via setting ranks to nodes
    // then ensure that the layer nodes have a certain distance to the boundary. Bigger than the offset applied


    void SetInsideProjection (ModelPart& rModelPart,
			      struct triangulateio &out,
			      std::vector<std::vector<int> > & list_of_neighbor_nodes)
    {

      KRATOS_TRY

      std::vector<double> nodes_ranks = SetRanks(rModelPart,out,list_of_neighbor_nodes);

      std::vector<double> nodes_layer = SetFirstLayer(rModelPart,out,list_of_neighbor_nodes);

      double movement_factor = 1.2;
      double contact_factor  = 2.0;
      const array_1d<double,3> ZeroPoint(3,0.0);

      std::vector<array_1d<double,3> > initial_nodes_distances (MesherUtilities::GetMaxNodeId(rModelPart)+1);
      //std::vector<array_1d<double,3> > initial_nodes_distances (rModelPart.NumberOfNodes()+1);
      std::fill( initial_nodes_distances.begin(), initial_nodes_distances.end(), ZeroPoint );

      for(int in = 0; in<out.numberofpoints; ++in)
	{

	  //if(nodes_ranks[in+1]<=1)
	  if(nodes_ranks[in+1]<1)
	    {
	      array_1d<double, 3 >  DeltaDisplacement  = rModelPart.Nodes()[in+1].FastGetSolutionStepValue(DISPLACEMENT) - rModelPart.Nodes()[in+1].FastGetSolutionStepValue(DISPLACEMENT,1);

	      array_1d<double, 3>&  Normal= rModelPart.Nodes()[in+1].FastGetSolutionStepValue(NORMAL); //BOUNDARY_NORMAL must be set as nodal variable

	      double projection=inner_prod(DeltaDisplacement,Normal);

	      bool contact_active = false;
	      if( rModelPart.Nodes()[in+1].SolutionStepsDataHas(CONTACT_FORCE) ){
		array_1d<double, 3 > & ContactForce = rModelPart.Nodes()[in+1].FastGetSolutionStepValue(CONTACT_FORCE);
		if(norm_2(ContactForce)!=0)
		  contact_active = true;
	      }

	      if(contact_active){
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
	      for(unsigned int i = 0; i < NumberOfNeighbours; ++i)
		{
		  shared_node = list_of_neighbor_nodes[in+1][i];

		  if(rModelPart.Nodes()[shared_node].Is(BOUNDARY)){

		    //neighbour position
		    Q[0] = out.pointlist[(list_of_neighbor_nodes[in+1][i]-1)*dimension];
		    Q[1] = out.pointlist[(list_of_neighbor_nodes[in+1][i]-1)*dimension+1];
		    Q[2] = 0;
		    array_1d<double, 3>&  Normal= rModelPart.Nodes()[shared_node].FastGetSolutionStepValue(NORMAL); //BOUNDARY_NORMAL must be set as nodal variable

		    D = Q-P;

		    double projection=inner_prod(D,Normal);

		    array_1d<double, 3>& Offset = rModelPart.Nodes()[shared_node].FastGetSolutionStepValue(OFFSET);
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
	for(int in = 0; in<out.numberofpoints; ++in)
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


		for(unsigned int i = 0; i < NumberOfNeighbours; ++i)
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
      for(int in = 0; in<out.numberofpoints; ++in)
	{

	  // array_1d<double, 3>&  Projection= rModelPart.Nodes()[in+1].FastGetSolutionStepValue(FORCE_EXTERNAL); //BOUNDARY_NORMAL must be set as nodal variable
	  // Projection = initial_nodes_distances[in+1];

	  if(nodes_ranks[in+1]>0)
	    {
	      //std::cout<<" Projection Set "<<initial_nodes_distances[in+1]<<std::endl;
	      out.pointlist[in*dimension]   += initial_nodes_distances[in+1][0];
	      out.pointlist[in*dimension+1] += initial_nodes_distances[in+1][1];
	    }


	}

      KRATOS_CATCH( "" )


    }


    //*******************************************************************************************
    //*******************************************************************************************


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


    //*******************************************************************************************
    //*******************************************************************************************


    inline void Clear(ModelPart::NodesContainerType::iterator node_it,  int step_data_size )
    {
      unsigned int buffer_size = node_it->GetBufferSize();

      for(unsigned int step = 0; step<buffer_size; ++step)
	{
	  //getting the data of the solution step
	  double* step_data = (node_it)->SolutionStepData().Data(step);

	  //copying this data in the position of the vector we are interested in
	  for(int j= 0; j< step_data_size; ++j)
	    {
	      step_data[j] = 0.0;
	    }
	}

    }

    //*******************************************************************************************
    //*******************************************************************************************


    inline void ClearVariables(ModelPart::NodesContainerType::iterator node_it , Variable<array_1d<double,3> >& rVariable)
    {
      array_1d<double, 3>& Aux_var = node_it->FastGetSolutionStepValue(rVariable, 0);

      noalias(Aux_var) = ZeroVector(3);

    }

    //*******************************************************************************************
    //*******************************************************************************************


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
