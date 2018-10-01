//
//   Project Name:        KratosPfemFluidApplication $
//   Created by:          $Author:           AFranci $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:           June 2017 $
//   Revision:            $Revision:             0.0 $
//
//

#if !defined(KRATOS_SPLIT_ELEMENTS_PROCESS_H_INCLUDED )
#define  KRATOS_SPLIT_ELEMENTS_PROCESS_H_INCLUDED


// System includes

// External includes

// Project includes


//#include "spatial_containers/spatial_containers.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "custom_utilities/mesher_utilities.hpp"

///VARIABLES used:
//Data:
//StepData:
//Flags:    (checked)
//          (set)
//          (modified)
//          (reset)


namespace Kratos
{

  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{
  typedef  ModelPart::NodesContainerType                      NodesContainerType;
  typedef  ModelPart::ElementsContainerType                ElementsContainerType;
  typedef  ModelPart::MeshType::GeometryType::PointsArrayType    PointsArrayType;


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
  /** Detail class definition.
   */
  class SplitElementsProcess
    : public Process
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SplitElementsProcess
    KRATOS_CLASS_POINTER_DEFINITION( SplitElementsProcess );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SplitElementsProcess(ModelPart& rModelPart,
				 int EchoLevel)
      : mrModelPart(rModelPart)
    {
      mEchoLevel = EchoLevel;
    }

    /// Destructor.
    virtual ~SplitElementsProcess()
    {
    }

  void operator()()
    {
      Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
      std::cout<<" Execute() in SplitElementsProcess"<<std::endl;

    };

    ///@}
    ///@name Operators
    ///@{



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
    std::string Info() const override
    {
      return "SplitElementsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "SplitElementsProcess";
    }

    void ExecuteInitialize() override
    {

      KRATOS_TRY

      for(ModelPart::SubModelPartIterator i_mp= mrModelPart.SubModelPartsBegin() ; i_mp!=mrModelPart.SubModelPartsEnd(); ++i_mp)
      	{
      	  if(  (i_mp->Is(ACTIVE) && i_mp->Is(SOLID)) || (i_mp->Is(ACTIVE) && i_mp->Is(FLUID)) ){

	    std::string ComputingModelPartName;
	    ComputingModelPartName = i_mp->Name();
	    ModelPart& rComputingModelPart = mrModelPart.GetSubModelPart(ComputingModelPartName);
	    const unsigned int dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
	    MesherUtilities MesherUtils;
	    double modelPartVolume=MesherUtils.ComputeModelPartVolume(rComputingModelPart);
	    double criticalVolume=0.5*modelPartVolume/double(rComputingModelPart.Elements().size());

      	    for(ElementsContainerType::iterator i_elem = rComputingModelPart.Elements().begin() ; i_elem != rComputingModelPart.Elements().end() ; ++i_elem)
      	      {
      		unsigned int numFreeSurface=0;
      		unsigned int numRigid=0;
      		PointsArrayType& vertices=i_elem->GetGeometry().Points();
      		const unsigned int numNodes = i_elem->GetGeometry().size();
		for(unsigned int i=0; i<numNodes; ++i)
		  {
		    if(vertices[i].Is(FREE_SURFACE) || vertices[i].Is(ISOLATED)) {
		      numFreeSurface++;
		    }else if(vertices[i].Is(RIGID)){
		      numRigid++;
		    }
		  }

		if(dimension==2){
		  double elementalVolume=i_elem->GetGeometry().Area();
		  if(numNodes==3){
		    // if((numRigid+numFreeSurface)==vertices.size() && vertices.size()>2 && elementalVolume>criticalVolume){
		    if(numRigid>0 && numFreeSurface==0 && vertices.size()>2 && elementalVolume>criticalVolume){
		      i_elem->Set(ACTIVE,false);
      		    }else{
		      i_elem->Set(ACTIVE,true);
		    }
		  }else{
		    std::cout<<"split_element_process not yet implemented for quadratic elements"<<std::endl;
		  }
		}
		else if(dimension==3){
		  double elementalVolume = 0;
		  if( i_elem->GetGeometry().Dimension() == 3 )
		    elementalVolume=i_elem->GetGeometry().Volume();
		  if(numNodes==4){
		    if(numRigid==0 && (numRigid+numFreeSurface)==vertices.size() && vertices.size()>3 && elementalVolume>criticalVolume){
		      i_elem->Set(ACTIVE,false);
		    }else{
		      i_elem->Set(ACTIVE,true);
		    }
		  }else{
		    std::cout<<"split_element_process not yet implemented for quadratic elements"<<std::endl;
		  }
		}

      	      }


   	    for(ElementsContainerType::iterator i_elem = rComputingModelPart.Elements().begin() ; i_elem != rComputingModelPart.Elements().end() ; ++i_elem)
      	      {
		if(i_elem->IsNot(ACTIVE)){
		  PointsArrayType& vertices=i_elem->GetGeometry().Points();
		  bool addedNode=false;
		  ModelPart::NodeType::Pointer newNode=this->CreateAndAddNewNodeToSubModelPart(i_mp, i_elem, vertices,addedNode);
		  if(addedNode==true){
		    unsigned int rElementId = 0;
		    ModelPart::PropertiesType::Pointer pProp = i_elem->pGetProperties();
		    if(dimension==2){

		      Triangle2D3<Node < 3 > > firstGeom(newNode,vertices(0),vertices(1));
		      rElementId = MesherUtilities::GetMaxElementId(mrModelPart)+1;
		      Element::Pointer newElement = i_elem->Create(rElementId,firstGeom,pProp);
		      newElement->Set(ACTIVE,true);
		      newElement->Set(FLUID);
		      newElement->Set(TO_ERASE);
		      newElement->Initialize();
		      rComputingModelPart.AddElement(newElement);

		      Triangle2D3<Node < 3 > > secondGeom(newNode,vertices(2),vertices(0));
		      rElementId = MesherUtilities::GetMaxElementId(mrModelPart)+1;
		      newElement = i_elem->Create(rElementId,secondGeom,pProp);
		      newElement->Set(ACTIVE,true);
		      newElement->Set(FLUID);
		      newElement->Set(TO_ERASE);
		      newElement->Initialize();
		      rComputingModelPart.AddElement(newElement);

		      Triangle2D3<Node < 3 > > thirdGeom(newNode,vertices(1),vertices(2));
		      rElementId = MesherUtilities::GetMaxElementId(mrModelPart)+1;
		      newElement = i_elem->Create(rElementId,thirdGeom,pProp);
		      newElement->Set(ACTIVE,true);
		      newElement->Set(FLUID);
		      newElement->Set(TO_ERASE);
		      newElement->Initialize();
		      rComputingModelPart.AddElement(newElement);

		    }else if(dimension==3){

		      Tetrahedra3D4<Node < 3 > > firstGeom(newNode,vertices(1),vertices(2),vertices(3));
		      rElementId = MesherUtilities::GetMaxElementId(mrModelPart)+1;
		      Element::Pointer newElement = i_elem->Create(rElementId,firstGeom,pProp);
		      newElement->Set(ACTIVE,true);
		      newElement->Set(FLUID);
		      newElement->Set(TO_ERASE);
		      newElement->Initialize();
		      rComputingModelPart.AddElement(newElement);

		      Tetrahedra3D4<Node < 3 > > secondGeom(vertices(0),newNode,vertices(2),vertices(3));
		      rElementId = MesherUtilities::GetMaxElementId(mrModelPart)+1;
		      newElement = i_elem->Create(rElementId,secondGeom,pProp);
		      newElement->Set(ACTIVE,true);
		      newElement->Set(FLUID);
		      newElement->Set(TO_ERASE);
		      newElement->Initialize();
		      rComputingModelPart.AddElement(newElement);

		      Tetrahedra3D4<Node < 3 > > thirdGeom(vertices(0),vertices(1),newNode,vertices(3));
		      rElementId = MesherUtilities::GetMaxElementId(mrModelPart)+1;
		      newElement = i_elem->Create(rElementId,thirdGeom,pProp);
		      newElement->Set(ACTIVE,true);
		      newElement->Set(FLUID);
		      newElement->Set(TO_ERASE);
		      newElement->Initialize();
		      rComputingModelPart.AddElement(newElement);

		      Tetrahedra3D4<Node < 3 > > fourthGeom(vertices(0),vertices(1),vertices(2),newNode);
		      rElementId = MesherUtilities::GetMaxElementId(mrModelPart)+1;
		      newElement = i_elem->Create(rElementId,fourthGeom,pProp);
		      newElement->Set(ACTIVE,true);
		      newElement->Set(FLUID);
		      newElement->Set(TO_ERASE);
		      newElement->Initialize();
		      rComputingModelPart.AddElement(newElement);

		    }
		  }else{
		    i_elem->Set(ACTIVE,true);
		  }
		}

	      }

 	    //Sort
	    rComputingModelPart.Nodes().Sort();
	    rComputingModelPart.Elements().Sort();

	    //Unique
	    rComputingModelPart.Nodes().Unique();
	    rComputingModelPart.Elements().Unique();

      	  }
      	}

      // std::cout<<""<<mrModelPart<<std::endl;


      KRATOS_CATCH(" ")
	}



    template< class TDataType > void  AddUniqueWeakPointer
    (WeakPointerVector< TDataType >& v, const typename TDataType::WeakPointer candidate)
    {
      typename WeakPointerVector< TDataType >::iterator i = v.begin();
      typename WeakPointerVector< TDataType >::iterator endit = v.end();
      while ( i != endit && (i)->Id() != (candidate.lock())->Id())
	{
	  ++i;
	}
      if( i == endit )
	{
	  v.push_back(candidate);
	}

    }


    ModelPart::NodeType::Pointer CreateAndAddNewNodeToSubModelPart(ModelPart::SubModelPartIterator i_mp,
								   ElementsContainerType::iterator i_elem,
								   PointsArrayType vertices,
								   bool &addedNode){

      unsigned int numNodes=vertices.size();
      std::vector<std::vector<double> > ElementPointCoordinates(numNodes);
      std::vector<double> PointCoordinates(3);
      array_1d<double, 3> rPoint=ZeroVector(3);
      unsigned int masterNode=0;
      for(unsigned int i=0; i<vertices.size(); ++i)
	{
	  rPoint+=vertices[i].Coordinates()/double(numNodes);
	  PointCoordinates[0] = vertices[i].X();
	  PointCoordinates[1] = vertices[i].Y();
	  PointCoordinates[2] = vertices[i].Z();
	  ElementPointCoordinates[i] = PointCoordinates;
	  if(vertices[i].IsNot(RIGID)){
	    masterNode=i;
	  }
	}
      PointCoordinates[0] =rPoint[0];
      PointCoordinates[1] =rPoint[1];
      PointCoordinates[2] =rPoint[2];
      unsigned int rNodeId = MesherUtilities::GetMaxNodeId(mrModelPart)+1;

      ModelPart::NodeType::Pointer newNode = i_mp->CreateNewNode( rNodeId, rPoint[0], rPoint[1], rPoint[2]);

      //generating the dofs
      ModelPart::NodeType::DofsContainerType& reference_dofs = vertices[masterNode].GetDofs();

      for(ModelPart::NodeType::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); ++iii)
	{
	  ModelPart::NodeType::DofType& rDof = *iii;
	  Node<3>::DofType::Pointer p_new_dof = newNode->pAddDof( rDof );
	  (p_new_dof)->FreeDof();
	}

      MeshDataTransferUtilities DataTransferUtilities;
      std::vector<double> ShapeFunctionsN;
      VariablesList&  variables_list = i_mp->GetNodalSolutionStepVariablesList();

      // //giving model part variables list to the node
      newNode->SetSolutionStepVariablesList(vertices[masterNode].pGetVariablesList());

      // //set buffer size
      newNode->SetBufferSize(vertices[masterNode].GetBufferSize());
      newNode->Set(FLUID);
      newNode->Reset(FREE_SURFACE);
      bool is_inside = MesherUtilities::CalculatePosition(ElementPointCoordinates,PointCoordinates,ShapeFunctionsN );
      if(is_inside == true)
	{
	  double alpha = 1; //1 to interpolate, 0 to leave the original data
	  DataTransferUtilities.Interpolate( i_elem->GetGeometry(), ShapeFunctionsN, variables_list, newNode, alpha );
	  newNode->Set(TO_ERASE);
	  newNode->Set(ACTIVE);
	  newNode->Set(FLUID);
	  const array_1d<double,3>& displacement = newNode->FastGetSolutionStepValue(DISPLACEMENT);
	  newNode->X0() = newNode->X() - displacement[0];
	  newNode->Y0() = newNode->Y() - displacement[1];
	  newNode->Z0() = newNode->Z() - displacement[2];

	  WeakPointerVector<Node<3> >& rN = newNode->GetValue(NEIGHBOUR_NODES);
	  for(unsigned int spn=0; spn<numNodes; ++spn)
	    {
	      ModelPart::NodeType::WeakPointer temp = i_elem->GetGeometry()(spn);
	      this->AddUniqueWeakPointer< Node<3> >(rN,temp);

	    }

	  i_mp->Nodes().push_back(newNode);
	  addedNode=true;
	  // std::cout<<"    NEW NODE "<<rNodeId<<")"<<rPoint[0]<<" "<<rPoint[1]<<" "<<rPoint[2]<<std::endl;

	}else{
	addedNode=false;
	// std::cout<<"NODE IS OUTSIDE "<<rNodeId<<")"<<rPoint[0]<<" "<<rPoint[1]<<" "<<rPoint[2]<<std::endl;
      }
      return newNode;

    }

   void ExecuteFinalize() override
    {

      KRATOS_TRY

	NodesContainerType temporal_nodes;
      temporal_nodes.reserve(mrModelPart.Nodes().size());
      temporal_nodes.swap(mrModelPart.Nodes());
      mrModelPart.Nodes().clear();

      for(NodesContainerType::iterator i_node = temporal_nodes.begin() ; i_node != temporal_nodes.end() ; ++i_node)
      	{
      	  if( i_node->IsNot(TO_ERASE) ){
      	    (mrModelPart.Nodes()).push_back(*(i_node.base()));
      	  }
      	}
      mrModelPart.Nodes().Sort();

      ElementsContainerType temporal_elements;
      temporal_elements.reserve(mrModelPart.Elements().size());
      temporal_elements.swap(mrModelPart.Elements());
      mrModelPart.Elements().clear();
      for(ElementsContainerType::iterator i_elem = temporal_elements.begin() ; i_elem != temporal_elements.end() ; ++i_elem)
      	{
      	  if( i_elem->IsNot(TO_ERASE) ){
      	    (mrModelPart.Elements()).push_back(*(i_elem.base()));
      	  }

      	}

      mrModelPart.Elements().Sort();


      KRATOS_CATCH(" ")
	}





  protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{



    //*******************************************************************************************
    //*******************************************************************************************


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
    ModelPart& mrModelPart;

    int mEchoLevel;


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


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
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    SplitElementsProcess& operator=(SplitElementsProcess const& rOther);

    /// Copy constructor.
    //SplitElementsProcess(SplitElementsProcess const& rOther);


    ///@}

  }; // Class SplitElementsProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    SplitElementsProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const SplitElementsProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif // KRATOS_SPLIT_ELEMENTS_PROCESS_H_INCLUDED  defined
