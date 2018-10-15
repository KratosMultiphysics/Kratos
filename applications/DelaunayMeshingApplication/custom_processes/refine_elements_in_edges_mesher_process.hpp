//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined(KRATOS_REFINE_ELEMENTS_IN_EDGES_MESHER_PROCESS_H_INCLUDED )
#define  KRATOS_REFINE_ELEMENTS_IN_EDGES_MESHER_PROCESS_H_INCLUDED


// External includes

// System includes

// Project includes
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"

#include "includes/model_part.h"
#include "custom_utilities/mesh_error_calculation_utilities.hpp"
#include "custom_utilities/mesher_utilities.hpp"
#include "custom_processes/mesher_process.hpp"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Refine Mesh Elements Process 2D and 3D
/** The process inserts nodes in the edge elements to be splitted in the remeshing process
*/

class RefineElementsInEdgesMesherProcess
    : public MesherProcess
{
 public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of Process
  KRATOS_CLASS_POINTER_DEFINITION( RefineElementsInEdgesMesherProcess );

  typedef ModelPart::ConditionType         ConditionType;
  typedef ModelPart::PropertiesType       PropertiesType;
  typedef ConditionType::GeometryType       GeometryType;

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  RefineElementsInEdgesMesherProcess(ModelPart& rModelPart,
				     MesherUtilities::MeshingParameters& rRemeshingParameters,
				     int EchoLevel)
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
  {
    mEchoLevel = EchoLevel;
  }


  /// Destructor.
  virtual ~RefineElementsInEdgesMesherProcess() {}


  ///@}
  ///@name Operators
  ///@{

  /// This operator is provided to call the process as a function and simply calls the Execute method.
  void operator()()
  {
    Execute();
  }


  ///@}
  ///@name Operations
  ///@{


  /// Execute method is used to execute the Process algorithms.
  void Execute() override
  {
    KRATOS_TRY

    if( ( mrRemesh.Refine->RefiningOptions.Is(MesherUtilities::REFINE_ADD_NODES) ||
        mrRemesh.Refine->RefiningOptions.Is(MesherUtilities::REFINE_INSERT_NODES) ) )
    {

      //1.- Select Elements to split (edge elements with all nodes as boundary-free surface)
      ModelPart::ElementsContainerType   BoundaryEdgedElements;
      ModelPart::ConditionsContainerType BoundaryEdgedConditions;
      this->SelectFullBoundaryEdgedElements(mrModelPart, BoundaryEdgedElements);

      //2.- Select Inside Faces to refine
      std::vector<Geometry< Node<3> > > ListOfFacesToSplit;
      this->SelectFacesToSplit(BoundaryEdgedElements,ListOfFacesToSplit);


      //3.- Create and insert new nodes
      std::vector<Node<3>::Pointer>  ListOfNewNodes;
      this->GenerateNewNodes(mrModelPart,ListOfNewNodes,ListOfFacesToSplit);

      //4.- Insert new nodes to model part
      this->SetNodesToModelPart(mrModelPart, ListOfNewNodes);
    }

    KRATOS_CATCH(" ")
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
  std::string Info() const override
  {
    return "RefineElementsInEdgesMesherProcess";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << "RefineElementsInEdgesMesherProcess";
  }

  /// Print object's data.
  void PrintData(std::ostream& rOStream) const override
  {
  }


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

  ModelPart& mrModelPart;

  MesherUtilities::MeshingParameters& mrRemesh;

  MesherUtilities mMesherUtilities;

  int mEchoLevel;

  ///@}
  ///@name Protected Operators
  ///@{
  ///@}
  ///@name Protected Operations
  ///@{
  ///@}


  //**************************************************************************
  //**************************************************************************

  virtual void SelectFullBoundaryEdgedElements(ModelPart& rModelPart,
                                               ModelPart::ElementsContainerType& rBoundaryEdgedElements)
  {
    KRATOS_TRY

    bool is_full_boundary = false;
    for(ModelPart::ElementsContainerType::iterator i_elem = rModelPart.ElementsBegin();
        i_elem != rModelPart.ElementsEnd(); ++i_elem)
    {
      Geometry< Node<3> >& rGeometry = i_elem->GetGeometry();

      is_full_boundary = true;
      for(unsigned int i=0; i<rGeometry.size(); ++i)
      {
        if( rGeometry[i].IsNot(BOUNDARY) ){
          is_full_boundary = false;
          break;
        }
      }

      if( is_full_boundary ){
        rBoundaryEdgedElements.push_back(*(i_elem.base()));
      }

    }

    KRATOS_CATCH( "" )
  }


  //**************************************************************************
  //**************************************************************************

  void SelectFacesToSplit(ModelPart::ElementsContainerType& rBoundaryEdgedElements,
                          std::vector<Geometry< Node<3> > >& rListOfFacesToSplit)

  {
    KRATOS_TRY

    DenseMatrix<unsigned int> lpofa; //connectivities of points defining faces
    DenseVector<unsigned int> lnofa; //number of points defining faces

    for(ModelPart::ElementsContainerType::iterator i_elem = rBoundaryEdgedElements.begin();
        i_elem != rBoundaryEdgedElements.end(); ++i_elem)
    {
      WeakPointerVector<Element>& neighb_elems = i_elem->GetValue(NEIGHBOUR_ELEMENTS);

      unsigned int face=0;
      bool accepted_face = false;
      for(WeakPointerVector< Element >::iterator ne = neighb_elems.begin(); ne!=neighb_elems.end(); ++ne)
      {

        if (ne->Id() != i_elem->Id())  // If there is a shared element in face nf
        {
          accepted_face = true;

          Geometry< Node<3> >& rGeometry = i_elem->GetGeometry();

          rGeometry.NodesInFaces(lpofa);
          rGeometry.NumberNodesInFaces(lnofa);

          unsigned int NumberNodesInFace = lnofa[face];
          Condition::NodesArrayType FaceNodes;
          FaceNodes.reserve(NumberNodesInFace);

          unsigned int split_counter=0;
          unsigned int counter=0;
          for(unsigned int j=1; j<=NumberNodesInFace; ++j)
          {
            if(rGeometry[lpofa(j,face)].Is(TO_SPLIT))
              ++split_counter;
            FaceNodes.push_back(rGeometry(lpofa(j,face)));

            //do not SPLIT small edges only big edges, else a lot of volume is added
            if(mrModelPart.Is(FLUID) && counter>0){
              if( 4.5 * this->mrRemesh.Refine->CriticalRadius > norm_2(FaceNodes[counter].Coordinates()-FaceNodes[counter-1].Coordinates()) ){
                accepted_face = false;
              }
            }
            ++counter;
          }

          if(split_counter==NumberNodesInFace)
            accepted_face = false;

          unsigned int free_surface = 0;
          if(accepted_face){

            if(mrModelPart.Is(FLUID) ){
              //set FLUID flag in RIGID nodes to false
              //and to not refine the edge in order to erase blocked edge elements
              for(unsigned int i=0; i<FaceNodes.size(); ++i)
              {
                if(FaceNodes[i].Is(RIGID) && FaceNodes[i].Is(FREE_SURFACE)){
                  ++free_surface;
                }
              }

              if(free_surface != 0){
                accepted_face = false;
              }
            }

          }

          if( accepted_face ){
            Geometry<Node<3> > InsideFace(FaceNodes);
            rListOfFacesToSplit.push_back(InsideFace);

            //set TO_SPLIT to make the insertion unique
            for(unsigned int i=0; i<FaceNodes.size(); ++i)
              FaceNodes[i].Set(TO_SPLIT);

            break;
          }

        }
        ++face;
      }
    }

    // reset TO_SPLIT
    for(std::vector<Geometry< Node<3> > >::iterator nf = rListOfFacesToSplit.begin(); nf != rListOfFacesToSplit.end(); ++nf)
    {
      for(unsigned int i=0; i<nf->size(); ++i)
      {
        (*nf)[i].Set(TO_SPLIT,false);
      }
    }

    // std::cout<<" rEdgeElements "<< rBoundaryEdgedElements.size()<<std::endl;
    // std::cout<<" FacesToSplit "<<rListOfFacesToSplit.size()<<std::endl;

    KRATOS_CATCH( "" )
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void GenerateNewNodes(ModelPart& rModelPart,
                        std::vector<Node<3>::Pointer>& rListOfNewNodes,
                        std::vector<Geometry<Node<3> > >& rListOfFacesToSplit)
  {
    KRATOS_TRY

    //ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

    MeshDataTransferUtilities DataTransferUtilities;

    Node<3>::Pointer pNode;

    //center
    double xc = 0;
    double yc = 0;
    double zc = 0;

    //radius
    double radius = 0;

    //assign data to dofs
    Node<3>::DofsContainerType& ReferenceDofs = rModelPart.Nodes().front().GetDofs();

    VariablesList& VariablesList = rModelPart.GetNodalSolutionStepVariablesList();


    std::vector<double> ShapeFunctionsN;

    unsigned int id = MesherUtilities::GetMaxNodeId(*(rModelPart.GetParentModelPart())) + 1;

    unsigned int size  = 0;
    //unsigned int count = 0;

    for (std::vector<Geometry<Node<3> > >::iterator i_face = rListOfFacesToSplit.begin() ; i_face != rListOfFacesToSplit.end(); ++i_face)
    {

      size = i_face->size();


      ShapeFunctionsN.resize(size);


      if( size == 2 )
        DataTransferUtilities.CalculateCenterAndSearchRadius( (*i_face)[0].X(), (*i_face)[0].Y(),
                                                              (*i_face)[1].X(), (*i_face)[1].Y(),
                                                              xc,yc,radius);


      if( size == 3 )
        DataTransferUtilities.CalculateCenterAndSearchRadius( (*i_face)[0].X(), (*i_face)[0].Y(), (*i_face)[0].Z(),
                                                              (*i_face)[1].X(), (*i_face)[1].Y(), (*i_face)[1].Z(),
                                                              (*i_face)[2].X(), (*i_face)[2].Y(), (*i_face)[2].Z(),
                                                              xc,yc,zc,radius);


      //create a new node
      pNode = Kratos::make_shared< Node<3> >( id, xc, yc, zc );

      //giving model part variables list to the node
      pNode->SetSolutionStepVariablesList(&VariablesList);

      //set buffer size
      pNode->SetBufferSize(rModelPart.GetBufferSize());

      //generating the dofs
      for(Node<3>::DofsContainerType::iterator i_dof = ReferenceDofs.begin(); i_dof != ReferenceDofs.end(); ++i_dof)
      {
        Node<3>::DofType& rDof = *i_dof;
        Node<3>::DofType::Pointer pNewDof = pNode->pAddDof( rDof );

        // in rigid edges set it fix has no sense:
        (pNewDof)->FreeDof();

        // count = 0;
        // for( unsigned int i = 0; i<size; ++i )
        // {
        //   if((*i_face)[i].IsFixed(rDof.GetVariable()))
        //     count++;
        // }

        // if( count == size )
        //   (pNewDof)->FixDof();
        // else
        //   (pNewDof)->FreeDof();
      }

      std::fill(ShapeFunctionsN.begin(), ShapeFunctionsN.end(), 1.0/double(size));

      double alpha = 1;
      DataTransferUtilities.Interpolate( (*i_face), ShapeFunctionsN, VariablesList, pNode, alpha );

      //set flags from one of the nodes
      pNode->AssignFlags((*i_face)[0]);

      if( rModelPart.Is(FLUID) )
        pNode->Set(FLUID);

      if( rModelPart.Is(SOLID) )
        pNode->Set(SOLID);

      pNode->Set(NEW_ENTITY,true);
      pNode->Set(ACTIVE,true);

      //unset boundary flags
      pNode->Set(BOUNDARY,false);
      pNode->Set(FREE_SURFACE,false);
      pNode->Set(RIGID,false);

      //set variables
      this->SetNewNodeVariables(rModelPart, pNode);

      rListOfNewNodes.push_back(pNode);

      id++;
    }

    KRATOS_CATCH( "" )
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void SetNewNodeVariables(ModelPart& rModelPart, Node<3>::Pointer& pNode)
  {
    KRATOS_TRY

    //set model part
    pNode->SetValue(MODEL_PART_NAME,rModelPart.Name());

    //set nodal_h
    pNode->FastGetSolutionStepValue(NODAL_H) = mrRemesh.Refine->CriticalSide*2.0;

    //set original position
    const array_1d<double,3>& Displacement = pNode->FastGetSolutionStepValue(DISPLACEMENT);
    pNode->X0() = pNode->X() - Displacement[0];
    pNode->Y0() = pNode->Y() - Displacement[1];
    pNode->Z0() = pNode->Z() - Displacement[2];

    //reset contact force
    if( pNode->SolutionStepsDataHas(CONTACT_FORCE) )
      pNode->FastGetSolutionStepValue(CONTACT_FORCE).clear();


    KRATOS_CATCH( "" )
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void SetNodesToModelPart(ModelPart& rModelPart,
                           std::vector<Node<3>::Pointer>& rListOfNewNodes)
  {
    KRATOS_TRY

        if(rListOfNewNodes.size()){

          //add new conditions: ( SOLID body model part )
          for(std::vector<Node<3>::Pointer>::iterator i_node = rListOfNewNodes.begin(); i_node!= rListOfNewNodes.end(); ++i_node)
	  {
	    rModelPart.Nodes().push_back(*(i_node));
	  }

        }

    KRATOS_CATCH( "" )
  }

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
  ///@name Private Static Member Variables
  ///@{
  ///@}
  ///@name Private Static Member Variables
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
  RefineElementsInEdgesMesherProcess& operator=(RefineElementsInEdgesMesherProcess const& rOther);


  /// this function is a private function


  /// Copy constructor.
  //Process(Process const& rOther);


  ///@}

}; // Class Process

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  RefineElementsInEdgesMesherProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const RefineElementsInEdgesMesherProcess& rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_REFINE_ELEMENTS_IN_EDGES_MESHER_PROCESS_H_INCLUDED  defined
