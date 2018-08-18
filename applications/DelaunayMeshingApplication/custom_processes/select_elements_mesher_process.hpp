//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined( KRATOS_SELECT_ELEMENTS_MESHER_PROCESS_H_INCLUDED )
#define KRATOS_SELECT_ELEMENTS_MESHER_PROCESS_H_INCLUDED


// External includes

// System includes

// Project includes
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"

#include "custom_utilities/mesher_utilities.hpp"
#include "custom_processes/mesher_process.hpp"

///VARIABLES used:
//Data:
//StepData: NODAL_H, CONTACT_FORCE ,VELOCITY, DISPLACEMENT
//Flags:    (checked) TO_ERASE, BOUNDARY, RIGID, SOLID, FREE_SURFACE, FLUID
//          (set)
//          (modified)
//          (reset)
//(set):=(set in this process)

namespace Kratos
{

///@name Kratos Classes
///@{

/// Refine Mesh Elements Process 2D and 3D
/** The process labels the elements to be refined in the mesher
    it applies a size constraint to elements that must be refined.

*/
class SelectElementsMesherProcess
    : public MesherProcess
{
 public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of Process
  KRATOS_CLASS_POINTER_DEFINITION( SelectElementsMesherProcess );

  typedef ModelPart::ConditionType         ConditionType;
  typedef ModelPart::PropertiesType       PropertiesType;
  typedef ConditionType::GeometryType       GeometryType;

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  SelectElementsMesherProcess(ModelPart& rModelPart,
			      MesherUtilities::MeshingParameters& rRemeshingParameters,
			      int EchoLevel)
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
  {
    mEchoLevel = EchoLevel;
  }


  /// Destructor.
  virtual ~SelectElementsMesherProcess() {}


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

    if( mEchoLevel > 0 )
      std::cout<<" [ SELECT MESH ELEMENTS: ("<<mrRemesh.OutMesh.GetNumberOfElements()<<") "<<std::endl;

    const int& OutNumberOfElements = mrRemesh.OutMesh.GetNumberOfElements();
    mrRemesh.PreservedElements.clear();
    mrRemesh.PreservedElements.resize(OutNumberOfElements);
    std::fill( mrRemesh.PreservedElements.begin(), mrRemesh.PreservedElements.end(), 0 );
    mrRemesh.MeshElementsSelectedFlag = true;

    mrRemesh.Info->NumberOfElements=0;

    bool box_side_element = false;
    bool wrong_added_node = false;

    unsigned int number_of_slivers = 0;

    unsigned int passed_alpha_shape = 0;
    unsigned int passed_inner_outer = 0;


    if(mrRemesh.ExecutionOptions.IsNot(MesherUtilities::SELECT_TESSELLATION_ELEMENTS))
    {
      for(int el=0; el<OutNumberOfElements; ++el)
      {
        mrRemesh.PreservedElements[el]=1;
        mrRemesh.Info->NumberOfElements+=1;
      }
    }
    else
    {
      if( mEchoLevel > 0 )
        std::cout<<"   Start Element Selection "<<OutNumberOfElements<<std::endl;

      //reassign fluid and rigid flags to eliminate full rigid elements
      //if( !IsFirstTimeStep() )
      this->LabelEdgeNodes(mrModelPart);

      unsigned int dimension = 0;
      unsigned int number_of_vertices = 0;

      this->GetElementDimension(dimension,number_of_vertices);

      const int* OutElementList = mrRemesh.OutMesh.GetElementList();

      ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes();

      int el = 0;
      int number = 0;

      //void CheckIds(OuterElementsList,OutNumberOfElements,number_of_vertices);

      //#pragma omp parallel for reduction(+:number) private(el)
      for(el=0; el<OutNumberOfElements; ++el)
      {
        GeometryType Vertices;

        // std::cout<<" selected vertices ["<<OutElementList[el*number_of_vertices];
        // for(unsigned int d=1; d<number_of_vertices; ++d)
        // 	{
        // 	  std::cout<<", "<<OutElementList[el*number_of_vertices+d];
        // 	}
        // std::cout<<"] "<<std::endl;

        wrong_added_node = false;
        box_side_element = false;

        NodalFlags VerticesFlags;
        for(unsigned int pn=0; pn<number_of_vertices; ++pn)
        {
          //std::cout<<" pn "<<pn<<" id "<<OutElementList[id]<<" size "<<rNodes.size()<<" IDS "<<mrRemesh.NodalPreIds.size()<<" preid "<<mrRemesh.NodalPreIds[OutElementList[id]]<<std::endl;
          unsigned int id = el*number_of_vertices+pn;

          if(OutElementList[id]<=0)
            std::cout<<" ERROR: something is wrong: nodal id < 0 "<<el<<std::endl;

          //check if the number of nodes are considered in the nodal pre ids
          if( (unsigned int)OutElementList[id] >= mrRemesh.NodalPreIds.size() ){
            if(mrRemesh.Options.Is(MesherUtilities::CONTACT_SEARCH))
              wrong_added_node = true;
            std::cout<<" ERROR: something is wrong: node out of bounds "<<std::endl;
            break;
          }

          //check if is a vertex of an artificial external bounding box
          if(mrRemesh.NodalPreIds[OutElementList[id]]<0){
            if(mrRemesh.Options.IsNot(MesherUtilities::CONTACT_SEARCH))
              std::cout<<" ERROR: something is wrong: nodal id < 0 "<<std::endl;
            box_side_element = true;
            break;
          }

          //get node from model part and set it as vertex
          Vertices.push_back(rNodes(OutElementList[id]));

          //vertices flags
          VerticesFlags.CountFlags(Vertices.back());
        }


        if(box_side_element || wrong_added_node){
          //std::cout<<" Box_Side_Element "<<std::endl;
          continue;
        }

        bool accepted=false;

        //monitoring bools
        // bool boundary_accepted = false;
        // bool alpha_accepted = false;
        // bool subdomain_accepted = false;
        // bool center_accepted = false;
        // bool shape_accepted = false;

        accepted = this->CheckElementBoundaries(Vertices,VerticesFlags);

        double Alpha = mrRemesh.AlphaParameter;

        this->GetAlphaParameter(Alpha,Vertices,VerticesFlags,dimension);

        //Alpha = 1.5;

        MesherUtilities MesherUtils;

        //2.- to control the alpha size
        if( accepted )
        {
          //boundary_accepted = true;
          if(mrRemesh.Options.Is(MesherUtilities::CONTACT_SEARCH))
          {
            accepted=MesherUtils.ShrankAlphaShape(Alpha,Vertices,mrRemesh.OffsetFactor,dimension);
          }
          else
          {
            if(mrModelPart.Is(FLUID)){
              accepted=MesherUtils.AlphaShape(Alpha,Vertices,dimension,4.0*mrRemesh.Refine->CriticalRadius);
            }
            else{ //SOLID
              accepted=MesherUtils.AlphaShape(Alpha,Vertices,dimension);
            }
          }
        }

        //3.- to control all nodes from the same subdomain (problem, domain is not already set for new inserted particles on mesher)
        bool self_contact = false;
        if(accepted)
        {
          //alpha_accepted = true;
          ++passed_alpha_shape;

          if(mrRemesh.Options.Is(MesherUtilities::CONTACT_SEARCH))
            self_contact = MesherUtils.CheckSubdomain(Vertices);
        }


        //4.- to control that the element is inside of the domain boundaries
        if(accepted)
        {
          //subdomain_accepted = true;
          if(mrRemesh.Options.Is(MesherUtilities::CONTACT_SEARCH))
          {
            //problems in 3D: be careful
            if(self_contact)
              accepted = MesherUtils.CheckOuterCentre(Vertices,mrRemesh.OffsetFactor, self_contact);
          }
          else
          {
            //accepted=MesherUtils.CheckInnerCentre(Vertices); //problems in 3D: when slivers are released, a boundary is created and the normals calculated, then elements that are inside suddently its center is calculated as outside... // some corrections are needded.
          }
        }

        //5.- to control that the element has a good shape
        if(accepted)
        {
          //center_accepted = true;
          ++passed_inner_outer;
          accepted = this->CheckElementShape(Vertices,VerticesFlags,dimension,number_of_slivers);
        }


        // if all checks have been passed, accept the element
        if(accepted)
        {
          //shape_accepted = true;
          number+=1;
          mrRemesh.PreservedElements[el] = number;
        }
        // else{
        //   std::cout<<" Element ["<<el<<"] with alpha "<<mrRemesh.AlphaParameter<<"("<<Alpha<<")"<<std::endl;
        //   for( unsigned int n=0; n<number_of_vertices; ++n)
        //   {
        //     std::cout<<" ("<<n+1<<"): Id["<<Vertices[n].Id()<<"] PreID["<<mrRemesh.NodalPreIds[Vertices[n].Id()]<<"] "<<Vertices[n].Coordinates()<<std::endl;
        //   }
        //   std::cout<<" (alpha:"<<alpha_accepted<<" subdomain:"<<subdomain_accepted<<" center:"<<center_accepted<<" shape:"<<shape_accepted<<") "<< std::endl;

        // }

      }

      mrRemesh.Info->NumberOfElements=number;

    }

    if( mEchoLevel > 0 ){
      std::cout<<"  [Preserved Elements "<<mrRemesh.Info->NumberOfElements<<"] ("<<mrModelPart.NumberOfElements() <<") :: (slivers detected: "<<number_of_slivers<<") "<<std::endl;
      std::cout<<"  (passed_alpha_shape: "<<passed_alpha_shape<<", passed_inner_outer: "<<passed_inner_outer<<") "<<std::endl;
    }

    if( mrModelPart.IsNot(CONTACT) )
      this->SelectNodesToErase();

    // AF for fluids ???
    // mrRemesh.InputInitializedFlag = false;
    // mMesherUtilities.SetNodes(mrModelPart,mrRemesh);
    // mrRemesh.InputInitializedFlag = true;


    if( mEchoLevel > 0 ){
      // std::cout<<"   Generated_Elements :"<<OutNumberOfElements<<std::endl;
      // std::cout<<"   Passed_AlphaShape  :"<<mrRemesh.Info->NumberOfElements<<std::endl;
      // if(OutNumberOfElements-mrRemesh.Info->NumberOfElements!=0)
      //   std::cout<<" DELETED ELEMENTS "<<std::endl;
      std::cout<<"   SELECT MESH ELEMENTS ("<<mrRemesh.Info->NumberOfElements<<") ]; "<<std::endl;

    }

    KRATOS_CATCH( "" )

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
    return "SelectElementsMesherProcess";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << "SelectElementsMesherProcess";
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

  struct NodalFlags{

    unsigned int  Solid;
    unsigned int  Fluid;
    unsigned int  Rigid;
    unsigned int  Boundary;
    unsigned int  FreeSurface;
    unsigned int  NoWallFreeSurface;
    unsigned int  Contact;
    unsigned int  Inlet;
    unsigned int  Isolated;
    unsigned int  Sliver;
    unsigned int  NewEntity;
    unsigned int  OldEntity;

    double Radius;

    //constructor
    NodalFlags()
    {
      Solid = 0;
      Fluid = 0;
      Rigid = 0;
      Boundary = 0;
      FreeSurface = 0;
      NoWallFreeSurface = 0;
      Contact = 0;
      Inlet = 0;
      Isolated = 0;
      Sliver = 0;
      NewEntity = 0;
      OldEntity = 0;
      Radius = 0;
    }

    //counter method
    void CountFlags(const Node<3>& rNode)
    {
      if(rNode.Is(SOLID))
        ++Solid;
      if(rNode.Is(FLUID))
        ++Fluid;
      if(rNode.Is(RIGID))
        ++Rigid;
      if(rNode.Is(BOUNDARY))
        ++Boundary;
      if(rNode.Is(CONTACT))
        ++Contact;
      if(rNode.Is(INLET))
        ++Inlet;
      if(rNode.Is(ISOLATED))
        ++Isolated;
      if(rNode.Is(FREE_SURFACE)){
        ++FreeSurface;
        if(rNode.IsNot(SOLID) && rNode.IsNot(RIGID))
          ++NoWallFreeSurface;
      }
      if(rNode.Is(SELECTED))
        ++Sliver;
      if(rNode.Is(NEW_ENTITY))
        ++NewEntity;
      if(rNode.Is(OLD_ENTITY))
        ++OldEntity;

      Radius+=rNode.FastGetSolutionStepValue(NODAL_H);
    }

  };

  ///@}
  ///@name Protected Operators
  ///@{
  ///@}
  ///@name Protected Operations
  ///@{

  //*******************************************************************************************
  //*******************************************************************************************

  bool IsFirstTimeStep()
  {
    const ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();
    const double& CurrentTime = rCurrentProcessInfo[TIME];
    const double& TimeStep = rCurrentProcessInfo[DELTA_TIME];
    if(CurrentTime<=TimeStep)
      return true;
    else
      return false;
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void LabelEdgeNodes(ModelPart& rModelPart)
  {

    //reset domain flags in nodes before new assignment
    if( rModelPart.Is(FLUID) ){
      //reset domain flags in nodes before new assignment
      unsigned int count_rigid = 0;
      unsigned int count_free_surface = 0;
      unsigned int count_positive_pressure = 0;
      for(ModelPart::ElementsContainerType::iterator i_elem = rModelPart.ElementsBegin() ; i_elem != rModelPart.ElementsEnd() ; ++i_elem)
      {
        PointsArrayType& vertices=i_elem->GetGeometry().Points();

        count_rigid = 0;
        count_free_surface = 0;
        count_positive_pressure = 0;
        for(unsigned int i=0; i<vertices.size(); ++i)
        {
          if( vertices[i].Is(RIGID) )
            ++count_rigid;
        }

        //to release full rigid elements with no fluid surrounding them
        if( count_rigid == vertices.size() ){

          for(unsigned int i=0; i<vertices.size(); ++i)
          {
            vertices[i].Set(OLD_ENTITY,true);
          }
        }
      }


      //set new element and domain flags to nodes
      for(ModelPart::ElementsContainerType::iterator i_elem = rModelPart.ElementsBegin() ; i_elem != rModelPart.ElementsEnd() ; ++i_elem)
      {
        PointsArrayType& vertices=i_elem->GetGeometry().Points();
        count_rigid = 0;
        for(unsigned int i=0; i<vertices.size(); ++i)
        {
          if( vertices[i].Is(RIGID) )
            ++count_rigid;
        }

        //to release full rigid elements with no fluid surrounding them
        if( count_rigid != vertices.size() ){

          for(unsigned int i=0; i<vertices.size(); ++i)
          {
            //set domain type to nodes
            vertices[i].Set(OLD_ENTITY,false);
          }

        }
      }
    }


  }

  //*******************************************************************************************
  //*******************************************************************************************

  void SelectNodesToErase()
  {

    // Set disconnected nodes to erase
    unsigned int dimension = 0;
    unsigned int number_of_vertices = 0;
    unsigned int isolated_nodes = 0;

    this->GetElementDimension(dimension,number_of_vertices);

    int* OutElementList = mrRemesh.OutMesh.GetElementList();

    ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes();

    const int& OutNumberOfElements = mrRemesh.OutMesh.GetNumberOfElements();

    //check engaged nodes
    for(int el=0; el<OutNumberOfElements; ++el)
    {
      if( mrRemesh.PreservedElements[el] ){
        for(unsigned int pn=0; pn<number_of_vertices; ++pn)
        {
          //set vertices
          rNodes[OutElementList[el*number_of_vertices+pn]].Set(BLOCKED);
        }
      }

    }

    int count_release = 0;
    for(ModelPart::NodesContainerType::iterator i_node = rNodes.begin() ; i_node != rNodes.end() ; ++i_node)
    {
      if( i_node->IsNot(BLOCKED) ){

        if(mrModelPart.Is(FLUID)){

          if( i_node->Is(RIGID) || i_node->Is(SOLID) ){
            if( i_node->Is(TO_ERASE) ){
              i_node->Set(TO_ERASE,false);
              std::cout<<" WARNING TRYING TO DELETE A WALL NODE (fluid): "<<i_node->Id()<<std::endl;
            }
          }
          else{

            if(mrRemesh.ExecutionOptions.Is(MesherUtilities::KEEP_ISOLATED_NODES)){

              if( i_node->IsNot(TO_ERASE) ){
                i_node->Set(ISOLATED,true);
                ++isolated_nodes;
              }
              // else{
              //   std::cout<<" ISOLATED TO ERASE "<<std::endl;
              // }

            }
            else{

              //std::cout<<" ISOLATED NODES NOT KEPT "<<std::endl;
              i_node->Set(TO_ERASE,true);
              if( mEchoLevel > 0 )
                std::cout<<" NODE "<<i_node->Id()<<" IS INSIDE RELEASE "<<std::endl;
              if( i_node->Is(BOUNDARY) )
                std::cout<<" NODE "<<i_node->Id()<<" IS BOUNDARY RELEASE "<<std::endl;
              ++count_release;

            }

          }

        }
        else{

          if( i_node->Is(RIGID) ){

            if( i_node->Is(TO_ERASE) ){
              i_node->Set(TO_ERASE,false);
              std::cout<<" WARNING TRYING TO DELETE A WALL NODE (solid): "<<i_node->Id()<<std::endl;
            }

          }
          else{

            if(mrRemesh.ExecutionOptions.Is(MesherUtilities::KEEP_ISOLATED_NODES)){

              if( i_node->IsNot(TO_ERASE) ){
                i_node->Set(ISOLATED,true);
                ++isolated_nodes;
              }

            }
            else{

              i_node->Set(TO_ERASE);
              if( mEchoLevel > 0 )
                std::cout<<" NODE "<<i_node->Id()<<" IS INSIDE RELEASE "<<std::endl;
              if( i_node->Is(BOUNDARY) )
                std::cout<<" NODE "<<i_node->Id()<<" IS BOUNDARY RELEASE "<<std::endl;
              ++count_release;
            }
          }
        }

      }
      else{

        i_node->Set(TO_ERASE,false);
        i_node->Set(ISOLATED,false);

      }

      i_node->Set(BLOCKED,false);
      i_node->Set(OLD_ENTITY,false);
      if(i_node->Is(VISITED))
        i_node->Set(SELECTED,true);
      else
        i_node->Set(SELECTED,false);
      i_node->Set(VISITED,false);

    }

    if( mEchoLevel > 0 ){
      std::cout<<"   NUMBER OF RELEASED NODES "<<count_release<<std::endl;
      std::cout<<"   NUMBER OF ISOLATED NODES "<<isolated_nodes<<std::endl;
    }
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void CheckIds(const int* OutElementList, const int& OutNumberOfElements, const unsigned int number_of_vertices)
  {

    int max_out_id = 0;
    for(int el=0; el<OutNumberOfElements; ++el)
    {
      for(unsigned int pn=0; pn<number_of_vertices; ++pn)
      {
        unsigned int id = el*number_of_vertices+pn;
        if( max_out_id < OutElementList[id] )
          max_out_id = OutElementList[id];
      }
    }

    if( max_out_id >= mrRemesh.NodalPreIds.size() )
      std::cout<<" ERROR ID PRE IDS "<<max_out_id<<" > "<<mrRemesh.NodalPreIds.size()<<" (nodes size:"<<mrModelPart.Nodes().size()<<")"<<std::endl;

  }

  //*******************************************************************************************
  //*******************************************************************************************

  void GetElementDimension(unsigned int& dimension, unsigned int& number_of_vertices)
  {
    if( mrModelPart.NumberOfElements() ){
      ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin();
      dimension = element_begin->GetGeometry().WorkingSpaceDimension();
      number_of_vertices = element_begin->GetGeometry().size();
    }
    else if ( mrModelPart.NumberOfConditions() ){
      ModelPart::ConditionsContainerType::iterator condition_begin = mrModelPart.ConditionsBegin();
      dimension = condition_begin->GetGeometry().WorkingSpaceDimension();
      if( dimension == 3 ) //number of nodes of a tetrahedron
        number_of_vertices = 4;
      else if( dimension == 2 ) //number of nodes of a triangle
        number_of_vertices = 3;
    }

  }



  //*******************************************************************************************
  //*******************************************************************************************

  //it can set TO_ERASE to new entities inserted in full rigid elements
  bool CheckElementBoundaries(GeometryType& rVertices,const NodalFlags& rVerticesFlags)
  {
    bool accepted = true;
    unsigned int NumberOfVertices = rVertices.size();

    if ( mrModelPart.Is(FLUID) ){

      //do not accept full rigid elements (no fluid)
      if( rVerticesFlags.Rigid == NumberOfVertices && rVerticesFlags.Fluid>0){
        //accept when it has more than two fluid nodes (2D) and more than 3 fluid nodes (3D)
        if( rVerticesFlags.Fluid < NumberOfVertices-1 ){
          //std::cout<<" RIGID EDGE DISCARDED (free_surface: "<<rVerticesFlags.FreeSurface<<" fluid: "<<rVerticesFlags.Fluid<<" rigid: "<<rVerticesFlags.Rigid<<")"<<std::endl;
          accepted=false;
        }
        else if( rVerticesFlags.Fluid == NumberOfVertices && rVerticesFlags.OldEntity >= 2  ){
          //std::cout<<" OLD RIGID EDGE DISCARDED (old_entity: "<<rVerticesFlags.OldEntity<<" fluid: "<<rVerticesFlags.Fluid<<" rigid: "<<rVerticesFlags.Rigid<<")"<<std::endl;
          accepted=false;
        }
      }

      //do not accept full rigid elements with a new inserted node (no fluid)
      if( (rVerticesFlags.Rigid + rVerticesFlags.NewEntity) == NumberOfVertices && rVerticesFlags.Fluid>0){
        if( rVerticesFlags.Fluid == NumberOfVertices && rVerticesFlags.OldEntity >= NumberOfVertices-1 ){
          std::cout<<" OLD RIGID NEW ENTITY EDGE DISCARDED (old_entity: "<<rVerticesFlags.OldEntity<<" fluid: "<<rVerticesFlags.Fluid<<" rigid: "<<rVerticesFlags.Rigid<<" free_surface: "<<rVerticesFlags.FreeSurface<<")"<<std::endl;
          accepted=false;
        }
      }

      //do not accept full rigid-solid elements (no fluid)
      if( (rVerticesFlags.Solid + rVerticesFlags.Rigid) >= NumberOfVertices && rVerticesFlags.Fluid == 0)
        accepted=false;

      //do not accept full solid elements
      if( rVerticesFlags.Solid == NumberOfVertices )
        accepted=false;
    }

    return accepted;
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void GetAlphaParameter(double &rAlpha,GeometryType& rVertices,const NodalFlags& rVerticesFlags,const unsigned int& rDimension)
  {
    unsigned int NumberOfVertices = rVertices.size();

    rAlpha = mrRemesh.AlphaParameter;

    if( mrModelPart.Is(SOLID) ){

      if( rVerticesFlags.Boundary >= NumberOfVertices )
        rAlpha*=1.2;

    }
    else if( mrModelPart.Is(FLUID) ){

      MesherUtilities MesherUtils;

      //avoid penetrating/scaping isolated nodes at large speeds and full free surface fluid nodes at large speeds
      if( (rVerticesFlags.Isolated+rVerticesFlags.FreeSurface) == NumberOfVertices && (rVerticesFlags.Isolated > 0 || rVerticesFlags.NoWallFreeSurface == NumberOfVertices ||  (rVerticesFlags.NoWallFreeSurface+rVerticesFlags.Isolated)== NumberOfVertices) ){
        const double MaxRelativeVelocity = 1.5; //arbitrary value, will depend on time step (AF)
        if( MesherUtils.CheckRelativeVelocities(rVertices, MaxRelativeVelocity) ){
          rAlpha=0;
          //set isolated nodes to erase if they approach to fast
          for(unsigned int i=0; i<rVertices.size(); ++i)
          {
            if( rVertices[i].Is(ISOLATED) )
              rVertices[i].Set(TO_ERASE);
          }
        }
      }

      if(rDimension==2){
        this->GetTriangleFluidElementAlpha(rAlpha,rVertices,rVerticesFlags,rDimension);
      }
      else if(rDimension==3){
        this->GetTetrahedronFluidElementAlpha(rAlpha,rVertices,rVerticesFlags,rDimension);
      }
    }

  }

  //*******************************************************************************************
  //*******************************************************************************************

  void GetTriangleFluidElementAlpha(double &rAlpha,GeometryType& rVertices,const NodalFlags& rVerticesFlags,const unsigned int& rDimension)
  {

    MesherUtilities MesherUtils;

    double VolumeChange = 0;
    double VolumeTolerance = 1.15e-4*pow(4.0*mrRemesh.Refine->CriticalRadius,rDimension);
    unsigned int NumberOfVertices = rVertices.size();

    //criterion for elements formed with new wall nodes
    if( rVerticesFlags.Fluid != NumberOfVertices ){

      //there are not fluid nodes
      if( rVerticesFlags.Fluid == 0 ){
        rAlpha = 0;
      }
      else{

        //if the element is decreasing its volume:
        VolumeChange = 0;
        if(MesherUtils.CheckVolumeDecrease(rVertices,rDimension,VolumeTolerance,VolumeChange)){

          //accept new elements formed with isolated nodes (here speedy test passed)
          if( rVerticesFlags.Isolated > 0 ){
            rAlpha*=0.80;
          }
          //one wall vertex
          else if( rVerticesFlags.Rigid == 1 || rVerticesFlags.Solid == 1 ){
            //to avoid new approaching wall elements with two non-wall free-surface nodes
            if( rVerticesFlags.NoWallFreeSurface == 2 )
              rAlpha*=0.60;
            else
              rAlpha*=0.80;
          }
          //two wall vertices
          else if( rVerticesFlags.Rigid == 2 || rVerticesFlags.Solid == 2 || (rVerticesFlags.Rigid+rVerticesFlags.Solid)==2 ){
            //to avoid new approaching wall elements with only one non-wall free-surface node
            if( rVerticesFlags.NoWallFreeSurface == 1 ){
              rAlpha*=0.60;
            }
            else{
              rAlpha*=0.80;
            }
          }
          //there are no wall vertices (impossible)
          else{
            rAlpha*=0.80;
            //std::cout<<" WARNING: new element with non-fluid particles and non wall-particles (rigid: "<<rVerticesFlags.Rigid<<" solid: "<<rVerticesFlags.Solid<<" fluid: "<<rVerticesFlags.Fluid<<" free-surface: "<<rVerticesFlags.FreeSurface<<")"<<std::endl;
          }

        }
        //if the element is not decreasing its volume:
        else{

          //if the element does not change the volume (all rigid or moving tangentially to the wall)
          if( VolumeChange > 0 && VolumeChange < VolumeTolerance ){
            rAlpha*=0.95;
            //std::cout<<" CONSIDERING VOLUME CHANGE "<<VolumeChange<<std::endl;
          }
          else{
            rAlpha=0;
          }
        }
      }
    }
    //all nodes are fluid (pre-existing elements) or new elements formed in the free-surface
    else{ //fluid element

      //all nodes in the free-surface
      if( rVerticesFlags.FreeSurface == 3 ){

        //all nodes in the non-wall free-surface
        if( rVerticesFlags.NoWallFreeSurface == 3 ){

          //if the element is decreasing its volume:
          VolumeChange = 0;
          //to avoid closing fluid voids with new elements
          if(MesherUtils.CheckVolumeDecrease(rVertices,rDimension,VolumeTolerance,VolumeChange)){
            //to avoid too avoid approaching free surfaces at very high speeds
            const double MaxRelativeVelocity = 1.5;
            if(MesherUtils.CheckRelativeVelocities(rVertices, MaxRelativeVelocity))
              rAlpha=0;
            else
              rAlpha*=0.80;
          }
          //to avoid increasing fluid surface with new elements
          else{
            rAlpha*=0.60;
          }
        }
        //two nodes in non-wall free-surface
        else if( rVerticesFlags.NoWallFreeSurface == 2 ){
          rAlpha*=0.80;
        }
        //one node in non-wall free-surface
        else if( rVerticesFlags.NoWallFreeSurface == 1 ){
          rAlpha*=0.80;
        }
        //one node in non-wall free-surface
        else{
          rAlpha*=1.20;
        }
      }
      //two nodes in the free-surface
      else if( rVerticesFlags.FreeSurface == 2 ){

        //to avoid closing fluid voids with new elements with a wall fluid node
        if( rVerticesFlags.NoWallFreeSurface == 2 && (rVerticesFlags.Rigid == 1 || rVerticesFlags.Solid == 1) )
          rAlpha*=0.80;
        else
          rAlpha*=1.20;

      }
      else{
        rAlpha*=1.20;
      }
    }

  }


  //*******************************************************************************************
  //*******************************************************************************************

  void GetTetrahedronFluidElementAlpha(double &rAlpha,GeometryType& rVertices,const NodalFlags& rVerticesFlags,const unsigned int& rDimension)
  {
    MesherUtilities MesherUtils;

    double VolumeChange = 0;
    double VolumeTolerance = 2.5e-3*pow(4.0*mrRemesh.Refine->CriticalRadius,rDimension);
    unsigned int NumberOfVertices = rVertices.size();

    //criterion for elements formed with new wall nodes
    if( rVerticesFlags.Fluid != NumberOfVertices ){ //element formed with a wall

      //there are not fluid nodes
      if( rVerticesFlags.Fluid == 0 ){
        rAlpha = 0;
      }
      else{

        //if the element is decreasing its volume:
        VolumeChange = 0;
        if(MesherUtils.CheckVolumeDecrease(rVertices,rDimension,VolumeTolerance,VolumeChange)){

          //accept new elements formed with isolated nodes (here speedy test passed)
          if( rVerticesFlags.Isolated > 0 ){
            rAlpha*=0.8;
          }
          //one wall vertex
          else if( (rVerticesFlags.Rigid == 1 || rVerticesFlags.Solid == 1) ){
            //to avoid new approaching wall elements with three non-wall free-surface nodes
            if( rVerticesFlags.NoWallFreeSurface == 3 )
              rAlpha *= 0.60;
            else
              rAlpha *= 0.70;
          }
          //two wall vertices
          else if( rVerticesFlags.Rigid == 2 || rVerticesFlags.Solid == 2 || (rVerticesFlags.Rigid+rVerticesFlags.Solid)==2 ){
            //to avoid new approaching wall elements with two non-wall free-surface nodes
            if( rVerticesFlags.NoWallFreeSurface == 2 )
              rAlpha *= 0.50;
            else
              rAlpha *= 0.80;
          }
          //three wall vertices
          else if( rVerticesFlags.Rigid == 3 || rVerticesFlags.Solid == 3 || (rVerticesFlags.Rigid+rVerticesFlags.Solid)==3 ){
            //to avoid new approaching wall elements with only one non-wall free-surface node
            if( rVerticesFlags.NoWallFreeSurface == 1 ){
              if( rVerticesFlags.Fluid == 1 )
                rAlpha *= 0.40;
              else
                rAlpha *= 0.60;
            }
            else{
              rAlpha *= 0.80;
            }
          }
          //there are no wall vertices (impossible)
          else{
            rAlpha*=0.80;
            std::cout<<" WARNING: new element with non-fluid particles and non wall-particles (rigid: "<<rVerticesFlags.Rigid<<" solid: "<<rVerticesFlags.Solid<<" fluid: "<<rVerticesFlags.Fluid<<" free-surface: "<<rVerticesFlags.FreeSurface<<" new_entity: "<<rVerticesFlags.NewEntity<<" isolated: "<<rVerticesFlags.Isolated<<" old_entity: "<<rVerticesFlags.OldEntity<<")"<<std::endl;
          }

        }
        //if the element is not decreasing its volume:
        else{

          //if the element does not change the volume (all rigid or moving tangentially to the wall)
          if( VolumeChange > 0 && VolumeChange < VolumeTolerance && rVerticesFlags.Rigid == NumberOfVertices ){
            rAlpha*=0.80;
            //std::cout<<" CONSIDERING VOLUME CHANGE "<<VolumeChange<<std::endl;
          }
          else{
            rAlpha=0;
          }

        }
      }

      //in 3D alpha must be larger
      rAlpha*=1.10;

    }
    //all nodes are fluid (pre-existing elements) or new elements formed inside of in the fulle free-surface
    else{ //fluid element

      //all nodes in the free-surface
      if( rVerticesFlags.FreeSurface == 4 && rVerticesFlags.Sliver == 0 ){

        //all nodes in the non-wall free-surface (free surface flat element can be artificial)
        if( rVerticesFlags.NoWallFreeSurface == 4){

          //if the element is decreasing its volume:
          VolumeChange = 0;
          //to avoid closing fluid voids with new elements
          if(MesherUtils.CheckVolumeDecrease(rVertices,rDimension,VolumeTolerance,VolumeChange))
            rAlpha*=0.70;
          //to avoid increasing fluid surface with new elements
          else
            rAlpha*=0.40;

        }
        //three nodes in non-wall free-surface
        else if( rVerticesFlags.NoWallFreeSurface == 3 ){
          if(MesherUtils.CheckVolumeDecrease(rVertices,rDimension,VolumeTolerance,VolumeChange))
            rAlpha*=0.70;
          else
            rAlpha*=0.50;
        }
        //two nodes in non-wall free-surface
        else if( rVerticesFlags.NoWallFreeSurface == 2 ){
          if(MesherUtils.CheckVolumeDecrease(rVertices,rDimension,VolumeTolerance,VolumeChange))
            rAlpha*=0.70;
          else
            rAlpha*=0.50;
        }
        //one node in non-wall free-surface
        else{
          rAlpha*=0.80;
        }
      }
      else{

        if( rVerticesFlags.FreeSurface > 0 && rVerticesFlags.Rigid > 0 ){
          if( MesherUtils.CheckVolumeDecrease(rVertices,rDimension,VolumeTolerance,VolumeChange) )
            rAlpha*=1.50;
          else
            rAlpha*=1.10;
        }
        else if(rVerticesFlags.FreeSurface == 0 &&  rVerticesFlags.Rigid == 3 ){
          rAlpha*=2.0; //accept large elements
        }
        else{
          rAlpha*=1.70;
        }
      }

      //in 3D alpha must be larger
      rAlpha*=1.10;

    }

  }

  //*******************************************************************************************
  //*******************************************************************************************

  bool CheckElementShape(GeometryType& rVertices, const NodalFlags& rVerticesFlags,const unsigned int& rDimension, unsigned int& number_of_slivers)
  {
    bool accepted = true;
    unsigned int NumberOfVertices = rVertices.size();

    if( mrModelPart.Is(SOLID) ){

      int sliver = 0;
      if(rDimension == 3 && NumberOfVertices==4){

        Tetrahedra3D4<Node<3> > Tetrahedron(rVertices);

        MesherUtilities MesherUtils;
        accepted = MesherUtils.CheckGeometryShape(Tetrahedron,sliver);

        if( sliver ){
          if(mrRemesh.Options.Is(MesherUtilities::CONTACT_SEARCH))
            accepted = true;
          else
            accepted = false;

          ++number_of_slivers;
        }
        else{

          if(mrRemesh.Options.Is(MesherUtilities::CONTACT_SEARCH))
            accepted = false;
          else
            accepted = true;
        }

      }

    }
    else if ( mrModelPart.Is(FLUID) ){

      if(rDimension == 2 && NumberOfVertices==3){

        //check 2D distorted elements with edges decreasing very fast
        if( rVerticesFlags.NoWallFreeSurface > 1 && (rVerticesFlags.Rigid+rVerticesFlags.Solid)==0 ){

          double MaxEdgeLength = std::numeric_limits<double>::min();
          double MinEdgeLength = std::numeric_limits<double>::max();

          for(unsigned int i=0; i<NumberOfVertices-1; ++i)
          {
            for(unsigned int j=i+1; j<NumberOfVertices; ++j)
            {
              double Length = norm_2(rVertices[j].Coordinates() - rVertices[i].Coordinates());
              if( Length < MinEdgeLength ){
                MinEdgeLength = Length;
              }
              if( Length > MaxEdgeLength ){
                MaxEdgeLength = Length;
              }

            }
          }
          MesherUtilities MesherUtils;
          const double MaxRelativeVelocity = 1.5; //arbitrary value, will depend on time step (AF)
          if( MinEdgeLength*5 < MaxEdgeLength ){
            if( MesherUtils.CheckRelativeVelocities(rVertices, MaxRelativeVelocity) ){
              std::cout<<" WARNING 2D sliver "<<std::endl;
              //commented to no erase but labeled to not calculate them
              //accepted = false;
              //only erase full free surface sliver elements (visualization purposes)
              if( rVerticesFlags.FreeSurface == NumberOfVertices)
                accepted = false;
              ++number_of_slivers;
              for(unsigned int i=0; i<NumberOfVertices; ++i)
              {
                if( rVertices[i].Is(VISITED) )
                  std::cout<<" WARNING Second sliver in the same node bis "<<std::endl;
                rVertices[i].Set(VISITED);
              }
            }
          }
        }

      }
      else if(rDimension == 3 && NumberOfVertices==4){

        if( rVerticesFlags.FreeSurface >=2 || (rVerticesFlags.Boundary == 4 && rVerticesFlags.NoWallFreeSurface !=4) ){

          Tetrahedra3D4<Node<3> > Tetrahedron(rVertices);
          double Volume = Tetrahedron.Volume();

          if( Volume < 0.01*pow(4.0*mrRemesh.Refine->CriticalRadius,rDimension) ){
            //KRATOS_INFO("SLIVER")<<" Volume="<<Volume<<" VS Critical Volume="<<0.01*pow(4.0*mrRemesh.Refine->CriticalRadius,rDimension)<<std::endl;
            //std::cout<<" SLIVER Volume="<<Volume<<" VS Critical Volume="<< 0.01*pow(4.0*mrRemesh.Refine->CriticalRadius,rDimension) <<std::endl;

            //commented to no erase but labeled to not calculate them
            //accepted = false;
            //only erase full free surface sliver elements (visualization purposes)
            if( rVerticesFlags.FreeSurface == NumberOfVertices)
              accepted = false;

            ++number_of_slivers;

            for(unsigned int i=0; i<NumberOfVertices; ++i)
            {
              // if( rVertices[i].Is(VISITED) )
              //   std::cout<<" WARNING Second sliver in the same node "<<std::endl;
              rVertices[i].Set(VISITED);
            }

          }
          // else if( Volume < 0.02*pow(4.0*mrRemesh.Refine->CriticalRadius,rDimension) ){
          //   int sliver = 0;
          //   MesherUtilities MesherUtils;
          //   bool distorted = MesherUtils.CheckGeometryShape(Tetrahedron,sliver);

          //   //std::cout << " SLIVER " << sliver <<" DISTORTED "<< distorted << std::endl;

          //   if( sliver ){
          //     std::cout<<" SLIVER SHAPE Volume="<<Volume<<" VS Critical Volume="<< 0.02*pow(4.0*mrRemesh.Refine->CriticalRadius,rDimension) <<std::endl;
          //     accepted = false;
          //     ++number_of_slivers;

          //     for(unsigned int i=0; i<rVertices.size(); ++i)
          //     {
          //       rVertices[i].Set(VISITED);
          //     }

          //   }

          // }
        }
        else if( ( rVerticesFlags.Rigid == 1 || rVerticesFlags.Rigid == 2 )  ){

          MesherUtilities MesherUtils;

          double VolumeIncrement = MesherUtils.GetDeformationGradientDeterminant(rVertices,rDimension);

          double MovedVolume = MesherUtils.GetMovedVolume(rVertices,rDimension,1.0);

          //check if the element is going to invert
          if( VolumeIncrement < 0.0 || MovedVolume < 0.0 ){

            //std::cout<<" SLIVER FOR INVERSION "<<VolumeIncrement<<" VS Critical Volume="<< 0.01*pow(4.0*mrRemesh.Refine->CriticalRadius,rDimension) <<"( rigid: "<<rVerticesFlags.Rigid<<" )"<<std::endl;

            //commented to no erase but labeled to not calculate them
            //accepted = false;

            ++number_of_slivers;

            for(unsigned int i=0; i<NumberOfVertices; ++i)
            {
              // if( rVertices[i].Is(VISITED) )
              //   std::cout<<" WARNING Second sliver in the same node bis "<<std::endl;
              rVertices[i].Set(VISITED);
            }
          }

        }

      }


    }
    return accepted;

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
  SelectElementsMesherProcess& operator=(SelectElementsMesherProcess const& rOther);

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
                                  SelectElementsMesherProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SelectElementsMesherProcess& rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_SELECT_ELEMENTS_MESHER_PROCESS_H_INCLUDED defined
