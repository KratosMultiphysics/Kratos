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

        accepted = this->CheckElementBoundaries(Vertices,VerticesFlags);

        //std::cout<<" ******** ELEMENT "<<el+1<<" ********** "<<std::endl;

        double Alpha = mrRemesh.AlphaParameter;

        this->GetAlphaParameter(Alpha,Vertices,VerticesFlags,dimension);

        // std::cout<<" vertices for the contact element "<<std::endl;
        // if(mrRemesh.Options.Is(MesherUtilities::CONTACT_SEARCH)){
        // 	for( unsigned int n=0; n<number_of_vertices; ++n)
        // 	  {
        // 	    std::cout<<" ("<<n+1<<"): ["<<mrRemesh.NodalPreIds[Vertices[n].Id()]<<"] "<<std::endl;
        // 	  }
        // }
        // std::cout<<" vertices for the subdomain element "<<std::endl;
        // for( unsigned int n=0; n<number_of_vertices; ++n)
        //  	{
        //  	  std::cout<<" ("<<n+1<<"): ["<<Vertices[n].Id()<<"]  NodalH "<<Vertices[n].FastGetSolutionStepValue(NODAL_H)<<std::endl;
        //  	}
        //std::cout<<" Element "<<el<<" with alpha "<<mrRemesh.AlphaParameter<<"("<<Alpha<<")"<<std::endl;


        MesherUtilities MesherUtils;

        if( accepted )
        {
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
        // if(accepted)
        // {
        //   std::cout<<" Element passed Alpha Shape "<<std::endl;
        //     if(mrRemesh.Refine->Options.IsNot(MesherUtilities::CONTACT_SEARCH))
        //   	accepted=MesherUtilities::CheckSubdomain(Vertices);
        // }

        //3.1.-
        bool self_contact = false;
        if(mrRemesh.Options.Is(MesherUtilities::CONTACT_SEARCH))
          self_contact = MesherUtils.CheckSubdomain(Vertices);

        //4.- to control that the element is inside of the domain boundaries
        if(accepted)
        {
          ++passed_alpha_shape;

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
        // else{

        //  	for( unsigned int n=0; n<number_of_vertices; ++n)
        // 	  {
        //  	    std::cout<<" ("<<n+1<<"): ["<<Vertices[n].Id()<<"]  NodalH "<<Vertices[n].FastGetSolutionStepValue(NODAL_H)<<std::endl;
        // 	  }

        // 	std::cout<<" Element "<<el<<" with alpha "<<mrRemesh.AlphaParameter<<"("<<Alpha<<")"<<std::endl;

        // }

        //5.- to control that the element has a good shape
        if(accepted)
        {
          ++passed_inner_outer;
          accepted = this->CheckElementShape(Vertices,VerticesFlags,dimension,number_of_slivers);
        }


        if(accepted)
        {
          //std::cout<<" Element ACCEPTED after cheking Center+Sliver "<<number<<std::endl;
          number+=1;
          mrRemesh.PreservedElements[el] = number;
        }
        // else{

        //   std::cout<<" Element DID NOT pass INNER/OUTER check "<<std::endl;
        // }


      }

      mrRemesh.Info->NumberOfElements=number;

    }

    std::cout<<"   [Preserved Elements "<<mrRemesh.Info->NumberOfElements<<"] :: (slivers detected: "<<number_of_slivers<<") "<<std::endl;
    std::cout<<"   (passed_alpha_shape: "<<passed_alpha_shape<<", passed_inner_outer: "<<passed_inner_outer<<") "<<std::endl;

    if(mrRemesh.ExecutionOptions.IsNot(MesherUtilities::KEEP_ISOLATED_NODES)){


      unsigned int dimension = 0;
      unsigned int number_of_vertices = 0;

      this->GetElementDimension(dimension,number_of_vertices);

      int* OutElementList = mrRemesh.OutMesh.GetElementList();

      ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes();

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
        if( i_node->IsNot(BLOCKED)  ){
          if(!(i_node->Is(FREE_SURFACE) || i_node->Is(RIGID))){
            i_node->Set(TO_ERASE);
            if( mEchoLevel > 0 )
              std::cout<<" NODE "<<i_node->Id()<<" RELEASE "<<std::endl;
            if( i_node->Is(BOUNDARY) )
              std::cout<<" ERROR: node "<<i_node->Id()<<" IS BOUNDARY RELEASE "<<std::endl;
            ++count_release;
          }
        }

        i_node->Reset(BLOCKED);
      }

      if( mEchoLevel > 0 )
        std::cout<<"   NUMBER OF RELEASED NODES "<<count_release<<std::endl;

    }
    else{

      ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes();

      for(ModelPart::NodesContainerType::iterator i_node = rNodes.begin() ; i_node != rNodes.end() ; ++i_node)
      {
        i_node->Reset(BLOCKED);
      }
    }

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
    unsigned int  Contact;
    unsigned int  Inlet;
    unsigned int  Isolated;

    double Radius;

    //constructor
    NodalFlags()
    {
      Solid = 0;
      Fluid = 0;
      Rigid = 0;
      Boundary = 0;
      FreeSurface = 0;
      Contact = 0;
      Inlet = 0;
      Isolated = 0;
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
      if(rNode.Is(FREE_SURFACE))
        ++FreeSurface;
      if(rNode.Is(CONTACT))
        ++Contact;
      if(rNode.Is(INLET))
        ++Inlet;
      if(rNode.Is(ISOLATED))
        ++Isolated;

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

      //avoid penetrating/scaping isolated nodes at large speeds
      if( (rVerticesFlags.Isolated+rVerticesFlags.FreeSurface) == NumberOfVertices && rVerticesFlags.Isolated > 0 ){
        const double MaxRelativeVelocity = 1.5; //arbitrary value, will depend on time step (AF)
        if( MesherUtils.CheckRelativeVelocities(rVertices, MaxRelativeVelocity) )
          rAlpha=0;
      }

      if(rDimension==2){

        if( rVerticesFlags.Fluid != NumberOfVertices ){ //element formed with a wall

          if( rVerticesFlags.Fluid == 0 ){
            rAlpha = 0;
          }
          else{

            if(MesherUtils.CheckVolumeDecrease(rVertices,rDimension)){

              //one wall vertex or only one free surface vertex for a new element
              if( (rVerticesFlags.Rigid == 1 || rVerticesFlags.Solid == 1) || ((rVerticesFlags.Rigid == 2 || rVerticesFlags.Solid == 2) && rVerticesFlags.FreeSurface == 1) ){
                rAlpha*=0.65;
              }
              else if( (rVerticesFlags.Rigid == 2 || rVerticesFlags.Solid == 2) && rVerticesFlags.Isolated == 1){
                rAlpha*=1.20;
              }
              else{
                rAlpha*=0.75;
              }

            }
            else{
              rAlpha=0;
              //std::cout<<" VOLUME release [isolated:"<<rVerticesFlags.Isolated<<" rigid:"<<rVerticesFlags.Rigid<<" freesurface:"<<rVerticesFlags.FreeSurface<<"]"<<std::endl;
            }
          }
        }
        else{ //fluid element

          if( rVerticesFlags.FreeSurface == 3 ){
           if( (rVerticesFlags.Solid+rVerticesFlags.Rigid)==0 ){
              if(MesherUtils.CheckVolumeDecrease(rVertices,rDimension))
                rAlpha*=0.65;
              else
                rAlpha=0;
            }
            else{
              rAlpha*=0.85;
            }
          }
          else if( (rVerticesFlags.Rigid == 1 || rVerticesFlags.Solid == 1) && rVerticesFlags.FreeSurface >= 2 ){
            rAlpha*=1.05;
          }
          else if( (rVerticesFlags.Rigid == 2 || rVerticesFlags.Solid == 2) && rVerticesFlags.FreeSurface >= 1 ){
            rAlpha*=1.10;
          }
          else{
            rAlpha*=1.75;
          }
        }

      }
      else if(rDimension==3){

        if( rVerticesFlags.Fluid != NumberOfVertices ){ //element formed with a wall

          if( rVerticesFlags.Fluid == 0 ){
            rAlpha = 0;
          }
          else{

            if(MesherUtils.CheckVolumeDecrease(rVertices,rDimension)){

              //one wall vertex or only one free surface vertex for a new element
              if( (rVerticesFlags.Rigid == 1 || rVerticesFlags.Solid == 1) || ((rVerticesFlags.Rigid == 3 || rVerticesFlags.Solid == 3) && rVerticesFlags.FreeSurface == 1) )
                rAlpha *= 0.75;
              else
                rAlpha*=0.95;

            }
            else{
              rAlpha = 0;
            }
          }
        }
        else{ //fluid element

          if( rVerticesFlags.FreeSurface == 4 ){
            if( (rVerticesFlags.Solid+rVerticesFlags.Rigid)==0 ){
              if(MesherUtils.CheckVolumeDecrease(rVertices,rDimension))
                rAlpha*=0.75;
              else
                rAlpha=0;
            }
            else{
              rAlpha*=0.85;
            }
          }
          else if( (rVerticesFlags.Rigid == 1 || rVerticesFlags.Solid == 1) && rVerticesFlags.FreeSurface == 3 ){
            rAlpha*=0.85;
          }
          else if( (rVerticesFlags.Rigid == 3 || rVerticesFlags.Solid == 3) && rVerticesFlags.FreeSurface == 3 ){
            rAlpha*=1.05;
          }
          else if( (rVerticesFlags.Rigid == 2 || rVerticesFlags.Solid == 2) && rVerticesFlags.FreeSurface == 2 ){
            rAlpha*=1.10;
          }
          else{
            rAlpha*=1.75;
          }
        }

      }

    }

  }


  //*******************************************************************************************
  //*******************************************************************************************

  bool CheckElementBoundaries(const GeometryType& rVertices,const NodalFlags& rVerticesFlags)
  {
    bool accepted = true;
    unsigned int NumberOfVertices = rVertices.size();

    if ( mrModelPart.Is(FLUID) ){

      //do not accept full rigid elements (no fluid)
      if( rVerticesFlags.Rigid == NumberOfVertices && rVerticesFlags.Fluid < NumberOfVertices )
        accepted=false;

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

  bool CheckElementShape(const GeometryType& rVertices,const NodalFlags& rVerticesFlags,const unsigned int& rDimension, unsigned int& number_of_slivers)
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

      if( rVerticesFlags.FreeSurface > 0 || rVerticesFlags.Boundary == 4 ){

        if(rDimension == 3 && NumberOfVertices==4){

          Tetrahedra3D4<Node<3> > Tetrahedron(rVertices);
          double Volume = Tetrahedron.Volume();

          if( Volume < 0.01*mrRemesh.Refine->MeanVolume ){
            KRATOS_INFO("SLIVER")<<" Volume="<<Volume<<" VS Critical Volume="<<0.01*mrRemesh.Refine->MeanVolume<<std::endl;
            accepted = false;
            ++number_of_slivers;
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
