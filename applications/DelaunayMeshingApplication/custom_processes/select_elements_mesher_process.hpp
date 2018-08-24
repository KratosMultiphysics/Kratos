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

#include "includes/model_part.h"
#include "custom_utilities/mesher_utilities.hpp"
#include "geometries/tetrahedra_3d_4.h"
#include "custom_processes/mesher_process.hpp"

///VARIABLES used:
//Data:
//StepData: NODAL_H, CONTACT_FORCE
//Flags:    (checked) TO_ERASE, BOUNDARY, NEW_ENTITY
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

      int& OutNumberOfElements = mrRemesh.OutMesh.GetNumberOfElements();
      mrRemesh.PreservedElements.clear();
      mrRemesh.PreservedElements.resize(OutNumberOfElements);
      std::fill( mrRemesh.PreservedElements.begin(), mrRemesh.PreservedElements.end(), 0 );
      mrRemesh.MeshElementsSelectedFlag = true;

      mrRemesh.Info->NumberOfElements=0;

      bool box_side_element = false;
      bool wrong_added_node = false;

      int number_of_slivers = 0;

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

	  unsigned int nds = 3;
	  unsigned int dimension = 2;
	  if( mrModelPart.NumberOfElements() ){
	    ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin();
	    nds = element_begin->GetGeometry().size();
	    dimension = element_begin->GetGeometry().WorkingSpaceDimension();
	  }
	  else if ( mrModelPart.NumberOfConditions() ){
	    ModelPart::ConditionsContainerType::iterator condition_begin = mrModelPart.ConditionsBegin();
	    dimension = condition_begin->GetGeometry().WorkingSpaceDimension();
	    if( dimension == 3 ) //number of nodes of a tetrahedron
	      nds = 4;
	    else if( dimension == 2 ) //number of nodes of a triangle
	      nds = 3;
	  }

	  int* OutElementList = mrRemesh.OutMesh.GetElementList();

	  ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes();

	  int el = 0;
	  int number = 0;

	  //CHECK:
	  // int max_out_id = 0;
	  // for(el=0; el<OutNumberOfElements; ++el)
	  //   {
	  //     for(unsigned int pn=0; pn<nds; ++pn)
	  // 	{
	  // 	  if( max_out_id < OutElementList[el*nds+pn] )
	  // 	    max_out_id = OutElementList[el*nds+pn];
	  // 	}
	  //   }

	  // std::cout<<"   MaxOutID "<<max_out_id<<std::endl;
	  // std::cout<<"   NumberOfNodes "<<rNodes.size()<<std::endl;
	  // std::cout<<"   NodalPreIdsSize "<<mrRemesh.NodalPreIds.size()<<std::endl;

	  // if( max_out_id >= mrRemesh.NodalPreIds.size() )
	  //   std::cout<<" ERROR ID PRE IDS "<<max_out_id<<" > "<<mrRemesh.NodalPreIds.size()<<std::endl;

	  //#pragma omp parallel for reduction(+:number) private(el)
	  for(el=0; el<OutNumberOfElements; ++el)
	    {
	      Geometry<Node<3> > vertices;
	      //double Alpha   = 0;
	      //double nodal_h = 0;
	      //double param   = 0.3333333;

	      // int  numflying=0;
	      // int  numlayer =0;
	      // int  numfixed =0;

	      unsigned int  numfreesurf =0;
	      unsigned int  numboundary =0;

	      // std::cout<<" selected vertices ["<<OutElementList[el*nds];
	      // for(unsigned int d=1; d<nds; ++d)
	      // 	{
	      // 	  std::cout<<", "<<OutElementList[el*nds+d];
	      // 	}

	      // std::cout<<"] "<<std::endl;
	      wrong_added_node = false;
	      box_side_element = false;
	      for(unsigned int pn=0; pn<nds; ++pn)
		{
		  //std::cout<<" pn "<<pn<<" id "<<OutElementList[el*nds+pn]<<" size "<<rNodes.size()<<" IDS "<<mrRemesh.NodalPreIds.size()<<" preid "<<mrRemesh.NodalPreIds[OutElementList[el*nds+pn]]<<std::endl;

		  if(OutElementList[el*nds+pn]<=0)
		    std::cout<<" ERROR: something is wrong: nodal id < 0 "<<el<<std::endl;

		  //check if the number of nodes are considered in the nodal pre ids
		  if( (unsigned int)OutElementList[el*nds+pn] >= mrRemesh.NodalPreIds.size() ){
		      if(mrRemesh.Options.Is(MesherUtilities::CONTACT_SEARCH))
			  wrong_added_node = true;
		      std::cout<<" ERROR: something is wrong: node out of bounds "<<std::endl;
		      break;
		  }

		  //check if is a vertex of an artificial external bounding box
		  if(mrRemesh.NodalPreIds[OutElementList[el*nds+pn]]<0){
		      if(mrRemesh.Options.IsNot(MesherUtilities::CONTACT_SEARCH))
			  std::cout<<" ERROR: something is wrong: nodal id < 0 "<<std::endl;
		      box_side_element = true;
		      break;
		  }

		  //vertices.push_bac( *((rNodes).find( OutElementList[el*nds+pn] ).base() ) );
		  vertices.push_back(rNodes(OutElementList[el*nds+pn]));

		  //check flags on nodes
		  if(vertices.back().Is(FREE_SURFACE))
		    numfreesurf++;

		  if(vertices.back().Is(BOUNDARY))
		    numboundary++;

		  // if(VertexPa[pn].match(_wall_))
		  // 	numfixed++;

		  // if(VertexPa[pn].match(_flying_))
		  // 	numflying++;

		  // if(VertexPa[pn].match(_layer_))
		  // 	numlayer++;

		  //nodal_h+=vertices.back().FastGetSolutionStepValue(NODAL_H);

		}


	      if(box_side_element || wrong_added_node){
		//std::cout<<" Box_Side_Element "<<std::endl;
		continue;
	      }

	      //1.- to not consider wall elements
	      // if(numfixed==3)
	      //   Alpha=0;

	      //2.- alpha shape:
	      //Alpha  = nodal_h * param;
	      //Alpha *= mrRemesh.AlphaParameter; //1.4; 1.35;

	      //2.1.- correction to avoid big elements on boundaries
	      // if(numflying>0){
	      //   Alpha*=0.8;
	      // }
	      // else{
	      //   if(numfixed+numsurf<=2){
	      //     //2.2.- correction to avoid voids in the fixed boundaries
	      //     if(numfixed>0)
	      // 	Alpha*=1.4;

	      //     //2.3.- correction to avoid voids on the free surface
	      //     if(numsurf>0)
	      // 	Alpha*=1.3;

	      //     //2.4.- correction to avoid voids in the next layer after fixed boundaries
	      //     if(numlayer>0 && !numsurf)
	      // 	Alpha*=1.2;
	      //   }

	      // }

	      //std::cout<<" ******** ELEMENT "<<el+1<<" ********** "<<std::endl;

	      double Alpha = mrRemesh.AlphaParameter; //*nds;
	      if(numboundary>=nds-1)
		Alpha*=1.2;


	      // std::cout<<" vertices for the contact element "<<std::endl;
	      // if(mrRemesh.Options.Is(MesherUtilities::CONTACT_SEARCH)){
	      // 	for( unsigned int n=0; n<nds; ++n)
	      // 	  {
	      // 	    std::cout<<" ("<<n+1<<"): ["<<mrRemesh.NodalPreIds[vertices[n].Id()]<<"] "<<std::endl;
	      // 	  }
	      // }

	      // std::cout<<" vertices for the subdomain element "<<std::endl;
	      // for( unsigned int n=0; n<nds; ++n)
	      //  	{
	      //  	  std::cout<<" ("<<n+1<<"): ["<<vertices[n].Id()<<"]  NodalH "<<vertices[n].FastGetSolutionStepValue(NODAL_H)<<std::endl;
	      //  	}

	      //std::cout<<" Element "<<el<<" with alpha "<<mrRemesh.AlphaParameter<<"("<<Alpha<<")"<<std::endl;

	      bool accepted=false;

	      MesherUtilities MesherUtils;

	      if(mrRemesh.Options.Is(MesherUtilities::CONTACT_SEARCH))
		{
		  accepted=MesherUtils.ShrankAlphaShape(Alpha,vertices,mrRemesh.OffsetFactor,dimension);
		}
	      else
		{
		  accepted=MesherUtils.AlphaShape(Alpha,vertices,dimension);
		}

	      //3.- to control all nodes from the same subdomain (problem, domain is not already set for new inserted particles on mesher)
	      // if(accepted)
	      // {
	      //   std::cout<<" Element passed Alpha Shape "<<std::endl;
	      //     if(mrRemesh.Refine->Options.IsNot(MesherUtilities::CONTACT_SEARCH))
	      //   	accepted=MesherUtilities::CheckSubdomain(vertices);
	      // }

	      //3.1.-
	      bool self_contact = false;
	      if(mrRemesh.Options.Is(MesherUtilities::CONTACT_SEARCH))
		self_contact = MesherUtils.CheckSubdomain(vertices);

	      //4.- to control that the element is inside of the domain boundaries
	      if(accepted)
		{
		  passed_alpha_shape++;

		  if(mrRemesh.Options.Is(MesherUtilities::CONTACT_SEARCH))
		    {
		      //problems in 3D: be careful
		      if(self_contact)
			  accepted = MesherUtils.CheckOuterCentre(vertices,mrRemesh.OffsetFactor, self_contact);
		    }
		  else
		    {
		      //accepted=MesherUtils.CheckInnerCentre(vertices); //problems in 3D: when slivers are released, a boundary is created and the normals calculated, then elements that are inside suddently its center is calculated as outside... // some corrections are needded.
		    }
		}
	      // else{

	      //  	for( unsigned int n=0; n<nds; ++n)
	      // 	  {
	      //  	    std::cout<<" ("<<n+1<<"): ["<<vertices[n].Id()<<"]  NodalH "<<vertices[n].FastGetSolutionStepValue(NODAL_H)<<std::endl;
	      // 	  }

	      // 	std::cout<<" Element "<<el<<" with alpha "<<mrRemesh.AlphaParameter<<"("<<Alpha<<")"<<std::endl;

	      // }

	      //5.- to control that the element has a good shape
	      int sliver = 0;
	      if(accepted)
		{
		  passed_inner_outer++;

		  if(nds==4){
		    Geometry<Node<3> >* tetrahedron = new Tetrahedra3D4<Node<3> > (vertices);

		    accepted = MesherUtils.CheckGeometryShape(*tetrahedron,sliver);

		    if( sliver ){

		      if(mrRemesh.Options.Is(MesherUtilities::CONTACT_SEARCH))
			accepted = true;
		      else
			accepted = false;

		      number_of_slivers++;
		    }
		    else{

		      if(mrRemesh.Options.Is(MesherUtilities::CONTACT_SEARCH))
			accepted = false;
		      else
			accepted = true;
		    }

		    delete tetrahedron;
		  }
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

	unsigned int nds = 3;
	if( mrModelPart.NumberOfElements() ){
	  ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin();
	  nds = element_begin->GetGeometry().size();
	}
	else if ( mrModelPart.NumberOfConditions() ){
	  ModelPart::ConditionsContainerType::iterator condition_begin = mrModelPart.ConditionsBegin();
	  unsigned int dimension = condition_begin->GetGeometry().WorkingSpaceDimension();
	  if( dimension == 3 ) //number of nodes of a tetrahedron
	    nds = 4;
	  else if( dimension == 2 ) //number of nodes of a triangle
	    nds = 3;
	}

	int* OutElementList = mrRemesh.OutMesh.GetElementList();

	ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes();

	//check engaged nodes
	for(int el=0; el<OutNumberOfElements; ++el)
	  {
	    if( mrRemesh.PreservedElements[el] ){
	      for(unsigned int pn=0; pn<nds; ++pn)
		{
		  //set vertices
		  rNodes[OutElementList[el*nds+pn]].Set(BLOCKED);
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
		count_release++;
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


private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Static Member Variables
    ///@{
    ModelPart& mrModelPart;

    MesherUtilities::MeshingParameters& mrRemesh;

    MesherUtilities mMesherUtilities;

    int mEchoLevel;

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
