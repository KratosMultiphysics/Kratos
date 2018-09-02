//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                August 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_GENERATE_NEW_CONTACT_CONDITIONS_MESHER_PROCESS_H_INCLUDED )
#define  KRATOS_GENERATE_NEW_CONTACT_CONDITIONS_MESHER_PROCESS_H_INCLUDED


// System includes
#include <string>
#include <iostream>

// External includes
#include <boost/timer.hpp>

// Project includes
#include "includes/model_part.h"
#include "custom_utilities/mesher_utilities.hpp"
#include "custom_processes/select_elements_mesher_process.hpp"

///VARIABLES used:
//Data:     MASTER_ELEMENTS(set), MASTER_CONDITION(set)
//StepData:
//Flags:    (checked)
//          (set)     CONTACT(set)
//          (modified)
//          (reset)
// (set):=(set in this process)

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
  /** Detail class definition.
   */
  class GenerateNewContactConditionsMesherProcess
    : public MesherProcess
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GenerateNewContactConditionsMesherProcess
    KRATOS_CLASS_POINTER_DEFINITION( GenerateNewContactConditionsMesherProcess );

    typedef ModelPart::ConditionType         ConditionType;
    typedef ModelPart::PropertiesType       PropertiesType;
    typedef ConditionType::GeometryType       GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GenerateNewContactConditionsMesherProcess(ModelPart& rModelPart,
			     MesherUtilities::MeshingParameters& rRemeshingParameters,
			     int EchoLevel)
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
    {
      mEchoLevel = EchoLevel;
    }

    /// Destructor.
    virtual ~GenerateNewContactConditionsMesherProcess()
    {
    }


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
      Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
      KRATOS_TRY

      if( mEchoLevel > 0 ){
	std::cout<<" [ GENERATE NEW CONTACT ELEMENTS: "<<std::endl;
        std::cout<<"   Total Conditions BEFORE: ["<<mrModelPart.Conditions().size()<<"] ];"<<std::endl;
      }
      
      if( mrModelPart.Name() != mrRemesh.SubModelPartName )
	std::cout<<" ModelPart Supplied do not corresponds to the Meshing Domain: ("<<mrModelPart.Name()<<" != "<<mrRemesh.SubModelPartName<<")"<<std::endl;

      //*******************************************************************
      //selecting elements
      if( !mrRemesh.MeshElementsSelectedFlag )  //Select Mesh Elements not performed  ... is needed to be done before building new elements
	{
	  std::cout<<" ERROR : no selection of elements performed before building the elements "<<std::endl;
	  SelectElementsMesherProcess SelectElements(mrModelPart,mrRemesh,mEchoLevel);
	  SelectElements.Execute();
	}


      //*******************************************************************
      //set consecutive ids for global conditions
      int id=1;
      for(ModelPart::ConditionsContainerType::iterator i_cond = mrModelPart.GetParentModelPart()->ConditionsBegin() ; i_cond != mrModelPart.GetParentModelPart()->ConditionsEnd(); ++i_cond)
	{
	  i_cond->SetId(id);
	  id++;
	}

      //*******************************************************************
      //properties to be used in the generation
      Properties::Pointer pProperties = mrRemesh.GetProperties();

      Condition const & rReferenceCondition=mrRemesh.GetReferenceCondition(); //contact element

      const unsigned int nds = rReferenceCondition.GetGeometry().size();

      if( mEchoLevel > 0 )
        std::cout<<"   [START contact Element Generation "<<std::endl;

      int& OutNumberOfElements = mrRemesh.OutMesh.GetNumberOfElements();
      int* OutElementList      = mrRemesh.OutMesh.GetElementList();

      int previous_id = mrModelPart.GetParentModelPart()->Conditions().back().Id();
      id = previous_id;

      ModelPart::NodesContainerType::iterator nodes_begin = mrModelPart.NodesBegin();

      for(int el = 0; el<OutNumberOfElements; ++el)
	{
	  if(mrRemesh.PreservedElements[el])
	    {
	      Geometry<Node<3> > Vertices;
	      for(unsigned int i=0; i<nds; ++i)
		{
		  //note that OutElementList, starts from node 1, not from node 0, it can be directly assigned to mrRemesh.NodalPreIds.

		  // detected problems in find() method for mesh nodes
		  // bool node_exists = mrModelPart.GetMesh().HasNode(OutElementList[el*nds+i]);
		  // if(!node_exists){

		  //   std::cout<<" ERROR node "<<mrRemesh.NodalPreIds[OutElementList[el*nds+i]]<<" is not in the modelpart "<<std::endl;
		  //   std::cout<<" Element["<<el<<"] lnode["<<i<<"]: ["<<el*nds+i<<"] ElList["<<OutElementList[el*nds+i]<<"]: NODE "<<mrRemesh.NodalPreIds[OutElementList[el*nds+i]]<<std::endl;
		  //   for(ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin() ; i_node != mrModelPart.NodesEnd() ; ++i_node)
		  //     {
		  // 	if( i_node->Id() == (unsigned int)mrRemesh.NodalPreIds[OutElementList[el*nds+i]] )
		  // 	  std::cout<<" Node "<<i_node->Id()<<" is in Model PART !!!! "<<std::endl;
		  //     }
		  // }

		  Vertices.push_back(*(nodes_begin + OutElementList[el*nds+i]-1).base());
		  //Vertices.push_back(mrModelPart.pGetNode(mrRemesh.NodalPreIds[OutElementList[el*nds+i]]));
		  Vertices.back().Set(CONTACT);
		}

	      id += 1;

	      Condition::Pointer pContactCondition = rReferenceCondition.Create(id, Vertices, pProperties);


	      //search the model part condition associated to this contact element MASTER CONDITION
	      //assign the MASTER ELEMENT and the MASTER NODE

	      bool condition_found=false;
	      MesherUtilities MesherUtils;
	      Condition::Pointer pMasterCondition = MesherUtils.FindMasterCondition(pContactCondition,mrModelPart.Conditions(),condition_found);

	      if(!condition_found){

		std::cout<<" MASTER CONDITION NOT FOUND: Contact Element Release "<<std::endl;
		id -= 1;
	      }
	      else{

		// std::cout<<" contact master "<<std::endl;
		// pMasterCondition->GetValue(MASTER_ELEMENTS)[0].PrintInfo(std::cout);
		// pMasterCondition->GetValue(MASTER_ELEMENTS)[0].PrintData(std::cout);
		// pMasterCondition->GetValue(MASTER_ELEMENTS)[0].GetProperties().PrintData(std::cout);
		// std::cout<<std::endl;

		pContactCondition->SetValue(MASTER_CONDITION, pMasterCondition );
		pContactCondition->SetValue(MASTER_ELEMENTS, pMasterCondition->GetValue(MASTER_ELEMENTS) );
		pContactCondition->SetValue(MASTER_NODES, pMasterCondition->GetValue(MASTER_NODES) );

		if( pContactCondition->Is(SELECTED) ){ //two master nodes needed

		  Element::ElementType& rMasterElement  = pMasterCondition->GetValue(MASTER_ELEMENTS).back();
		  Geometry< Node<3> >&  rMasterGeometry = rMasterElement.GetGeometry();
		  Element::NodeType&    rMasterNode     = pContactCondition->GetValue(MASTER_NODES).back();
		  Geometry< Node<3> >&  rGeometry       = pContactCondition->GetGeometry();

		  std::vector<bool> edge_nodes(4);
		  std::fill(edge_nodes.begin(), edge_nodes.end(), false);

		  for(unsigned int i=0; i<rMasterGeometry.PointsNumber(); ++i)
		    {
		      for(unsigned int j=0; j<rGeometry.PointsNumber(); ++j)
			{
			  if(rGeometry[j].Id()==rMasterGeometry[i].Id()){
			    edge_nodes[i] = true;
			    break;
			  }
			}
		    }

		  for(unsigned int i=0; i<4; ++i)
		    {
		      if(!edge_nodes[i] && rMasterGeometry[i].Id() != rMasterNode.Id())
			pContactCondition->GetValue(MASTER_NODES).push_back( Node<3>::WeakPointer(rMasterGeometry(i)) );
		    }
		}


		pContactCondition->SetValue(NORMAL, pMasterCondition->GetValue(NORMAL) );
		pContactCondition->Set(CONTACT);

		//set ACTIVE if it is going to be considered in computation:
		//here one can check the geometrical gap to dismiss some contact elements here
		//it will be done later in the contact element calculation
		pContactCondition->Set(ACTIVE);

		//setting new elements
		mrModelPart.AddCondition(pContactCondition);

	      }
	    }
	}


      //Restore global ID's
      for(ModelPart::NodesContainerType::iterator in = mrModelPart.NodesBegin() ; in != mrModelPart.NodesEnd(); ++in)
	{
	  in->SetId( mrRemesh.NodalPreIds[ in->Id() ] );
	}

      if( mEchoLevel > 0 ){
        std::cout<<"   [END   contact Elements Generation ["<<id-previous_id<<"] ]"<<std::endl;

        std::cout<<"   Total Conditions AFTER: ["<<mrModelPart.Conditions().size()<<"] ];"<<std::endl;
      }


      std::cout<<"  [Contact Candidates:"<<id-previous_id<<"]"<<std::endl;


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
      return "GenerateNewContactConditionsMesherProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "GenerateNewContactConditionsMesherProcess";
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
    ModelPart& mrModelPart;

    MesherUtilities::MeshingParameters& mrRemesh;

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
    GenerateNewContactConditionsMesherProcess& operator=(GenerateNewContactConditionsMesherProcess const& rOther);

    /// Copy constructor.
    //GenerateNewContactConditionsMesherProcess(GenerateNewContactConditionsMesherProcess const& rOther);


    ///@}

  }; // Class GenerateNewContactConditionsMesherProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    GenerateNewContactConditionsMesherProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const GenerateNewContactConditionsMesherProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif // KRATOS_GENERATE_NEW_CONTACT_CONDITIONS_MESHER_PROCESS_H_INCLUDED  defined
