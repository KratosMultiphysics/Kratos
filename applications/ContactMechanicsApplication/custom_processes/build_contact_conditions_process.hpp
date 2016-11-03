//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                August 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_BUILD_CONTACT_CONDITIONS_PROCESS_H_INCLUDED )
#define  KRATOS_BUILD_CONTACT_CONDITIONS_PROCESS_H_INCLUDED


// System includes
#include <string>
#include <iostream>

// External includes
#include <boost/timer.hpp>

// Project includes
#include "includes/model_part.h"
#include "custom_utilities/modeler_utilities.hpp"
#include "custom_processes/select_mesh_elements_process.hpp"

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
  class BuildContactConditionsProcess
    : public Process
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BuildContactConditionsProcess
    KRATOS_CLASS_POINTER_DEFINITION( BuildContactConditionsProcess );

    typedef ModelPart::ConditionType         ConditionType;
    typedef ModelPart::PropertiesType       PropertiesType;
    typedef ConditionType::GeometryType       GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BuildContactConditionsProcess(ModelPart& rModelPart,
			     ModelerUtilities::MeshingParameters& rRemeshingParameters,
			     int EchoLevel)
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
    { 
      mMeshId = mrRemesh.MeshId;
      mEchoLevel = EchoLevel;
    }

    /// Destructor.
    virtual ~BuildContactConditionsProcess()
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

    virtual void Execute()
    {
      KRATOS_TRY

      if( mEchoLevel > 0 )
	std::cout<<" [ GENERATE NEW CONTACT ELEMENTS: "<<std::endl;

      if( mrModelPart.Name() != mrRemesh.SubModelPartName )
	std::cout<<" ModelPart Supplied do not corresponds to the Meshing Domain: ("<<mrModelPart.Name()<<" != "<<mrRemesh.SubModelPartName<<")"<<std::endl;
      
      //*******************************************************************
      //selecting elements
      
	if( !mrRemesh.MeshElementsSelectedFlag )  //Select Mesh Elements not performed  ... is needed to be done before building new elements
	{  
	  std::cout<<" ERROR : no selection of elements performed before building the elements "<<std::endl;
	  SelectMeshElementsProcess SelectElements(mrModelPart,mrRemesh,mEchoLevel);
	  SelectElements.Execute();
	}


      //*******************************************************************
      //restore global ID's
      int id=0;
      for(ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin(mMeshId) ; i_node != mrModelPart.NodesEnd(mMeshId) ; i_node++)
	{
	  id= i_node->Id();
	  i_node->SetId( mrRemesh.NodalPreIds[ id ] );
	}


      //*******************************************************************
      //set consecutive ids for global conditions
      id=1;
      for(ModelPart::ConditionsContainerType::iterator i_cond = mrModelPart.GetParentModelPart()->ConditionsBegin(mMeshId) ; i_cond != mrModelPart.GetParentModelPart()->ConditionsEnd(mMeshId) ; i_cond++)
	{
	  i_cond->SetId(id);
	  id++;
	}

      //*******************************************************************
      //properties to be used in the generation
      Properties::Pointer pProperties = mrRemesh.GetProperties();

      Condition const & rReferenceCondition=mrRemesh.GetReferenceCondition(); //contact element
      
      const unsigned int nds = rReferenceCondition.GetGeometry().size();

      std::cout<<"   [START contact Element Generation "<<std::endl;

      int& OutNumberOfElements = mrRemesh.OutMesh.GetNumberOfElements();
      int* OutElementList      = mrRemesh.OutMesh.GetElementList();

      int previous_id = mrModelPart.GetParentModelPart()->Conditions().back().Id();
      id = previous_id;
      for(int el = 0; el<OutNumberOfElements; el++)
	{
	  if(mrRemesh.PreservedElements[el])
	    {
	      Geometry<Node<3> > Vertices;
	      for(unsigned int i=0; i<nds; i++)
		{
		  Vertices.push_back(mrModelPart.pGetNode(mrRemesh.NodalPreIds[OutElementList[el*nds+i]],mMeshId));
		  Vertices.back().Set(CONTACT);
		}
	      
	      id += 1;
	      
	      Condition::Pointer pContactCondition = rReferenceCondition.Create(id, Vertices, pProperties);
	      

	      //search the model part condition associated to this contact element MASTER CONDITION
	      //assign the MASTER ELEMENT and the MASTER NODE
				
	      bool condition_found=false;
	      ModelerUtilities ModelerUtils;
	      Condition::Pointer pMasterCondition = ModelerUtils.FindMasterCondition(pContactCondition,mrModelPart.Conditions(),condition_found);

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
		pContactCondition->SetValue(NORMAL, pMasterCondition->GetValue(NORMAL) );
		pContactCondition->Set(CONTACT);

		//setting new elements
		mrModelPart.AddCondition(pContactCondition);

	      }
	    }
	}
      
      std::cout<<"   [END   contact Elements Generation ["<<id-previous_id<<"] ]"<<std::endl;
      
      std::cout<<"   Total Conditions AFTER: ["<<mrModelPart.Conditions().size()<<"] ];"<<std::endl;


      if( mEchoLevel > 0 )
	std::cout<<"   GENERATE NEW CONTACT ELEMENTS ]; "<<std::endl;


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
    virtual std::string Info() const
    {
      return "BuildContactConditionsProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
      rOStream << "BuildContactConditionsProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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

    ModelerUtilities::MeshingParameters& mrRemesh;

    int mMeshId;
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
    BuildContactConditionsProcess& operator=(BuildContactConditionsProcess const& rOther);

    /// Copy constructor.
    //BuildContactConditionsProcess(BuildContactConditionsProcess const& rOther);


    ///@}

  }; // Class BuildContactConditionsProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    BuildContactConditionsProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const BuildContactConditionsProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif // KRATOS_BUILD_CONTACT_CONDITIONS_PROCESS_H_INCLUDED  defined 
