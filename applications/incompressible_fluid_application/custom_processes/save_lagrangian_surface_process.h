//
//   Project Name:        Kratos
//   Last Modified by:    $Author: pavel $
//   Date:                $Date: 2008-10-02 15:38:16 $
//   Revision:            $Revision: 1.1 $
//
//  this process saves the boundary of a Lagrangian part for EMBEDDED technique

#if !defined(KRATOS_SAVE_LAGRANGIAN_SURFACE_PROCESS_INCLUDED )
#define  KRATOS_SAVE_LAGRANGIAN_SURFACE_PROCESS_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/condition.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"



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
	Update the PRESSURE_FORCE on the nodes


*/

class SaveLagrangianSurfaceProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PushStructureProcess
    KRATOS_CLASS_POINTER_DEFINITION(SaveLagrangianSurfaceProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SaveLagrangianSurfaceProcess()
    {
    }

    /// Destructor.
    virtual ~SaveLagrangianSurfaceProcess()
    {
    }


    ///@}
    ///@name Operators
    ///@{

//		void operator()()
//		{
//			SaveStructure();
//		}


    ///@}
    ///@name Operations
    ///@{

	//FLUID_MODEL_PART is the whole Lagrangian fluid domain, surface model_part is its boundary
    void  SaveSurfaceConditions(ModelPart& lagrangian_model_part, ModelPart& surface_model_part)
    {
        KRATOS_TRY
        surface_model_part.Elements().clear();
        surface_model_part.Conditions().clear();
        surface_model_part.Nodes().clear();
	
	KRATOS_WATCH(lagrangian_model_part.Conditions().size())
       
        KRATOS_WATCH("SAVING LAGRANGIAN SURFACE MADE BY CONDITIONS FOR EMBEDDED")
        for(ModelPart::ConditionsContainerType::iterator ic = lagrangian_model_part.ConditionsBegin() ;
                ic != lagrangian_model_part.ConditionsEnd() ; ++ic)
        {          
	    surface_model_part.Conditions().push_back(*(ic.base()));
        }

	for(ModelPart::NodesContainerType::iterator in = lagrangian_model_part.NodesBegin() ;
                in != lagrangian_model_part.NodesEnd() ; ++in)
        {          
	    if (in->FastGetSolutionStepValue(IS_BOUNDARY)==1.0 && in->GetValue(NEIGHBOUR_CONDITIONS).size()!=0)
		    surface_model_part.Nodes().push_back(*(in.base()));
        }
        KRATOS_CATCH("")
    }
	///////////////////////////////////////////////////////////////////////
    void  SaveSurfaceElements(ModelPart& lagrangian_model_part, ModelPart& structure_model_part, ModelPart& surface_model_part)
    {
        KRATOS_TRY
        surface_model_part.Elements().clear();
        surface_model_part.Conditions().clear();
        surface_model_part.Nodes().clear();
	
	surface_model_part.Elements()=structure_model_part.Elements();
	KRATOS_WATCH(surface_model_part)
       
        KRATOS_WATCH("SAVING LAGRANGIAN SURFACE MADE BY STRUCTURE ELEMENTS FOR EMBEDDED")
	//this one is for treating the membrane problems
       

	for(ModelPart::NodesContainerType::iterator in = lagrangian_model_part.NodesBegin() ;
                in != lagrangian_model_part.NodesEnd() ; ++in)
        {          
	    //DO IT IN A BETTER WAY AFTERWAEDS:: THERE MIGHT BE A MEMBRANE ELEMENT THAT TOUCHES THE LAGRANGIAN FLUID
	    if (in->FastGetSolutionStepValue(IS_STRUCTURE)==1.0 && in->FastGetSolutionStepValue(IS_FLUID)==0.0 && in->GetValue(NEIGHBOUR_ELEMENTS).size()!=0)
		    surface_model_part.Nodes().push_back(*(in.base()));
        }
	KRATOS_WATCH(surface_model_part)
        KRATOS_CATCH("")
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
        return "SaveLagrangianSurfaceProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "SaveLagrangianSurfaceProcess";
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
    //ModelPart& mr_fluid_model_part;
    //ModelPart& mr_structure_model_part;

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
//		SaveLagrangianSurfaceProcess& operator=(SaveLagrangianSurfaceProcess const& rOther);

    /// Copy constructor.
//		SaveLagrangianSurfaceProcess(SaveLagrangianSurfaceProcess const& rOther);


    ///@}

}; // Class SaveLagrangianSurfaceProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  SaveLagrangianSurfaceProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SaveLagrangianSurfaceProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_SAVE_LAGRANGIAN_SURFACE_PROCESS_INCLUDED  defined 


