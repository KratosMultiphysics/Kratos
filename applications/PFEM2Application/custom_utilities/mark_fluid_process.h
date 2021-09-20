//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Julio Marti.
//

//
//  this process puts a IS_IS_FLUID flag on the nodes of fluid elements

#if !defined(KRATOS_MARK_IS_FLUID_PROCESS_INCLUDED )
#define  KRATOS_MARK_IS_FLUID_PROCESS_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
//#include "custom_utilities/geometry_utilities2D.h"
//#include "custom_elements/updated_lagrangian_fluid.h"
//#include "custom_elements/updated_lagrangian_fluid3D.h"
//#include "custom_elements/updated_lagrangian_fluid_inc.h"
//#include "custom_elements/updated_lagrangian_fluid3D_inc.h"


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

class MarkFluidProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PushStructureProcess
    KRATOS_CLASS_POINTER_DEFINITION(MarkFluidProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MarkFluidProcess(ModelPart& model_part)
        : mr_model_part(model_part)
    {
    }

    /// Destructor.
    virtual ~MarkFluidProcess()
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

    virtual void Execute() override
    {
        KRATOS_TRY

        ProcessInfo& proc_info = mr_model_part.GetProcessInfo();
        double temp;

        for(ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin(); in!=mr_model_part.NodesEnd(); in++)
        {
            in->FastGetSolutionStepValue(IS_FLUID) = 0.0;
        }
        //set to 1 all the nodes surrounded by at least one fluid element
        for(ModelPart::ElementsContainerType::iterator ie = mr_model_part.ElementsBegin(); ie!=mr_model_part.ElementsEnd(); ie++)
        {
            ie->Calculate(IS_FLUID,temp,proc_info);
        }

        /*		for(ModelPart::NodesContainerType::const_iterator in = mr_model_part.NodesBegin(); in!=mr_model_part.NodesEnd(); in++)
        			{
        				if( (in->GetValue(NEIGHBOUR_ELEMENTS)).size() != 0)
        					{
        					in->FastGetSolutionStepValue(IS_IS_FLUID) = 1.0;

        					}
        				else
        							in->FastGetSolutionStepValue(IS_IS_FLUID) = 0.0;
        			}*/
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
    virtual std::string Info() const override
    {
        return "MarkFluidProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MarkFluidProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
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
    ModelPart& mr_model_part;

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
//		MarkFluidProcess& operator=(MarkFluidProcess const& rOther);

    /// Copy constructor.
//		MarkFluidProcess(MarkFluidProcess const& rOther);


    ///@}

}; // Class MarkFluidProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MarkFluidProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MarkFluidProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_MARK_IS_FLUID_PROCESS_INCLUDED  defined
