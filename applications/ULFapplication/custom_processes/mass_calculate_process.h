//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pavel Ryzhakov
//
// THIS PROCESS is INVENTED in order to CALCULATE THE MASS and
//STORE it node-wise.
//this is necessary for the projection step of ulf-frac method

#if !defined(KRATOS_MASS_CALCULATE_PROCESS_INCLUDED )
#define  KRATOS_MASS_CALCULATE_PROCESS_INCLUDED



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
#include "utilities/geometry_utilities.h"
#include "ULF_application.h"


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
	Update the MASS_FORCE on the nodes


*/

class MassCalculateProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MassCalculateProcess
    KRATOS_CLASS_POINTER_DEFINITION(MassCalculateProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MassCalculateProcess(ModelPart& model_part)
        : mr_model_part(model_part)
    {
    }

    /// Destructor. 
    ~MassCalculateProcess() override
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

        ProcessInfo& proc_info = mr_model_part.GetProcessInfo();
        double dummy;
        //first initialize the Mass force to the old value

        //set the Mass to the old value
        for(ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin() ;
                in != mr_model_part.NodesEnd() ; ++in)
        {
            in->FastGetSolutionStepValue(NODAL_MASS)=0.0;
        }

        for(ModelPart::ElementsContainerType::iterator im = mr_model_part.ElementsBegin() ;
                im != mr_model_part.ElementsEnd() ; ++im)
        {
            im->Calculate(NODAL_MASS,dummy,proc_info);
        }

        KRATOS_WATCH("Execute of Mass Calulate Process");



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
    std::string Info() const override
    {
        return "MassCalculateProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MassCalculateProcess";
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
//		MassCalculateProcess& operator=(MassCalculateProcess const& rOther);

    /// Copy constructor.
//		MassCalculateProcess(MassCalculateProcess const& rOther);


    ///@}

}; // Class MassCalculateProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MassCalculateProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MassCalculateProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_MASS_CALCULATE_PROCESS_INCLUDED  defined 


