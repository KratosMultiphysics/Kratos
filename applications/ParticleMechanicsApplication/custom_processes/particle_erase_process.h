//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta
//


#if !defined(KRATOS_PARTICLE_ERASE_PROCESS_INCLUDED )
#define  KRATOS_PARTICLE_ERASE_PROCESS_INCLUDED

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
#include "includes/kratos_flags.h"
#include "utilities/math_utils.h"

#include "particle_mechanics_application.h"

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
//erases the nodes marked as
/** Detail class definition.
*/

class ParticleEraseProcess
        : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NodeAndElementEraseProcess
    KRATOS_CLASS_POINTER_DEFINITION(ParticleEraseProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ParticleEraseProcess(ModelPart& model_part)
        : mr_model_part(model_part)
    {
        KRATOS_TRY
                KRATOS_CATCH("");
    }

    /// Destructor.
    virtual ~ParticleEraseProcess()
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
        KRATOS_TRY;

        int initial_num_element = mr_model_part.NumberOfElements();
        mr_model_part.RemoveElements( TO_ERASE );
        int num_element = mr_model_part.NumberOfElements();

        KRATOS_INFO("ParticleEraseProcess") << "WARNING: " << num_element - initial_num_element << " particle elements has been erased.";

        KRATOS_CATCH("");
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
        return "ParticleEraseProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ParticleEraseProcess";
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
    PointerVector<Node<3> > mTrashedNodes;


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
    ParticleEraseProcess& operator=(ParticleEraseProcess const& rOther);

    /// Copy constructor.
    //NodeAndElementEraseProcess(NodeAndElementEraseProcess const& rOther);


    ///@}

}; // Class NodeAndElementEraseProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ParticleEraseProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ParticleEraseProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_NODE_AND_ELEMENT_ERASE_PROCESS_INCLUDED   defined

