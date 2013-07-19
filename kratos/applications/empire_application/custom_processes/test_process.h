//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Pavel Ryzhakov, Daniel Baumgärtner, Johannes Wolf $
//   Date:                $Date: 2013-07-15 17:20:00 $
//   Revision:            $Revision: 1.0 $
//
//  this process is a test process which reads a model part and returns the element IDs

#if !defined(KRATOS_TEST_PROCESS_INCLUDED )
#define  KRATOS_TEST_PROCESS_INCLUDED



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

class TestProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TestProcess
    KRATOS_CLASS_POINTER_DEFINITION(TestProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TestProcess(ModelPart& model_part)
        : mr_model_part(model_part)
    {
    }

    /// Destructor.
    virtual ~TestProcess()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    void TestMemberFunction()
    {
        KRATOS_TRY

        for(ModelPart::NodeIterator i_node = mr_model_part.NodesBegin() ;
            i_node != mr_model_part.NodesEnd() ;
            i_node++)
        {
            KRATOS_WATCH(i_node->Id());
        }

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
        return "TestProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "TestProcess";
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
//		TestProcess& operator=(TestProcess const& rOther);

    /// Copy constructor.
//		TestProcess(TestProcess const& rOther);


    ///@}

}; // Class TestProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/*
/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  TestProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const TestProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
*/
///@}


}  // namespace Kratos.

#endif // KRATOS_TEST_PROCESS_INCLUDED  defined 


