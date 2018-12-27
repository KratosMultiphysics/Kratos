/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 2.0

Copyright 2010
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

/*
 * File:   set_communicator_process.h
 * Author: jcotela
 *
 * Created on December 20, 2010, 7:43 PM
 */

#ifndef KRATOS_SET_MPI_COMMUNICATOR_PROCESS_H_INCLUDED
#define	KRATOS_SET_MPI_COMMUNICATOR_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <mpi.h>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "mpi/includes/mpi_communicator.h"
#include "processes/process.h"

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

/// Set up an MPI communicator for a model part
class SetMPICommunicatorProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SetMPICommunicatorProcess
    KRATOS_CLASS_POINTER_DEFINITION(SetMPICommunicatorProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SetMPICommunicatorProcess(ModelPart& rModelPart):
        mrModelPart(rModelPart)
    {}

    /// Destructor.
    virtual ~SetMPICommunicatorProcess() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
        KRATOS_WARNING("SetMPICommunicatorProcess")
            << "Calling deprecated process SetMPICommunicatorProcess.\n"
            << "Please use ModelPartCommunicatorUtilities::SetMPICommunicator(ModelPart) instead."
            << std::endl;

        VariablesList * mVariables_List = &mrModelPart.GetNodalSolutionStepVariablesList();
        mrModelPart.SetCommunicator(Communicator::Pointer(new MPICommunicator(mVariables_List)));
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
        std::stringstream buffer;
        buffer << "SetMPICommunicatorProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SetMPICommunicatorProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


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
    SetMPICommunicatorProcess& operator=(SetMPICommunicatorProcess const& rOther)
    {
        if(&rOther == this)
            return *this;
        else
        {
            *this = rOther;
            return *this;
        }
    }

//      /// Copy constructor.
//      SetMPICommunicatorProcess(SetMPICommunicatorProcess const& rOther) :
//        mrModelPart(rOther.mrModelPart)
//      {}


    ///@}

}; // Class SetMPICommunicatorProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  SetMPICommunicatorProcess& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SetMPICommunicatorProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif	/* KRATOS_SET_MPI_COMMUNICATOR_PROCESS_H_INCLUDED */

