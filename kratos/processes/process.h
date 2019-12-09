//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

#if !defined(KRATOS_PROCESS_H_INCLUDED )
#define  KRATOS_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{

///@name Kratos Classes
///@{

// Some forward declarations to avoid increase a lot the compilation time
class Model;

/**
 * @class Process
 * @ingroup KratosCore
 * @brief The base class for all processes in Kratos.
 * @details The process is the base class for all processes and defines a simple interface for them.
    Execute method is used to execute the Process algorithms. While the parameters of this method
  can be very different from one Process to other there is no way to create enough overridden
  versions of it. For this reason this method takes no argument and all Process parameters must
  be passed at construction time. The reason is that each constructor can take different set of
  argument without any dependency to other processes or the base Process class.
  @author Pooyan Dadvand
  @author Riccardo Rossi
*/
class Process : public Flags
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(Process);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Process() : Flags() {}
    explicit Process(Flags options) : Flags( options ) {}

    /// Destructor.
    ~Process() override {}


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

    /**
     * @brief This method creates an pointer of the process
     * @details We consider as input a Mmodel and a set of Parameters for the sake of generality
     * @warning Must be overrided in each process implementation
     * @param rModel The model to be consider
     * @param ThisParameters The configuration parameters
     */
    virtual Process::Pointer Create(
        Model& rModel,
        Parameters ThisParameters
        )
    {
        KRATOS_ERROR << "Calling base class create. Please override this method in the corresonding Process" << std::endl;
        return nullptr;
    }

    /**
     * @brief Execute method is used to execute the Process algorithms.
     */
    virtual void Execute() {}

    /**
     * @brief This function is designed for being called at the beginning of the computations
     * right after reading the model and the groups
     */
    virtual void ExecuteInitialize()
    {
    }

    /**
     * @brief This function is designed for being execute once before the solution loop but after
     * all of the solvers where built
     */
    virtual void ExecuteBeforeSolutionLoop()
    {
    }


    /**
     * @brief This function will be executed at every time step BEFORE performing the solve phase
     */
    virtual void ExecuteInitializeSolutionStep()
    {
    }

    /**
     * @brief This function will be executed at every time step AFTER performing the solve phase
     */
    virtual void ExecuteFinalizeSolutionStep()
    {
    }


    /**
     * @brief This function will be executed at every time step BEFORE  writing the output
     */
    virtual void ExecuteBeforeOutputStep()
    {
    }


    /**
     * @brief This function will be executed at every time step AFTER writing the output
     */
    virtual void ExecuteAfterOutputStep()
    {
    }


    /**
     * @brief This function is designed for being called at the end of the computations
     */
    virtual void ExecuteFinalize()
    {
    }

    /**
     * @brief This function is designed for being called after ExecuteInitialize ONCE
     * to verify that the input is correct.
     */
    virtual int Check()
    {
        return 0;
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
        return "Process";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Process";
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
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    Process& operator=(Process const& rOther);

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
                                  Process& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Process& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_PROCESS_H_INCLUDED  defined


