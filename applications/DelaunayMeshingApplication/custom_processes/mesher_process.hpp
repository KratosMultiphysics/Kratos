//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:               August 2016 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined(KRATOS_MESHER_PROCESS_H_INCLUDED)
#define  KRATOS_MESHER_PROCESS_H_INCLUDED



// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// The base class for processes passed to the solution scheme

/**
   They are used for the configuration of the solver
*/

class MesherProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MesherProcess
    KRATOS_CLASS_POINTER_DEFINITION(MesherProcess);

    ///@}
    ///@name Life Cycle
    ///@{ç

    /// Default constructor.
    MesherProcess() : Process() {}

    /// Constructor.
    MesherProcess(Flags options) : Process( options ) {}

    /// Destructor.
    virtual ~MesherProcess() {}


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


    /// Execute method is used to execute the MesherProcess algorithms.
    void Execute()  override
    {
      KRATOS_WARNING(" MesherProcess ") << " this method must be implemented and override " << std::endl;
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() final override
    {
      KRATOS_WARNING(" MesherProcess ") << " method not available " << std::endl;
    }

    /// this function is designed for being execute once before the solution loop but after all of the
    /// solvers where built
    void ExecuteBeforeSolutionLoop() final override
    {
      KRATOS_WARNING(" MesherProcess ") << " method not available " << std::endl;
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() final override
    {
      KRATOS_WARNING(" MesherProcess ") << " method not available " << std::endl;
    }


    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteFinalizeSolutionStep() final override
    {
      KRATOS_WARNING(" MesherProcess ") << " method not available " << std::endl;
    }


    /// this function will be executed at every time step BEFORE  writing the output
    void ExecuteBeforeOutputStep() final override
    {
      KRATOS_WARNING(" MesherProcess ") << " method not available " << std::endl;
    }


    /// this function will be executed at every time step AFTER writing the output
    void ExecuteAfterOutputStep() final override
    {
      KRATOS_WARNING(" MesherProcess ") << " method not available " << std::endl;
    }


    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    void ExecuteFinalize() final override
    {
      KRATOS_WARNING(" MesherProcess ") << " method not available " << std::endl;
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
        return "MesherProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MesherProcess";
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

    /// Copy constructor.
    MesherProcess(MesherProcess const& rOther);

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
    ///@}
    ///@name Private Operators
    ///@{
    ///@}
    ///@name Private Operations
    ///@{
    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    MesherProcess& operator=(MesherProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class MesherProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MesherProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MesherProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_MESHER_PROCESS_H_INCLUDED  defined
