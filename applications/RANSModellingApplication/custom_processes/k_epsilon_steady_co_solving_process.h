//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(K_EPSILON_STEADY_CO_SOLVING_PROCESS_H_INCLUDED)
#define K_EPSILON_STEADY_CO_SOLVING_PROCESS_H_INCLUDED

// System includes
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "processes/process.h"

// Application includes
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.h"
#include "custom_processes/k_epsilon_co_solving_process.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_variable_utils.h"
#include "rans_modelling_application_variables.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// The base class for all KEpsilonSteadyCoSolvingProcess in Kratos.
/** The KEpsilonSteadyCoSolvingProcess is the base class for all KEpsilonSteadyCoSolvingProcess and defines a simple interface for them.
    Execute method is used to execute the KEpsilonSteadyCoSolvingProcess algorithms. While the parameters of this method
  can be very different from one KEpsilonSteadyCoSolvingProcess to other there is no way to create enough overridden
  versions of it. For this reason this method takes no argument and all KEpsilonSteadyCoSolvingProcess parameters must
  be passed at construction time. The reason is that each constructor can take different set of
  argument without any dependency to other KEpsilonSteadyCoSolvingProcess or the base KEpsilonSteadyCoSolvingProcess class.
*/
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class KEpsilonSteadyCoSolvingProcess
    : public KEpsilonCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    typedef KEpsilonCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef ModelPart::NodeType NodeType;

    typedef ModelPart::NodesContainerType NodesContainerType;

    /// Pointer definition of KEpsilonSteadyCoSolvingProcess
    KRATOS_CLASS_POINTER_DEFINITION(KEpsilonSteadyCoSolvingProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    KEpsilonSteadyCoSolvingProcess(ModelPart& rModelPart, Parameters& rParameters, Process& rYPlusModelProcess)
        : BaseType(rModelPart, rParameters, rYPlusModelProcess)
    {
    }

    /// Destructor.
    ~KEpsilonSteadyCoSolvingProcess() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Execute method is used to execute the ScalarCoSolvingProcess algorithms.
    void Execute() override
    {
        if (this->mrModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] > 5)
            this->SolveSolutionStep();
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
        return "KEpsilonSteadyCoSolvingProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "KEpsilonSteadyCoSolvingProcess";
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
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KEpsilonSteadyCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>& operator=(
        KEpsilonSteadyCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver> const& rOther);

    /// Copy constructor.
    // KEpsilonSteadyCoSolvingProcess(KEpsilonSteadyCoSolvingProcess const& rOther);

    ///@}

}; // Class KEpsilonSteadyCoSolvingProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
inline std::istream& operator>>(std::istream& rIStream,
                                KEpsilonSteadyCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>& rThis);

/// output stream function
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
inline std::ostream& operator<<(std::ostream& rOStream,
                                const KEpsilonSteadyCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // K_EPSILON_STEADY_CO_SOLVING_PROCESS_H_INCLUDED defined
