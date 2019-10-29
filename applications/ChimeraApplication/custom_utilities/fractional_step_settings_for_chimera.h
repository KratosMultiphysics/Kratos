//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Authors:        Aditya Ghantasala, https://github.com/adityaghantasala
// 					Navaneeth K Narayanan
//					Rishith Ellath Meethal
// ==============================================================================
//


#ifndef KRATOS_FRACTIONAL_STEP_SETTINGS_FOR_CHIMERA_H
#define KRATOS_FRACTIONAL_STEP_SETTINGS_FOR_CHIMERA_H

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"

#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "processes/process.h"
#include "processes/fast_transfer_between_model_parts_process.h"
#include "containers/model.h"

// Application includes
#include "custom_utilities/solver_settings.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "custom_strategies/custom_builder_and_solver/residualbased_block_builder_and_solver_with_constraints_for_chimera.h"


namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

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

/// Helper class to define solution strategies for FS_Strategy.
template< class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver
          >
class KRATOS_API(CHIMERA_APPLICATION) FractionalStepSettingsForChimera: public SolverSettings<TSparseSpace,TDenseSpace,TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FractionalStepSettingsForChimera
    KRATOS_CLASS_POINTER_DEFINITION(FractionalStepSettingsForChimera);

    typedef SolverSettings<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;
    typedef typename BaseType::StrategyType StrategyType;
    typedef typename BaseType::StrategyPointerType StrategyPointerType;
    typedef typename BaseType::ProcessPointerType ProcessPointerType;
    typedef typename BaseType::StrategyLabel StrategyLabel;
    typedef ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> ResidualBasedBlockBuilderAndSolverType;
    typedef typename ResidualBasedBlockBuilderAndSolverType::Pointer ResidualBasedBlockBuilderAndSolverPointerType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    FractionalStepSettingsForChimera(ModelPart& rModelPart,
                   const std::size_t ThisDomainSize,
                   const std::size_t ThisTimeOrder,
                   const bool UseSlip,
                   const bool MoveMeshFlag,
                   const bool ReformDofSet);

    /// Destructor.
    ~FractionalStepSettingsForChimera() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    void SetStrategy(StrategyLabel const& rStrategyLabel,
                             typename TLinearSolver::Pointer pLinearSolver,
                             const double Tolerance,
                             const unsigned int MaxIter) override ;


    bool FindStrategy(StrategyLabel const& rStrategyLabel,
                              StrategyPointerType& pThisStrategy) override;

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;


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

    /// Default constructor.
    FractionalStepSettingsForChimera() = delete;

    /// Assignment operator.
    FractionalStepSettingsForChimera& operator=(FractionalStepSettingsForChimera const& rOther) = delete;

    /// Copy constructor.
    FractionalStepSettingsForChimera(FractionalStepSettingsForChimera const& rOther) = delete;


    ///@}

}; // Class FractionalStepSettingsForChimera

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class TDenseSpace, class TSparseSpace, class TLinearSolver >
inline std::istream& operator >> (std::istream& rIStream,
                                  FractionalStepSettingsForChimera<TSparseSpace,TDenseSpace,TLinearSolver>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TDenseSpace, class TSparseSpace, class TLinearSolver >
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FractionalStepSettingsForChimera<TSparseSpace,TDenseSpace,TLinearSolver>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FRACTIONAL_STEP_SETTINGS_FOR_CHIMERA_H
