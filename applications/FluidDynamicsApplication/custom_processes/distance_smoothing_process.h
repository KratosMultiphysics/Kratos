//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohammad R. Hashemi
//

#ifndef KRATOS_DISTANCE_SMOOTHING_H
#define KRATOS_DISTANCE_SMOOTHING_H

// System includes
#include <string>
#include <iostream>

// Project includes
#include "processes/process.h"

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

/// Utility for distance smoothing
/// Based on Tornberg, Anna-Karin, and Björn Engquist. "A finite element based level-set method
/// for multiphase flow applications." Computing and Visualization in Science 3, no. 1-2 (2000): 93-101.
/// The algorithm is improved by imposing a boundary condition and a correction step.
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) DistanceSmoothingProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DistanceSmoothingProcess
    KRATOS_CLASS_POINTER_DEFINITION(DistanceSmoothingProcess);

    typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;
    typedef typename SchemeType::Pointer SchemePointerType;
    typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverPointerType;
    typedef SolvingStrategy< TSparseSpace, TDenseSpace, TLinearSolver > SolvingStrategyType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    DistanceSmoothingProcess(
        ModelPart& rModelPart,
        typename TLinearSolver::Pointer);

    /// Constructor with Kratos parameters.
    DistanceSmoothingProcess(
        ModelPart& rModelPart,
        Parameters Parameters);

    /// Constructor with Kratos model
    DistanceSmoothingProcess(
        Model& rModel,
        Parameters Parameters);

    /// Destructor.
    ~DistanceSmoothingProcess() override
    {
        if(mrModel.HasModelPart( mAuxModelPartName ))
            mrModel.DeleteModelPart( mAuxModelPartName );
    }

    ///@}
    ///@name Operators
    ///@{

    void Execute() override;

    void Clear() override;

    ///@}
    ///@name Operations
    ///@{

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
        buffer << "Distance Smoothing Process" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "Distance Smoothing Process";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}

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

    ModelPart& mrModelPart;
    Model& mrModel;
    std::string mAuxModelPartName;

    typename SolvingStrategyType::UniquePointer mp_solving_strategy;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Create a Solution Strategy object
     * This method creates the linear solution strategy
     * @param pBuilderAndSolver Builder and solver pointer
     */
    void CreateSolutionStrategy(typename TLinearSolver::Pointer pLinearSolver);

    /**
     * @brief Initialize the process
     * Create a Model Part with the DistanceSmoothingElement
     */
    void CreateAuxModelPart();

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
    DistanceSmoothingProcess() = delete;

    /// Assignment operator.
    DistanceSmoothingProcess& operator=(DistanceSmoothingProcess const& rOther) = delete;

    /// Copy constructor.
    DistanceSmoothingProcess(DistanceSmoothingProcess const& rOther) = delete;

    ///@}

}; // Class DistanceSmoothingProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_DISTANCE_SMOOTHING__H
