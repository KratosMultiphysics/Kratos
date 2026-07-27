//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Riccardo Rossi
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#ifdef KRATOS_USE_FUTURE
#include "future/solving_strategies/strategies/implicit_strategy_data.h"
#endif

namespace Kratos::Future
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

/**
 * @class ConvergenceCriteria
 * @ingroup KratosCore
 * @brief This is the base class to define the different convergence criterion considered
 * @tparam TLinearAlgebra The linear algebra type
 * @author Ruben Zorrilla
 * @author Riccardo Rossi
*/
template<class TLinearAlgebra>
class ConvergenceCriteria
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ConvergenceCriteria
    KRATOS_CLASS_POINTER_DEFINITION(ConvergenceCriteria);

    /// The definition of the current class
    typedef ConvergenceCriteria< TLinearAlgebra > ClassType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    explicit ConvergenceCriteria()
    {
        SetEchoLevel(0);
    }

    /// Constructor with Parameters
    explicit ConvergenceCriteria(Kratos::Parameters ThisParameters)
    {
        // Validate and assign defaults
        ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /// Copy constructor
    explicit ConvergenceCriteria( ConvergenceCriteria const& rOther) = delete;

    /// Destructor
    virtual ~ConvergenceCriteria() = default;

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Create method
    /// @param ThisParameters The configuration parameters
    virtual typename ClassType::Pointer Create(Parameters ThisParameters) const
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
    }

    /**
     * @brief It sets the level of echo for the solving strategy
     * @param Level The level to set
     * @details The different levels of echo are:
     * - 0: Mute... no echo at all
     * - 1: Printing time and basic information
     * - 2: Printing linear solver data
     * - 3: Print of debug information: Echo of stiffness matrix, Dx, b...
     */
    void SetEchoLevel(const int Level)
    {
        mEchoLevel = Level;
    }

    /**
     * @brief Criterias that need to be called before getting the solution
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param rImplicitStrategyData Data container of the implicit strategy
     * @return true if convergence is achieved, false otherwise
     */
    //FIXME: I THINK WE SHOULD REMOVE THIS!
    virtual bool PreCriteria(
        ModelPart& rModelPart,
        ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData)
    {
        return true;
    }

    /**
     * @brief Criterias that need to be called after getting the solution
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param rImplicitStrategyData Data container of the implicit strategy
     * @return true if convergence is achieved, false otherwise
     */
    //FIXME: I THINK WE SHOULD REMOVE THIS AND JUST CALL IT IsConverged()!
    virtual bool PostCriteria(
        ModelPart& rModelPart,
        ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData)
    {
        return true;
    }

    /**
     * @brief Checks if the solution is converged
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param rImplicitStrategyData Data container of the implicit strategy
     * @return true if the solution is converged, false otherwise
     */
    virtual bool IsConverged(
        ModelPart& rModelPart,
        ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData) = 0;

    /**
     * @brief This function initialize the convergence criteria
     * @param rImplicitStrategyData Data container of the implicit strategy
     */
    virtual void Initialize(ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData)
    {
    }

    /**
     * @brief This function initializes the solution step
     * @warning Must be defined on the derived classes
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param rImplicitStrategyData Data container of the implicit strategy
     */
    virtual void InitializeSolutionStep(
        ModelPart& rModelPart,
        const ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData)
    {
    }

    /**
     * @brief This function initializes the non-linear iteration
     * @warning Must be defined on the derived classes
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param rImplicitStrategyData Data container of the implicit strategy
     */
    virtual void InitializeNonLinearIteration(
        ModelPart& rModelPart,
        const ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData)
    {
    }

    /**
     * @brief This function finalizes the non-linear iteration
     * @warning Must be defined on the derived classes
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param rImplicitStrategyData Data container of the implicit strategy
     */
    virtual void FinalizeNonLinearIteration(
        ModelPart& rModelPart,
        const ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData)
    {
    }

    /**
     * @brief This function finalizes the solution step
     * @warning Must be defined on the derived classes
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param rImplicitStrategyData Data container of the implicit strategy
     */
    virtual void FinalizeSolutionStep(
        ModelPart& rModelPart,
        const ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData)
    {
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed on the input provided. Checks can be "expensive" as the function is designed to catch user's errors.
     * @warning Must be defined on the derived classes
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @return 0 all OK, 1 otherwise
     */
    virtual int Check(ModelPart& rModelPart)
    {
        KRATOS_TRY

        return 0;
        KRATOS_CATCH("");
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    virtual Parameters GetDefaultParameters() const
    {
        const Parameters default_parameters = Parameters(R"({
            "name"       : "convergence_criteria",
            "echo_level" : 1
        })");
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "convergence_criteria";
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief This returns the level of echo for the solving strategy
     * @details The different levels of echo are:
     * - 0: Mute... no echo at all
     * - 1: Printing time and basic information
     * - 2: Printing linear solver data
     * - 3: Print of debug information: Echo of stiffness matrix, Dx, b...
     * @return Level of echo for the solving strategy
     */
    int GetEchoLevel() const
    {
        return mEchoLevel;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ConvergenceCriteria";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << Info();
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

    int mEchoLevel; /// The echo level

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    virtual void AssignSettings(const Parameters ThisParameters)
    {
        mEchoLevel = ThisParameters["echo_level"].GetInt();
    }


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

    ///@}

}; /// Class ConvergenceCriteria
} // namespace Kratos::Future.

