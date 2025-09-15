//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla, Eduard Gómez
//
//

#if !defined(KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_SOLVING_STRATEGY_RUNGE_KUTTA)
#define KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_SOLVING_STRATEGY_RUNGE_KUTTA

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "factories/factory.h"
#include "solving_strategies/strategies/explicit_solving_strategy_runge_kutta.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "compressible_navier_stokes_explicit_solving_strategy.h"

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

/** @brief Explicit solving strategy base class
 * @details This is the base class from which we will derive all the explicit strategies (FE, RK4, ...)
 */
template <class TSparseSpace, class TDenseSpace, class TButcherTableau>
class CompressibleNavierStokesExplicitSolvingStrategyRungeKutta
: public CompressibleNavierStokesExplicitSolvingStrategy<ExplicitSolvingStrategyRungeKutta<TSparseSpace, TDenseSpace, TButcherTableau>>
{
public:
    ///@name Type Definitions
    ///@{

    /// The base class definition
    typedef CompressibleNavierStokesExplicitSolvingStrategy<ExplicitSolvingStrategyRungeKutta<TSparseSpace, TDenseSpace, TButcherTableau>> BaseType;

    /// The explicit builder and solver definition
    typedef typename BaseType::ExplicitBuilderType ExplicitBuilderType;

    /// The local vector definition
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    /// Pointer definition of CompressibleNavierStokesExplicitSolvingStrategyRungeKutta
    KRATOS_CLASS_POINTER_DEFINITION(CompressibleNavierStokesExplicitSolvingStrategyRungeKutta);

    /// Local Flags
    KRATOS_DEFINE_LOCAL_FLAG(SHOCK_CAPTURING);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit CompressibleNavierStokesExplicitSolvingStrategyRungeKutta(
        ModelPart &rModelPart,
        Parameters ThisParameters)
        : BaseType(rModelPart)
    {
        KRATOS_TRY

        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);

        KRATOS_CATCH("")
    }

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param pExplicitBuilder The pointer to the explicit builder and solver
     * @param MoveMeshFlag The flag to set if the mesh is moved or not
     */
    explicit CompressibleNavierStokesExplicitSolvingStrategyRungeKutta(
        ModelPart &rModelPart,
        typename ExplicitBuilderType::Pointer pExplicitBuilder,
        bool MoveMeshFlag = false,
        int RebuildLevel = 0)
        : BaseType(rModelPart, pExplicitBuilder, MoveMeshFlag, RebuildLevel)
    {
    }

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param MoveMeshFlag The flag to set if the mesh is moved or not
     */
    explicit CompressibleNavierStokesExplicitSolvingStrategyRungeKutta(
        ModelPart &rModelPart,
        bool MoveMeshFlag = false,
        int RebuildLevel = 0)
        : BaseType(rModelPart, MoveMeshFlag, RebuildLevel)
    {
    }

    /** Copy constructor.
     */
    CompressibleNavierStokesExplicitSolvingStrategyRungeKutta(const CompressibleNavierStokesExplicitSolvingStrategyRungeKutta &Other) = delete;

    /** Destructor.
     */
    virtual ~CompressibleNavierStokesExplicitSolvingStrategyRungeKutta() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        KRATOS_TRY

        Parameters default_parameters {};
        default_parameters.AddString("explicit_solving_strategy", Name());

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;

        KRATOS_CATCH("")
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "compressible_navier_stokes_explicit_solving_strategy_runge_kutta_" + TButcherTableau::Name();
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "CompressibleNavierStokesExplicitSolvingStrategyRungeKutta";
    }

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

    void InitializeRungeKuttaIntermediateSubStep() override
    {
        BaseType::InitializeRungeKuttaIntermediateSubStep();
        this->InitializeEverySubstep();
    }

    void InitializeRungeKuttaLastSubStep() override
    {
        BaseType::InitializeRungeKuttaLastSubStep();
        this->InitializeEverySubstep();
    }

    /**
     * @brief Finalize the Runge-Kutta intermediate substep
     * In this method we calculate the linearised time derivatives after the intermediate substep
     */
    void FinalizeRungeKuttaIntermediateSubStep() override
    {
        BaseType::FinalizeRungeKuttaIntermediateSubStep();
        this->FinalizeEverySubstep();
    }

    void FinalizeRungeKuttaLastSubStep() override
    {
        BaseType::FinalizeRungeKuttaLastSubStep();
        this->FinalizeEverySubstep();
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
}; /* Class CompressibleNavierStokesExplicitSolvingStrategyRungeKutta */

///@}

///@name Type Definitions
///@{

///@}

template<class TSparseSpace, class TDenseSpace>
using CompressibleNavierStokesExplicitSolvingStrategyForwardEuler = CompressibleNavierStokesExplicitSolvingStrategyRungeKutta<TSparseSpace, TDenseSpace, ButcherTableauForwardEuler>;

template<class TSparseSpace, class TDenseSpace>
using CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4 = CompressibleNavierStokesExplicitSolvingStrategyRungeKutta<TSparseSpace, TDenseSpace, ButcherTableauRK4>;

template<class TSparseSpace, class TDenseSpace>
using CompressibleNavierStokesExplicitSolvingStrategyRungeKutta3TVD = CompressibleNavierStokesExplicitSolvingStrategyRungeKutta<TSparseSpace, TDenseSpace, ButcherTableauRK3TVD>;

} /* namespace Kratos.*/

#endif /* KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_SOLVING_STRATEGY_RUNGE_KUTTA  defined */
