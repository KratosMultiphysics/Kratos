//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Eduard Gómez
//
//

#if !defined(KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_SOLVING_STRATEGY_BFECC)
#define KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_SOLVING_STRATEGY_BFECC

/* System includes */
#if __cplusplus >= 201703L
#include <optional>
#endif

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "factories/factory.h"
#include "solving_strategies/strategies/explicit_solving_strategy_bfecc.h"

/* Appication includes */
#include "fluid_dynamics_application_variables.h"
#include "compressible_navier_stokes_explicit_solving_strategy.h"

namespace Kratos
{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

template <class TSparseSpace, class TDenseSpace>
class CompressibleNavierStokesExplicitSolvingStrategyBFECC
: public CompressibleNavierStokesExplicitSolvingStrategy<ExplicitSolvingStrategyBFECC<TSparseSpace, TDenseSpace>>
{
public:
    ///@name Type Definitions
    ///@{

    // The base solving strategy class definition
    typedef SolvingStrategy<TSparseSpace, TDenseSpace> SolvingStrategyType;

    // The base class definition
    typedef CompressibleNavierStokesExplicitSolvingStrategy<ExplicitSolvingStrategyBFECC<TSparseSpace, TDenseSpace>> BaseType;

    /// The definition of the current class
    typedef CompressibleNavierStokesExplicitSolvingStrategyBFECC<TSparseSpace, TDenseSpace> ClassType;

    // The explicit builder and solver definition
    typedef typename BaseType::ExplicitBuilderType ExplicitBuilderType;

    /// The DOF type
    typedef typename BaseType::DofType DofType;

    /// The local vector definition
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;
    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;

    // Replace with proper std::optional when upgrading to c++17
#if __cplusplus >= 201703L
    // using std::optional;
#else
    /** Naive implementation of std::optional, used as a replacement for c++11 support.
     * Use only with simple types.
     */
    template<typename T>
    struct optional {
        optional()           noexcept : mHasValue(false) { }
        optional(const T& V) noexcept : mHasValue(true), mValue(V) { }

        void reset() noexcept { mHasValue = false; }

        T& operator*() noexcept { return mValue; }
        T operator*() const noexcept { return mValue; }
        
        bool has_value() const noexcept { return mHasValue;}
    private:
        bool mHasValue;
        T mValue;
    };
#endif

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(CompressibleNavierStokesExplicitSolvingStrategyBFECC);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor. (empty)
     */
    explicit CompressibleNavierStokesExplicitSolvingStrategyBFECC()
        : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit CompressibleNavierStokesExplicitSolvingStrategyBFECC(
        ModelPart &rModelPart,
        Parameters ThisParameters)
        : BaseType(rModelPart)
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param pExplicitBuilder The pointer to the explicit builder and solver
     * @param MoveMeshFlag The flag to set if the mesh is moved or not
     */
    explicit CompressibleNavierStokesExplicitSolvingStrategyBFECC(
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
    explicit CompressibleNavierStokesExplicitSolvingStrategyBFECC(
        ModelPart &rModelPart,
        bool MoveMeshFlag = false,
        int RebuildLevel = 0)
        : BaseType(rModelPart, MoveMeshFlag, RebuildLevel)
    {
    }

    /**
     * @brief Create method
     * @param rModelPart The model part to be computed
     * @param ThisParameters The configuration parameters
     */
    typename SolvingStrategyType::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters
        ) const override
    {
        return Kratos::make_shared<ClassType>(rModelPart, ThisParameters);
    }

    /** Copy constructor.
     */
    CompressibleNavierStokesExplicitSolvingStrategyBFECC(const CompressibleNavierStokesExplicitSolvingStrategyBFECC &Other) = delete;

    /** Destructor.
     */
    ~CompressibleNavierStokesExplicitSolvingStrategyBFECC() override = default;

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        std::stringstream s;
        s << "compressible_navier_stokes_explicit_solving_strategy_bfecc";
        return s.str();
    }

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


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream ss;
        ss << "CompressibleNavierStokesExplicitSolvingStrategyBFECC";
        return ss.str();
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

    /**
     * @brief Initialize the BFECC forward substep
     * This method is intended to implement all the operations required before each BFECC initial forward substep
     */
    void InitializeBFECCForwardSubstep() override
    {
        BaseType::InitializeBFECCForwardSubstep();
        StashDiffusiveConstants();
        this->InitializeEverySubstep();
    };

    /**
     * @brief Initialize the BFECC backward substep
     * This method is intended to implement all the operations required before each BFECC backward substep
     */
    void InitializeBFECCBackwardSubstep() override
    {
        BaseType::InitializeBFECCBackwardSubstep();
        this->InitializeEverySubstep();
    };

    /**
     * @brief Initialize the BFECC final substep
     * This method is intended to implement all the operations required before each BFECC final substep
     */
    void InitializeBFECCFinalSubstep() override
    {
        BaseType::InitializeBFECCFinalSubstep();
        PopDiffusiveConstants();
        this->InitializeEverySubstep();
    };

    /**
     * @brief Initialize the BFECC forward substep
     * This method is intended to implement all the operations required before each BFECC initial forward substep
     */
    void FinalizeBFECCForwardSubstep() override
    {
        BaseType::FinalizeBFECCForwardSubstep();
        this->FinalizeEverySubstep();
    };

    /**
     * @brief Finalize the BFECC backward substep
     * This method is intended to implement all the operations required before each BFECC backward substep
     */
    void FinalizeBFECCBackwardSubstep() override
    {
        BaseType::FinalizeBFECCBackwardSubstep();
        this->FinalizeEverySubstep();
    };

    /**
     * @brief Finalize the BFECC final substep
     * This method is intended to implement all the operations required before each BFECC final substep
     */
    void FinalizeBFECCFinalSubstep() override
    {
        BaseType::FinalizeBFECCFinalSubstep();
        PopDiffusiveConstants();
        this->FinalizeEverySubstep();
    };

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

#if __cplusplus >= 201703L
    struct Stash {
        std::optional<double> conductivity = {};
        std::optional<double> dynamic_viscosity = {};
    } mDiffusionStash;
#else
    struct Stash {
        optional<double> conductivity = {};
        optional<double> dynamic_viscosity = {};
    } mDiffusionStash;
#endif

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void StashDiffusiveConstants()
    {
        KRATOS_TRY

        auto& model_part = BaseType::GetModelPart();
        
        KRATOS_ERROR_IF(model_part.NumberOfElements() == 0) << "Model part is devoid of elements!"; // This error must be removed when MPI is implemented

        auto& properties = model_part.ElementsBegin()->GetProperties();

        if(properties.Has(CONDUCTIVITY))
        {
            auto& r_conductivity = properties.GetValue(CONDUCTIVITY);
            mDiffusionStash.conductivity = r_conductivity;
            r_conductivity = 0;
        }

        if(properties.Has(DYNAMIC_VISCOSITY))
        {
            auto& r_dynamic_viscosity = properties.GetValue(DYNAMIC_VISCOSITY);
            mDiffusionStash.dynamic_viscosity = r_dynamic_viscosity;
            r_dynamic_viscosity = 0;
        }

        KRATOS_CATCH("")
    }

    void PopDiffusiveConstants()
    {
        KRATOS_TRY

        auto& model_part = BaseType::GetModelPart();
        
        KRATOS_ERROR_IF(model_part.NumberOfElements() == 0) << "Model part is devoid of elements!"; // This error must be removed when MPI is implemented

        auto& properties = model_part.ElementsBegin()->GetProperties();

        if(mDiffusionStash.conductivity.has_value())
        {
            properties.SetValue(CONDUCTIVITY, *mDiffusionStash.conductivity);
        }

        if(mDiffusionStash.dynamic_viscosity.has_value())
        {
            properties.SetValue(DYNAMIC_VISCOSITY, *mDiffusionStash.dynamic_viscosity);
        }

        mDiffusionStash.conductivity.reset();
        mDiffusionStash.dynamic_viscosity.reset();

        KRATOS_CATCH("")
    }

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
};

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_SOLVING_STRATEGY_BFECC  defined */
