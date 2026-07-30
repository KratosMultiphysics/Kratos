//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Riccardo Rossi
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#ifdef KRATOS_USE_FUTURE
#include "future/solving_strategies/convergence_criteria/convergence_criteria.h"
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
 * @class ResidualCriteria
 * @ingroup KratosCore
 * @brief This is a convergence criteria that considers the residual as criteria
 * @details The reactions from the RHS are not computed in the residual
 * @author Ruben Zorrilla
 * @author Riccardo Rossi
*/
template<class TLinearAlgebra>
class ResidualCriteria : public ConvergenceCriteria< TLinearAlgebra >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ResidualCriteria
    KRATOS_CLASS_POINTER_DEFINITION( ResidualCriteria );

    /// The definition of the base ConvergenceCriteria
    using BaseType = ConvergenceCriteria< TLinearAlgebra >;

    /// The definition of the current class
    using ClassType = ResidualCriteria< TLinearAlgebra >;

    /// The data type
    using DataType = typename TLinearAlgebra::DataType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    explicit ResidualCriteria()
        : BaseType()
    {
    }

    /// Constructor with Parameters
    explicit ResidualCriteria(Kratos::Parameters ThisParameters)
        : BaseType()
    {
        // Validate and assign defaults
        ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /// Copy constructor
    explicit ResidualCriteria( ResidualCriteria const& rOther ) = delete;

    /// Destructor
    ~ResidualCriteria() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method creates a new instance of the convergence criteria
     * @param ThisParameters The configuration parameters
     * @return A pointer to the new instance
     */
    typename BaseType::Pointer Create(Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
    }

    /**
     * @brief Checks if the solution is converged
     * @param rImplicitStrategyData Data container of the implicit strategy
     * @return true if the solution is converged, false otherwise
     */
    bool IsConverged(ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData) override
    {
        // Get effective solution residual vector
        // Note that this already accounts for MPC constraints
        const auto p_eff_dof_set = rImplicitStrategyData.pGetEffectiveDofSet();
        const auto p_eff_lin_sys = rImplicitStrategyData.pGetEffectiveLinearSystem();
        const auto p_eff_rhs = p_eff_lin_sys->pGetVector(Future::LinearSystemTags::DenseVectorTag::RHS);

        // Check if we are solving forsomething (i.e., not all DOFs are fixed)
        if (p_eff_rhs->size() == 0) {
            return true;
        }

        // Calculate the norm of the current residual (RHS)
        const auto [n_free_dofs, eff_rhs_norm] = this->CalculateVectorNorm(*p_eff_dof_set, *p_eff_rhs);

        // Calculate convergence ratio
        auto& r_model_part = this->GetModelPart();
        const DataType ratio = mInitialResidualNorm < std::numeric_limits<DataType>::epsilon() ? 0.0 : eff_rhs_norm / mInitialResidualNorm;
        r_model_part.GetProcessInfo()[CONVERGENCE_RATIO] = ratio;

        // Calculate absolute residual norm (RHS / sqrt(ndof))
        const DataType absolute_norm = (eff_rhs_norm / std::sqrt(static_cast<DataType>(n_free_dofs)));
        r_model_part.GetProcessInfo()[RESIDUAL_NORM] = absolute_norm;

        // Print current iteration information
        const int rank = r_model_part.GetCommunicator().GetDataCommunicator().Rank();
        KRATOS_INFO_IF("ResidualCriterion", this->GetEchoLevel() > 1 && rank == 0) << " :: [Initial residual norm = " << mInitialResidualNorm << "; Current residual norm =  " << eff_rhs_norm << "]" << std::endl;
        KRATOS_INFO_IF("ResidualCriterion", this->GetEchoLevel() > 0 && rank == 0) << " :: [Obtained ratio = " << ratio << "; Expected ratio = " << mRelativeTolerance << "; Absolute norm = " << absolute_norm << "; Expected norm =  " << mAbsoluteTolerance << "]" << std::endl;

        // Check convergence
        const bool is_converged = ratio <= mRelativeTolerance || absolute_norm < mAbsoluteTolerance;
        KRATOS_INFO_IF("ResidualCriterion", is_converged && this->GetEchoLevel() > 0 && rank == 0) << "Convergence achieved" << std::endl;

        return is_converged;
    }

    /**
     * @brief This function initializes the solution step
     * @param rImplicitStrategyData Data container of the implicit strategy
     */
    void InitializeSolutionStep(const ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData) override
    {
        BaseType::InitializeSolutionStep(rImplicitStrategyData);

        // Get effective solution residual vector
        // Note that this already accounts for MPC constraints
        const auto p_eff_dof_set = rImplicitStrategyData.pGetEffectiveDofSet();
        const auto p_eff_lin_sys = rImplicitStrategyData.pGetEffectiveLinearSystem();
        const auto p_eff_rhs = p_eff_lin_sys->pGetVector(Future::LinearSystemTags::DenseVectorTag::RHS);

        // Calculate the residual norm at the beginning of the solution step
        // Note that this will be used as reference for the convergence ratio
        const auto output = this->CalculateVectorNorm(*p_eff_dof_set, *p_eff_rhs);
        mInitialResidualNorm = std::get<1>(output);
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name" : "residual_criteria",
            "absolute_tolerance" : 1.0e-4,
            "relative_tolerance" : 1.0e-9
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "residual_criteria";
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
        return "ResidualCriteria";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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

    DataType mRelativeTolerance = DataType(); /// The ratio threshold for the norm of the residual

    DataType mAbsoluteTolerance = DataType(); /// The absolute value threshold for the norm of the residual

    DataType mInitialResidualNorm = DataType(); /// The reference norm of the residual computed at the beginning of the solution step

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
    void AssignSettings(const Parameters ThisParameters) override
    {
        // Assign base class settings (e.g., echo level)
        BaseType::AssignSettings(ThisParameters);

        // Assign tolerances
        mRelativeTolerance = ThisParameters["relative_tolerance"].GetDouble();
        mAbsoluteTolerance = ThisParameters["absolute_tolerance"].GetDouble();
    }

    ///@}

}; // Class ResidualCriteria

///@}

///@name Type Definitions
///@{

///@}

}  // namespace Kratos::Future
