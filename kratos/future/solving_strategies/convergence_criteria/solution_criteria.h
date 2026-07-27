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
 * @class SolutionCriteria
 * @ingroup KratosCore
 * @brief This is a convergence criteria that considers the increment on the solution as criteria
 * @details The reactions from the RHS are not computed in the solution
 * @author Riccardo Rossi
*/
template<class TLinearAlgebra>
class SolutionCriteria : public ConvergenceCriteria< TLinearAlgebra >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SolutionCriteria
    KRATOS_CLASS_POINTER_DEFINITION( SolutionCriteria );

    /// The definition of the base ConvergenceCriteria
    using BaseType = ConvergenceCriteria< TLinearAlgebra >;

    /// The definition of the current class
    using ClassType = SolutionCriteria< TLinearAlgebra >;

    /// The data type
    using DataType = typename TLinearAlgebra::DataType;

    /// The DOFs array type
    using DofsArrayType = typename ModelPart::DofsArrayType;

    /// The definition of the DOF type
    using DofType = typename ModelPart::DofType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    explicit SolutionCriteria()
        : BaseType()
    {
    }

    /// Constructor with Parameters
    explicit SolutionCriteria(Kratos::Parameters ThisParameters)
        : BaseType()
    {
        // Validate and assign defaults
        ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /// Copy constructor
    explicit SolutionCriteria( SolutionCriteria const& rOther ) = delete;

    /// Destructor
    ~SolutionCriteria() override = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create method
     * @param ThisParameters The configuration parameters
     */
    typename BaseType::Pointer Create(Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
    }

    /**
     * @brief Checks if the solution is converged
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param rImplicitStrategyData Data container of the implicit strategy
     * @return true if the solution is converged, false otherwise
     */
    bool IsConverged(
        ModelPart& rModelPart,
        ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData) override
    {
        // Get effective solution vector
        const auto p_eff_dof_set = rImplicitStrategyData.pGetEffectiveDofSet();
        const auto p_eff_lin_sys = rImplicitStrategyData.pGetEffectiveLinearSystem();
        const auto p_eff_dx = p_eff_lin_sys->pGetVector(Future::LinearSystemTags::DenseVectorTag::Dx);

        // Check if we are solving forsomething (i.e., not all DOFs are fixed)
        if (p_eff_dx->size() == 0) {
            return true;
        }

        // Calculate the norm of the solution increment (dx)
        const auto [n_free_dofs, dx_norm] = CalculateDxNorm(rModelPart, *p_eff_dof_set, *p_eff_dx);

        // Calculate the norm of the solution (reference norm)
        DataType sol_norm = CalculateReferenceNorm(rModelPart, *p_eff_dof_set);
        if (sol_norm < std::numeric_limits<DataType>::epsilon()) {
            KRATOS_WARNING("SolutionCriteria") << "Zero solution norm detected. Setting reference norm to dx norm" << std::endl;
            sol_norm = dx_norm;
        }

        // Calculate convergence ratio
        const DataType ratio = dx_norm < std::numeric_limits<DataType>::epsilon() ? 0.0 : dx_norm / sol_norm;
        rModelPart.GetProcessInfo()[CONVERGENCE_RATIO] = ratio;

        // Calculate absolute solution increment norm (dx / ndof)
        const DataType absolute_norm = (dx_norm / std::sqrt(static_cast<DataType>(n_free_dofs)));
        rModelPart.GetProcessInfo()[RESIDUAL_NORM] = absolute_norm;

        // Print current iteration information
        const int rank = rModelPart.GetCommunicator().GetDataCommunicator().Rank();
        KRATOS_INFO_IF("SolutionCriteria", this->GetEchoLevel() > 0 && rank == 0) << "[Obtained ratio = " << ratio << "; Expected ratio = " << mRelativeTolerance << "; Absolute norm = " << absolute_norm << "; Expected norm =  " << mAbsoluteTolerance << "]" << std::endl;

        // Check convergence
        const bool is_converged = ratio <= mRelativeTolerance || absolute_norm < mAbsoluteTolerance;
        KRATOS_INFO_IF("SolutionCriteria", is_converged && this->GetEchoLevel() > 0 && rank == 0) << "Convergence achieved" << std::endl;

        return is_converged;
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name" : "solution_criteria",
            "variable_name" : "DISPLACEMENT",
            "relative_tolerance" : 1.0e-4,
            "absolute_tolerance" : 1.0e-9
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
        return "solution_criteria";
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
        return "SolutionCriteria";
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

    const VariableData* mpSolutionVariable = nullptr; /// Pointer to the variable to be checked

    DataType mRelativeTolerance = 0.0; /// The ratio threshold for the norm of the solution increment (dx/x)

    DataType mAbsoluteTolerance = 0.0; /// The absolute value threshold for the norm of the solution (dx/ndof)

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

        // Check variable type (scalar and array) and set solution variable pointer
        const std::string variable_name = ThisParameters["variable_name"].GetString();
        if (KratosComponents<Variable<double>>::Has(variable_name)) {
            mpSolutionVariable = dynamic_cast<const VariableData*>(&KratosComponents<Variable<double>>::Get(variable_name));
        } else if (KratosComponents<Variable<array_1d<double,3>>>::Has(variable_name)) {
            mpSolutionVariable = dynamic_cast<const VariableData*>(&KratosComponents<Variable<array_1d<double, 3>>>::Get(variable_name));
        } else {
            KRATOS_ERROR << "Variable " << variable_name << " is not a valid solution variable. Note that only scalar and array variables are allowed." << std::endl;
        }

        // Assign tolerances
        mRelativeTolerance = ThisParameters["relative_tolerance"].GetDouble();
        mAbsoluteTolerance = ThisParameters["absolute_tolerance"].GetDouble();
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

    // /**
    //  * @brief Check if a Degree of Freedom (DOF) is free.
    //  * @details This function checks if a given Degree of Freedom (DOF) is free.
    //  * The reason why PARTITION_INDEX is considered in distributed runs is to avoid adding twice (or even more times) the same value into the norm
    //  * @param rDof The Degree of Freedom to check
    //  * @param Rank The rank of the DOF
    //  * @return True if the DOF is free, false otherwise
    //  */
    // bool IsFreeAndLocalDof(
    //     const DofType& rDof,
    //     const int Rank)
    // {
    //     if constexpr (!TSparseSpace::IsDistributed()) {
    //         return rDof.IsFree();
    //     } else {
    //         return (rDof.IsFree() && (rDof.GetSolutionStepValue(PARTITION_INDEX) == Rank));
    //     }
    // }

    /**
     * @brief This method computes the reference norm
     * @details It checks if the dof is fixed
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @todo We should doo as in the residual criteria, and consider the active DoFs (not just free), taking into account the MPC in addition to fixed DoFs
     */
    DataType CalculateReferenceNorm(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet)
    {
        // Retrieve the data communicator
        const auto& r_data_communicator = rModelPart.GetCommunicator().GetDataCommunicator();

        // Check if the problem is distributed
        DataType sol_norm;
        if (r_data_communicator.IsDistributed()) {

            // // The current MPI rank
            // const int rank = r_data_communicator.Rank();

            // // Auxiliary struct
            // struct TLS {TDataType dof_value{};};

            // const TDataType reference_disp_norm = block_for_each<SumReduction<TDataType>>(rDofSet, TLS(), [this, &rank](auto& rDof, TLS& rTLS) {
            //     if (ClassType::IsFreeAndLocalDof(rDof, rank)) {
            //         rTLS.dof_value = rDof.GetSolutionStepValue();
            //         return std::pow(rTLS.dof_value, 2);
            //     } else {
            //         return TDataType();
            //     }
            // });
            // mReferenceDispNorm = std::sqrt(r_data_communicator.SumAll(reference_disp_norm));
        } else {
            // Loop over Dofs and add the contribution of each free DOF to the norm
            sol_norm = block_for_each<SumReduction<DataType>>(rDofSet, [](auto& rDof) {
                if (rDof.IsFree()) {
                    return std::pow(rDof.GetSolutionStepValue(), 2);
                } else {
                    return DataType();
                }
            });
        }

        // Synchronize the norm among all processes and return
        return std::sqrt(r_data_communicator.SumAll(sol_norm));
    }

    /**
     * @brief This method computes the final solution increment norm
     * @param rDofSet Reference to the container of the problem's DOFs
     * @param rDx Vector of solution increment values (variations on nodal variables)
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @return A pair containing the number of free DOFs and the norm of the solution increment
     */
    std::pair<std::size_t, DataType> CalculateDxNorm(
        const ModelPart& rModelPart,
        const DofsArrayType& rDofSet,
        const TLinearAlgebra::VectorType& rDx)
    {
        // Retrieve the data communicator
        const auto& r_data_communicator = rModelPart.GetCommunicator().GetDataCommunicator();

        // Define custom reduction for parallel computation
        // First item in the reduction tuple: sum of the squared norm of the variation of the DOFs
        // Second item in the reduction tuple: number of free DOFs
        using CustomReductionType = CombinedReduction<SumReduction<DataType>,SumReduction<std::size_t>>;

        // Check if the problem is distributed
        DataType dx_norm;
        std::size_t n_free_dofs;
        if (r_data_communicator.IsDistributed()) {
            // // The current MPI rank
            // const int rank = r_data_communicator.Rank();

            // // Loop over Dofs and add the contribution of each free DOF to the norm
            // // Note that the PARTITION_INDEX is considered in distributed runs to avoid adding more than once the same value into the norm
            // std::tie(final_correction_norm, n_free_dofs) = block_for_each<CustomReductionType>(rDofSet, TLS(), [&rDx, rank](auto& rDof, TLS& rTLS) {
            //     if (rDof.IsFree() && (rDof.GetSolutionStepValue(PARTITION_INDEX) == rank)) {

            //         DataType dof_dx_value = rDx[]





            //         rTLS.variation_dof_value = SparseSpaceType::GetValue(rDx, rDof.EquationId());
            //         return std::make_tuple(std::pow(rTLS.variation_dof_value, 2), 1);
            //     } else {
            //         return std::make_tuple(DataType(), 0);
            //     }
            // });

        } else {

            // Loop over Dofs and add the contribution of each free DOF to the norm
            std::tie(dx_norm, n_free_dofs) = block_for_each<CustomReductionType>(rDofSet, [&rDx](auto& rDof) {
                if (rDof.IsFree()) {
                    return std::make_tuple(std::pow(rDx[rDof.EquationId()], 2), 1);
                } else {
                    return std::make_tuple(DataType(), 0);
                }
            });
        }

        // Communicator reduction
        n_free_dofs = r_data_communicator.SumAll(n_free_dofs);
        dx_norm = std::sqrt(r_data_communicator.SumAll(dx_norm));

        return std::make_pair(n_free_dofs, dx_norm);
    }

    ///@}
}; // class SolutionCriteria

///@}
///@name Type Definitions
///@{

///@}

}  /* namespace Kratos.*/