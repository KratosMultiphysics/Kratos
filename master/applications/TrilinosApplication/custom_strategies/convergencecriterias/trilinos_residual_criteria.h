//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#pragma once

// System includes

// External includes

// Project includes
#include "solving_strategies/convergencecriterias/residual_criteria.h"

namespace Kratos
{
///@addtogroup TrilinosApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class TrilinosResidualCriteria
 * @ingroup TrilinosApplication
 * @brief MPI version of the ResidualCriteria.
 * @details Implements a convergence criteria based on the norm of the (free rows of) the RHS vector.
 * @see ResidualCriteria
 * @author Jordi Cotela
 */
template<class TSparseSpace,
         class TDenseSpace
         >
class TrilinosResidualCriteria
    : public ResidualCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TrilinosResidualCriteria
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosResidualCriteria);

    /// The definition of the base ConvergenceCriteria
    using BaseConvergenceCriteriaType = ConvergenceCriteria< TSparseSpace, TDenseSpace >;

    /// The definition of the base ResidualCriteria
    using BaseType = ResidualCriteria< TSparseSpace, TDenseSpace >;

    /// The definition of the current class
    using ClassType = TrilinosResidualCriteria< TSparseSpace, TDenseSpace >;

    /// The data type
    using TDataType = typename BaseType::TDataType;

    /// The dofs array type
    using DofsArrayType = typename BaseType::DofsArrayType;

    /// The definition of the DoF data type
    using DofType = typename Node::DofType;

    /// The sparse matrix type
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;

    /// The dense vector type
    using TSystemVectorType = typename BaseType::TSystemVectorType;

    /// Definition of the IndexType
    using IndexType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    explicit TrilinosResidualCriteria()
        : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param ThisParameters The configuration parameters
     */
    explicit TrilinosResidualCriteria(Kratos::Parameters ThisParameters)
        : BaseType(ThisParameters)
    {
    }

    /**
     * @brief Constructor 2 arguments
     * @param NewRatioTolerance The ratio tolerance for the convergence.
     * @param AlwaysConvergedNorm The absolute tolerance for the convergence.
     */
    explicit TrilinosResidualCriteria(TDataType NewRatioTolerance,TDataType AlwaysConvergedNorm):
        BaseType(NewRatioTolerance, AlwaysConvergedNorm)
    {
    }

    /**
     * @brief Copy constructor
     * @param rOther The criteria to be copied
     */
    explicit TrilinosResidualCriteria(const TrilinosResidualCriteria& rOther):
        BaseType(rOther)
    {
    }

    /// Destructor.
    ~TrilinosResidualCriteria() override {}

    ///@}
    ///@name Operators
    ///@{

    /// Deleted assignment operator.
    TrilinosResidualCriteria& operator=(TrilinosResidualCriteria const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create method
     * @param ThisParameters The configuration parameters
     */
    typename BaseConvergenceCriteriaType::Pointer Create(Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
    }

    /**
     * @brief This function initializes the solution step
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual + reactions)
     */
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {
        // We set the initial Id
        const auto& r_data_communicator = rModelPart.GetCommunicator().GetDataCommunicator();
        const int rank = r_data_communicator.Rank();
        for (auto& r_dof : rDofSet) {
            if (r_dof.GetSolutionStepValue(PARTITION_INDEX) == rank) {
                BaseType::mInitialDoFId = r_dof.EquationId();
                break;
            }
        }

        // Calling base class
        BaseType::InitializeSolutionStep(rModelPart, rDofSet, rA, rDx, rb);
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                        : "trilinos_residual_criteria",
            "residual_absolute_tolerance" : 1.0e-4,
            "residual_relative_tolerance" : 1.0e-9
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
        return "trilinos_residual_criteria";
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
        return "TrilionosResidualCriteria";
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
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method computes the active dofs
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param rDofSet The whole set of dofs
     */
    void ComputeActiveDofs(
        ModelPart& rModelPart,
        const ModelPart::DofsArrayType& rDofSet
        ) override
    {
        // Filling mActiveDofs when MPC exist
        ConstraintUtilities::DistributedComputeActiveDofs(rModelPart, BaseType::mActiveDofs, rDofSet, BaseType::mInitialDoFId);
    }

    ///@}
}; // Class TrilinosResidualCriteria

///@}

///@} addtogroup block

}  // namespace Kratos.
