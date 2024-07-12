//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#pragma once

// System includes

// External includes
#include <limits>

// Project includes
#include "includes/model_part.h"
#include "utilities/constraint_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

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

/**
 * @class ResidualCriteria
 * @ingroup KratosCore
 * @brief This is a convergence criteria that considers the residual as criteria
 * @details The reactions from the RHS are not computed in the residual
 * @author Riccardo Rossi
*/
template<class TSparseSpace,
         class TDenseSpace
         >
class ResidualCriteria
    : public  ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ResidualCriteria
    KRATOS_CLASS_POINTER_DEFINITION( ResidualCriteria );

    /// The definition of the base ConvergenceCriteria
    using BaseType = ConvergenceCriteria< TSparseSpace, TDenseSpace >;

    /// The definition of the current class
    using ClassType = ResidualCriteria< TSparseSpace, TDenseSpace >;

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

    /// Definition of the size type
    using SizeType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    explicit ResidualCriteria()
        : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param ThisParameters The configuration parameters
     */
    explicit ResidualCriteria(Kratos::Parameters ThisParameters)
        : BaseType()
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);

        this->mActualizeRHSIsNeeded = true;
    }

    /**
    * @brief Constructor 2 arguments
    * @param NewRatioTolerance The ratio tolerance for the convergence.
    * @param AlwaysConvergedNorm The absolute tolerance for the convergence.
    */
    explicit ResidualCriteria(
        TDataType NewRatioTolerance,
        TDataType AlwaysConvergedNorm)
        : BaseType(),
          mRatioTolerance(NewRatioTolerance),
          mAlwaysConvergedNorm(AlwaysConvergedNorm)
    {
        this->mActualizeRHSIsNeeded = true;
    }

    /**
    * @brief Copy constructor
    * @param rOther The criteria to be copied
    */
    explicit ResidualCriteria( ResidualCriteria const& rOther )
      :BaseType(rOther)
      ,mRatioTolerance(rOther.mRatioTolerance)
      ,mInitialResidualNorm(rOther.mInitialResidualNorm)
      ,mCurrentResidualNorm(rOther.mCurrentResidualNorm)
      ,mAlwaysConvergedNorm(rOther.mAlwaysConvergedNorm)
      ,mReferenceDispNorm(rOther.mReferenceDispNorm)
      ,mActiveDofs(rOther.mActiveDofs)
      ,mInitialDoFId(rOther.mInitialDoFId)
    {
        this->mActualizeRHSIsNeeded = true;
    }

    /**
     * @brief Destructor.
     */
    ~ResidualCriteria() override {}

    ///@}
    ///@name Operators
    ///@{

    /// Deleted assignment operator.
    ResidualCriteria& operator=(ResidualCriteria const& rOther) = delete;

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
     * @brief Criterion that need to be called after getting the solution
     * @details Compute relative and absolute error.
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual + reactions)
     * @return true if convergence is achieved, false otherwise
     */
    bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {
        if (TSparseSpace::Size(rb) != 0) { //if we are solving for something
            // Some values
            const int rank = rModelPart.GetCommunicator().GetDataCommunicator().Rank();

            // Calculate the residual norm
            SizeType size_residual;
            CalculateResidualNorm(rModelPart, mCurrentResidualNorm, size_residual, rDofSet, rb);

            TDataType ratio{};
            if(mInitialResidualNorm < std::numeric_limits<TDataType>::epsilon()) {
                ratio = 0.0;
            } else {
                ratio = mCurrentResidualNorm/mInitialResidualNorm;
            }

            const TDataType float_size_residual = static_cast<TDataType>(size_residual);
            const TDataType absolute_norm = (mCurrentResidualNorm/float_size_residual);

            KRATOS_INFO_IF("RESIDUAL CRITERION", this->GetEchoLevel() > 1 && rank == 0) << " :: [ Initial residual norm = " << mInitialResidualNorm << "; Current residual norm =  " << mCurrentResidualNorm << "]" << std::endl;
            KRATOS_INFO_IF("RESIDUAL CRITERION", this->GetEchoLevel() > 0 && rank == 0) << " :: [ Obtained ratio = " << ratio << "; Expected ratio = " << mRatioTolerance << "; Absolute norm = " << absolute_norm << "; Expected norm =  " << mAlwaysConvergedNorm << "]" << std::endl;

            rModelPart.GetProcessInfo()[CONVERGENCE_RATIO] = ratio;
            rModelPart.GetProcessInfo()[RESIDUAL_NORM] = absolute_norm;

            const bool has_achieved_convergence = ratio <= mRatioTolerance || absolute_norm < mAlwaysConvergedNorm;
            KRATOS_INFO_IF("RESIDUAL CRITERION", has_achieved_convergence && this->GetEchoLevel() > 0 && rank == 0) << "Convergence is achieved" << std::endl;
            return has_achieved_convergence;
        } else {
            return true;
        }
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
        BaseType::InitializeSolutionStep(rModelPart, rDofSet, rA, rDx, rb);

        // Filling mActiveDofs when MPC exist
        if (rModelPart.NumberOfMasterSlaveConstraints() > 0) {
            ComputeActiveDofs(rModelPart, rDofSet);
        }

        SizeType size_residual;
        CalculateResidualNorm(rModelPart, mInitialResidualNorm, size_residual, rDofSet, rb);
    }

    /**
     * @brief This function finalizes the solution step
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual + reactions)
     */
    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {
        BaseType::FinalizeSolutionStep(rModelPart, rDofSet, rA, rDx, rb);
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                        : "residual_criteria",
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

    TDataType mRatioTolerance{};                                     /// The ratio threshold for the norm of the residual

    TDataType mInitialResidualNorm{};                                /// The reference norm of the residual

    TDataType mCurrentResidualNorm{};                                /// The current norm of the residual

    TDataType mAlwaysConvergedNorm{};                                /// The absolute value threshold for the norm of the residual

    TDataType mReferenceDispNorm{};                                  /// The norm at the beginning of the iterations

    std::vector<int> mActiveDofs;                                    /// This vector contains the dofs that are active

    IndexType mInitialDoFId = std::numeric_limits<IndexType>::max(); /// The initial DoF Id

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method computes the active dofs
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param rDofSet The whole set of dofs
     */
    virtual void ComputeActiveDofs(
        ModelPart& rModelPart,
        const ModelPart::DofsArrayType& rDofSet
        )
    {
        // Filling mActiveDofs when MPC exist
        ConstraintUtilities::ComputeActiveDofs(rModelPart, mActiveDofs, rDofSet);
    }

    /**
     * @brief Check if a Degree of Freedom (Dof) is active
     * @details This function checks if a given Degree of Freedom (Dof) is active.
     * The reason why PARTITION_INDEX is considered in distributed runs is to avoid adding twice (or even more times) the same value into the norm
     * @param rDof The Degree of Freedom to check.
     * @param Rank The rank of the Dof.
     * @return True if the Dof is free, false otherwise.
     */
    bool IsActiveAndLocalDof(
        const DofType& rDof,
        const int Rank
        )
    {
        const IndexType dof_id = rDof.EquationId();
        if constexpr (!TSparseSpace::IsDistributedSpace()) {
            return mActiveDofs[dof_id] == 1;
        } else {
            KRATOS_DEBUG_ERROR_IF((dof_id - mInitialDoFId) >= mActiveDofs.size() && (rDof.GetSolutionStepValue(PARTITION_INDEX) == Rank)) << "DofId is greater than the size of the active Dofs vector. DofId: " << dof_id << "\tInitialDoFId: " << mInitialDoFId << "\tActiveDofs size: " << mActiveDofs.size() << std::endl;
            return (mActiveDofs[dof_id - mInitialDoFId] == 1 && (rDof.GetSolutionStepValue(PARTITION_INDEX) == Rank));
        }
    }

    /**
     * @brief Check if a Degree of Freedom (Dof) is free.
     * @details This function checks if a given Degree of Freedom (Dof) is free.
     * The reason why PARTITION_INDEX is considered in distributed runs is to avoid adding twice (or even more times) the same value into the norm
     * @param rDof The Degree of Freedom to check.
     * @param Rank The rank of the Dof.
     * @return True if the Dof is free, false otherwise.
     */
    bool IsFreeAndLocalDof(
        const DofType& rDof,
        const int Rank
        )
    {
        if constexpr (!TSparseSpace::IsDistributedSpace()) {
            return rDof.IsFree();
        } else {
            return (rDof.IsFree() && (rDof.GetSolutionStepValue(PARTITION_INDEX) == Rank));
        }
    }

    /**
     * @brief This method computes the norm of the residual
     * @details It checks if the dof is fixed
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rResidualSolutionNorm The norm of the residual
     * @param rDofNum The number of DoFs
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rb RHS vector (residual + reactions)
     */
    virtual void CalculateResidualNorm(
        ModelPart& rModelPart,
        TDataType& rResidualSolutionNorm,
        SizeType& rDofNum,
        DofsArrayType& rDofSet,
        const TSystemVectorType& rb
        )
    {
        // Retrieve the data communicator
        const auto& r_data_communicator = rModelPart.GetCommunicator().GetDataCommunicator();

        // The current MPI rank
        const int rank = r_data_communicator.Rank();

        // Initialize
        TDataType residual_solution_norm = TDataType();
        unsigned int dof_num = 0;

        // Custom reduction
        using CustomReduction = CombinedReduction<SumReduction<TDataType>,SumReduction<unsigned int>>;

        // Auxiliary struct
        struct TLS {TDataType residual_dof_value{};};

        // Loop over Dofs
        if (rModelPart.NumberOfMasterSlaveConstraints() > 0) {
            std::tie(residual_solution_norm, dof_num) = block_for_each<CustomReduction>(rDofSet, TLS(), [this, &rb, &rank](auto& rDof, TLS& rTLS) {
                if (this->IsActiveAndLocalDof(rDof, rank)) {
                    rTLS.residual_dof_value = TSparseSpace::GetValue(rb, rDof.EquationId());
                    return std::make_tuple(std::pow(rTLS.residual_dof_value, 2), 1);
                } else {
                    return std::make_tuple(TDataType(), 0);
                }
            });
        } else {
            std::tie(residual_solution_norm, dof_num) = block_for_each<CustomReduction>(rDofSet, TLS(), [this, &rb, &rank](auto& rDof, TLS& rTLS) {
                if (this->IsFreeAndLocalDof(rDof, rank)) {
                    rTLS.residual_dof_value = TSparseSpace::GetValue(rb, rDof.EquationId());
                    return std::make_tuple(std::pow(rTLS.residual_dof_value, 2), 1);
                } else {
                    return std::make_tuple(TDataType(), 0);
                }
            });
        }

        rDofNum = static_cast<SizeType>(r_data_communicator.SumAll(dof_num));
        rResidualSolutionNorm = std::sqrt(r_data_communicator.SumAll(residual_solution_norm));
    }

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);
        mAlwaysConvergedNorm = ThisParameters["residual_absolute_tolerance"].GetDouble();
        mRatioTolerance = ThisParameters["residual_relative_tolerance"].GetDouble();
    }

    ///@}

}; // Class ClassName

///@}

///@name Type Definitions
///@{

///@}

}  // namespace Kratos.
