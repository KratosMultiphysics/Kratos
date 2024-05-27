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
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
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
 * @class DisplacementCriteria
 * @ingroup KratosCore
 * @brief This is a convergence criteria that considers the increment on the solution as criteria
 * @details The reactions from the RHS are not computed in the solution
 * @author Riccardo Rossi
*/
template<class TSparseSpace,
         class TDenseSpace
         >
class DisplacementCriteria
    : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DisplacementCriteria
    KRATOS_CLASS_POINTER_DEFINITION( DisplacementCriteria );

    /// The definition of the base ConvergenceCriteria
    using BaseType = ConvergenceCriteria< TSparseSpace, TDenseSpace >;

    /// The definition of the current class
    using ClassType = DisplacementCriteria< TSparseSpace, TDenseSpace >;

    /// The definition of the sparse space type
    using SparseSpaceType = TSparseSpace;

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
    explicit DisplacementCriteria()
        : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param ThisParameters The configuration parameters
     */
    explicit DisplacementCriteria(Kratos::Parameters ThisParameters)
        : BaseType()
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /**
     * @brief Constructor 2 arguments
     * @param NewRatioTolerance The ratio tolerance for the convergence.
     * @param AlwaysConvergedNorm The absolute tolerance for the convergence.
     */
    explicit DisplacementCriteria(
        TDataType NewRatioTolerance,
        TDataType AlwaysConvergedNorm)
        : BaseType(),
          mRatioTolerance(NewRatioTolerance),
          mAlwaysConvergedNorm(AlwaysConvergedNorm)
    {
    }

    /**
     * @brief Copy constructor
     * @param rOther The criteria to be copied
     */
    explicit DisplacementCriteria( DisplacementCriteria const& rOther )
      :BaseType(rOther)
      ,mRatioTolerance(rOther.mRatioTolerance)
      ,mAlwaysConvergedNorm(rOther.mAlwaysConvergedNorm)
      ,mReferenceDispNorm(rOther.mReferenceDispNorm)
    {
    }

    /**
     * @brief Destructor.
     */
    ~DisplacementCriteria() override {}

    ///@}
    ///@name Operators
    ///@{

    /// Deleted assignment operator.
    DisplacementCriteria& operator=(DisplacementCriteria const& rOther) = delete;

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
     * @brief Compute relative and absolute error.
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
        if (SparseSpaceType::Size(rDx) != 0) { //if we are solving for something
            // Some values
            const int rank = rModelPart.GetCommunicator().GetDataCommunicator().Rank();
            const TDataType approx_zero_tolerance = std::numeric_limits<TDataType>::epsilon();

            // The final correction norm calculation
            SizeType size_solution = 0;
            const TDataType final_correction_norm = CalculateFinalCorrectionNorm(size_solution, rDofSet, rDx, rModelPart);

            TDataType ratio = 0.0;

            CalculateReferenceNorm(rDofSet, rModelPart);
            if (mReferenceDispNorm < approx_zero_tolerance) {
                KRATOS_WARNING("DisplacementCriteria") << "NaN norm is detected. Setting reference to convergence criteria" << std::endl;
                mReferenceDispNorm = final_correction_norm;
            }

            if(final_correction_norm < approx_zero_tolerance) {
                ratio = 0.0;
            } else {
                ratio = final_correction_norm/mReferenceDispNorm;
            }

            const TDataType float_size_solution = static_cast<TDataType>(size_solution);

            const TDataType absolute_norm = (final_correction_norm/std::sqrt(float_size_solution));

            KRATOS_INFO_IF("DISPLACEMENT CRITERION", this->GetEchoLevel() > 0 && rank == 0) << " :: [ Obtained ratio = " << ratio << "; Expected ratio = " << mRatioTolerance << "; Absolute norm = " << absolute_norm << "; Expected norm =  " << mAlwaysConvergedNorm << "]" << std::endl;

            rModelPart.GetProcessInfo()[CONVERGENCE_RATIO] = ratio;
            rModelPart.GetProcessInfo()[RESIDUAL_NORM] = absolute_norm;

            const bool has_achieved_convergence = ratio <= mRatioTolerance || absolute_norm < mAlwaysConvergedNorm; //  || (final_correction_norm/x.size())<=1e-7)
            KRATOS_INFO_IF("DISPLACEMENT CRITERION", has_achieved_convergence && this->GetEchoLevel() > 0 && rank == 0) << "Convergence is achieved" << std::endl;
            return has_achieved_convergence;
        } else { //in this case all the displacements are imposed!
            return true;
        }
    }

    /**
     * @brief This function initialize the convergence criteria
     * @param rModelPart Reference to the ModelPart containing the problem. (unused)
     */
    void Initialize(
        ModelPart& rModelPart
        ) override
    {
        BaseType::mConvergenceCriteriaIsInitialized = true;
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
    }

    /**
     * @brief This function finalizes the solution step
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual + reactions)
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
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "displacement_criteria";
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                            : "displacement_criteria",
            "displacement_relative_tolerance" : 1.0e-4,
            "displacement_absolute_tolerance" : 1.0e-9
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
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
        return "DisplacementCriteria";
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

    TDataType mRatioTolerance = 0.0;      /// The ratio threshold for the norm of the residual

    TDataType mAlwaysConvergedNorm = 0.0; /// The absolute value threshold for the norm of the residual

    TDataType mReferenceDispNorm = 0.0;   /// The norm at the beginning of the iterations

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
        BaseType::AssignSettings(ThisParameters);
        mAlwaysConvergedNorm = ThisParameters["displacement_absolute_tolerance"].GetDouble();
        mRatioTolerance = ThisParameters["displacement_relative_tolerance"].GetDouble();
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

    /**
     * @brief Check if a Degree of Freedom (Dof) is free.
     * @details This function checks if a given Degree of Freedom (Dof) is free.
     * The reason why PARTITION_INDEX is considered in distributed runs is to avoid adding twice (or even more times) the same value into the norm
     * @param rDof The Degree of Freedom to check.
     * @param Rank The rank of the Dof.
     * @return True if the Dof is free, false otherwise.
     * @todo We should doo as in the residual criteria, and consider the active DoFs (not just free), taking into account the MPC in addition to fixed DoFs
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
     * @brief This method computes the reference norm
     * @details It checks if the dof is fixed
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @todo We should doo as in the residual criteria, and consider the active DoFs (not just free), taking into account the MPC in addition to fixed DoFs
     */
    void CalculateReferenceNorm(
        DofsArrayType& rDofSet,
        ModelPart& rModelPart
        )
    {
        // Retrieve the data communicator
        const auto& r_data_communicator = rModelPart.GetCommunicator().GetDataCommunicator();

        // The current MPI rank
        const int rank = r_data_communicator.Rank();

        // Auxiliary struct
        struct TLS {TDataType dof_value{};};

        const TDataType reference_disp_norm = block_for_each<SumReduction<TDataType>>(rDofSet, TLS(), [this, &rank](auto& rDof, TLS& rTLS) {
            if (ClassType::IsFreeAndLocalDof(rDof, rank)) {
                rTLS.dof_value = rDof.GetSolutionStepValue();
                return std::pow(rTLS.dof_value, 2);
            } else {
                return TDataType();
            }
        });
        mReferenceDispNorm = std::sqrt(r_data_communicator.SumAll(reference_disp_norm));
    }

    /**
     * @brief This method computes the final norm
     * @details It checks if the dof is fixed
     * @param rDofNum The number of DoFs
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rDx Vector of results (variations on nodal variables)
     * @todo We should doo as in the residual criteria, and consider the active DoFs (not just free), taking into account the MPC in addition to fixed DoFs
     */
    TDataType CalculateFinalCorrectionNorm(
        SizeType& rDofNum,
        DofsArrayType& rDofSet,
        const TSystemVectorType& rDx,
        ModelPart& rModelPart
        )
    {
        // Retrieve the data communicator
        const auto& r_data_communicator = rModelPart.GetCommunicator().GetDataCommunicator();

        // The current MPI rank
        const int rank = r_data_communicator.Rank();

        // Initialize
        TDataType final_correction_norm = TDataType();
        unsigned int dof_num = 0;

        // Custom reduction
        using CustomReduction = CombinedReduction<SumReduction<TDataType>,SumReduction<unsigned int>>;

        // Auxiliary struct
        struct TLS {TDataType dof_value{}; TDataType variation_dof_value{};};

        // Loop over Dofs
        std::tie(final_correction_norm, dof_num) = block_for_each<CustomReduction>(rDofSet, TLS(), [this, &rDx, &rank](auto& rDof, TLS& rTLS) {
            if (this->IsFreeAndLocalDof(rDof, rank)) {
                rTLS.variation_dof_value = SparseSpaceType::GetValue(rDx, rDof.EquationId());
                return std::make_tuple(std::pow(rTLS.variation_dof_value, 2), 1);
            } else {
                return std::make_tuple(TDataType(), 0);
            }
        });

        rDofNum = static_cast<SizeType>(r_data_communicator.SumAll(dof_num));
        return std::sqrt(r_data_communicator.SumAll(final_correction_norm));
    }

    ///@}
}; /* Class DisplacementCriteria */

///@}
///@name Type Definitions
///@{

///@}

}  /* namespace Kratos.*/