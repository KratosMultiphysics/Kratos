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

//  System includes

//  External includes

//  Project includes
#include "solving_strategies/convergencecriterias/displacement_criteria.h"

namespace Kratos
{

///@addtogroup TrilinosApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class TrilinosDisplacementCriteria
 * @ingroup TrilinosApplication
 * @brief MPI version of the DisplacementCriteria.
 * @details This is a convergence criteria that considers the increment on the solution as criteria. The reactions from the RHS are not computed in the solution
 * @see DisplacementCriteria
 * @author Jordi Cotela
 */
template<class TSparseSpace,
         class TDenseSpace
         >
class TrilinosDisplacementCriteria
    : public DisplacementCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TrilinosDisplacementCriteria
    KRATOS_CLASS_POINTER_DEFINITION( TrilinosDisplacementCriteria );

    /// The definition of the base ConvergenceCriteria
    using BaseConvergenceCriteriaType = ConvergenceCriteria< TSparseSpace, TDenseSpace >;

    /// The definition of the base DisplacementCriteria
    using BaseType = DisplacementCriteria<TSparseSpace, TDenseSpace>;

    /// The definition of the class type
    using ClassType = TrilinosDisplacementCriteria<TSparseSpace, TDenseSpace>;

    /// The definition of the data type
    using TDataType = typename BaseType::TDataType;

    /// The definition of the DoF data type
    using DofType = typename Node::DofType;

    /// The definition of the dofs array type
    using DofsArrayType = typename BaseType::DofsArrayType;

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
    explicit TrilinosDisplacementCriteria()
        : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param ThisParameters The configuration parameters
     */
    explicit TrilinosDisplacementCriteria(Kratos::Parameters ThisParameters)
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
    explicit TrilinosDisplacementCriteria(
        TDataType NewRatioTolerance,
        TDataType AlwaysConvergedNorm)
        : BaseType(NewRatioTolerance, AlwaysConvergedNorm)
    {
    }

    /**
     * @brief Copy constructor
     * @param rOther The criteria to be copied
     */
    explicit TrilinosDisplacementCriteria( TrilinosDisplacementCriteria const& rOther )
      :BaseType(rOther)
    {
    }

    /**
     * @brief Destructor.
     */
    ~TrilinosDisplacementCriteria() override {}

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
    typename BaseConvergenceCriteriaType::Pointer Create(Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "trilinos_displacement_criteria";
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                            : "trilinos_displacement_criteria",
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
        return "TrilinosDisplacementCriteria";
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
        ) override
    {
        if constexpr (TSparseSpace::LinearAlgebraLibrary() == TrilinosLinearAlgebraLibrary::EPETRA) {
            BaseType::CalculateReferenceNorm(rDofSet, rModelPart);
        } else if constexpr (TSparseSpace::LinearAlgebraLibrary() == TrilinosLinearAlgebraLibrary::TPETRA) { // Kokkos conflicts with OpenMP in the block_for_each, so we need to specialize the code for TPETRA
        #if (HAVE_TPETRA)
            // Retrieve the data communicator
            const auto& r_data_communicator = rModelPart.GetCommunicator().GetDataCommunicator();

            // The current MPI rank
            const int rank = r_data_communicator.Rank();

            TDataType reference_disp_norm = TDataType();
            for (auto& rDof : rDofSet) {
                if (this->IsFreeAndLocalDof(rDof, rank)) {
                    const TDataType dof_value = rDof.GetSolutionStepValue();
                    reference_disp_norm += std::pow(dof_value, 2);
                }
            }
            BaseType::mReferenceDispNorm = std::sqrt(r_data_communicator.SumAll(reference_disp_norm));
        #else
            KRATOS_ERROR << "You must compile Kratos with TPETRA support" << std::endl;
        #endif
        } else {
            KRATOS_ERROR << "Only EPETRA and TPETRA are supported for now" << std::endl;
        }
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
        ) override
    {
        if constexpr (TSparseSpace::LinearAlgebraLibrary() == TrilinosLinearAlgebraLibrary::EPETRA) {
            return BaseType::CalculateFinalCorrectionNorm(rDofNum, rDofSet, rDx, rModelPart);
        } else if constexpr (TSparseSpace::LinearAlgebraLibrary() == TrilinosLinearAlgebraLibrary::TPETRA) { // Kokkos conflicts with OpenMP in the block_for_each, so we need to specialize the code for TPETRA
        #if (HAVE_TPETRA)
            // Retrieve the data communicator
            const auto& r_data_communicator = rModelPart.GetCommunicator().GetDataCommunicator();

            // The current MPI rank
            const int rank = r_data_communicator.Rank();

            // Initialize
            TDataType final_correction_norm = TDataType();
            unsigned int dof_num = 0;

            // Loop over Dofs
            for (auto& rDof : rDofSet) {
                if (this->IsFreeAndLocalDof(rDof, rank)) {
                    const TDataType variation_dof_value = TSparseSpace::GetValue(rDx, rDof.EquationId());
                    final_correction_norm += std::pow(variation_dof_value, 2);
                    ++dof_num;
                }
            }

            rDofNum = static_cast<SizeType>(r_data_communicator.SumAll(dof_num));
            return std::sqrt(r_data_communicator.SumAll(final_correction_norm));
        #else
            KRATOS_ERROR << "You must compile Kratos with TPETRA support" << std::endl;
        #endif
        } else {
            KRATOS_ERROR << "Only EPETRA and TPETRA are supported for now" << std::endl;
        }
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
}; //  Class ClassName

///@}
///@name Type Definitions
///@{

///@}

}  //  namespace Kratos.
