// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Anna Rehr
//                   Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/color_utilities.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

// Processes
#include "custom_processes/spr_error_process.h"

namespace Kratos
{
///@addtogroup StructuralMechanicsApplication
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

/**
 * @class ErrorMeshCriteria
 * @ingroup StructuralMechanicsApplication
 * @brief Custom convergence for used to check the convergence in the mesh error
 * @tparam TSparseSpace The sparse space considered
 * @tparam TDenseSpace The dense space considered
 * @author Anna Rehr
 * @author Vicente Mataix Ferrandiz
 */
template<class TSparseSpace, class TDenseSpace>
class ErrorMeshCriteria
    : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ErrorMeshCriteria );

    /// The definition of the base ConvergenceCriteria
    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace >        BaseType;

    /// The definition of the current class
    typedef ErrorMeshCriteria< TSparseSpace, TDenseSpace >         ClassType;

    /// The data type
    typedef typename BaseType::TDataType TDataType;

    /// The dofs array type
    typedef typename BaseType::DofsArrayType DofsArrayType;

    /// The sparse matrix type
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    /// The dense vector type
    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    /// Definition of the IndexType
    typedef std::size_t IndexType;

    /// Definition of the size type
    typedef std::size_t SizeType;

    ///@}
    ///@name Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    explicit ErrorMeshCriteria(Parameters ThisParameters = Parameters(R"({})"))
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >()
    {
        // Validate and assign defaults
        mThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(mThisParameters);
    }

    ///Copy constructor
    ErrorMeshCriteria( ErrorMeshCriteria const& rOther )
      :BaseType(rOther)
      ,mErrorTolerance(rOther.mErrorTolerance)
      ,mConstantError(rOther.mConstantError)
    {
    }

    /// Destructor
    ~ErrorMeshCriteria() override = default;

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
     * @brief This function initialize the convergence criteria
     * @param rModelPart The model part of interest
     */
    void Initialize(ModelPart& rModelPart) override
    {
        BaseType::Initialize(rModelPart);
    }

    /**
     * @brief Compute relative and absolute error.
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */
    bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
        // The process info
        const ProcessInfo& process_info = rModelPart.GetProcessInfo();

        // Computing error
        if (process_info[DOMAIN_SIZE] == 2) {
            SPRErrorProcess<2> compute_error_process = SPRErrorProcess<2>(rModelPart, mThisParameters["compute_error_extra_parameters"]);
            compute_error_process.Execute();
        } else {
            SPRErrorProcess<3> compute_error_process = SPRErrorProcess<3>(rModelPart, mThisParameters["compute_error_extra_parameters"]);
            compute_error_process.Execute();
        }

        // We get the estimated error
        const double estimated_error = process_info[ERROR_RATIO];

        // We check if converged
        const bool converged_error = (estimated_error > mErrorTolerance) ? false : true;

        if (converged_error) {
            KRATOS_INFO_IF("ErrorMeshCriteria", rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) << "NL ITERATION: " << process_info[NL_ITERATION_NUMBER] << "\tThe error due to the mesh size: " << estimated_error << " is under the tolerance prescribed: " << mErrorTolerance << ". " << BOLDFONT(FGRN("No remeshing required")) << std::endl;
        } else {
            KRATOS_INFO_IF("ErrorMeshCriteria", rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            << "NL ITERATION: " << process_info[NL_ITERATION_NUMBER] << "\tThe error due to the mesh size: " << estimated_error << " is bigger than the tolerance prescribed: " << mErrorTolerance << ". "<< BOLDFONT(FRED("Remeshing required")) << std::endl;
        }

        return converged_error;
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        /**
         * error_mesh_tolerance: The tolerance in the convergence criteria of the error
         * error_mesh_constant: The constant considered in the remeshing process
         * penalty_normal: The penalty used in the normal direction (for the contact patch)
         * penalty_tangential: The penalty used in the tangent direction (for the contact patch)
         * echo_level: The verbosity
         */
        Parameters default_parameters = Parameters(R"(
        {
            "name"                 : "error_mesh_criteria",
            "error_mesh_tolerance" : 5.0e-3,
            "error_mesh_constant"  : 5.0e-3,
            "compute_error_extra_parameters":
            {
                "echo_level"                          : 0
            }
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
        return "error_mesh_criteria";
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

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
        mErrorTolerance = ThisParameters["error_mesh_tolerance"].GetDouble();
        mConstantError = ThisParameters["error_mesh_constant"].GetDouble();
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

    Parameters mThisParameters; /// The parameters

    double mErrorTolerance;     /// The error tolerance considered
    double mConstantError;      /// The constant considered in the remeshing process

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

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}

}; // Class ErrorMeshCriteria

///@name Explicit Specializations
///@{

}  // namespace Kratos

