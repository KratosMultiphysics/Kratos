// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Anna Rehr
//                   Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_CONTACT_ERROR_MESH_CRITERIA_H)
#define  KRATOS_CONTACT_ERROR_MESH_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/color_utilities.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

// Processes
#include "custom_processes/contact_spr_error_process.h"

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
 * @class ContactErrorMeshCriteria
 * @ingroup ContactStructuralMechanicsApplication
 * @brief Custom convergence for used to check the convergence in the mesh error
 * @tparam TSparseSpace The sparse space considered
 * @tparam TDenseSpace The dense space considered
 * @author Anna Rehr
 * @author Vicente Mataix Ferrandiz
 */
template<class TSparseSpace, class TDenseSpace>
class ContactErrorMeshCriteria
    : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ContactErrorMeshCriteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace >        BaseType;

    typedef TSparseSpace                                     SparseSpaceType;

    typedef typename BaseType::TDataType                           TDataType;

    typedef typename BaseType::DofsArrayType                   DofsArrayType;

    typedef typename BaseType::TSystemMatrixType           TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType           TSystemVectorType;

    typedef ModelPart::ConditionsContainerType           ConditionsArrayType;

    typedef ModelPart::NodesContainerType                     NodesArrayType;

    typedef std::size_t                                              KeyType;

    typedef std::size_t                                             SizeType;

    ///@}
    ///@name Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    explicit ContactErrorMeshCriteria(Parameters ThisParameters = Parameters(R"({})"))
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >(),
          mThisParameters(ThisParameters)
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
            "error_mesh_tolerance" : 5.0e-3,
            "error_mesh_constant"  : 5.0e-3,
            "compute_error_extra_parameters":
            {
                "penalty_normal"                      : 1.0e4,
                "penalty_tangential"                  : 1.0e4,
                "echo_level"                          : 0
            }
        })" );

        mThisParameters.ValidateAndAssignDefaults(default_parameters);

        mErrorTolerance = mThisParameters["error_mesh_tolerance"].GetDouble();
        mConstantError = mThisParameters["error_mesh_constant"].GetDouble();

    }

    ///Copy constructor
    ContactErrorMeshCriteria( ContactErrorMeshCriteria const& rOther )
      :BaseType(rOther)
      ,mErrorTolerance(rOther.mErrorTolerance)
      ,mConstantError(rOther.mConstantError)
    {
    }

    /// Destructor
    ~ContactErrorMeshCriteria() override = default;

    ///@}
    ///@name Operators
    ///@{

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
            SPRErrorProcess<2> compute_error_process = ContactSPRErrorProcess<2>(rModelPart, mThisParameters["compute_error_extra_parameters"]);
            compute_error_process.Execute();
        } else {
            SPRErrorProcess<3> compute_error_process = ContactSPRErrorProcess<3>(rModelPart, mThisParameters["compute_error_extra_parameters"]);
            compute_error_process.Execute();
        }

        // We get the estimated error
        const double estimated_error = process_info[ERROR_RATIO];

        // We check if converged
        const bool converged_error = (estimated_error > mErrorTolerance) ? false : true;

        if (converged_error) {
            KRATOS_INFO_IF("ContactErrorMeshCriteria", rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) << "NL ITERATION: " << process_info[NL_ITERATION_NUMBER] << "\tThe error due to the mesh size: " << estimated_error << " is under the tolerance prescribed: " << mErrorTolerance << ". " << BOLDFONT(FGRN("No remeshing required")) << std::endl;
        } else {
            KRATOS_INFO_IF("ContactErrorMeshCriteria", rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            << "NL ITERATION: " << process_info[NL_ITERATION_NUMBER] << "\tThe error due to the mesh size: " << estimated_error << " is bigger than the tolerance prescribed: " << mErrorTolerance << ". "<< BOLDFONT(FRED("Remeshing required")) << std::endl;
        }

        return converged_error;
    }

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Acces
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

}; // Class ContactErrorMeshCriteria

///@name Explicit Specializations
///@{

}  // namespace Kratos

#endif /* KRATOS_CONTACT_ERROR_MESH_CRITERIA_H  defined */

