// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____ 
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _ 
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Anna Rehr
//

#if !defined(KRATOS_ERROR_MESH_CRITERIA_H)
#define  KRATOS_ERROR_MESH_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "meshing_application.h"
#include "includes/kratos_parameters.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "utilities/process_factory_utility.h"
// Processes
#include "processes/find_nodal_h_process.h"
#include "custom_processes/metric_fast_init_process.h"
#include "custom_processes/metrics_spr_error_process.h"
#ifdef INCLUDE_MMG
    #include "custom_processes/mmg_process.h"
#endif

namespace Kratos
{
///@addtogroup MeshingApplication
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
 * @ingroup MeshingApplication
 * @brief Custom convergence criteria for the mortar condition
 * @tparam TSparseSpace The sparse space considered
 * @tparam TDenseSpace The dense space considered
 * @author Vicente Mataix Ferrandiz
 * @author Anna Rehr
 */
template<class TSparseSpace, class TDenseSpace>
class ErrorMeshCriteria : public virtual  ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ErrorMeshCriteria );

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
    
    typedef ProcessFactoryUtility::Pointer                 ProcessesListType;

    ///@}
    ///@name Enum's
    ///@{

    enum class RemeshingUtilities {MMG = 0};

    ///@}
    ///@name Life Cycle
    ///@{
    
    /// Default constructors
    ErrorMeshCriteria(
        Parameters ThisParameters = Parameters(R"({})"),
        ProcessesListType pMyProcesses = nullptr
        )
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >(),
          mThisParameters(ThisParameters),
          mpMyProcesses(pMyProcesses)
    {
        Parameters default_parameters = Parameters(R"(
        {
            "error_mesh_tolerance" : 1.0e-3,
            "error_mesh_constant"  : 1.0e-3,
            "remeshing_utility"    : "MMG",
            "strategy"             : "Error",
            "remeshing_parameters" :
            {
                "filename"                             : "out",
                "framework"                            : "Lagrangian",
                "internal_variables_parameters"        :
                {
                    "allocation_size"                      : 1000, 
                    "bucket_size"                          : 4, 
                    "search_factor"                        : 2, 
                    "interpolation_type"                   : "LST",
                    "internal_variable_interpolation_list" :[]
                },
                "save_external_files"              : false,
                "max_number_of_searchs"            : 1000,
                "echo_level"                       : 3
            },
            "error_strategy_parameters": 
            {
                "minimal_size"                        : 0.1,
                "maximal_size"                        : 10.0, 
                "error"                               : 0.05
            }
        })" );
        
        mThisParameters.ValidateAndAssignDefaults(default_parameters);
        
        mErrorTolerance = mThisParameters["error_mesh_tolerance"].GetDouble();
        mConstantError = mThisParameters["error_mesh_constant"].GetDouble();
        mRemeshingUtilities = ConvertRemeshUtil(mThisParameters["remeshing_utility"].GetString());
        
#if !defined(INCLUDE_MMG)
        KRATOS_ERROR_IF(mRemeshingUtilities == RemeshingUtilities::MMG) << "YOU CANNOT USE MMG LIBRARY. CHECK YOUR COMPILATION" << std::endl;
#endif
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

    /**
     * @brief This function initialize the convergence criteria
     * @param rModelPart The model part of interest
     */
    void Initialize(ModelPart& rModelPart) override
    {
        BaseType::Initialize(rModelPart);

        // Initialize metrics
        if (rModelPart.GetProcessInfo()[DOMAIN_SIZE] == 2) {
            MetricFastInit<2> MetricInit = MetricFastInit<2>(rModelPart);
            MetricInit.Execute();
        } else {
            MetricFastInit<3> MetricInit = MetricFastInit<3>(rModelPart);
            MetricInit.Execute();
        }
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
        const bool converged_error = ComputeRemesh(rModelPart);
        
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

    /**
     * @brief This converts the remehing utility string to an enum
     * @param str The string that you want to comvert in the equivalent enum
     * @return RemeshingUtilities: The equivalent enum (this requires less memmory than a std::string)
     */
    RemeshingUtilities ConvertRemeshUtil(const std::string& str)
    {
        if(str == "MMG") {
            return RemeshingUtilities::MMG;
        } else {
            return RemeshingUtilities::MMG;
        }
    }
    
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

    Parameters mThisParameters;                /// The parameters
    
    RemeshingUtilities mRemeshingUtilities;    /// The remeshing utilities to use
    
    double mErrorTolerance;                    /// The error tolerance considered
    double mConstantError;                     /// The constant considered in the remeshing process
    
    ProcessesListType mpMyProcesses;           /// The processes list
    
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This process computes the remeshing depending of the error
     * @param rModelPart The model part of interest
     * @return True if converged, false otherwise
     */
    bool ComputeRemesh(ModelPart& rModelPart)
    {
        // Computing metric
        double estimated_error = 0;
        if (rModelPart.GetProcessInfo()[DOMAIN_SIZE] == 2) {
            SPRMetricProcess<2> ComputeMetric = SPRMetricProcess<2>(rModelPart, mThisParameters["error_strategy_parameters"]);
            ComputeMetric.Execute();
        } else {
            SPRMetricProcess<3> ComputeMetric = SPRMetricProcess<3>(rModelPart, mThisParameters["error_strategy_parameters"]);
            ComputeMetric.Execute();
        }

        // We get the estimated error
        estimated_error = rModelPart.GetProcessInfo()[ERROR_ESTIMATE];

        // We check if converged
        const bool converged_error = (estimated_error > mErrorTolerance) ? false : true;

        if (converged_error) {
            KRATOS_INFO_IF("ErrorMeshCriteria", rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) << "The error due to the mesh size: " << estimated_error << " is under the tolerance prescribed: " << mErrorTolerance << ". No remeshing required" << std::endl;
        } else {
            KRATOS_INFO_IF("ErrorMeshCriteria", rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            << "The error due to the mesh size: " << estimated_error << " is bigger than the tolerance prescribed: " << mErrorTolerance << ". Remeshing required" << std::endl
            << "AVERAGE_NODAL_ERROR: " << rModelPart.GetProcessInfo()[AVERAGE_NODAL_ERROR] << std::endl;

            // Remeshing
            if (mRemeshingUtilities == RemeshingUtilities::MMG) {
            #ifdef INCLUDE_MMG
                if (rModelPart.GetProcessInfo()[DOMAIN_SIZE] == 2) {
                    MmgProcess<2> MmgRemesh = MmgProcess<2>(rModelPart, mThisParameters["remeshing_parameters"]);
                    MmgRemesh.Execute();
                } else {
                    MmgProcess<3> MmgRemesh = MmgProcess<3>(rModelPart, mThisParameters["remeshing_parameters"]);
                    MmgRemesh.Execute();
                }
            #else
                KRATOS_ERROR << "Please compile with MMG to use this utility" << std::endl;
            #endif
            } else {
                KRATOS_ERROR << "Not an alternative utility" << std::endl;
            }

            // We set the model part as modified
            rModelPart.Set(MODIFIED, true);

            FindNodalHProcess find_nodal_h_process = FindNodalHProcess(rModelPart);
            find_nodal_h_process.Execute();

            // Processes initialization
            mpMyProcesses->ExecuteInitialize();
            // Processes before the loop
            mpMyProcesses->ExecuteBeforeSolutionLoop();
            // Processes of initialize the solution step
            mpMyProcesses->ExecuteInitializeSolutionStep();

            // We reset the model part as modified
            rModelPart.Set(MODIFIED, false);
        }

        return converged_error;
    }

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

#endif /* KRATOS_ERROR_MESH_CRITERIA_H  defined */

