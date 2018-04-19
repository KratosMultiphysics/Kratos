// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____ 
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _ 
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
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
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
// Processes
#include "processes/find_nodal_h_process.h"
#include "custom_processes/metric_fast_init_process.h"
#include "custom_processes/metrics_error_process.h"
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
    
    #if !defined(REMESHING_UTILITIES)
    #define REMESHING_UTILITIES
        enum RemeshingUtilities {MMG = 0};
    #endif
    
///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/** @brief Custom convergence criteria for the mortar condition 
 */
template<class TSparseSpace, class TDenseSpace>
class ErrorMeshCriteria : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ErrorMeshCriteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

    typedef TSparseSpace                              SparseSpaceType;

    typedef typename BaseType::TDataType                    TDataType;

    typedef typename BaseType::DofsArrayType            DofsArrayType;

    typedef typename BaseType::TSystemMatrixType    TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType    TSystemVectorType;
    
    typedef ModelPart::ConditionsContainerType    ConditionsArrayType;
    
    typedef ModelPart::NodesContainerType              NodesArrayType;
    
    typedef std::size_t KeyType;

    ///@}
    ///@name Life Cycle
    ///@{
    
    /// Default constructors
    ErrorMeshCriteria(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        )
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >(),
          mThisModelPart(rThisModelPart),
          mDimension(rThisModelPart.GetProcessInfo()[DOMAIN_SIZE]),
          mThisParameters(ThisParameters),
          mFindNodalH(FindNodalHProcess(mThisModelPart))
    {
        Parameters DefaultParameters = Parameters(R"(
        {
            "error_mesh_tolerance" : 1.0e-3,
            "error_mesh_constant" : 1.0e-3,
            "remeshing_utility"   : "MMG",
            "remeshing_parameters": 
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
                "enforce_current"                     : true
            }
        })" );
        
        mThisParameters.ValidateAndAssignDefaults(DefaultParameters);
        
        mErrorTolerance = mThisParameters["error_mesh_tolerance"].GetDouble();
        mConstantError = mThisParameters["error_mesh_constant"].GetDouble();
        mRemeshingUtilities = ConvertRemeshUtil(mThisParameters["remeshing_utility"].GetString());
        
        #if !defined(INCLUDE_MMG)
            if (mRemeshingUtilities == MMG)
            {
                KRATOS_ERROR << "YOU CAN USE MMG LIBRARY. CHECK YOUR COMPILATION" << std::endl;
            }
        #endif
    }

    ///Copy constructor 
    ErrorMeshCriteria( ErrorMeshCriteria const& rOther )
      :BaseType(rOther)
      ,mThisModelPart(rOther.mThisModelPart)
      ,mDimension(rOther.mDimension)
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
     * Compute relative and absolute error.
     * @param rModelPart Reference to the ModelPart containing the contact problem.
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
        // We recompute the NODAL_H
        mFindNodalH.Execute();
                
        // We initialize the check
        bool ConvergedError = true;
        
        // We initialize the total error
        double TotalErrorPow2 = 0.0;
        double CurrentSolPow2 = 0.0;
        
        // Iterate in the nodes
        NodesArrayType& NodesArray = rModelPart.Nodes();
        int numNodes = NodesArray.end() - NodesArray.begin();
        
//         #pragma omp parallel for // FIXME: Parallel not working
        for(int i = 0; i < numNodes; i++) 
        {
            auto itNode = NodesArray.begin() + i;
            
            double MainDoFNodalError = 0.0;           
            double OtherDoFNodalError = 0.0;        
            
            Node<3>::DofsContainerType& NodeDofs = (itNode)->GetDofs();

            std::size_t DoFId;
//             TDataType DofValue;
            
            for (typename Node<3>::DofsContainerType::const_iterator itDof = NodeDofs.begin(); itDof != NodeDofs.end(); itDof++)
            {
                if (itDof->IsFree())
                {
                    DoFId = itDof->EquationId();
//                     DofValue = itDof->GetSolutionStepValue(0);
                    
                    KeyType CurrVar = itDof->GetVariable().Key();
                    if ((CurrVar == DISPLACEMENT_X) || (CurrVar == DISPLACEMENT_Y) || (CurrVar == DISPLACEMENT_Z) || (CurrVar == VELOCITY_X) || (CurrVar == VELOCITY_Y) || (CurrVar == VELOCITY_Z))
                    {
                        MainDoFNodalError += b[DoFId] * b[DoFId];
                    }
                    else
                    {
                        OtherDoFNodalError += b[DoFId] * b[DoFId];
                    }
                }
//                 else
//                 {
// //                     const double Reaction = itDof->GetSolutionStepReactionValue();
// // //                     #pragma omp atomic
// //                     CurrentSolPow2 += Reaction * Reaction;
//                     
//                     DoFId = itDof->EquationId();
//             
// //                     #pragma omp atomic
//                     CurrentSolPow2 += b[DoFId] * b[DoFId];
//                 }
            }
        
            const double NodalH = itNode->FastGetSolutionStepValue(NODAL_H);
            
            const double NodalError = NodalH * NodalH * std::sqrt(MainDoFNodalError) + NodalH * std::sqrt(OtherDoFNodalError);
            itNode->SetValue(NODAL_ERROR, NodalError);
            
//             #pragma omp atomic
            TotalErrorPow2 += (NodalError * NodalError);
        }
        
        // Setting the average nodal error
        rModelPart.GetProcessInfo()[AVERAGE_NODAL_ERROR] = mErrorTolerance * std::sqrt((CurrentSolPow2 + TotalErrorPow2)/numNodes);
        
//         // Debug
//         KRATOS_WATCH(CurrentSolPow2)
//         KRATOS_WATCH(TotalErrorPow2)
        
        // Final check
        const double MeshError = mConstantError * std::sqrt(TotalErrorPow2);
//         const double MeshError = std::sqrt(TotalErrorPow2/(CurrentSolPow2 + TotalErrorPow2));
        if (MeshError > mErrorTolerance)
        {
            ConvergedError = false;
        }
    
        if (ConvergedError == true)
        {
            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            {
                std::cout << "The error due to the mesh size: " << MeshError << " is under the tolerance prescribed: " << mErrorTolerance << ". No remeshing required" << std::endl;
            }
        }
        else
        {
            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            {
                std::cout << "The error due to the mesh size: " << MeshError << " is bigger than the tolerance prescribed: " << mErrorTolerance << ". Remeshing required" << std::endl;
                std::cout << "AVERAGE_NODAL_ERROR: " << rModelPart.GetProcessInfo()[AVERAGE_NODAL_ERROR] << std::endl;
            }
            
            // Computing metric
            if (mDimension == 2)
            {
                MetricFastInit<2> MetricInit = MetricFastInit<2>(mThisModelPart);
                MetricInit.Execute();
                ComputeErrorSolMetricProcess<2> ComputeMetric = ComputeErrorSolMetricProcess<2>(mThisModelPart, mThisParameters["error_strategy_parameters"]);
                ComputeMetric.Execute();
            }
            else
            {
                MetricFastInit<3> MetricInit = MetricFastInit<3>(mThisModelPart);
                MetricInit.Execute();
                ComputeErrorSolMetricProcess<3> ComputeMetric = ComputeErrorSolMetricProcess<3>(mThisModelPart, mThisParameters["error_strategy_parameters"]);
                ComputeMetric.Execute();
            }
            
            // Remeshing
            if (mRemeshingUtilities == MMG)
            {
                #ifdef INCLUDE_MMG
                    if (mDimension == 2)
                    {
                        MmgProcess<2> MmgRemesh = MmgProcess<2>(mThisModelPart, mThisParameters["remeshing_parameters"]); 
                        MmgRemesh.Execute();
                    }
                    else
                    {
                        MmgProcess<3> MmgRemesh = MmgProcess<3>(mThisModelPart, mThisParameters["remeshing_parameters"]); 
                        MmgRemesh.Execute();
                    }
                #else 
                    KRATOS_ERROR << "Please compile with MMG to use this utility" << std::endl;
                #endif
            }
            else
            {
                KRATOS_ERROR << "Not an alternative utility" << std::endl;
            }
            
            mFindNodalH.Execute();
        }
        
        return ConvergedError;
    }
    
    /**
     * This function initialize the convergence criteria
     * @param rModelPart: The model part of interest
     */ 
    
    void Initialize(ModelPart& rModelPart) override
    {
        BaseType::mConvergenceCriteriaIsInitialized = true;
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
     * This converts the remehing utility string to an enum
     * @param str: The string that you want to comvert in the equivalent enum
     * @return RemeshingUtilities: The equivalent enum (this requires less memmory than a std::string)
     */
        
    RemeshingUtilities ConvertRemeshUtil(const std::string& str)
    {
        if(str == "MMG") 
        {
            return MMG;
        }
        else
        {
            return MMG;
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

    ModelPart& mThisModelPart;              // The model part where the refinement is computed
    const unsigned int mDimension;          // The dimension of the problem
    Parameters mThisParameters;             // The parameters
    
    FindNodalHProcess mFindNodalH;          // The process to copmpute NODAL_H
    RemeshingUtilities mRemeshingUtilities; // The remeshing utilities to use
    
    double mErrorTolerance;                 // The error tolerance considered
    double mConstantError;                  // The constant considered in the remeshing process
    
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

#endif /* KRATOS_ERROR_MESH_CRITERIA_H  defined */

