// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_BASE_MORTAR_CRITERIA_H)
#define  KRATOS_BASE_MORTAR_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "contact_structural_mechanics_application_variables.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/process_info.h"
#include "utilities/mortar_utilities.h"
#include "utilities/variable_utils.h"
#include "custom_processes/aalm_adapt_penalty_value_process.h"
#include "custom_processes/compute_dynamic_factor_process.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

// DEBUG
#include "includes/gid_io.h"

namespace Kratos
{
///@addtogroup ContactStructuralMechanicsApplication
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
 * @class DisplacementLagrangeMultiplierContactCriteria 
 * @ingroup ContactStructuralMechanicsApplication 
 * @brief Custom convergence criteria for the mortar condition
 * @author Vicente Mataix Ferrandiz 
 */
template<class TSparseSpace, class TDenseSpace>
class BaseMortarConvergenceCriteria 
    : public  ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    /// Creating the corresponding pointer
    KRATOS_CLASS_POINTER_DEFINITION( BaseMortarConvergenceCriteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

    typedef TSparseSpace                              SparseSpaceType;

    typedef typename BaseType::TDataType                    TDataType;

    typedef typename BaseType::DofsArrayType            DofsArrayType;

    typedef typename BaseType::TSystemMatrixType    TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType    TSystemVectorType;
    
    typedef ModelPart::ConditionsContainerType    ConditionsArrayType;
    
    typedef ModelPart::NodesContainerType              NodesArrayType;
    
    typedef GidIO<> GidIOBaseType;

    ///@}
    ///@name Life Cycle
    ///@{
    
    /// Default constructors
    BaseMortarConvergenceCriteria(const bool IODebug = false)
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >(),
          mIODebug(IODebug),
          mpGidIO(nullptr)
    {
        if (mIODebug)
            mpGidIO = Kratos::make_shared<GidIOBaseType>("POST_LINEAR_ITER", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
    }

    ///Copy constructor 
    BaseMortarConvergenceCriteria( BaseMortarConvergenceCriteria const& rOther )
      :BaseType(rOther)
    {
    }

    /// Destructor
    ~BaseMortarConvergenceCriteria() override = default;

    ///@}
    ///@name Operators
    ///@{
    
    /**
     * @brief Criterias that need to be called before getting the solution
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */
    
    bool PreCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
        // The current process info
        ProcessInfo& process_info = rModelPart.GetProcessInfo();

        // We update the normals if necessary
        const auto normal_variation = process_info.Has(CONSIDER_NORMAL_VARIATION) ? static_cast<NormalDerivativesComputation>(process_info.GetValue(CONSIDER_NORMAL_VARIATION)) : NO_DERIVATIVES_COMPUTATION;
        if (normal_variation != NO_DERIVATIVES_COMPUTATION)
            MortarUtilities::ComputeNodesMeanNormalModelPart( rModelPart.GetSubModelPart("Contact") ); // Update normal of the conditions
        
        const bool adapt_penalty = process_info.Has(ADAPT_PENALTY) ? process_info.GetValue(ADAPT_PENALTY) : false;
        const bool dynamic_case = rModelPart.NodesBegin()->SolutionStepsDataHas(VELOCITY_X);
        
        /* Compute weighthed gap */
        if (adapt_penalty || dynamic_case) {
            // Set to zero the weighted gap
            ResetWeightedGap(rModelPart);
            
            ConditionsArrayType& conditions_array = rModelPart.GetSubModelPart("ComputingContact").Conditions();
        
            KRATOS_TRACE_IF("Empty model part", conditions_array.size() == 0) << "YOUR COMPUTING CONTACT MODEL PART IS EMPTY" << std::endl;
            
            #pragma omp parallel for
            for(int i = 0; i < static_cast<int>(conditions_array.size()); ++i)
                (conditions_array.begin() + i)->AddExplicitContribution(process_info);
        }
         
//         // In dynamic case
//         if ( dynamic_case ) {
//             ComputeDynamicFactorProcess compute_dynamic_factor_process = ComputeDynamicFactorProcess( rModelPart.GetSubModelPart("Contact") );
//             compute_dynamic_factor_process.Execute();
//         }
        
        // We recalculate the penalty parameter
        if ( adapt_penalty ) {
            AALMAdaptPenaltyValueProcess aalm_adaptation_of_penalty = AALMAdaptPenaltyValueProcess( rModelPart.GetSubModelPart("Contact") );
            aalm_adaptation_of_penalty.Execute();
        }
        
        return true;
    }
    
    /**
     * @brief Compute relative and absolute error.
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
        // We save the current WEIGHTED_GAP in the buffer 
        NodesArrayType& nodes_array = rModelPart.GetSubModelPart("Contact").Nodes();
        
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
            auto it_node = nodes_array.begin() + i;
            it_node->FastGetSolutionStepValue(WEIGHTED_GAP, 1) = it_node->FastGetSolutionStepValue(WEIGHTED_GAP);
        }
        
        // Set to zero the weighted gap
        ResetWeightedGap(rModelPart);
        
        ConditionsArrayType& conditions_array = rModelPart.GetSubModelPart("ComputingContact").Conditions();
        
        if (conditions_array.size() == 0) 
            KRATOS_TRACE("Empty model part") << "WARNING:: YOUR COMPUTING CONTACT MODEL PART IS EMPTY" << std::endl;
        
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(conditions_array.size()); ++i)
            (conditions_array.begin() + i)->AddExplicitContribution(rModelPart.GetProcessInfo());
        
        // GiD IO for debugging
        if (mIODebug == true) {            
            const int nl_iter = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
            const double label = static_cast<double>(nl_iter);
            
            if (nl_iter == 1) {
                mpGidIO->InitializeMesh(label);
                mpGidIO->WriteMesh(rModelPart.GetMesh());
                mpGidIO->FinalizeMesh();
                mpGidIO->InitializeResults(label, rModelPart.GetMesh());
            }
            
            mpGidIO->WriteNodalFlags(INTERFACE, "INTERFACE", rModelPart.Nodes(), label);
            mpGidIO->WriteNodalFlags(ACTIVE, "ACTIVE", rModelPart.Nodes(), label);
            mpGidIO->WriteNodalFlags(SLAVE, "SLAVE", rModelPart.Nodes(), label);
            mpGidIO->WriteNodalFlags(ISOLATED, "ISOLATED", rModelPart.Nodes(), label);
            mpGidIO->WriteNodalResults(NORMAL, rModelPart.Nodes(), label, 0);
            mpGidIO->WriteNodalResultsNonHistorical(DYNAMIC_FACTOR, rModelPart.Nodes(), label);
            mpGidIO->WriteNodalResultsNonHistorical(AUGMENTED_NORMAL_CONTACT_PRESSURE, rModelPart.Nodes(), label);
            mpGidIO->WriteNodalResults(DISPLACEMENT, rModelPart.Nodes(), label, 0);
            if (rModelPart.Nodes().begin()->SolutionStepsDataHas(VELOCITY_X) == true) {
                mpGidIO->WriteNodalResults(VELOCITY, rModelPart.Nodes(), label, 0);
                mpGidIO->WriteNodalResults(ACCELERATION, rModelPart.Nodes(), label, 0);
            }
            if (nodes_array.begin()->SolutionStepsDataHas(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE))
                mpGidIO->WriteNodalResults(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE, rModelPart.Nodes(), label, 0);
            else
                mpGidIO->WriteNodalResults(VECTOR_LAGRANGE_MULTIPLIER, rModelPart.Nodes(), label, 0);
            mpGidIO->WriteNodalResults(WEIGHTED_GAP, rModelPart.Nodes(), label, 0);
        }
        
        return true;
    }
    
    /**
     * @brief This function initialize the convergence criteria
     * @param rModelPart The model part of interest
     */ 
    
    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_ERROR << "YOUR ARE CALLING THE BASE MORTAR CRITERIA" << std::endl;
    }
    
    /**
     * @brief This function initializes the solution step
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
    
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    { 
        // Update normal of the conditions
        MortarUtilities::ComputeNodesMeanNormalModelPart( rModelPart.GetSubModelPart("Contact") );
        
        // GiD IO for debugging
        if (mIODebug == true) {
            mpGidIO->CloseResultFile();
            std::ostringstream new_name ;
            new_name << "POST_LINEAR_ITER_STEP=""POST_LINEAR_ITER_STEP=" << rModelPart.GetProcessInfo()[STEP];
            mpGidIO->ChangeOutputName(new_name.str());
        }
    }
    
    /**
     * @brief This function finalizes the solution step
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
    
    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    { 
        // GiD IO for debugging
        if (mIODebug == true)
            mpGidIO->FinalizeResults();
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
    
    /**
     * @brief This method resets the weighted gap in the nodes of the problem
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     */
    
    virtual void ResetWeightedGap(ModelPart& rModelPart)
    {       
        NodesArrayType& nodes_array = rModelPart.GetSubModelPart("Contact").Nodes();
        VariableUtils().SetScalarVar<Variable<double>>(WEIGHTED_GAP, 0.0, nodes_array);
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
    
    bool mIODebug;                  /// If we generate an output gid file in order to debug         
    GidIOBaseType::Pointer mpGidIO; /// The pointer to the debugging GidIO    
    
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
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}

}; // Class BaseMortarConvergenceCriteria 

///@name Explicit Specializations
///@{

}  // namespace Kratos 

#endif /* KRATOS_BASE_MORTAR_CRITERIA_H  defined */

