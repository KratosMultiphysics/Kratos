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
#include "custom_utilities/contact_utilities.h"
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
 * @class BaseMortarConvergenceCriteria
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

    /// Pointer definition of BaseMortarConvergenceCriteria
    KRATOS_CLASS_POINTER_DEFINITION( BaseMortarConvergenceCriteria );

    /// Local Flags
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_DYNAMIC_FACTOR );
    KRATOS_DEFINE_LOCAL_FLAG( IO_DEBUG );
    KRATOS_DEFINE_LOCAL_FLAG( PURE_SLIP );

    /// The base class definition (and it subclasses)
    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;
    typedef typename BaseType::TDataType                    TDataType;
    typedef typename BaseType::DofsArrayType            DofsArrayType;
    typedef typename BaseType::TSystemMatrixType    TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType    TSystemVectorType;

    /// The sparse space used
    typedef TSparseSpace                              SparseSpaceType;

    /// The components containers
    typedef ModelPart::ConditionsContainerType    ConditionsArrayType;
    typedef ModelPart::NodesContainerType              NodesArrayType;

    typedef GidIO<> GidIOBaseType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    explicit BaseMortarConvergenceCriteria(
        const bool ComputeDynamicFactor = false,
        const bool IODebug = false,
        const bool PureSlip = false
        )
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >(),
          mpIO(nullptr)
    {
        // Set local flags
        mOptions.Set(BaseMortarConvergenceCriteria::COMPUTE_DYNAMIC_FACTOR, ComputeDynamicFactor);
        mOptions.Set(BaseMortarConvergenceCriteria::IO_DEBUG, IODebug);
        mOptions.Set(BaseMortarConvergenceCriteria::PURE_SLIP, PureSlip);

        if (mOptions.Is(BaseMortarConvergenceCriteria::IO_DEBUG)) {
            mpIO = Kratos::make_shared<GidIOBaseType>("POST_LINEAR_ITER", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
        }
    }

    ///Copy constructor
    BaseMortarConvergenceCriteria( BaseMortarConvergenceCriteria const& rOther )
      :BaseType(rOther),
       mOptions(rOther.mOptions),
       mpIO(rOther.mpIO)
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
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */
    bool PreCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {
        // The current process info
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

        // The contact model part
        ModelPart& r_contact_model_part = rModelPart.GetSubModelPart("Contact");

        // We update the normals if necessary
        const auto normal_variation = r_process_info.Has(CONSIDER_NORMAL_VARIATION) ? static_cast<NormalDerivativesComputation>(r_process_info.GetValue(CONSIDER_NORMAL_VARIATION)) : NO_DERIVATIVES_COMPUTATION;
        if (normal_variation != NO_DERIVATIVES_COMPUTATION) {
            ComputeNodesMeanNormalModelPartWithPairedNormal(rModelPart); // Update normal of the conditions
        }

        // Update tangent (must be updated even for constant normal)
        const bool frictional_problem = rModelPart.IsDefined(SLIP) ? rModelPart.Is(SLIP) : false;
        if (frictional_problem) {
            const bool has_lm = rModelPart.HasNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER);
            if (has_lm && mOptions.IsNot(BaseMortarConvergenceCriteria::PURE_SLIP)) {
                MortarUtilities::ComputeNodesTangentModelPart(r_contact_model_part);
            } else {
                MortarUtilities::ComputeNodesTangentModelPart(r_contact_model_part, &WEIGHTED_SLIP, 1.0, true);
            }
        }

        const bool adapt_penalty = r_process_info.Has(ADAPT_PENALTY) ? r_process_info.GetValue(ADAPT_PENALTY) : false;
        const bool dynamic_case = rModelPart.HasNodalSolutionStepVariable(VELOCITY);

        /* Compute weighthed gap */
        if (adapt_penalty || dynamic_case) {
            // Set to zero the weighted gap
            ResetWeightedGap(rModelPart);

            // Compute the contribution
            ContactUtilities::ComputeExplicitContributionConditions(rModelPart.GetSubModelPart("ComputingContact"));
        }

        // In dynamic case
        if ( dynamic_case && mOptions.Is(BaseMortarConvergenceCriteria::COMPUTE_DYNAMIC_FACTOR)) {
            ComputeDynamicFactorProcess compute_dynamic_factor_process( r_contact_model_part );
            compute_dynamic_factor_process.Execute();
        }

        // We recalculate the penalty parameter
        if ( adapt_penalty ) {
            AALMAdaptPenaltyValueProcess aalm_adaptation_of_penalty( r_contact_model_part );
            aalm_adaptation_of_penalty.Execute();
        }

        return true;
    }

    /**
     * @brief Compute relative and absolute error.
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual)
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
        // We save the current WEIGHTED_GAP in the buffer
        NodesArrayType& r_nodes_array = rModelPart.GetSubModelPart("Contact").Nodes();
        const auto it_node_begin = r_nodes_array.begin();

        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
            auto it_node = it_node_begin + i;
            it_node->FastGetSolutionStepValue(WEIGHTED_GAP, 1) = it_node->FastGetSolutionStepValue(WEIGHTED_GAP);
        }

        // Set to zero the weighted gap
        ResetWeightedGap(rModelPart);

        // Compute the contribution
        ContactUtilities::ComputeExplicitContributionConditions(rModelPart.GetSubModelPart("ComputingContact"));

        // GiD IO for debugging
        if (mOptions.Is(BaseMortarConvergenceCriteria::IO_DEBUG)) {
            const bool frictional_problem = rModelPart.IsDefined(SLIP) ? rModelPart.Is(SLIP) : false;
            const int nl_iter = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
            const double label = static_cast<double>(nl_iter);

            if (nl_iter == 1) {
                mpIO->InitializeMesh(label);
                mpIO->WriteMesh(rModelPart.GetMesh());
                mpIO->FinalizeMesh();
                mpIO->InitializeResults(label, rModelPart.GetMesh());
            }

            mpIO->WriteNodalFlags(INTERFACE, "INTERFACE", rModelPart.Nodes(), label);
            mpIO->WriteNodalFlags(ACTIVE, "ACTIVE", rModelPart.Nodes(), label);
            mpIO->WriteNodalFlags(SLAVE, "SLAVE", rModelPart.Nodes(), label);
            mpIO->WriteNodalFlags(ISOLATED, "ISOLATED", rModelPart.Nodes(), label);
            mpIO->WriteNodalResults(NORMAL, rModelPart.Nodes(), label, 0);
            mpIO->WriteNodalResultsNonHistorical(DYNAMIC_FACTOR, rModelPart.Nodes(), label);
            mpIO->WriteNodalResultsNonHistorical(AUGMENTED_NORMAL_CONTACT_PRESSURE, rModelPart.Nodes(), label);
            mpIO->WriteNodalResults(DISPLACEMENT, rModelPart.Nodes(), label, 0);
            if (rModelPart.Nodes().begin()->SolutionStepsDataHas(VELOCITY_X)) {
                mpIO->WriteNodalResults(VELOCITY, rModelPart.Nodes(), label, 0);
                mpIO->WriteNodalResults(ACCELERATION, rModelPart.Nodes(), label, 0);
            }
            if (r_nodes_array.begin()->SolutionStepsDataHas(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE))
                mpIO->WriteNodalResults(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE, rModelPart.Nodes(), label, 0);
            else if (r_nodes_array.begin()->SolutionStepsDataHas(VECTOR_LAGRANGE_MULTIPLIER_X))
                mpIO->WriteNodalResults(VECTOR_LAGRANGE_MULTIPLIER, rModelPart.Nodes(), label, 0);
            mpIO->WriteNodalResults(WEIGHTED_GAP, rModelPart.Nodes(), label, 0);
            if (frictional_problem) {
                mpIO->WriteNodalFlags(SLIP, "SLIP", rModelPart.Nodes(), label);
                mpIO->WriteNodalResults(WEIGHTED_SLIP, rModelPart.Nodes(), label, 0);
                mpIO->WriteNodalResultsNonHistorical(AUGMENTED_TANGENT_CONTACT_PRESSURE, rModelPart.Nodes(), label);
            }
        }

        return true;
    }

    /**
     * @brief This function initialize the convergence criteria
     * @param rModelPart The model part of interest
     */
    void Initialize(ModelPart& rModelPart) override
    {
        // Calling base criteria
        BaseType::Initialize(rModelPart);

        // The current process info
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        r_process_info.SetValue(ACTIVE_SET_COMPUTED, false);
    }

    /**
     * @brief This function initializes the solution step
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual)
     */
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {
        // Update normal of the conditions
        ModelPart& r_contact_model_part = rModelPart.GetSubModelPart("Contact");
        MortarUtilities::ComputeNodesMeanNormalModelPart(r_contact_model_part);
        const bool frictional_problem = rModelPart.IsDefined(SLIP) ? rModelPart.Is(SLIP) : false;
        if (frictional_problem) {
            const bool has_lm = rModelPart.HasNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER);
            if (has_lm && mOptions.IsNot(BaseMortarConvergenceCriteria::PURE_SLIP)) {
                MortarUtilities::ComputeNodesTangentModelPart(r_contact_model_part);
            } else {
                MortarUtilities::ComputeNodesTangentModelPart(r_contact_model_part, &WEIGHTED_SLIP, 1.0, true);
            }
        }

        // IO for debugging
        if (mOptions.Is(BaseMortarConvergenceCriteria::IO_DEBUG)) {
            mpIO->CloseResultFile();
            std::ostringstream new_name ;
            new_name << "POST_LINEAR_ITER_STEP=""POST_LINEAR_ITER_STEP=" << rModelPart.GetProcessInfo()[STEP];
            mpIO->ChangeOutputName(new_name.str());
        }
    }

    /**
     * @brief This function finalizes the solution step
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual)
     */
    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {
        // IO for debugging
        if (mOptions.Is(BaseMortarConvergenceCriteria::IO_DEBUG)) {
            mpIO->FinalizeResults();
        }
    }

    /**
     * @brief This function finalizes the non-linear iteration
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual + reactions)
     */
    void FinalizeNonLinearIteration(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {
        // Calling base criteria
        BaseType::FinalizeNonLinearIteration(rModelPart, rDofSet, rA, rDx, rb);

        // The current process info
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        r_process_info.SetValue(ACTIVE_SET_COMPUTED, false);
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

    Flags mOptions; /// Local flags

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
        NodesArrayType& r_nodes_array = rModelPart.GetSubModelPart("Contact").Nodes();
        VariableUtils().SetVariable(WEIGHTED_GAP, 0.0, r_nodes_array);
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

    GidIOBaseType::Pointer mpIO; /// The pointer to the debugging IO

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief It computes the mean of the normal in the condition in all the nodes
     * @param rModelPart The model part to compute
     */
    inline void ComputeNodesMeanNormalModelPartWithPairedNormal(ModelPart& rModelPart)
    {
        // Compute normal and tangent
        ModelPart& r_contact_model_part = rModelPart.GetSubModelPart("Contact");
        MortarUtilities::ComputeNodesMeanNormalModelPart(r_contact_model_part);

        // Iterate over the computing conditions
        ModelPart& r_computing_contact_model_part = rModelPart.GetSubModelPart("ComputingContact");
        ConditionsArrayType& r_conditions_array = r_computing_contact_model_part.Conditions();
        const auto it_cond_begin = r_conditions_array.begin();

        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i) {
            auto it_cond = it_cond_begin + i;

            // Aux coordinates
            Point::CoordinatesArrayType aux_coords;

            // We update the paired normal
            GeometryType& r_parent_geometry = it_cond->GetGeometry().GetGeometryPart(0);
            aux_coords = r_parent_geometry.PointLocalCoordinates(aux_coords, r_parent_geometry.Center());
            it_cond->SetValue(NORMAL, r_parent_geometry.UnitNormal(aux_coords));
        }
    }

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

///@name Local flags creation
///@{

/// Local Flags
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags BaseMortarConvergenceCriteria<TSparseSpace, TDenseSpace>::COMPUTE_DYNAMIC_FACTOR(Kratos::Flags::Create(0));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags BaseMortarConvergenceCriteria<TSparseSpace, TDenseSpace>::IO_DEBUG(Kratos::Flags::Create(1));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags BaseMortarConvergenceCriteria<TSparseSpace, TDenseSpace>::PURE_SLIP(Kratos::Flags::Create(2));

}  // namespace Kratos

#endif /* KRATOS_BASE_MORTAR_CRITERIA_H  defined */
