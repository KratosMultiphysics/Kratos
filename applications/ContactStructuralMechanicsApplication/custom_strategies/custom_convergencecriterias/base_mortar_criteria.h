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

    /// Pointer definition of BaseMortarConvergenceCriteria
    KRATOS_CLASS_POINTER_DEFINITION( BaseMortarConvergenceCriteria );

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
        const bool IODebug = false
        )
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >(),
          mComputeDynamicFactor(ComputeDynamicFactor),
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
        ProcessInfo& process_info = rModelPart.GetProcessInfo();

        // The contact model part
        ModelPart& r_contact_model_part = rModelPart.GetSubModelPart("Contact");

        // We update the normals if necessary
        const auto normal_variation = process_info.Has(CONSIDER_NORMAL_VARIATION) ? static_cast<NormalDerivativesComputation>(process_info.GetValue(CONSIDER_NORMAL_VARIATION)) : NO_DERIVATIVES_COMPUTATION;
        if (normal_variation != NO_DERIVATIVES_COMPUTATION)
            ComputeNodesMeanNormalModelPartWithPairedNormal(rModelPart); // Update normal of the conditions

        const bool adapt_penalty = process_info.Has(ADAPT_PENALTY) ? process_info.GetValue(ADAPT_PENALTY) : false;
        const bool dynamic_case = rModelPart.HasNodalSolutionStepVariable(VELOCITY);

        /* Compute weighthed gap */
        if (adapt_penalty || dynamic_case) {
            // Set to zero the weighted gap
            ResetWeightedGap(rModelPart);

            // Compute the contribution
            ContactUtilities::ComputeExplicitContributionConditions(rModelPart.GetSubModelPart("ComputingContact"));
        }

        // In dynamic case
        if ( dynamic_case && mComputeDynamicFactor) {
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
        if (mIODebug) {
            const bool frictional_problem = rModelPart.IsDefined(SLIP) ? rModelPart.Is(SLIP) : false;
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
            if (rModelPart.Nodes().begin()->SolutionStepsDataHas(VELOCITY_X)) {
                mpGidIO->WriteNodalResults(VELOCITY, rModelPart.Nodes(), label, 0);
                mpGidIO->WriteNodalResults(ACCELERATION, rModelPart.Nodes(), label, 0);
            }
            if (r_nodes_array.begin()->SolutionStepsDataHas(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE))
                mpGidIO->WriteNodalResults(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE, rModelPart.Nodes(), label, 0);
            else if (r_nodes_array.begin()->SolutionStepsDataHas(VECTOR_LAGRANGE_MULTIPLIER_X))
                mpGidIO->WriteNodalResults(VECTOR_LAGRANGE_MULTIPLIER, rModelPart.Nodes(), label, 0);
            mpGidIO->WriteNodalResults(WEIGHTED_GAP, rModelPart.Nodes(), label, 0);
            if (frictional_problem) {
                mpGidIO->WriteNodalFlags(SLIP, "SLIP", rModelPart.Nodes(), label);
                mpGidIO->WriteNodalResults(WEIGHTED_SLIP, rModelPart.Nodes(), label, 0);
                mpGidIO->WriteNodalResultsNonHistorical(AUGMENTED_TANGENT_CONTACT_PRESSURE, rModelPart.Nodes(), label);
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
        KRATOS_ERROR << "YOUR ARE CALLING THE BASE MORTAR CRITERIA" << std::endl;
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
        MortarUtilities::ComputeNodesMeanNormalModelPart(rModelPart.GetSubModelPart("Contact"));

        // GiD IO for debugging
        if (mIODebug) {
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
        // GiD IO for debugging
        if (mIODebug)
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
        NodesArrayType& r_nodes_array = rModelPart.GetSubModelPart("Contact").Nodes();
        VariableUtils().SetScalarVar<Variable<double>>(WEIGHTED_GAP, 0.0, r_nodes_array);
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

    bool mComputeDynamicFactor;     /// If we compute the dynamic factor
    bool mIODebug;                  /// If we generate an output gid file in order to debug
    GidIOBaseType::Pointer mpGidIO; /// The pointer to the debugging GidIO

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
    static inline void ComputeNodesMeanNormalModelPartWithPairedNormal(ModelPart& rModelPart) {
        MortarUtilities::ComputeNodesMeanNormalModelPart(rModelPart.GetSubModelPart("Contact"));

        // Iterate over the computing conditions
        ConditionsArrayType& r_conditions_array = rModelPart.GetSubModelPart("ComputingContact").Conditions();
        const auto it_cond_begin = r_conditions_array.begin();

        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i) {
            auto it_cond = it_cond_begin + i;

            // Aux coordinates
            Point::CoordinatesArrayType aux_coords;

            // We update the paired normal
            GeometryType& this_geometry = it_cond->GetGeometry();
            aux_coords = this_geometry.PointLocalCoordinates(aux_coords, this_geometry.Center());
            it_cond->SetValue(NORMAL, this_geometry.UnitNormal(aux_coords));

            // We update the paired normal
            GeometryType::Pointer p_paired_geometry = it_cond->GetValue(PAIRED_GEOMETRY);
            aux_coords = p_paired_geometry->PointLocalCoordinates(aux_coords, p_paired_geometry->Center());
            it_cond->SetValue(PAIRED_NORMAL, p_paired_geometry->UnitNormal(aux_coords));
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

///@name Explicit Specializations
///@{

}  // namespace Kratos

#endif /* KRATOS_BASE_MORTAR_CRITERIA_H  defined */

