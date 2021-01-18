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

#if !defined(KRATOS_MPC_CONTACT_CRITERIA_H)
#define  KRATOS_MPC_CONTACT_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "utilities/color_utilities.h"
#include "utilities/variable_utils.h"
#include "custom_utilities/contact_utilities.h"
#include "processes/simple_mortar_mapper_wrapper_process.h"

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
 * @class MPCContactCriteria
 * @ingroup ContactStructuralMechanicsApplication
 * @brief Custom convergence criteria for the contact problem
 * @author Vicente Mataix Ferrandiz
 */
template<class TSparseSpace, class TDenseSpace>
class MPCContactCriteria
    : public  ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MPCContactCriteria
    KRATOS_CLASS_POINTER_DEFINITION( MPCContactCriteria );

    /// The base class definition (and it subclasses)
    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace >                    BaseType;
    typedef typename BaseType::TDataType                                       TDataType;
    typedef typename BaseType::DofsArrayType                               DofsArrayType;
    typedef typename BaseType::TSystemMatrixType                       TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType                       TSystemVectorType;

    /// The sparse space used
    typedef TSparseSpace                                                 SparseSpaceType;

    /// The components containers
    typedef ModelPart::NodesContainerType                                 NodesArrayType;
    typedef ModelPart::ElementsContainerType                           ElementsArrayType;
    typedef ModelPart::ConditionsContainerType                       ConditionsArrayType;
    typedef ModelPart::MasterSlaveConstraintContainerType            ConstraintArrayType;

    /// The table stream definition TODO: Replace by logger
    typedef TableStreamUtility::Pointer                          TablePrinterPointerType;

    /// The index type definition
    typedef std::size_t                                                        IndexType;

    // Geometry definition
    typedef Node<3>                                                             NodeType;
    typedef Geometry<NodeType>                                              GeometryType;
    typedef CouplingGeometry<NodeType>                              CouplingGeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    explicit MPCContactCriteria()
        : BaseType()
    {
    }

    ///Copy constructor
    MPCContactCriteria( MPCContactCriteria const& rOther )
      : BaseType(rOther)
    {
    }

    /// Destructor
    ~MPCContactCriteria() override = default;

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
        BaseType::PreCriteria(rModelPart, rDofSet, rA, rDx, rb);

        // Auxiliar zero array
        const array_1d<double, 3> zero_array = ZeroVector(3);

        // We initailize the contact force
        NodesArrayType& r_nodes_array = rModelPart.GetSubModelPart("Contact").Nodes();
        const auto it_node_begin = r_nodes_array.begin();

        // We save the current WEIGHTED_GAP in the buffer and reset the CONTACT_FORCE
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
            auto it_node = it_node_begin + i;
            it_node->SetValue(CONTACT_FORCE, zero_array);
            it_node->FastGetSolutionStepValue(WEIGHTED_GAP, 1) = it_node->FastGetSolutionStepValue(WEIGHTED_GAP);
        }

        // Compute weighted gap
        ComputeWeightedGap(rModelPart);

        // Reset the NODAL_AREA
        VariableUtils().SetNonHistoricalVariableToZero(NODAL_AREA, r_nodes_array);

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
        // We call the base class
        BaseType::PostCriteria(rModelPart, rDofSet, rA, rDx, rb);

        // Getting process info
        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        if (r_process_info[NL_ITERATION_NUMBER] > 0) {
            // Getting REACTION_CHECK_STIFFNESS_FACTOR
            const double reaction_check_stiffness_factor = r_process_info.Has(REACTION_CHECK_STIFFNESS_FACTOR) ?  r_process_info.GetValue(REACTION_CHECK_STIFFNESS_FACTOR) : 1.0e-12;

            // Compute weighted gap
            ComputeWeightedGap(rModelPart);

            // Transfer reaction from master to slave
            std::size_t sub_contact_counter = 0;
            CounterContactModelParts(rModelPart, sub_contact_counter);

            // Mapping reaction
            Parameters mapping_parameters = Parameters(R"({"distance_threshold" : 1.0e24, "update_interface" : false, "origin_variable" : "REACTION", "mapping_coefficient" : -1.0e0})" );
            if (r_process_info.Has(DISTANCE_THRESHOLD)) {
                mapping_parameters["distance_threshold"].SetDouble(r_process_info[DISTANCE_THRESHOLD]);
            }
            auto& r_contact_model_part = rModelPart.GetSubModelPart("Contact");
            for (std::size_t i_contact = 0; i_contact < sub_contact_counter; ++i_contact) {
                auto& r_sub = r_contact_model_part.GetSubModelPart("ContactSub" + std::to_string(i_contact));
                auto& r_sub_master = r_sub.GetSubModelPart("MasterSubModelPart" + std::to_string(i_contact));
                auto& r_sub_slave = r_sub.GetSubModelPart("SlaveSubModelPart" + std::to_string(i_contact));
                SimpleMortarMapperProcessWrapper(r_sub_master, r_sub_slave, mapping_parameters).Execute();
            }

            // TODO: Add frictional check

            // Getting process info
            Properties::Pointer p_properties = rModelPart.Elements().begin()->pGetProperties();
            for (auto& r_elements : rModelPart.Elements()) {
                if (r_elements.pGetProperties()->Has(YOUNG_MODULUS)) {
                    p_properties = r_elements.pGetProperties();
                }
            }

            // Defining the convergence
            IndexType is_active_set_converged = 0, is_slip_converged = 0;

            // Checking just after first iteration
            // We get the YOUNG_MODULUS
            const double young_modulus = p_properties->Has(YOUNG_MODULUS) ? p_properties->GetValue(YOUNG_MODULUS) : 0.0;
            const double auxiliar_check = young_modulus > 0.0 ? -(reaction_check_stiffness_factor * young_modulus) : 0.0;

            // We check the active/inactive set during the first non-linear iteration or for the general semi-smooth case
            NodesArrayType& r_nodes_array = r_contact_model_part.Nodes();
            const auto it_node_begin = r_nodes_array.begin();

            // If frictionaless or mesh tying
            if (rModelPart.IsNot(SLIP)) {
                #pragma omp parallel for reduction(+:is_active_set_converged)
                for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
                    auto it_node = it_node_begin + i;

                    if (it_node->Is(SLAVE)) {
                        // The contact force corresponds with the reaction in the normal direction
                        const array_1d<double, 3>& r_total_force = it_node->FastGetSolutionStepValue(REACTION);

                        const double nodal_area = it_node->Has(NODAL_AREA) ? it_node->GetValue(NODAL_AREA) : 1.0;
                        const double gap = it_node->FastGetSolutionStepValue(WEIGHTED_GAP)/nodal_area;
                        const array_1d<double, 3>& r_normal = it_node->FastGetSolutionStepValue(NORMAL);
                        const double contact_force = inner_prod(r_total_force, r_normal);
                        const double contact_pressure = contact_force/it_node->GetValue(NODAL_MAUX);

                        if (contact_pressure < auxiliar_check || gap < 0.0) { // NOTE: This could be conflictive (< or <=)
                            // We save the contact force
                            it_node->SetValue(CONTACT_FORCE, contact_force/it_node->GetValue(NODAL_PAUX) * r_normal);
                            it_node->SetValue(NORMAL_CONTACT_STRESS, contact_pressure);
                            if (it_node->IsNot(ACTIVE)) {
                                it_node->Set(ACTIVE, true);
                                is_active_set_converged += 1;
                            }
                        } else {
                            if (it_node->Is(ACTIVE)) {
                                it_node->Set(ACTIVE, false);
                                is_active_set_converged += 1;
                            }
                        }
                    }
                }
            } else { // If frictional
                #pragma omp parallel for reduction(+:is_active_set_converged, is_slip_converged)
                for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
                    auto it_node = it_node_begin + i;

                    if (it_node->Is(SLAVE)) {
                        const double auxiliar_check = young_modulus > 0.0 ? -(reaction_check_stiffness_factor * young_modulus) : 0.0;
                        // The contact force corresponds with the reaction in the normal direction
                        const array_1d<double, 3>& r_total_force = it_node->FastGetSolutionStepValue(REACTION);

                        const double nodal_area = it_node->Has(NODAL_AREA) ? it_node->GetValue(NODAL_AREA) : 1.0;
                        const double gap = it_node->FastGetSolutionStepValue(WEIGHTED_GAP)/nodal_area;
                        const array_1d<double, 3>& r_normal = it_node->FastGetSolutionStepValue(NORMAL);
                        const double contact_force = inner_prod(r_total_force, r_normal);
                        const double normal_contact_pressure = contact_force/it_node->GetValue(NODAL_MAUX);

                        if (normal_contact_pressure < auxiliar_check || gap < 0.0) { // NOTE: This could be conflictive (< or <=)
                            // We save the contact force
                            it_node->SetValue(CONTACT_FORCE, r_total_force/it_node->GetValue(NODAL_PAUX));
                            it_node->SetValue(NORMAL_CONTACT_STRESS, normal_contact_pressure);
                            if (it_node->IsNot(ACTIVE)) {
                                it_node->Set(ACTIVE, true);
                                is_active_set_converged += 1;
                            }

                            // The friction coefficient
                            const double tangential_contact_pressure = norm_2(r_total_force - contact_force * r_normal)/it_node->GetValue(NODAL_MAUX);
                            const bool is_slip = it_node->Is(SLIP);
                            const double mu = it_node->GetValue(FRICTION_COEFFICIENT);

                            if (tangential_contact_pressure <= - mu * contact_force) { // STICK CASE // TODO: Check the <=
                                it_node->SetValue(TANGENTIAL_CONTACT_STRESS, tangential_contact_pressure);
                                if (is_slip) {
                                    it_node->Set(SLIP, false);
                                    is_slip_converged += 1;
                                }
                            } else { // SLIP CASE
                                it_node->SetValue(TANGENTIAL_CONTACT_STRESS, - mu * contact_force);
                                if (!is_slip) {
                                    it_node->Set(SLIP, true);
                                    is_slip_converged += 1;
                                }
                            }
                        } else {
                            if (it_node->Is(ACTIVE)) {
                                it_node->Set(ACTIVE, false);
                                it_node->Reset(SLIP);
                                is_active_set_converged += 1;
                            }
                        }
                    }
                }
            }

            // We set the constraints active and inactive in function of the active set
            ConditionsArrayType& r_conditions_array = rModelPart.GetSubModelPart("ComputingContact").Conditions();
            auto it_cond_begin = r_conditions_array.begin();

            #pragma omp parallel for
            for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i) {
                auto it_cond = it_cond_begin + i;

                const auto& r_slave_geometry = it_cond->GetGeometry().GetGeometryPart(CouplingGeometryType::Master);
                std::size_t counter = 0;
                for (auto& r_node : r_slave_geometry) {
                    if (r_node.IsNot(ACTIVE)) {
                        ++counter;
                    }
                }

                // In case of traction we deactivate
                if (counter == r_slave_geometry.size()) {
                    it_cond->Set(ACTIVE, false);
                    // We deactivate the constraints on inactive conditions
                    if (it_cond->Has(CONSTRAINT_POINTER)) {
                        auto p_const = it_cond->GetValue(CONSTRAINT_POINTER);

                        // In case of traction we deactivate
                        p_const->Set(ACTIVE, false);
                    } else {
                        KRATOS_ERROR << "Contact conditions must have defined CONSTRAINT_POINTER" << std::endl;
                    }
                }
            }

            // We save to the process info if the active set has converged
            const bool active_set_converged = (is_active_set_converged == 0 ? true : false);
            const bool slip_set_converged = (is_slip_converged == 0 ? true : false);

            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
                if (active_set_converged) {
                    KRATOS_INFO("MPCContactCriteria")  << BOLDFONT("\tActive set") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                } else {
                    KRATOS_INFO("MPCContactCriteria")  << BOLDFONT("\tActive set") << " convergence is " << BOLDFONT(FRED("not achieved")) << std::endl;
                }
                if (slip_set_converged) {
                    KRATOS_INFO("MPCContactCriteria")  << BOLDFONT("\tSlip set") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                } else {
                    KRATOS_INFO("MPCContactCriteria")  << BOLDFONT("\tSlip set") << " convergence is " << BOLDFONT(FRED("not achieved")) << std::endl;
                }
            }

            return (active_set_converged && slip_set_converged);
        }

        return true;
    }

    /**
     * @brief This function initialize the convergence criteria
     * @param rModelPart The model part of interest
     */
    void Initialize(ModelPart& rModelPart) override
    {
        BaseType::Initialize(rModelPart);
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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method computes the weighted gap in the nodes of the problem
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     */
    void ComputeWeightedGap(ModelPart& rModelPart)
    {
        NodesArrayType& r_nodes_array = rModelPart.GetSubModelPart("Contact").Nodes();
        // Set to zero the weighted gap
        if (rModelPart.Is(SLIP)) {
            // Reset
            VariableUtils().SetHistoricalVariableToZero(WEIGHTED_GAP, r_nodes_array);
            VariableUtils().SetHistoricalVariableToZero(WEIGHTED_SLIP, r_nodes_array);
        } else {
            VariableUtils().SetHistoricalVariableToZero(WEIGHTED_GAP, r_nodes_array);
        }

        // Compute the contribution
        ContactUtilities::ComputeExplicitContributionConditions(rModelPart.GetSubModelPart("ComputingContact"));
    }

    /**
     * @brief This method computes the weighted gap in the nodes of the problem
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rCounter Reference to the counter
     */
    void CounterContactModelParts(
        ModelPart& rModelPart,
        std::size_t& rCounter
        )
    {
        for (auto& r_name : rModelPart.GetSubModelPartNames()) {
            if (r_name.find("ContactSub") != std::string::npos && r_name.find("ComputingContactSub") == std::string::npos) {
                ++rCounter;
            }
            auto& r_sub = rModelPart.GetSubModelPart(r_name);
            if (r_sub.IsSubModelPart()) {
                CounterContactModelParts(r_sub, rCounter);
            }
        }
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Unaccessible methods
    ///@{

    ///@}

}; // Class MPCContactCriteria

///@name Explicit Specializations
///@{

}  // namespace Kratos

#endif /* KRATOS_MPC_CONTACT_CRITERIA_H  defined */

