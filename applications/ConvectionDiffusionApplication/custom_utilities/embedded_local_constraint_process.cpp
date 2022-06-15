// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Franziska Wahl
//

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "processes/find_nodal_neighbours_process.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "embedded_local_constraint_process.h"

namespace Kratos
{

/* Public functions *******************************************************/

EmbeddedLocalConstraintProcess::EmbeddedLocalConstraintProcess(
    Model& rModel,
    Parameters ThisParameters)
    : Process()
{
    // Validate input settings with defaults
    ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    // Retrieve the required model parts
    const std::string model_part_name = ThisParameters["model_part_name"].GetString();
    mpModelPart = &rModel.GetModelPart(model_part_name);

    // Set whether nodal distances will be modified to avoid levelset zeros
    mAvoidZeroDistances = ThisParameters["avoid_zero_distances"].GetBool();

    // Set which elements will not be active
    mDeactivateNegativeElements = ThisParameters["deactivate_negative_elements"].GetBool();
    mDeactivateIntersectedElements = ThisParameters["deactivate_intersected_elements"].GetBool();
}

void EmbeddedLocalConstraintProcess::Execute()
{
    if (mAvoidZeroDistances) { ModifyDistances(); }

    // Deactivate intersected and negative elements as well as their nodes depending on the process settings
    DeactivateElementsAndNodes();

    // Apply a constraint to negative nodes of split elements
    ApplyConstraints();
}

/* Protected functions ****************************************************/

/* Private functions ****************************************************/

void EmbeddedLocalConstraintProcess::ApplyConstraints()
{
    // Initialize counter of master slave constraints
    ModelPart::IndexType id = mpModelPart->NumberOfMasterSlaveConstraints()+1;
    // Get variable to constrain  //TODO: different variable??
    const auto& r_var = KratosComponents<Variable<double>>::Get("TEMPERATURE");

    // Loop through all elements to get negative nodes of split elements (slave nodes)
    // TODO: only add constraint for small cut??
    for (auto& rElement : mpModelPart->Elements()) {
    auto& r_geom = rElement.GetGeometry();
        if (IsSplit(r_geom)) {  // && IsSmallCut(r_geom)) {

            //TODO: Calculate using boundary condition
            //const double temp_bc = rElement.GetValue(EMBEDDED_SCALAR);
            // Intersection point coordinates
            // TODO: get shape function values for intersection point
            // array_1d<double,3> xg_coords = ZeroVector(3);
            // for (std::size_t i = 0; i < TDim+1; ++i) {
            //     noalias(xg_coords) += N(i) * r_geom[i].Coordinates();
            // const double aux_temp_bc = std::pow(xg_coords[0],2) + std::pow(xg_coords[1],2);

            // Containers for negative and positive side nodes of split element
            std::vector<NodeType::Pointer> neg_nodes = {};
            std::vector<NodeType::Pointer> pos_nodes = {};
            // .clear().reserve(n_nodes);  .size()

            for (auto& rNode : r_geom) {
                if (rNode.FastGetSolutionStepValue(DISTANCE) > 0.0) {
                    pos_nodes.push_back(&rNode);
                } else {
                    neg_nodes.push_back(&rNode);
                }
            }
            const double support_node_weight = 0.0; //1.0 / pos_nodes.size();

            //TODO: Add constraints for a negative node only for one split element
            // --> dictonary of cut elements negative nodes, first add elements with 1 negative node, then 2, then 3 (3D), only add if node is not already there
            // Add one master-slave constraint for every positive distance node of the split element for all negative distance nodes
            // The contributions of each master will be summed up in the BuilderAndSolver to give an equation for the slave dof
            for(auto p_slave_node : neg_nodes) {
                for (auto p_support_node : pos_nodes) {
                    mpModelPart->CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", id++,
                    *p_support_node, r_var, *p_slave_node, r_var,
                    support_node_weight, 0.5);
                }
            }
        }
    }
}

void EmbeddedLocalConstraintProcess::DeactivateElementsAndNodes()
{
    // Initialize flags to true \\TODO: necessary??
    block_for_each(mpModelPart->Nodes(), [](NodeType& rNode){
        rNode.Set(ACTIVE, true);  // Nodes that belong to the elements to be assembled
    });
    block_for_each(mpModelPart->Elements(), [](Element& rElement){
        rElement.Set(ACTIVE, true);  // Elements in the positive distance region (the ones to be assembled)
    });

    // Deactivate intersected elements and their nodes  \\TODO: move to where the constraints are applied??
    if ( mDeactivateIntersectedElements ) {
        for (auto& rElement : mpModelPart->Elements()) {
            // Check if the element is split
            const auto& r_geom = rElement.GetGeometry();
            if (IsSplit(r_geom)) {
                rElement.Set(ACTIVE, false);
                for (auto& rNode : r_geom) {
                    rNode.Set(ACTIVE, false);
                }
            }
        }
    }
    // Deactivate negative elements and their nodes
    if ( mDeactivateNegativeElements ) {
        for (auto& rElement : mpModelPart->Elements()) {
            // Check if the element is negative
            const auto& r_geom = rElement.GetGeometry();
            if (IsNegative(r_geom)) {
                rElement.Set(ACTIVE, false);
                for (auto& rNode : r_geom) {
                    rNode.Set(ACTIVE, false);
                }
            }
        }
    }
}

void EmbeddedLocalConstraintProcess::ModifyDistances()
{
    auto& r_nodes = mpModelPart->Nodes();
    const double tol_d = 1.0e-12;

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
        auto it_node = r_nodes.begin() + i;
        double& d = it_node->FastGetSolutionStepValue(DISTANCE);

        // Check if the distance values are close to zero, if so set the tolerance as distance value
        if (std::abs(d) < tol_d) {
            d = (d > 0.0) ? tol_d : -tol_d;
        }
    }
}

bool EmbeddedLocalConstraintProcess::IsSplit(const GeometryType& rGeometry)
{
    std::size_t n_neg = 0;
    std::size_t n_pos = 0;
    for (const auto& r_node : rGeometry) {
        if (r_node.FastGetSolutionStepValue(DISTANCE) < 0.0) {
            n_neg++;
        } else {
            n_pos++;
        }
    }
    return (n_pos != 0 && n_neg != 0);
}

bool EmbeddedLocalConstraintProcess::IsSmallCut(const GeometryType& rGeometry)
{
    const double tol_d = 0.001;
    for (const auto& r_node : rGeometry) {
        if (abs(r_node.FastGetSolutionStepValue(DISTANCE)) < tol_d) {
            true;
        }
    }
    return false;
}

bool EmbeddedLocalConstraintProcess::IsNegative(const GeometryType& rGeometry)
{
    std::size_t n_neg = 0;
    for (const auto& r_node : rGeometry) {
        if (r_node.FastGetSolutionStepValue(DISTANCE) < 0.0) {
            n_neg++;
        }
    }
    return (n_neg == rGeometry.PointsNumber());
}

}; // namespace Kratos.
