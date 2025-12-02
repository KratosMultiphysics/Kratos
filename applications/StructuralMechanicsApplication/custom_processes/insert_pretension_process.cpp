// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Máté Kelemen
//

// --- Structural Includes ---
#include "custom_processes/insert_pretension_process.hpp" // InsertPretensionProcess

// --- Core Includes ---
#include "includes/master_slave_constraint.h" // MasterSlaveConstraint
#include "includes/model_part.h" // ModelPart
#include "processes/find_global_nodal_elemental_neighbours_process.h"
#include "utilities/comparison.h" // Comparison
#include "includes/element.h" // Element
//#include "geometries/point_3d.h" // Point3D

// --- STL Includes ---
#include <array> // std::array
#include <sstream> // std::stringstream


namespace Kratos {


//class PretensionElement : public Element {
//
//}; // class PretensionElement


struct InsertPretensionProcess::Impl {
    ModelPart* mpModelPart;
    std::string mImpositionName;
    double mPretensionValue;
}; // struct InsertPretensionProcess::Impl


InsertPretensionProcess::InsertPretensionProcess(Model& rModel, Parameters Settings)
    : Process(),
      mpImpl(new Impl) {
    Settings.ValidateAndAssignDefaults(this->GetDefaultParameters());

    KRATOS_TRY
    mpImpl->mpModelPart = &rModel.GetModelPart(Settings["model_part_name"].GetString());
    mpImpl->mImpositionName = Settings["imposition"].GetString();
    mpImpl->mPretensionValue = Settings["pretension_value"].GetDouble();
    KRATOS_CATCH("")

    const std::array<std::string,3> valid_imposition_names {
        std::string {"lagrange_element"},
        std::string {"penalty_element"},
        std::string {"constraint"}
    };

    if (std::find(valid_imposition_names.begin(),
                  valid_imposition_names.end(),
                  mpImpl->mImpositionName) == valid_imposition_names.end()) {
        std::stringstream message;
        message << "Invalid setting for \"imposition\". Options are: ";
        for (const auto& r_imposition_name : valid_imposition_names)
            message << '"' << r_imposition_name << "\" ";
        KRATOS_ERROR << message.str();
    }

    // This process expects the pre-tension surface to be defined
    // by a set of Conditions in the input ModelPart. It requires
    // at least one Condition.
    KRATOS_ERROR_IF(mpImpl->mpModelPart->Conditions().empty())
        << "InsertPretensionProcess requires at least one Condition in the provided ModelPart, "
        << "but there aren't any in '" << mpImpl->mpModelPart->Name() << "'.";
}


InsertPretensionProcess::InsertPretensionProcess(InsertPretensionProcess&&) noexcept = default;


InsertPretensionProcess::~InsertPretensionProcess() = default;


struct PretensionSurfacePartition {
    array_1d<double,3> normal;
    std::vector<decltype(NEIGHBOUR_ELEMENTS)::Type::pointer> positive_side, negative_side;
}; // struct SurfacePartition


void FindElementsAdjacentToConditions(ModelPart& rModelPart) {
    using NeighborElements = decltype(NEIGHBOUR_ELEMENTS)::Type;

    // Make sure every node on the pre-tension surface is in the model part.
    KRATOS_TRY
        ModelPart::NodesContainerType nodes;
        for (Condition& r_condition : rModelPart.Conditions()) {
            auto& r_geometry = r_condition.GetGeometry();
            nodes.insert(r_geometry.begin(), r_geometry.end());
        } // for r_condition in rModelPart.Conditions
        rModelPart.Nodes().insert(nodes);
    KRATOS_CATCH("")

    KRATOS_TRY
        // Find elements containing each node in the model part.
        FindGlobalNodalElementalNeighboursProcess(rModelPart.GetRootModelPart()).Execute();

        block_for_each(
            rModelPart.Conditions(),
            [] (Condition& r_condition) -> void {
                NeighborElements neighbor_elements;

                // Collect all neighbor elements from the condition's nodes.
                for (Node& r_node : r_condition.GetGeometry()) {
                    NeighborElements& r_containing_elements = r_node.GetValue(NEIGHBOUR_ELEMENTS);
                    for (unsigned i_element=0u; i_element<r_containing_elements.size(); ++i_element) {
                        const auto& rp_element = *(r_containing_elements.ptr_begin() + i_element);
                        const auto it_element = std::find_if(
                            neighbor_elements.ptr_begin(),
                            neighbor_elements.ptr_end(),
                            [&rp_element](const auto& rp_other_element) -> bool {
                                return rp_other_element->Id() == rp_element->Id();
                            });
                        if (it_element == neighbor_elements.ptr_end()) {
                            neighbor_elements.push_back(rp_element);
                        }
                    }
                }

                // Filter elements that don't contain all nodes of the condition.
                neighbor_elements.erase(std::remove_if(
                    neighbor_elements.begin(),
                    neighbor_elements.end(),
                    [&r_condition] (const Element& r_element) -> bool {
                        for (const Node& r_condition_node : r_condition.GetGeometry()) {
                            const auto it_node = std::find_if(
                                r_element.GetGeometry().begin(),
                                r_element.GetGeometry().end(),
                                [&r_condition_node] (const Node& r_element_node) -> bool {
                                    return r_element_node.Id() == r_condition_node.Id();
                                });
                            if (it_node == r_element.GetGeometry().end())
                                return true;
                        }
                        return false;
                    }
                ), neighbor_elements.end());

                // Set neighbor elements on the condition.
                r_condition.SetValue(NEIGHBOUR_ELEMENTS, neighbor_elements);
            });
    KRATOS_CATCH("")
}


template <unsigned SurfaceDimension>
PretensionSurfacePartition PretensionPlanePartition(const ModelPart::ConditionsContainerType::iterator itConditionBegin,
                                                    const ModelPart::ConditionsContainerType::iterator itConditionEnd) {
    PretensionSurfacePartition pretension_surface;

    // There are three valid possibilities here.
    // 1) All surfaces are 0-dimensional (corners)
    //    => expecting exactly 1 condition
    //    => expecting exactly 2 connected elements, both
    //       must have line geometries
    //    => the normal is not well-defined, so it is chosen
    //       as the average of the connected elements' directions
    // 2) All surfaces are 1-dimensional (edges)
    //    => expecting at least 1 condition
    //    => expecting at least 2 connected elements, all of which
    //       must be 2-dimensional
    //    => the normal is not well-defined, so it is chosen as the
    //       average of the rotated endpoints in 2D, in positive
    //       direction. This assumes coordinates in z-direction
    //       to vanish. This assumption is checked.
    // 3) All surfaces are 2-dimensional (surfaces)
    //    => expecting at least 1 condition
    //    => expecting at least 2 connected elements, all
    //       of which must be 3-dimensional
    //    => normals are well-defined
    // In all cases, at least one connected element is required on
    // either side of the average pre-tension plane.

    // Sanity checks.
    if constexpr (SurfaceDimension == 0u) {
        KRATOS_ERROR_IF_NOT(std::distance(itConditionBegin, itConditionEnd) == 1)
            << "Expecting exactly 1 condition on the 0-dimensional pre-tension surface , but got "
            << std::distance(itConditionBegin, itConditionEnd) << ".";
    } else if constexpr (SurfaceDimension == 1u) {
        KRATOS_ERROR_IF_NOT(1 <= std::distance(itConditionBegin, itConditionEnd))
            << "Expecting at least 1 condition on the 1-dimensional pre-tension surface , but got "
            << std::distance(itConditionBegin, itConditionEnd) << ".";
    } else if constexpr (SurfaceDimension == 2u) {
        KRATOS_ERROR_IF_NOT(1 <= std::distance(itConditionBegin, itConditionEnd))
            << "Expecting at least 1 condition on the 2-dimensional pre-tension surface , but got "
            << std::distance(itConditionBegin, itConditionEnd) << ".";
    } else {
        static_assert(float(SurfaceDimension) == 0.5f, "Invalid surface dimension");
    }

    for (auto it_condition=itConditionBegin; it_condition!=itConditionEnd; ++it_condition) {
        Condition& r_condition = *it_condition;

        // Sanity checks on the condition.
        KRATOS_ERROR_IF_NOT(r_condition.GetGeometry().LocalSpaceDimension() == SurfaceDimension)
            << "Expecting all geometries defining the pre-tension surface to be " << SurfaceDimension << " dimensional, "
            << "but condition " << r_condition.Id() << " on geometry " << r_condition.GetGeometry().Id() << " (" << r_condition.GetGeometry().Name() << ") "
            << "is " << r_condition.GetGeometry().LocalSpaceDimension() << "-dimensional.";

        KRATOS_ERROR_IF_NOT(r_condition.Has(NEIGHBOUR_ELEMENTS))
            << "Expecting condition " << r_condition.Id() << " to have exactly "
            << "2 neighboring elements, but found none.";
        auto& r_neighbor_elements = r_condition.GetValue(NEIGHBOUR_ELEMENTS);
        if (r_neighbor_elements.size() != 2) {
            std::stringstream message;
            message << "Expecting condition " << r_condition.Id() << " to have exactly "
                    << "2 neighboring elements, but found " << r_neighbor_elements.size() << " [";
            for (const Element& r_element : r_neighbor_elements) message << r_element.Id() << " ";
            message << "].";
            KRATOS_ERROR << message.str();
        }

        // Check whether the 2 connected elements lie on opposite sides of the surface.
        const array_1d<double,3>
            first_connected_direction  = r_neighbor_elements[0].GetGeometry().Center() - r_condition.GetGeometry().Center(),
            second_connected_direction = r_neighbor_elements[1].GetGeometry().Center() - r_condition.GetGeometry().Center();
        const double direction_inner_product = std::inner_product(
            first_connected_direction.begin(),
            first_connected_direction.end(),
            second_connected_direction.begin(),
            0.0);
        KRATOS_ERROR_IF_NOT(direction_inner_product < 0.0)
            << "Element " << r_neighbor_elements[0].Id() << " and " << r_neighbor_elements[1].Id() << " "
            << "do not lie on opposite sides of condition " << r_condition.Id() << " "
            << "defining a pre-tension surface.";

        if constexpr (SurfaceDimension == 0u) {
            for (unsigned i_element=0u; i_element<r_neighbor_elements.size(); ++i_element) {
                const Element& r_element = r_neighbor_elements[i_element];
                KRATOS_ERROR_IF_NOT(r_element.GetGeometry().GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Linear)
                    << "Expecting only " << SurfaceDimension + 1 << "-dimensional elements connected to the " << SurfaceDimension << "-dimensional pre-tension surface, "
                    << "but element " << r_element.Id() << " has a geometry of type " << r_element.GetGeometry().Name() << ".";

                const std::size_t surface_node_id = itConditionBegin->GetGeometry()[0].Id();
                const array_1d<double,3> surface_node_position = itConditionBegin->GetGeometry()[0];
                array_1d<double,3> opposite_node_position;
                if (r_element.GetGeometry()[0].Id() == surface_node_id) {
                    opposite_node_position = r_element.GetGeometry()[1];
                } else if (r_element.GetGeometry()[1].Id() == surface_node_id) {
                    opposite_node_position = r_element.GetGeometry()[0];
                } else {
                    KRATOS_ERROR << "Element " << r_element.Id() << " is assumed to be connected to node "
                                 << surface_node_id << " but is not.";
                }

                if (i_element % 2) {
                    pretension_surface.positive_side.push_back(*(r_neighbor_elements.ptr_begin() + i_element));
                    pretension_surface.normal = opposite_node_position - surface_node_position;
                } else {
                    pretension_surface.negative_side.push_back(*(r_neighbor_elements.ptr_begin() + i_element));
                }
            } // for i_element in range(r_neighbor_elements.size())
        } /*if SurfaceDimension == 0*/ else if constexpr (SurfaceDimension == 1u) {
            // Define the line condition's plane.
            const array_1d<double,3> surface_begin = r_condition.GetGeometry()[0],
                                     surface_end   = *(r_condition.GetGeometry().end() - 1);
            KRATOS_ERROR_IF_NOT(Comparison<double>::Equal(1e-12, 1e-9)(surface_begin[2], 0.0))
                << "Assuming that 1-dimensional pre-tension surfaces lie on the XY plane, but node "
                << r_condition.GetGeometry()[0].Id() << " of condition " << r_condition.Id() << " "
                << "has a Z-coordinate of " << surface_begin[2] << ".";
            KRATOS_ERROR_IF_NOT(Comparison<double>::Equal(1e-12, 1e-9)(surface_end[2], 0.0))
                << "Assuming that 1-dimensional pre-tension surfaces lie on the XY plane, but node "
                << r_condition.GetGeometry()[1].Id() << " of condition " << r_condition.Id() << " "
                << "has a Z-coordinate of " << surface_end[2] << ".";
            const array_1d<double,3>
                surface_normal {
                    -(surface_end[1] - surface_begin[1]),
                    surface_end[0] - surface_begin[0],
                    0.0},
                surface_center = r_condition.GetGeometry().Center();
            pretension_surface.normal += surface_normal;

            // Sort connected elements.
            for (unsigned i_element=0u; i_element<r_neighbor_elements.size(); ++i_element) {
                const Element& r_element = r_neighbor_elements[i_element];
                const array_1d<double,3> element_direction = r_element.GetGeometry().Center() - surface_center;
                const double direction_inner_product = std::inner_product(
                    surface_normal.begin(),
                    surface_normal.end(),
                    element_direction.begin(),
                    0.0);
                if (0 < direction_inner_product) pretension_surface.negative_side.push_back(*(r_neighbor_elements.ptr_begin() + i_element));
                else pretension_surface.positive_side.push_back(*(r_neighbor_elements.ptr_begin() + i_element));
            } // for i_element in range(r_neighbor_elements.size())

            KRATOS_ERROR_IF_NOT(pretension_surface.positive_side.size() == pretension_surface.negative_side.size())
                << "Elements " << r_neighbor_elements[0].Id() << " and " << r_neighbor_elements[1].Id() << " "
                << "lie on the same side of condition " << r_condition.Id() << " that defines a pre-tension surface.";
        } /*else if SurfaceDimension == 1*/ else if constexpr (SurfaceDimension == 2u) {
            // Define the surface's plane.
            const array_1d<double,3> surface_center = r_condition.GetGeometry().Center();
            const array_1d<double,3> surface_normal = r_condition.GetGeometry().Normal(surface_center);
            pretension_surface.normal += surface_normal;

            // Sort connected elements.
            for (unsigned i_element=0u; i_element<r_neighbor_elements.size(); ++i_element) {
                const Element& r_element = r_neighbor_elements[i_element];
                const array_1d<double,3> element_direction = r_element.GetGeometry().Center() - surface_center;
                const double direction_inner_product = std::inner_product(
                    surface_normal.begin(),
                    surface_normal.end(),
                    element_direction.begin(),
                    0.0);
                if (0 < direction_inner_product) pretension_surface.negative_side.push_back(*(r_neighbor_elements.ptr_begin() + i_element));
                else pretension_surface.positive_side.push_back(*(r_neighbor_elements.ptr_begin() + i_element));
            } // for i_element in range(r_neighbor_elements.size())

            KRATOS_ERROR_IF_NOT(pretension_surface.positive_side.size() == pretension_surface.negative_side.size())
                << "Elements " << r_neighbor_elements[0].Id() << " and " << r_neighbor_elements[1].Id() << " "
                << "lie on the same side of condition " << r_condition.Id() << " that defines a pre-tension surface.";
        } /*else if SurfaceDimension == 2*/
    } // for r_condition in r_model_part.Conditions

    return pretension_surface;
}


void InsertPretensionProcess::ExecuteBeforeSolutionLoop() {
    ModelPart& r_model_part = *mpImpl->mpModelPart;

    // Find elements connected to the pre-tension surface.
    // => each node will have a NEIGHBOR_ELEMENTS variable
    //    in its non-historical storage that lists elements
    //    containing them.
    KRATOS_TRY
        FindElementsAdjacentToConditions(r_model_part);
    KRATOS_CATCH("")

    // Sort connected elements into 2 categories, lying on either side of
    // the pre-tension surface.
    PretensionSurfacePartition pretension_surface;
    KRATOS_TRY
        // Find the dimension of the pre-tension surface,
        // and make sure all conditions share it.
        KRATOS_ERROR_IF(r_model_part.Conditions().empty());
        const int pre_tension_surface_dimension = r_model_part.Conditions().front().GetGeometry().LocalSpaceDimension();
        switch (pre_tension_surface_dimension) {
            case 0u: {
                pretension_surface = PretensionPlanePartition<0>(r_model_part.Conditions().begin(), r_model_part.Conditions().end());
                break;
            }
            case 1u: {
                pretension_surface = PretensionPlanePartition<1>(r_model_part.Conditions().begin(), r_model_part.Conditions().end());
                break;
            }
            case 2u: {
                pretension_surface = PretensionPlanePartition<2>(r_model_part.Conditions().begin(), r_model_part.Conditions().end());
                break;
            }
            default: KRATOS_ERROR << "Invalid pre-tension surface dimension: " << pre_tension_surface_dimension << ".";
        }
    KRATOS_CATCH("")

    // Duplicate nodes on the pre-tension surface and collect relevant displacement components.
    std::unordered_map<Node::IndexType,Node::Pointer> duplicated_nodes;
    KRATOS_TRY
        std::unordered_set<Node*> surface_nodes;

        // Collect nodes on the surface.
        for (Condition& r_condition : r_model_part.Conditions()) {
            for (Node& r_node : r_condition.GetGeometry()) {
                surface_nodes.insert(&r_node);
            } // for r_node in r_condition.GetGeometry()
        } // for r_condition in r_model_part.Conditions()

        // Find the first available node ID that can safely be
        // used to insert new ones.
        Node::IndexType node_id = r_model_part.GetRootModelPart().Nodes().back().Id() + 1;

        // Clone surface nodes and associate them with the original ones.
        for (Node* p_node : surface_nodes) {
            const auto emplace_result = duplicated_nodes.emplace(
                p_node->Id(),
                p_node->Clone(node_id++));
            r_model_part.AddNode(emplace_result.first->second);
        } // for p_node in surface_nodes

        // Replace nodes on the negative side elements with duplicated ones.
        for (auto& rp_element : pretension_surface.negative_side) {
            for (auto itp_node=rp_element->GetGeometry().ptr_begin(); itp_node!=rp_element->GetGeometry().ptr_end(); ++itp_node) {
                const auto it_duplicate_pair = duplicated_nodes.find((*itp_node)->Id());
                if (it_duplicate_pair != duplicated_nodes.end()) {
                    (*itp_node) = it_duplicate_pair->second;
                }
            } // for itp_node in rp_element->GetGeometry()
        } // for rp_element in pretension_surface.negative_side
    KRATOS_CATCH("")

    // Construct the constraint equation.
    struct PretensionConstraintData {
        std::vector<Dof<double>::Pointer> dofs;
        std::vector<double> gradient;
        double gap;
    } constraint_data;
    KRATOS_TRY
        const std::array<const Variable<double>*,3> all_displacement_components {
            &DISPLACEMENT_X,
            &DISPLACEMENT_Y,
            &DISPLACEMENT_Z
        };

        Condition::DofsVectorType dofs;
        for (Condition& r_condition : r_model_part.Conditions()) {
            r_condition.GetDofList(dofs, r_model_part.GetProcessInfo());
            for (Dof<double>::Pointer p_negative_side_dof : dofs) {
                const auto it_displacement_component = std::find_if(
                    all_displacement_components.begin(),
                    all_displacement_components.end(),
                    [p_negative_side_dof](const Variable<double>* p_variable) {
                        return p_variable->Key() == p_negative_side_dof->GetVariable().Key();
                    });
                if (it_displacement_component != all_displacement_components.end()) {
                    // Insert the negative side component.
                    const std::size_t i_displacement_component = std::distance(all_displacement_components.begin(), it_displacement_component);

                    {
                        const double gradient_component = pretension_surface.normal[i_displacement_component];
                        constraint_data.dofs.push_back(p_negative_side_dof);
                        constraint_data.gradient.push_back(gradient_component);
                    }

                    // Find the duplicate node on the positive side.
                    {
                        const Node::IndexType negative_side_node_id = p_negative_side_dof->Id();
                        const auto it_positive_side_node = duplicated_nodes.find(negative_side_node_id);
                        KRATOS_ERROR_IF(it_positive_side_node == duplicated_nodes.end());
                        Node& r_positive_side_node = *it_positive_side_node->second;
                        const auto& r_positive_side_dofs = r_positive_side_node.GetDofs();
                        const auto it_positive_side_dof = std::find_if(
                            r_positive_side_dofs.begin(),
                            r_positive_side_dofs.end(),
                            [it_displacement_component](const auto& rp_dof) {
                                return rp_dof->GetVariable().Key() == (*it_displacement_component)->Key();
                            });
                        KRATOS_ERROR_IF(it_positive_side_dof == r_positive_side_dofs.end());
                        const double gradient_component = -pretension_surface.normal[i_displacement_component];
                        constraint_data.dofs.push_back(it_positive_side_dof->get());
                        constraint_data.gradient.push_back(gradient_component);
                    }
                } // if p_dof in all_displacement_components
            } // for p_negative_side_dof in dofs
        } // for r_condition in r_model_part.Conditions

        constraint_data.gap = mpImpl->mPretensionValue;
    KRATOS_CATCH("")

    // Insert the constraint.
    const std::size_t constraint_id = r_model_part.MasterSlaveConstraints().empty() ?
                                      1ul :
                                      r_model_part.MasterSlaveConstraints().back().Id() + 1;

    MasterSlaveConstraint::DofPointerVectorType masters, slaves;
    MasterSlaveConstraint::MatrixType relation_matrix;
    MasterSlaveConstraint::VectorType constants;

    slaves.push_back(constraint_data.dofs.back());
    masters.reserve(constraint_data.dofs.size() - 1);
    std::copy(
        constraint_data.dofs.begin(),
        constraint_data.dofs.end() - 1,
        std::back_inserter(masters));

    relation_matrix.resize(1, constraint_data.dofs.size() - 1, false);
    std::transform(
        constraint_data.gradient.begin(),
        constraint_data.gradient.end() - 1,
        relation_matrix.begin2(),
        [&constraint_data](double component){return -component / constraint_data.gradient.back();});

    constants.resize(1, false);
    constants[0] = -constraint_data.gap;

    KRATOS_ERROR_IF_NOT(KratosComponents<MasterSlaveConstraint>::Has("LinearMasterSlaveConstraint"));
    MasterSlaveConstraint::Pointer p_constraint;

    KRATOS_TRY
        p_constraint = KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint").Create(
            constraint_id,
            masters,
            slaves,
            relation_matrix,
            constants);
        r_model_part.AddMasterSlaveConstraint(p_constraint);
    KRATOS_CATCH("")
}


const Parameters InsertPretensionProcess::GetDefaultParameters() const {
    return Parameters(R"({
        "model_part_name" : "",
        "pretension_value" : 0.0,
        "imposition" : "lagrange_element"
    })");
}


} // namespace Kratos
