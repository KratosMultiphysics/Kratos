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
#include "custom_operations/insert_pre_tension_operation.hpp" // InsertPreTensionOperation
#include "custom_conditions/point_load_condition_1d_1n.h" // PointLoadCondition1D1N
#include "structural_mechanics_application_variables.h" // POINT_LOAD_X

// --- Core Includes ---
#include "includes/master_slave_constraint.h" // MasterSlaveConstraint
#include "includes/model_part.h" // ModelPart
#include "processes/find_global_nodal_elemental_neighbours_process.h"
#include "utilities/comparison.h" // Comparison
#include "includes/element.h" // Element
#include "input_output/logger.h" // KRATOS_INFO
#include "solving_strategies/builder_and_solvers/p_multigrid/linear_multifreedom_constraint.hpp" // LinearMultifreedomConstraint

// --- STL Includes ---
#include <array> // std::array
#include <sstream> // std::stringstream
#include <optional> // std::optional
#include <memory> // std::shared_ptr
#include <array> // std::array
#include <vector> // std::vector
#include <algorithm>
#include <numeric>
#include <string>
#include <unordered_set>
#include <unordered_map>


namespace Kratos {


struct InsertPreTensionOperation::Impl {
    ModelPart* mpModelPart;
    int mVerbosity;

    // Setting NEIGHBOUR_ELEMENTS on data value containers is fucked,
    // so geometry-element adjacency has to be stored in another manner.
    std::shared_ptr<std::unordered_map<
        Geometry<Node>::IndexType,
        std::vector<Element::Pointer>>
    > mpAdjacencyMap;
}; // struct InsertPreTensionOperation::Impl


InsertPreTensionOperation::InsertPreTensionOperation(Model& rModel, Parameters Settings)
    : Operation(),
      mpImpl(new Impl) {
    Settings.ValidateAndAssignDefaults(this->GetDefaultParameters());

    KRATOS_TRY
        mpImpl->mpModelPart         = &rModel.GetModelPart(Settings["model_part_name"].GetString());
        mpImpl->mVerbosity          = Settings["verbosity"].GetInt();
    KRATOS_CATCH("")

    // This operation expects the pre-tension surface to be defined
    // by a set of Geometries in the input ModelPart. It requires
    // at least one Geometry.
    KRATOS_ERROR_IF(mpImpl->mpModelPart->Geometries().empty())
        << "InsertPreTensionOperation requires at least one Geometry in the provided ModelPart, "
        << "but there aren't any in '" << mpImpl->mpModelPart->Name() << "'.";

    mpImpl->mpAdjacencyMap = std::make_shared<std::unordered_map<
        Geometry<Node>::IndexType,
        std::vector<Element::Pointer>>>();
}


InsertPreTensionOperation::InsertPreTensionOperation(InsertPreTensionOperation&&) noexcept = default;


InsertPreTensionOperation::~InsertPreTensionOperation() = default;


struct PreTensionSurfacePartition {
    array_1d<double,3> normal;
    std::vector<Element::Pointer> positive_side, negative_side;
}; // struct SurfacePartition


template <bool Prune = true>
void FindElementsAdjacentToGeometries(
    ModelPart& rModelPart,
    std::unordered_map<Geometry<Node>::IndexType,std::vector<Element::Pointer>>& rAdjacencyMap) {
        for (const auto& r_geometry : rModelPart.Geometries())
            rAdjacencyMap.emplace(r_geometry.Id(), std::vector<Element::Pointer> {});

        // Find elements containing geometries.
        KRATOS_TRY
            block_for_each(
                rModelPart.Geometries(),
                [&rAdjacencyMap] (Geometry<Node>& r_geometry) -> void {
                    std::vector<Element::Pointer> adjacent_elements;

                    // Collect all neighbor elements from the geometry's nodes.
                    for (Node& r_node : r_geometry) {
                        auto& r_containing_elements = r_node.GetValue(NEIGHBOUR_ELEMENTS);
                        for (unsigned i_element=0u; i_element<r_containing_elements.size(); ++i_element) {
                            auto& rp_element = *(r_containing_elements.ptr_begin() + i_element);
                            const auto it_element = std::find_if(
                                adjacent_elements.begin(),
                                adjacent_elements.end(),
                                [&rp_element](const auto& rp_other_element) -> bool {
                                    return rp_other_element->Id() == rp_element->Id();
                                });
                            if (it_element == adjacent_elements.end()) {
                                adjacent_elements.emplace_back(Element::Pointer(rp_element.get()));
                            }
                        } // for i_element in r_containing_elements.size()
                    } // for node in geometry

                    // Set neighbor elements on the geometry.
                    rAdjacencyMap[r_geometry.Id()] = std::move(adjacent_elements);
                });
        KRATOS_CATCH("")
}


template <unsigned SurfaceDimension>
PreTensionSurfacePartition PartitionPreTensionSurface(
    const ModelPart::GeometryContainerType::iterator itGeometryBegin,
    const ModelPart::GeometryContainerType::iterator itGeometryEnd,
    const std::unordered_map<
        Geometry<Node>::IndexType,
        std::vector<Element::Pointer>
    >& rAdjacencyMap) {
        PreTensionSurfacePartition surface_partition;
        std::fill(
            surface_partition.normal.begin(),
            surface_partition.normal.end(),
            0.0);

        // There are three valid possibilities here, separately for split and unsplit cases.
        // 1) All surfaces are 0-dimensional (corners)
        //    => expecting exactly 1 geometry
        //    => all connected elements must have line geometries
        //    => the normal is not well-defined, so it is chosen
        //       as the average of the connected elements' directions
        // 2) All surfaces are 1-dimensional (edges)
        //    => expecting at least 1 geometry
        //    => all connected elements must be 2-dimensional
        //    => the normal is not well-defined, so it is chosen as the
        //       average of the rotated endpoints in 2D, in positive
        //       direction. This assumes coordinates in z-direction vanish.
        //       This assumption is checked.
        // 3) All surfaces are 2-dimensional (surfaces)
        //    => expecting at least 1 geometry
        //    => normals are well-defined
        //    => all connected elements must be 3-dimensional
        //
        // In all unsplit cases, at least one connected element is required on
        // either side of the average pre-tension plane.

        // Sanity checks.
        if constexpr (SurfaceDimension == 0u) {
            KRATOS_ERROR_IF_NOT(std::distance(itGeometryBegin, itGeometryEnd) == 1)
                << "Expecting exactly 1 geometry on the 0-dimensional pre-tension surface , but got "
                << std::distance(itGeometryBegin, itGeometryEnd) << ".";
        } else if constexpr (SurfaceDimension == 1u) {
            KRATOS_ERROR_IF_NOT(1 <= std::distance(itGeometryBegin, itGeometryEnd))
                << "Expecting at least 1 geometry on the 1-dimensional pre-tension surface , but got "
                << std::distance(itGeometryBegin, itGeometryEnd) << ".";
        } else if constexpr (SurfaceDimension == 2u) {
            KRATOS_ERROR_IF_NOT(1 <= std::distance(itGeometryBegin, itGeometryEnd))
                << "Expecting at least 1 geometry on the 2-dimensional pre-tension surface , but got "
                << std::distance(itGeometryBegin, itGeometryEnd) << ".";
        } else {
            static_assert(float(SurfaceDimension) == 0.5f, "Invalid surface dimension");
        }

        for (auto it_geometry=itGeometryBegin; it_geometry!=itGeometryEnd; ++it_geometry) {
            Geometry<Node>& r_geometry = *it_geometry;

            // Sanity checks on the geometry.
            KRATOS_ERROR_IF_NOT(r_geometry.LocalSpaceDimension() == SurfaceDimension)
                << "Expecting all geometries defining the pre-tension surface to be " << SurfaceDimension << " dimensional, "
                << "but geometry " << r_geometry.Id() << " on geometry " << r_geometry.Id() << " (" << r_geometry.Name() << ") "
                << "is " << r_geometry.LocalSpaceDimension() << "-dimensional.";

            const auto it_adjacency = rAdjacencyMap.find(r_geometry.Id());
            KRATOS_ERROR_IF(it_adjacency == rAdjacencyMap.end())
                << "no adjacency map is stored for geometry " << r_geometry.Id();

            const std::vector<Element::Pointer>& r_adjacent_elements = it_adjacency->second;

            if constexpr (SurfaceDimension == 0u) {
                for (unsigned i_element=0u; i_element<r_adjacent_elements.size(); ++i_element) {
                    const Element& r_element = *r_adjacent_elements[i_element];
                    KRATOS_ERROR_IF_NOT(r_element.GetGeometry().GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Linear)
                        << "Expecting only " << SurfaceDimension + 1 << "-dimensional elements connected to the " << SurfaceDimension << "-dimensional pre-tension surface, "
                        << "but element " << r_element.Id() << " has a geometry of type " << r_element.GetGeometry().Name() << ".";

                    const std::size_t surface_node_id = itGeometryBegin->operator[](0).Id();
                    const array_1d<double,3> surface_node_position = itGeometryBegin->operator[](0);
                    array_1d<double,3> opposite_node_position;
                    if (r_element.GetGeometry()[0].Id() == surface_node_id) {
                        opposite_node_position = r_element.GetGeometry()[1];
                    } else if (r_element.GetGeometry()[1].Id() == surface_node_id) {
                        opposite_node_position = r_element.GetGeometry()[0];
                    } else {
                        KRATOS_ERROR
                            << "Element " << r_element.Id() << " is assumed to be connected to node "
                            << surface_node_id << " but is not.";
                    }

                    if (i_element % 2) {
                        surface_partition.positive_side.push_back(*(r_adjacent_elements.begin() + i_element));
                        surface_partition.normal = opposite_node_position - surface_node_position;
                    } else {
                        surface_partition.negative_side.push_back(*(r_adjacent_elements.begin() + i_element));
                        surface_partition.normal = surface_node_position - opposite_node_position;
                    }
                } // for i_element in range(r_adjacent_elements.size())
            } /*if SurfaceDimension == 0*/ else if constexpr (SurfaceDimension == 1u) {
                // Define the line geometry's plane.
                const array_1d<double,3>
                    surface_begin = r_geometry[0],
                    surface_end   = *(r_geometry.end() - 1);
                KRATOS_ERROR_IF_NOT(Comparison<double>::Equal(1e-12, 1e-9)(surface_begin[2], 0.0))
                    << "Assuming that 1-dimensional pre-tension surfaces lie on the XY plane, but node "
                    << r_geometry[0].Id() << " of geometry " << r_geometry.Id() << " "
                    << "has a Z-coordinate of " << surface_begin[2] << ".";
                KRATOS_ERROR_IF_NOT(Comparison<double>::Equal(1e-12, 1e-9)(surface_end[2], 0.0))
                    << "Assuming that 1-dimensional pre-tension surfaces lie on the XY plane, but node "
                    << r_geometry[1].Id() << " of geometry " << r_geometry.Id() << " "
                    << "has a Z-coordinate of " << surface_end[2] << ".";
                const array_1d<double,3>
                    surface_normal {
                        -(surface_end[1] - surface_begin[1]),
                        surface_end[0] - surface_begin[0],
                        0.0},
                    surface_center = r_geometry.Center();
                surface_partition.normal += surface_normal;

                // Sort connected elements.
                for (unsigned i_element=0u; i_element<r_adjacent_elements.size(); ++i_element) {
                    const Element& r_element = *r_adjacent_elements[i_element];
                    const array_1d<double,3> element_direction = r_element.GetGeometry().Center() - surface_center;
                    const double direction_inner_product = std::inner_product(
                        surface_normal.begin(),
                        surface_normal.end(),
                        element_direction.begin(),
                        0.0);
                    if (0 < direction_inner_product) surface_partition.negative_side.push_back(*(r_adjacent_elements.begin() + i_element));
                    else surface_partition.positive_side.push_back(*(r_adjacent_elements.begin() + i_element));
                } // for i_element in range(r_adjacent_elements.size())
            } /*else if SurfaceDimension == 1*/ else if constexpr (SurfaceDimension == 2u) {
                // Define the surface's plane.
                const array_1d<double,3> surface_center = r_geometry.Center();
                const array_1d<double,3> surface_normal = r_geometry.Normal(surface_center);
                surface_partition.normal += surface_normal;

                // Sort connected elements.
                for (unsigned i_element=0u; i_element<r_adjacent_elements.size(); ++i_element) {
                    const Element& r_element = *r_adjacent_elements[i_element];
                    const array_1d<double,3> element_direction = r_element.GetGeometry().Center() - surface_center;
                    const double direction_inner_product = std::inner_product(
                        surface_normal.begin(),
                        surface_normal.end(),
                        element_direction.begin(),
                        0.0);
                    if (0 < direction_inner_product) surface_partition.negative_side.push_back(*(r_adjacent_elements.begin() + i_element));
                    else surface_partition.positive_side.push_back(*(r_adjacent_elements.begin() + i_element));
                } // for i_element in range(r_adjacent_elements.size())
            } /*else if SurfaceDimension == 2*/
        } // for r_condition in r_model_part.Conditions

        // Normalize the surface's normal.
        const double norm = std::sqrt(std::accumulate(
            surface_partition.normal.begin(),
            surface_partition.normal.end(),
            0.0,
            [] (double sum, double component) -> double {
                return sum + component * component;
            }));
        KRATOS_ERROR_IF_NOT(norm) << "normal vanished on the pretension surface";
        std::transform(
            surface_partition.normal.begin(),
            surface_partition.normal.end(),
            surface_partition.normal.begin(),
            [norm] (double component) -> double {
                return component / norm;
            });

        return surface_partition;
}


/// @brief Base class for @ref InPlaneRelativeDisplacementConstraint and @ref OutOfPlaneDisplacementConstraint.
class PlaneDisplacementConstraint : public MultifreedomConstraint {
public:
    using ProtectedNormal = std::pair<
        std::optional<array_1d<double,3>>,
        LockObject>;

    PlaneDisplacementConstraint() noexcept = default;

    PlaneDisplacementConstraint(
        MasterSlaveConstraint::IndexType Id,
        const ModelPart::GeometryContainerType& rSurfaceGeometries,
        DofPointerVectorType&& rDofs,
        const std::vector<std::size_t>& rConstraintLabels,
        std::shared_ptr<ProtectedNormal> pSharedSurfaceNormal,
        std::shared_ptr<const std::unordered_map<Geometry<Node>::IndexType,std::vector<Element::Pointer>>> pAdjacencyMap,
        int Verbosity)
            :   MultifreedomConstraint(
                    Id,
                    std::move(rDofs),
                    rConstraintLabels),
                mVerbosity(Verbosity),
                mSurfaceGeometries(rSurfaceGeometries),
                mpAdjacencyMap(pAdjacencyMap),
                mpSharedSurfaceNormal(pSharedSurfaceNormal) {
        KRATOS_ERROR_IF(rSurfaceGeometries.empty());
    }

    void InitializeSolutionStep(const ProcessInfo& rProcessInfo) override {
        this->InitializeNonLinearIteration(rProcessInfo);
    }

    /// @brief Linearize the constraint equation.
    void InitializeNonLinearIteration(const ProcessInfo&) override {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFileName() << " is virtual";
    }

    void FinalizeNonLinearIteration(const ProcessInfo&) override {
        KRATOS_TRY
            std::scoped_lock<LockObject> lock(mpSharedSurfaceNormal->second);
            mpSharedSurfaceNormal->first.reset();
        KRATOS_CATCH("")
    }

    void FinalizeSolutionStep(const ProcessInfo& rProcessInfo) override {
        this->FinalizeNonLinearIteration(rProcessInfo);
    }

    void CalculateLocalSystem(
        MatrixType& rConstraintGradient,
        VectorType& rConstraintGaps,
        const ProcessInfo&) const override {
            rConstraintGradient = this->GetData().GetValue(CONSTITUTIVE_MATRIX);
            rConstraintGaps = this->GetData().GetValue(INTERNAL_FORCES_VECTOR);
    }

    std::string Info() const override {
        return "PlaneDisplacementConstraint";
    }

    void save(Serializer& rSerializer) const override {
        KRATOS_TRY
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, MultifreedomConstraint);
            rSerializer.save("SurfaceGeometries", mSurfaceGeometries);
            rSerializer.save("Verbosity", mVerbosity);
        KRATOS_CATCH("")
    }

    void load(Serializer& rDeserializer) override {
        KRATOS_TRY
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rDeserializer, MultifreedomConstraint);
            rDeserializer.load("SurfaceGeometries", mSurfaceGeometries);
            rDeserializer.load("Verbosity", mVerbosity);
        KRATOS_CATCH("")
    }

protected:
    array_1d<double,3> GetSurfaceNormal() {
        KRATOS_TRY
            std::scoped_lock<LockObject> lock(mpSharedSurfaceNormal->second);
            if (mpSharedSurfaceNormal->first.has_value()) {
                return mpSharedSurfaceNormal->first.value();
            } else {
                PreTensionSurfacePartition surface_partition;
                switch (mSurfaceGeometries.front().LocalSpaceDimension()) {
                    case 0: {
                        surface_partition = PartitionPreTensionSurface<0>(
                            mSurfaceGeometries.begin(),
                            mSurfaceGeometries.end(),
                            *mpAdjacencyMap);
                        break;
                    }
                    case 1: {
                        surface_partition = PartitionPreTensionSurface<1>(
                            mSurfaceGeometries.begin(),
                            mSurfaceGeometries.end(),
                            *mpAdjacencyMap);
                        break;
                    }
                    case 2: {
                        surface_partition = PartitionPreTensionSurface<2>(
                            mSurfaceGeometries.begin(),
                            mSurfaceGeometries.end(),
                            *mpAdjacencyMap);
                        break;
                    }
                    default: KRATOS_ERROR << "Invalid pre-tension surface dimension: " << mSurfaceGeometries.front().LocalSpaceDimension() << ".";
                } // switch dimension

                KRATOS_INFO_IF(this->Info() + " " + std::to_string(this->Id()), 3 <= mVerbosity)
                    << "pre-tension surface normal: ["
                    << surface_partition.normal[0] << ' '
                    << surface_partition.normal[1] << ' '
                    << surface_partition.normal[2] << "]\n";
                mpSharedSurfaceNormal->first = surface_partition.normal;
                return surface_partition.normal;
            }
        KRATOS_CATCH("")
    }

    int mVerbosity;

    ModelPart::GeometryContainerType mSurfaceGeometries;

    std::shared_ptr<const std::unordered_map<
        Geometry<Node>::IndexType,
        std::vector<Element::Pointer>>
    > mpAdjacencyMap;

    /// @details The surface's normal is potentially shared by a set of constraints,
    ///          since each constraint only restricts one pair of nodes. If the surface
    ///          contains more than a single pair, computing its normal should not be
    ///          unnecessarily repeated.
    std::shared_ptr<std::pair<
        std::optional<array_1d<double,3>>,
        LockObject
    >> mpSharedSurfaceNormal;
}; // class PlaneDisplacementConstraint


/// @brief Forbid the relative in-plane displacement of a pair of nodes.
/// @details The surface's normal is computed in each non-linear iteration
///          as the weighted average of a set of geometries' normals.
class InPlaneRelativeDisplacementConstraint : public PlaneDisplacementConstraint {
public:
    using PlaneDisplacementConstraint::PlaneDisplacementConstraint;

    MasterSlaveConstraint::Pointer Clone(IndexType Id) const override {
        const auto& r_constraint_labels = this->GetValue(CONSTRAINT_LABELS);
        std::vector<std::size_t> constraint_labels(
            r_constraint_labels.begin(),
            r_constraint_labels.end());
        return MasterSlaveConstraint::Pointer(new InPlaneRelativeDisplacementConstraint(
            Id,
            mSurfaceGeometries,
            DofPointerVectorType(this->GetDofs()),
            constraint_labels,
            mpSharedSurfaceNormal,
            mpAdjacencyMap,
            mVerbosity));
    }

    void InitializeNonLinearIteration(const ProcessInfo&) override {
        KRATOS_TRY
            const array_1d<double,3> normal = this->GetSurfaceNormal();

            // DoFs are assumed to be stored in a specific layout.
            // - contained DoFs are related to two (an original and its duplicate) nodes
            // - the same set (n) of variables are referenced
            // - the first half of the array contains DISPLACEMENT_X, DISPLACEMENT_Y ... of the original node in that order
            // - the second half is identical, but for the duplicate node
            auto& r_dofs = this->GetDofs();
            KRATOS_ERROR_IF(r_dofs.size() % 2);
            const std::size_t dimension_count = r_dofs.size() / 2;

            Matrix constraint_gradients(dimension_count, 2 * dimension_count);
            Vector constraint_gaps, displacements(2 * dimension_count);

            // Build the constraint equations.
            for (std::size_t i_dimension=0ul; i_dimension<dimension_count; ++i_dimension) {
                for (std::size_t j_dimension=0ul; j_dimension<dimension_count; ++j_dimension) {
                    double gradient_component = -normal[i_dimension] * normal[j_dimension];
                    if (i_dimension == j_dimension) gradient_component += 1.0;
                    constraint_gradients(i_dimension, j_dimension) = gradient_component;
                    constraint_gradients(i_dimension, j_dimension + dimension_count) = -gradient_component;
                } // for j_dimension in range(dimension_count)
                displacements[i_dimension] = r_dofs[i_dimension]->GetSolutionStepValue();
                displacements[i_dimension + dimension_count] = r_dofs[i_dimension + dimension_count]->GetSolutionStepValue();
            } // for i_dimension in range(dimension_count)

            constraint_gaps = prod(constraint_gradients, displacements);

            // Store constraint gradients and constraint gaps.
            this->SetValue(CONSTITUTIVE_MATRIX, constraint_gradients);
            this->SetValue(INTERNAL_FORCES_VECTOR, constraint_gaps);
        KRATOS_CATCH("")
    }

    std::string Info() const override {
        return "InPlaneRelativeDisplacementConstraint";
    }
}; // class InPlaneRelativeDisplacementConstraint


Dof<double>* FindDofInNode(Node& rNode, const Variable<double>& rVariable) noexcept {
    auto& r_dofs = rNode.GetDofs();
    const auto it_dof = std::find_if(
        r_dofs.begin(),
        r_dofs.end(),
        [&rVariable] (const auto& rp_dof) -> bool {
            return rp_dof->GetVariable() == rVariable;
        });
    return it_dof == r_dofs.end()
        ? nullptr
        : it_dof->get();
}


void InsertPreTensionOperation::Execute() {
    ModelPart& r_model_part = *mpImpl->mpModelPart;

    // Find elements connected to the pre-tension surface.
    // => each node will have a NEIGHBOUR_ELEMENTS variable
    //    in its non-historical storage that lists elements
    //    containing them.
    KRATOS_TRY
        FindGlobalNodalElementalNeighboursProcess(r_model_part.GetRootModelPart()).Execute();
        FindElementsAdjacentToGeometries<false>(r_model_part, *mpImpl->mpAdjacencyMap);
    KRATOS_CATCH("")

    // Sort connected elements into 2 categories, lying on either side of
    // the pre-tension surface.
    PreTensionSurfacePartition surface_partition;
    KRATOS_TRY
        // Find the dimension of the pre-tension surface,
        // and make sure all conditions share it.
        KRATOS_ERROR_IF(r_model_part.Geometries().empty());
        const int pre_tension_surface_dimension = r_model_part.Geometries().front().LocalSpaceDimension();
        switch (pre_tension_surface_dimension) {
            case 0u: {
                surface_partition = PartitionPreTensionSurface<0>(
                    r_model_part.Geometries().begin(),
                    r_model_part.Geometries().end(),
                    *mpImpl->mpAdjacencyMap);
                break;
            }
            case 1u: {
                surface_partition = PartitionPreTensionSurface<1>(
                    r_model_part.Geometries().begin(),
                    r_model_part.Geometries().end(),
                    *mpImpl->mpAdjacencyMap);
                break;
            }
            case 2u: {
                surface_partition = PartitionPreTensionSurface<2>(
                    r_model_part.Geometries().begin(),
                    r_model_part.Geometries().end(),
                    *mpImpl->mpAdjacencyMap);
                break;
            }
            default: KRATOS_ERROR << "Invalid pre-tension surface dimension: " << pre_tension_surface_dimension << ".";
        }
    KRATOS_CATCH("")

    if (3 <= mpImpl->mVerbosity) {
        KRATOS_INFO(this->Info())
            << "pre-tension surface normal: ["
            << surface_partition.normal[0] << ' '
            << surface_partition.normal[1] << ' '
            << surface_partition.normal[2] << "]\n";

        std::stringstream message;
        if (3 <= mpImpl->mVerbosity) {
            message << "elements on the negative side: [";
            for (const auto& rp_element : surface_partition.negative_side)
                message << rp_element->Id() << ' ';
            message << "]\n";
            KRATOS_INFO(this->Info()) << message.view();

            message = std::stringstream();
            message << "elements on the positive side: [";
            for (const auto& rp_element : surface_partition.positive_side)
                message << rp_element->Id() << ' ';
            message << "]\n";
            KRATOS_INFO(this->Info()) << message.view();
        } // if 3 <= verbosity
    } // if 3 <= verbosity

    // Insert a new control node acting as a master for all
    // other nodes on the pre-tension surface.
    Node::Pointer p_control_node = mpImpl->mpModelPart->CreateNewNode(
        /*Id=*/mpImpl->mpModelPart->GetRootModelPart().Nodes().empty()
            ? 1
            : mpImpl->mpModelPart->GetRootModelPart().Nodes().back().Id() + 1,
        /*x=*/0.0,
        /*y=*/0.0,
        /*z=*/0.0);
    p_control_node->AddDof(DISPLACEMENT_X);
    KRATOS_INFO_IF(this->Info(), 2 <= mpImpl->mVerbosity)
        << "insert control node " << p_control_node->Id() << ' '
        << "at [" << p_control_node->X() << ' ' << p_control_node->Y() << ' ' << p_control_node->Z() << "]\n";

    // Duplicate nodes on the pre-tension surface.
    std::unordered_map<Node*,Node::Pointer> duplicated_nodes;
    KRATOS_TRY
        std::unordered_set<Node*> surface_nodes;

        // Collect nodes on the surface.
        for (Geometry<Node>& r_geometry : r_model_part.Geometries()) {
            for (Node& r_node : r_geometry) {
                surface_nodes.insert(&r_node);
            } // for r_node in r_geometry
        } // for r_geometry in r_model_part.Conditions()

        // Find the first available node ID that can safely be
        // used to insert new ones.
        Node::IndexType node_id = r_model_part.GetRootModelPart().Nodes().back().Id() + 1;

        // Clone surface nodes and associate them with the original ones.
        for (Node* p_node : surface_nodes) {
            const auto emplace_result = duplicated_nodes.emplace(
                p_node,
                p_node->Clone(node_id++));
            KRATOS_INFO_IF(this->Info(), 3 <= mpImpl->mVerbosity)
                << "duplicate node " << p_node->Id() << " => " << emplace_result.first->second->Id() << "\n";
            r_model_part.AddNode(emplace_result.first->second);
        } // for p_node in surface_nodes

        // Replace nodes on the negative side elements with duplicated ones.
        for (auto& rp_element : surface_partition.negative_side) {
            for (auto itp_node=rp_element->GetGeometry().ptr_begin(); itp_node!=rp_element->GetGeometry().ptr_end(); ++itp_node) {
                const auto it_duplicate_pair = duplicated_nodes.find(&**itp_node);
                if (it_duplicate_pair != duplicated_nodes.end()) {
                    (*itp_node) = it_duplicate_pair->second;
                }
            } // for itp_node in rp_element->GetGeometry()
        } // for rp_element in surface_partition.negative_side
    KRATOS_CATCH("")

    // Prune elements that don't have at least one face on the pre-tension surface.
    KRATOS_TRY
        for (auto& [id_geometry, r_adjacent_elements] : *mpImpl->mpAdjacencyMap) {
            const Geometry<Node>& r_geometry = r_model_part.GetGeometry(id_geometry);
            r_adjacent_elements.erase(std::remove_if(
                r_adjacent_elements.begin(),
                r_adjacent_elements.end(),
                [&r_geometry] (const Element::Pointer& rp_element) -> bool {
                    for (const Node& r_geometry_node : r_geometry) {
                        const auto it_node = std::find_if(
                            rp_element->GetGeometry().begin(),
                            rp_element->GetGeometry().end(),
                            [&r_geometry_node] (const Node& r_element_node) -> bool {
                                return r_element_node.Id() == r_geometry_node.Id();
                            });
                        if (it_node == rp_element->GetGeometry().end())
                            return true;
                    } // for r_geometry_node in r_geometry
                    return false;
                }
            ), r_adjacent_elements.end());
        } // for id_geometry, r_adjacency_elements
    KRATOS_CATCH("")

    KRATOS_TRY
        FindGlobalNodalElementalNeighboursProcess(r_model_part.GetRootModelPart()).Execute();
    KRATOS_CATCH("")

    // Collect DoFs from elements on the positive side of the pre-tension surface.
    // This is necessary because elements on the pre-tension surface may not have
    // DoFs for all displacement components. In this case, constraints should not
    // be inserted for these components.
    std::unordered_set<const Dof<double>*> positive_side_dofs;
    {
        Element::DofsVectorType dof_buffer;
        for (const auto& rp_element : surface_partition.positive_side) {
            rp_element->GetDofList(dof_buffer, mpImpl->mpModelPart->GetProcessInfo());
            positive_side_dofs.insert(
                dof_buffer.begin(),
                dof_buffer.end());
        }
    }

    LinearMultifreedomConstraint::IndexType id_constraint = mpImpl->mpModelPart->GetRootModelPart().MasterSlaveConstraints().empty()
        ? 1
        : mpImpl->mpModelPart->GetRootModelPart().MasterSlaveConstraints().back().Id() + 1;

    // Forbid relative in-plane displacements for each duplicated node pair.
    {
        const std::array<const Variable<double>*,3> all_displacement_components {
            &DISPLACEMENT_X,
            &DISPLACEMENT_Y,
            &DISPLACEMENT_Z};

        std::vector<std::array<Dof<double>*,2>> dof_pairs;
        auto p_protected_surface_normal = std::make_shared<std::pair<
            std::optional<array_1d<double,3>>,
            LockObject>>();

        for (auto& [rp_positive_side_node, rp_negative_side_node] : duplicated_nodes) {
            dof_pairs.clear();

            // Collect relevant DoF pairs.
            // dof_pairs will contain [{u_x_positive, u_x_negative}, {u_y_positive, u_y_negative}, ...]
            for (const auto& rp_variable : all_displacement_components) {
                Dof<double>* p_positive_side_dof = FindDofInNode(*rp_positive_side_node, *rp_variable);
                if (positive_side_dofs.find(p_positive_side_dof) == positive_side_dofs.end()) continue;

                Dof<double>* p_negative_side_dof = FindDofInNode(*rp_negative_side_node, *rp_variable);
                if (p_positive_side_dof && p_negative_side_dof)
                    dof_pairs.push_back({p_positive_side_dof, p_negative_side_dof});
            } // for rp_variable in all_displacement_components

            MasterSlaveConstraint::DofPointerVectorType dofs(2 * dof_pairs.size());
            for (std::size_t i_pair=0ul; i_pair<dof_pairs.size(); ++i_pair) {
                dofs[i_pair] = dof_pairs[i_pair].front();
                dofs[i_pair + dof_pairs.size()] = dof_pairs[i_pair].back();
            }

            std::vector<std::size_t> constraint_labels(dof_pairs.size());
            std::iota(
                constraint_labels.begin(),
                constraint_labels.end(),
                id_constraint);
            MasterSlaveConstraint::Pointer p_constraint(new InPlaneRelativeDisplacementConstraint(
                id_constraint,
                mpImpl->mpModelPart->Geometries(),
                std::move(dofs),
                constraint_labels,
                p_protected_surface_normal,
                mpImpl->mpAdjacencyMap,
                mpImpl->mVerbosity));

            KRATOS_INFO_IF(this->Info(), 3 <= mpImpl->mVerbosity)
                << "insert constraint " << p_constraint->Id() << ' '
                << "prohibiting the in-plane relative displacement of nodes "
                << rp_positive_side_node->Id() << " and " << rp_negative_side_node->Id() << '\n';
            mpImpl->mpModelPart->AddMasterSlaveConstraint(p_constraint);

            id_constraint += constraint_labels.size();
        } // for id_positive_side_node, rp_negative_side_node : duplicated_nodes
    }

    // Tie rotation DoFs between duplicate nodes.
    {
        const std::array<const Variable<double>*,3> all_rotation_components {
            &ROTATION_X,
            &ROTATION_Y,
            &ROTATION_Z};

        for (auto& [rp_positive_side_node, rp_negative_side_node] : duplicated_nodes) {
            for (const auto& rp_variable : all_rotation_components) {
                Dof<double>* p_positive_side_dof = FindDofInNode(*rp_positive_side_node, *rp_variable);
                if (positive_side_dofs.find(p_positive_side_dof) == positive_side_dofs.end()) continue;
                Dof<double>* p_negative_side_dof = FindDofInNode(*rp_negative_side_node, *rp_variable);
                if (p_positive_side_dof && p_negative_side_dof) {
                    LinearMultifreedomConstraint::DofPointerVectorType dofs(2);
                    LinearMultifreedomConstraint::MatrixType constraint_gradient(1, 2);
                    LinearMultifreedomConstraint::VectorType constraint_gap(1);
                    dofs[0] = p_positive_side_dof;
                    dofs[1] = p_negative_side_dof;
                    constraint_gradient(0, 0) =  1.0;
                    constraint_gradient(0, 1) = -1.0;
                    constraint_gap[0] = 0.0;
                    MasterSlaveConstraint::Pointer p_constraint(new LinearMultifreedomConstraint(
                        id_constraint,
                        std::move(dofs),
                        {id_constraint},
                        constraint_gradient,
                        constraint_gap));
                    KRATOS_INFO_IF(this->Info(), 3 <= mpImpl->mVerbosity)
                        << "tie " << p_positive_side_dof->GetVariable().Name() << ' '
                        << "of node " << rp_positive_side_node->Id() << ' '
                        << "to " << p_negative_side_dof->GetVariable().Name() << ' '
                        << "of node " << rp_negative_side_node->Id() << '\n';
                    mpImpl->mpModelPart->AddMasterSlaveConstraint(p_constraint);
                    ++id_constraint;
                } // if p_positive_side_dof && p_negative_side_dof
            } // for rp_variable in all_rotation_components
        } // for rp_positive_side_node, rp_negative_side_node
    }

    // Decide what to do with the control node in derived classes.
    this->InsertControlNodeConstraints(
        *mpImpl->mpModelPart,
        surface_partition.normal,
        duplicated_nodes,
        p_control_node,
        positive_side_dofs);
}


const Parameters InsertPreTensionOperation::GetDefaultParameters() const {
    return Parameters(R"({
        "model_part_name" : "",
        "verbosity" : 1
    })");
}


std::string InsertPreTensionOperation::Info() const {
    return "InsertPreTensionOperation";
}


/// @brief Tie the average [relative] out-of-plane displacement to a common control DoF.
template <bool IsRelative>
class OutOfPlaneDisplacementConstraint : public PlaneDisplacementConstraint {
public:
    using PlaneDisplacementConstraint::PlaneDisplacementConstraint;

    MasterSlaveConstraint::Pointer Clone(IndexType Id) const override {
        const auto& r_constraint_labels = this->GetValue(CONSTRAINT_LABELS);
        std::vector<std::size_t> constraint_labels(
            r_constraint_labels.begin(),
            r_constraint_labels.end());
        return MasterSlaveConstraint::Pointer(new OutOfPlaneDisplacementConstraint(
            Id,
            mSurfaceGeometries,
            DofPointerVectorType(this->GetDofs()),
            constraint_labels,
            mpSharedSurfaceNormal,
            mpAdjacencyMap,
            mVerbosity));
    }

    void InitializeNonLinearIteration(const ProcessInfo&) override {
        KRATOS_TRY
            // Compute the plane's normal.
            const array_1d<double,3> normal = this->GetSurfaceNormal();

            // DoFs are assumed in a specific layout.
            // - DISPLACEMENT_X of node 0
            // - DISPLACEMENT_X of node 0's pair (only if relative)
            // - DISPLACEMENT_Y of node 0
            // - DISPLACEMENT_Y of node 0's pair (only if relative)
            // ...
            // - DISPLACEMENT_Z (depending on the plane's dimension) of node n
            // - DISPLACEMENT_Z (depending on the plane's dimension) of node n's pair (only if relative)
            // - DoF of the control node
            auto& r_dofs = this->GetDofs();
            KRATOS_ERROR_IF(r_dofs.empty());

            const std::size_t controlled_dof_count = r_dofs.size() - 1;
            KRATOS_ERROR_IF(IsRelative && controlled_dof_count % 2);

            const std::size_t positive_side_dof_count = IsRelative
                ? controlled_dof_count / 2
                : controlled_dof_count;

            // Find the plane's dimension by checking how many unique variables the DoFs reference.
            std::size_t dimension_count = 0ul;
            {
                std::unordered_set<std::size_t> keys;
                for (std::size_t i_dof=0ul; i_dof<r_dofs.size()-1; ++i_dof) {
                    keys.insert(r_dofs[i_dof]->GetVariable().Key());
                } // for i_dof in r_dofs.size()-1
                dimension_count = keys.size();
                KRATOS_ERROR_IF_NOT(0 < dimension_count && dimension_count < 4)
                    << "unexpected surface dimension " << dimension_count;
            }

            // Build the linearized constraint equation.
            Matrix constraint_gradient(1, r_dofs.size());
            Vector constraint_gap(1), displacements(r_dofs.size());

            constexpr std::size_t dof_stride = IsRelative ? 2ul : 1ul;
            for (std::size_t i_dof=0ul; i_dof<positive_side_dof_count; ++i_dof) {
                constraint_gradient(0, dof_stride * i_dof) = normal[i_dof % dimension_count];
                if constexpr (IsRelative)
                    constraint_gradient(0, dof_stride * i_dof + 1) = -normal[i_dof % dimension_count];
            }

            constraint_gap[0] = 0.0;

            // Add the control DoF's coefficient.
            constraint_gradient(0, controlled_dof_count) = -1.0;

            // Adjust for current the current state.
            std::transform(
                r_dofs.begin(),
                r_dofs.end(),
                displacements.begin(),
                [] (const Dof<double>* p_dof) -> double {
                    return p_dof->GetSolutionStepValue();
                });

            constraint_gap += prod(constraint_gradient, displacements);

            // Register the linearized constraint equation.
            this->SetValue(CONSTITUTIVE_MATRIX, constraint_gradient);
            this->SetValue(INTERNAL_FORCES_VECTOR, constraint_gap);
        KRATOS_CATCH("")
    }

    std::string Info() const override {
        if constexpr (IsRelative)
            return "OutOfPlaneRelativeDisplacementConstraint";
        else
            return "OutOfPlaneDisplacementConstraint";
    }
}; // class OutOfPlaneDisplacementConstraint


void InsertDirichletPreTensionOperation::InsertControlNodeConstraints(
    ModelPart& rModelPart,
    array_1d<double,3> SurfaceNormal,
    const std::unordered_map<Node*,Node::Pointer> rDuplicateNodeMap,
    Node::Pointer pControlNode,
    const std::unordered_set<const Dof<double>*> rPositiveSideDofs) {
        mpControlNode = pControlNode;
        mPreTensionSurfaceSize = rDuplicateNodeMap.size();

        KRATOS_TRY
            // Insert a constraint that ties the average out-of-plane relative displacement
            // to a prescribed value.
            const std::array<const Variable<double>*,3> all_displacement_components {
                &DISPLACEMENT_X,
                &DISPLACEMENT_Y,
                &DISPLACEMENT_Z};
            LinearMultifreedomConstraint::DofPointerVectorType dofs;

            for (const auto& [rp_positive_side_node, rp_negative_side_node] : rDuplicateNodeMap) {
                for (std::size_t i_variable=0ul; i_variable<all_displacement_components.size(); ++i_variable) {
                    const Variable<double>& r_variable = *all_displacement_components[i_variable];
                    Dof<double>* p_positive_side_dof = FindDofInNode(*rp_positive_side_node, r_variable);
                    if (rPositiveSideDofs.find(p_positive_side_dof) == rPositiveSideDofs.end()) continue;

                    Dof<double>* p_negative_side_dof = FindDofInNode(*rp_negative_side_node, r_variable);
                    if (p_positive_side_dof && p_negative_side_dof) {
                        dofs.push_back(p_positive_side_dof);
                        dofs.push_back(p_negative_side_dof);
                    } // if p_positive_side_dof && p_negative_side_dof
                } // for rp_variable in all_displacement_components
            } // for rp_positive_side_node, rp_negative_side_node

            // Add the control node's DoF to the constraint equation.
            Dof<double>* p_control_dof = mpControlNode->GetDofs()[0].get();
            dofs.push_back(p_control_dof);

            // Construct and register the constraint equation.
            const LinearMultifreedomConstraint::IndexType constraint_id = mpImpl->mpModelPart->GetRootModelPart().MasterSlaveConstraints().empty()
                ? 1
                : mpImpl->mpModelPart->GetRootModelPart().MasterSlaveConstraints().back().Id() + 1;
            auto p_protected_surface_normal = std::make_shared<std::pair<
                std::optional<array_1d<double,3>>,
                LockObject>>();
            MasterSlaveConstraint::Pointer p_constraint(new OutOfPlaneDisplacementConstraint<true>(
                constraint_id,
                mpImpl->mpModelPart->Geometries(),
                std::move(dofs),
                {constraint_id},
                p_protected_surface_normal,
                mpImpl->mpAdjacencyMap,
                mpImpl->mVerbosity));
            mpImpl->mpModelPart->AddMasterSlaveConstraint(p_constraint);
KRATOS_CATCH("")
}


void InsertDirichletPreTensionOperation::Apply(double Magnitude) {
    KRATOS_ERROR_IF_NOT(mpControlNode);
    KRATOS_ERROR_IF_NOT(mpControlNode->GetDofs().size() == 1);
    KRATOS_TRY
        Dof<double>& r_control_dof = *mpControlNode->GetDofs()[0].get();
        r_control_dof.GetSolutionStepValue() = 0.5 * mPreTensionSurfaceSize * Magnitude - r_control_dof.GetSolutionStepValue();
        r_control_dof.FixDof();
    KRATOS_CATCH("")
}


std::string InsertDirichletPreTensionOperation::Info() const {
    return "InsertDirichletPreTensionOperation";
}


void InsertNeumannPreTensionOperation::InsertControlNodeConstraints(
    ModelPart& rModelPart,
    array_1d<double,3> SurfaceNormal,
    const std::unordered_map<Node*,Node::Pointer> rDuplicateNodeMap,
    Node::Pointer pPositiveSideControlNode,
    const std::unordered_set<const Dof<double>*> rPositiveSideDofs) {
        KRATOS_TRY
            const std::array<const Variable<double>*,3> all_displacement_components {
                &DISPLACEMENT_X,
                &DISPLACEMENT_Y,
                &DISPLACEMENT_Z};
            LinearMultifreedomConstraint::DofPointerVectorType positive_side_dofs, negative_side_dofs;

            for (const auto& [rp_positive_side_node, rp_negative_side_node] : rDuplicateNodeMap) {
                for (std::size_t i_variable=0ul; i_variable<all_displacement_components.size(); ++i_variable) {
                    const Variable<double>& r_variable = *all_displacement_components[i_variable];
                    Dof<double>* p_positive_side_dof = FindDofInNode(*rp_positive_side_node, r_variable);
                    if (rPositiveSideDofs.find(p_positive_side_dof) == rPositiveSideDofs.end()) continue;

                    Dof<double>* p_negative_side_dof = FindDofInNode(*rp_negative_side_node, r_variable);
                    if (p_positive_side_dof && p_negative_side_dof) {
                        positive_side_dofs.push_back(p_positive_side_dof);
                        negative_side_dofs.push_back(p_negative_side_dof);
                    } // if p_positive_side_dof && p_negative_side_dof
                } // for rp_variable in all_displacement_components
            } // for rp_positive_side_node, rp_negative_side_node

            // Insert a new node acting as the control node on the negative side.
            Node::Pointer p_negative_side_control_node = pPositiveSideControlNode->Clone(
                mpImpl->mpModelPart->Nodes().empty()
                    ? 1
                    : mpImpl->mpModelPart->Nodes().back().Id() + 1);
            KRATOS_INFO_IF(this->Info(), 2 <= mpImpl->mVerbosity)
                << "insert negative side control node " << p_negative_side_control_node->Id() << ' '
                << "at [" << p_negative_side_control_node->X() << ' '
                << p_negative_side_control_node->Y() << ' '
                << p_negative_side_control_node->Z() << "]\n";
            mpImpl->mpModelPart->AddNode(p_negative_side_control_node);

            // Add control DoFs.
            positive_side_dofs.push_back(pPositiveSideControlNode->GetDofs()[0].get());
            negative_side_dofs.push_back(p_negative_side_control_node->GetDofs()[0].get());

            LinearMultifreedomConstraint::IndexType constraint_id = mpImpl->mpModelPart->GetRootModelPart().MasterSlaveConstraints().empty()
                ? 1
                : mpImpl->mpModelPart->GetRootModelPart().MasterSlaveConstraints().back().Id() + 1;
            Condition::IndexType condition_id = mpImpl->mpModelPart->GetRootModelPart().Conditions().empty()
                ? 1
                : mpImpl->mpModelPart->GetRootModelPart().Conditions().back().Id() + 1;
            auto p_protected_surface_normal = std::make_shared<std::pair<
                    std::optional<array_1d<double,3>>,
                    LockObject>>();

            Properties::Pointer p_properties = mpImpl->mpModelPart->PropertiesArray().empty()
                ? Properties::Pointer(new Properties(1))
                : *mpImpl->mpModelPart->PropertiesArray().begin();
            // Insert the positive side constraint and condition.
            {
                MasterSlaveConstraint::Pointer p_constraint(new OutOfPlaneDisplacementConstraint<false>(
                    constraint_id,
                    mpImpl->mpModelPart->Geometries(),
                    std::move(positive_side_dofs),
                    {constraint_id},
                    p_protected_surface_normal,
                    mpImpl->mpAdjacencyMap,
                    mpImpl->mVerbosity));
                mpPositiveSideLoad = Condition::Pointer(new PointLoadCondition1D1N(
                    condition_id,
                    Geometry<Node>::Pointer(new Point2D<Node>(pPositiveSideControlNode)),
                    p_properties));

                KRATOS_INFO_IF(this->Info(), 3 <= mpImpl->mVerbosity)
                    << "insert constraint " << p_constraint->Id() <<' '
                    << "tying the average out-of-plane displacement on the positive side to "
                    << pPositiveSideControlNode->GetDofs()[0]->GetVariable().Name() << ' '
                    << "of node " << pPositiveSideControlNode->Id() << "\n";
                mpImpl->mpModelPart->AddMasterSlaveConstraint(p_constraint);

                KRATOS_INFO_IF(this->Info(), 3 <= mpImpl->mVerbosity)
                    << "insert condition " << mpPositiveSideLoad->Id() <<' '
                    << "loading " << pPositiveSideControlNode->GetDofs()[0]->GetVariable().Name() << ' '
                    << "of node " << pPositiveSideControlNode->Id() << "\n";
                mpImpl->mpModelPart->AddCondition(mpPositiveSideLoad);

                ++constraint_id;
                ++condition_id;
            }

            // Insert the negative side constraint and condition.
            {
                MasterSlaveConstraint::Pointer p_constraint(new OutOfPlaneDisplacementConstraint<false>(
                    constraint_id,
                    mpImpl->mpModelPart->Geometries(),
                    std::move(negative_side_dofs),
                    {constraint_id},
                    p_protected_surface_normal,
                    mpImpl->mpAdjacencyMap,
                    mpImpl->mVerbosity));
                mpNegativeSideLoad = Condition::Pointer(new PointLoadCondition1D1N(
                    condition_id,
                    Geometry<Node>::Pointer(new Point2D<Node>(p_negative_side_control_node)),
                    p_properties));

                KRATOS_INFO_IF(this->Info(), 3 <= mpImpl->mVerbosity)
                    << "insert constraint " << p_constraint->Id() <<' '
                    << "tying the average out-of-plane displacement on the positive side to "
                    << p_negative_side_control_node->GetDofs()[0]->GetVariable().Name() << ' '
                    << "of node " << p_negative_side_control_node->Id() << "\n";
                mpImpl->mpModelPart->AddMasterSlaveConstraint(p_constraint);

                KRATOS_INFO_IF(this->Info(), 3 <= mpImpl->mVerbosity)
                    << "insert condition " << mpNegativeSideLoad->Id() <<' '
                    << "loading " << p_negative_side_control_node->GetDofs()[0]->GetVariable().Name() << ' '
                    << "of node " << p_negative_side_control_node->Id() << "\n";
                mpImpl->mpModelPart->AddCondition(mpNegativeSideLoad);

                ++constraint_id;
                ++condition_id;
            }
        KRATOS_CATCH("")
}


void InsertNeumannPreTensionOperation::Apply(double Magnitude) {
    KRATOS_ERROR_IF_NOT(mpPositiveSideLoad && mpNegativeSideLoad);
    KRATOS_TRY
        mpPositiveSideLoad->SetValue(POINT_LOAD_X, Magnitude);
        mpNegativeSideLoad->SetValue(POINT_LOAD_X, -Magnitude);
    KRATOS_CATCH("")
}


std::string InsertNeumannPreTensionOperation::Info() const {
    return "InsertNeumannPreTensionOperation";
}


} // namespace Kratos
