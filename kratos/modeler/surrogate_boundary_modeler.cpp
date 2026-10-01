//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//


// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_map>
#include <utility>

// Project includes
#include "geometries/line_2d_2.h"
#include "includes/kratos_components.h"
#include "includes/variables.h"
#include "modeler/surrogate_boundary_modeler.h"
#include "spatial_containers/bins_dynamic.h"
#include "utilities/geometry_utilities/nearest_point_utilities.h"
#include "utilities/structured_surrogate_boundary_utility.h"

namespace Kratos {

namespace {

using TopologyUtility2D = StructuredSurrogateBoundaryUtility<2>;
using IndexArray2D = TopologyUtility2D::IndexArrayType;
using SkinSearchPointType = Point;
using SkinSearchPointPointerType = SkinSearchPointType::Pointer;
using SkinSearchPointVectorType = std::vector<SkinSearchPointPointerType>;
using SkinSearchPointIteratorType = SkinSearchPointVectorType::iterator;
using SkinSearchDistanceVectorType = std::vector<double>;
using SkinSearchDistanceIteratorType = SkinSearchDistanceVectorType::iterator;
using SkinBinsType = BinsDynamic<
    2,
    SkinSearchPointType,
    SkinSearchPointVectorType,
    SkinSearchPointPointerType,
    SkinSearchPointIteratorType,
    SkinSearchDistanceIteratorType>;
using SkinConditionLookupType = std::unordered_map<const Point*, Condition*>;

struct CellClassificationResult
{
    std::vector<bool> ActiveCells;
    std::size_t NumberOfCutCells = 0;
    std::size_t NumberOfUncutRegions = 0;
};

void CheckLineSkin(const ModelPart::ConditionsContainerType& rSkinConditions)
{
    for (const auto& r_condition : rSkinConditions) {
        const auto& r_geometry = r_condition.GetGeometry();
        KRATOS_ERROR_IF(r_geometry.LocalSpaceDimension() != 1 || r_geometry.size() != 2)
            << "SurrogateBoundaryModeler 2D expects a skin made of Line2 conditions. "
            << "Condition " << r_condition.Id() << " has geometry "
            << r_geometry.Info() << "." << std::endl;
    }
}

bool IsPointInsideClosedLineSkin(
    const ModelPart::ConditionsContainerType& rSkinConditions,
    const Point& rPoint)
{
    bool is_inside = false;
    const double x = rPoint.X();
    const double y = rPoint.Y();

    for (const auto& r_condition : rSkinConditions) {
        const auto& r_geometry = r_condition.GetGeometry();
        const double x_0 = r_geometry[0].X();
        const double y_0 = r_geometry[0].Y();
        const double x_1 = r_geometry[1].X();
        const double y_1 = r_geometry[1].Y();

        const bool crosses_y = (y_0 > y) != (y_1 > y);
        if (crosses_y) {
            const double intersection_x = x_0 + (y - y_0) * (x_1 - x_0) / (y_1 - y_0);
            if (x < intersection_x) {
                is_inside = !is_inside;
            }
        }
    }
    return is_inside;
}

bool CellIntersectsSkin(
    const SkinSearchPointVectorType& rCandidates,
    const std::size_t NumberOfCandidates,
    const SkinConditionLookupType& rConditionBySearchPoint,
    const Point& rLowPoint,
    const Point& rHighPoint)
{
    for (std::size_t candidate_index = 0;
         candidate_index < NumberOfCandidates;
         ++candidate_index) {
        const auto p_condition = rConditionBySearchPoint.at(
            rCandidates[candidate_index].get());
        const auto& r_geometry = p_condition->GetGeometry();
        double entry_parameter = 0.0;
        double exit_parameter = 1.0;
        bool intersects = true;
        for (std::size_t d = 0; d < 2; ++d) {
            const double origin = r_geometry[0][d];
            const double direction = r_geometry[1][d] - origin;
            if (std::abs(direction) <= std::numeric_limits<double>::epsilon()) {
                if (origin < rLowPoint[d] || origin > rHighPoint[d]) {
                    intersects = false;
                    break;
                }
            } else {
                double first = (rLowPoint[d] - origin) / direction;
                double second = (rHighPoint[d] - origin) / direction;
                if (first > second) {
                    std::swap(first, second);
                }
                entry_parameter = std::max(entry_parameter, first);
                exit_parameter = std::min(exit_parameter, second);
                if (entry_parameter > exit_parameter) {
                    intersects = false;
                    break;
                }
            }
        }
        if (intersects) {
            return true;
        }
    }
    return false;
}

Point ProjectOnNearestSkinCondition(
    const Point& rPoint,
    SkinBinsType& rSkinBins,
    const SkinConditionLookupType& rConditionBySearchPoint,
    SkinSearchPointVectorType& rCandidates,
    SkinSearchDistanceVectorType& rCandidateDistances,
    const double MaximumHalfSegmentLength,
    const double InitialSearchRadius,
    const double MaximumSearchRadius,
    const double SearchTolerance)
{
    double search_radius = InitialSearchRadius;
    while (true) {
        const std::size_t number_of_candidates = rSkinBins.SearchInRadius(
            rPoint,
            search_radius,
            rCandidates.begin(),
            rCandidateDistances.begin(),
            rCandidates.size());

        Point closest_point;
        double closest_distance = std::numeric_limits<double>::max();
        for (std::size_t candidate_index = 0;
             candidate_index < number_of_candidates;
             ++candidate_index) {
            const auto p_condition = rConditionBySearchPoint.at(
                rCandidates[candidate_index].get());
            const auto& r_geometry = p_condition->GetGeometry();
            const Point projected_point = NearestPointUtilities::LineNearestPoint(
                rPoint.Coordinates(),
                r_geometry[0].Coordinates(),
                r_geometry[1].Coordinates());
            const double distance = std::hypot(
                projected_point.X() - rPoint.X(),
                projected_point.Y() - rPoint.Y());
            if (distance < closest_distance) {
                closest_distance = distance;
                closest_point = projected_point;
            }
        }

        // A segment whose center lies outside the search circle cannot be
        // closer than search_radius minus its half length. This makes the
        // center-based bin search exact for segments of different lengths.
        const bool closest_is_global =
            closest_distance <= search_radius - MaximumHalfSegmentLength + SearchTolerance;
        if (number_of_candidates > 0 &&
            (closest_is_global || search_radius >= MaximumSearchRadius)) {
            return closest_point;
        }

        KRATOS_ERROR_IF(search_radius >= MaximumSearchRadius)
            << "Could not find a true-skin condition within radius "
            << MaximumSearchRadius << " for surrogate point " << rPoint << "."
            << std::endl;
        search_radius = std::min(2.0 * search_radius, MaximumSearchRadius);
    }
}

CellClassificationResult ClassifyCellsWithSpatialSearch(
    const IndexArray2D& rNumberOfCells,
    const Vector& rLower,
    const double Dx,
    const double Dy,
    const bool PhysicalDomainIsInside,
    const double Lambda,
    const int SamplesPerDirection,
    const ModelPart::ConditionsContainerType& rSkinConditions,
    SkinBinsType& rSkinBins,
    const SkinConditionLookupType& rConditionBySearchPoint,
    const std::size_t NumberOfSkinConditions,
    const double MaximumHalfSegmentLength,
    const double SearchTolerance)
{
    const std::size_t number_of_cells = rNumberOfCells[0] * rNumberOfCells[1];
    std::vector<bool> cut_cells(number_of_cells, false);
    std::vector<int> region_inside_state(number_of_cells, -1);
    SkinSearchPointVectorType candidates(NumberOfSkinConditions);
    SkinSearchDistanceVectorType candidate_distances(NumberOfSkinConditions);
    const double cell_search_radius =
        0.5 * std::hypot(Dx, Dy) + MaximumHalfSegmentLength + SearchTolerance;

    for (std::size_t j = 0; j < rNumberOfCells[1]; ++j) {
        for (std::size_t i = 0; i < rNumberOfCells[0]; ++i) {
            const IndexArray2D index{i, j};
            const std::size_t linear_index = TopologyUtility2D::LinearIndex(index, rNumberOfCells);
            const Point low_point(rLower[0] + i * Dx, rLower[1] + j * Dy, rLower[2]);
            const Point high_point(low_point.X() + Dx, low_point.Y() + Dy, rLower[2]);
            const Point center(
                0.5 * (low_point.X() + high_point.X()),
                0.5 * (low_point.Y() + high_point.Y()),
                rLower[2]);
            const std::size_t number_of_candidates = rSkinBins.SearchInRadius(
                center,
                cell_search_radius,
                candidates.begin(),
                candidate_distances.begin(),
                candidates.size());
            if (CellIntersectsSkin(
                    candidates,
                    number_of_candidates,
                    rConditionBySearchPoint,
                    low_point,
                    high_point)) {
                cut_cells[linear_index] = true;
                region_inside_state[linear_index] = -2;
            }
        }
    }

    std::size_t number_of_regions = 0;
    constexpr std::array<std::array<int, 2>, 4> neighbour_offsets{{
        {{-1, 0}}, {{1, 0}}, {{0, -1}}, {{0, 1}}
    }};
    std::queue<std::size_t> cells_to_visit;

    for (std::size_t seed = 0; seed < number_of_cells; ++seed) {
        if (region_inside_state[seed] != -1) {
            continue;
        }

        const std::size_t seed_i = seed % rNumberOfCells[0];
        const std::size_t seed_j = seed / rNumberOfCells[0];
        const Point representative(
            rLower[0] + (static_cast<double>(seed_i) + 0.5) * Dx,
            rLower[1] + (static_cast<double>(seed_j) + 0.5) * Dy,
            rLower[2]);
        const int is_inside = IsPointInsideClosedLineSkin(
            rSkinConditions, representative) ? 1 : 0;
        region_inside_state[seed] = is_inside;
        cells_to_visit.push(seed);
        ++number_of_regions;

        while (!cells_to_visit.empty()) {
            const std::size_t current = cells_to_visit.front();
            cells_to_visit.pop();
            const int current_i = static_cast<int>(current % rNumberOfCells[0]);
            const int current_j = static_cast<int>(current / rNumberOfCells[0]);

            for (const auto& r_offset : neighbour_offsets) {
                const int neighbour_i = current_i + r_offset[0];
                const int neighbour_j = current_j + r_offset[1];
                if (neighbour_i < 0 || neighbour_j < 0 ||
                    neighbour_i >= static_cast<int>(rNumberOfCells[0]) ||
                    neighbour_j >= static_cast<int>(rNumberOfCells[1])) {
                    continue;
                }
                const IndexArray2D neighbour_index{
                    static_cast<std::size_t>(neighbour_i),
                    static_cast<std::size_t>(neighbour_j)};
                const std::size_t neighbour = TopologyUtility2D::LinearIndex(
                    neighbour_index, rNumberOfCells);
                if (region_inside_state[neighbour] == -1) {
                    region_inside_state[neighbour] = is_inside;
                    cells_to_visit.push(neighbour);
                }
            }
        }
    }

    CellClassificationResult result;
    result.ActiveCells.resize(number_of_cells, false);
    result.NumberOfCutCells = std::count(cut_cells.begin(), cut_cells.end(), true);
    result.NumberOfUncutRegions = number_of_regions;
    const double epsilon = std::numeric_limits<double>::epsilon();
    const bool retain_all_cut_cells = Lambda <= epsilon;
    const bool discard_all_cut_cells = Lambda >= 1.0 - epsilon;
    const std::size_t total_samples = static_cast<std::size_t>(
        SamplesPerDirection * SamplesPerDirection);

    for (std::size_t linear_index = 0; linear_index < number_of_cells; ++linear_index) {
        if (!cut_cells[linear_index]) {
            result.ActiveCells[linear_index] =
                (region_inside_state[linear_index] == 1) == PhysicalDomainIsInside;
            continue;
        }
        if (retain_all_cut_cells) {
            result.ActiveCells[linear_index] = true;
            continue;
        }
        if (discard_all_cut_cells) {
            continue;
        }

        const std::size_t i = linear_index % rNumberOfCells[0];
        const std::size_t j = linear_index / rNumberOfCells[0];
        std::size_t physical_samples = 0;
        for (int sample_j = 0; sample_j < SamplesPerDirection; ++sample_j) {
            for (int sample_i = 0; sample_i < SamplesPerDirection; ++sample_i) {
                const Point sample(
                    rLower[0] + (static_cast<double>(i) +
                        (static_cast<double>(sample_i) + 0.5) / SamplesPerDirection) * Dx,
                    rLower[1] + (static_cast<double>(j) +
                        (static_cast<double>(sample_j) + 0.5) / SamplesPerDirection) * Dy,
                    rLower[2]);
                if (IsPointInsideClosedLineSkin(rSkinConditions, sample) == PhysicalDomainIsInside) {
                    ++physical_samples;
                }
            }
        }
        const double physical_fraction =
            static_cast<double>(physical_samples) / total_samples;
        result.ActiveCells[linear_index] = physical_fraction >= Lambda;
    }

    return result;
}

std::size_t FindMaximumNodeId(const ModelPart& rModelPart)
{
    std::size_t maximum_id = 0;
    for (const auto& r_node : rModelPart.Nodes()) {
        maximum_id = std::max(maximum_id, r_node.Id());
    }
    return maximum_id;
}

std::size_t FindMaximumElementId(const ModelPart& rModelPart)
{
    std::size_t maximum_id = 0;
    for (const auto& r_element : rModelPart.Elements()) {
        maximum_id = std::max(maximum_id, r_element.Id());
    }
    return maximum_id;
}

std::size_t FindMaximumConditionId(const ModelPart& rModelPart)
{
    std::size_t maximum_id = 0;
    for (const auto& r_condition : rModelPart.Conditions()) {
        maximum_id = std::max(maximum_id, r_condition.Id());
    }
    return maximum_id;
}

} // unnamed namespace

SurrogateBoundaryModeler::SurrogateBoundaryModeler(
    Model& rModel,
    Parameters ModelerParameters)
    : Modeler(rModel, ModelerParameters),
      mpModel(&rModel)
{
    ModelerParameters.ValidateAndAssignDefaults(GetDefaultParameters());
    mParameters = ModelerParameters;
}

const Parameters SurrogateBoundaryModeler::GetDefaultParameters() const
{
    return Parameters(R"({
        "echo_level" : 0,
        "domain_size" : 2,
        "input_model_part_name" : "SkinModelPart",
        "output_model_part_name" : "SurrogateModelPart",
        "volume_sub_model_part_name" : "FluidDomain",
        "surrogate_sub_model_part_name" : "SurrogateBoundary",
        "lower_point" : [0.0, 0.0, 0.0],
        "upper_point" : [1.0, 1.0, 0.0],
        "number_of_elements" : [10, 10, 1],
        "physical_domain" : "outside",
        "lambda" : 1.0,
        "classification_samples_per_direction" : 3,
        "element_name" : "Element2D4N",
        "condition_name" : "LineCondition2D2N",
        "properties_id" : 1,
        "number_of_boundary_gauss_points" : 2
    })");
}

void SurrogateBoundaryModeler::SetupModelPart()
{
    KRATOS_ERROR_IF(mpModel == nullptr) << "SurrogateBoundaryModeler has no Model container." << std::endl;
    const int domain_size = mParameters["domain_size"].GetInt();
    if (domain_size == 2) {
        SetupModelPart2D();
    } else {
        KRATOS_ERROR << "SurrogateBoundaryModeler domain_size=" << domain_size
                     << " is not implemented yet. The structured topology kernel supports 3D; "
                     << "the triangle-skin classifier and Hex8/Quad4 entity adapter are the remaining extension." << std::endl;
    }
}

void SurrogateBoundaryModeler::SetupModelPart2D()
{
    const std::string skin_name = mParameters["input_model_part_name"].GetString();
    const std::string output_name = mParameters["output_model_part_name"].GetString();
    KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(skin_name))
        << "Input skin ModelPart '" << skin_name << "' does not exist." << std::endl;

    ModelPart& r_skin = mpModel->GetModelPart(skin_name);
    KRATOS_ERROR_IF(r_skin.NumberOfConditions() == 0)
        << "Input skin ModelPart '" << skin_name << "' has no conditions." << std::endl;

    ModelPart& r_output = mpModel->HasModelPart(output_name)
        ? mpModel->GetModelPart(output_name)
        : mpModel->CreateModelPart(output_name);
    r_output.GetProcessInfo()[DOMAIN_SIZE] = 2;

    ModelPart& r_volume = r_output.HasSubModelPart(mParameters["volume_sub_model_part_name"].GetString())
        ? r_output.GetSubModelPart(mParameters["volume_sub_model_part_name"].GetString())
        : r_output.CreateSubModelPart(mParameters["volume_sub_model_part_name"].GetString());
    ModelPart& r_surrogate = r_output.HasSubModelPart(mParameters["surrogate_sub_model_part_name"].GetString())
        ? r_output.GetSubModelPart(mParameters["surrogate_sub_model_part_name"].GetString())
        : r_output.CreateSubModelPart(mParameters["surrogate_sub_model_part_name"].GetString());

    const auto lower = mParameters["lower_point"].GetVector();
    const auto upper = mParameters["upper_point"].GetVector();
    const auto divisions_vector = mParameters["number_of_elements"].GetVector();
    KRATOS_ERROR_IF(lower.size() != 3 || upper.size() != 3 || divisions_vector.size() != 3)
        << "lower_point, upper_point and number_of_elements must have three entries." << std::endl;

    const IndexArray2D number_of_cells{
        static_cast<std::size_t>(std::lround(divisions_vector[0])),
        static_cast<std::size_t>(std::lround(divisions_vector[1]))};
    KRATOS_ERROR_IF(number_of_cells[0] == 0 || number_of_cells[1] == 0)
        << "number_of_elements must be positive in x and y." << std::endl;
    KRATOS_ERROR_IF(upper[0] <= lower[0] || upper[1] <= lower[1])
        << "upper_point must be greater than lower_point in x and y." << std::endl;

    const double dx = (upper[0] - lower[0]) / number_of_cells[0];
    const double dy = (upper[1] - lower[1]) / number_of_cells[1];
    const bool physical_domain_is_inside = mParameters["physical_domain"].GetString() == "inside";
    KRATOS_ERROR_IF(!physical_domain_is_inside && mParameters["physical_domain"].GetString() != "outside")
        << "physical_domain must be either 'inside' or 'outside'." << std::endl;
    const double lambda = mParameters["lambda"].GetDouble();
    KRATOS_ERROR_IF(lambda < 0.0 || lambda > 1.0) << "lambda must lie in [0,1]." << std::endl;
    const int samples_per_direction = mParameters["classification_samples_per_direction"].GetInt();
    KRATOS_ERROR_IF(samples_per_direction < 1)
        << "classification_samples_per_direction must be positive." << std::endl;

    auto& r_skin_conditions = r_skin.Conditions();
    CheckLineSkin(r_skin_conditions);
    const double domain_scale = std::max({
        upper[0] - lower[0], upper[1] - lower[1], 1.0});
    const double search_tolerance = 1.0e-12 * domain_scale;
    SkinSearchPointVectorType skin_search_points;
    skin_search_points.reserve(r_skin_conditions.size());
    SkinConditionLookupType condition_by_search_point;
    condition_by_search_point.reserve(r_skin_conditions.size());
    double maximum_half_segment_length = 0.0;
    for (auto& r_condition : r_skin_conditions) {
        const auto& r_geometry = r_condition.GetGeometry();
        auto p_search_point = Kratos::make_shared<Point>(
            0.5 * (r_geometry[0].X() + r_geometry[1].X()),
            0.5 * (r_geometry[0].Y() + r_geometry[1].Y()),
            0.5 * (r_geometry[0].Z() + r_geometry[1].Z()));
        maximum_half_segment_length = std::max(
            maximum_half_segment_length,
            0.5 * std::hypot(
                r_geometry[1].X() - r_geometry[0].X(),
                r_geometry[1].Y() - r_geometry[0].Y()));
        condition_by_search_point.emplace(p_search_point.get(), &r_condition);
        skin_search_points.push_back(std::move(p_search_point));
    }
    SkinBinsType skin_bins(skin_search_points.begin(), skin_search_points.end());
    const auto classification = ClassifyCellsWithSpatialSearch(
        number_of_cells,
        lower,
        dx,
        dy,
        physical_domain_is_inside,
        lambda,
        samples_per_direction,
        r_skin_conditions,
        skin_bins,
        condition_by_search_point,
        skin_search_points.size(),
        maximum_half_segment_length,
        search_tolerance);
    const auto& active_cells = classification.ActiveCells;
    const double initial_projection_search_radius =
        std::hypot(dx, dy) + maximum_half_segment_length + search_tolerance;
    const double maximum_projection_search_radius = std::hypot(
        upper[0] - lower[0], upper[1] - lower[1]) +
        maximum_half_segment_length + initial_projection_search_radius;

    auto p_properties = r_output.HasProperties(mParameters["properties_id"].GetInt())
        ? r_output.pGetProperties(mParameters["properties_id"].GetInt())
        : r_output.CreateNewProperties(mParameters["properties_id"].GetInt());
    const auto& r_element_prototype = KratosComponents<Element>::Get(mParameters["element_name"].GetString());
    const auto& r_condition_prototype = KratosComponents<Condition>::Get(mParameters["condition_name"].GetString());

    std::size_t next_node_id = FindMaximumNodeId(r_output) + 1;
    std::size_t next_element_id = FindMaximumElementId(r_output) + 1;
    std::size_t next_condition_id = FindMaximumConditionId(r_output) + 1;
    const std::size_t number_of_grid_nodes_x = number_of_cells[0] + 1;
    std::unordered_map<std::size_t, Node::Pointer> grid_nodes;
    std::vector<Element::Pointer> cell_elements(active_cells.size());

    auto get_node = [&](const std::size_t i, const std::size_t j) {
        const std::size_t key = i + j * number_of_grid_nodes_x;
        const auto found = grid_nodes.find(key);
        if (found != grid_nodes.end()) {
            return found->second;
        }
        auto p_node = r_output.CreateNewNode(next_node_id++, lower[0] + i * dx, lower[1] + j * dy, lower[2]);
        grid_nodes.emplace(key, p_node);
        r_volume.AddNode(p_node);
        return p_node;
    };

    for (std::size_t j = 0; j < number_of_cells[1]; ++j) {
        for (std::size_t i = 0; i < number_of_cells[0]; ++i) {
            const IndexArray2D index{i, j};
            const std::size_t linear_index = TopologyUtility2D::LinearIndex(index, number_of_cells);
            if (!active_cells[linear_index]) {
                continue;
            }

            Element::NodesArrayType nodes;
            nodes.push_back(get_node(i, j));
            nodes.push_back(get_node(i + 1, j));
            nodes.push_back(get_node(i + 1, j + 1));
            nodes.push_back(get_node(i, j + 1));
            auto p_element = r_element_prototype.Create(next_element_id++, nodes, p_properties);
            r_output.AddElement(p_element);
            r_volume.AddElement(p_element);
            cell_elements[linear_index] = p_element;
        }
    }

    const auto boundary_faces = TopologyUtility2D::ExtractBoundaryFaces(active_cells, number_of_cells, false);
    const int number_of_gauss_points = mParameters["number_of_boundary_gauss_points"].GetInt();
    KRATOS_ERROR_IF(number_of_gauss_points != 2)
        << "The initial 2D face adapter currently supports exactly two Gauss points." << std::endl;
    SkinSearchPointVectorType projection_candidates(skin_search_points.size());
    SkinSearchDistanceVectorType projection_candidate_distances(skin_search_points.size());

    for (const auto& r_face : boundary_faces) {
        const std::size_t i = r_face.OwnerCell[0];
        const std::size_t j = r_face.OwnerCell[1];
        Condition::NodesArrayType face_nodes;

        if (r_face.Direction == 0 && !r_face.IsPositiveSide) {
            face_nodes.push_back(get_node(i, j + 1));
            face_nodes.push_back(get_node(i, j));
        } else if (r_face.Direction == 0) {
            face_nodes.push_back(get_node(i + 1, j));
            face_nodes.push_back(get_node(i + 1, j + 1));
        } else if (!r_face.IsPositiveSide) {
            face_nodes.push_back(get_node(i, j));
            face_nodes.push_back(get_node(i + 1, j));
        } else {
            face_nodes.push_back(get_node(i + 1, j + 1));
            face_nodes.push_back(get_node(i, j + 1));
        }

        const std::size_t owner_index = TopologyUtility2D::LinearIndex(r_face.OwnerCell, number_of_cells);
        auto p_line_geometry = Kratos::make_shared<Line2D2<Node>>(face_nodes);
        Condition::NodesArrayType condition_nodes;
        const auto& r_owner_geometry = cell_elements[owner_index]->GetGeometry();
        for (std::size_t n = 0; n < r_owner_geometry.size(); ++n) {
            condition_nodes.push_back(r_owner_geometry.pGetPoint(n));
        }
        auto p_condition = r_condition_prototype.Create(
            next_condition_id++, condition_nodes, p_properties);
        p_condition->SetValue(NEIGHBOUR_ELEMENTS, GlobalPointersVector<Element>({cell_elements[owner_index]}));

        Matrix face_coordinates(2, 3);
        for (std::size_t n = 0; n < 2; ++n) {
            for (std::size_t d = 0; d < 3; ++d) {
                face_coordinates(n, d) = face_nodes[n][d];
            }
        }
        p_condition->SetValue(SURROGATE_BOUNDARY_FACE_COORDINATES, face_coordinates);

        const auto integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
        const auto& r_integration_points = p_line_geometry->IntegrationPoints(integration_method);
        const Matrix& r_N = p_line_geometry->ShapeFunctionsValues(integration_method);
        Matrix projections(r_integration_points.size(), 3);
        for (std::size_t g = 0; g < r_integration_points.size(); ++g) {
            Point surrogate_point(0.0, 0.0, 0.0);
            for (std::size_t n = 0; n < face_nodes.size(); ++n) {
                surrogate_point.Coordinates() += r_N(g, n) * face_nodes[n].Coordinates();
            }
            const Point projected_point = ProjectOnNearestSkinCondition(
                surrogate_point,
                skin_bins,
                condition_by_search_point,
                projection_candidates,
                projection_candidate_distances,
                maximum_half_segment_length,
                initial_projection_search_radius,
                maximum_projection_search_radius,
                search_tolerance);
            for (std::size_t d = 0; d < 3; ++d) {
                projections(g, d) = projected_point[d];
            }
        }
        p_condition->SetValue(SURROGATE_BOUNDARY_PROJECTION, projections);
        KRATOS_INFO_IF("SurrogateBoundaryModeler", mParameters["echo_level"].GetInt() > 1)
            << "condition=" << p_condition->Id()
            << ", owner=(" << i << "," << j << ")"
            << ", direction=" << r_face.Direction
            << ", positive_side=" << r_face.IsPositiveSide
            << ", face=[(" << face_coordinates(0, 0) << "," << face_coordinates(0, 1)
            << "),(" << face_coordinates(1, 0) << "," << face_coordinates(1, 1) << ")]"
            << std::endl;
        r_output.AddCondition(p_condition);
        r_surrogate.AddCondition(p_condition);
        for (std::size_t n = 0; n < p_line_geometry->size(); ++n) {
            r_surrogate.AddNode(p_line_geometry->pGetPoint(n));
        }
    }

    KRATOS_INFO_IF("SurrogateBoundaryModeler", mParameters["echo_level"].GetInt() > 0)
        << "Classified " << number_of_cells[0] * number_of_cells[1]
        << " cells (" << classification.NumberOfCutCells << " cut, "
        << classification.NumberOfUncutRegions << " uncut regions); created "
        << r_volume.NumberOfElements() << " active elements and "
        << r_surrogate.NumberOfConditions() << " surrogate faces." << std::endl;
}

} // namespace Kratos
