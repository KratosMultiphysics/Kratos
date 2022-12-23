//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/intersection_utilities.h"
#include "geometries/line_3d_2.h"
#include "utilities/geometrical_projection_utilities.h"
#include "processes/calculate_distance_to_path_process.h"
#include "utilities/variable_utils.h"
#include "spatial_containers/spatial_containers.h" // kd-tree

namespace Kratos
{

template<bool THistorical>
CalculateDistanceToPathProcess<THistorical>::CalculateDistanceToPathProcess(
    Model& rModel, 
    Parameters ThisParameters
    ) : mrModel(rModel),
        mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
    mpDistanceVariable = &KratosComponents<Variable<double>>::Get(mThisParameters["distance_variable_name"].GetString());
}

/***********************************************************************************/
/***********************************************************************************/

template<bool THistorical>
void CalculateDistanceToPathProcess<THistorical>::Execute()
{
    // Getting the model parts
    const std::string& r_distance_model_part_name = mThisParameters["distance_model_part_name"].GetString();
    auto& r_distance_model_part = mrModel.GetModelPart(r_distance_model_part_name);
    const std::string& r_path_model_part_name = mThisParameters["path_model_part_name"].GetString();
    auto& r_path_model_part = mrModel.GetModelPart(r_path_model_part_name);

    // Getting communicator
    auto& r_comm = r_distance_model_part.GetCommunicator();
    auto& r_data_comm = r_comm.GetDataCommunicator();
    
    // Initialize distance variable
    if constexpr ( THistorical) {
        VariableUtils().SetHistoricalVariableToZero(*mpDistanceVariable, r_distance_model_part.Nodes());
        r_comm.SynchronizeVariable(*mpDistanceVariable);
    } else {
        VariableUtils().SetNonHistoricalVariableToZero(*mpDistanceVariable, r_distance_model_part.Nodes());
        r_comm.SynchronizeNonHistoricalVariable(*mpDistanceVariable);
    }

    /* Getting a global vector of Line3D2N geometries */
    // First we check that every element in the modelpart is a line
    const auto& r_elements_array = r_path_model_part.Elements();
    // By default we consider that the model part is composed of elements
    if (r_elements_array.size() > 0) {
        for (auto& r_elem : r_elements_array) {
            const auto& r_geom = r_elem.GetGeometry();
            const auto geom_type = r_geom.GetGeometryType();
            KRATOS_ERROR_IF((geom_type != GeometryData::KratosGeometryType::Kratos_Line2D2) && geom_type != GeometryData::KratosGeometryType::Kratos_Line3D2) << "The geometry type is not correct, it is suppossed to be a linear line" << std::endl;
        }
    } else {
        KRATOS_ERROR << "The model part is empty. Please check the model part provided" << std::endl;
    }

    // Now we create the vector of geometries
    const unsigned int size_vector = r_elements_array.size();
    if (r_distance_model_part.IsDistributed()) {
        r_data_comm.SumAll(size_vector);
    }
    std::vector<Geometry<NodeType>::Pointer> geometries_vector(size_vector);
    if (r_distance_model_part.IsDistributed()) {
        std::size_t counter = 0;
        for (int i=0; i<r_data_comm.Size(); ++i) {
            auto& r_path_local_mesh = r_comm.LocalMesh(i);
            auto& r_local_elements_array = r_path_local_mesh.Elements();
            const auto it_elem_begin = r_local_elements_array.begin();
            const std::size_t number_elements = r_local_elements_array.size();
            IndexPartition<std::size_t>(number_elements).for_each([&](std::size_t i){
                geometries_vector[counter + i] = (it_elem_begin + i)->pGetGeometry();
            });
            counter += number_elements;
        }
    } else {
            const auto it_elem_begin = r_elements_array.begin();
            const std::size_t number_elements = r_elements_array.size();
            IndexPartition<std::size_t>(number_elements).for_each([&](std::size_t i){
                geometries_vector[i] = (it_elem_begin + i)->pGetGeometry();
            });
    }

    // Gettings if compute by brute force or not
    const bool brute_force_calculation = mThisParameters["brute_force_calculation"].GetBool();
    if (brute_force_calculation) {
        this->CalculateDistanceByBruteForce(r_distance_model_part, geometries_vector);
    } else {
        this->CalculateDistance(r_distance_model_part, geometries_vector);
    }

    // Synchronize variables
    if constexpr ( THistorical) {
        r_comm.SynchronizeCurrentDataToAbsMin(*mpDistanceVariable);
    } else {
        r_comm.SynchronizeNonHistoricalDataToAbsMin(*mpDistanceVariable);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<bool THistorical>
void CalculateDistanceToPathProcess<THistorical>::CalculateDistance(
    ModelPart& rModelPart,
    std::vector<Geometry<NodeType>::Pointer>& rVectorSegments
    )
{
    // Max length of the segments considered
    double max_length = 0.0;
    max_length = block_for_each<MaxReduction<double>>(rVectorSegments, [&](Geometry<NodeType>::Pointer pGeometry) {
        return pGeometry->Length();
    });

    /// Type definitions for the tree
    using PointType = PointGeometry;
    using PointTypePointer = PointType::Pointer;
    using PointVector = std::vector<PointTypePointer>;
    using PointIterator = PointVector::iterator;
    using DistanceVector = std::vector<double>;
    using DistanceIterator = DistanceVector::iterator;

    /// KDtree definitions
    using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator >;
    using KDTree = Tree< KDTreePartition<BucketType> >;

    // Some auxiliary values
    const IndexType allocation_size = mThisParameters["search_parameters"]["allocation_size"].GetInt();                 // Allocation size for the vectors and max number of potential results
    const double search_factor = mThisParameters["search_parameters"]["search_factor"].GetDouble();                     // The search factor to be considered
    const double search_increment_factor = mThisParameters["search_parameters"]["search_increment_factor"].GetDouble(); // The search increment factor to be considered
    IndexType bucket_size = mThisParameters["search_parameters"]["bucket_size"].GetInt();                               // Bucket size for kd-tree

    KRATOS_ERROR_IF(rVectorSegments.size() == 0) << "Path not initialized" << std::endl;
    PointVector points_vector;
    points_vector.reserve(rVectorSegments.size());
    for (auto& p_geom : rVectorSegments) {
        points_vector.push_back(PointTypePointer(new PointType(p_geom)));
    }
    KDTree tree_points(points_vector.begin(), points_vector.end(), bucket_size);

    const double radius_path = mThisParameters["radius_path"].GetDouble();
    const double distance_tolerance = mThisParameters["distance_tolerance"].GetDouble();
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode) {
        double search_radius = search_factor * max_length;

        // Initialize values
        PointVector points_found(allocation_size);
        IndexType number_points_found = 0;
        while (number_points_found == 0) {
            search_radius *= search_increment_factor;
            const PointGeometry point(rNode.X(), rNode.Y(), rNode.Z());
            number_points_found = tree_points.SearchInRadius(point, search_radius, points_found.begin(), allocation_size);
        }

        double min_value = std::numeric_limits<double>::max();
        Geometry<NodeType>::Pointer p_closest_geometry = nullptr;
        for (IndexType i = 0; i < number_points_found; ++i) {
            auto p_point = points_found[i];
            auto p_segment = p_point->pGetGeometry();
            const double potential_min = GeometricalProjectionUtilities::FastMinimalDistanceOnLine(*p_segment, rNode, distance_tolerance);
            if (std::abs(potential_min) < std::abs(min_value)) {
                min_value = potential_min;
                p_closest_geometry = p_segment;
            }
        }

        FastMinimalDistanceOnLineWithRadius(min_value, *p_closest_geometry, rNode, radius_path, distance_tolerance);
        if constexpr (THistorical) {
            rNode.FastGetSolutionStepValue(*mpDistanceVariable) = min_value;
        } else {
            rNode.GetValue(*mpDistanceVariable) = min_value;
        }
    });
}

/***********************************************************************************/
/***********************************************************************************/

template<bool THistorical>
void CalculateDistanceToPathProcess<THistorical>::CalculateDistanceByBruteForce(
    ModelPart& rModelPart,
    std::vector<Geometry<NodeType>::Pointer>& rVectorSegments
    )
{
    KRATOS_ERROR_IF(rVectorSegments.size() == 0) << "Path not initialized" << std::endl;

    const double radius_path = mThisParameters["radius_path"].GetDouble();
    const double distance_tolerance = mThisParameters["distance_tolerance"].GetDouble();
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode) {
        double min_value = std::numeric_limits<double>::max();
        Geometry<NodeType>::Pointer p_closest_geometry = nullptr;
        for (auto& p_segment : rVectorSegments) {
            const double potential_min = GeometricalProjectionUtilities::FastMinimalDistanceOnLine(*p_segment, rNode, distance_tolerance);
            if (std::abs(potential_min) < std::abs(min_value)) {
                min_value = potential_min;
                p_closest_geometry = p_segment;
            }
        }
        
        FastMinimalDistanceOnLineWithRadius(min_value, *p_closest_geometry, rNode, radius_path, distance_tolerance);
        if constexpr (THistorical) {
            rNode.FastGetSolutionStepValue(*mpDistanceVariable) = min_value;
        } else {
            rNode.GetValue(*mpDistanceVariable) = min_value;
        }
    });
}

/***********************************************************************************/
/***********************************************************************************/

template<bool THistorical>
DistanceComputed CalculateDistanceToPathProcess<THistorical>::FastMinimalDistanceOnLineWithRadius(
    double& rDistance,
    const Geometry<NodeType>& rSegment,
    const Point& rPoint,
    const double Radius,
    const double Tolerance
    )
{
    // If radius is zero, we keep the distance as it is
    if (Radius < std::numeric_limits<double>::epsilon()) {
        return DistanceComputed::NO_RADIUS;
    } else {
        Line3D2<NodeType> line(rSegment.Points()); // NOTE: TO ENSURE THAT IT IS 3D IN CASE IS DECLARED AS 2D
        Point projected_point;
        const double projected_distance = GeometricalProjectionUtilities::FastProjectOnLine(line, rPoint, projected_point);
        typename Geometry<NodeType>::CoordinatesArrayType projected_local;
        // If projection is inside, just remove the radius
        if (line.IsInside(projected_point.Coordinates(), projected_local, Tolerance)) {
            rDistance = projected_distance - Radius;
            return DistanceComputed::RADIUS_PROJECTED;
        } else { // Othwerise we compute the distance to the closest node and compute the difference with the "radius cylinder"
            // Distances to the nodes
            const Point::Pointer point = Kratos::make_shared<Point>(rPoint.Coordinates());
            const Point::Pointer point_a = Kratos::make_shared<Point>(line[0].Coordinates());
            const Point::Pointer point_b = Kratos::make_shared<Point>(line[1].Coordinates());
            const double distance_a = rPoint.Distance(*point_a);
            const double distance_b = rPoint.Distance(*point_b);

            // Positive distance. Remove distance to parallel line in Radius
            if (projected_distance > Radius) {
                array_1d<double, 3> vector_line = line[1].Coordinates() - line[0].Coordinates();
                vector_line /= norm_2(vector_line);
                const double N_line = Radius/projected_distance;
                const Point::Pointer point_distance_r_projection = Kratos::make_shared<Point>((1.0 - N_line) * projected_point.Coordinates() + N_line * rPoint.Coordinates());
                const Point::Pointer aux_point_parallel = Kratos::make_shared<Point>(point_distance_r_projection->Coordinates() + vector_line);

                // Parallel line to the segment
                Geometry<Point>::PointsArrayType points_array_parallel_line;
                points_array_parallel_line.reserve(2);
                points_array_parallel_line.push_back(point_distance_r_projection);
                points_array_parallel_line.push_back(aux_point_parallel);
                Line3D2<Point> parallel_line(points_array_parallel_line);

                // Line from the point to the segment
                Geometry<Point>::PointsArrayType points_array_distance;
                points_array_distance.reserve(2);
                points_array_distance.push_back(point);
                if (distance_a < distance_b) {
                    points_array_distance.push_back(point_a);
                } else {
                    points_array_distance.push_back(point_b);
                }
                Line3D2<Point> distance_line(points_array_distance);

                // Compute intersection
                const auto line_intersection = Line3D2<Point>(IntersectionUtilities::ComputeShortestLineBetweenTwoLines(parallel_line, distance_line));
                const Point intersection_point(line_intersection.Center().Coordinates());

                // Compute distance
                if (distance_a < distance_b) {
                    rDistance = distance_a - point_a->Distance(intersection_point);
                } else {
                    rDistance = distance_b - point_b->Distance(intersection_point);
                }
                return DistanceComputed::RADIUS_NOT_PROJECTED_OUTSIDE;
            } else { // Negative distance
                array_1d<double, 3> projection_vector = rPoint.Coordinates() - projected_point.Coordinates();
                projection_vector /= norm_2(projection_vector);

                // Parallel projected point
                const array_1d<double, 3> aux_parallel_projection = (distance_a < distance_b) ? point_a->Coordinates() + projection_vector * Radius : point_b->Coordinates() + projection_vector * Radius;
                Point parallel_projection_point(aux_parallel_projection);

                // Compute distance
                rDistance = - point->Distance(parallel_projection_point);
                return DistanceComputed::RADIUS_NOT_PROJECTED_INSIDE;
            }
        }
    }
    return DistanceComputed::ERROR;
}

/***********************************************************************************/
/***********************************************************************************/

template<bool THistorical>
const Parameters CalculateDistanceToPathProcess<THistorical>::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"({
        "distance_model_part_name" :  "",
        "path_model_part_name"     :  "",
        "distance_variable_name"   : "DISTANCE",
        "brute_force_calculation"  : false,
        "radius_path"              : 0.0,
        "distance_tolerance"       : 1.0e-9,
        "search_parameters"        :  {
            "allocation_size"         : 100,
            "bucket_size"             : 4,
            "search_factor"           : 2.0,
            "search_increment_factor" : 1.5
        }
    })" );

    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template class CalculateDistanceToPathProcess<true>;
template class CalculateDistanceToPathProcess<false>;

} /// namespace Kratos