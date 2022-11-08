//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
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
            KRATOS_ERROR_IF((geom_type != GeometryData::KratosGeometryType::Kratos_Line2D2) && geom_type != GeometryData::KratosGeometryType::Kratos_Line3D2) << "The geometry type is not correct, it is suppossed to be a line" << std::endl;
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

    // Synchronize variables (TODO: Should be absolute minimum)
    if constexpr ( THistorical) {
        r_comm.SynchronizeCurrentDataToMin(*mpDistanceVariable);
    } else {
        r_comm.SynchronizeNonHistoricalDataToMin(*mpDistanceVariable);
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
    // TODO: Use BinUtils
    const double radius_path = mThisParameters["radius_path"].GetDouble();
    const double distance_tolerance = mThisParameters["distance_tolerance"].GetDouble();
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode) {
        double min_value = std::numeric_limits<double>::max();
        for (auto& p_segment : rVectorSegments) {
            const double potential_min = GeometricalProjectionUtilities::FastMinimalDistanceOnLine(*p_segment, rNode, distance_tolerance);
            // const auto distance_computed_type = FastMinimalDistanceOnLineWithRadius(potential_min, *p_segment, rNode, radius_path, distance_tolerance);
            min_value = std::abs(potential_min) < std::abs(min_value) ? potential_min : min_value;
        }
        if constexpr (THistorical) {
            rNode.FastGetSolutionStepValue(*mpDistanceVariable) = min_value - radius_path;
        } else {
            rNode.GetValue(*mpDistanceVariable) = min_value - radius_path;
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
    const double radius_path = mThisParameters["radius_path"].GetDouble();
    const double distance_tolerance = mThisParameters["distance_tolerance"].GetDouble();
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode) {
        double min_value = std::numeric_limits<double>::max();
        for (auto& p_segment : rVectorSegments) {
            const double potential_min = GeometricalProjectionUtilities::FastMinimalDistanceOnLine(*p_segment, rNode, distance_tolerance);
            // const auto distance_computed_type = FastMinimalDistanceOnLineWithRadius(potential_min, *p_segment, rNode, radius_path, distance_tolerance);
            min_value = std::abs(potential_min) < std::abs(min_value) ? potential_min : min_value;
        }
        if constexpr (THistorical) {
            rNode.FastGetSolutionStepValue(*mpDistanceVariable) = min_value - radius_path;
        } else {
            rNode.GetValue(*mpDistanceVariable) = min_value - radius_path;
        }
    });
}

/***********************************************************************************/
/***********************************************************************************/

// template<bool THistorical>
// DistanceComputed CalculateDistanceToPathProcess<THistorical>::FastMinimalDistanceOnLineWithRadius(
//     double& rDistance,
//     const Geometry<NodeType>& rSegment,
//     const Point& rPoint,
//     const double Radius,
//     const double Tolerance
//     )
// {
//     // If radius is zero, we compute the distance to the line
//     Line3D2<NodeType> line(rSegment.Points()); // NOTE: TO ENSURE THAT IT IS 3D IN CASE IS DECLARED AS 2D
//     if (Radius < std::numeric_limits<double>::epsilon()) {
//         rDistance = GeometricalProjectionUtilities::FastMinimalDistanceOnLine(line, rPoint, Tolerance);
//         return DistanceComputed::NO_RADIUS;
//     } else {
//         Point projected_point;
//         const double projected_distance = GeometricalProjectionUtilities::FastProjectOnLine(line, rPoint, projected_point);
//         array_1d<double, 3> vector_line = line[1].Coordinates() - line[0].Coordinates();
//         vector_line /= norm_2(vector_line);
//         typename Geometry<NodeType>::CoordinatesArrayType projected_local;
//         // If projection is inside, just remove the radius
//         if (line.IsInside(projected_point.Coordinates(), projected_local, Tolerance)) {
//             rDistance = projected_distance - Radius;
//             return DistanceComputed::RADIUS_PROJECTED;
//         } else { // Othwerise we compute the distance to the closest node and compute the difference with the "radius cylinder"
//             // Distances to the nodes
//             const Point::Pointer point = Kratos::make_shared<Point>(rPoint.Coordinates());
//             const Point::Pointer aux_projected_point = Kratos::make_shared<Point>(projected_point.Coordinates());
//             const Point::Pointer point_a = Kratos::make_shared<Point>(line[0].Coordinates());
//             const Point::Pointer point_b = Kratos::make_shared<Point>(line[1].Coordinates());
//             const double distance_a = rPoint.Distance(*point_a);
//             const double distance_b = rPoint.Distance(*point_b);

//             // Positive distance. Remove distance to parallel line in Radius
//             if (projected_distance > Radius) {
//                 const double N_line = Radius/projected_distance;
//                 const Point::Pointer point_distance_r_projection = Kratos::make_shared<Point>(N_line * projected_point.Coordinates() + (1.0 - N_line) * rPoint.Coordinates());
//                 const Point::Pointer aux_point_parallel = Kratos::make_shared<Point>(point_distance_r_projection->Coordinates() + vector_line);
                
//                 // Parallel line to the segment
//                 Geometry<Point>::PointsArrayType points_array_parallel_line;
//                 points_array_parallel_line.reserve(2);
//                 points_array_parallel_line.push_back(point_distance_r_projection);
//                 points_array_parallel_line.push_back(aux_point_parallel);
//                 Line3D2<Point> parallel_line(points_array_parallel_line);

//                 // Line from the point to the segment
//                 Geometry<Point>::PointsArrayType points_array_distance;
//                 points_array_distance.reserve(2);
//                 points_array_distance.push_back(point);
//                 if (distance_a < distance_b) {
//                     points_array_distance.push_back(point_a);
//                 } else {
//                     points_array_distance.push_back(point_b);
//                 }
//                 Line3D2<Point> distance_line(points_array_distance);

//                 // Compute intersection
//                 const auto line_intersection = Line3D2<Point>(IntersectionUtilities::ComputeClosestLineLineIntersection(parallel_line, distance_line));
//                 const Point intersection_point(line_intersection.Center().Coordinates());

//                 // Compute distance
//                 if (distance_a < distance_b) {
//                     rDistance = distance_a - point_a->Distance(intersection_point);
//                 } else {
//                     rDistance = distance_b - point_b->Distance(intersection_point);
//                 }
//                 return DistanceComputed::RADIUS_NOT_PROJECTED_OUTSIDE;
//             } else { // Negative distance
//                 array_1d<double, 3> projection_vector = rPoint.Coordinates() - projected_point.Coordinates();
//                 projection_vector /= norm_2(projection_vector);

//                 // Parallel projected point
//                 const array_1d<double, 3> aux_parallel_projection = (distance_a < distance_b) ? point_a->Coordinates() + projection_vector * Radius : point_b->Coordinates() + projection_vector * Radius;
//                 Point parallel_projection_point(aux_parallel_projection);

//                 // Compute distance
//                 rDistance = - point->Distance(parallel_projection_point);
//                 return DistanceComputed::RADIUS_NOT_PROJECTED_OUTSIDE;
//             }
//         }
//     }
//     return DistanceComputed::ERROR;
// }

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
        "distance_tolerance"       : 1.0e-9
    })" );

    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template class CalculateDistanceToPathProcess<true>;
template class CalculateDistanceToPathProcess<false>;

} /// namespace Kratos