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
#include "spatial_containers/point_object.h"
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

    // MPI not supported for the moment
    KRATOS_ERROR_IF(r_data_comm.IsDistributed()) << "MPI not supported for the moment" << std::endl;

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

    /// Type definitions for the KDtree
    using KDtreePointType = PointObject<Geometry<Node>>;
    using KDtreePointTypePointer = KDtreePointType::Pointer;
    using KDtreePointVector = std::vector<KDtreePointTypePointer>;
    using KDtreePointIterator = KDtreePointVector::iterator;
    using DistanceVector = std::vector<double>;
    using DistanceIterator = DistanceVector::iterator;

    /// KDtree definitions
    using BucketType = Bucket< 3ul, KDtreePointType, KDtreePointVector, KDtreePointTypePointer, KDtreePointIterator, DistanceIterator >;
    using KDTree = Tree< KDTreePartition<BucketType> >;

    // Some auxiliary values
    const IndexType allocation_size = mThisParameters["search_parameters"]["allocation_size"].GetInt();                 // Allocation size for the vectors and max number of potential results
    const double search_factor = mThisParameters["search_parameters"]["search_factor"].GetDouble();                     // The search factor to be considered
    const double search_increment_factor = mThisParameters["search_parameters"]["search_increment_factor"].GetDouble(); // The search increment factor to be considered
    IndexType bucket_size = mThisParameters["search_parameters"]["bucket_size"].GetInt();                               // Bucket size for kd-tree

    KRATOS_ERROR_IF(rVectorSegments.size() == 0) << "Path not initialized" << std::endl;
    KDtreePointVector points_vector;
    points_vector.reserve(rVectorSegments.size());
    for (auto& p_geom : rVectorSegments) {
        points_vector.push_back(KDtreePointTypePointer(new KDtreePointType(p_geom)));
    }
    KDTree tree_points(points_vector.begin(), points_vector.end(), bucket_size);

    const double radius_path = mThisParameters["radius_path"].GetDouble();
    const double distance_tolerance = mThisParameters["distance_tolerance"].GetDouble();
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode) {
        double search_radius = search_factor * max_length;

        // Initialize values
        KDtreePointVector points_found(allocation_size);
        IndexType number_points_found = 0;
        while (number_points_found == 0) {
            search_radius *= search_increment_factor;
            const KDtreePointType point(rNode.Coordinates());
            number_points_found = tree_points.SearchInRadius(point, search_radius, points_found.begin(), allocation_size);
        }

        double min_value = std::numeric_limits<double>::max();
        Geometry<NodeType>::Pointer p_closest_geometry = nullptr;
        for (IndexType i = 0; i < number_points_found; ++i) {
            auto p_point = points_found[i];
            auto p_segment = p_point->pGetObject();
            const double potential_min = GeometricalProjectionUtilities::FastMinimalDistanceOnLine(*p_segment, rNode, distance_tolerance);
            if (std::abs(potential_min) < std::abs(min_value)) {
                min_value = potential_min;
                p_closest_geometry = p_segment;
            }
        }

        GeometricalProjectionUtilities::FastMinimalDistanceOnLineWithRadius(min_value, *p_closest_geometry, rNode, radius_path, distance_tolerance);
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

        GeometricalProjectionUtilities::FastMinimalDistanceOnLineWithRadius(min_value, *p_closest_geometry, rNode, radius_path, distance_tolerance);
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