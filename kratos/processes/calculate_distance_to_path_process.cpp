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

    // Local meshes
    auto& r_distance_local_mesh = r_comm.LocalMesh();
    auto& r_path_local_mesh = r_path_model_part.GetCommunicator().LocalMesh();
    
    // Initialize distance variable
    if constexpr ( THistorical) {
        VariableUtils().SetHistoricalVariableToZero(*mpDistanceVariable, r_distance_local_mesh.Nodes());
        r_comm.SynchronizeVariable(*mpDistanceVariable);
    } else {
        VariableUtils().SetNonHistoricalVariableToZero(*mpDistanceVariable, r_distance_local_mesh.Nodes());
        r_comm.SynchronizeNonHistoricalVariable(*mpDistanceVariable);
    }

    /* Getting a global vector of Line3D2N geometries */
    // First we check that every element in the modelpart is a line
    const auto& r_elements_array = r_distance_local_mesh.Elements();
    // By default we consider that the model part is composed of elements
    if (r_elements_array.size() > 0) {
        for (auto& r_elem : r_elements_array) {
            const auto& r_geom = r_elem.GetGeometry();
            const auto geom_type = r_geom.GetGeometryType();
            KRATOS_ERROR_IF((geom_type != GeometryData::KratosGeometryType::Kratos_Line2D2) || geom_type != GeometryData::KratosGeometryType::Kratos_Line3D2) << "The geometry type is not correct, it is suppossed to be a line" << std::endl;
        }
    } else {
        KRATOS_ERROR << "The model part is empty. Please check the model part provided" << std::endl;
    }

    // Now we create the vector of geometries
    const unsigned int size_vector = r_elements_array.size();
    if (r_distance_model_part.IsDistributed()) {
        r_data_comm.SumAll(size_vector);
    }
    std::vector<Geometry<NodeType>> geometries_vector(size_vector);
    // const auto it_elem_begin = r_elements_array.begin();
    // IndexPartition<std::size_t>(r_elements_array.size()).for_each([&](std::size_t i){
    //     const auto& r_geom = (it_elem_begin + i)->GetGeometry();
    //     geometries_vector[i] = Line3D2<NodeType>(r_geom);
    // });

    // Gettings if compute by brute force or not
    const bool brute_force_calculation = mThisParameters["brute_force_calculation"].GetBool();
    if (brute_force_calculation) {
        this->CalculateDistanceByBruteForce(r_distance_model_part, geometries_vector);
    } else {
        this->CalculateDistance(r_distance_model_part, geometries_vector);
    }

    // Synchronize variables
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
    std::vector<Geometry<NodeType>>& rVectorSegments
    )
{
    /// TODO
    const double radius_path = mThisParameters["radius_path"].GetDouble();
}

/***********************************************************************************/
/***********************************************************************************/

template<bool THistorical>
void CalculateDistanceToPathProcess<THistorical>::CalculateDistanceByBruteForce(
    ModelPart& rModelPart,
    std::vector<Geometry<NodeType>>& rVectorSegments
    )
{
    /// TODO
    const double radius_path = mThisParameters["radius_path"].GetDouble();
}

/***********************************************************************************/
/***********************************************************************************/

template<bool THistorical>
const Parameters CalculateDistanceToPathProcess<THistorical>::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"({
        "distance_model_part_name" :  "",
        "path_model_part_name"     :  "",
        "distance_variable"        : "DISTANCE",
        "brute_force_calculation"  : false,
        "radius_path"              : 0.0
    })" );

    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template class CalculateDistanceToPathProcess<true>;
template class CalculateDistanceToPathProcess<false>;

} /// namespace Kratos