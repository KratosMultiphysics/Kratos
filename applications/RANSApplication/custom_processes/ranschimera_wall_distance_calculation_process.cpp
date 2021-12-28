//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//                   Rahul Kikkeri Nagaraja
//

// System includes
#include <limits>

// External includes

// Project includes
#include "includes/communicator.h"
#include "includes/variables.h"
#include "processes/find_global_nodal_neighbours_process.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/parallel_levelset_distance_calculator.h"
#include "utilities/variable_utils.h"

// Application includes

// Include base h
#include "ranschimera_wall_distance_calculation_process.h"

namespace Kratos
{
namespace RansChimeraWallDistanceCalculationUtilities
{
/**
 * @brief Calculates wall distance and updates nodal value
 *
 * This method calculates wall distance. Wall distance is the projection of distance vector
 * from node to wall point in the direction of normal
 *
 * @param rNode                 Node, where wall distance is required
 * @param rWallLocation         Nearest wall location
 * @param rUnitNormal           Wall normal (this should be always outward pointing normal)
 * @param rDistanceVariable     Distance variable to store nodal distance
 */
void CalculateAndUpdateNodalMinimumWallDistance(
    ModelPart::NodeType& rNode,
    const array_1d<double, 3>& rWallLocation,
    const array_1d<double, 3>& rUnitNormal,
    const Variable<double>& rDistanceVariable)
{
    // rUnitNormal is assumed to be outward pointing, hence wall_distance will be always positive.
    const double wall_distance =
        inner_prod(rWallLocation - rNode.Coordinates(), rUnitNormal);

    rNode.SetLock();
    double& current_distance = rNode.FastGetSolutionStepValue(rDistanceVariable);
    current_distance = std::min(current_distance, wall_distance);
    rNode.Set(VISITED, true);
    rNode.UnSetLock();
}
} // namespace RansChimeraWallDistanceCalculationUtilities

template <int TDim>
RansChimeraWallDistanceCalculationProcess<TDim>::RansChimeraWallDistanceCalculationProcess(
    Model& rModel,
    Parameters rParameters,
    ChimeraProcessType& rChimeraProcess)
    : mrModel(rModel), mrChimeraProcess(rChimeraProcess)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mMaxLevels = rParameters["max_levels"].GetInt();
    mEchoLevel = rParameters["echo_level"].GetInt();
    mMaxDistance = rParameters["max_distance"].GetDouble();
    mMainModelPartName = rParameters["main_model_part_name"].GetString();
    mWallModelPartName = rParameters["wall_model_part_name"].GetString();
    mDistanceVariableName = rParameters["distance_variable_name"].GetString();
    mNodalAreaVariableName = rParameters["nodal_area_variable_name"].GetString();
    mRecalculateAtEachTimeStep = rParameters["re_calculate_at_each_time_step"].GetBool();
    mIsFormulated = false;

    KRATOS_CATCH("");
}

template <int TDim>
int RansChimeraWallDistanceCalculationProcess<TDim>::Check()
{
    KRATOS_TRY

    const auto& r_model_part = mrModel.GetModelPart(mMainModelPartName);

    const auto& r_distance_variable =
        KratosComponents<Variable<double>>::Get(mDistanceVariableName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(r_distance_variable))
        << r_distance_variable.Name() << " is not found in nodal solution step variables list of "
        << mMainModelPartName << " used to store nodal distances.";

    const auto& r_nodal_area_variable = KratosComponents<Variable<double>>::Get(mNodalAreaVariableName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(r_nodal_area_variable))
        << r_nodal_area_variable.Name() << " is not found in nodal solution step variables list of "
        << mMainModelPartName << " used to store nodal areas.";

    return RansFormulationProcess::Check();

    KRATOS_CATCH("");
}

// void RansChimeraWallDistanceCalculationProcess::ExecuteInitialize()
// {
//     CalculateWallDistances();
// }

template <int TDim>
void RansChimeraWallDistanceCalculationProcess<TDim>::ExecuteInitializeSolutionStep()
{
    if (!mIsFormulated) {
        CalculateWallDistances();
        // CalculateAnalyticalWallDistances();
        mIsFormulated = true;
    }
}

template <int TDim>
void RansChimeraWallDistanceCalculationProcess<TDim>::ExecuteFinalizeSolutionStep()
{
    if (mRecalculateAtEachTimeStep){
        mIsFormulated = false;
    }
}

template <int TDim>
std::string RansChimeraWallDistanceCalculationProcess<TDim>::Info() const
{
    return std::string("RansChimeraWallDistanceCalculationProcess");
}

template <int TDim>
void RansChimeraWallDistanceCalculationProcess<TDim>::CalculateWallDistances()
{
    KRATOS_TRY

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Chimera wall distances calculation started for " << mMainModelPartName << std::endl;

    using namespace RansChimeraWallDistanceCalculationUtilities;

    auto& r_main_model_part = mrModel.GetModelPart(mMainModelPartName);
    auto& r_wall_model_part = mrModel.GetModelPart(mWallModelPartName);
    auto& r_communicator = r_main_model_part.GetCommunicator();

    const auto& r_distance_variable = KratosComponents<Variable<double>>::Get(mDistanceVariableName);
    const auto& r_nodal_area_variable = KratosComponents<Variable<double>>::Get(mNodalAreaVariableName);

    BuiltinTimer wall_distance_calculation_time;

    // variable initialization
    block_for_each(r_main_model_part.Nodes(), [&](ModelPart::NodeType& rNode){
        rNode.SetValue(NORMAL, NORMAL.Zero());
        rNode.Set(VISITED, false);
        rNode.FastGetSolutionStepValue(r_distance_variable) = mMaxDistance;
    });

    // identify all wall nodes to negative small distance values
    block_for_each(r_wall_model_part.Nodes(), [&](ModelPart::NodeType& rNode){
        rNode.FastGetSolutionStepValue(r_distance_variable) = -1e-12;
        rNode.Set(VISITED, true);
    });

    // calculate wall distances based on conditions
    block_for_each(r_wall_model_part.Conditions(), [&](ModelPart::ConditionType& rCondition) {
        array_1d<double, 3> normal = rCondition.GetValue(NORMAL);
        const double normal_magnitude = norm_2(normal);
        KRATOS_ERROR_IF(normal_magnitude == 0.0)
            << "NORMAL is not properly initialized in condition with id "
            << rCondition.Id() << " at " << rCondition.GetGeometry().Center() << ".\n";
        normal /= normal_magnitude;

        auto& parent_geometry = rCondition.GetValue(NEIGHBOUR_ELEMENTS)[0].GetGeometry();
        for (auto& r_parent_node : parent_geometry) {
            if (r_parent_node.FastGetSolutionStepValue(r_distance_variable) > 0.0) {
                CalculateAndUpdateNodalMinimumWallDistance(
                    r_parent_node, rCondition.GetGeometry().Center(), normal, r_distance_variable);
            }
        }

        // update nodal non-historical NORMAL to be used in case where
        // some nodes needs updating of nodal distances via Elements.
        auto& r_geometry = rCondition.GetGeometry();
        for (auto& r_node : r_geometry) {
            r_node.SetLock();
            r_node.GetValue(NORMAL) += normal;
            r_node.UnSetLock();
        }
    });

    // communication for all ranks
    r_communicator.AssembleNonHistoricalData(NORMAL);
    r_communicator.SynchronizeCurrentDataToMin(r_distance_variable);
    r_communicator.SynchronizeOrNodalFlags(VISITED);

    // calculate distances based on elements (on wall adjacent nodes which are not covered in NEIGHBOURING_ELEMENTS of conditions)
    block_for_each(r_main_model_part.Elements(), [&](ModelPart::ElementType& rElement) {
        auto& r_geometry = rElement.GetGeometry();

        array_1d<double, 3> normal = ZeroVector(3);
        array_1d<double, 3> normal_center = ZeroVector(3);
        double normals_count = 0.0;
        std::vector<int> nodal_indices_to_update;

        // identify nodes' distances which are not calculated by conditions, but are adjacent to wall nodes.
        // compute average normal and average wall location from nodal normals and nodal coordinates which
        // lies on the wall.
        for (std::size_t i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node) {
            const auto& r_node = r_geometry[i_node];
            if (r_node.FastGetSolutionStepValue(r_distance_variable) < 0.0) {
                noalias(normal) += r_node.GetValue(NORMAL);
                noalias(normal_center) += r_node.Coordinates();
                normals_count += 1.0;
            } else if (!r_node.Is(VISITED)) {
                nodal_indices_to_update.push_back(i_node);
            }
        }

        if (normals_count > 0.0 && nodal_indices_to_update.size() > 0) {
            const double normal_magnitude = norm_2(normal);
            KRATOS_ERROR_IF(normal_magnitude == 0.0)
                << "NORMAL is not properly initialized in adjacent wall nodes "
                   "in element with id "
                << rElement.Id() << " at " << rElement.GetGeometry().Center() << ".\n";

            normal /= normal_magnitude;
            normal_center /= normals_count;

            for (const auto i_node : nodal_indices_to_update) {
                auto& r_node = r_geometry[i_node];
                CalculateAndUpdateNodalMinimumWallDistance(
                    r_node, normal_center, normal, r_distance_variable);
            }
        }
    });

    // Transfer the negative distance values of patch wall nodes to the background subdomain
    VariableUtils().SetVariable(r_distance_variable, 0.0, r_main_model_part.Nodes(),
                                VISITED, false);

    // get the pointLocatorsMap from the chimera process
    PointLocatorsMapType const& mPointLocatorsMap = mrChimeraProcess.GetPointLocatorsMap();
    std::cout << "Point Locators size before start of Wall distance calculation: " << mPointLocatorsMap.size() << std::endl;

    // check if the modelparts in the mPointLocatorsMap are the submodelparts of r_main_model_part.
    for (auto const& rPointLocatorMap : mPointLocatorsMap){
        std::string modelpart_name = rPointLocatorMap.first;
        KRATOS_ERROR_IF(!(r_main_model_part.HasSubModelPart(modelpart_name)))
            << "Error: The patch modelpart " << modelpart_name << " in the mPointLocatorsMap was not found as the SubModelPart of "
            << mMainModelPartName << "." << std::endl;
    }

    for (auto& rBackgroundPointLocatorMap: mPointLocatorsMap){
        std::string background_mp_name = rBackgroundPointLocatorMap.first;
        auto& r_background_mp = r_main_model_part.GetSubModelPart(background_mp_name);
        auto const& p_point_locator_on_background = rBackgroundPointLocatorMap.second;

        for (auto& rPatchPointLocatorMap: mPointLocatorsMap){
            int found = 0;
            std::string patch_mp_name = rPatchPointLocatorMap.first;
            auto& r_patch_mp = r_main_model_part.GetSubModelPart(patch_mp_name);

            // if same modelpart pairs are encountered, move to the next pair of patch-background
            if (patch_mp_name == background_mp_name){
                continue;
            }

            // Transfer the negative distances of the patch nodes to their corresponding host elements in the background subdomain 
            for (auto& r_node_on_patch: r_patch_mp.Nodes()){
                if (r_node_on_patch.IsNot(BLOCKED)){
                    if (r_node_on_patch.Is(VISITED) && (r_node_on_patch.FastGetSolutionStepValue(r_distance_variable) < 0.0)){
                        Element::Pointer p_host_element;
                        Vector weights;
                        bool is_found = SearchNode(*p_point_locator_on_background, r_node_on_patch, p_host_element, weights);

                        if (is_found){
                            found += 1;
                            // data from the p_host_element
                            auto& rGeometry_of_host_element = p_host_element->GetGeometry();
                            for (auto& rNode_of_host_element: rGeometry_of_host_element){
                                CalculateAndUpdateNodalMinimumWallDistance( rNode_of_host_element,
                                                                            r_node_on_patch.Coordinates(),
                                                                            r_node_on_patch.GetValue(NORMAL), 
                                                                            r_distance_variable);
                                rNode_of_host_element.Set(BLOCKED, true); // so these nodes if they have negative dist are not considered as wall nodes henceforth
                                // not considering the minimum values of distances
                            }                
                        }
                    }                    
                }
            }
            std::cout << "bkg: " << background_mp_name << " and patch: " << patch_mp_name << std::endl;
            std::cout << "patch nodes: " << r_patch_mp.Nodes().size() << std::endl;
            std::cout << "overlapping nodes found: " << found << std::endl;
        }
    }

    // communication for all ranks
    r_communicator.SynchronizeCurrentDataToMin(r_distance_variable);
    r_communicator.SynchronizeOrNodalFlags(VISITED);

    // update rest of the domain
    const int domain_size = r_main_model_part.GetProcessInfo()[DOMAIN_SIZE];
    if (domain_size == 2) {
        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Parallel level set method started for " << r_main_model_part.Name() << std::endl;
        ParallelDistanceCalculator<2>().CalculateDistances(
            r_main_model_part, r_distance_variable, r_nodal_area_variable, mMaxLevels, mMaxDistance);
        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Parallel level set method completed for " << r_main_model_part.Name() << std::endl;
    } else if (domain_size == 3) {
        ParallelDistanceCalculator<3>().CalculateDistances(
            r_main_model_part, r_distance_variable, r_nodal_area_variable, mMaxLevels, mMaxDistance);
    } else {
        KRATOS_ERROR << "Unsupported domain size in " << r_main_model_part.Name()
                     << ". [ DOMAIN_SIZE = " << domain_size << " ]\n";
    }

    // revert boundary negative distances to zero
    VariableUtils().SetVariable(r_distance_variable, 0.0, r_wall_model_part.Nodes());
    // Expensive process: change the method.
    block_for_each(r_main_model_part.Nodes(), [&](ModelPart::NodeType& rNode){
        if (rNode.FastGetSolutionStepValue(r_distance_variable) < 0.0){
            rNode.FastGetSolutionStepValue(r_distance_variable) = 0.0;
        }
    });


    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Chimera wall distances calculation in " << mMainModelPartName << " took " 
        << wall_distance_calculation_time.ElapsedSeconds() << " seconds" << std::endl;

    KRATOS_CATCH("");
}

template <int TDim>
bool RansChimeraWallDistanceCalculationProcess<TDim>::SearchNode(PointLocatorType& rBinLocator,
                                                           NodeType& rNodeToFind,
                                                           Element::Pointer& rpHostElement,
                                                           Vector& rWeights)
{
    const int max_results = 10000;
    typename PointLocatorType::ResultContainerType results(max_results);
    typename PointLocatorType::ResultIteratorType result_begin = results.begin();

    bool is_found = false;
    is_found = rBinLocator.FindPointOnMesh(rNodeToFind.Coordinates(), rWeights,
                                           rpHostElement, result_begin, max_results);

    return is_found;
}

template <int TDim>
const Parameters RansChimeraWallDistanceCalculationProcess<TDim>::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "main_model_part_name"             : "PLEASE_SPECIFY_MAIN_MODEL_PART_NAME",
            "wall_model_part_name"             : "ALL_WALL_MODEL_PART",
            "max_levels"                       : 100,
            "max_distance"                     : 1e+30,
            "echo_level"                       : 0,
            "distance_variable_name"           : "WALL_DISTANCE",
            "nodal_area_variable_name"         : "NODAL_AREA",
            "re_calculate_at_each_time_step"   : false
        })");

    return default_parameters;
}


template <int TDim>
void RansChimeraWallDistanceCalculationProcess<TDim>::CalculateAnalyticalWallDistances(){

    KRATOS_TRY;

    auto& r_main_model_part = mrModel.GetModelPart(mMainModelPartName);
    auto& r_wall_model_part = mrModel.GetModelPart(mWallModelPartName);
    auto& r_communicator = r_main_model_part.GetCommunicator();

    // get the pointLocatorsMap from the chimera process
    PointLocatorsMapType const& mPointLocatorsMap = mrChimeraProcess.GetPointLocatorsMap();
    std::cout << "Point Locators size before start of Wall distance calculation: " << mPointLocatorsMap.size() << std::endl;

    const auto& r_distance_variable = KratosComponents<Variable<double>>::Get(mDistanceVariableName);
    // set radius and center
    const double radius = 0.05;
    Vector center(3);
    center[0] = 0.5;
    center[1] = 0.6;
    center[2] = 0.0;

    BuiltinTimer wall_distance_calculation_time;

    for (auto rPointLocatorMap : mPointLocatorsMap){
        std::string modelpart_name = rPointLocatorMap.first;
        auto& submodelpart = r_main_model_part.GetSubModelPart(modelpart_name);

        block_for_each(submodelpart.Nodes(), [&](ModelPart::NodeType& rNode){
            double distance_from_center = fabs(sqrt(
                std::pow((rNode.Coordinates()[0] - center[0]), 2) + 
                std::pow((rNode.Coordinates()[1] - center[1]), 2) +
                std::pow((rNode.Coordinates()[2] - center[2]), 2)            
            ));
            // std::cout << rNode.Id() << ": " << fabs(distance_from_center - radius) << std::endl;
            rNode.FastGetSolutionStepValue(r_distance_variable) = fabs(distance_from_center - radius);
            rNode.Set(VISITED, true);
        });
    }

    block_for_each(r_wall_model_part.Nodes(), [&](ModelPart::NodeType& rNode){
        rNode.FastGetSolutionStepValue(r_distance_variable) = 0.0;
    });

    r_communicator.SynchronizeCurrentDataToMin(r_distance_variable);
    r_communicator.SynchronizeOrNodalFlags(VISITED);

    KRATOS_INFO ("") << "Analytical Wall distances calculation took " 
        << wall_distance_calculation_time.ElapsedSeconds() << " seconds." << std::endl;

    KRATOS_CATCH("");
}

// Template declarations
template class RansChimeraWallDistanceCalculationProcess<2>;
template class RansChimeraWallDistanceCalculationProcess<3>;

} // namespace Kratos.
