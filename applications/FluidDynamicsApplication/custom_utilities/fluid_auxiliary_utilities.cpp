//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Suneth Warnakulasuriya
//

// System includes


// External includes


// Project includes
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "fluid_auxiliary_utilities.h"

namespace Kratos
{

bool FluidAuxiliaryUtilities::IsSplit(const Vector& rElementDistancesVector)
{
    std::size_t n_pos(0);
    std::size_t n_neg(0);
    const std::size_t pts_number = rElementDistancesVector.size();
    for (std::size_t i_node = 0; i_node < pts_number; ++i_node){
        if (rElementDistancesVector[i_node] > 0.0)
            n_pos++;
        else
            n_neg++;
    }
    return (n_pos > 0 && n_neg > 0) ? true : false;
}

bool FluidAuxiliaryUtilities::IsPositive(const Vector &rElementDistancesVector)
{
    std::size_t n_pos (0);
    const std::size_t pts_number = rElementDistancesVector.size();
    for (std::size_t i_node = 0; i_node < pts_number; ++i_node){
        if (rElementDistancesVector[i_node] > 0.0)
            n_pos++;
    }
    return (n_pos == pts_number) ? true : false;
}

bool FluidAuxiliaryUtilities::IsNegative(const Vector &rElementDistancesVector)
{
    return !IsPositive(rElementDistancesVector);
}

double FluidAuxiliaryUtilities::CalculateFluidVolume(const ModelPart& rModelPart)
{
    KRATOS_ERROR_IF(rModelPart.GetCommunicator().GlobalNumberOfElements() == 0) << "There are no elements in the provided model part. Fluid volume cannot be computed." << std::endl;

    double fluid_volume = 0.0;
    const auto& r_communicator = rModelPart.GetCommunicator();
    if (r_communicator.LocalMesh().NumberOfElements() != 0) {
        fluid_volume = block_for_each<SumReduction<double>>(r_communicator.LocalMesh().Elements(), [](Element& rElement){
            const auto& r_geom = rElement.GetGeometry();
            return r_geom.DomainSize();
        });
    }
    r_communicator.GetDataCommunicator().SumAll(fluid_volume);

    return fluid_volume;
}

double FluidAuxiliaryUtilities::CalculateFluidPositiveVolume(const ModelPart& rModelPart)
{
    // Check that there are elements and distance variable in the nodal database
    KRATOS_ERROR_IF(rModelPart.GetCommunicator().GlobalNumberOfElements() == 0) << "There are no elements in the provided model part. Fluid volume cannot be computed." << std::endl;
    const auto& r_communicator = rModelPart.GetCommunicator();
    if (r_communicator.LocalMesh().NumberOfNodes() !=0) {
        KRATOS_ERROR_IF_NOT(r_communicator.LocalMesh().NodesBegin()->SolutionStepsDataHas(DISTANCE)) << "Nodal solution step data has no \'DISTANCE\' variable. Positive volume cannot be computed" << std::endl;
    }

    double fluid_volume = 0.0;
    if (r_communicator.LocalMesh().NumberOfElements() != 0) {
        // Set the modified shape functions fatory
        const auto& r_geom_begin = r_communicator.LocalMesh().ElementsBegin()->GetGeometry();
        auto mod_sh_func_factory = GetStandardModifiedShapeFunctionsFactory(r_geom_begin);

        // Calculate the positive volume
        Vector nodal_distances(r_geom_begin.PointsNumber());
        fluid_volume = block_for_each<SumReduction<double>>(rModelPart.GetCommunicator().LocalMesh().Elements(), nodal_distances, [&mod_sh_func_factory](Element& rElement, Vector& rNodalDistancesTLS){
            // Set the distances vector to check if the element is split
            const auto& r_geom = rElement.GetGeometry();
            const std::size_t n_nodes = r_geom.PointsNumber();
            for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                rNodalDistancesTLS[i_node] = r_geom[i_node].FastGetSolutionStepValue(DISTANCE);
            }
            // Split element check
            double elem_volume = 0.0;
            if (IsSplit(rNodalDistancesTLS)) {
                // Compute positive volume fraction with the modified shape functions
                auto p_mod_sh_func = mod_sh_func_factory(rElement.pGetGeometry(), rNodalDistancesTLS);
                elem_volume = p_mod_sh_func->ComputePositiveSideDomainSize();
            } else if (IsPositive(rNodalDistancesTLS)) {
                // If the element is positive, compute the volume from geometry
                elem_volume = r_geom.DomainSize();
            }

            // Return the value to be reduced
            return elem_volume;
        });
    }

    // Synchronize among processors
    r_communicator.GetDataCommunicator().SumAll(fluid_volume);

    return fluid_volume;
}

double FluidAuxiliaryUtilities::CalculateFluidNegativeVolume(const ModelPart& rModelPart)
{
    // Check that there are elements and distance variable in the nodal database
    KRATOS_ERROR_IF(rModelPart.GetCommunicator().GlobalNumberOfElements() == 0) << "There are no elements in the provided model part. Fluid volume cannot be computed." << std::endl;
    const auto& r_communicator = rModelPart.GetCommunicator();
    if (r_communicator.LocalMesh().NumberOfNodes() !=0) {
        KRATOS_ERROR_IF_NOT(r_communicator.LocalMesh().NodesBegin()->SolutionStepsDataHas(DISTANCE)) << "Nodal solution step data has no \'DISTANCE\' variable. Negative volume cannot be computed" << std::endl;
    }

    double fluid_volume = 0.0;
    if (r_communicator.LocalMesh().NumberOfElements() != 0) {
        // Set the modified shape functions fatory
        const auto& r_geom_begin = r_communicator.LocalMesh().ElementsBegin()->GetGeometry();
        auto mod_sh_func_factory = GetStandardModifiedShapeFunctionsFactory(r_geom_begin);

        // Calculate the negative volume
        Vector nodal_distances(r_geom_begin.PointsNumber());
        fluid_volume = block_for_each<SumReduction<double>>(rModelPart.GetCommunicator().LocalMesh().Elements(), nodal_distances, [&mod_sh_func_factory](Element& rElement, Vector& rNodalDistancesTLS){
            // Set the distances vector to check if the element is split
            const auto& r_geom = rElement.GetGeometry();
            const std::size_t n_nodes = r_geom.PointsNumber();
            for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                rNodalDistancesTLS[i_node] = r_geom[i_node].FastGetSolutionStepValue(DISTANCE);
            }
            // Split element check
            double elem_volume = 0.0;
            if (IsSplit(rNodalDistancesTLS)) {
                // Compute negative volume fraction with the modified shape functions
                auto p_mod_sh_func = mod_sh_func_factory(rElement.pGetGeometry(), rNodalDistancesTLS);
                elem_volume = p_mod_sh_func->ComputeNegativeSideDomainSize();
            } else if (IsNegative(rNodalDistancesTLS)) {
                // If the element is negative, compute the volume from geometry
                elem_volume = r_geom.DomainSize();
            }

            // Return the value to be reduced
            return elem_volume;
        });
    }

    // Synchronize among processors
    r_communicator.GetDataCommunicator().SumAll(fluid_volume);

    return fluid_volume;
}

double FluidAuxiliaryUtilities::CalculateFlowRate(const ModelPart& rModelPart)
{
    // Check that there are conditions and distance variable in the nodal database
    KRATOS_ERROR_IF(rModelPart.GetCommunicator().GlobalNumberOfConditions() == 0) << "There are no conditions in the provided model part. Flow rate cannot be computed." << std::endl;
    const auto& r_communicator = rModelPart.GetCommunicator();
    if (r_communicator.LocalMesh().NumberOfNodes() !=0) {
        KRATOS_ERROR_IF_NOT(r_communicator.LocalMesh().NodesBegin()->SolutionStepsDataHas(VELOCITY)) << "Nodal solution step data has no \'VELOCITY\' variable. Flow rate cannot be computed" << std::endl;
    }

    double flow_rate = 0.0;
    if (r_communicator.LocalMesh().NumberOfConditions() != 0) {
        flow_rate = block_for_each<SumReduction<double>>(r_communicator.LocalMesh().Conditions(), [](Condition& rCondition){
            // Calculate the condition area (length in 2D)
            const auto& r_geom = rCondition.GetGeometry();
            // Calculate the condition area normal
            GeometryType::CoordinatesArrayType point_local;
            r_geom.PointLocalCoordinates(point_local, r_geom.Center()) ;
            const array_1d<double,3> area_normal = r_geom.Normal(point_local);
            // Check condition area and calculate the condition average flow rate
            double condition_flow_rate = 0.0;
            if (norm_2(area_normal) > std::numeric_limits<double>::epsilon()) {
                for (auto& r_node : r_geom) {
                    condition_flow_rate += MathUtils<double>::Dot(r_node.FastGetSolutionStepValue(VELOCITY), area_normal);
                }
                condition_flow_rate /= static_cast<double>(r_geom.PointsNumber());
            } else {
                KRATOS_WARNING("CalculateFlowRate") << "Condition " << rCondition.Id() << " area is close to zero. Flow rate not considered." << std::endl;
            }
            return condition_flow_rate;
        });
    }

    // Synchronize among processors
    flow_rate = r_communicator.GetDataCommunicator().SumAll(flow_rate);

    return flow_rate;
}

FluidAuxiliaryUtilities::ModifiedShapeFunctionsFactoryType FluidAuxiliaryUtilities::GetStandardModifiedShapeFunctionsFactory(const GeometryType& rGeometry)
{
    switch (rGeometry.GetGeometryType()) {
        case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
            return [](const GeometryType::Pointer pGeometry, const Vector& rNodalDistances)->ModifiedShapeFunctions::UniquePointer{
                return Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(pGeometry, rNodalDistances);};
        case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
            return [](const GeometryType::Pointer pGeometry, const Vector& rNodalDistances)->ModifiedShapeFunctions::UniquePointer{
                return Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(pGeometry, rNodalDistances);};
        default:
            KRATOS_ERROR << "Asking for a non-implemented modified shape functions geometry.";
    }
}

} // namespace Kratos
