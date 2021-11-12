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
    std::size_t n_neg (0);
    const std::size_t pts_number = rElementDistancesVector.size();
    for (std::size_t i_node = 0; i_node < pts_number; ++i_node){
        if (rElementDistancesVector[i_node] < 0.0)
            n_neg++;
    }
    return n_neg == pts_number;
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
            const auto& r_geom = rCondition.GetGeometry();
            return CalculateConditionFlowRate(r_geom);
        });
    }

    // Synchronize among processors
    flow_rate = r_communicator.GetDataCommunicator().SumAll(flow_rate);

    return flow_rate;
}

double FluidAuxiliaryUtilities::CalculateFlowRatePositiveSkin(const ModelPart& rModelPart)
{
    return CalculateFlowRateAuxiliary<true, false>(rModelPart);
}

double FluidAuxiliaryUtilities::CalculateFlowRateNegativeSkin(const ModelPart& rModelPart)
{
    return CalculateFlowRateAuxiliary<false, false>(rModelPart);
}

double FluidAuxiliaryUtilities::CalculateFlowRatePositiveSkin(
    const ModelPart& rModelPart,
    const Flags& rSkinFlag)
{
    return CalculateFlowRateAuxiliary<true, true>(rModelPart, rSkinFlag);
}

double FluidAuxiliaryUtilities::CalculateFlowRateNegativeSkin(
    const ModelPart& rModelPart,
    const Flags& rSkinFlag)
{
    return CalculateFlowRateAuxiliary<false, true>(rModelPart, rSkinFlag);
}

template<bool IsPositiveSubdomain, bool CheckConditionFlag>
double FluidAuxiliaryUtilities::CalculateFlowRateAuxiliary(
    const ModelPart& rModelPart,
    const Flags& rSkinFlag)
{
    // Check that there are conditions and distance variable in the nodal database
    KRATOS_ERROR_IF(rModelPart.GetCommunicator().GlobalNumberOfConditions() == 0) << "There are no conditions in the provided model part. Flow rate cannot be computed." << std::endl;
    const auto& r_communicator = rModelPart.GetCommunicator();
    if (r_communicator.LocalMesh().NumberOfNodes() !=0) {
        KRATOS_ERROR_IF_NOT(r_communicator.LocalMesh().NodesBegin()->SolutionStepsDataHas(DISTANCE)) << "Nodal solution step data has no \'DISTANCE\' variable. Flow rate cannot be computed" << std::endl;
        KRATOS_ERROR_IF_NOT(r_communicator.LocalMesh().NodesBegin()->SolutionStepsDataHas(VELOCITY)) << "Nodal solution step data has no \'VELOCITY\' variable. Flow rate cannot be computed" << std::endl;
    }

    double flow_rate = 0.0;
    if (r_communicator.LocalMesh().NumberOfConditions() != 0) {
        // Create the modified shape functions factory with the first condition parent as prototype
        const auto& r_cond_begin = r_communicator.LocalMesh().ConditionsBegin();
        const auto& r_parent_cond_begin = r_cond_begin->GetValue(NEIGHBOUR_ELEMENTS)[0];
        auto mod_sh_func_factory = GetStandardModifiedShapeFunctionsFactory(r_parent_cond_begin.GetGeometry());

        std::size_t n_dim = rModelPart.GetProcessInfo().GetValue(DOMAIN_SIZE);
        Vector nodal_distances(r_cond_begin->GetGeometry().PointsNumber());
        flow_rate = block_for_each<SumReduction<double>>(r_communicator.LocalMesh().Conditions(), nodal_distances, [&](Condition& rCondition, Vector& rNodalDistancesTLS){
            // Check if the condition is to be added to the flow contribution
            double cond_flow_rate = 0.0;
            if (CheckConditionFlagAuxiliary<CheckConditionFlag>(rCondition, rSkinFlag)) {
                // Get geometry data
                const auto& r_geom = rCondition.GetGeometry();
                const std::size_t n_nodes = r_geom.PointsNumber();

                // Set up distances vector
                for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                    rNodalDistancesTLS(i_node) = r_geom[i_node].FastGetSolutionStepValue(DISTANCE);
                }
                // Check if the condition is in the positive subdomain or intersected
                if (CheckNonSplitConditionSubdomain<IsPositiveSubdomain>(rNodalDistancesTLS)) {
                    cond_flow_rate = CalculateConditionFlowRate(r_geom);
                } else if (IsSplit(rNodalDistancesTLS)){
                    // Get the current condition parent
                    const auto p_parent_element = rCondition.GetValue(NEIGHBOUR_ELEMENTS)(0);
                    KRATOS_ERROR_IF_NOT(p_parent_element.get()) << "Condition " << rCondition.Id() << " has no parent element assigned." << std::endl;
                    const auto& r_parent_geom = p_parent_element->GetGeometry();

                    // Get the corresponding face id of the current condition
                    const std::size_t n_parent_faces = r_parent_geom.FacesNumber();
                    DenseMatrix<unsigned int> nodes_in_faces(n_parent_faces, n_parent_faces);
                    r_parent_geom.NodesInFaces(nodes_in_faces);
                    std::size_t face_id = 0;
                    for (std::size_t i_face = 0; i_face < n_parent_faces; ++i_face) {
                        std::size_t match_nodes = 0;
                        for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                            std::size_t parent_local_id = nodes_in_faces(i_node + 1, i_face);
                            for (std::size_t j_node = 0; j_node < n_nodes; ++j_node) {
                                if (r_geom[j_node].Id() == r_parent_geom[parent_local_id].Id()) {
                                    match_nodes++;
                                    break;
                                }
                            }
                        }
                        if (match_nodes == n_nodes) {
                            face_id = i_face;
                            break;
                        }
                    }

                    // Calculate the modified shape functions in the face of interest
                    const std::size_t n_nodes_parent = r_parent_geom.PointsNumber();
                    Vector parent_distances(n_nodes_parent);
                    for (std::size_t i_node = 0; i_node < n_nodes_parent; ++i_node) {
                        parent_distances(i_node) = r_parent_geom[i_node].FastGetSolutionStepValue(DISTANCE);
                    }
                    auto p_mod_sh_func = mod_sh_func_factory(p_parent_element->pGetGeometry(), parent_distances);

                    Vector w_vect;
                    Matrix N_container;
                    std::vector<array_1d<double,3>> normals_vect;
                    CalculateSplitConditionGeometryData<IsPositiveSubdomain>(p_mod_sh_func, face_id, N_container, normals_vect, w_vect);

                    // Interpolate the flow rate in the positive subdomain
                    Vector i_normal(n_dim);
                    Vector i_N(n_nodes_parent);
                    array_1d<double,3> aux_vel;
                    const std::size_t n_gauss = w_vect.size();
                    for (std::size_t i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
                        aux_vel = ZeroVector(3);
                        i_N = row(N_container, i_gauss);
                        i_normal = normals_vect[i_gauss];
                        const double i_normal_norm = norm_2(i_normal);
                        KRATOS_WARNING_IF("CalculateFlowRatePositiveSkin", i_normal_norm < 1.0e-12) << "Condition " << rCondition.Id() << " normal close to zero." << std::endl;
                        i_normal /= i_normal_norm;
                        for (std::size_t i_parent_node = 0; i_parent_node < n_nodes_parent; ++i_parent_node) {
                            noalias(aux_vel) += i_N[i_parent_node] * r_parent_geom[i_parent_node].FastGetSolutionStepValue(VELOCITY);
                        }
                        cond_flow_rate += w_vect[i_gauss] * inner_prod(aux_vel, i_normal);
                    }
                }
            }

            return cond_flow_rate;
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

double FluidAuxiliaryUtilities::CalculateConditionFlowRate(const GeometryType& rGeometry)
{
    // Calculate the condition area normal
    GeometryType::CoordinatesArrayType point_local;
    rGeometry.PointLocalCoordinates(point_local, rGeometry.Center()) ;
    const array_1d<double,3> area_normal = rGeometry.Normal(point_local);
    // Check condition area and calculate the condition average flow rate
    double condition_flow_rate = 0.0;
    if (norm_2(area_normal) > std::numeric_limits<double>::epsilon()) {
        for (auto& r_node : rGeometry) {
            condition_flow_rate += MathUtils<double>::Dot(r_node.FastGetSolutionStepValue(VELOCITY), area_normal);
        }
        condition_flow_rate /= static_cast<double>(rGeometry.PointsNumber());
    } else {
        KRATOS_WARNING("CalculateFlowRate") << "Condition area is close to zero. Flow rate not considered." << std::endl;
    }
    return condition_flow_rate;
}

template<>
bool FluidAuxiliaryUtilities::CheckConditionFlagAuxiliary<true>(
    const Condition& rCondition,
    const Flags& rSkinFlag)
{
    return rCondition.Is(rSkinFlag);
}

template<>
bool FluidAuxiliaryUtilities::CheckConditionFlagAuxiliary<false>(
    const Condition& rCondition,
    const Flags& rSkinFlag)
{
    return true;
}

template<>
bool FluidAuxiliaryUtilities::CheckNonSplitConditionSubdomain<true>(const Vector &rElementDistancesVector)
{
    return IsPositive(rElementDistancesVector);
}

template<>
bool FluidAuxiliaryUtilities::CheckNonSplitConditionSubdomain<false>(const Vector &rElementDistancesVector)
{
    return IsNegative(rElementDistancesVector);
}

template<>
void FluidAuxiliaryUtilities::CalculateSplitConditionGeometryData<true>(
    const ModifiedShapeFunctions::UniquePointer& rpModShapeFunc,
    const std::size_t FaceId,
    Matrix& rShapeFunctions,
    std::vector<array_1d<double,3>>& rNormals,
    Vector& rWeights)
{
    //TODO: Use a method without gradients when we implement it
    ModifiedShapeFunctions::ShapeFunctionsGradientsType n_pos_DN_DX;
    rpModShapeFunc->ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(rShapeFunctions, n_pos_DN_DX, rWeights, FaceId, GeometryData::IntegrationMethod::GI_GAUSS_2);
    rpModShapeFunc->ComputePositiveExteriorFaceAreaNormals(rNormals, FaceId, GeometryData::IntegrationMethod::GI_GAUSS_2);
}

template<>
void FluidAuxiliaryUtilities::CalculateSplitConditionGeometryData<false>(
    const ModifiedShapeFunctions::UniquePointer& rpModShapeFunc,
    const std::size_t FaceId,
    Matrix& rShapeFunctions,
    std::vector<array_1d<double,3>>& rNormals,
    Vector& rWeights)
{
    //TODO: Use a method without gradients when we implement it
    ModifiedShapeFunctions::ShapeFunctionsGradientsType n_pos_DN_DX;
    rpModShapeFunc->ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(rShapeFunctions, n_pos_DN_DX, rWeights, FaceId, GeometryData::IntegrationMethod::GI_GAUSS_2);
    rpModShapeFunc->ComputeNegativeExteriorFaceAreaNormals(rNormals, FaceId, GeometryData::IntegrationMethod::GI_GAUSS_2);
}

} // namespace Kratos
