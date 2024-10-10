//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//                   Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "utilities/normal_calculation_utils.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
template <class TContainerType>
void NormalCalculationUtils::CalculateNormalsInContainer(
    ModelPart& rModelPart,
    const Variable<array_1d<double,3>>& rNormalVariable
    )
{
    KRATOS_TRY;

    // Declare auxiliar coordinates
    Point::CoordinatesArrayType aux_coords;

    auto& r_entity_array = GetContainer<TContainerType>(rModelPart);
    const auto it_entity_begin = r_entity_array.begin();

    // TODO: Change to OMP parallel utilities
#pragma omp parallel for firstprivate(aux_coords)
    for (int i = 0; i < static_cast<int>(r_entity_array.size()); ++i) {
        auto it_entity = it_entity_begin + i;
        const auto& r_geometry = it_entity->GetGeometry();

        // Avoid not "flat" conditions
        if (r_geometry.WorkingSpaceDimension() != r_geometry.LocalSpaceDimension() + 1) {
            continue;
        }

        // Set entity normal
        r_geometry.PointLocalCoordinates(aux_coords, r_geometry.Center());
        it_entity->SetValue(rNormalVariable, r_geometry.UnitNormal(aux_coords));
    }

    KRATOS_CATCH("Error in CalculateNormalsInContainer");
}

/***********************************************************************************/
/***********************************************************************************/

template <class TContainerType, bool TIsHistorical>
void NormalCalculationUtils::InitializeNormals(
    ModelPart& rModelPart,
    const Variable<array_1d<double,3>>& rNormalVariable
    )
{
    KRATOS_TRY

    // Resetting the normals
    const array_1d<double, 3> zero = ZeroVector(3);

    auto& r_entity_array = GetContainer<TContainerType>(rModelPart);

    if (rModelPart.GetCommunicator().GetDataCommunicator().IsDistributed()) {
        // If Parallel make sure normals are reset in all partitions
        VariableUtils().SetFlag(VISITED, false, rModelPart.Nodes());

        block_for_each(r_entity_array, [](typename TContainerType::value_type& rEntity) {
            for (auto& r_node : rEntity.GetGeometry()) {
                r_node.SetLock();
                r_node.Set(VISITED, true);
                r_node.UnSetLock();
            }
        });

        rModelPart.GetCommunicator().SynchronizeOrNodalFlags(VISITED);
        if (TIsHistorical) {
            VariableUtils().SetVariable(rNormalVariable, zero, rModelPart.Nodes(), VISITED);
        } else {
            VariableUtils().SetNonHistoricalVariable(rNormalVariable, zero, rModelPart.Nodes(), VISITED);
        }
    } else {
        // In serial iterate normally over the condition nodes
        block_for_each(r_entity_array, [zero, &rNormalVariable, this](typename TContainerType::value_type& rEntity) {
                for (auto& r_node : rEntity.GetGeometry()) {
                    r_node.SetLock();
                    SetNormalValue<TIsHistorical>(r_node, rNormalVariable, zero);
                    r_node.UnSetLock();
                }
        });
    }

    KRATOS_CATCH("Error in InitializeNormals");
}

/***********************************************************************************/
/***********************************************************************************/

template <>
KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::CalculateNormals<ModelPart::ConditionsContainerType, true>(
    ModelPart& rModelPart,
    const bool EnforceGenericGeometryAlgorithm,
    const bool ConsiderUnitNormal,
    const Variable<array_1d<double,3>>& rNormalVariable
    )
{
    KRATOS_TRY

    if (CheckUseSimplex(rModelPart, EnforceGenericGeometryAlgorithm)) {
        const auto& r_process_info = rModelPart.GetProcessInfo();
        const SizeType dimension = r_process_info.GetValue(DOMAIN_SIZE);
        CalculateOnSimplex(rModelPart, dimension, rNormalVariable);
    } else {
        CalculateNormalsUsingGenericAlgorithm<ModelPart::ConditionsContainerType, true>(rModelPart, ConsiderUnitNormal, rNormalVariable);
    }

    KRATOS_CATCH("Error in CalculateNormals");
}

/***********************************************************************************/
/***********************************************************************************/

template <>
KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::CalculateNormals<ModelPart::ConditionsContainerType, false>(
    ModelPart& rModelPart,
    const bool EnforceGenericGeometryAlgorithm,
    const bool ConsiderUnitNormal,
    const Variable<array_1d<double,3>>& rNormalVariable
    )
{
    if (CheckUseSimplex(rModelPart, EnforceGenericGeometryAlgorithm)) {
        const auto& r_process_info = rModelPart.GetProcessInfo();
        const SizeType dimension = r_process_info.GetValue(DOMAIN_SIZE);
        CalculateOnSimplexNonHistorical(rModelPart, dimension, rNormalVariable);
    } else {
        CalculateNormalsUsingGenericAlgorithm<ModelPart::ConditionsContainerType, false>(rModelPart, ConsiderUnitNormal, rNormalVariable);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::CalculateNormals<ModelPart::ElementsContainerType, true>(
    ModelPart& rModelPart,
    const bool EnforceGenericGeometryAlgorithm,
    const bool ConsiderUnitNormal,
    const Variable<array_1d<double,3>>& rNormalVariable
    )
{
    CalculateNormalsUsingGenericAlgorithm<ModelPart::ElementsContainerType, true>(rModelPart, ConsiderUnitNormal, rNormalVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::CalculateNormals<ModelPart::ElementsContainerType, false>(
    ModelPart& rModelPart,
    const bool EnforceGenericGeometryAlgorithm,
    const bool ConsiderUnitNormal,
    const Variable<array_1d<double,3>>& rNormalVariable)
{
    CalculateNormalsUsingGenericAlgorithm<ModelPart::ElementsContainerType, false>(rModelPart, ConsiderUnitNormal, rNormalVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TContainerType, bool TIsHistorical>
void NormalCalculationUtils::CalculateUnitNormals(
    ModelPart& rModelPart,
    const bool EnforceGenericGeometryAlgorithm,
    const Variable<array_1d<double,3>>& rNormalVariable
    )
{
    // Compute area normals
    CalculateNormals<TContainerType, TIsHistorical>(rModelPart, EnforceGenericGeometryAlgorithm, true, rNormalVariable);

    // Compute unit normals
    ComputeUnitNormalsFromAreaNormals<TContainerType, TIsHistorical>(rModelPart, rNormalVariable);
}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::CalculateOnSimplex(
    ConditionsArrayType& rConditions,
    const std::size_t Dimension,
    const Variable<array_1d<double,3>>& rNormalVariable
    )
{
    KRATOS_TRY

    AuxiliaryCalculateOnSimplex<true>(rConditions, Dimension, rNormalVariable);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::CalculateOnSimplexNonHistorical(
    ConditionsArrayType& rConditions,
    const std::size_t Dimension,
    const Variable<array_1d<double,3>>& rNormalVariable
    )
{
    KRATOS_TRY

    AuxiliaryCalculateOnSimplex<false>(rConditions, Dimension, rNormalVariable);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::CalculateNormalShapeDerivativesOnSimplex(
    ConditionsArrayType& rConditions,
    const std::size_t Dimension
    )
{
    KRATOS_TRY

    if (Dimension == 2) {
        // reset NORMAL_SHAPE_DERIVATIVE variable for all conditions
        VariableUtils().SetNonHistoricalVariable(NORMAL_SHAPE_DERIVATIVE,
                                                 ZeroMatrix(4, 2), rConditions);

        // calculate condition normal shape derivatives
        block_for_each(rConditions, [](ConditionType& rCondition) {
            if (rCondition.GetGeometry().GetGeometryType() ==
                GeometryData::KratosGeometryType::Kratos_Line2D2) {
                CalculateNormalShapeDerivative2D(rCondition);
            }
        });
    } else {
        // reset NORMAL_SHAPE_DERIVATIVE variable for all conditions
        VariableUtils().SetNonHistoricalVariable(NORMAL_SHAPE_DERIVATIVE,
                                                 ZeroMatrix(9, 3), rConditions);

        // calculate condition normal shape derivatives
        block_for_each(rConditions, [](ConditionType& rCondition) {
            if (rCondition.GetGeometry().GetGeometryType() ==
                GeometryData::KratosGeometryType::Kratos_Triangle3D3) {
                CalculateNormalShapeDerivative3D(rCondition);
            }
        });
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::CalculateOnSimplex(
    ModelPart& rModelPart,
    const Variable<array_1d<double,3>>& rNormalVariable
    )
{
    KRATOS_TRY

    const auto& r_process_info = rModelPart.GetProcessInfo();

    KRATOS_ERROR_IF(!r_process_info.Has(DOMAIN_SIZE))
        << "DOMAIN_SIZE not found in process info of " << rModelPart.Name() << ".";
    CalculateOnSimplex(rModelPart, r_process_info[DOMAIN_SIZE], rNormalVariable);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::CalculateOnSimplexNonHistorical(
    ModelPart& rModelPart,
    const Variable<array_1d<double,3>>& rNormalVariable)
{
    KRATOS_TRY

    const auto& r_process_info = rModelPart.GetProcessInfo();

    KRATOS_ERROR_IF(!r_process_info.Has(DOMAIN_SIZE))
        << "DOMAIN_SIZE not found in process info of " << rModelPart.Name() << ".";
    CalculateOnSimplexNonHistorical(rModelPart, r_process_info[DOMAIN_SIZE], rNormalVariable);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::CalculateOnSimplex(
    ModelPart& rModelPart,
    const std::size_t Dimension,
    const Variable<array_1d<double,3>>& rNormalVariable
    )
{
    KRATOS_TRY

    // Initialize the normals
    InitializeNormals<ModelPart::ConditionsContainerType,true>(rModelPart, rNormalVariable);

    // Calling CalculateOnSimplex for conditions
    AuxiliaryCalculateOnSimplex<true>(rModelPart.Conditions(), Dimension, rNormalVariable);

    // Synchronize the normal
    rModelPart.GetCommunicator().AssembleCurrentData(rNormalVariable);

    KRATOS_CATCH("Error in CalculateOnSimplex");
}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::CalculateOnSimplexNonHistorical(
    ModelPart& rModelPart,
    const std::size_t Dimension,
    const Variable<array_1d<double,3>>& rNormalVariable)
{
    // Initialize the normals
    InitializeNormals<ModelPart::ConditionsContainerType,false>(rModelPart, rNormalVariable);

    // Calling CalculateOnSimplex for conditions
    AuxiliaryCalculateOnSimplex<false>(rModelPart.Conditions(), Dimension, rNormalVariable);

    // Synchronize the normal
    rModelPart.GetCommunicator().AssembleNonHistoricalData(rNormalVariable);
}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::SwapNormals(ModelPart& rModelPart)
{
    KRATOS_TRY

    for(auto& r_cond : rModelPart.Conditions()) {
        GeometryType& r_geometry = r_cond.GetGeometry();
        Node::Pointer paux = r_geometry(0);
        r_geometry(0) = r_geometry(1);
        r_geometry(1) = paux;
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<class TContainerType, bool TIsHistorical>
void NormalCalculationUtils::ComputeUnitNormalsFromAreaNormals(
    ModelPart& rModelPart,
    const Variable<array_1d<double,3>>& rNormalVariable)
{
    KRATOS_TRY

    // Auxiliary lambda to calculate the unit normal from the area one
    auto unit_normal_func = [&rNormalVariable, this](NodeType& rNode){
        auto& r_normal = GetNormalValue<TIsHistorical>(rNode, rNormalVariable);
        const double norm_normal = norm_2(r_normal);
        KRATOS_ERROR_IF(rNode.Is(INTERFACE) && norm_normal <= std::numeric_limits<double>::epsilon())
            << "Zero norm normal in INTERFACE flagged node: " << rNode.Id() << std::endl;
        r_normal /= norm_normal;
    };

    // We iterate over the nodes of the corresponding container
    // Note that we need to iterate the corresponding entities in case there are
    // nodes in the model part that do not belong to any entity.
    auto& r_entity_container = GetContainer<TContainerType>(rModelPart);
    if (rModelPart.GetCommunicator().GetDataCommunicator().IsDistributed()) {
        // Flag the nodes to which we will compute the unit normal in all the partitions
        VariableUtils().SetFlag(VISITED, false, rModelPart.Nodes());
        block_for_each(r_entity_container, [](typename TContainerType::value_type& rEntity) {
            for (auto& r_node : rEntity.GetGeometry()) {
                r_node.SetLock();
                r_node.Set(VISITED, true);
                r_node.UnSetLock();
            }
        });
        rModelPart.GetCommunicator().SynchronizeOrNodalFlags(VISITED);

        // For the flagged nodes in the local mesh calculate the unit normal
        auto& r_local_nodes = rModelPart.GetCommunicator().LocalMesh().Nodes();
        block_for_each(r_local_nodes, [&unit_normal_func](NodeType& rNode) {
            if (rNode.Is(VISITED)) {
                unit_normal_func(rNode);
            }
        });

        // Synchronize unit normals to ghost nodes
        if (TIsHistorical) {
            rModelPart.GetCommunicator().SynchronizeVariable(rNormalVariable);
        } else {
            rModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rNormalVariable);
        }
    } else {
        // In serial iterate normally over the container nodes
        // Note that we only calculate the unit normal in the nodes belonging to the entities of interest
        block_for_each(r_entity_container, [&unit_normal_func](typename TContainerType::value_type& rEntity) {
            for (auto& r_node : rEntity.GetGeometry()) {
                r_node.SetLock();
                unit_normal_func(r_node);
                r_node.UnSetLock();
            }
        });
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::CalculateNormal2D(
    Condition& rCondition,
    array_1d<double,3>& rAn,
    const Variable<array_1d<double,3>>& rNormalVariable
    )
{
    const GeometryType& r_geometry = rCondition.GetGeometry();

    rAn[0] =    r_geometry[1].Y() - r_geometry[0].Y();
    rAn[1] = - (r_geometry[1].X() - r_geometry[0].X());
    rAn[2] =    0.00;

    rCondition.SetValue(rNormalVariable, rAn);
}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::CalculateNormalShapeDerivative2D(
    ConditionType& rCondition
    )
{
    KRATOS_TRY

    // since this is a private static method (no external calls),
    // this method assumes matrix is properly sized always before this call.
    Matrix& r_normal_shape_derivatives = rCondition.GetValue(NORMAL_SHAPE_DERIVATIVE);

    KRATOS_DEBUG_ERROR_IF(r_normal_shape_derivatives.size1() != 4 || r_normal_shape_derivatives.size2() != 2)
        << "Condition (id: " << rCondition.Id() << ") does not have properly assigned \"NORMAL_SHAPE_DERIVATIVE\" "
        << "variable. Please assigned this variable with proper matrix size (required matrix size: (4, 2) != given matrix size ("
        << r_normal_shape_derivatives.size1() << ", " << r_normal_shape_derivatives.size2() << ")).";

    r_normal_shape_derivatives(0, 0) = 0.0;  // n_x derivative w.r.t. node_0_x
    r_normal_shape_derivatives(1, 0) = -1.0; // n_x derivative w.r.t. node_0_y
    r_normal_shape_derivatives(2, 0) = 0.0;  // n_x derivative w.r.t. node_1_x
    r_normal_shape_derivatives(3, 0) = 1.0;  // n_x derivative w.r.t. node_1_y

    r_normal_shape_derivatives(0, 1) = 1.0;  // n_y derivative w.r.t. node_0_x
    r_normal_shape_derivatives(1, 1) = 0.0;  // n_y derivative w.r.t. node_0_y
    r_normal_shape_derivatives(2, 1) = -1.0; // n_y derivative w.r.t. node_1_x
    r_normal_shape_derivatives(3, 1) = 0.0;  // n_y derivative w.r.t. node_1_y

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::CalculateNormal3D(
    Condition& rCondition,
    array_1d<double,3>& rAn,
    array_1d<double,3>& rv1,
    array_1d<double,3>& rv2,
    const Variable<array_1d<double,3>>& rNormalVariable
    )
{
    const GeometryType& r_geometry = rCondition.GetGeometry();

    rv1[0] = r_geometry[1].X() - r_geometry[0].X();
    rv1[1] = r_geometry[1].Y() - r_geometry[0].Y();
    rv1[2] = r_geometry[1].Z() - r_geometry[0].Z();

    rv2[0] = r_geometry[2].X() - r_geometry[0].X();
    rv2[1] = r_geometry[2].Y() - r_geometry[0].Y();
    rv2[2] = r_geometry[2].Z() - r_geometry[0].Z();

    MathUtils<double>::CrossProduct(rAn,rv1,rv2);
    rAn *= 0.5;

    rCondition.SetValue(rNormalVariable, rAn);
}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::CalculateNormalShapeDerivative3D(
    ConditionType& rCondition
    )
{
    KRATOS_TRY

    const auto& rGeometry = rCondition.GetGeometry();

    // since this is a private static method (no external calls),
    // this method assumes matrix is properly sized always before this call.
    Matrix& r_normal_shape_derivatives = rCondition.GetValue(NORMAL_SHAPE_DERIVATIVE);

    KRATOS_DEBUG_ERROR_IF(r_normal_shape_derivatives.size1() != 9 || r_normal_shape_derivatives.size2() != 3)
        << "Condition (id: " << rCondition.Id() << ") does not have properly assigned \"NORMAL_SHAPE_DERIVATIVE\" "
        << "variable. Please assigned this variable with proper matrix size (required matrix size: (9, 3) != given matrix size ("
        << r_normal_shape_derivatives.size1() << ", " << r_normal_shape_derivatives.size2() << ")).";

    r_normal_shape_derivatives(0, 0) = 0;                                    // normal_x derivative w.r.t. node_0_x
    r_normal_shape_derivatives(1, 0) = rGeometry[1].Z() - rGeometry[2].Z();  // normal_x derivative w.r.t. node_0_y
    r_normal_shape_derivatives(2, 0) = -rGeometry[1].Y() + rGeometry[2].Y(); // normal_x derivative w.r.t. node_0_z
    r_normal_shape_derivatives(3, 0) = 0;                                    // normal_x derivative w.r.t. node_1_x
    r_normal_shape_derivatives(4, 0) = -rGeometry[0].Z() + rGeometry[2].Z(); // normal_x derivative w.r.t. node_1_y
    r_normal_shape_derivatives(5, 0) = rGeometry[0].Y() - rGeometry[2].Y();  // normal_x derivative w.r.t. node_1_z
    r_normal_shape_derivatives(6, 0) = 0;                                    // normal_x derivative w.r.t. node_2_x
    r_normal_shape_derivatives(7, 0) = rGeometry[0].Z() - rGeometry[1].Z();  // normal_x derivative w.r.t. node_2_y
    r_normal_shape_derivatives(8, 0) = -rGeometry[0].Y() + rGeometry[1].Y(); // normal_x derivative w.r.t. node_2_z

    r_normal_shape_derivatives(0, 1) = -rGeometry[1].Z() + rGeometry[2].Z(); // normal_y derivative w.r.t. node_0_x
    r_normal_shape_derivatives(1, 1) = 0;                                    // normal_y derivative w.r.t. node_0_y
    r_normal_shape_derivatives(2, 1) = rGeometry[1].X() - rGeometry[2].X();  // normal_y derivative w.r.t. node_0_z
    r_normal_shape_derivatives(3, 1) = rGeometry[0].Z() - rGeometry[2].Z();  // normal_y derivative w.r.t. node_1_x
    r_normal_shape_derivatives(4, 1) = 0;                                    // normal_y derivative w.r.t. node_1_y
    r_normal_shape_derivatives(5, 1) = -rGeometry[0].X() + rGeometry[2].X(); // normal_y derivative w.r.t. node_1_z
    r_normal_shape_derivatives(6, 1) = -rGeometry[0].Z() + rGeometry[1].Z(); // normal_y derivative w.r.t. node_2_x
    r_normal_shape_derivatives(7, 1) = 0;                                    // normal_y derivative w.r.t. node_2_y
    r_normal_shape_derivatives(8, 1) = rGeometry[0].X() - rGeometry[1].X();  // normal_y derivative w.r.t. node_2_z

    r_normal_shape_derivatives(0, 2) = rGeometry[1].Y() - rGeometry[2].Y();  // normal_z derivative w.r.t. node_0_x
    r_normal_shape_derivatives(1, 2) = -rGeometry[1].X() + rGeometry[2].X(); // normal_z derivative w.r.t. node_0_y
    r_normal_shape_derivatives(2, 2) = 0;                                    // normal_z derivative w.r.t. node_0_z
    r_normal_shape_derivatives(3, 2) = -rGeometry[0].Y() + rGeometry[2].Y(); // normal_z derivative w.r.t. node_1_x
    r_normal_shape_derivatives(4, 2) = rGeometry[0].X() - rGeometry[2].X();  // normal_z derivative w.r.t. node_1_y
    r_normal_shape_derivatives(5, 2) = 0;                                    // normal_z derivative w.r.t. node_1_z
    r_normal_shape_derivatives(6, 2) = rGeometry[0].Y() - rGeometry[1].Y();  // normal_z derivative w.r.t. node_2_x
    r_normal_shape_derivatives(7, 2) = -rGeometry[0].X() + rGeometry[1].X(); // normal_z derivative w.r.t. node_2_y
    r_normal_shape_derivatives(8, 2) = 0;                                    // normal_z derivative w.r.t. node_2_z

    noalias(r_normal_shape_derivatives) = r_normal_shape_derivatives * 0.5;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template <>
KRATOS_API(KRATOS_CORE) ModelPart::ConditionsContainerType& NormalCalculationUtils::GetContainer(
    ModelPart& rModelPart)
{
    return rModelPart.Conditions();
}

/***********************************************************************************/
/***********************************************************************************/

template <>
KRATOS_API(KRATOS_CORE) ModelPart::ElementsContainerType& NormalCalculationUtils::GetContainer(
    ModelPart& rModelPart)
{
    return rModelPart.Elements();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TContainerType, bool TIsHistorical>
void NormalCalculationUtils::CalculateNormalsUsingGenericAlgorithm(
    ModelPart& rModelPart,
    const bool ConsiderUnitNormal,
    const NormalVariableType& rNormalVariable
    )
{
    KRATOS_TRY

    // Initialize the normals
    InitializeNormals<TContainerType,TIsHistorical>(rModelPart, rNormalVariable);

    // Calculate normals in entities
    CalculateNormalsInContainer<TContainerType>(rModelPart, rNormalVariable);

    // Declare auxiliar coordinates
    Point::CoordinatesArrayType aux_coords;

    auto& r_entities_array = GetContainer<TContainerType>(rModelPart);
    const auto it_entity_begin = r_entities_array.begin();

    const std::function<array_1d<double, 3>(const GeometryType&, const GeometryType::CoordinatesArrayType&, const double)> retrieve_normal_unit_normal =
    [](const GeometryType& rGeometry, const GeometryType::CoordinatesArrayType& rLocalCoordinates, const double Coefficient) -> array_1d<double, 3> {return rGeometry.UnitNormal(rLocalCoordinates);};
    const std::function<array_1d<double, 3>(const GeometryType&, const GeometryType::CoordinatesArrayType&, const double)> retrieve_normal_area_normal =
    [](const GeometryType& rGeometry, const GeometryType::CoordinatesArrayType& rLocalCoordinates, const double Coefficient) -> array_1d<double, 3> {return Coefficient * rGeometry.Normal(rLocalCoordinates);};

    const auto* p_retrieve_normal = ConsiderUnitNormal ? &retrieve_normal_unit_normal : &retrieve_normal_area_normal;

    // TODO: Use TLS in ParallelUtilities
    #pragma omp parallel for firstprivate(aux_coords)
    for (int i = 0; i < static_cast<int>(r_entities_array.size()); ++i) {
        auto it_entity = it_entity_begin + i;
        auto& r_geometry = it_entity->GetGeometry();

        // Avoid not "flat" elements
        if (r_geometry.WorkingSpaceDimension() != r_geometry.LocalSpaceDimension() + 1) {
            continue;
        }

        // Iterate over nodes
        const double coefficient = 1.0 / static_cast<double>(r_geometry.PointsNumber());
        for (auto& r_node : r_geometry) {
            r_geometry.PointLocalCoordinates(aux_coords, r_node.Coordinates());
            r_node.SetLock();
            noalias(GetNormalValue<TIsHistorical>(r_node, rNormalVariable)) += (*p_retrieve_normal)(r_geometry, aux_coords, coefficient);
            r_node.UnSetLock();
        }
    }

    // For MPI: correct values on partition boundaries
    if (TIsHistorical) {
        rModelPart.GetCommunicator().AssembleCurrentData(rNormalVariable);
    } else {
        rModelPart.GetCommunicator().AssembleNonHistoricalData(rNormalVariable);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
KRATOS_API(KRATOS_CORE) array_1d<double,3>& NormalCalculationUtils::GetNormalValue<true>(
    NodeType& rNode,
    const NormalVariableType& rNormalVariable
    )
{
    return rNode.FastGetSolutionStepValue(rNormalVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
KRATOS_API(KRATOS_CORE) array_1d<double,3>& NormalCalculationUtils::GetNormalValue<false>(
    NodeType& rNode,
    const NormalVariableType& rNormalVariable
    )
{
    return rNode.GetValue(rNormalVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::SetNormalValue<true>(
    NodeType& rNode,
    const NormalVariableType& rNormalVariable,
    const array_1d<double,3>& rNormalValue
    )
{
    noalias(rNode.FastGetSolutionStepValue(rNormalVariable)) = rNormalValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::SetNormalValue<false>(
    NodeType& rNode,
    const NormalVariableType& rNormalVariable,
    const array_1d<double,3>& rNormalValue
    )
{
    rNode.SetValue(rNormalVariable, rNormalValue);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_API(KRATOS_CORE) bool NormalCalculationUtils::CheckUseSimplex(
    const ModelPart& rModelPart,
    const bool EnforceGenericGeometryAlgorithm
    )
{
    // Getting dimension
    const auto& r_process_info = rModelPart.GetProcessInfo();
    const SizeType dimension = r_process_info.GetValue(DOMAIN_SIZE);

    // Check that the conditions are simplex
    // Note that in here we assume unique condition type in the model part
    bool use_simplex = false;
    auto& r_conditions_array = rModelPart.Conditions();
    if (r_conditions_array.size() != 0) {
        // Get geometry type from first condition
        const auto it_cond_begin = r_conditions_array.begin();
        const auto geometry_type = it_cond_begin->GetGeometry().GetGeometryType();

        // Checking if we can compute with simplex
        const bool is_condition_simplex = dimension == 2
            ? geometry_type == GeometryData::KratosGeometryType::Kratos_Line2D2
            : geometry_type == GeometryData::KratosGeometryType::Kratos_Triangle3D3;
        use_simplex = EnforceGenericGeometryAlgorithm ? false : is_condition_simplex;
    }

    // Synchronize obtained results
    use_simplex = rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(use_simplex);

    return use_simplex;
}

/***********************************************************************************/
/***********************************************************************************/

template<bool TIsHistorical>
void NormalCalculationUtils::AuxiliaryCalculateOnSimplex(
    ConditionsArrayType& rConditions,
    const std::size_t Dimension,
    const NormalVariableType& rNormalVariable
    )
{
    KRATOS_TRY

    // Calculating the normals and storing on the conditions
    if (Dimension == 2) {
        block_for_each(rConditions, array_1d<double,3>(), [&rNormalVariable](Condition& rCondition, array_1d<double,3>& rAn){
            if (rCondition.GetGeometry().PointsNumber() == 2) {
                CalculateNormal2D(rCondition, rAn, rNormalVariable);
            }
        });
    } else if (Dimension == 3) {
        array_1d<double,3> An;
        array_1d<double,3> v1;
        array_1d<double,3> v2;
        auto aux_tls = std::make_tuple(An, v1, v2);
        using TLSType = std::tuple<array_1d<double,3>,array_1d<double,3>,array_1d<double,3>>;
        block_for_each(rConditions, aux_tls, [&rNormalVariable](Condition& rCondition, TLSType& rTLS){
            if (rCondition.GetGeometry().PointsNumber() == 3) {
                auto& r_An = std::get<0>(rTLS);
                auto& r_v1 = std::get<1>(rTLS);
                auto& r_v2 = std::get<2>(rTLS);
                CalculateNormal3D(rCondition, r_An, r_v1, r_v2, rNormalVariable);
            }
        });
    }

    const auto condition_simplex_type = Dimension == 2
            ? GeometryData::KratosGeometryType::Kratos_Line2D2
            : GeometryData::KratosGeometryType::Kratos_Triangle3D3;

    // Adding the normals to the nodes
    block_for_each(rConditions, [&rNormalVariable, &condition_simplex_type, this](Condition& rCondition) {
        auto& r_geometry = rCondition.GetGeometry();
        if (rCondition.GetGeometry().GetGeometryType() == condition_simplex_type) {
            const double coeff = 1.0 / r_geometry.size();
            const auto& r_normal = rCondition.GetValue(rNormalVariable);
            for (unsigned int i = 0; i < r_geometry.size(); ++i) {
                auto& r_node = r_geometry[i];
                r_node.SetLock();
                noalias(GetNormalValue<TIsHistorical>(r_node, rNormalVariable)) += coeff * r_normal;
                r_node.UnSetLock();
            }
        }
    });

    KRATOS_CATCH("AuxiliaryCalculateOnSimplex");
}

// template instantiations

template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::CalculateNormalsUsingGenericAlgorithm<ModelPart::ConditionsContainerType, true>(ModelPart&, const bool, const NormalVariableType&);
template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::CalculateNormalsUsingGenericAlgorithm<ModelPart::ConditionsContainerType, false>(ModelPart&, const bool, const NormalVariableType&);
template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::CalculateNormalsUsingGenericAlgorithm<ModelPart::ElementsContainerType, true>(ModelPart&, const bool, const NormalVariableType&);
template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::CalculateNormalsUsingGenericAlgorithm<ModelPart::ElementsContainerType, false>(ModelPart&, const bool, const NormalVariableType&);
template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::InitializeNormals<ModelPart::ConditionsContainerType, true>(ModelPart&, const NormalVariableType&);
template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::InitializeNormals<ModelPart::ConditionsContainerType, false>(ModelPart&, const NormalVariableType&);
template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::InitializeNormals<ModelPart::ElementsContainerType, true>(ModelPart&, const NormalVariableType&);
template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::InitializeNormals<ModelPart::ElementsContainerType, false>(ModelPart&, const NormalVariableType&);
template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::CalculateNormalsInContainer<ModelPart::ConditionsContainerType>(ModelPart&, const Variable<array_1d<double,3>>&);
template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::CalculateNormalsInContainer<ModelPart::ElementsContainerType>(ModelPart&, const Variable<array_1d<double,3>>&);
template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::CalculateUnitNormals<ModelPart::ConditionsContainerType, true>(ModelPart&, const bool, const NormalVariableType&);
template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::CalculateUnitNormals<ModelPart::ConditionsContainerType, false>(ModelPart&, const bool, const NormalVariableType&);
template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::CalculateUnitNormals<ModelPart::ElementsContainerType, true>(ModelPart&, const bool, const NormalVariableType&);
template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::CalculateUnitNormals<ModelPart::ElementsContainerType, false>(ModelPart&, const bool, const NormalVariableType&);
template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::AuxiliaryCalculateOnSimplex<true>(ConditionsArrayType&, const std::size_t, const NormalVariableType&);
template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::AuxiliaryCalculateOnSimplex<false>(ConditionsArrayType&, const std::size_t, const NormalVariableType&);

} // namespace Kratos
