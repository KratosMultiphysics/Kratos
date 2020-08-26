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
//
//

// System includes

// External includes

// Project includes
#include "utilities/normal_calculation_utils.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
template <class TContainerType>
void NormalCalculationUtils::CalculateNormalsInContainer(ModelPart& rModelPart)
{
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
        it_entity->SetValue(NORMAL, r_geometry.UnitNormal(aux_coords));
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TContainerType>
void NormalCalculationUtils::InitializeNormals(ModelPart& rModelPart)
{
    // Resetting the normals
    const array_1d<double, 3> zero = ZeroVector(3);

    auto& r_entity_array = GetContainer<TContainerType>(rModelPart);

    if (rModelPart.GetCommunicator().GetDataCommunicator().IsDistributed()) {
        // If Parallel make sure normals are reset in all partitions
        VariableUtils().SetFlag(VISITED, false, rModelPart.Nodes());

        BlockPartition<TContainerType>(r_entity_array).for_each([](typename TContainerType::value_type& rEntity) {
            for (auto& r_node : rEntity.GetGeometry()) {
                r_node.SetLock();
                r_node.Set(VISITED, true);
                r_node.UnSetLock();
            }
        });

        rModelPart.GetCommunicator().SynchronizeOrNodalFlags(VISITED);
        VariableUtils().SetVariable(NORMAL, zero, rModelPart.Nodes(), VISITED);
    } else {
        // In serial iteratre normally over the condition nodes
        BlockPartition<TContainerType>(r_entity_array)
            .for_each([zero](typename TContainerType::value_type& rEntity) {
                for (auto& r_node : rEntity.GetGeometry()) {
                    r_node.SetLock();
                    noalias(r_node.FastGetSolutionStepValue(NORMAL)) = zero;
                    r_node.UnSetLock();
                }
            });
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void NormalCalculationUtils::CalculateNormals<Condition>(
    ModelPart& rModelPart,
    const bool EnforceGenericGeometryAlgorithm)
{
    // Getting process info
    const auto& r_process_info = rModelPart.GetProcessInfo();

    // Getting dimension
    const SizeType dimension = r_process_info.GetValue(DOMAIN_SIZE);

    // Sum all the nodes normals
    auto& r_conditions_array = rModelPart.Conditions();
    const auto it_cond_begin = r_conditions_array.begin();

    // Checking if we can compute with simplex
    const GeometryData::KratosGeometryType geometry_type =
        it_cond_begin->GetGeometry().GetGeometryType();
    const bool use_simplex =
        EnforceGenericGeometryAlgorithm
            ? false
            : dimension == 2
                  ? geometry_type == GeometryData::KratosGeometryType::Kratos_Line2D2
                  : geometry_type == GeometryData::KratosGeometryType::Kratos_Triangle3D3;

    if (use_simplex) {
        CalculateOnSimplex(rModelPart, dimension);
    } else {
        CalculateNormalsUsingGenericAlgorithm<ModelPart::ConditionsContainerType>(rModelPart);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void NormalCalculationUtils::CalculateNormals<Element>(
    ModelPart& rModelPart,
    const bool
    )
{
    CalculateNormalsUsingGenericAlgorithm<ModelPart::ElementsContainerType>(rModelPart);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntityType>
void NormalCalculationUtils::CalculateUnitNormals(
    ModelPart& rModelPart,
    const bool EnforceGenericGeometryAlgorithm
    )
{
    // Compute area normals
    CalculateNormals<TEntityType>(rModelPart, EnforceGenericGeometryAlgorithm);

    // Compute unit normals
    ComputeUnitNormalsFromAreaNormals(rModelPart);
}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::CalculateOnSimplex(
    ConditionsArrayType& rConditions,
    const std::size_t Dimension
    )
{
    KRATOS_TRY

    // Calculating the normals and storing on the conditions
    array_1d<double, 3> An;
    if (Dimension == 2) {
        // TODO: Use TLS OpenMP Parallel utils
#pragma omp parallel for firstprivate(An)
        for (int i_cond = 0; i_cond < static_cast<int>(rConditions.size()); ++i_cond) {
            auto it = rConditions.begin() + i_cond;
            if (it->GetGeometry().PointsNumber() == 2)
                CalculateNormal2D(*it, An);
        }
    } else if (Dimension == 3) {
        // TODO: Use TLS OpenMP Parallel utils
        array_1d<double, 3> v1, v2;
#pragma omp parallel for firstprivate(An, v1, v2)
        for (int i_cond = 0; i_cond < static_cast<int>(rConditions.size()); ++i_cond) {
            auto it = rConditions.begin() + i_cond;
            if (it->GetGeometry().PointsNumber() == 3)
                CalculateNormal3D(*it, An, v1, v2);
        }
    }

    // Adding the normals to the nodes
    BlockPartition<ConditionsArrayType>(rConditions).for_each([](Condition& rCondition) {
        auto& r_geometry = rCondition.GetGeometry();
        const double coeff = 1.0 / r_geometry.size();
        const auto& r_normal = rCondition.GetValue(NORMAL);
        for (unsigned int i = 0; i < r_geometry.size(); ++i) {
            auto& r_node = r_geometry[i];
            r_node.SetLock();
            noalias(r_node.FastGetSolutionStepValue(NORMAL)) += coeff * r_normal;
            r_node.UnSetLock();
        }
    });

    KRATOS_CATCH("")

}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::CalculateOnSimplex(
    ModelPart& rModelPart)
{
    KRATOS_TRY

    const auto& r_process_info = rModelPart.GetProcessInfo();

    KRATOS_ERROR_IF(!r_process_info.Has(DOMAIN_SIZE))
        << "DOMAIN_SIZE not found in process info of " << rModelPart.Name() << ".";
    CalculateOnSimplex(rModelPart, r_process_info[DOMAIN_SIZE]);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::CalculateOnSimplex(
    ModelPart& rModelPart,
    const std::size_t Dimension
    )
{
    // Initialize the normals
    InitializeNormals<ModelPart::ConditionsContainerType>(rModelPart);

    // Calling CalculateOnSimplex for conditions
    this->CalculateOnSimplex(rModelPart.Conditions(), Dimension);

    // Synchronize the normal
    rModelPart.GetCommunicator().AssembleCurrentData(NORMAL);
}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::SwapNormals(ModelPart& rModelPart)
{
    KRATOS_TRY

    for(auto& r_cond : rModelPart.Conditions()) {
        GeometryType& r_geometry = r_cond.GetGeometry();
        Node<3>::Pointer paux = r_geometry(0);
        r_geometry(0) = r_geometry(1);
        r_geometry(1) = paux;
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::ComputeUnitNormalsFromAreaNormals(ModelPart& rModelPart)
{
    KRATOS_TRY

    // We iterate over nodes
    auto& r_nodes_array = rModelPart.GetCommunicator().LocalMesh().Nodes();

    BlockPartition<ModelPart::NodesContainerType>(r_nodes_array).for_each([](NodeType& rNode) {
        auto& r_normal = rNode.FastGetSolutionStepValue(NORMAL);
        const double norm_normal = norm_2(r_normal);

        KRATOS_ERROR_IF(rNode.Is(INTERFACE) && norm_normal <= std::numeric_limits<double>::epsilon())
            << "ERROR:: ZERO NORM NORMAL IN NODE: " << rNode.Id() << std::endl;
        r_normal /= norm_normal;
    });

    // For MPI: correct values on partition boundaries
    rModelPart.GetCommunicator().SynchronizeVariable(NORMAL);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::CalculateNormal2D(
    Condition& rCondition,
    array_1d<double,3>& rAn
    )
{
    const GeometryType& r_geometry = rCondition.GetGeometry();

    rAn[0] =    r_geometry[1].Y() - r_geometry[0].Y();
    rAn[1] = - (r_geometry[1].X() - r_geometry[0].X());
    rAn[2] =    0.00;

    rCondition.SetValue(NORMAL, rAn);
}

/***********************************************************************************/
/***********************************************************************************/

void NormalCalculationUtils::CalculateNormal3D(
    Condition& rCondition,
    array_1d<double,3>& rAn,
    array_1d<double,3>& rv1,
    array_1d<double,3>& rv2
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

    rCondition.SetValue(NORMAL, rAn);
}

template <>
ModelPart::ConditionsContainerType& NormalCalculationUtils::GetContainer(
    ModelPart& rModelPart)
{
    return rModelPart.Conditions();
}

template <>
ModelPart::ElementsContainerType& NormalCalculationUtils::GetContainer(
    ModelPart& rModelPart)
{
    return rModelPart.Elements();
}


template<class TContainerType>
void NormalCalculationUtils::CalculateNormalsUsingGenericAlgorithm(
    ModelPart& rModelPart)
{
    KRATOS_TRY

    // Initialize the normals
    InitializeNormals<TContainerType>(rModelPart);

    // Calculate normals in entities
    CalculateNormalsInContainer<TContainerType>(rModelPart);

    // Declare auxiliar coordinates
    Point::CoordinatesArrayType aux_coords;

    auto& r_entities_array = GetContainer<TContainerType>(rModelPart);
    const auto it_entity_begin = r_entities_array.begin();

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
            const auto normal = r_geometry.Normal(aux_coords);
            r_node.SetLock();
            noalias(r_node.FastGetSolutionStepValue(NORMAL)) += normal * coefficient;
            r_node.UnSetLock();
        }
    }

    // For MPI: correct values on partition boundaries
    rModelPart.GetCommunicator().AssembleCurrentData(NORMAL);

    KRATOS_CATCH("");
}


// template instantiations

template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::CalculateNormalsUsingGenericAlgorithm<ModelPart::ConditionsContainerType>(ModelPart&);
template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::CalculateNormalsUsingGenericAlgorithm<ModelPart::ElementsContainerType>(ModelPart&);
template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::InitializeNormals<ModelPart::ConditionsContainerType>(ModelPart&);
template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::InitializeNormals<ModelPart::ElementsContainerType>(ModelPart&);
template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::CalculateNormalsInContainer<ModelPart::ConditionsContainerType>(ModelPart&);
template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::CalculateNormalsInContainer<ModelPart::ElementsContainerType>(ModelPart&);
template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::CalculateUnitNormals<Condition>(ModelPart&, const bool);
template KRATOS_API(KRATOS_CORE) void NormalCalculationUtils::CalculateUnitNormals<Element>(ModelPart&, const bool);

} // namespace Kratos
