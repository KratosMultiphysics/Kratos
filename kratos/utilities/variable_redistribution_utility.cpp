//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//	                 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Jordi Cotela
//

#include "includes/fsi_variables.h"
#include "utilities/atomic_utilities.h"
#include "utilities/parallel_utilities.h"
#include "variable_redistribution_utility.h"

namespace Kratos
{

// Historical API
// ConvertDistributedValuesToPoint for double variables
void VariableRedistributionUtility::ConvertDistributedValuesToPoint(
    ModelPart& rModelPart,
    const Variable<double>& rDistributedVariable,
    const Variable<double>& rPointVariable)
{
    KRATOS_WARNING("VariableRedistributionUtility") << "This ConvertDistributedValuesToPoint() signature is deprecated. Use the one with Elements/Conditions container." << std::endl;
    CallSpecializedConvertDistributedValuesToPoint<true>(rModelPart,rModelPart.Conditions(),rDistributedVariable,rPointVariable);
}

void VariableRedistributionUtility::ConvertDistributedValuesToPoint(
    ModelPart& rModelPart,
    ModelPart::ConditionsContainerType& rConditions,
    const Variable<double>& rDistributedVariable,
    const Variable<double>& rPointVariable)
{
    CallSpecializedConvertDistributedValuesToPoint<true>(rModelPart,rConditions,rDistributedVariable,rPointVariable);
}

void VariableRedistributionUtility::ConvertDistributedValuesToPoint(
    ModelPart& rModelPart,
    ModelPart::ElementsContainerType& rElements,
    const Variable<double>& rDistributedVariable,
    const Variable<double>& rPointVariable)
{
    CallSpecializedConvertDistributedValuesToPoint<true>(rModelPart,rElements,rDistributedVariable,rPointVariable);
}

// ConvertDistributedValuesToPoint for array variables
void VariableRedistributionUtility::ConvertDistributedValuesToPoint(
    ModelPart& rModelPart,
    const Variable< array_1d<double,3> >& rDistributedVariable,
    const Variable< array_1d<double,3> >& rPointVariable)
{
    KRATOS_WARNING("VariableRedistributionUtility") << "This ConvertDistributedValuesToPoint() signature is deprecated. Use the one with Elements/Conditions container." << std::endl;
    CallSpecializedConvertDistributedValuesToPoint<true>(rModelPart,rModelPart.Conditions(),rDistributedVariable,rPointVariable);
}

void VariableRedistributionUtility::ConvertDistributedValuesToPoint(
    ModelPart& rModelPart,
    ModelPart::ConditionsContainerType& rConditions,
    const Variable< array_1d<double,3> >& rDistributedVariable,
    const Variable< array_1d<double,3> >& rPointVariable)
{
    CallSpecializedConvertDistributedValuesToPoint<true>(rModelPart,rConditions,rDistributedVariable,rPointVariable);
}

void VariableRedistributionUtility::ConvertDistributedValuesToPoint(
    ModelPart& rModelPart,
    ModelPart::ElementsContainerType& rElements,
    const Variable< array_1d<double,3> >& rDistributedVariable,
    const Variable< array_1d<double,3> >& rPointVariable)
{
    CallSpecializedConvertDistributedValuesToPoint<true>(rModelPart,rElements,rDistributedVariable,rPointVariable);
}

// DistributePointValues for double variables
void VariableRedistributionUtility::DistributePointValues(
    ModelPart& rModelPart,
    const Variable<double>& rPointVariable,
    const Variable<double>& rDistributedVariable,
    double Tolerance,
    unsigned int MaximumIterations)
{
    KRATOS_WARNING("VariableRedistributionUtility") << "This DistributePointValues() signature is deprecated. Use the one with Elements/Conditions container." << std::endl;
    CallSpecializedDistributePointValues<true>(rModelPart, rModelPart.Conditions(), rPointVariable, rDistributedVariable, Tolerance, MaximumIterations);
}

void VariableRedistributionUtility::DistributePointValues(
    ModelPart& rModelPart,
    ModelPart::ConditionsContainerType& rConditions,
    const Variable<double>& rPointVariable,
    const Variable<double>& rDistributedVariable,
    double Tolerance,
    unsigned int MaximumIterations)
{
    CallSpecializedDistributePointValues<true>(rModelPart, rConditions, rPointVariable, rDistributedVariable, Tolerance, MaximumIterations);
}

void VariableRedistributionUtility::DistributePointValues(
    ModelPart& rModelPart,
    ModelPart::ElementsContainerType& rElements,
    const Variable<double>& rPointVariable,
    const Variable<double>& rDistributedVariable,
    double Tolerance,
    unsigned int MaximumIterations)
{
    CallSpecializedDistributePointValues<true>(rModelPart, rElements, rPointVariable, rDistributedVariable, Tolerance, MaximumIterations);
}

// DistributePointValues for array variables
void VariableRedistributionUtility::DistributePointValues(
    ModelPart& rModelPart,
    const Variable< array_1d<double,3> >& rPointVariable,
    const Variable< array_1d<double,3> >& rDistributedVariable,
    double Tolerance,
    unsigned int MaximumIterations)
{
    KRATOS_WARNING("VariableRedistributionUtility") << "This DistributePointValues() signature is deprecated. Use the one with Elements/Conditions container." << std::endl;
    CallSpecializedDistributePointValues<true>(rModelPart, rModelPart.Conditions(), rPointVariable, rDistributedVariable, Tolerance, MaximumIterations);
}

void VariableRedistributionUtility::DistributePointValues(
    ModelPart& rModelPart,
    ModelPart::ConditionsContainerType& rConditions,
    const Variable< array_1d<double,3> >& rPointVariable,
    const Variable< array_1d<double,3> >& rDistributedVariable,
    double Tolerance,
    unsigned int MaximumIterations)
{
    CallSpecializedDistributePointValues<true>(rModelPart, rConditions, rPointVariable, rDistributedVariable, Tolerance, MaximumIterations);
}

void VariableRedistributionUtility::DistributePointValues(
    ModelPart& rModelPart,
    ModelPart::ElementsContainerType& rElements,
    const Variable< array_1d<double,3> >& rPointVariable,
    const Variable< array_1d<double,3> >& rDistributedVariable,
    double Tolerance,
    unsigned int MaximumIterations)
{
    CallSpecializedDistributePointValues<true>(rModelPart, rElements, rPointVariable, rDistributedVariable, Tolerance, MaximumIterations);
}

// Non-historical API
// ConvertDistributedValuesToPointNonHistorical for double variables
void VariableRedistributionUtility::ConvertDistributedValuesToPointNonHistorical(
    ModelPart& rModelPart,
    ModelPart::ConditionsContainerType& rConditions,
    const Variable<double>& rDistributedVariable,
    const Variable<double>& rPointVariable)
{
    CallSpecializedConvertDistributedValuesToPoint<false>(rModelPart,rConditions,rDistributedVariable,rPointVariable);
}

void VariableRedistributionUtility::ConvertDistributedValuesToPointNonHistorical(
    ModelPart& rModelPart,
    ModelPart::ElementsContainerType& rElements,
    const Variable<double>& rDistributedVariable,
    const Variable<double>& rPointVariable)
{
    CallSpecializedConvertDistributedValuesToPoint<false>(rModelPart,rElements,rDistributedVariable,rPointVariable);
}

// ConvertDistributedValuesToPointNonHistorical for array variables
void VariableRedistributionUtility::ConvertDistributedValuesToPointNonHistorical(
    ModelPart& rModelPart,
    ModelPart::ConditionsContainerType& rConditions,
    const Variable< array_1d<double,3> >& rDistributedVariable,
    const Variable< array_1d<double,3> >& rPointVariable)
{
    CallSpecializedConvertDistributedValuesToPoint<false>(rModelPart,rConditions,rDistributedVariable,rPointVariable);
}

void VariableRedistributionUtility::ConvertDistributedValuesToPointNonHistorical(
    ModelPart& rModelPart,
    ModelPart::ElementsContainerType& rElements,
    const Variable< array_1d<double,3> >& rDistributedVariable,
    const Variable< array_1d<double,3> >& rPointVariable)
{
    CallSpecializedConvertDistributedValuesToPoint<false>(rModelPart,rElements,rDistributedVariable,rPointVariable);
}

// DistributePointValuesNonHistorical for double variables
void VariableRedistributionUtility::DistributePointValuesNonHistorical(
    ModelPart& rModelPart,
    ModelPart::ConditionsContainerType& rConditions,
    const Variable<double>& rPointVariable,
    const Variable<double>& rDistributedVariable,
    double Tolerance,
    unsigned int MaximumIterations)
{
    CallSpecializedDistributePointValues<false>(rModelPart, rConditions, rPointVariable, rDistributedVariable, Tolerance, MaximumIterations);
}

void VariableRedistributionUtility::DistributePointValuesNonHistorical(
    ModelPart& rModelPart,
    ModelPart::ElementsContainerType& rElements,
    const Variable<double>& rPointVariable,
    const Variable<double>& rDistributedVariable,
    double Tolerance,
    unsigned int MaximumIterations)
{
    CallSpecializedDistributePointValues<false>(rModelPart, rElements, rPointVariable, rDistributedVariable, Tolerance, MaximumIterations);
}

// DistributePointValuesNonHistorical for array variables
void VariableRedistributionUtility::DistributePointValuesNonHistorical(
    ModelPart& rModelPart,
    ModelPart::ConditionsContainerType& rConditions,
    const Variable< array_1d<double,3> >& rPointVariable,
    const Variable< array_1d<double,3> >& rDistributedVariable,
    double Tolerance,
    unsigned int MaximumIterations)
{
    CallSpecializedDistributePointValues<false>(rModelPart, rConditions, rPointVariable, rDistributedVariable, Tolerance, MaximumIterations);
}

void VariableRedistributionUtility::DistributePointValuesNonHistorical(
    ModelPart& rModelPart,
    ModelPart::ElementsContainerType& rElements,
    const Variable< array_1d<double,3> >& rPointVariable,
    const Variable< array_1d<double,3> >& rDistributedVariable,
    double Tolerance,
    unsigned int MaximumIterations)
{
    CallSpecializedDistributePointValues<false>(rModelPart, rElements, rPointVariable, rDistributedVariable, Tolerance, MaximumIterations);
}

// Private methods
template< const bool TIsHistorical, class TContainerType, class TValueType >
void VariableRedistributionUtility::CallSpecializedConvertDistributedValuesToPoint(
    ModelPart& rModelPart,
    TContainerType& rEntitiesContainer,
    const Variable<TValueType>& rDistributedVariable,
    const Variable<TValueType>& rPointVariable)
{
    // Check if there is any entity in the current partition
    const int n_loc_entities = NumberOfLocalEntities(rModelPart, rEntitiesContainer);
    const int n_tot_entities = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(n_loc_entities);

    // If there is conditions, this function dispatches the call to the correct specialization
    if (n_tot_entities != 0){
        if (n_loc_entities != 0){
            Geometry< Node<3> >& rReferenceGeometry = rEntitiesContainer.begin()->GetGeometry();
            const GeometryData::KratosGeometryFamily GeometryFamily = rReferenceGeometry.GetGeometryFamily();
            const unsigned int PointsNumber = rReferenceGeometry.PointsNumber();

            if (GeometryFamily == GeometryData::Kratos_Linear && PointsNumber == 2){
                VariableRedistributionUtility::SpecializedConvertDistributedValuesToPoint<TIsHistorical,TContainerType,GeometryData::Kratos_Linear,2,TValueType>(
                    rModelPart,
                    rEntitiesContainer,
                    rDistributedVariable,
                    rPointVariable);
            } else if (GeometryFamily == GeometryData::Kratos_Triangle && PointsNumber == 3) {
                VariableRedistributionUtility::SpecializedConvertDistributedValuesToPoint<TIsHistorical,TContainerType,GeometryData::Kratos_Triangle,3,TValueType>(
                    rModelPart,
                    rEntitiesContainer,
                    rDistributedVariable,
                    rPointVariable);
            } else if (GeometryFamily == GeometryData::Kratos_Quadrilateral && PointsNumber == 4) {
                VariableRedistributionUtility::SpecializedConvertDistributedValuesToPoint<TIsHistorical,TContainerType,GeometryData::Kratos_Quadrilateral,4,TValueType>(
                    rModelPart,
                    rEntitiesContainer,
                    rDistributedVariable,
                    rPointVariable);
            } else {
                KRATOS_ERROR << "Unsupported geometry type with " << PointsNumber << " points." << std::endl;
            }
        } else {
            VariableRedistributionUtility::DummySpecializedConvertDistributedValuesToPoint<TIsHistorical, TValueType>(
                rModelPart,
                rDistributedVariable,
                rPointVariable);
        }
    }
}

template< const bool TIsHistorical, class TContainerType, class TValueType >
void VariableRedistributionUtility::CallSpecializedDistributePointValues(
    ModelPart& rModelPart,
    TContainerType& rEntitiesContainer,
    const Variable<TValueType>& rPointVariable,
    const Variable<TValueType>& rDistributedVariable,
    double Tolerance,
    unsigned int MaximumIterations)
{
    // Check if there is any entity in the current partition
    const int n_loc_entities = NumberOfLocalEntities(rModelPart, rEntitiesContainer);
    const int n_tot_entities = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(n_loc_entities);

    // If there is conditions, this function dispatches the call to the correct specialization
    if (n_tot_entities != 0){
        if (n_loc_entities != 0){
            Geometry< Node<3> >& rReferenceGeometry = rEntitiesContainer.begin()->GetGeometry();
            const GeometryData::KratosGeometryFamily GeometryFamily = rReferenceGeometry.GetGeometryFamily();
            const unsigned int PointsNumber = rReferenceGeometry.PointsNumber();

            if (GeometryFamily == GeometryData::Kratos_Linear && PointsNumber == 2){
                VariableRedistributionUtility::SpecializedDistributePointValues<TIsHistorical, TContainerType,GeometryData::Kratos_Linear,2,TValueType>(
                    rModelPart,
                    rEntitiesContainer,
                    rPointVariable,
                    rDistributedVariable,
                    Tolerance,
                    MaximumIterations);
            } else if (GeometryFamily == GeometryData::Kratos_Triangle && PointsNumber == 3){
                VariableRedistributionUtility::SpecializedDistributePointValues<TIsHistorical, TContainerType,GeometryData::Kratos_Triangle,3,TValueType>(
                    rModelPart,
                    rEntitiesContainer,
                    rPointVariable,
                    rDistributedVariable,
                    Tolerance,
                    MaximumIterations);
            } else if (GeometryFamily == GeometryData::Kratos_Quadrilateral && PointsNumber == 4){
                VariableRedistributionUtility::SpecializedDistributePointValues<TIsHistorical, TContainerType,GeometryData::Kratos_Quadrilateral,4,TValueType>(
                    rModelPart,
                    rEntitiesContainer,
                    rPointVariable,
                    rDistributedVariable,
                    Tolerance,
                    MaximumIterations);
            } else {
                KRATOS_ERROR << "Unsupported geometry type with " << PointsNumber << " points." << std::endl;
            }
        } else {
            VariableRedistributionUtility::DummySpecializedDistributePointValues<TIsHistorical, TValueType>(
                rModelPart,
                rDistributedVariable,
                Tolerance,
                MaximumIterations);
        }
    }

    // Synchronize the obtained resultas
    rModelPart.GetCommunicator().SynchronizeVariable(rDistributedVariable);
}

///////////////////////////////////////////////////////////////////////////////

template< const bool TIsHistorical, class TContainerType, GeometryData::KratosGeometryFamily TFamily, unsigned int TPointNumber, class TValueType >
void VariableRedistributionUtility::SpecializedConvertDistributedValuesToPoint(
    ModelPart& rModelPart,
    TContainerType& rEntitiesContainer,
    const Variable< TValueType >& rDistributedVariable,
    const Variable< TValueType >& rPointVariable)
{
    // Get the corresponding mass matrix
    BoundedMatrix<double,TPointNumber,TPointNumber> MassMatrix;
    ConsistentMassMatrix<TFamily,TPointNumber>(MassMatrix);

    // Initialize result to zero
    const TValueType zero = rPointVariable.Zero();
    block_for_each(rModelPart.Nodes(), zero, [&](NodeType& rNode, const TValueType& rZero){
        AuxiliarySet<TIsHistorical>(rPointVariable, rZero, rNode);
    });

    // Make sure that the distributed values are equal between processors
    if (TIsHistorical) {
        rModelPart.GetCommunicator().SynchronizeVariable(rDistributedVariable);
    } else {
        rModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rDistributedVariable);
    }

    TValueType value_j;
    block_for_each(rEntitiesContainer, value_j, [&](typename TContainerType::value_type& rEntity, TValueType& rValueJ){
        auto& r_geometry = rEntity.GetGeometry();
        const double size = r_geometry.DomainSize();
        for (unsigned int j = 0; j < TPointNumber; ++j) {
            rValueJ = rDistributedVariable.Zero();
            for (unsigned int k = 0; k < TPointNumber; ++k) {
                rValueJ += size * MassMatrix(j,k) * r_geometry[k].FastGetSolutionStepValue(rDistributedVariable);
            }
            ThreadsafeAdd(AuxiliaryGet<TIsHistorical>(rPointVariable, r_geometry[j]), rValueJ);
        }
    });

    // Add the contributions between processors
    if (TIsHistorical) {
        rModelPart.GetCommunicator().AssembleCurrentData(rPointVariable);
    } else {
        rModelPart.GetCommunicator().AssembleNonHistoricalData(rPointVariable);
    }
}

template< const bool TIsHistorical, class TValueType >
void VariableRedistributionUtility::DummySpecializedConvertDistributedValuesToPoint(
    ModelPart& rModelPart,
    const Variable< TValueType >& rDistributedVariable,
    const Variable< TValueType >& rPointVariable)
{
    // Make sure that the distributed values are equal between processors
    if (TIsHistorical) {
        rModelPart.GetCommunicator().SynchronizeVariable(rDistributedVariable);
    } else {
        rModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rDistributedVariable);
    }

    // Add the contributions between processors
    if (TIsHistorical) {
        rModelPart.GetCommunicator().AssembleCurrentData(rPointVariable);
    } else {
        rModelPart.GetCommunicator().AssembleNonHistoricalData(rPointVariable);
    }
}

///////////////////////////////////////////////////////////////////////////////

template< const bool TIsHistorical, class TContainerType, GeometryData::KratosGeometryFamily TFamily, unsigned int TPointNumber, class TValueType >
void VariableRedistributionUtility::SpecializedDistributePointValues(
    ModelPart& rModelPart,
    TContainerType& rEntitiesContainer,
    const Variable< TValueType >& rPointVariable,
    const Variable< TValueType >& rDistributedVariable,
    double Tolerance,
    unsigned int MaximumIterations)
{
    // Get the corresponding geometry mass matrix
    BoundedMatrix< double, TPointNumber, TPointNumber > mass_matrix;
    ConsistentMassMatrix<TFamily,TPointNumber>(mass_matrix);

    // Initial guess (NodalValue / NodalSize)
    ComputeNodalSizes(rModelPart, rEntitiesContainer);
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        const double nodal_maux = rNode.GetValue(NODAL_MAUX);
        const TValueType& r_point_value = AuxiliaryGet<TIsHistorical>(rPointVariable,rNode);
        AuxiliarySet<TIsHistorical>(rDistributedVariable,  TValueType(r_point_value / nodal_maux), rNode);
    });

    // Make sure that the initial approximation is the same between processes
    if (TIsHistorical) {
        rModelPart.GetCommunicator().SynchronizeVariable(rDistributedVariable);
    } else {
        rModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rDistributedVariable);
    }

    // Iteration: LumpedMass * delta_distributed = point_value - ConsistentMass * distributed_old
    double error_l2_norm = 0.0;
    unsigned int iteration = 0;
    for (std::size_t it = 0; it < MaximumIterations; ++it) {
        UpdateDistributionRHS<TIsHistorical,TContainerType,TPointNumber,TValueType>(rModelPart, rEntitiesContainer, rPointVariable, rDistributedVariable, mass_matrix);
        error_l2_norm = SolveDistributionIteration<TIsHistorical>(rModelPart,rDistributedVariable);
        error_l2_norm = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(error_l2_norm);
        error_l2_norm = std::sqrt(error_l2_norm);
        // Check convergence
        if (error_l2_norm <= Tolerance) {
            break;
        }
    }

    KRATOS_WARNING_IF("VariableRedistributionUtility", iteration == MaximumIterations)
        << "DistributePointValues did not converge in " << iteration << " iterations. L2 error norm: " << error_l2_norm << std::endl;
}

template< const bool TIsHistorical, class TValueType >
void VariableRedistributionUtility::DummySpecializedDistributePointValues(
    ModelPart& rModelPart,
    const Variable< TValueType >& rDistributedVariable,
    double Tolerance,
    unsigned int MaximumIterations)
{
    // Make sure that the initial approximation is the same between processes
    if (TIsHistorical) {
        rModelPart.GetCommunicator().SynchronizeVariable(rDistributedVariable);
    } else {
        rModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rDistributedVariable);
    }

    // Iteration: LumpedMass * delta_distributed = point_value - ConsistentMass * distributed_old
    double error_l2_norm = 0.0;
    unsigned int iteration = 0;
    while (iteration < MaximumIterations) {
        DummyUpdateDistributionRHS<TValueType>(rModelPart, rDistributedVariable);

        error_l2_norm = DummySolveDistributionIteration(rModelPart,rDistributedVariable);;
        error_l2_norm = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(error_l2_norm);

        // Check convergence
        iteration++;
        if (error_l2_norm <= Tolerance*Tolerance) {
            break;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

template<class TContainerType>
void VariableRedistributionUtility::ComputeNodalSizes(
    ModelPart& rModelPart,
    TContainerType& rEntitiesContainer)
{
    // Initialize NOAL_MAUX
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        rNode.SetValue(NODAL_MAUX, 0.0);
    });

    // Add the nodal contributions to NODAL_MAUX from each entity
    block_for_each(rEntitiesContainer, [&](typename TContainerType::value_type& rEntity){
        auto& r_geometry = rEntity.GetGeometry();
        const std::size_t n_nodes = r_geometry.PointsNumber();
        const double nodal_weight = r_geometry.DomainSize() / n_nodes;
        for (unsigned int j = 0; j < n_nodes; ++j) {
            double& r_lumped_mass = r_geometry[j].GetValue(NODAL_MAUX);
            AtomicAdd(r_lumped_mass, nodal_weight);
        }
    });

    // Assembe NODAL_MAUX values among processes
    rModelPart.GetCommunicator().AssembleNonHistoricalData(NODAL_MAUX);
}

template<>
void VariableRedistributionUtility::ConsistentMassMatrix< GeometryData::Kratos_Linear, 2 >(
    BoundedMatrix<double, 2, 2>& rMassMatrix)
{
    // First row
    rMassMatrix(0, 0) = 2.0/6.0;
    rMassMatrix(0, 1) = 1.0/6.0;

    // Second row
    rMassMatrix(1, 0) = 1.0/6.0;
    rMassMatrix(1, 1) = 2.0/6.0;
}


template<>
void VariableRedistributionUtility::ConsistentMassMatrix< GeometryData::Kratos_Triangle, 3 >(
    BoundedMatrix<double, 3, 3>& rMassMatrix)
{
    // First row
    rMassMatrix(0, 0) = 1.0/ 6.0;
    rMassMatrix(0, 1) = 1.0/12.0;
    rMassMatrix(0, 2) = 1.0/12.0;

    // Second row
    rMassMatrix(1, 0) = 1.0/12.0;
    rMassMatrix(1, 1) = 1.0/ 6.0;
    rMassMatrix(1, 2) = 1.0/12.0;

    // Third row
    rMassMatrix(2, 0) = 1.0/12.0;
    rMassMatrix(2, 1) = 1.0/12.0;
    rMassMatrix(2, 2) = 1.0/ 6.0;
}

template<>
void VariableRedistributionUtility::ConsistentMassMatrix< GeometryData::Kratos_Quadrilateral, 4 >(
    BoundedMatrix<double, 4, 4>& rMassMatrix)
{
    const double one_div_36 = 1.0 / 36.0;
    const double two_div_36 = 2.0 / 36.0;
    const double four_div_36 = 4.0 / 36.0;

    // First row
    rMassMatrix(0, 0) = four_div_36;
    rMassMatrix(0, 1) = two_div_36;
    rMassMatrix(0, 2) = one_div_36;
    rMassMatrix(0, 3) = two_div_36;

    // Second row
    rMassMatrix(1, 0) = two_div_36;
    rMassMatrix(1, 1) = four_div_36;
    rMassMatrix(1, 2) = two_div_36;
    rMassMatrix(1, 3) = one_div_36;

    // Third row
    rMassMatrix(2, 0) = one_div_36;
    rMassMatrix(2, 1) = two_div_36;
    rMassMatrix(2, 2) = four_div_36;
    rMassMatrix(2, 3) = two_div_36;

    // Forth row
    rMassMatrix(3, 0) = two_div_36;
    rMassMatrix(3, 1) = one_div_36;
    rMassMatrix(3, 2) = two_div_36;
    rMassMatrix(3, 3) = four_div_36;
}

template< const bool TIsHistorical, class TContainerType, unsigned int TNumNodes, class TValueType >
void VariableRedistributionUtility::UpdateDistributionRHS(
    ModelPart& rModelPart,
    TContainerType& rEntitiesContainer,
    const Variable< TValueType >& rPointVariable,
    const Variable< TValueType >& rDistributedVariable,
    BoundedMatrix<double, TNumNodes, TNumNodes>& rMassMatrix)
{
    const Variable<TValueType>& rhs_variable = GetRHSVariable(rDistributedVariable);
    const TValueType rhs_zero = rhs_variable.Zero(); // something of the correct type to initialize our values to zero

    // Reset the RHS container to zero
    block_for_each(rModelPart.Nodes(), rhs_zero, [&](NodeType& rNode, const TValueType& rRHSZero){
        rNode.SetValue(rhs_variable, rRHSZero);
    });

    // Calculate updated RHS
    block_for_each(rEntitiesContainer, [&](typename TContainerType::data_type& rEntity){
        auto& r_geometry = rEntity.GetGeometry();
        const double size = r_geometry.DomainSize();
        for (unsigned int j = 0; j < TNumNodes; j++) {
            TValueType rhs_j = rhs_zero;
            for (unsigned int k = 0; k < TNumNodes; k++) {
                rhs_j -= size * rMassMatrix(j,k) * AuxiliaryGet<TIsHistorical>(rDistributedVariable, r_geometry[k]);
            }
            ThreadsafeAdd(r_geometry[j].GetValue(rhs_variable), rhs_j);
        }
    });

    // Assemble distributed contributions
    rModelPart.GetCommunicator().AssembleNonHistoricalData(rhs_variable);

    // Add the nodal part of the RHS (the point-wise values)
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        rNode.GetValue(rhs_variable) += AuxiliaryGet<TIsHistorical>(rPointVariable, rNode);
    });
}

template< class TValueType >
void VariableRedistributionUtility::DummyUpdateDistributionRHS(
    ModelPart& rModelPart,
    const Variable< TValueType >& rDistributedVariable)
{
    // Assemble distributed contributions
    const auto& r_rhs_variable = GetRHSVariable(rDistributedVariable);
    rModelPart.GetCommunicator().AssembleNonHistoricalData(r_rhs_variable);
}

template< const bool TIsHistorical, class TValueType >
double VariableRedistributionUtility::SolveDistributionIteration(
    ModelPart& rModelPart,
    const Variable< TValueType >& rDistributedVariable)
{
    TValueType delta;
    double domain_size, error_l2_norm;
    const auto& r_rhs_variable = GetRHSVariable(rDistributedVariable);
    typedef CombinedReduction<SumReduction<double>,SumReduction<double>> TwoSumReduction;
    std::tie(domain_size, error_l2_norm) = block_for_each<TwoSumReduction>(rModelPart.Nodes(), delta, [&](NodeType& rNode, TValueType& rDelta){
        // Add correction to the current distributed nodal values
        const double size = rNode.GetValue(NODAL_MAUX);
        rDelta = rNode.GetValue(r_rhs_variable) / size;
        AuxiliaryGet<TIsHistorical>(rDistributedVariable, rNode) += rDelta;
        // Reduce error norm and
        return std::make_tuple(size, AddToNorm(rDelta, size));
    });

    return error_l2_norm /= domain_size;
}

template< class TValueType >
double VariableRedistributionUtility::DummySolveDistributionIteration(
    ModelPart& rModelPart,
    const Variable< TValueType >& rDistributedVariable)
{
    return 0.0;
}

template<>
const Variable<double>& VariableRedistributionUtility::GetRHSVariable<double>(const Variable<double>& rVariable)
{
    return MAPPER_SCALAR_PROJECTION_RHS;
}

template<>
const Variable< array_1d<double,3> >& VariableRedistributionUtility::GetRHSVariable< array_1d<double,3> >(const Variable< array_1d<double,3> >& rVariable)
{
    return MAPPER_VECTOR_PROJECTION_RHS;
}


template<>
double VariableRedistributionUtility::AddToNorm<double>(double NodalValue, double NodalSize)
{
    return NodalSize*NodalValue*NodalValue;
}


template<>
double VariableRedistributionUtility::AddToNorm< array_1d<double,3> >(array_1d<double,3> NodalValue, double NodalSize)
{
    return NodalSize*( NodalValue[0]*NodalValue[0]+NodalValue[1]*NodalValue[1]+NodalValue[2]*NodalValue[2] );
}

template<>
void VariableRedistributionUtility::ThreadsafeAdd<double>(
    double& rLHS,
    const double& rRHS)
{
    AtomicAdd<double>(rLHS,rRHS);
}

template<>
void VariableRedistributionUtility::ThreadsafeAdd<array_1d<double,3>>(
    array_1d<double,3>& rLHS,
    const array_1d<double,3>& rRHS)
{
    AtomicAdd<array_1d<double,3>, array_1d<double,3>>(rLHS, rRHS);
}

template<>
std::size_t VariableRedistributionUtility::NumberOfLocalEntities(
    const ModelPart& rModelPart,
    const ModelPart::ConditionsContainerType& rEntitiesContainer)
{
    return rModelPart.GetCommunicator().LocalMesh().NumberOfConditions();
}

template<>
std::size_t VariableRedistributionUtility::NumberOfLocalEntities(
    const ModelPart& rModelPart,
    const ModelPart::ElementsContainerType& rEntitiesContainer)
{
    return rModelPart.GetCommunicator().LocalMesh().NumberOfElements();
}

template<>
double& VariableRedistributionUtility::AuxiliaryGet<true, double>(
    const Variable<double>& rVariable,
    VariableRedistributionUtility::NodeType& rNode)
{
    return rNode.FastGetSolutionStepValue(rVariable);
}

template<>
array_1d<double,3>& VariableRedistributionUtility::AuxiliaryGet<true, array_1d<double,3>>(
    const Variable<array_1d<double,3>>& rVariable,
    VariableRedistributionUtility::NodeType& rNode)
{
    return rNode.FastGetSolutionStepValue(rVariable);
}

template<>
double& VariableRedistributionUtility::AuxiliaryGet<false, double>(
    const Variable<double>& rVariable,
    VariableRedistributionUtility::NodeType& rNode)
{
    return rNode.GetValue(rVariable);
}

template<>
array_1d<double,3>& VariableRedistributionUtility::AuxiliaryGet<false, array_1d<double,3>>(
    const Variable<array_1d<double,3>>& rVariable,
    VariableRedistributionUtility::NodeType& rNode)
{
    return rNode.GetValue(rVariable);
}

template<>
void VariableRedistributionUtility::AuxiliarySet<true, double>(
    const Variable<double>& rVariable,
    const double& rData,
    VariableRedistributionUtility::NodeType& rNode)
{
    rNode.FastGetSolutionStepValue(rVariable) = rData;
}

template<>
void VariableRedistributionUtility::AuxiliarySet<true, array_1d<double,3>>(
    const Variable<array_1d<double,3>>& rVariable,
    const array_1d<double,3>& rData,
    VariableRedistributionUtility::NodeType& rNode)
{
    rNode.FastGetSolutionStepValue(rVariable) = rData;
}

template<>
void VariableRedistributionUtility::AuxiliarySet<false, double>(
    const Variable<double>& rVariable,
    const double& rData,
    VariableRedistributionUtility::NodeType& rNode)
{
    rNode.SetValue(rVariable, rData);
}

template<>
void VariableRedistributionUtility::AuxiliarySet<false, array_1d<double,3>>(
    const Variable<array_1d<double,3>>& rVariable,
    const array_1d<double,3>& rData,
    VariableRedistributionUtility::NodeType& rNode)
{
    rNode.SetValue(rVariable, rData);
}

}
