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
#include "variable_redistribution_utility.h"

namespace Kratos
{

void VariableRedistributionUtility::ConvertDistributedValuesToPoint(
    ModelPart& rModelPart,
    const Variable<double>& rDistributedVariable,
    const Variable<double>& rPointVariable)
{
    CallSpecializedConvertDistributedValuesToPoint(rModelPart,rDistributedVariable,rPointVariable);
}

void VariableRedistributionUtility::DistributePointValues(
    ModelPart& rModelPart,
    const Variable<double>& rPointVariable,
    const Variable<double>& rDistributedVariable,
    double Tolerance,
    unsigned int MaximumIterations)
{
    CallSpecializedDistributePointValues(rModelPart,rPointVariable,rDistributedVariable, Tolerance, MaximumIterations);
}

void VariableRedistributionUtility::ConvertDistributedValuesToPoint(
    ModelPart& rModelPart,
    const Variable< array_1d<double,3> >& rDistributedVariable,
    const Variable< array_1d<double,3> >& rPointVariable)
{
    CallSpecializedConvertDistributedValuesToPoint(rModelPart,rDistributedVariable,rPointVariable);
}

void VariableRedistributionUtility::DistributePointValues(
    ModelPart& rModelPart,
    const Variable< array_1d<double,3> >& rPointVariable,
    const Variable< array_1d<double,3> >& rDistributedVariable,
    double Tolerance,
    unsigned int MaximumIterations)
{
    CallSpecializedDistributePointValues(rModelPart,rPointVariable,rDistributedVariable, Tolerance, MaximumIterations);
}

template< class TValueType >
void VariableRedistributionUtility::CallSpecializedConvertDistributedValuesToPoint(
    ModelPart& rModelPart,
    const Variable<TValueType>& rDistributedVariable,
    const Variable<TValueType>& rPointVariable)
{
    // Check if there is any condition in the current partition
    const unsigned int n_loc_conds = rModelPart.GetCommunicator().LocalMesh().NumberOfConditions();
    int n_tot_conds = n_loc_conds;
    rModelPart.GetCommunicator().SumAll(n_tot_conds);

    // If there is conditions, this function dispatches the call to the correct specialization
    if (n_tot_conds != 0){
        if (n_loc_conds != 0){
            Geometry< Node<3> >& rReferenceGeometry = rModelPart.ConditionsBegin()->GetGeometry();
            const GeometryData::KratosGeometryFamily GeometryFamily = rReferenceGeometry.GetGeometryFamily();
            const unsigned int PointsNumber = rReferenceGeometry.PointsNumber();

            if (GeometryFamily == GeometryData::Kratos_Linear && PointsNumber == 2){
                VariableRedistributionUtility::SpecializedConvertDistributedValuesToPoint<GeometryData::Kratos_Linear,2,TValueType>(
                    rModelPart,
                    rDistributedVariable,
                    rPointVariable);
            } else if (GeometryFamily == GeometryData::Kratos_Triangle && PointsNumber == 3) {
                VariableRedistributionUtility::SpecializedConvertDistributedValuesToPoint<GeometryData::Kratos_Triangle,3,TValueType>(
                    rModelPart,
                    rDistributedVariable,
                    rPointVariable);
            } else {
                KRATOS_ERROR << "Unsupported geometry type with " << PointsNumber << " points." << std::endl;
            }
        } else {
            VariableRedistributionUtility::DummySpecializedConvertDistributedValuesToPoint<TValueType>(
                rModelPart,
                rDistributedVariable,
                rPointVariable);
        }
    }
}

template< class TValueType >
void VariableRedistributionUtility::CallSpecializedDistributePointValues(
    ModelPart& rModelPart,
    const Variable<TValueType>& rPointVariable,
    const Variable<TValueType>& rDistributedVariable,
    double Tolerance,
    unsigned int MaximumIterations)
{
    // Check if there is any condition in the current partition
    const unsigned int n_loc_conds = rModelPart.GetCommunicator().LocalMesh().NumberOfConditions();
    int n_tot_conds = n_loc_conds;
    rModelPart.GetCommunicator().SumAll(n_tot_conds);

    // If there is conditions, this function dispatches the call to the correct specialization
    if (n_tot_conds != 0){
        if (n_loc_conds != 0){
            Geometry< Node<3> >& rReferenceGeometry = rModelPart.ConditionsBegin()->GetGeometry();
            const GeometryData::KratosGeometryFamily GeometryFamily = rReferenceGeometry.GetGeometryFamily();
            const unsigned int PointsNumber = rReferenceGeometry.PointsNumber();

            if (GeometryFamily == GeometryData::Kratos_Linear && PointsNumber == 2){
                VariableRedistributionUtility::SpecializedDistributePointValues<GeometryData::Kratos_Linear,2,TValueType>(
                    rModelPart,
                    rPointVariable,
                    rDistributedVariable,
                    Tolerance,
                    MaximumIterations);
            } else if (GeometryFamily == GeometryData::Kratos_Triangle && PointsNumber == 3){
                VariableRedistributionUtility::SpecializedDistributePointValues<GeometryData::Kratos_Triangle,3,TValueType>(
                    rModelPart,
                    rPointVariable,
                    rDistributedVariable,
                    Tolerance,
                    MaximumIterations);
            } else {
                KRATOS_ERROR << "Unsupported geometry type with " << PointsNumber << " points." << std::endl;
            }
        } else {
            VariableRedistributionUtility::DummySpecializedDistributePointValues<TValueType>(
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

template< GeometryData::KratosGeometryFamily TFamily, unsigned int TPointNumber, class TValueType >
void VariableRedistributionUtility::SpecializedConvertDistributedValuesToPoint(
    ModelPart& rModelPart,
    const Variable< TValueType >& rDistributedVariable,
    const Variable< TValueType >& rPointVariable)
{
    const int number_of_nodes_in_model_part = rModelPart.NumberOfNodes();
    TValueType Zero = rPointVariable.Zero();

    // Initialize result to zero
    #pragma omp parallel for
    for (int i = 0; i < number_of_nodes_in_model_part; i++)
    {
        ModelPart::NodeIterator node_iter = rModelPart.NodesBegin() + i;
        node_iter->FastGetSolutionStepValue(rPointVariable) = Zero;
    }

    BoundedMatrix<double,TPointNumber,TPointNumber> MassMatrix;
    ConsistentMassMatrix<TFamily,TPointNumber>(MassMatrix);

    // Make sure that the distributed values are equal between processors
    rModelPart.GetCommunicator().SynchronizeVariable(rDistributedVariable);

    const int number_of_conditions = rModelPart.NumberOfConditions();
    #pragma omp for
    for (int i = 0; i < number_of_conditions; i++)
    {
        auto condition_iterator = rModelPart.ConditionsBegin()+i;
        Geometry< Node<3> >& r_geometry = condition_iterator->GetGeometry();
        const double size = r_geometry.DomainSize();

        for (unsigned int j = 0; j < TPointNumber; j++)
        {
            TValueType value_j = rDistributedVariable.Zero();

            for (unsigned int k = 0; k < TPointNumber; k++)
            {
                value_j += size * MassMatrix(j,k) * r_geometry[k].FastGetSolutionStepValue(rDistributedVariable);
            }

            ThreadsafeAdd(r_geometry[j].FastGetSolutionStepValue(rPointVariable), value_j);
        }
    }

    // Add the contributions between processors
    rModelPart.GetCommunicator().AssembleCurrentData(rPointVariable);
}

template< class TValueType >
void VariableRedistributionUtility::DummySpecializedConvertDistributedValuesToPoint(
    ModelPart& rModelPart,
    const Variable< TValueType >& rDistributedVariable,
    const Variable< TValueType >& rPointVariable)
{
    // Make sure that the distributed values are equal between processors
    rModelPart.GetCommunicator().SynchronizeVariable(rDistributedVariable);

    // Add the contributions between processors
    rModelPart.GetCommunicator().AssembleCurrentData(rPointVariable);
}

///////////////////////////////////////////////////////////////////////////////

template< GeometryData::KratosGeometryFamily TFamily, unsigned int TPointNumber, class TValueType >
void VariableRedistributionUtility::SpecializedDistributePointValues(
    ModelPart& rModelPart,
    const Variable< TValueType >& rPointVariable,
    const Variable< TValueType >& rDistributedVariable,
    double Tolerance,
    unsigned int MaximumIterations)
{
    ComputeNodalSizes(rModelPart);

    BoundedMatrix< double, TPointNumber, TPointNumber > mass_matrix;
    ConsistentMassMatrix<TFamily,TPointNumber>(mass_matrix);

    const int number_of_nodes_in_model_part = rModelPart.NumberOfNodes();

    // Initial guess (NodalValue / NodalSize)
    #pragma omp parallel for
    for (int i = 0; i < number_of_nodes_in_model_part; i++)
    {
        ModelPart::NodeIterator node_iter = rModelPart.NodesBegin() + i;
        node_iter->FastGetSolutionStepValue(rDistributedVariable) = node_iter->FastGetSolutionStepValue(rPointVariable) / node_iter->GetValue(NODAL_MAUX);
    }

    // Make sure that the initial approximation is the same between processes
    rModelPart.GetCommunicator().SynchronizeVariable(rDistributedVariable);

    // Iteration: LumpedMass * delta_distributed = point_value - ConsistentMass * distributed_old
    double error_l2_norm = 0.0;
    unsigned int iteration = 0;
    while ( iteration < MaximumIterations )
    {
        UpdateDistributionRHS<TPointNumber,TValueType>(rModelPart,rPointVariable, rDistributedVariable, mass_matrix);

        error_l2_norm = SolveDistributionIteration(rModelPart,rDistributedVariable);
        rModelPart.GetCommunicator().SumAll(error_l2_norm);

        // Check convergence
        iteration++;
        if (error_l2_norm <= Tolerance*Tolerance)
        {
            break;
        }
    }

    KRATOS_WARNING_IF("VariableRedistributionUtility", iteration == MaximumIterations)
        << "DistributePointValues did not converge in " << iteration << " iterations. L2 error norm: " << error_l2_norm << std::endl;
}

template< class TValueType >
void VariableRedistributionUtility::DummySpecializedDistributePointValues(
    ModelPart& rModelPart,
    const Variable< TValueType >& rDistributedVariable,
    double Tolerance,
    unsigned int MaximumIterations)
{
    // Make sure that the initial approximation is the same between processes
    rModelPart.GetCommunicator().SynchronizeVariable(rDistributedVariable);

    // Iteration: LumpedMass * delta_distributed = point_value - ConsistentMass * distributed_old
    double error_l2_norm = 0.0;
    unsigned int iteration = 0;
    while ( iteration < MaximumIterations )
    {
        DummyUpdateDistributionRHS<TValueType>(rModelPart, rDistributedVariable);

        rModelPart.GetCommunicator().SumAll(error_l2_norm);

        // Check convergence
        iteration++;
        if (error_l2_norm <= Tolerance*Tolerance){
            break;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void VariableRedistributionUtility::ComputeNodalSizes(ModelPart& rModelPart)
{
    const int number_of_nodes = rModelPart.NumberOfNodes();
    #pragma omp for
    for (int i = 0; i < number_of_nodes; i++)
    {
        auto node_iterator = rModelPart.NodesBegin()+i;
        node_iterator->SetValue(NODAL_MAUX,0.0);
    }

    const int number_of_conditions = rModelPart.NumberOfConditions();
    #pragma omp for
    for (int i = 0; i < number_of_conditions; i++)
    {
        auto condition_iterator = rModelPart.ConditionsBegin()+i;
        Geometry< Node<3> >& r_geometry = condition_iterator->GetGeometry();
        const double condition_size = r_geometry.DomainSize();
        const unsigned int nodes_in_condition = r_geometry.PointsNumber();
        const double nodal_factor = 1.0 / double(nodes_in_condition);
        for (unsigned int j = 0; j < nodes_in_condition; j++)
        {
            double& lumped_mass = r_geometry[j].GetValue(NODAL_MAUX);
            #pragma omp atomic
            lumped_mass += nodal_factor * condition_size;
        }
    }

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

template< unsigned int TNumNodes, class TValueType >
void VariableRedistributionUtility::UpdateDistributionRHS(
    ModelPart& rModelPart,
    const Variable< TValueType >& rPointVariable,
    const Variable< TValueType >& rDistributedVariable,
    BoundedMatrix<double, TNumNodes, TNumNodes>& rMassMatrix)
{
    const Variable<TValueType>& rhs_variable = GetRHSVariable(rDistributedVariable);
    const TValueType rhs_zero = rhs_variable.Zero(); // something of the correct type to initialize our values to zero

    const int number_of_nodes_in_model_part = rModelPart.NumberOfNodes();
    const int number_of_conditions_in_model_part = rModelPart.NumberOfConditions();

    // Reset the RHS container to zero
    #pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes_in_model_part; i_node++ )
    {
        ModelPart::NodesContainerType::iterator node_iter = rModelPart.NodesBegin() + i_node;
        node_iter->SetValue(rhs_variable, rhs_zero );
    }

    // Calculate updated RHS
    #pragma omp parallel for
    for (int i_condition = 0; i_condition < number_of_conditions_in_model_part; i_condition++)
    {
        auto condition_iterator = rModelPart.ConditionsBegin() + i_condition;
        Geometry< Node<3> >& r_geometry = condition_iterator->GetGeometry();

        const double size = r_geometry.DomainSize();

        for (unsigned int j = 0; j < TNumNodes; j++)
        {
            TValueType rhs_j = rhs_zero;

            for (unsigned int k = 0; k < TNumNodes; k++)
            {
                rhs_j -= size * rMassMatrix(j,k) * r_geometry[k].FastGetSolutionStepValue(rDistributedVariable);
            }

            ThreadsafeAdd(r_geometry[j].GetValue(rhs_variable), rhs_j);
        }
    }

    // Assemble distributed contributions
    rModelPart.GetCommunicator().AssembleNonHistoricalData(rhs_variable);

    // Add the nodal part of the RHS (the point-wise values)
    #pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes_in_model_part; i_node++ )
    {
        ModelPart::NodesContainerType::iterator node_iter = rModelPart.NodesBegin() + i_node;
        node_iter->GetValue(rhs_variable) += node_iter->FastGetSolutionStepValue(rPointVariable);
    }
}

template< class TValueType >
void VariableRedistributionUtility::DummyUpdateDistributionRHS(
    ModelPart& rModelPart,
    const Variable< TValueType >& rDistributedVariable)
{
    const Variable<TValueType>& rhs_variable = GetRHSVariable(rDistributedVariable);
    // Assemble distributed contributions
    rModelPart.GetCommunicator().AssembleNonHistoricalData(rhs_variable);
}

template< class TValueType >
double VariableRedistributionUtility::SolveDistributionIteration(
    ModelPart& rModelPart,
    const Variable< TValueType >& rDistributedVariable)
{

    const int number_of_nodes_in_model_part = rModelPart.NumberOfNodes();

    const Variable<TValueType>& rhs_variable = GetRHSVariable(rDistributedVariable);
    TValueType delta = rDistributedVariable.Zero();
    double error_l2_norm = 0.0;
    double domain_size = 0.0;

    #pragma omp parallel for reduction(+: error_l2_norm, domain_size) private(delta)
    for (int i_node = 0; i_node < number_of_nodes_in_model_part; i_node++)
    {
        auto node_iter = rModelPart.NodesBegin() + i_node;

        const double size = node_iter->GetValue(NODAL_MAUX);
        delta = node_iter->GetValue(rhs_variable) / size;

        TValueType& r_updated_value = node_iter->FastGetSolutionStepValue(rDistributedVariable);
        r_updated_value += delta;

        error_l2_norm += AddToNorm(delta, size);
        domain_size += size;
    }

    return error_l2_norm /= domain_size;
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
void VariableRedistributionUtility::ThreadsafeAdd<double>(double& rLHS, const double& rRHS)
{
    #pragma omp atomic
    rLHS += rRHS;
}

template<>
void VariableRedistributionUtility::ThreadsafeAdd< array_1d<double,3> >(array_1d<double,3>& rLHS, const array_1d<double,3>& rRHS)
{
    for (unsigned int i = 0; i < 3; i++)
    {
        #pragma omp atomic
        rLHS[i] += rRHS[i];
    }
}

}
