//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "wave_condition.h"
#include "includes/checks.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "shallow_water_application_variables.h"
#include "custom_utilities/shallow_water_utilities.h"

namespace Kratos
{

template<std::size_t TNumNodes>
int WaveCondition<TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int err = Condition::Check(rCurrentProcessInfo);
    if (err != 0) return err;

    // Check that the condition's nodes contain all required SolutionStepData and Degrees of freedom
    for (const auto& r_node : this->GetGeometry())
    {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MOMENTUM, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HEIGHT, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TOPOGRAPHY, r_node)

        KRATOS_CHECK_DOF_IN_NODE(this->GetUnknownComponent(0), r_node)
        KRATOS_CHECK_DOF_IN_NODE(this->GetUnknownComponent(1), r_node)
        KRATOS_CHECK_DOF_IN_NODE(this->GetUnknownComponent(2), r_node)
    }

    return err;

    KRATOS_CATCH("")
}

template<std::size_t TNumNodes>
void WaveCondition<TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if(rResult.size() != mLocalSize)
        rResult.resize(mLocalSize, false); // False says not to preserve existing storage!!

    const GeometryType& r_geom = this->GetGeometry();
    int counter = 0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rResult[counter++] = r_geom[i].GetDof(this->GetUnknownComponent(0)).EquationId();
        rResult[counter++] = r_geom[i].GetDof(this->GetUnknownComponent(1)).EquationId();
        rResult[counter++] = r_geom[i].GetDof(this->GetUnknownComponent(2)).EquationId();
    }

    KRATOS_CATCH("")
}

template<std::size_t TNumNodes>
void WaveCondition<TNumNodes>::GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if(rConditionDofList.size() != mLocalSize)
        rConditionDofList.resize(mLocalSize);

    const GeometryType& r_geom = this->GetGeometry();
    int counter=0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rConditionDofList[counter++] = r_geom[i].pGetDof(this->GetUnknownComponent(0));
        rConditionDofList[counter++] = r_geom[i].pGetDof(this->GetUnknownComponent(1));
        rConditionDofList[counter++] = r_geom[i].pGetDof(this->GetUnknownComponent(2));
    }

    KRATOS_CATCH("")
}

template<std::size_t TNumNodes>
void WaveCondition<TNumNodes>::GetValuesVector(Vector& rValues, int Step) const
{
    if (rValues.size() != mLocalSize)
        rValues.resize(mLocalSize, false);

    const GeometryType& r_geom = this->GetGeometry();
    IndexType counter = 0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(this->GetUnknownComponent(0), Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(this->GetUnknownComponent(1), Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(this->GetUnknownComponent(2), Step);
    }
}

template<std::size_t TNumNodes>
void WaveCondition<TNumNodes>::GetFirstDerivativesVector(Vector& rValues, int Step) const
{
    if (rValues.size() != mLocalSize)
        rValues.resize(mLocalSize, false);

    const GeometryType& r_geom = this->GetGeometry();
    IndexType counter = 0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(ACCELERATION_X, Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(ACCELERATION_Y, Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(VERTICAL_VELOCITY, Step);
    }
}

template<std::size_t TNumNodes>
void WaveCondition<TNumNodes>::GetSecondDerivativesVector(Vector& rValues, int Step) const
{
    KRATOS_ERROR << "WaveCondition::GetSecondDerivativesVector This method is not supported by the formulation" << std::endl;
}

template<std::size_t TNumNodes>
const Variable<double>& WaveCondition<TNumNodes>::GetUnknownComponent(int Index) const
{
    switch (Index) {
        case 0: return VELOCITY_X;
        case 1: return VELOCITY_Y;
        case 2: return HEIGHT;
        default: KRATOS_ERROR << "WaveCondition::GetUnknownComponent index out of bounds." << std::endl;
    }
}

template<std::size_t TNumNodes>
typename WaveCondition<TNumNodes>::LocalVectorType WaveCondition<TNumNodes>::GetUnknownVector(ConditionData& rData)
{
    std::size_t index = 0;
    array_1d<double,mLocalSize> unknown;
    for (std::size_t i = 0; i < TNumNodes; ++i) {
        unknown[index++] = rData.nodal_v[i][0];
        unknown[index++] = rData.nodal_v[i][1];
        unknown[index++] = rData.nodal_h[i];
    }
    return unknown;
}

template<std::size_t TNumNodes>
void WaveCondition<TNumNodes>::CalculateGeometryData(Vector &rGaussWeights, Matrix &rNContainer) const
{
    Vector det_j_vector;
    const auto& r_geom = this->GetGeometry();
    const auto integration_method = r_geom.GetDefaultIntegrationMethod();

    rNContainer = r_geom.ShapeFunctionsValues(integration_method);

    const unsigned int number_of_gauss_points = r_geom.IntegrationPointsNumber(integration_method);
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints(integration_method);
    r_geom.DeterminantOfJacobian(det_j_vector, integration_method);

    if (rGaussWeights.size() != number_of_gauss_points)
        rGaussWeights.resize(number_of_gauss_points, false);

    for (unsigned int g = 0; g < number_of_gauss_points; ++g)
        rGaussWeights[g] = det_j_vector[g] * integration_points[g].Weight();
}

template<std::size_t TNumNodes>
void WaveCondition<TNumNodes>::InitializeData(
    ConditionData& rData,
    const ProcessInfo& rProcessInfo)
{
    rData.integrate_by_parts = rProcessInfo[INTEGRATE_BY_PARTS]; //since it is passed as const it will return false if it doesn't have INTEGRATE_BY_PARTS
    rData.gravity = rProcessInfo[GRAVITY_Z];
    rData.stab_factor = rProcessInfo[STABILIZATION_FACTOR];
    rData.relative_dry_height = rProcessInfo[RELATIVE_DRY_HEIGHT];

    auto& r_geom = this->GetGeometry();
    rData.length = r_geom.Length();

    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rData.nodal_h[i] = r_geom[i].FastGetSolutionStepValue(HEIGHT);
        rData.nodal_z[i] = r_geom[i].FastGetSolutionStepValue(TOPOGRAPHY);
        rData.nodal_v[i] = r_geom[i].FastGetSolutionStepValue(VELOCITY);
        rData.nodal_q[i] = r_geom[i].FastGetSolutionStepValue(MOMENTUM);
    }
}

template<std::size_t TNumNodes>
void WaveCondition<TNumNodes>::CalculateGaussPointData(
    ConditionData& rData,
    const IndexType PointIndex,
    const array_1d<double,TNumNodes>& rN)
{
    const double h = inner_prod(rData.nodal_h, rN);
    const array_1d<double,3> v = VectorProduct(rData.nodal_v, rN);

    rData.height = h;
    rData.velocity = v;

    rData.A1 = ZeroMatrix(3, 3);
    rData.A1(0,2) = rData.gravity;
    rData.A1(2,0) = h;

    rData.A2 = ZeroMatrix(3, 3);
    rData.A2(1,2) = rData.gravity;
    rData.A2(2,1) = h;

    rData.b1 = ZeroVector(3);
    rData.b1[0] = rData.gravity;

    rData.b2 = ZeroVector(3);
    rData.b2[1] = rData.gravity;

    auto integration_point = this->GetGeometry().IntegrationPoints()[PointIndex];
    rData.normal = this->GetGeometry().Normal(integration_point);
    rData.normal /= norm_2(rData.normal);
}

template<std::size_t TNumNodes>
const array_1d<double,3> WaveCondition<TNumNodes>::VectorProduct(
    const array_1d<array_1d<double,3>,TNumNodes>& rV,
    const array_1d<double,TNumNodes>& rN) const
{
    array_1d<double,3> result = ZeroVector(3);
    for (std::size_t i = 0; i < TNumNodes; ++i)
    {
        result += rV[i] * rN[i];
    }
    return result;
}

template<std::size_t TNumNodes>
void WaveCondition<TNumNodes>::AddWaveTerms(
    LocalMatrixType& rMatrix,
    LocalVectorType& rVector,
    const ConditionData& rData,
    const array_1d<double,TNumNodes>& rN,
    const double Weight)
{
    const auto z = rData.nodal_z;
    const auto n = rData.normal;

    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        for (IndexType j = 0; j < TNumNodes; ++j)
        {
            double n_ij;
            if (rData.integrate_by_parts) {
                n_ij = rN[i] * rN[j];
            } else {
                n_ij = 0.0;
            }

            /// First component
            MathUtils<double>::AddMatrix(rMatrix,  Weight*n_ij*rData.A1*n[0], 3*i, 3*j);
            MathUtils<double>::AddVector(rVector, -Weight*n_ij*rData.b1*n[0]*z[j], 3*i);

            /// Second component
            MathUtils<double>::AddMatrix(rMatrix,  Weight*n_ij*rData.A2*n[1], 3*i, 3*j);
            MathUtils<double>::AddVector(rVector, -Weight*n_ij*rData.b2*n[1]*z[j], 3*i);
        }
    }
}

template<std::size_t TNumNodes>
void WaveCondition<TNumNodes>::AddFluxTerms(
    LocalVectorType& rVector,
    const ConditionData& rData,
    const array_1d<double,TNumNodes>& rN,
    const double Weight)
{
}

template<std::size_t TNumNodes>
void WaveCondition<TNumNodes>::AddMassTerms(
    LocalMatrixType& rMatrix,
    const ConditionData& rData,
    const array_1d<double,TNumNodes>& rN,
    const double Weight)
{
}

template<std::size_t TNumNodes>
void WaveCondition<TNumNodes>::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    if(rRightHandSideVector.size() != mLocalSize)
        rRightHandSideVector.resize(mLocalSize, false);

    MatrixType lhs = ZeroMatrix(mLocalSize, mLocalSize);

    CalculateLocalSystem(lhs, rRightHandSideVector, rCurrentProcessInfo);
}

template<std::size_t TNumNodes>
void WaveCondition<TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    if(rLeftHandSideMatrix.size1() != mLocalSize)
        rLeftHandSideMatrix.resize(mLocalSize, mLocalSize, false);

    if(rRightHandSideVector.size() != mLocalSize)
        rRightHandSideVector.resize(mLocalSize, false);
    
    LocalMatrixType lhs = ZeroMatrix(mLocalSize, mLocalSize);
    LocalVectorType rhs = ZeroVector(mLocalSize);

    ConditionData data;
    InitializeData(data, rCurrentProcessInfo);

    Vector weights;
    Matrix N_container;
    CalculateGeometryData(weights, N_container);
    const IndexType num_gauss_points = weights.size();

    for (IndexType g = 0; g < num_gauss_points; ++g)
    {
        const double weight = weights[g];
        const array_1d<double,TNumNodes> N = row(N_container, g);
        CalculateGaussPointData(data, g, N);
        AddWaveTerms(lhs, rhs, data, N, weight);
        AddFluxTerms(rhs, data, N, weight);
    }
    noalias(rhs) -= prod(lhs, this->GetUnknownVector(data));

    noalias(rLeftHandSideMatrix) = lhs;
    noalias(rRightHandSideVector) = rhs;
}

template<std::size_t TNumNodes>
void WaveCondition<TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    if(rMassMatrix.size1() != mLocalSize)
        rMassMatrix.resize(mLocalSize, mLocalSize, false);

    LocalMatrixType m = ZeroMatrix(mLocalSize, mLocalSize);

    ConditionData data;
    InitializeData(data, rCurrentProcessInfo);

    Vector weights;
    Matrix N_container;
    CalculateGeometryData(weights, N_container);
    const IndexType num_gauss_points = weights.size();

    for (IndexType g = 0; g < num_gauss_points; ++g)
    {
        const double weight = weights[g];
        const array_1d<double,TNumNodes> N = row(N_container, g);
        CalculateGaussPointData(data, g, N);
        AddMassTerms(m, data, N, weight);
    }
    noalias(rMassMatrix) = m;
}

template class WaveCondition<2>;

} // namespace Kratos
