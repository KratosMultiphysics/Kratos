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
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HEIGHT, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TOPOGRAPHY, r_node)

        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X, r_node)
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y, r_node)
        KRATOS_CHECK_DOF_IN_NODE(HEIGHT, r_node)
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

    const GeometryType& r_geom = GetGeometry();
    int counter = 0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rResult[counter++] = r_geom[i].GetDof(VELOCITY_X).EquationId();
        rResult[counter++] = r_geom[i].GetDof(VELOCITY_Y).EquationId();
        rResult[counter++] = r_geom[i].GetDof(HEIGHT).EquationId();
    }

    KRATOS_CATCH("")
}

template<std::size_t TNumNodes>
void WaveCondition<TNumNodes>::GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if(rConditionDofList.size() != mLocalSize)
        rConditionDofList.resize(mLocalSize);

    const GeometryType& r_geom = GetGeometry();
    int counter=0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rConditionDofList[counter++] = r_geom[i].pGetDof(VELOCITY_X);
        rConditionDofList[counter++] = r_geom[i].pGetDof(VELOCITY_Y);
        rConditionDofList[counter++] = r_geom[i].pGetDof(HEIGHT);
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
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(VELOCITY_X, Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(VELOCITY_Y, Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(HEIGHT, Step);
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
    KRATOS_ERROR << "WaveCondition: This method is not supported by the formulation" << std::endl;
}

template<std::size_t TNumNodes>
void WaveCondition<TNumNodes>::CalculateGeometryData(
    Vector &rGaussWeights,
    Matrix &rNContainer,
    ShapeFunctionsGradientsType &rDN_DXContainer) const
{
    Vector det_j_vector;
    const auto& r_geom = this->GetGeometry();
    const auto integration_method = r_geom.GetDefaultIntegrationMethod();
    r_geom.ShapeFunctionsIntegrationPointsGradients(rDN_DXContainer, det_j_vector, integration_method, rNContainer);

    const unsigned int number_of_gauss_points = r_geom.IntegrationPointsNumber(integration_method);
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints(integration_method);

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
    rData.gravity = rProcessInfo[GRAVITY_Z];
    rData.stab_factor = rProcessInfo[STABILIZATION_FACTOR];
    rData.relative_dry_height = rProcessInfo[RELATIVE_DRY_HEIGHT];

    auto& r_geom = GetGeometry();
    rData.length = r_geom.Length();

    for (IndexType i = 0; i < TNumNodes; i++)
    {
        const IndexType block = 3 * i;

        const auto h = r_geom[i].FastGetSolutionStepValue(HEIGHT);
        const array_1d<double,3> v = r_geom[i].FastGetSolutionStepValue(VELOCITY);

        rData.topography[i] = r_geom[i].FastGetSolutionStepValue(TOPOGRAPHY);

        rData.unknown[block]     = v[0];
        rData.unknown[block + 1] = v[1];
        rData.unknown[block + 2] = h;
    }
}

template<std::size_t TNumNodes>
void WaveCondition<TNumNodes>::CalculateGaussPointData(
    ConditionData& rData,
    const IndexType PointIndex,
    const array_1d<double,TNumNodes>& rN)
{
    rData.height = 0.0;
    rData.velocity = ZeroVector(3);

    for (IndexType i = 0; i < TNumNodes; i++)
    {
        const IndexType block = 3 * i;

        const auto h = rData.unknown[block + 2];
        array_1d<double,3> v;
        v[0] = rData.unknown[block];
        v[1] = rData.unknown[block + 1];
        v[2] = 0.0;

        rData.height += rN[i] * h;
        rData.velocity += rN[i] * v;
    }
    auto integration_point = GetGeometry().IntegrationPoints()[PointIndex];
    rData.normal = GetGeometry().Normal(integration_point);
    rData.normal /= norm_2(rData.normal);
}

template<std::size_t TNumNodes>
void WaveCondition<TNumNodes>::AddWaveTerms(
    LocalMatrixType& rMatrix,
    LocalVectorType& rVector,
    const ConditionData& rData,
    const array_1d<double,TNumNodes>& rN,
    const double Weight)
{
    const bool integrate_by_parts = false;
    const auto h = rData.height;
    const auto g = rData.gravity;
    const auto z = rData.topography;
    const auto n = rData.normal;

    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        const IndexType i_block = 3 * i;
        for (IndexType j = 0; j < TNumNodes; ++j)
        {
            const IndexType j_block = 3 * j;

            double n_ij;
            if (integrate_by_parts) {
                n_ij = rN[i] * rN[j];
            } else {
                n_ij = 0.0;
            }

            /* First component
            * A_1 = {{ 0   0   g },
            *        { 0   0   0 },
            *        { h   0   0 }}
            */
            rMatrix(i_block,     j_block + 2) -= Weight * n_ij * g * n[0];
            rMatrix(i_block + 2, j_block)     -= Weight * n_ij * h * n[0];
            rVector(i_block)                  += Weight * n_ij * g * n[0] * z[j];

            /* Second component
            * A_2 = {{ 0   0   0 },
            *        { 0   0   g },
            *        { 0   h   0 }}
            */
            rMatrix(i_block + 1, j_block + 2) -= Weight * n_ij * g * n[1];
            rMatrix(i_block + 2, j_block + 1) -= Weight * n_ij * h * n[1];
            rVector(i_block + 1)              += Weight * n_ij * g * n[1] * z[j];
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
    const double h = rData.height;
    const double g = rData.gravity;
    const auto n = rData.normal;
    const auto inv_c = std::sqrt(InverseHeight(rData) / g);
    const auto l = StabilizationParameter(rData);

    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        const IndexType i_block = 3 * i;
        for (IndexType j = 0; j < TNumNodes; ++j)
        {
            const IndexType j_block = 3 * j;
            const double n_ij = rN[i] * rN[j];

            /* Stabilization x
             * l / sqrt(gh) * A1 * n1
             */
            rMatrix(i_block,     j_block + 2) -= Weight * n_ij * l * g * inv_c * n[0];
            rMatrix(i_block + 2, j_block)     -= Weight * n_ij * l * h * inv_c * n[0];

            /* Stabilization y
             * l / sqrt(gh) * A2 * n2
             */
            rMatrix(i_block + 1, j_block + 2) -= Weight * n_ij * l * g * inv_c * n[1];
            rMatrix(i_block + 2, j_block + 1) -= Weight * n_ij * l * h * inv_c * n[1];
        }
    }
}

template<std::size_t TNumNodes>
double WaveCondition<TNumNodes>::StabilizationParameter(const ConditionData& rData) const
{
    return rData.length * rData.stab_factor;
}

template<std::size_t TNumNodes>
double WaveCondition<TNumNodes>::InverseHeight(const ConditionData& rData) const
{
    const double height = rData.height;
    const double epsilon = rData.relative_dry_height * rData.length;
    return ShallowWaterUtilities().InverseHeight(height, epsilon);
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
    ShapeFunctionsGradientsType DN_DX_container;
    CalculateGeometryData(weights, N_container, DN_DX_container);
    const IndexType num_gauss_points = weights.size();

    for (IndexType g = 0; g < num_gauss_points; ++g)
    {
        const double weight = weights[g];
        const array_1d<double,TNumNodes> N = row(N_container, g);
        CalculateGaussPointData(data, g, N);
        AddWaveTerms(lhs, rhs, data, N, weight);
        AddFluxTerms(rhs, data, N, weight);
    }
    noalias(rhs) -= prod(lhs, data.unknown);

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
    ShapeFunctionsGradientsType DN_DX_container;
    CalculateGeometryData(weights, N_container, DN_DX_container);
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
