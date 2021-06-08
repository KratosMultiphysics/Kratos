#include "stokes_wall_condition.h"

namespace Kratos
{

///@name Specialized implementation of VMS for functions that depend on TDim
///@{

/**
 * @see StokesWallCondition::EquationIdVector
 */
template <>
void StokesWallCondition<2,2>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const unsigned int NumNodes = 2;
    const unsigned int LocalSize = 6;
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
    }
}

/**
 * @see StokesWallCondition::EquationIdVector
 */
template <>
void StokesWallCondition<3,3>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 12;
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
    }
}

/**
 * @see StokesWallCondition::EquationIdVector
 */
template <>
void StokesWallCondition<3,4>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const SizeType NumNodes = 4;
    const SizeType LocalSize = 16;
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
    }
}

/**
 * @see StokesWallCondition::GetDofList
 */
template <>
void StokesWallCondition<2,2>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const SizeType NumNodes = 2;
    const SizeType LocalSize = 6;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
    }
}

/**
 * @see StokesWallCondition::GetDofList
 */
template <>
void StokesWallCondition<3,3>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 12;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Z);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
    }
}


/**
 * @see StokesWallCondition::GetDofList
 */
template <>
void StokesWallCondition<3,4>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const SizeType NumNodes = 4;
    const SizeType LocalSize = 16;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Z);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
    }
}

template<unsigned int TDim, unsigned int TNumNodes>
void StokesWallCondition<TDim,TNumNodes>::ApplyNeumannCondition(MatrixType &rLocalMatrix, VectorType &rLocalVector)
{
    const unsigned int LocalSize = TDim+1;
    const GeometryType& rGeom = this->GetGeometry();
    const GeometryData::IntegrationMethod integration_method = static_cast<GeometryData::IntegrationMethod> (GeometryData::GI_GAUSS_2);
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(integration_method);
    const unsigned int NumGauss = IntegrationPoints.size();
    Vector gauss_pts_J_det = ZeroVector(NumGauss);
    rGeom.DeterminantOfJacobian(gauss_pts_J_det, integration_method);
    MatrixType NContainer = rGeom.ShapeFunctionsValues(integration_method);

    for (unsigned int g = 0; g < NumGauss; g++)
    {
        Vector N = row(NContainer,g);
        double Weight = gauss_pts_J_det[g] * IntegrationPoints[g].Weight();

        // Compute Unit normal
        array_1d<double, 3> Normal;
        Normal = rGeom.Normal(IntegrationPoints[g].Coordinates());
        double A = norm_2(Normal);
        Normal /= A;

        // Neumann boundary condition
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            //unsigned int row = i*LocalSize;
            const NodeType& rConstNode = this->GetGeometry()[i];
            const double pext = rConstNode.FastGetSolutionStepValue(EXTERNAL_PRESSURE);
            for (unsigned int j = 0; j < TNumNodes; j++)
            {
                unsigned int row = j*LocalSize;
                for (unsigned int d = 0; d < TDim;d++)
                    rLocalVector[row+d] -= Weight*N[j]*N[i]*pext*Normal[d];
            }
        }
    }
}

template class StokesWallCondition<2,2>;
template class StokesWallCondition<3,3>;
template class StokesWallCondition<3,4>;

} // namespace Kratos
