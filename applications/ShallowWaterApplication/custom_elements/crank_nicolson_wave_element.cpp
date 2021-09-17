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
#include "utilities/geometry_utilities.h"
#include "crank_nicolson_wave_element.h"

namespace Kratos
{

template<>
void CrankNicolsonWaveElement<3>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    if(rLeftHandSideMatrix.size1() != mLocalSize)
        rLeftHandSideMatrix.resize(mLocalSize, mLocalSize, false);

    if(rRightHandSideVector.size() != mLocalSize)
        rRightHandSideVector.resize(mLocalSize, false);

    // Integrated system
    LocalMatrixType lhs;
    LocalVectorType rhs;

    // Mass matrix
    LocalMatrixType m = ZeroMatrix(mLocalSize, mLocalSize);

    // Current time step
    LocalMatrixType k_0 = ZeroMatrix(mLocalSize, mLocalSize);
    LocalVectorType f_0 = ZeroVector(mLocalSize);
    LocalVectorType u_0;

    // Previous time step
    LocalMatrixType k_1 = ZeroMatrix(mLocalSize, mLocalSize);
    LocalVectorType f_1 = ZeroVector(mLocalSize);
    LocalVectorType u_1;

    // Struct to pass around the data
    ElementData data;
    InitializeData(data, rCurrentProcessInfo);
    const double dt_inv = 1.0 / rCurrentProcessInfo[DELTA_TIME];

    // Geometry data
    BoundedMatrix<double,3,2> DN_DX; // Gradients matrix
    array_1d<double,3> N;            // Position of the gauss point
    double area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, area);

    // Evaluate at the current time step
    GetNodalData(data, GetGeometry(), 0);
    CalculateGaussPointData(data, N);
    u_0 = this->GetUnknownVector(data);
    AddMassTerms(m, data, N, DN_DX);
    AddWaveTerms(k_0, f_0, data, N, DN_DX);
    AddFrictionTerms(k_0, f_0, data, N, DN_DX);

    // Evaluate at the previous time step
    GetNodalData(data, GetGeometry(), 1);
    CalculateGaussPointData(data, N);
    u_1 = this->GetUnknownVector(data);
    AddWaveTerms(k_1, f_1, data, N, DN_DX);
    AddFrictionTerms(k_1, f_1, data, N, DN_DX);

    // Assembly the system
    noalias(lhs) = dt_inv * m + 0.5 * k_0;
    noalias(rhs) = dt_inv * prod(m, u_1) - 0.5 * prod(k_1, u_1) + 0.5 * f_0 + 0.5 * f_1;

    // Substracting the Dirichlet term (since we use a residualbased approach)
    noalias(rhs) -= prod(lhs, u_0);

    noalias(rLeftHandSideMatrix) = area * lhs;
    noalias(rRightHandSideVector) = area * rhs;
}

template<std::size_t TNumNodes>
void CrankNicolsonWaveElement<TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR << "CrankNicolsonWaveElement: This element integrates in time and CalculateMassMatrix is not supported. Call CalculateLocalSystem instead." << std::endl;
}

template class CrankNicolsonWaveElement<3>;

} // namespace Kratos
