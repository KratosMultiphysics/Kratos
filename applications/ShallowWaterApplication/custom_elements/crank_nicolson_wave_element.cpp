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
#include "custom_friction_laws/manning_law.h"
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

    LocalMatrixType lhs = ZeroMatrix(mLocalSize, mLocalSize);
    LocalVectorType rhs = ZeroVector(mLocalSize);

    // Struct to pass around the data
    ElementData data;
    InitializeData(data, rCurrentProcessInfo);
    GetNodalData(data, GetGeometry());

    BoundedMatrix<double,3,2> DN_DX; // Gradients matrix
    array_1d<double,3> N;            // Position of the gauss point
    double area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, area);

    CalculateGaussPointData(data, N);

    AddWaveTerms(lhs, rhs, data, N, DN_DX);
    AddFrictionTerms(lhs, rhs, data, N, DN_DX);

    // Substracting the Dirichlet term (since we use a residualbased approach)
    noalias(rhs) -= prod(lhs, data.unknown);

    noalias(rLeftHandSideMatrix) = area * lhs;
    noalias(rRightHandSideVector) = area * rhs;
}

template<std::size_t TNumNodes>
void CrankNicolsonWaveElement<TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR << "CrankNicolsonWaveElement: This element integrates in time and CalculateMassMatrix is no needed. Call CalculateLocalSystem instead." << std::endl;
}

template class CrankNicolsonWaveElement<3>;

} // namespace Kratos
