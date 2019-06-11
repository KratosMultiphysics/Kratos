// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Jordi Cotela
//

#include "adjoint_heat_diffusion_element.h"
#include "laplacian_element.h"

namespace Kratos
{

template<class PrimalElement>
AdjointHeatDiffusionElement<PrimalElement>::AdjointHeatDiffusionElement(
    IndexType NewId, typename GeometryType::Pointer pGeometry, Properties::Pointer pProperties):
    PrimalElement(NewId, pGeometry, pProperties)
{}

template<class PrimalElement>
AdjointHeatDiffusionElement<PrimalElement>::~AdjointHeatDiffusionElement() {}

template<class PrimalElement>
Element::Pointer AdjointHeatDiffusionElement<PrimalElement>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<AdjointHeatDiffusionElement<PrimalElement>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template<class PrimalElement>
Element::Pointer AdjointHeatDiffusionElement<PrimalElement>::Create(
    IndexType NewId,
    typename GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<AdjointHeatDiffusionElement<PrimalElement>>(NewId, pGeometry, pProperties);
}

template<class PrimalElement>
void AdjointHeatDiffusionElement<PrimalElement>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    const Geometry<Node<3>>& r_geom = this->GetGeometry();
    const unsigned int num_nodes = r_geom.PointsNumber();

    if (rLeftHandSideMatrix.size1() != num_nodes || rLeftHandSideMatrix.size2() != num_nodes)
    {
        rLeftHandSideMatrix.resize(num_nodes,num_nodes,false);
    }

    noalias(rLeftHandSideMatrix) = ZeroMatrix(num_nodes,num_nodes);
}

template class AdjointHeatDiffusionElement<LaplacianElement>;

}