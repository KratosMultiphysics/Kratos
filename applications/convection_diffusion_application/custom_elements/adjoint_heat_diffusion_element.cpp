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

#include "convection_diffusion_application_variables.h"

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

template<class PrimalElement>
void AdjointHeatDiffusionElement<PrimalElement>::GetValuesVector(Vector& rValues, int Step)
{
    const GeometryType& r_geom = this->GetGeometry();
    const unsigned int num_nodes = r_geom.PointsNumber();

    if (rValues.size() != num_nodes)
    {
        rValues.resize(num_nodes,false);
    }

    for (unsigned int i = 0; i < num_nodes; i++)
    {
        rValues[i] = r_geom[i].FastGetSolutionStepValue(ADJOINT_HEAT_TRANSFER, Step);
    }
}

template<class PrimalElement>
void AdjointHeatDiffusionElement<PrimalElement>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& r_geom = this->GetGeometry();
    const unsigned int num_nodes = r_geom.PointsNumber();

    if (rResult.size() != num_nodes)
    {
        rResult.resize(num_nodes,false);
    }

    for (unsigned int i = 0; i < num_nodes; i++)
    {
        rResult[i] = r_geom[i].GetDof(ADJOINT_HEAT_TRANSFER).EquationId();
    }
}

template<class PrimalElement>
void AdjointHeatDiffusionElement<PrimalElement>::GetDofList(
    DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& r_geom = this->GetGeometry();
    const unsigned int num_nodes = r_geom.PointsNumber();

    if (rElementalDofList.size() != num_nodes)
    {
        rElementalDofList.resize(num_nodes);
    }

    for (unsigned int i = 0; i < num_nodes; i++)
    {
        rElementalDofList.push_back(r_geom[i].pGetDof(ADJOINT_HEAT_TRANSFER));
    }
}

template class AdjointHeatDiffusionElement<LaplacianElement>;

}