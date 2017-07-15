//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#include "flux_condition.h"
#include "includes/convection_diffusion_settings.h"

namespace Kratos
{

// Public Life Cycle //////////////////////////////////////////////////////////

template< unsigned int TNodeNumber >
FluxCondition<TNodeNumber>::FluxCondition(IndexType NewId, Geometry< Node<3> >::Pointer pGeometry):
    Condition(NewId,pGeometry)
{
}

template< unsigned int TNodeNumber >
FluxCondition<TNodeNumber>::FluxCondition(
    IndexType NewId,
    Geometry< Node<3> >::Pointer pGeometry,
    Properties::Pointer pProperties):
    Condition(NewId,pGeometry,pProperties)
{
}

template< unsigned int TNodeNumber >
FluxCondition<TNodeNumber>::~FluxCondition()
{
}

// Public Operations //////////////////////////////////////////////////////////

template< unsigned int TNodeNumber >
Condition::Pointer FluxCondition<TNodeNumber>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Condition::Pointer(new FluxCondition(NewId, GetGeometry().Create(ThisNodes), pProperties) );
}

template< unsigned int TNodeNumber >
void FluxCondition<TNodeNumber>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rLeftHandSideMatrix.size1() != TNodeNumber)
    {
        rLeftHandSideMatrix.resize(TNodeNumber,TNodeNumber,false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(TNodeNumber,TNodeNumber);

    this->CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);

    KRATOS_CATCH("")
}

template< unsigned int TNodeNumber >
void FluxCondition<TNodeNumber>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rRightHandSideVector.size() != TNodeNumber)
    {
        rRightHandSideVector.resize(TNodeNumber,false);
    }
    noalias(rRightHandSideVector) = ZeroVector(TNodeNumber);

    const ProcessInfo& rProcessInfo = rCurrentProcessInfo; // to ensure that Gets are threadsafe
    ConvectionDiffusionSettings& rSettings = *(rProcessInfo[CONVECTION_DIFFUSION_SETTINGS]);

    const Variable<double>& rFluxVar = rSettings.GetSurfaceSourceVariable();

    FluxConditionInternals::IntegrationData<TNodeNumber> Data(this->GetGeometry(), rFluxVar);

    for (unsigned int g = 0; g < Data.NumGauss; g++)
    {
        Data.SetCurrentGaussPoint(g);

        this->AddIntegrationPointRHSContribution(rRightHandSideVector, Data);
    }

    KRATOS_CATCH("")
}

template< unsigned int TNodeNumber >
void FluxCondition<TNodeNumber>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const ProcessInfo& rProcessInfo = rCurrentProcessInfo; // to ensure that Gets are threadsafe
    ConvectionDiffusionSettings& rSettings = *(rProcessInfo[CONVECTION_DIFFUSION_SETTINGS]);

    const Variable<double>& rUnknownVar = rSettings.GetUnknownVariable();
    
    if (rResult.size() != TNodeNumber)
    {
        rResult.resize(TNodeNumber,false);
    }

    Geometry< Node<3> >& rGeometry = this->GetGeometry();

    for (unsigned int i = 0; i < TNodeNumber; i++)
    {
        rResult[i] = rGeometry[i].GetDof(rUnknownVar).EquationId();
    }

    KRATOS_CATCH("")
}

template< unsigned int TNodeNumber >
void FluxCondition<TNodeNumber>::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const ProcessInfo& rProcessInfo = rCurrentProcessInfo; // to ensure that Gets are threadsafe
    ConvectionDiffusionSettings& rSettings = *(rProcessInfo[CONVECTION_DIFFUSION_SETTINGS]);

    const Variable<double>& rUnknownVar = rSettings.GetUnknownVariable();

    if (rConditionalDofList.size() != TNodeNumber)
    {
        rConditionalDofList.resize(TNodeNumber);
    }

    Geometry< Node<3> >& rGeometry = this->GetGeometry();

    for (unsigned int i = 0; i < TNodeNumber; i++)
    {
        rConditionalDofList[i] = rGeometry[i].pGetDof(rUnknownVar);
    }

    KRATOS_CATCH("")
}

template< unsigned int TNodeNumber >
GeometryData::IntegrationMethod FluxCondition<TNodeNumber>::GetIntegrationMethod()
{
    return GeometryData::GI_GAUSS_2;
}

// Input and Output ///////////////////////////////////////////////////////////

template< unsigned int TNodeNumber >
std::string FluxCondition<TNodeNumber>::Info() const
{
    std::stringstream buffer;
    buffer << "FluxCondition #" << Id();
    return buffer.str();
}

template< unsigned int TNodeNumber >
void FluxCondition<TNodeNumber>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "FluxCondition #" << Id();
}

template< unsigned int TNodeNumber >
void FluxCondition<TNodeNumber>::PrintData(std::ostream& rOStream) const
{
    rOStream << "FluxCondition #" << Id() << std::endl;
    this->GetGeometry().PrintData(rOStream);
}


// Finite element functions ///////////////////////////////////////////////////

template< unsigned int TNodeNumber >
void FluxCondition<TNodeNumber>::AddIntegrationPointRHSContribution(
    VectorType& rRightHandSideVector,
    const FluxConditionInternals::IntegrationData<TNodeNumber>& rData)
{
    double InterpolatedFlux = rData.GaussPointFlux();
    for (unsigned int i = 0; i < TNodeNumber; i++)
    {
        rRightHandSideVector[i] += rData.N(i) * InterpolatedFlux * rData.IntegrationWeight();
    }
}


// Serialization //////////////////////////////////////////////////////////////

template< unsigned int TNodeNumber >
FluxCondition<TNodeNumber>::FluxCondition():
    Condition()
{
}

template< unsigned int TNodeNumber >
void FluxCondition<TNodeNumber>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

template< unsigned int TNodeNumber >
void FluxCondition<TNodeNumber>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

template class FluxCondition<2>;
template class FluxCondition<3>;

}