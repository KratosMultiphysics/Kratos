// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___ 
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _| 
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//


#include "flux_condition.h"
#include "includes/convection_diffusion_settings.h"
#include "includes/variables.h"

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
    return Kratos::make_shared<FluxCondition<TNodeNumber>>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

template< unsigned int TNodeNumber >
Condition::Pointer FluxCondition<TNodeNumber>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_shared<FluxCondition<TNodeNumber>>(NewId, pGeom, pProperties);
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

template< unsigned int TNodeNumber >
void FluxCondition<TNodeNumber>::GetValueOnIntegrationPoints(
    const Variable<array_1d<double,3> > &rVariable,
    std::vector<array_1d<double,3> > &rValues,
    const ProcessInfo &rCurrentProcessInfo)
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(NumGauss);
    if (rVariable == NORMAL)
    {
        this->CalculateNormal(rValues[0]);
    }
    else
    {
        /* The cast is done to avoid modification of the element's data. Data modification
         * would happen if rVariable is not stored now (would initialize a pointer to &rVariable
         * with associated value of 0.0). This is catastrophic if the variable referenced
         * goes out of scope.
         */
        const FluxCondition* const_this = static_cast< const FluxCondition* >(this);
        rValues[0] = const_this->GetValue(rVariable);
    }

    // Copy the values to the different gauss points
    for (unsigned int g = 1; g < NumGauss; g++)
    {
        noalias(rValues[g]) = rValues[0];
    }
}

template< unsigned int TNodeNumber >
void FluxCondition<TNodeNumber>::GetValueOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(NumGauss);
    const FluxCondition* const_this = static_cast< const FluxCondition* >(this);
    rValues[0] = const_this->GetValue(rVariable);
    for (unsigned int g = 1; g < NumGauss; g++)
    {
        rValues[g] = rValues[0];
    }
}

template< unsigned int TNodeNumber >
void FluxCondition<TNodeNumber>::GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 6 > >& rVariable,
    std::vector<array_1d<double, 6 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(NumGauss);
    const FluxCondition* const_this = static_cast< const FluxCondition* >(this);
    rValues[0] = const_this->GetValue(rVariable);
    for (unsigned int g = 1; g < NumGauss; g++)
    {
        noalias(rValues[g]) = rValues[0];
    }
}

template< unsigned int TNodeNumber >
void FluxCondition<TNodeNumber>::GetValueOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(NumGauss);
    const FluxCondition* const_this = static_cast< const FluxCondition* >(this);
    rValues[0] = const_this->GetValue(rVariable);
    for (unsigned int g = 1; g < NumGauss; g++)
    {
        noalias(rValues[g]) = rValues[0];
    }
}

template< unsigned int TNodeNumber >
void FluxCondition<TNodeNumber>::GetValueOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(NumGauss);
    const FluxCondition* const_this = static_cast< const FluxCondition* >(this);
    rValues[0] = const_this->GetValue(rVariable);
    for (unsigned int g = 1; g < NumGauss; g++)
    {
        noalias(rValues[g]) = rValues[0];
    }
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


template <>
void FluxCondition<2>::CalculateNormal(array_1d<double,3>& An)
{
    Geometry<Node<3> >& pGeometry = this->GetGeometry();

    An[0] =   pGeometry[1].Y() - pGeometry[0].Y();
    An[1] = - (pGeometry[1].X() - pGeometry[0].X());
    An[2] =    0.00;
}


template <>
void FluxCondition<3>::CalculateNormal(array_1d<double,3>& An )
{
    Geometry<Node<3> >& pGeometry = this->GetGeometry();

    array_1d<double,3> v1,v2;
    v1[0] = pGeometry[1].X() - pGeometry[0].X();
    v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
    v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

    v2[0] = pGeometry[2].X() - pGeometry[0].X();
    v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
    v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

    MathUtils<double>::CrossProduct(An,v1,v2);
    An *= 0.5;
}

template <>
void FluxCondition<4>::CalculateNormal(array_1d<double,3>& An )
{
    KRATOS_ERROR << "This function is not yet implemented" << std::endl;
  
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
template class FluxCondition<4>;

}
