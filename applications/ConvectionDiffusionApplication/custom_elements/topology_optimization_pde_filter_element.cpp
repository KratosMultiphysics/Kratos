//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Gianmarco Boscolo
//

#include "topology_optimization_pde_filter_element.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"

#include "utilities/element_size_calculator.h"

// include Fluid Topology Optimization Data
#include "custom_utilities/topology_optimization_pde_filter_element_data.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TElementData >
TopologyOptimizationPdeFilterElement<TElementData>::TopologyOptimizationPdeFilterElement(IndexType NewId):
    Element(NewId)
{}

template< class TElementData >
TopologyOptimizationPdeFilterElement<TElementData>::TopologyOptimizationPdeFilterElement(IndexType NewId, const NodesArrayType& ThisNodes):
    Element(NewId,ThisNodes)
{}


template< class TElementData >
TopologyOptimizationPdeFilterElement<TElementData>::TopologyOptimizationPdeFilterElement(IndexType NewId, GeometryType::Pointer pGeometry):
    Element(NewId,pGeometry)
{}


template< class TElementData >
TopologyOptimizationPdeFilterElement<TElementData>::TopologyOptimizationPdeFilterElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
    Element(NewId,pGeometry,pProperties)
{}


template< class TElementData >
TopologyOptimizationPdeFilterElement<TElementData>::~TopologyOptimizationPdeFilterElement()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <class TElementData>
Element::Pointer TopologyOptimizationPdeFilterElement<TElementData>::Create(
    IndexType NewId,
    NodesArrayType const &ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<TopologyOptimizationPdeFilterElement>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template <class TElementData>
Element::Pointer TopologyOptimizationPdeFilterElement<TElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<TopologyOptimizationPdeFilterElement>(NewId, pGeom, pProperties);
}

template< class TElementData >
void TopologyOptimizationPdeFilterElement<TElementData>::Initialize(const ProcessInfo& rCurrentProcessInfo) {
    KRATOS_TRY;

    // DO NOTHING

    KRATOS_CATCH("");
}

template <class TElementData>
void TopologyOptimizationPdeFilterElement<TElementData>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                      VectorType& rRightHandSideVector,
                                                      const ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    if constexpr (TElementData::ElementManagesTimeIntegration) {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(
            gauss_weights, shape_functions, shape_derivatives);
        const unsigned int number_of_gauss_points = gauss_weights.size();

        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        for (unsigned int g = 0; g < number_of_gauss_points; g++) 
        {

            this->UpdateIntegrationPointData(
                data, g, gauss_weights[g],
                row(shape_functions, g),shape_derivatives[g]);

            this->AddTimeIntegratedSystem(
                data, rLeftHandSideMatrix, rRightHandSideVector);
        }
    }
}

template <class TElementData>
void TopologyOptimizationPdeFilterElement<TElementData>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                                       const ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);

    if constexpr (TElementData::ElementManagesTimeIntegration) {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(
            gauss_weights, shape_functions, shape_derivatives);
        const unsigned int number_of_gauss_points = gauss_weights.size();

        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        for (unsigned int g = 0; g < number_of_gauss_points; g++) 
        {
            this->UpdateIntegrationPointData(
                data, g, gauss_weights[g],
                row(shape_functions, g),shape_derivatives[g]);

            this->AddTimeIntegratedLHS(data, rLeftHandSideMatrix);
        }
    }
}

template <class TElementData>
void TopologyOptimizationPdeFilterElement<TElementData>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                        const ProcessInfo& rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);

    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    if constexpr (TElementData::ElementManagesTimeIntegration) 
    {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(
            gauss_weights, shape_functions, shape_derivatives);
        const unsigned int number_of_gauss_points = gauss_weights.size();

        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);
        
        for (unsigned int g = 0; g < number_of_gauss_points; g++) 
        {
            this->UpdateIntegrationPointData(
                data, g, gauss_weights[g],
                row(shape_functions, g),shape_derivatives[g]);

            this->AddTimeIntegratedRHS(data, rRightHandSideVector);
        }
    }
}

template <class TElementData>
void TopologyOptimizationPdeFilterElement<TElementData>::CalculateLocalVelocityContribution(
    MatrixType& rDampMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if( rDampMatrix.size1() != LocalSize )
        rDampMatrix.resize(LocalSize,LocalSize,false);

    if( rRightHandSideVector.size() != LocalSize )
        rRightHandSideVector.resize(LocalSize,false);

    noalias(rDampMatrix) = ZeroMatrix(LocalSize,LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    if (!TElementData::ElementManagesTimeIntegration) {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(
            gauss_weights, shape_functions, shape_derivatives);
        const unsigned int number_of_gauss_points = gauss_weights.size();

        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        // Iterate over integration points to evaluate local contribution
        for (unsigned int g = 0; g < number_of_gauss_points; g++) {
            const auto& r_dndx = shape_derivatives[g];
            this->UpdateIntegrationPointData(
                data, g, gauss_weights[g],
                row(shape_functions, g),r_dndx);

            this->AddVelocitySystem(data, rDampMatrix, rRightHandSideVector);
        }
    }
}

template <class TElementData>
void TopologyOptimizationPdeFilterElement<TElementData>::CalculateMassMatrix(MatrixType& rMassMatrix,
                                                     const ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if (rMassMatrix.size1() != LocalSize)
        rMassMatrix.resize(LocalSize, LocalSize, false);

    noalias(rMassMatrix) = ZeroMatrix(LocalSize, LocalSize);

    if (!TElementData::ElementManagesTimeIntegration) {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(
            gauss_weights, shape_functions, shape_derivatives);
        const unsigned int number_of_gauss_points = gauss_weights.size();

        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        // Iterate over integration points to evaluate local contribution
        for (unsigned int g = 0; g < number_of_gauss_points; g++) {
            this->UpdateIntegrationPointData(
                data, g, gauss_weights[g],
                row(shape_functions, g),shape_derivatives[g]);

            this->AddMassLHS(data, rMassMatrix);
        }
    }
}

template< class TElementData >
void TopologyOptimizationPdeFilterElement< TElementData >::EquationIdVector(EquationIdVectorType &rResult, const ProcessInfo &rCurrentProcessInfo) const
{
    const GeometryType& r_geometry = this->GetGeometry();

    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

        const unsigned int cpos = this->GetGeometry()[0].GetDofPosition(PDE_FILTER_RESULT);

        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            rResult[LocalIndex++] = r_geometry[i].GetDof(PDE_FILTER_RESULT, cpos).EquationId();
        }
}

template< class TElementData >
void TopologyOptimizationPdeFilterElement< TElementData >::GetDofList(DofsVectorType &rElementalDofList, const ProcessInfo &rCurrentProcessInfo) const
{
    const GeometryType& r_geometry = this->GetGeometry();

    if (rElementalDofList.size() != LocalSize)
         rElementalDofList.resize(LocalSize);

    const unsigned int cpos = this->GetGeometry()[0].GetDofPosition(PDE_FILTER_RESULT);

    unsigned int LocalIndex = 0;
    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(PDE_FILTER_RESULT,cpos);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////


template< class TElementData >
GeometryData::IntegrationMethod TopologyOptimizationPdeFilterElement<TElementData>::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_2;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template< class TElementData >
int TopologyOptimizationPdeFilterElement<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    // THE CHECK FOR THIS ELEMENT DOES NOT CHECK THE EXISTENCE OF A CONSITUTITVE LAW SINCE IT HAS BEEN IMPLEMENTED WITHOUT IT
    // IN RERALITY THIS ELEMENT SOLVES A CONVECTION-DIFFUSION-REACTION EQUATION

    // Generic geometry check
    int out = Element::Check(rCurrentProcessInfo);
    if (out != 0) {
        return out;
    }

    // Check variables used by TElementData
    out = TElementData::Check(*this, rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Something is wrong with the elemental data of Element "
        << this->Info() << std::endl;

    const GeometryType& r_geometry = this->GetGeometry();

    for(unsigned int i=0; i<NumNodes; ++i)
    {
        const Node& rNode = r_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PDE_FILTER_RESULT,rNode);
    }

    // If this is a 2D problem, check that nodes are in XY plane
    if ( Dim == 2)
    {
        for (unsigned int i=0; i<NumNodes; ++i) {
            if (r_geometry[i].Z() != 0.0)
                KRATOS_ERROR << "Node " << r_geometry[i].Id() << "has non-zero Z coordinate." << std::endl;
        }
    }

    return out;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void TopologyOptimizationPdeFilterElement<TElementData>::CalculateOnIntegrationPoints(
    Variable<array_1d<double, 3 > > const& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{}

template< class TElementData >
void TopologyOptimizationPdeFilterElement<TElementData>::CalculateOnIntegrationPoints(
    Variable<double> const& rVariable,
    std::vector<double>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{}

template <class TElementData>
void TopologyOptimizationPdeFilterElement<TElementData>::CalculateOnIntegrationPoints(
    Variable<array_1d<double, 6>> const& rVariable,
    std::vector<array_1d<double, 6>>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{}

template <class TElementData>
void TopologyOptimizationPdeFilterElement<TElementData>::CalculateOnIntegrationPoints(
    Variable<Vector> const& rVariable,
    std::vector<Vector>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{}

template <class TElementData>
void TopologyOptimizationPdeFilterElement<TElementData>::CalculateOnIntegrationPoints(
    Variable<Matrix> const& rVariable,
    std::vector<Matrix>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output
template <class TElementData>
void TopologyOptimizationPdeFilterElement<TElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info() << std::endl;
}

template <class TElementData>
std::string TopologyOptimizationPdeFilterElement<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "TopologyOptimizationPdeFilterElement" << Dim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template <class TElementData>
double TopologyOptimizationPdeFilterElement<TElementData>::GetAtCoordinate(
    const typename TElementData::NodalScalarData& rValues,
    const typename TElementData::ShapeFunctionsType& rN) const
{
    double result = 0.0;

    for (size_t i = 0; i < NumNodes; i++) {
        result += rN[i] * rValues[i];
    }

    return result;
}

template <class TElementData>
array_1d<double, 3> TopologyOptimizationPdeFilterElement<TElementData>::GetAtCoordinate(
    const typename TElementData::NodalVectorData& rValues,
    const typename TElementData::ShapeFunctionsType& rN) const
{
    array_1d<double, 3> result = ZeroVector(3);

    for (size_t i = 0; i < NumNodes; i++) {
        for (size_t j = 0; j < Dim; j++) {
            result[j] += rN[i] * rValues(i, j);
        }
    }

    return result;
}

template <class TElementData>
BoundedMatrix<double, TElementData::Dim, TElementData::Dim> TopologyOptimizationPdeFilterElement<TElementData>::GetAtCoordinate(
    const typename TElementData::NodalTensorData &rValues,
    const typename TElementData::ShapeFunctionsType &rN) const
{
    BoundedMatrix<double,Dim,Dim> result = ZeroMatrix(Dim,Dim);

    for (size_t i = 0; i < NumNodes; i++) {
        noalias(result) += rN[i] * rValues[i];
    }

    return result;
}

template <class TElementData>
double TopologyOptimizationPdeFilterElement<TElementData>::GetAtCoordinate(
    const double Value,
    const typename TElementData::ShapeFunctionsType& rN) const
{
    return Value;
}

template <class TElementData>
void TopologyOptimizationPdeFilterElement<TElementData>::UpdateIntegrationPointData(
    TElementData& rData,
    unsigned int IntegrationPointIndex,
    double Weight,
    const typename TElementData::MatrixRowType& rN,
    const typename TElementData::ShapeDerivativesType& rDN_DX) const
{
    rData.UpdateGeometryValues(IntegrationPointIndex, Weight, rN, rDN_DX);
}

template< class TElementData >
void TopologyOptimizationPdeFilterElement<TElementData>::CalculateGeometryData(Vector &rGaussWeights,
                                      Matrix &rNContainer,
                                      ShapeFunctionDerivativesArrayType &rDN_DX) const
{
    const GeometryData::IntegrationMethod integration_method = this->GetIntegrationMethod();
    const GeometryType& r_geometry = this->GetGeometry();
    const unsigned int number_of_gauss_points = r_geometry.IntegrationPointsNumber(integration_method);

    Vector DetJ;
    r_geometry.ShapeFunctionsIntegrationPointsGradients(rDN_DX,DetJ,integration_method);

    if (rNContainer.size1() != number_of_gauss_points || rNContainer.size2() != NumNodes) {
        rNContainer.resize(number_of_gauss_points,NumNodes,false);
    }
    rNContainer = r_geometry.ShapeFunctionsValues(integration_method);

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = r_geometry.IntegrationPoints(integration_method);

    if (rGaussWeights.size() != number_of_gauss_points) {
        rGaussWeights.resize(number_of_gauss_points,false);
    }

    for (unsigned int g = 0; g < number_of_gauss_points; g++)
        rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
}

template <class TElementData>
void TopologyOptimizationPdeFilterElement<TElementData>::Calculate(
    const Variable<double> &rVariable,
    double &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    Element::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TElementData>
void TopologyOptimizationPdeFilterElement<TElementData>::Calculate(
    const Variable<array_1d<double, 3>> &rVariable,
    array_1d<double, 3> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    Element::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TElementData>
void TopologyOptimizationPdeFilterElement<TElementData>::Calculate(
    const Variable<Vector >& rVariable,
    Vector& rOutput,
    const ProcessInfo& rCurrentProcessInfo )
{
    Element::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TElementData>
void TopologyOptimizationPdeFilterElement<TElementData>::Calculate(
    const Variable<Matrix> &rVariable,
    Matrix &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    Element::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void TopologyOptimizationPdeFilterElement<TElementData>::ConvectionOperator(Vector &rResult,
                                   const array_1d<double,3> &rConvVel,
                                   const ShapeFunctionDerivativesType &DN_DX) const
{
    if(rResult.size() != NumNodes) rResult.resize(NumNodes,false);

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        rResult[i] = rConvVel[0]*DN_DX(i,0);
        for(unsigned int k = 1; k < Dim; k++)
            rResult[i] += rConvVel[k]*DN_DX(i,k);
    }
}

template <class TElementData>
void TopologyOptimizationPdeFilterElement<TElementData>::AddTimeIntegratedSystem(
    TElementData& rData, MatrixType& rLHS, VectorType& rRHS) 
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void TopologyOptimizationPdeFilterElement<TElementData>::AddTimeIntegratedLHS(
    TElementData& rData, MatrixType& rLHS) 
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
}

template <class TElementData>
void TopologyOptimizationPdeFilterElement<TElementData>::AddTimeIntegratedRHS(
    TElementData& rData, VectorType& rRHS) 
{
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

// BUILD NS SYSTEM
template <>
void TopologyOptimizationPdeFilterElement< TopologyOptimizationPdeFilterElementData<2,3,true>>::ComputeGaussPointLHSContribution(
    TopologyOptimizationPdeFilterElementData<2,3,true> & rData,
    MatrixType& rLHS)
{
    const array_1d<double,3> D = rData.PdeFilterDiffusion;
    const array_1d<double,3> r = rData.PdeFilterReaction;

    const double h = rData.ElementSize;
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;

    // TRANSPORT ELEMENTAL LHS MATRIX
    const double crLHS0 = pow(N[0], 2);
const double crLHS1 = N[0]*r[0] + N[1]*r[1] + N[2]*r[2];
const double crLHS2 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2];
const double crLHS3 = pow(crLHS1, 2);
const double crLHS4 = 1.0/(crLHS1*stab_c3 + crLHS2*stab_c1/pow(h, 2));
const double crLHS5 = N[0]*crLHS1;
const double crLHS6 = -gauss_weight*(1.0*N[0]*N[1]*crLHS3*crLHS4 - N[1]*crLHS5 - crLHS2*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1)));
const double crLHS7 = -gauss_weight*(1.0*N[0]*N[2]*crLHS3*crLHS4 - N[2]*crLHS5 - crLHS2*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1)));
const double crLHS8 = pow(N[1], 2);
const double crLHS9 = -gauss_weight*(-N[1]*N[2]*crLHS1 + 1.0*N[1]*N[2]*crLHS3*crLHS4 - crLHS2*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1)));
const double crLHS10 = pow(N[2], 2);
rLHS(0,0)+=-gauss_weight*(-crLHS0*crLHS1 + 1.0*crLHS0*crLHS3*crLHS4 - crLHS2*(pow(DN(0,0), 2) + pow(DN(0,1), 2)));
rLHS(0,1)+=crLHS6;
rLHS(0,2)+=crLHS7;
rLHS(1,0)+=crLHS6;
rLHS(1,1)+=-gauss_weight*(-crLHS1*crLHS8 - crLHS2*(pow(DN(1,0), 2) + pow(DN(1,1), 2)) + 1.0*crLHS3*crLHS4*crLHS8);
rLHS(1,2)+=crLHS9;
rLHS(2,0)+=crLHS7;
rLHS(2,1)+=crLHS9;
rLHS(2,2)+=-gauss_weight*(-crLHS1*crLHS10 + 1.0*crLHS10*crLHS3*crLHS4 - crLHS2*(pow(DN(2,0), 2) + pow(DN(2,1), 2)));
 
    
}

template <>
void TopologyOptimizationPdeFilterElement< TopologyOptimizationPdeFilterElementData<2,4,true>>::ComputeGaussPointLHSContribution(
    TopologyOptimizationPdeFilterElementData<2,4,true> & rData,
    MatrixType& rLHS)
{
    const array_1d<double,4> D = rData.PdeFilterDiffusion;
    const array_1d<double,4> r = rData.PdeFilterReaction;

    const double h = rData.ElementSize;
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,2>& DN = rData.DN_DX;

    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;

    // TRANSPORT ELEMENTAL LHS MATRIX
    const double crLHS0 = pow(N[0], 2);
const double crLHS1 = N[0]*r[0] + N[1]*r[1] + N[2]*r[2] + N[3]*r[3];
const double crLHS2 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2] + D[3]*N[3];
const double crLHS3 = pow(crLHS1, 2);
const double crLHS4 = 1.0/(crLHS1*stab_c3 + crLHS2*stab_c1/pow(h, 2));
const double crLHS5 = N[0]*crLHS1;
const double crLHS6 = -gauss_weight*(1.0*N[0]*N[1]*crLHS3*crLHS4 - N[1]*crLHS5 - crLHS2*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1)));
const double crLHS7 = -gauss_weight*(1.0*N[0]*N[2]*crLHS3*crLHS4 - N[2]*crLHS5 - crLHS2*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1)));
const double crLHS8 = -gauss_weight*(1.0*N[0]*N[3]*crLHS3*crLHS4 - N[3]*crLHS5 - crLHS2*(DN(0,0)*DN(3,0) + DN(0,1)*DN(3,1)));
const double crLHS9 = pow(N[1], 2);
const double crLHS10 = N[1]*crLHS1;
const double crLHS11 = -gauss_weight*(1.0*N[1]*N[2]*crLHS3*crLHS4 - N[2]*crLHS10 - crLHS2*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1)));
const double crLHS12 = -gauss_weight*(1.0*N[1]*N[3]*crLHS3*crLHS4 - N[3]*crLHS10 - crLHS2*(DN(1,0)*DN(3,0) + DN(1,1)*DN(3,1)));
const double crLHS13 = pow(N[2], 2);
const double crLHS14 = -gauss_weight*(-N[2]*N[3]*crLHS1 + 1.0*N[2]*N[3]*crLHS3*crLHS4 - crLHS2*(DN(2,0)*DN(3,0) + DN(2,1)*DN(3,1)));
const double crLHS15 = pow(N[3], 2);
rLHS(0,0)+=-gauss_weight*(-crLHS0*crLHS1 + 1.0*crLHS0*crLHS3*crLHS4 - crLHS2*(pow(DN(0,0), 2) + pow(DN(0,1), 2)));
rLHS(0,1)+=crLHS6;
rLHS(0,2)+=crLHS7;
rLHS(0,3)+=crLHS8;
rLHS(1,0)+=crLHS6;
rLHS(1,1)+=-gauss_weight*(-crLHS1*crLHS9 - crLHS2*(pow(DN(1,0), 2) + pow(DN(1,1), 2)) + 1.0*crLHS3*crLHS4*crLHS9);
rLHS(1,2)+=crLHS11;
rLHS(1,3)+=crLHS12;
rLHS(2,0)+=crLHS7;
rLHS(2,1)+=crLHS11;
rLHS(2,2)+=-gauss_weight*(-crLHS1*crLHS13 + 1.0*crLHS13*crLHS3*crLHS4 - crLHS2*(pow(DN(2,0), 2) + pow(DN(2,1), 2)));
rLHS(2,3)+=crLHS14;
rLHS(3,0)+=crLHS8;
rLHS(3,1)+=crLHS12;
rLHS(3,2)+=crLHS14;
rLHS(3,3)+=-gauss_weight*(-crLHS1*crLHS15 + 1.0*crLHS15*crLHS3*crLHS4 - crLHS2*(pow(DN(3,0), 2) + pow(DN(3,1), 2)));
 
    
}

template <>
void TopologyOptimizationPdeFilterElement<TopologyOptimizationPdeFilterElementData<3,4,true>>::ComputeGaussPointLHSContribution(
    TopologyOptimizationPdeFilterElementData<3,4,true>& rData,
    MatrixType& rLHS)
{
    const array_1d<double,4> D = rData.PdeFilterDiffusion;
    const array_1d<double,4> r = rData.PdeFilterReaction;

    const double h = rData.ElementSize;
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,3>& DN = rData.DN_DX;

    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;

    // TRANSPORT ELEMENTAL LHS MATRIX
    const double crLHS0 = pow(N[0], 2);
const double crLHS1 = N[0]*r[0] + N[1]*r[1] + N[2]*r[2] + N[3]*r[3];
const double crLHS2 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2] + D[3]*N[3];
const double crLHS3 = pow(crLHS1, 2);
const double crLHS4 = 1.0/(crLHS1*stab_c3 + crLHS2*stab_c1/pow(h, 2));
const double crLHS5 = N[0]*crLHS1;
const double crLHS6 = -gauss_weight*(1.0*N[0]*N[1]*crLHS3*crLHS4 - N[1]*crLHS5 - crLHS2*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1) + DN(0,2)*DN(1,2)));
const double crLHS7 = -gauss_weight*(1.0*N[0]*N[2]*crLHS3*crLHS4 - N[2]*crLHS5 - crLHS2*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1) + DN(0,2)*DN(2,2)));
const double crLHS8 = -gauss_weight*(1.0*N[0]*N[3]*crLHS3*crLHS4 - N[3]*crLHS5 - crLHS2*(DN(0,0)*DN(3,0) + DN(0,1)*DN(3,1) + DN(0,2)*DN(3,2)));
const double crLHS9 = pow(N[1], 2);
const double crLHS10 = N[1]*crLHS1;
const double crLHS11 = -gauss_weight*(1.0*N[1]*N[2]*crLHS3*crLHS4 - N[2]*crLHS10 - crLHS2*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1) + DN(1,2)*DN(2,2)));
const double crLHS12 = -gauss_weight*(1.0*N[1]*N[3]*crLHS3*crLHS4 - N[3]*crLHS10 - crLHS2*(DN(1,0)*DN(3,0) + DN(1,1)*DN(3,1) + DN(1,2)*DN(3,2)));
const double crLHS13 = pow(N[2], 2);
const double crLHS14 = -gauss_weight*(-N[2]*N[3]*crLHS1 + 1.0*N[2]*N[3]*crLHS3*crLHS4 - crLHS2*(DN(2,0)*DN(3,0) + DN(2,1)*DN(3,1) + DN(2,2)*DN(3,2)));
const double crLHS15 = pow(N[3], 2);
rLHS(0,0)+=-gauss_weight*(-crLHS0*crLHS1 + 1.0*crLHS0*crLHS3*crLHS4 - crLHS2*(pow(DN(0,0), 2) + pow(DN(0,1), 2) + pow(DN(0,2), 2)));
rLHS(0,1)+=crLHS6;
rLHS(0,2)+=crLHS7;
rLHS(0,3)+=crLHS8;
rLHS(1,0)+=crLHS6;
rLHS(1,1)+=-gauss_weight*(-crLHS1*crLHS9 - crLHS2*(pow(DN(1,0), 2) + pow(DN(1,1), 2) + pow(DN(1,2), 2)) + 1.0*crLHS3*crLHS4*crLHS9);
rLHS(1,2)+=crLHS11;
rLHS(1,3)+=crLHS12;
rLHS(2,0)+=crLHS7;
rLHS(2,1)+=crLHS11;
rLHS(2,2)+=-gauss_weight*(-crLHS1*crLHS13 + 1.0*crLHS13*crLHS3*crLHS4 - crLHS2*(pow(DN(2,0), 2) + pow(DN(2,1), 2) + pow(DN(2,2), 2)));
rLHS(2,3)+=crLHS14;
rLHS(3,0)+=crLHS8;
rLHS(3,1)+=crLHS12;
rLHS(3,2)+=crLHS14;
rLHS(3,3)+=-gauss_weight*(-crLHS1*crLHS15 + 1.0*crLHS15*crLHS3*crLHS4 - crLHS2*(pow(DN(3,0), 2) + pow(DN(3,1), 2) + pow(DN(3,2), 2)));

}

template <>
void TopologyOptimizationPdeFilterElement<TopologyOptimizationPdeFilterElementData<2,3,true>>::ComputeGaussPointRHSContribution(
    TopologyOptimizationPdeFilterElementData<2,3,true>& rData,
    VectorType& rRHS)
{
    const array_1d<double,3> D = rData.PdeFilterDiffusion;
    const array_1d<double,3> r = rData.PdeFilterReaction;

    const double h = rData.ElementSize;
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const array_1d<double,3>& filter = rData.PdeFilterResult;
    const array_1d<double,3>& source = rData.PdeFilterForcing;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // TRANSPORT ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*source[0] + N[1]*source[1] + N[2]*source[2];
const double crRHS1 = N[0]*r[0] + N[1]*r[1] + N[2]*r[2];
const double crRHS2 = crRHS1*(N[0]*filter[0] + N[1]*filter[1] + N[2]*filter[2]);
const double crRHS3 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2];
const double crRHS4 = DN(0,0)*filter[0] + DN(1,0)*filter[1] + DN(2,0)*filter[2];
const double crRHS5 = DN(0,1)*filter[0] + DN(1,1)*filter[1] + DN(2,1)*filter[2];
const double crRHS6 = 1.0*crRHS1*(crRHS0 - crRHS2)/(crRHS1*stab_c3 + crRHS3*stab_c1/pow(h, 2));
rRHS[0]+=-gauss_weight*(-N[0]*crRHS0 + N[0]*crRHS2 + N[0]*crRHS6 + crRHS3*(DN(0,0)*crRHS4 + DN(0,1)*crRHS5));
rRHS[1]+=-gauss_weight*(-N[1]*crRHS0 + N[1]*crRHS2 + N[1]*crRHS6 + crRHS3*(DN(1,0)*crRHS4 + DN(1,1)*crRHS5));
rRHS[2]+=-gauss_weight*(-N[2]*crRHS0 + N[2]*crRHS2 + N[2]*crRHS6 + crRHS3*(DN(2,0)*crRHS4 + DN(2,1)*crRHS5));

}

template <>
void TopologyOptimizationPdeFilterElement<TopologyOptimizationPdeFilterElementData<2,4,true>>::ComputeGaussPointRHSContribution(
    TopologyOptimizationPdeFilterElementData<2,4,true>& rData,
    VectorType& rRHS)
{
    const array_1d<double,4> D = rData.PdeFilterDiffusion;
    const array_1d<double,4> r = rData.PdeFilterReaction;

    const double h = rData.ElementSize;
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,2>& DN = rData.DN_DX;

    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const array_1d<double,4>& filter = rData.PdeFilterResult;
    const array_1d<double,4>& source = rData.PdeFilterForcing;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // TRANSPORT ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*source[0] + N[1]*source[1] + N[2]*source[2] + N[3]*source[3];
const double crRHS1 = N[0]*r[0] + N[1]*r[1] + N[2]*r[2] + N[3]*r[3];
const double crRHS2 = crRHS1*(N[0]*filter[0] + N[1]*filter[1] + N[2]*filter[2] + N[3]*filter[3]);
const double crRHS3 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2] + D[3]*N[3];
const double crRHS4 = DN(0,0)*filter[0] + DN(1,0)*filter[1] + DN(2,0)*filter[2] + DN(3,0)*filter[3];
const double crRHS5 = DN(0,1)*filter[0] + DN(1,1)*filter[1] + DN(2,1)*filter[2] + DN(3,1)*filter[3];
const double crRHS6 = 1.0*crRHS1*(crRHS0 - crRHS2)/(crRHS1*stab_c3 + crRHS3*stab_c1/pow(h, 2));
rRHS[0]+=-gauss_weight*(-N[0]*crRHS0 + N[0]*crRHS2 + N[0]*crRHS6 + crRHS3*(DN(0,0)*crRHS4 + DN(0,1)*crRHS5));
rRHS[1]+=-gauss_weight*(-N[1]*crRHS0 + N[1]*crRHS2 + N[1]*crRHS6 + crRHS3*(DN(1,0)*crRHS4 + DN(1,1)*crRHS5));
rRHS[2]+=-gauss_weight*(-N[2]*crRHS0 + N[2]*crRHS2 + N[2]*crRHS6 + crRHS3*(DN(2,0)*crRHS4 + DN(2,1)*crRHS5));
rRHS[3]+=-gauss_weight*(-N[3]*crRHS0 + N[3]*crRHS2 + N[3]*crRHS6 + crRHS3*(DN(3,0)*crRHS4 + DN(3,1)*crRHS5));

}

template <>
void TopologyOptimizationPdeFilterElement<TopologyOptimizationPdeFilterElementData<3,4,true>>::ComputeGaussPointRHSContribution(
    TopologyOptimizationPdeFilterElementData<3,4,true>& rData,
    VectorType& rRHS)
{
    const array_1d<double,4> D = rData.PdeFilterDiffusion;
    const array_1d<double,4> r = rData.PdeFilterReaction;

    const double h = rData.ElementSize;
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,3>& DN = rData.DN_DX;

    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const array_1d<double,4>& filter = rData.PdeFilterResult;
    const array_1d<double,4>& source = rData.PdeFilterForcing;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // TRANSPORT ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*source[0] + N[1]*source[1] + N[2]*source[2] + N[3]*source[3];
const double crRHS1 = N[0]*r[0] + N[1]*r[1] + N[2]*r[2] + N[3]*r[3];
const double crRHS2 = crRHS1*(N[0]*filter[0] + N[1]*filter[1] + N[2]*filter[2] + N[3]*filter[3]);
const double crRHS3 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2] + D[3]*N[3];
const double crRHS4 = DN(0,0)*filter[0] + DN(1,0)*filter[1] + DN(2,0)*filter[2] + DN(3,0)*filter[3];
const double crRHS5 = DN(0,1)*filter[0] + DN(1,1)*filter[1] + DN(2,1)*filter[2] + DN(3,1)*filter[3];
const double crRHS6 = DN(0,2)*filter[0] + DN(1,2)*filter[1] + DN(2,2)*filter[2] + DN(3,2)*filter[3];
const double crRHS7 = 1.0*crRHS1*(crRHS0 - crRHS2)/(crRHS1*stab_c3 + crRHS3*stab_c1/pow(h, 2));
rRHS[0]+=-gauss_weight*(-N[0]*crRHS0 + N[0]*crRHS2 + N[0]*crRHS7 + crRHS3*(DN(0,0)*crRHS4 + DN(0,1)*crRHS5 + DN(0,2)*crRHS6));
rRHS[1]+=-gauss_weight*(-N[1]*crRHS0 + N[1]*crRHS2 + N[1]*crRHS7 + crRHS3*(DN(1,0)*crRHS4 + DN(1,1)*crRHS5 + DN(1,2)*crRHS6));
rRHS[2]+=-gauss_weight*(-N[2]*crRHS0 + N[2]*crRHS2 + N[2]*crRHS7 + crRHS3*(DN(2,0)*crRHS4 + DN(2,1)*crRHS5 + DN(2,2)*crRHS6));
rRHS[3]+=-gauss_weight*(-N[3]*crRHS0 + N[3]*crRHS2 + N[3]*crRHS7 + crRHS3*(DN(3,0)*crRHS4 + DN(3,1)*crRHS5 + DN(3,2)*crRHS6));

}

template <class TElementData>
void TopologyOptimizationPdeFilterElement<TElementData>::AddVelocitySystem(
    TElementData& rData,
    MatrixType& rLocalLHS,
    VectorType& rLocalRHS) 
    {
    KRATOS_TRY;

    KRATOS_ERROR << "Calling base TopologyOptimizationPdeFilterElement::AddVelocitySystem "
                    "implementation. This method is not supported by your "
                    "element."
                 << std::endl;

    KRATOS_CATCH("");
}

template <class TElementData>
void TopologyOptimizationPdeFilterElement<TElementData>::AddMassLHS(
    TElementData& rData, MatrixType& rMassMatrix) {
    KRATOS_TRY;

    KRATOS_ERROR << "Calling base TopologyOptimizationPdeFilterElement::AddMassLHS "
                    "implementation. This method is not supported by your "
                    "element."
                 << std::endl;

    KRATOS_CATCH("");
}

template <class TElementData>
void TopologyOptimizationPdeFilterElement<TElementData>::GetCurrentValuesVector(
    const TElementData& rData,
    array_1d<double,LocalSize>& rValues) const {

    int local_index = 0;

    const auto& r_filters = rData.PdeFilterResult;

    for (unsigned int i = 0; i < NumNodes; ++i) {
        rValues[local_index++] = r_filters[i];  // PdeFilterResult Dof
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template< class TElementData >
void TopologyOptimizationPdeFilterElement<TElementData>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    rSerializer.save("mpConstitutiveLaw",this->mpConstitutiveLaw);
}


template< class TElementData >
void TopologyOptimizationPdeFilterElement<TElementData>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    rSerializer.load("mpConstitutiveLaw",this->mpConstitutiveLaw);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class TopologyOptimizationPdeFilterElement< TopologyOptimizationPdeFilterElementData<2,3,true> >;
template class TopologyOptimizationPdeFilterElement< TopologyOptimizationPdeFilterElementData<2,4,true> >;
template class TopologyOptimizationPdeFilterElement< TopologyOptimizationPdeFilterElementData<3,4,true> >;

///////////////////////////////////////////////////////////////////////////////////////////////////
}