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

#include "transport_topology_optimization_element.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"

#include "utilities/element_size_calculator.h"

// include Fluid Topology Optimization Data
#include "custom_utilities/transport_topology_optimization_element_data.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TElementData >
TransportTopologyOptimizationElement<TElementData>::TransportTopologyOptimizationElement(IndexType NewId):
    Element(NewId)
{}

template< class TElementData >
TransportTopologyOptimizationElement<TElementData>::TransportTopologyOptimizationElement(IndexType NewId, const NodesArrayType& ThisNodes):
    Element(NewId,ThisNodes)
{}


template< class TElementData >
TransportTopologyOptimizationElement<TElementData>::TransportTopologyOptimizationElement(IndexType NewId, GeometryType::Pointer pGeometry):
    Element(NewId,pGeometry)
{}


template< class TElementData >
TransportTopologyOptimizationElement<TElementData>::TransportTopologyOptimizationElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
    Element(NewId,pGeometry,pProperties)
{}


template< class TElementData >
TransportTopologyOptimizationElement<TElementData>::~TransportTopologyOptimizationElement()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <class TElementData>
Element::Pointer TransportTopologyOptimizationElement<TElementData>::Create(
    IndexType NewId,
    NodesArrayType const &ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<TransportTopologyOptimizationElement>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template <class TElementData>
Element::Pointer TransportTopologyOptimizationElement<TElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<TransportTopologyOptimizationElement>(NewId, pGeom, pProperties);
}

template< class TElementData >
void TransportTopologyOptimizationElement<TElementData>::Initialize(const ProcessInfo& rCurrentProcessInfo) {
    KRATOS_TRY;

    // DO NOTHING

    KRATOS_CATCH("");
}

template <class TElementData>
void TransportTopologyOptimizationElement<TElementData>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
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

        // Now the distinction of what to do for the primal or adjoint problem must be done calling the relative Kratos property
        unsigned int problem_physics = data.TopOptProblemStage;
        // problem_physics = 1: NS equations
        // problem_physics = 2: ADJOINT NS equations
        if (problem_physics == 1)
        {
            // PRIMAL TRANSPORT
            // Iterate over integration points to evaluate local contribution
            for (unsigned int g = 0; g < number_of_gauss_points; g++) 
            {

                this->UpdateIntegrationPointData(
                    data, g, gauss_weights[g],
                    row(shape_functions, g),shape_derivatives[g]);

                this->AddTimeIntegratedSystem(
                    data, rLeftHandSideMatrix, rRightHandSideVector);
            }
        }
        else if (problem_physics == 2)
        {
            // ADJOINT TRANSPORT
            // Iterate over integration points to evaluate local contribution
            for (unsigned int g = 0; g < number_of_gauss_points; g++) {

                this->UpdateIntegrationPointData(
                    data, g, gauss_weights[g],
                    row(shape_functions, g),shape_derivatives[g]);

                this->AddTimeIntegratedSystemAdjoint(
                    data, rLeftHandSideMatrix, rRightHandSideVector);
            }
        }
        else
        {
            if (problem_physics != 0)
            {
                KRATOS_ERROR << "\nInvalid value for the variable TRANSPORT_TOP_OPT_PROBLEM_STAGE for the Fluid Topology Optimization aplication. |\t TopOptProblemStage = " << problem_physics << " |\t Accepted values are 1->NS, 2->ADJ_NS.\n";
            }
        }
    }
}

template <class TElementData>
void TransportTopologyOptimizationElement<TElementData>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
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

        // Now the distinction of what to do for the primal or adjoint problem must be done calling the relative Kratos property
        unsigned int problem_physics = data.TopOptProblemStage;
        // problem_physics = 1: NS equations
        // problem_physics = 2: ADJOINT NS equations
        if (problem_physics == 1)
        {
            // PRIMAL TRANSPORT
            // Iterate over integration points to evaluate local contribution
            for (unsigned int g = 0; g < number_of_gauss_points; g++) 
            {
                this->UpdateIntegrationPointData(
                    data, g, gauss_weights[g],
                    row(shape_functions, g),shape_derivatives[g]);

                this->AddTimeIntegratedLHS(data, rLeftHandSideMatrix);
            }
        }
        else if (problem_physics == 2)
        {
            // ADJOINT TRANSPORT
            // Iterate over integration points to evaluate local contribution
            // PRIMAL TRANSPORT
            // Iterate over integration points to evaluate local contribution
            for (unsigned int g = 0; g < number_of_gauss_points; g++) 
            {
                this->UpdateIntegrationPointData(
                    data, g, gauss_weights[g],
                    row(shape_functions, g),shape_derivatives[g]);

                this->AddTimeIntegratedLHSAdjoint(data, rLeftHandSideMatrix);
            }
        }
        else
        {
            if (problem_physics != 0)
            {
                KRATOS_ERROR << "\nInvalid value for the variable TRANSPORT_TOP_OPT_PROBLEM_STAGE for the Fluid Topology Optimization aplication. |\t TopOptProblemStage = " << problem_physics << " |\t Accepted values are 1->NS, 2->ADJ_NS.\n";
        
            }
        }
    }
}

template <class TElementData>
void TransportTopologyOptimizationElement<TElementData>::CalculateRightHandSide(VectorType& rRightHandSideVector,
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

        // Now the distinction of what to do for the primal or adjoint problem must be done calling the relative Kratos property
        unsigned int problem_physics = data.TopOptProblemStage;
        // problem_physics = 1: NS equations
        // problem_physics = 2: ADJOINT NS equations
        if (problem_physics == 1)
        {
            // PRIMAL TRANSPORT
            // Iterate over integration points to evaluate local contribution
            for (unsigned int g = 0; g < number_of_gauss_points; g++) 
            {
                this->UpdateIntegrationPointData(
                    data, g, gauss_weights[g],
                    row(shape_functions, g),shape_derivatives[g]);

                this->AddTimeIntegratedRHS(data, rRightHandSideVector);
            }
        }
        else if (problem_physics == 2)
        {
            // ADJOINT TRANSPORT
            // Iterate over integration points to evaluate local contribution
            for (unsigned int g = 0; g < number_of_gauss_points; g++) 
            {
                this->UpdateIntegrationPointData(
                    data, g, gauss_weights[g],
                    row(shape_functions, g),shape_derivatives[g]);

                this->AddTimeIntegratedRHSAdjoint(data, rRightHandSideVector);
            }
        }
        else
        {
            if (problem_physics != 0)
            {
                KRATOS_ERROR << "\nInvalid value for the variable TRANSPORT_TOP_OPT_PROBLEM_STAGE for the Fluid Topology Optimization aplication. |\t TopOptProblemStage = " << problem_physics << " |\t Accepted values are 1->NS, 2->ADJ_NS.\n";
        
            }
        }
    }
}

template <class TElementData>
void TransportTopologyOptimizationElement<TElementData>::CalculateLocalVelocityContribution(
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
void TransportTopologyOptimizationElement<TElementData>::CalculateMassMatrix(MatrixType& rMassMatrix,
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
void TransportTopologyOptimizationElement< TElementData >::EquationIdVector(EquationIdVectorType &rResult, const ProcessInfo &rCurrentProcessInfo) const
{
    const GeometryType& r_geometry = this->GetGeometry();

    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    // Now the distinction of what to do for the primal or adjoint problem must be done calling the relative Kratos property
    unsigned int problem_physics = rCurrentProcessInfo[TRANSPORT_TOP_OPT_PROBLEM_STAGE];
    // problem_physics = 1: NS equations
    // problem_physics = 2: ADJOINT NS equations
    if (problem_physics == 1)
    {
        // PRIMAL TRANSPORT
        const unsigned int cpos = this->GetGeometry()[0].GetDofPosition(TEMPERATURE);

        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            rResult[LocalIndex++] = r_geometry[i].GetDof(TEMPERATURE, cpos).EquationId();
        }
    }
    else if (problem_physics == 2)
    {
        // ADJOINT TRANSPORT
        const unsigned int cpos = this->GetGeometry()[0].GetDofPosition(TEMPERATURE_ADJ);

        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            rResult[LocalIndex++] = r_geometry[i].GetDof(TEMPERATURE_ADJ, cpos).EquationId();
        }
    }
    else
    {
        if (problem_physics != 0)
        {
            KRATOS_ERROR << "\nInvalid value for the variable TRANSPORT_TOP_OPT_PROBLEM_STAGE for the Fluid Topology Optimization aplication. |\t TopOptProblemStage = " << problem_physics << " |\t Accepted values are 1->T, 2->ADJ_T.\n";
        }
    }
}
template< class TElementData >
void TransportTopologyOptimizationElement< TElementData >::GetDofList(DofsVectorType &rElementalDofList, const ProcessInfo &rCurrentProcessInfo) const
{
    const GeometryType& r_geometry = this->GetGeometry();

    if (rElementalDofList.size() != LocalSize)
         rElementalDofList.resize(LocalSize);

    // Now the distinction of what to do for the primal or adjoint problem must be done calling the relative Kratos property
    unsigned int problem_physics = rCurrentProcessInfo[TRANSPORT_TOP_OPT_PROBLEM_STAGE];
    // problem_physics = 1: NS equations
    // problem_physics = 2: ADJOINT NS equations
    if (problem_physics == 1)
    {
        // PRIMAL TRANSPORT
        const unsigned int cpos = this->GetGeometry()[0].GetDofPosition(TEMPERATURE);

        unsigned int LocalIndex = 0;
        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(TEMPERATURE,cpos);
        }
    }
    else if (problem_physics == 2)
    {
        // ADJOINT TRANSPORT
        const unsigned int cpos = this->GetGeometry()[0].GetDofPosition(TEMPERATURE_ADJ);

        unsigned int LocalIndex = 0;
        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(TEMPERATURE_ADJ,cpos);
        }
    }
    else
    {
        if (problem_physics != 0)
        {
            KRATOS_ERROR << "\nInvalid value for the variable TRANSPORT_TOP_OPT_PROBLEM_STAGE for the Fluid Topology Optimization aplication. |\t TopOptProblemStage = " << problem_physics << " |\t Accepted values are 1->NS, 2->ADJ_NS.\n";
    
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void TransportTopologyOptimizationElement<TElementData>::GetFirstDerivativesVector(Vector &rValues, int Step) const
{
    const GeometryType& r_geometry = this->GetGeometry();

    if (rValues.size() != LocalSize)
        rValues.resize(LocalSize,false);

    unsigned int Index = 0;

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        rValues[Index++] = r_geometry[i].FastGetSolutionStepValue(TEMPERATURE,Step);
    }
}


template< class TElementData >
void TransportTopologyOptimizationElement<TElementData>::GetSecondDerivativesVector(Vector &rValues, int Step) const
{
    const GeometryType& r_geometry = this->GetGeometry();

    if (rValues.size() != LocalSize)
        rValues.resize(LocalSize,false);

    unsigned int index = 0;

    for (unsigned int i = 0; i < NumNodes; i++) {
        const array_1d<double,3>& r_acceleration = r_geometry[i].FastGetSolutionStepValue(ACCELERATION,Step);
        for (unsigned int d = 0; d < Dim; d++)
            rValues[index++] = r_acceleration[d];
        rValues[index++] = 0.0; // skip pressure Dof
    }
}


template< class TElementData >
GeometryData::IntegrationMethod TransportTopologyOptimizationElement<TElementData>::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_2;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template< class TElementData >
int TransportTopologyOptimizationElement<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
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

    unsigned int problem_physics = rCurrentProcessInfo[TRANSPORT_TOP_OPT_PROBLEM_STAGE];
    if (problem_physics == 1)
    {
        for(unsigned int i=0; i<NumNodes; ++i)
        {
            const Node& rNode = r_geometry[i];
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TEMPERATURE,rNode);
        }
    }
    else if (problem_physics == 2)
    {
        for(unsigned int i=0; i<NumNodes; ++i)
        {
            const Node& rNode = r_geometry[i];
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TEMPERATURE_ADJ,rNode);
        }
    }
    else
    {
        if (problem_physics != 0)
        {
            KRATOS_ERROR << "\nInvalid value for the variable TRANSPORT_TOP_OPT_PROBLEM_STAGE for the Fluid Topology Optimization aplication. |\t TopOptProblemStage = " << problem_physics << " |\t Accepted values are 1->NS, 2->ADJ_NS.\n";
    
        }    
    }

    // If this is a 2D problem, check that nodes are in XY plane
    if ( Dim == 2)
    {
        for (unsigned int i=0; i<NumNodes; ++i) {
            if (r_geometry[i].Z() != 0.0)
                KRATOS_ERROR << "Node " << r_geometry[i].Id() << "has non-zero Z coordinate." << std::endl;
        }
    }

    // Check the constitutive law
    // KRATOS_ERROR_IF(mpConstitutiveLaw == nullptr) << "Constitutive Law not initialized for Element " << this->Info() << std::endl;

    // constexpr auto dimension = Dim;  // I need to set this here otherwise it gives me a linking error when attempting to '<<' Dim.

    // KRATOS_ERROR_IF(mpConstitutiveLaw->WorkingSpaceDimension() != Dim)
    //     << "Wrong dimension: The " << mpConstitutiveLaw->WorkingSpaceDimension()
    //     << "D constitutive law " << mpConstitutiveLaw->Info()
    //     << " is not compatible with " << dimension << "D element " << this->Info()
    //     << "." << std::endl;

    // out = mpConstitutiveLaw->Check(this->GetProperties(),r_geometry,rCurrentProcessInfo);
    // KRATOS_ERROR_IF_NOT( out == 0) << "The Constitutive Law provided for Element " << this->Info() << " is not correct." << std::endl;

    return out;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void TransportTopologyOptimizationElement<TElementData>::CalculateOnIntegrationPoints(
    Variable<array_1d<double, 3 > > const& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    // if (rVariable == VORTICITY)
    // {
    //     // Get Shape function data
    //     Vector gauss_weights;
    //     Matrix shape_functions;
    //     ShapeFunctionDerivativesArrayType shape_function_gradients;
    //     this->CalculateGeometryData(gauss_weights,shape_functions,shape_function_gradients);

    //     VorticityUtilities<Dim>::CalculateVorticityVector(this->GetGeometry(),shape_function_gradients,rValues);
    // }
}


template< class TElementData >
void TransportTopologyOptimizationElement<TElementData>::CalculateOnIntegrationPoints(
    Variable<double> const& rVariable,
    std::vector<double>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    // if (rVariable == Q_VALUE)
    // {
    //     // Get Shape function data
    //     Vector gauss_weights;
    //     Matrix shape_functions;
    //     ShapeFunctionDerivativesArrayType shape_function_gradients;
    //     this->CalculateGeometryData(gauss_weights,shape_functions,shape_function_gradients);

    //     VorticityUtilities<Dim>::CalculateQValue(this->GetGeometry(),shape_function_gradients,rValues);
    // }
    // else if (rVariable == VORTICITY_MAGNITUDE)
    // {
    //     // Get Shape function data
    //     Vector gauss_weights;
    //     Matrix shape_functions;
    //     ShapeFunctionDerivativesArrayType shape_function_gradients;
    //     this->CalculateGeometryData(gauss_weights,shape_functions,shape_function_gradients);

    //     VorticityUtilities<Dim>::CalculateVorticityMagnitude(this->GetGeometry(),shape_function_gradients,rValues);
    // }
    // else if (rVariable == UPDATE_STATISTICS)
    // {
    //     KRATOS_DEBUG_ERROR_IF(!rCurrentProcessInfo.Has(STATISTICS_CONTAINER)) << "Trying to compute turbulent statistics, but ProcessInfo does not have STATISTICS_CONTAINER defined." << std::endl;
    //     rCurrentProcessInfo.GetValue(STATISTICS_CONTAINER)->UpdateStatistics(this);
    // }
}

template <class TElementData>
void TransportTopologyOptimizationElement<TElementData>::CalculateOnIntegrationPoints(
    Variable<array_1d<double, 6>> const& rVariable,
    std::vector<array_1d<double, 6>>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{}

template <class TElementData>
void TransportTopologyOptimizationElement<TElementData>::CalculateOnIntegrationPoints(
    Variable<Vector> const& rVariable,
    std::vector<Vector>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{}

template <class TElementData>
void TransportTopologyOptimizationElement<TElementData>::CalculateOnIntegrationPoints(
    Variable<Matrix> const& rVariable,
    std::vector<Matrix>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output
template <class TElementData>
void TransportTopologyOptimizationElement<TElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info() << std::endl;

    // if (this->GetConstitutiveLaw() != nullptr) {
    //     rOStream << "with constitutive law " << std::endl;
    //     this->GetConstitutiveLaw()->PrintInfo(rOStream);
    // }
}

template <class TElementData>
std::string TransportTopologyOptimizationElement<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "TransportTopologyOptimizationElement" << Dim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template <class TElementData>
double TransportTopologyOptimizationElement<TElementData>::GetAtCoordinate(
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
array_1d<double, 3> TransportTopologyOptimizationElement<TElementData>::GetAtCoordinate(
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
BoundedMatrix<double, TElementData::Dim, TElementData::Dim> TransportTopologyOptimizationElement<TElementData>::GetAtCoordinate(
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
double TransportTopologyOptimizationElement<TElementData>::GetAtCoordinate(
    const double Value,
    const typename TElementData::ShapeFunctionsType& rN) const
{
    return Value;
}

template <class TElementData>
void TransportTopologyOptimizationElement<TElementData>::UpdateIntegrationPointData(
    TElementData& rData,
    unsigned int IntegrationPointIndex,
    double Weight,
    const typename TElementData::MatrixRowType& rN,
    const typename TElementData::ShapeDerivativesType& rDN_DX) const
{
    rData.UpdateGeometryValues(IntegrationPointIndex, Weight, rN, rDN_DX);
}

template< class TElementData >
void TransportTopologyOptimizationElement<TElementData>::CalculateGeometryData(Vector &rGaussWeights,
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
void TransportTopologyOptimizationElement<TElementData>::Calculate(
    const Variable<double> &rVariable,
    double &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    Element::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TElementData>
void TransportTopologyOptimizationElement<TElementData>::Calculate(
    const Variable<array_1d<double, 3>> &rVariable,
    array_1d<double, 3> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    Element::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TElementData>
void TransportTopologyOptimizationElement<TElementData>::Calculate(
    const Variable<Vector >& rVariable,
    Vector& rOutput,
    const ProcessInfo& rCurrentProcessInfo )
{
    Element::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TElementData>
void TransportTopologyOptimizationElement<TElementData>::Calculate(
    const Variable<Matrix> &rVariable,
    Matrix &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    Element::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void TransportTopologyOptimizationElement<TElementData>::ConvectionOperator(Vector &rResult,
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
void TransportTopologyOptimizationElement<TElementData>::AddTimeIntegratedSystem(
    TElementData& rData, MatrixType& rLHS, VectorType& rRHS) 
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void TransportTopologyOptimizationElement<TElementData>::AddTimeIntegratedLHS(
    TElementData& rData, MatrixType& rLHS) 
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
}

template <class TElementData>
void TransportTopologyOptimizationElement<TElementData>::AddTimeIntegratedRHS(
    TElementData& rData, VectorType& rRHS) 
{
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void TransportTopologyOptimizationElement<TElementData>::AddTimeIntegratedSystemAdjoint(
    TElementData& rData, MatrixType& rLHS, VectorType& rRHS) 
{
    this->ComputeGaussPointLHSContributionAdjoint(rData, rLHS);
    this->ComputeGaussPointRHSContributionAdjoint(rData, rRHS);
}

template <class TElementData>
void TransportTopologyOptimizationElement<TElementData>::AddTimeIntegratedLHSAdjoint(
    TElementData& rData, MatrixType& rLHS) 
{
    this->ComputeGaussPointLHSContributionAdjoint(rData, rLHS);
}

template <class TElementData>
void TransportTopologyOptimizationElement<TElementData>::AddTimeIntegratedRHSAdjoint(
    TElementData& rData, VectorType& rRHS) 
{
    this->ComputeGaussPointRHSContributionAdjoint(rData, rRHS);
}

// BUILD NS SYSTEM
template <>
void TransportTopologyOptimizationElement< TransportTopologyOptimizationElementData<2,3,true>>::ComputeGaussPointLHSContribution(
    TransportTopologyOptimizationElementData<2,3,true> & rData,
    MatrixType& rLHS)
{
    const array_1d<double,3> D = rData.Conductivity;
    const array_1d<double,3> k = rData.Decay;
    const array_1d<double,3> c = rData.ConvectionCoefficient;

    const double h = rData.ElementSize;
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // const double dyn_tau = rData.DynamicTau;
    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const BoundedMatrix<double,2,3> vconv = rData.ConvectiveVelocity - rData.MeshVelocity;

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;

    // TRANSPORT ELEMENTAL LHS MATRIX
    const double crLHS0 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2];
const double crLHS1 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2];
const double crLHS2 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crLHS3 = DN(0,0)*crLHS2;
const double crLHS4 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crLHS5 = DN(0,1)*crLHS4;
const double crLHS6 = crLHS3 + crLHS5;
const double crLHS7 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2];
const double crLHS8 = N[0]*crLHS7;
const double crLHS9 = N[0]*crLHS0;
const double crLHS10 = crLHS3*crLHS7;
const double crLHS11 = crLHS5*crLHS7;
const double crLHS12 = crLHS10 + crLHS11 + crLHS9;
const double crLHS13 = 1.0/(crLHS0*stab_c3 + crLHS1*stab_c1*1.0/(h*h) + stab_c2*fabs(crLHS7)*1.0/h*sqrt(crLHS2*crLHS2 + crLHS4*crLHS4));
const double crLHS14 = 1.0*crLHS13;
const double crLHS15 = crLHS12*crLHS14;
const double crLHS16 = DN(1,0)*crLHS2;
const double crLHS17 = DN(1,1)*crLHS4;
const double crLHS18 = crLHS16 + crLHS17;
const double crLHS19 = N[1]*crLHS0;
const double crLHS20 = crLHS16*crLHS7;
const double crLHS21 = crLHS17*crLHS7;
const double crLHS22 = crLHS19 + crLHS20 + crLHS21;
const double crLHS23 = crLHS14*crLHS22;
const double crLHS24 = N[1]*crLHS9 + crLHS1*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1));
const double crLHS25 = DN(2,0)*crLHS2;
const double crLHS26 = DN(2,1)*crLHS4;
const double crLHS27 = crLHS25 + crLHS26;
const double crLHS28 = crLHS25*crLHS7;
const double crLHS29 = crLHS26*crLHS7;
const double crLHS30 = N[2]*crLHS0 + crLHS28 + crLHS29;
const double crLHS31 = crLHS14*crLHS30;
const double crLHS32 = N[2]*crLHS9 + crLHS1*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1));
const double crLHS33 = N[1]*crLHS7;
const double crLHS34 = N[2]*crLHS19 + crLHS1*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1));
const double crLHS35 = N[2]*crLHS7;
rLHS(0,0)+=-gauss_weight*(1.0*N[0]*crLHS0*crLHS12*crLHS13 - crLHS0*N[0]*N[0] - crLHS1*(DN(0,0)*DN(0,0) + DN(0,1)*DN(0,1)) - crLHS10*crLHS15 - crLHS11*crLHS15 - crLHS6*crLHS8);
rLHS(0,1)+=-gauss_weight*(1.0*N[0]*crLHS0*crLHS13*crLHS22 - crLHS10*crLHS23 - crLHS11*crLHS23 - crLHS18*crLHS8 - crLHS24);
rLHS(0,2)+=-gauss_weight*(1.0*N[0]*crLHS0*crLHS13*crLHS30 - crLHS10*crLHS31 - crLHS11*crLHS31 - crLHS27*crLHS8 - crLHS32);
rLHS(1,0)+=-gauss_weight*(1.0*N[1]*crLHS0*crLHS12*crLHS13 - crLHS15*crLHS20 - crLHS15*crLHS21 - crLHS24 - crLHS33*crLHS6);
rLHS(1,1)+=-gauss_weight*(1.0*N[1]*crLHS0*crLHS13*crLHS22 - crLHS0*N[1]*N[1] - crLHS1*(DN(1,0)*DN(1,0) + DN(1,1)*DN(1,1)) - crLHS18*crLHS33 - crLHS20*crLHS23 - crLHS21*crLHS23);
rLHS(1,2)+=-gauss_weight*(1.0*N[1]*crLHS0*crLHS13*crLHS30 - crLHS20*crLHS31 - crLHS21*crLHS31 - crLHS27*crLHS33 - crLHS34);
rLHS(2,0)+=-gauss_weight*(1.0*N[2]*crLHS0*crLHS12*crLHS13 - crLHS15*crLHS28 - crLHS15*crLHS29 - crLHS32 - crLHS35*crLHS6);
rLHS(2,1)+=-gauss_weight*(1.0*N[2]*crLHS0*crLHS13*crLHS22 - crLHS18*crLHS35 - crLHS23*crLHS28 - crLHS23*crLHS29 - crLHS34);
rLHS(2,2)+=-gauss_weight*(1.0*N[2]*crLHS0*crLHS13*crLHS30 - crLHS0*N[2]*N[2] - crLHS1*(DN(2,0)*DN(2,0) + DN(2,1)*DN(2,1)) - crLHS27*crLHS35 - crLHS28*crLHS31 - crLHS29*crLHS31);
 
    
}

template <>
void TransportTopologyOptimizationElement< TransportTopologyOptimizationElementData<2,4,true>>::ComputeGaussPointLHSContribution(
    TransportTopologyOptimizationElementData<2,4,true> & rData,
    MatrixType& rLHS)
{
    const array_1d<double,4> D = rData.Conductivity;
    const array_1d<double,4> k = rData.Decay;
    const array_1d<double,4> c = rData.ConvectionCoefficient;

    const double h = rData.ElementSize;
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,2>& DN = rData.DN_DX;

    // const double dyn_tau = rData.DynamicTau;
    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const BoundedMatrix<double,2,4> vconv = rData.ConvectiveVelocity - rData.MeshVelocity;

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;

    // TRANSPORT ELEMENTAL LHS MATRIX
    const double crLHS0 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2] + N[3]*k[3];
const double crLHS1 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2] + D[3]*N[3];
const double crLHS2 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crLHS3 = DN(0,0)*crLHS2;
const double crLHS4 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crLHS5 = DN(0,1)*crLHS4;
const double crLHS6 = crLHS3 + crLHS5;
const double crLHS7 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3];
const double crLHS8 = N[0]*crLHS7;
const double crLHS9 = N[0]*crLHS0;
const double crLHS10 = crLHS3*crLHS7;
const double crLHS11 = crLHS5*crLHS7;
const double crLHS12 = crLHS10 + crLHS11 + crLHS9;
const double crLHS13 = 1.0/(crLHS0*stab_c3 + crLHS1*stab_c1*1.0/(h*h) + stab_c2*fabs(crLHS7)*1.0/h*sqrt(crLHS2*crLHS2 + crLHS4*crLHS4));
const double crLHS14 = 1.0*crLHS13;
const double crLHS15 = crLHS12*crLHS14;
const double crLHS16 = DN(1,0)*crLHS2;
const double crLHS17 = DN(1,1)*crLHS4;
const double crLHS18 = crLHS16 + crLHS17;
const double crLHS19 = N[1]*crLHS0;
const double crLHS20 = crLHS16*crLHS7;
const double crLHS21 = crLHS17*crLHS7;
const double crLHS22 = crLHS19 + crLHS20 + crLHS21;
const double crLHS23 = crLHS14*crLHS22;
const double crLHS24 = N[1]*crLHS9 + crLHS1*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1));
const double crLHS25 = DN(2,0)*crLHS2;
const double crLHS26 = DN(2,1)*crLHS4;
const double crLHS27 = crLHS25 + crLHS26;
const double crLHS28 = N[2]*crLHS0;
const double crLHS29 = crLHS25*crLHS7;
const double crLHS30 = crLHS26*crLHS7;
const double crLHS31 = crLHS28 + crLHS29 + crLHS30;
const double crLHS32 = crLHS14*crLHS31;
const double crLHS33 = N[2]*crLHS9 + crLHS1*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1));
const double crLHS34 = DN(3,0)*crLHS2;
const double crLHS35 = DN(3,1)*crLHS4;
const double crLHS36 = crLHS34 + crLHS35;
const double crLHS37 = crLHS34*crLHS7;
const double crLHS38 = crLHS35*crLHS7;
const double crLHS39 = N[3]*crLHS0 + crLHS37 + crLHS38;
const double crLHS40 = crLHS14*crLHS39;
const double crLHS41 = N[3]*crLHS9 + crLHS1*(DN(0,0)*DN(3,0) + DN(0,1)*DN(3,1));
const double crLHS42 = N[1]*crLHS7;
const double crLHS43 = N[2]*crLHS19 + crLHS1*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1));
const double crLHS44 = N[3]*crLHS19 + crLHS1*(DN(1,0)*DN(3,0) + DN(1,1)*DN(3,1));
const double crLHS45 = N[2]*crLHS7;
const double crLHS46 = N[3]*crLHS28 + crLHS1*(DN(2,0)*DN(3,0) + DN(2,1)*DN(3,1));
const double crLHS47 = N[3]*crLHS7;
rLHS(0,0)+=-gauss_weight*(1.0*N[0]*crLHS0*crLHS12*crLHS13 - crLHS0*N[0]*N[0] - crLHS1*(DN(0,0)*DN(0,0) + DN(0,1)*DN(0,1)) - crLHS10*crLHS15 - crLHS11*crLHS15 - crLHS6*crLHS8);
rLHS(0,1)+=-gauss_weight*(1.0*N[0]*crLHS0*crLHS13*crLHS22 - crLHS10*crLHS23 - crLHS11*crLHS23 - crLHS18*crLHS8 - crLHS24);
rLHS(0,2)+=-gauss_weight*(1.0*N[0]*crLHS0*crLHS13*crLHS31 - crLHS10*crLHS32 - crLHS11*crLHS32 - crLHS27*crLHS8 - crLHS33);
rLHS(0,3)+=-gauss_weight*(1.0*N[0]*crLHS0*crLHS13*crLHS39 - crLHS10*crLHS40 - crLHS11*crLHS40 - crLHS36*crLHS8 - crLHS41);
rLHS(1,0)+=-gauss_weight*(1.0*N[1]*crLHS0*crLHS12*crLHS13 - crLHS15*crLHS20 - crLHS15*crLHS21 - crLHS24 - crLHS42*crLHS6);
rLHS(1,1)+=-gauss_weight*(1.0*N[1]*crLHS0*crLHS13*crLHS22 - crLHS0*N[1]*N[1] - crLHS1*(DN(1,0)*DN(1,0) + DN(1,1)*DN(1,1)) - crLHS18*crLHS42 - crLHS20*crLHS23 - crLHS21*crLHS23);
rLHS(1,2)+=-gauss_weight*(1.0*N[1]*crLHS0*crLHS13*crLHS31 - crLHS20*crLHS32 - crLHS21*crLHS32 - crLHS27*crLHS42 - crLHS43);
rLHS(1,3)+=-gauss_weight*(1.0*N[1]*crLHS0*crLHS13*crLHS39 - crLHS20*crLHS40 - crLHS21*crLHS40 - crLHS36*crLHS42 - crLHS44);
rLHS(2,0)+=-gauss_weight*(1.0*N[2]*crLHS0*crLHS12*crLHS13 - crLHS15*crLHS29 - crLHS15*crLHS30 - crLHS33 - crLHS45*crLHS6);
rLHS(2,1)+=-gauss_weight*(1.0*N[2]*crLHS0*crLHS13*crLHS22 - crLHS18*crLHS45 - crLHS23*crLHS29 - crLHS23*crLHS30 - crLHS43);
rLHS(2,2)+=-gauss_weight*(1.0*N[2]*crLHS0*crLHS13*crLHS31 - crLHS0*N[2]*N[2] - crLHS1*(DN(2,0)*DN(2,0) + DN(2,1)*DN(2,1)) - crLHS27*crLHS45 - crLHS29*crLHS32 - crLHS30*crLHS32);
rLHS(2,3)+=-gauss_weight*(1.0*N[2]*crLHS0*crLHS13*crLHS39 - crLHS29*crLHS40 - crLHS30*crLHS40 - crLHS36*crLHS45 - crLHS46);
rLHS(3,0)+=-gauss_weight*(1.0*N[3]*crLHS0*crLHS12*crLHS13 - crLHS15*crLHS37 - crLHS15*crLHS38 - crLHS41 - crLHS47*crLHS6);
rLHS(3,1)+=-gauss_weight*(1.0*N[3]*crLHS0*crLHS13*crLHS22 - crLHS18*crLHS47 - crLHS23*crLHS37 - crLHS23*crLHS38 - crLHS44);
rLHS(3,2)+=-gauss_weight*(1.0*N[3]*crLHS0*crLHS13*crLHS31 - crLHS27*crLHS47 - crLHS32*crLHS37 - crLHS32*crLHS38 - crLHS46);
rLHS(3,3)+=-gauss_weight*(1.0*N[3]*crLHS0*crLHS13*crLHS39 - crLHS0*N[3]*N[3] - crLHS1*(DN(3,0)*DN(3,0) + DN(3,1)*DN(3,1)) - crLHS36*crLHS47 - crLHS37*crLHS40 - crLHS38*crLHS40);
 
    
}

template <>
void TransportTopologyOptimizationElement<TransportTopologyOptimizationElementData<3,4,true>>::ComputeGaussPointLHSContribution(
    TransportTopologyOptimizationElementData<3,4,true>& rData,
    MatrixType& rLHS)
{
    const array_1d<double,4> D = rData.Conductivity;
    const array_1d<double,4> k = rData.Decay;
    const array_1d<double,4> c = rData.ConvectionCoefficient;

    const double h = rData.ElementSize;
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,3>& DN = rData.DN_DX;

    // const double dyn_tau = rData.DynamicTau;
    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const BoundedMatrix<double,3,4> vconv = rData.ConvectiveVelocity - rData.MeshVelocity;

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;

    // TRANSPORT ELEMENTAL LHS MATRIX
    const double crLHS0 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2] + N[3]*k[3];
const double crLHS1 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2] + D[3]*N[3];
const double crLHS2 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crLHS3 = DN(0,0)*crLHS2;
const double crLHS4 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crLHS5 = DN(0,1)*crLHS4;
const double crLHS6 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crLHS7 = DN(0,2)*crLHS6;
const double crLHS8 = crLHS3 + crLHS5 + crLHS7;
const double crLHS9 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3];
const double crLHS10 = N[0]*crLHS9;
const double crLHS11 = N[0]*crLHS0;
const double crLHS12 = crLHS3*crLHS9;
const double crLHS13 = crLHS5*crLHS9;
const double crLHS14 = crLHS7*crLHS9;
const double crLHS15 = crLHS11 + crLHS12 + crLHS13 + crLHS14;
const double crLHS16 = 1.0/(crLHS0*stab_c3 + crLHS1*stab_c1*1.0/(h*h) + stab_c2*fabs(crLHS9)*1.0/h*sqrt(crLHS2*crLHS2 + crLHS4*crLHS4 + crLHS6*crLHS6));
const double crLHS17 = 1.0*crLHS16;
const double crLHS18 = crLHS15*crLHS17;
const double crLHS19 = DN(1,0)*crLHS2;
const double crLHS20 = DN(1,1)*crLHS4;
const double crLHS21 = DN(1,2)*crLHS6;
const double crLHS22 = crLHS19 + crLHS20 + crLHS21;
const double crLHS23 = N[1]*crLHS0;
const double crLHS24 = crLHS19*crLHS9;
const double crLHS25 = crLHS20*crLHS9;
const double crLHS26 = crLHS21*crLHS9;
const double crLHS27 = crLHS23 + crLHS24 + crLHS25 + crLHS26;
const double crLHS28 = crLHS17*crLHS27;
const double crLHS29 = N[1]*crLHS11 + crLHS1*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1) + DN(0,2)*DN(1,2));
const double crLHS30 = DN(2,0)*crLHS2;
const double crLHS31 = DN(2,1)*crLHS4;
const double crLHS32 = DN(2,2)*crLHS6;
const double crLHS33 = crLHS30 + crLHS31 + crLHS32;
const double crLHS34 = N[2]*crLHS0;
const double crLHS35 = crLHS30*crLHS9;
const double crLHS36 = crLHS31*crLHS9;
const double crLHS37 = crLHS32*crLHS9;
const double crLHS38 = crLHS34 + crLHS35 + crLHS36 + crLHS37;
const double crLHS39 = crLHS17*crLHS38;
const double crLHS40 = N[2]*crLHS11 + crLHS1*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1) + DN(0,2)*DN(2,2));
const double crLHS41 = DN(3,0)*crLHS2;
const double crLHS42 = DN(3,1)*crLHS4;
const double crLHS43 = DN(3,2)*crLHS6;
const double crLHS44 = crLHS41 + crLHS42 + crLHS43;
const double crLHS45 = crLHS41*crLHS9;
const double crLHS46 = crLHS42*crLHS9;
const double crLHS47 = crLHS43*crLHS9;
const double crLHS48 = N[3]*crLHS0 + crLHS45 + crLHS46 + crLHS47;
const double crLHS49 = crLHS17*crLHS48;
const double crLHS50 = N[3]*crLHS11 + crLHS1*(DN(0,0)*DN(3,0) + DN(0,1)*DN(3,1) + DN(0,2)*DN(3,2));
const double crLHS51 = N[1]*crLHS9;
const double crLHS52 = N[2]*crLHS23 + crLHS1*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1) + DN(1,2)*DN(2,2));
const double crLHS53 = N[3]*crLHS23 + crLHS1*(DN(1,0)*DN(3,0) + DN(1,1)*DN(3,1) + DN(1,2)*DN(3,2));
const double crLHS54 = N[2]*crLHS9;
const double crLHS55 = N[3]*crLHS34 + crLHS1*(DN(2,0)*DN(3,0) + DN(2,1)*DN(3,1) + DN(2,2)*DN(3,2));
const double crLHS56 = N[3]*crLHS9;
rLHS(0,0)+=-gauss_weight*(1.0*N[0]*crLHS0*crLHS15*crLHS16 - crLHS0*N[0]*N[0] - crLHS1*(DN(0,0)*DN(0,0) + DN(0,1)*DN(0,1) + DN(0,2)*DN(0,2)) - crLHS10*crLHS8 - crLHS12*crLHS18 - crLHS13*crLHS18 - crLHS14*crLHS18);
rLHS(0,1)+=-gauss_weight*(1.0*N[0]*crLHS0*crLHS16*crLHS27 - crLHS10*crLHS22 - crLHS12*crLHS28 - crLHS13*crLHS28 - crLHS14*crLHS28 - crLHS29);
rLHS(0,2)+=-gauss_weight*(1.0*N[0]*crLHS0*crLHS16*crLHS38 - crLHS10*crLHS33 - crLHS12*crLHS39 - crLHS13*crLHS39 - crLHS14*crLHS39 - crLHS40);
rLHS(0,3)+=-gauss_weight*(1.0*N[0]*crLHS0*crLHS16*crLHS48 - crLHS10*crLHS44 - crLHS12*crLHS49 - crLHS13*crLHS49 - crLHS14*crLHS49 - crLHS50);
rLHS(1,0)+=-gauss_weight*(1.0*N[1]*crLHS0*crLHS15*crLHS16 - crLHS18*crLHS24 - crLHS18*crLHS25 - crLHS18*crLHS26 - crLHS29 - crLHS51*crLHS8);
rLHS(1,1)+=-gauss_weight*(1.0*N[1]*crLHS0*crLHS16*crLHS27 - crLHS0*N[1]*N[1] - crLHS1*(DN(1,0)*DN(1,0) + DN(1,1)*DN(1,1) + DN(1,2)*DN(1,2)) - crLHS22*crLHS51 - crLHS24*crLHS28 - crLHS25*crLHS28 - crLHS26*crLHS28);
rLHS(1,2)+=-gauss_weight*(1.0*N[1]*crLHS0*crLHS16*crLHS38 - crLHS24*crLHS39 - crLHS25*crLHS39 - crLHS26*crLHS39 - crLHS33*crLHS51 - crLHS52);
rLHS(1,3)+=-gauss_weight*(1.0*N[1]*crLHS0*crLHS16*crLHS48 - crLHS24*crLHS49 - crLHS25*crLHS49 - crLHS26*crLHS49 - crLHS44*crLHS51 - crLHS53);
rLHS(2,0)+=-gauss_weight*(1.0*N[2]*crLHS0*crLHS15*crLHS16 - crLHS18*crLHS35 - crLHS18*crLHS36 - crLHS18*crLHS37 - crLHS40 - crLHS54*crLHS8);
rLHS(2,1)+=-gauss_weight*(1.0*N[2]*crLHS0*crLHS16*crLHS27 - crLHS22*crLHS54 - crLHS28*crLHS35 - crLHS28*crLHS36 - crLHS28*crLHS37 - crLHS52);
rLHS(2,2)+=-gauss_weight*(1.0*N[2]*crLHS0*crLHS16*crLHS38 - crLHS0*N[2]*N[2] - crLHS1*(DN(2,0)*DN(2,0) + DN(2,1)*DN(2,1) + DN(2,2)*DN(2,2)) - crLHS33*crLHS54 - crLHS35*crLHS39 - crLHS36*crLHS39 - crLHS37*crLHS39);
rLHS(2,3)+=-gauss_weight*(1.0*N[2]*crLHS0*crLHS16*crLHS48 - crLHS35*crLHS49 - crLHS36*crLHS49 - crLHS37*crLHS49 - crLHS44*crLHS54 - crLHS55);
rLHS(3,0)+=-gauss_weight*(1.0*N[3]*crLHS0*crLHS15*crLHS16 - crLHS18*crLHS45 - crLHS18*crLHS46 - crLHS18*crLHS47 - crLHS50 - crLHS56*crLHS8);
rLHS(3,1)+=-gauss_weight*(1.0*N[3]*crLHS0*crLHS16*crLHS27 - crLHS22*crLHS56 - crLHS28*crLHS45 - crLHS28*crLHS46 - crLHS28*crLHS47 - crLHS53);
rLHS(3,2)+=-gauss_weight*(1.0*N[3]*crLHS0*crLHS16*crLHS38 - crLHS33*crLHS56 - crLHS39*crLHS45 - crLHS39*crLHS46 - crLHS39*crLHS47 - crLHS55);
rLHS(3,3)+=-gauss_weight*(1.0*N[3]*crLHS0*crLHS16*crLHS48 - crLHS0*N[3]*N[3] - crLHS1*(DN(3,0)*DN(3,0) + DN(3,1)*DN(3,1) + DN(3,2)*DN(3,2)) - crLHS44*crLHS56 - crLHS45*crLHS49 - crLHS46*crLHS49 - crLHS47*crLHS49);

}

template <>
void TransportTopologyOptimizationElement<TransportTopologyOptimizationElementData<2,3,true>>::ComputeGaussPointRHSContribution(
    TransportTopologyOptimizationElementData<2,3,true>& rData,
    VectorType& rRHS)
{
    const array_1d<double,3> D = rData.Conductivity;
    const array_1d<double,3> k = rData.Decay;
    const array_1d<double,3> c = rData.ConvectionCoefficient;

    const double h = rData.ElementSize;
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // const double dyn_tau = rData.DynamicTau;
    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const array_1d<double,3>& t = rData.Temperature;
    const array_1d<double,3>& tn = rData.Temperature_OldStep1;
    const array_1d<double,3>& tnn = rData.Temperature_OldStep2;
    const BoundedMatrix<double,2,3>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,2,3> vconv = rData.ConvectiveVelocity - rData.MeshVelocity;
    const array_1d<double,3>& source = rData.SourceTerm;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // TRANSPORT ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*source[0] + N[1]*source[1] + N[2]*source[2];
const double crRHS1 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2];
const double crRHS2 = crRHS1*(N[0]*t[0] + N[1]*t[1] + N[2]*t[2]);
const double crRHS3 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2];
const double crRHS4 = DN(0,0)*t[0] + DN(1,0)*t[1] + DN(2,0)*t[2];
const double crRHS5 = DN(0,1)*t[0] + DN(1,1)*t[1] + DN(2,1)*t[2];
const double crRHS6 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crRHS7 = crRHS4*crRHS6;
const double crRHS8 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crRHS9 = crRHS5*crRHS8;
const double crRHS10 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2];
const double crRHS11 = crRHS10*(crRHS7 + crRHS9);
const double crRHS12 = 1.0*(crRHS0 - crRHS10*crRHS7 - crRHS10*crRHS9 - crRHS2)*1.0/(crRHS1*stab_c3 + crRHS3*stab_c1*1.0/(h*h) + stab_c2*fabs(crRHS10)*1.0/h*sqrt(crRHS6*crRHS6 + crRHS8*crRHS8));
const double crRHS13 = crRHS1*crRHS12;
const double crRHS14 = crRHS10*crRHS12;
const double crRHS15 = crRHS14*crRHS6;
const double crRHS16 = crRHS14*crRHS8;
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS15 - DN(0,1)*crRHS16 - N[0]*crRHS0 + N[0]*crRHS11 + N[0]*crRHS13 + N[0]*crRHS2 + crRHS3*(DN(0,0)*crRHS4 + DN(0,1)*crRHS5));
rRHS[1]+=-gauss_weight*(-DN(1,0)*crRHS15 - DN(1,1)*crRHS16 - N[1]*crRHS0 + N[1]*crRHS11 + N[1]*crRHS13 + N[1]*crRHS2 + crRHS3*(DN(1,0)*crRHS4 + DN(1,1)*crRHS5));
rRHS[2]+=-gauss_weight*(-DN(2,0)*crRHS15 - DN(2,1)*crRHS16 - N[2]*crRHS0 + N[2]*crRHS11 + N[2]*crRHS13 + N[2]*crRHS2 + crRHS3*(DN(2,0)*crRHS4 + DN(2,1)*crRHS5));

}

template <>
void TransportTopologyOptimizationElement<TransportTopologyOptimizationElementData<2,4,true>>::ComputeGaussPointRHSContribution(
    TransportTopologyOptimizationElementData<2,4,true>& rData,
    VectorType& rRHS)
{
    const array_1d<double,4> D = rData.Conductivity;
    const array_1d<double,4> k = rData.Decay;
    const array_1d<double,4> c = rData.ConvectionCoefficient;

    const double h = rData.ElementSize;
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,2>& DN = rData.DN_DX;

    // const double dyn_tau = rData.DynamicTau;
    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const array_1d<double,4>& t = rData.Temperature;
    const array_1d<double,4>& tn = rData.Temperature_OldStep1;
    const array_1d<double,4>& tnn = rData.Temperature_OldStep2;
    const BoundedMatrix<double,2,4>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,2,4> vconv = rData.ConvectiveVelocity - rData.MeshVelocity;
    const array_1d<double,4>& source = rData.SourceTerm;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // TRANSPORT ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*source[0] + N[1]*source[1] + N[2]*source[2] + N[3]*source[3];
const double crRHS1 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2] + N[3]*k[3];
const double crRHS2 = crRHS1*(N[0]*t[0] + N[1]*t[1] + N[2]*t[2] + N[3]*t[3]);
const double crRHS3 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2] + D[3]*N[3];
const double crRHS4 = DN(0,0)*t[0] + DN(1,0)*t[1] + DN(2,0)*t[2] + DN(3,0)*t[3];
const double crRHS5 = DN(0,1)*t[0] + DN(1,1)*t[1] + DN(2,1)*t[2] + DN(3,1)*t[3];
const double crRHS6 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crRHS7 = crRHS4*crRHS6;
const double crRHS8 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crRHS9 = crRHS5*crRHS8;
const double crRHS10 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3];
const double crRHS11 = crRHS10*(crRHS7 + crRHS9);
const double crRHS12 = 1.0*(crRHS0 - crRHS10*crRHS7 - crRHS10*crRHS9 - crRHS2)*1.0/(crRHS1*stab_c3 + crRHS3*stab_c1*1.0/(h*h) + stab_c2*fabs(crRHS10)*1.0/h*sqrt(crRHS6*crRHS6 + crRHS8*crRHS8));
const double crRHS13 = crRHS1*crRHS12;
const double crRHS14 = crRHS10*crRHS12;
const double crRHS15 = crRHS14*crRHS6;
const double crRHS16 = crRHS14*crRHS8;
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS15 - DN(0,1)*crRHS16 - N[0]*crRHS0 + N[0]*crRHS11 + N[0]*crRHS13 + N[0]*crRHS2 + crRHS3*(DN(0,0)*crRHS4 + DN(0,1)*crRHS5));
rRHS[1]+=-gauss_weight*(-DN(1,0)*crRHS15 - DN(1,1)*crRHS16 - N[1]*crRHS0 + N[1]*crRHS11 + N[1]*crRHS13 + N[1]*crRHS2 + crRHS3*(DN(1,0)*crRHS4 + DN(1,1)*crRHS5));
rRHS[2]+=-gauss_weight*(-DN(2,0)*crRHS15 - DN(2,1)*crRHS16 - N[2]*crRHS0 + N[2]*crRHS11 + N[2]*crRHS13 + N[2]*crRHS2 + crRHS3*(DN(2,0)*crRHS4 + DN(2,1)*crRHS5));
rRHS[3]+=-gauss_weight*(-DN(3,0)*crRHS15 - DN(3,1)*crRHS16 - N[3]*crRHS0 + N[3]*crRHS11 + N[3]*crRHS13 + N[3]*crRHS2 + crRHS3*(DN(3,0)*crRHS4 + DN(3,1)*crRHS5));

}

template <>
void TransportTopologyOptimizationElement<TransportTopologyOptimizationElementData<3,4,true>>::ComputeGaussPointRHSContribution(
    TransportTopologyOptimizationElementData<3,4,true>& rData,
    VectorType& rRHS)
{
    const array_1d<double,4> D = rData.Conductivity;
    const array_1d<double,4> k = rData.Decay;
    const array_1d<double,4> c = rData.ConvectionCoefficient;

    const double h = rData.ElementSize;
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // const double dyn_tau = rData.DynamicTau;
    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const array_1d<double,4>& t = rData.Temperature;
    const array_1d<double,4>& tn = rData.Temperature_OldStep1;
    const array_1d<double,4>& tnn = rData.Temperature_OldStep2;
    const BoundedMatrix<double,3,4>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,3,4> vconv = rData.ConvectiveVelocity - rData.MeshVelocity;
    const array_1d<double,4>& source = rData.SourceTerm;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // TRANSPORT ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*source[0] + N[1]*source[1] + N[2]*source[2] + N[3]*source[3];
const double crRHS1 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2] + N[3]*k[3];
const double crRHS2 = crRHS1*(N[0]*t[0] + N[1]*t[1] + N[2]*t[2] + N[3]*t[3]);
const double crRHS3 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2] + D[3]*N[3];
const double crRHS4 = DN(0,0)*t[0] + DN(1,0)*t[1] + DN(2,0)*t[2] + DN(3,0)*t[3];
const double crRHS5 = DN(0,1)*t[0] + DN(1,1)*t[1] + DN(2,1)*t[2] + DN(3,1)*t[3];
const double crRHS6 = DN(0,2)*t[0] + DN(1,2)*t[1] + DN(2,2)*t[2] + DN(3,2)*t[3];
const double crRHS7 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crRHS8 = crRHS4*crRHS7;
const double crRHS9 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crRHS10 = crRHS5*crRHS9;
const double crRHS11 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crRHS12 = crRHS11*crRHS6;
const double crRHS13 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3];
const double crRHS14 = crRHS13*(crRHS10 + crRHS12 + crRHS8);
const double crRHS15 = 1.0*(crRHS0 - crRHS10*crRHS13 - crRHS12*crRHS13 - crRHS13*crRHS8 - crRHS2)*1.0/(crRHS1*stab_c3 + crRHS3*stab_c1*1.0/(h*h) + stab_c2*fabs(crRHS13)*1.0/h*sqrt(crRHS11*crRHS11 + crRHS7*crRHS7 + crRHS9*crRHS9));
const double crRHS16 = crRHS1*crRHS15;
const double crRHS17 = crRHS13*crRHS15;
const double crRHS18 = crRHS17*crRHS7;
const double crRHS19 = crRHS17*crRHS9;
const double crRHS20 = crRHS11*crRHS17;
rRHS[0]+=gauss_weight*(DN(0,0)*crRHS18 + DN(0,1)*crRHS19 + DN(0,2)*crRHS20 + N[0]*crRHS0 - N[0]*crRHS14 - N[0]*crRHS16 - N[0]*crRHS2 - crRHS3*(DN(0,0)*crRHS4 + DN(0,1)*crRHS5 + DN(0,2)*crRHS6));
rRHS[1]+=gauss_weight*(DN(1,0)*crRHS18 + DN(1,1)*crRHS19 + DN(1,2)*crRHS20 + N[1]*crRHS0 - N[1]*crRHS14 - N[1]*crRHS16 - N[1]*crRHS2 - crRHS3*(DN(1,0)*crRHS4 + DN(1,1)*crRHS5 + DN(1,2)*crRHS6));
rRHS[2]+=gauss_weight*(DN(2,0)*crRHS18 + DN(2,1)*crRHS19 + DN(2,2)*crRHS20 + N[2]*crRHS0 - N[2]*crRHS14 - N[2]*crRHS16 - N[2]*crRHS2 - crRHS3*(DN(2,0)*crRHS4 + DN(2,1)*crRHS5 + DN(2,2)*crRHS6));
rRHS[3]+=gauss_weight*(DN(3,0)*crRHS18 + DN(3,1)*crRHS19 + DN(3,2)*crRHS20 + N[3]*crRHS0 - N[3]*crRHS14 - N[3]*crRHS16 - N[3]*crRHS2 - crRHS3*(DN(3,0)*crRHS4 + DN(3,1)*crRHS5 + DN(3,2)*crRHS6));

}

// BUILD ADJOINT NS SYSTEM
template <>
void TransportTopologyOptimizationElement< TransportTopologyOptimizationElementData<2,3,true>>::ComputeGaussPointLHSContributionAdjoint(
    TransportTopologyOptimizationElementData<2,3,true> & rData,
    MatrixType& rLHS)
{
    const array_1d<double,3> D = rData.Conductivity;
    const array_1d<double,3> k = rData.Decay;
    const array_1d<double,3> c = rData.ConvectionCoefficient;

    const double h = rData.ElementSize;
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // const double dyn_tau = rData.DynamicTau;
    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const BoundedMatrix<double,2,3> vconv_adj = rData.ConvectiveVelocity - rData.MeshVelocity;

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT TRANSPORT ELEMENTAL LHS MATRIX
    const double crLHS0 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2];
const double crLHS1 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2];
const double crLHS2 = N[0]*vconv_adj(0,0) + N[1]*vconv_adj(1,0) + N[2]*vconv_adj(2,0);
const double crLHS3 = DN(0,0)*crLHS2;
const double crLHS4 = N[0]*vconv_adj(0,1) + N[1]*vconv_adj(1,1) + N[2]*vconv_adj(2,1);
const double crLHS5 = DN(0,1)*crLHS4;
const double crLHS6 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2];
const double crLHS7 = crLHS6*(crLHS3 + crLHS5);
const double crLHS8 = N[0]*crLHS0;
const double crLHS9 = crLHS3*crLHS6;
const double crLHS10 = crLHS5*crLHS6;
const double crLHS11 = crLHS10 - crLHS8 + crLHS9;
const double crLHS12 = 1.0*1.0/(crLHS0*stab_c3 + crLHS1*stab_c1*1.0/(h*h) + stab_c2*fabs(crLHS6)*1.0/h*sqrt(crLHS2*crLHS2 + crLHS4*crLHS4));
const double crLHS13 = crLHS12*crLHS8;
const double crLHS14 = crLHS11*crLHS12;
const double crLHS15 = N[1]*crLHS0;
const double crLHS16 = DN(1,0)*crLHS2;
const double crLHS17 = crLHS16*crLHS6;
const double crLHS18 = DN(1,1)*crLHS4;
const double crLHS19 = crLHS18*crLHS6;
const double crLHS20 = -crLHS15 + crLHS17 + crLHS19;
const double crLHS21 = crLHS12*crLHS20;
const double crLHS22 = N[1]*crLHS8 + crLHS1*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1));
const double crLHS23 = N[2]*crLHS0;
const double crLHS24 = DN(2,0)*crLHS2;
const double crLHS25 = crLHS24*crLHS6;
const double crLHS26 = DN(2,1)*crLHS4;
const double crLHS27 = crLHS26*crLHS6;
const double crLHS28 = -crLHS23 + crLHS25 + crLHS27;
const double crLHS29 = crLHS12*crLHS28;
const double crLHS30 = N[2]*crLHS8 + crLHS1*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1));
const double crLHS31 = crLHS6*(crLHS16 + crLHS18);
const double crLHS32 = crLHS12*crLHS15;
const double crLHS33 = N[2]*crLHS15 + crLHS1*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1));
const double crLHS34 = crLHS6*(crLHS24 + crLHS26);
const double crLHS35 = crLHS12*crLHS23;
rLHS(0,0)+=gauss_weight*(N[0]*crLHS7 + crLHS0*(N[0]*N[0]) + crLHS1*(DN(0,0)*DN(0,0) + DN(0,1)*DN(0,1)) + crLHS10*crLHS14 + crLHS11*crLHS13 + crLHS14*crLHS9);
rLHS(0,1)+=gauss_weight*(N[1]*crLHS7 + crLHS10*crLHS21 + crLHS13*crLHS20 + crLHS21*crLHS9 + crLHS22);
rLHS(0,2)+=gauss_weight*(N[2]*crLHS7 + crLHS10*crLHS29 + crLHS13*crLHS28 + crLHS29*crLHS9 + crLHS30);
rLHS(1,0)+=gauss_weight*(N[0]*crLHS31 + crLHS11*crLHS32 + crLHS14*crLHS17 + crLHS14*crLHS19 + crLHS22);
rLHS(1,1)+=gauss_weight*(N[1]*crLHS31 + crLHS0*(N[1]*N[1]) + crLHS1*(DN(1,0)*DN(1,0) + DN(1,1)*DN(1,1)) + crLHS17*crLHS21 + crLHS19*crLHS21 + crLHS20*crLHS32);
rLHS(1,2)+=gauss_weight*(N[2]*crLHS31 + crLHS17*crLHS29 + crLHS19*crLHS29 + crLHS28*crLHS32 + crLHS33);
rLHS(2,0)+=gauss_weight*(N[0]*crLHS34 + crLHS11*crLHS35 + crLHS14*crLHS25 + crLHS14*crLHS27 + crLHS30);
rLHS(2,1)+=gauss_weight*(N[1]*crLHS34 + crLHS20*crLHS35 + crLHS21*crLHS25 + crLHS21*crLHS27 + crLHS33);
rLHS(2,2)+=gauss_weight*(N[2]*crLHS34 + crLHS0*(N[2]*N[2]) + crLHS1*(DN(2,0)*DN(2,0) + DN(2,1)*DN(2,1)) + crLHS25*crLHS29 + crLHS27*crLHS29 + crLHS28*crLHS35);

}

template <>
void TransportTopologyOptimizationElement< TransportTopologyOptimizationElementData<2,4,true>>::ComputeGaussPointLHSContributionAdjoint(
    TransportTopologyOptimizationElementData<2,4,true> & rData,
    MatrixType& rLHS)
{
    const array_1d<double,4> D = rData.Conductivity;
    const array_1d<double,4> k = rData.Decay;
    const array_1d<double,4> c = rData.ConvectionCoefficient;

    const double h = rData.ElementSize;
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,2>& DN = rData.DN_DX;

    // const double dyn_tau = rData.DynamicTau;
    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const BoundedMatrix<double,2,4> vconv_adj = rData.ConvectiveVelocity - rData.MeshVelocity;

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT TRANSPORT ELEMENTAL LHS MATRIX
    const double crLHS0 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2] + N[3]*k[3];
const double crLHS1 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2] + D[3]*N[3];
const double crLHS2 = N[0]*vconv_adj(0,0) + N[1]*vconv_adj(1,0) + N[2]*vconv_adj(2,0) + N[3]*vconv_adj(3,0);
const double crLHS3 = DN(0,0)*crLHS2;
const double crLHS4 = N[0]*vconv_adj(0,1) + N[1]*vconv_adj(1,1) + N[2]*vconv_adj(2,1) + N[3]*vconv_adj(3,1);
const double crLHS5 = DN(0,1)*crLHS4;
const double crLHS6 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3];
const double crLHS7 = crLHS6*(crLHS3 + crLHS5);
const double crLHS8 = N[0]*crLHS0;
const double crLHS9 = crLHS3*crLHS6;
const double crLHS10 = crLHS5*crLHS6;
const double crLHS11 = crLHS10 - crLHS8 + crLHS9;
const double crLHS12 = 1.0*1.0/(crLHS0*stab_c3 + crLHS1*stab_c1*1.0/(h*h) + stab_c2*fabs(crLHS6)*1.0/h*sqrt(crLHS2*crLHS2 + crLHS4*crLHS4));
const double crLHS13 = crLHS12*crLHS8;
const double crLHS14 = crLHS11*crLHS12;
const double crLHS15 = N[1]*crLHS0;
const double crLHS16 = DN(1,0)*crLHS2;
const double crLHS17 = crLHS16*crLHS6;
const double crLHS18 = DN(1,1)*crLHS4;
const double crLHS19 = crLHS18*crLHS6;
const double crLHS20 = -crLHS15 + crLHS17 + crLHS19;
const double crLHS21 = crLHS12*crLHS20;
const double crLHS22 = N[1]*crLHS8 + crLHS1*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1));
const double crLHS23 = N[2]*crLHS0;
const double crLHS24 = DN(2,0)*crLHS2;
const double crLHS25 = crLHS24*crLHS6;
const double crLHS26 = DN(2,1)*crLHS4;
const double crLHS27 = crLHS26*crLHS6;
const double crLHS28 = -crLHS23 + crLHS25 + crLHS27;
const double crLHS29 = crLHS12*crLHS28;
const double crLHS30 = N[2]*crLHS8 + crLHS1*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1));
const double crLHS31 = N[3]*crLHS0;
const double crLHS32 = DN(3,0)*crLHS2;
const double crLHS33 = crLHS32*crLHS6;
const double crLHS34 = DN(3,1)*crLHS4;
const double crLHS35 = crLHS34*crLHS6;
const double crLHS36 = -crLHS31 + crLHS33 + crLHS35;
const double crLHS37 = crLHS12*crLHS36;
const double crLHS38 = N[3]*crLHS8 + crLHS1*(DN(0,0)*DN(3,0) + DN(0,1)*DN(3,1));
const double crLHS39 = crLHS6*(crLHS16 + crLHS18);
const double crLHS40 = crLHS12*crLHS15;
const double crLHS41 = N[2]*crLHS15 + crLHS1*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1));
const double crLHS42 = N[3]*crLHS15 + crLHS1*(DN(1,0)*DN(3,0) + DN(1,1)*DN(3,1));
const double crLHS43 = crLHS6*(crLHS24 + crLHS26);
const double crLHS44 = crLHS12*crLHS23;
const double crLHS45 = N[3]*crLHS23 + crLHS1*(DN(2,0)*DN(3,0) + DN(2,1)*DN(3,1));
const double crLHS46 = crLHS6*(crLHS32 + crLHS34);
const double crLHS47 = crLHS12*crLHS31;
rLHS(0,0)+=gauss_weight*(N[0]*crLHS7 + crLHS0*(N[0]*N[0]) + crLHS1*(DN(0,0)*DN(0,0) + DN(0,1)*DN(0,1)) + crLHS10*crLHS14 + crLHS11*crLHS13 + crLHS14*crLHS9);
rLHS(0,1)+=gauss_weight*(N[1]*crLHS7 + crLHS10*crLHS21 + crLHS13*crLHS20 + crLHS21*crLHS9 + crLHS22);
rLHS(0,2)+=gauss_weight*(N[2]*crLHS7 + crLHS10*crLHS29 + crLHS13*crLHS28 + crLHS29*crLHS9 + crLHS30);
rLHS(0,3)+=gauss_weight*(N[3]*crLHS7 + crLHS10*crLHS37 + crLHS13*crLHS36 + crLHS37*crLHS9 + crLHS38);
rLHS(1,0)+=gauss_weight*(N[0]*crLHS39 + crLHS11*crLHS40 + crLHS14*crLHS17 + crLHS14*crLHS19 + crLHS22);
rLHS(1,1)+=gauss_weight*(N[1]*crLHS39 + crLHS0*(N[1]*N[1]) + crLHS1*(DN(1,0)*DN(1,0) + DN(1,1)*DN(1,1)) + crLHS17*crLHS21 + crLHS19*crLHS21 + crLHS20*crLHS40);
rLHS(1,2)+=gauss_weight*(N[2]*crLHS39 + crLHS17*crLHS29 + crLHS19*crLHS29 + crLHS28*crLHS40 + crLHS41);
rLHS(1,3)+=gauss_weight*(N[3]*crLHS39 + crLHS17*crLHS37 + crLHS19*crLHS37 + crLHS36*crLHS40 + crLHS42);
rLHS(2,0)+=gauss_weight*(N[0]*crLHS43 + crLHS11*crLHS44 + crLHS14*crLHS25 + crLHS14*crLHS27 + crLHS30);
rLHS(2,1)+=gauss_weight*(N[1]*crLHS43 + crLHS20*crLHS44 + crLHS21*crLHS25 + crLHS21*crLHS27 + crLHS41);
rLHS(2,2)+=gauss_weight*(N[2]*crLHS43 + crLHS0*(N[2]*N[2]) + crLHS1*(DN(2,0)*DN(2,0) + DN(2,1)*DN(2,1)) + crLHS25*crLHS29 + crLHS27*crLHS29 + crLHS28*crLHS44);
rLHS(2,3)+=gauss_weight*(N[3]*crLHS43 + crLHS25*crLHS37 + crLHS27*crLHS37 + crLHS36*crLHS44 + crLHS45);
rLHS(3,0)+=gauss_weight*(N[0]*crLHS46 + crLHS11*crLHS47 + crLHS14*crLHS33 + crLHS14*crLHS35 + crLHS38);
rLHS(3,1)+=gauss_weight*(N[1]*crLHS46 + crLHS20*crLHS47 + crLHS21*crLHS33 + crLHS21*crLHS35 + crLHS42);
rLHS(3,2)+=gauss_weight*(N[2]*crLHS46 + crLHS28*crLHS47 + crLHS29*crLHS33 + crLHS29*crLHS35 + crLHS45);
rLHS(3,3)+=gauss_weight*(N[3]*crLHS46 + crLHS0*(N[3]*N[3]) + crLHS1*(DN(3,0)*DN(3,0) + DN(3,1)*DN(3,1)) + crLHS33*crLHS37 + crLHS35*crLHS37 + crLHS36*crLHS47);

}

template <>
void TransportTopologyOptimizationElement<TransportTopologyOptimizationElementData<3,4,true>>::ComputeGaussPointLHSContributionAdjoint(
    TransportTopologyOptimizationElementData<3,4,true>& rData,
    MatrixType& rLHS)
{
    const array_1d<double,4> D = rData.Conductivity;
    const array_1d<double,4> k = rData.Decay;
    const array_1d<double,4> c = rData.ConvectionCoefficient;

    const double h = rData.ElementSize;
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,3>& DN = rData.DN_DX;

    // const double dyn_tau = rData.DynamicTau;
    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const BoundedMatrix<double,3,4> vconv_adj = rData.ConvectiveVelocity - rData.MeshVelocity;

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT TRANSPORT ELEMENTAL LHS MATRIX
    const double crLHS0 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2] + N[3]*k[3];
const double crLHS1 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2] + D[3]*N[3];
const double crLHS2 = N[0]*vconv_adj(0,0) + N[1]*vconv_adj(1,0) + N[2]*vconv_adj(2,0) + N[3]*vconv_adj(3,0);
const double crLHS3 = DN(0,0)*crLHS2;
const double crLHS4 = N[0]*vconv_adj(0,1) + N[1]*vconv_adj(1,1) + N[2]*vconv_adj(2,1) + N[3]*vconv_adj(3,1);
const double crLHS5 = DN(0,1)*crLHS4;
const double crLHS6 = N[0]*vconv_adj(0,2) + N[1]*vconv_adj(1,2) + N[2]*vconv_adj(2,2) + N[3]*vconv_adj(3,2);
const double crLHS7 = DN(0,2)*crLHS6;
const double crLHS8 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3];
const double crLHS9 = crLHS8*(crLHS3 + crLHS5 + crLHS7);
const double crLHS10 = N[0]*crLHS0;
const double crLHS11 = crLHS3*crLHS8;
const double crLHS12 = crLHS5*crLHS8;
const double crLHS13 = crLHS7*crLHS8;
const double crLHS14 = -crLHS10 + crLHS11 + crLHS12 + crLHS13;
const double crLHS15 = 1.0*1.0/(crLHS0*stab_c3 + crLHS1*stab_c1*1.0/(h*h) + stab_c2*fabs(crLHS8)*1.0/h*sqrt(crLHS2*crLHS2 + crLHS4*crLHS4 + crLHS6*crLHS6));
const double crLHS16 = crLHS10*crLHS15;
const double crLHS17 = crLHS14*crLHS15;
const double crLHS18 = N[1]*crLHS0;
const double crLHS19 = DN(1,0)*crLHS2;
const double crLHS20 = crLHS19*crLHS8;
const double crLHS21 = DN(1,1)*crLHS4;
const double crLHS22 = crLHS21*crLHS8;
const double crLHS23 = DN(1,2)*crLHS6;
const double crLHS24 = crLHS23*crLHS8;
const double crLHS25 = -crLHS18 + crLHS20 + crLHS22 + crLHS24;
const double crLHS26 = crLHS15*crLHS25;
const double crLHS27 = N[1]*crLHS10 + crLHS1*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1) + DN(0,2)*DN(1,2));
const double crLHS28 = N[2]*crLHS0;
const double crLHS29 = DN(2,0)*crLHS2;
const double crLHS30 = crLHS29*crLHS8;
const double crLHS31 = DN(2,1)*crLHS4;
const double crLHS32 = crLHS31*crLHS8;
const double crLHS33 = DN(2,2)*crLHS6;
const double crLHS34 = crLHS33*crLHS8;
const double crLHS35 = -crLHS28 + crLHS30 + crLHS32 + crLHS34;
const double crLHS36 = crLHS15*crLHS35;
const double crLHS37 = N[2]*crLHS10 + crLHS1*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1) + DN(0,2)*DN(2,2));
const double crLHS38 = N[3]*crLHS0;
const double crLHS39 = DN(3,0)*crLHS2;
const double crLHS40 = crLHS39*crLHS8;
const double crLHS41 = DN(3,1)*crLHS4;
const double crLHS42 = crLHS41*crLHS8;
const double crLHS43 = DN(3,2)*crLHS6;
const double crLHS44 = crLHS43*crLHS8;
const double crLHS45 = -crLHS38 + crLHS40 + crLHS42 + crLHS44;
const double crLHS46 = crLHS15*crLHS45;
const double crLHS47 = N[3]*crLHS10 + crLHS1*(DN(0,0)*DN(3,0) + DN(0,1)*DN(3,1) + DN(0,2)*DN(3,2));
const double crLHS48 = crLHS8*(crLHS19 + crLHS21 + crLHS23);
const double crLHS49 = crLHS15*crLHS18;
const double crLHS50 = N[2]*crLHS18 + crLHS1*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1) + DN(1,2)*DN(2,2));
const double crLHS51 = N[3]*crLHS18 + crLHS1*(DN(1,0)*DN(3,0) + DN(1,1)*DN(3,1) + DN(1,2)*DN(3,2));
const double crLHS52 = crLHS8*(crLHS29 + crLHS31 + crLHS33);
const double crLHS53 = crLHS15*crLHS28;
const double crLHS54 = N[3]*crLHS28 + crLHS1*(DN(2,0)*DN(3,0) + DN(2,1)*DN(3,1) + DN(2,2)*DN(3,2));
const double crLHS55 = crLHS8*(crLHS39 + crLHS41 + crLHS43);
const double crLHS56 = crLHS15*crLHS38;
rLHS(0,0)+=gauss_weight*(N[0]*crLHS9 + crLHS0*(N[0]*N[0]) + crLHS1*(DN(0,0)*DN(0,0) + DN(0,1)*DN(0,1) + DN(0,2)*DN(0,2)) + crLHS11*crLHS17 + crLHS12*crLHS17 + crLHS13*crLHS17 + crLHS14*crLHS16);
rLHS(0,1)+=gauss_weight*(N[1]*crLHS9 + crLHS11*crLHS26 + crLHS12*crLHS26 + crLHS13*crLHS26 + crLHS16*crLHS25 + crLHS27);
rLHS(0,2)+=gauss_weight*(N[2]*crLHS9 + crLHS11*crLHS36 + crLHS12*crLHS36 + crLHS13*crLHS36 + crLHS16*crLHS35 + crLHS37);
rLHS(0,3)+=gauss_weight*(N[3]*crLHS9 + crLHS11*crLHS46 + crLHS12*crLHS46 + crLHS13*crLHS46 + crLHS16*crLHS45 + crLHS47);
rLHS(1,0)+=gauss_weight*(N[0]*crLHS48 + crLHS14*crLHS49 + crLHS17*crLHS20 + crLHS17*crLHS22 + crLHS17*crLHS24 + crLHS27);
rLHS(1,1)+=gauss_weight*(N[1]*crLHS48 + crLHS0*(N[1]*N[1]) + crLHS1*(DN(1,0)*DN(1,0) + DN(1,1)*DN(1,1) + DN(1,2)*DN(1,2)) + crLHS20*crLHS26 + crLHS22*crLHS26 + crLHS24*crLHS26 + crLHS25*crLHS49);
rLHS(1,2)+=gauss_weight*(N[2]*crLHS48 + crLHS20*crLHS36 + crLHS22*crLHS36 + crLHS24*crLHS36 + crLHS35*crLHS49 + crLHS50);
rLHS(1,3)+=gauss_weight*(N[3]*crLHS48 + crLHS20*crLHS46 + crLHS22*crLHS46 + crLHS24*crLHS46 + crLHS45*crLHS49 + crLHS51);
rLHS(2,0)+=gauss_weight*(N[0]*crLHS52 + crLHS14*crLHS53 + crLHS17*crLHS30 + crLHS17*crLHS32 + crLHS17*crLHS34 + crLHS37);
rLHS(2,1)+=gauss_weight*(N[1]*crLHS52 + crLHS25*crLHS53 + crLHS26*crLHS30 + crLHS26*crLHS32 + crLHS26*crLHS34 + crLHS50);
rLHS(2,2)+=gauss_weight*(N[2]*crLHS52 + crLHS0*(N[2]*N[2]) + crLHS1*(DN(2,0)*DN(2,0) + DN(2,1)*DN(2,1) + DN(2,2)*DN(2,2)) + crLHS30*crLHS36 + crLHS32*crLHS36 + crLHS34*crLHS36 + crLHS35*crLHS53);
rLHS(2,3)+=gauss_weight*(N[3]*crLHS52 + crLHS30*crLHS46 + crLHS32*crLHS46 + crLHS34*crLHS46 + crLHS45*crLHS53 + crLHS54);
rLHS(3,0)+=gauss_weight*(N[0]*crLHS55 + crLHS14*crLHS56 + crLHS17*crLHS40 + crLHS17*crLHS42 + crLHS17*crLHS44 + crLHS47);
rLHS(3,1)+=gauss_weight*(N[1]*crLHS55 + crLHS25*crLHS56 + crLHS26*crLHS40 + crLHS26*crLHS42 + crLHS26*crLHS44 + crLHS51);
rLHS(3,2)+=gauss_weight*(N[2]*crLHS55 + crLHS35*crLHS56 + crLHS36*crLHS40 + crLHS36*crLHS42 + crLHS36*crLHS44 + crLHS54);
rLHS(3,3)+=gauss_weight*(N[3]*crLHS55 + crLHS0*(N[3]*N[3]) + crLHS1*(DN(3,0)*DN(3,0) + DN(3,1)*DN(3,1) + DN(3,2)*DN(3,2)) + crLHS40*crLHS46 + crLHS42*crLHS46 + crLHS44*crLHS46 + crLHS45*crLHS56);

}

template <>
void TransportTopologyOptimizationElement<TransportTopologyOptimizationElementData<2,3,true>>::ComputeGaussPointRHSContributionAdjoint(
    TransportTopologyOptimizationElementData<2,3,true>& rData,
    VectorType& rRHS)
{
    const array_1d<double,3> D = rData.Conductivity;
    const array_1d<double,3> k = rData.Decay;
    const array_1d<double,3> c = rData.ConvectionCoefficient;

    const double h = rData.ElementSize;
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // const double dyn_tau = rData.DynamicTau;
    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    Vector functional_weights = rData.Functional_Weights; //  functional terms weights

    const array_1d<double,3>& t_adj = rData.Temperature_adj;
    const array_1d<double,3>& tn_adj = rData.Temperature_adj_OldStep1;
    const array_1d<double,3>& tnn_adj = rData.Temperature_adj_OldStep2;
    const BoundedMatrix<double,2,3>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,2,3> vconv_adj = rData.ConvectiveVelocity - rData.MeshVelocity;
    const array_1d<double,3>& source_adj = rData.SourceTerm_adj;
    const array_1d<double,3>& opt_t = rData.Optimization_Temperature;
    const array_1d<double,3>& t = rData.Temperature;
    const BoundedMatrix<double,2,3> vconv = rData.ConvectiveVelocity - rData.MeshVelocity;
    const array_1d<double,3>& source = rData.SourceTerm;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT TRANSPORT ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*source_adj[0] + N[1]*source_adj[1] + N[2]*source_adj[2];
const double crRHS1 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2];
const double crRHS2 = crRHS1*functional_weights[9];
const double crRHS3 = 2.0*functional_weights[4]*(N[0]*opt_t[0] + N[1]*opt_t[1] + N[2]*opt_t[2]);
const double crRHS4 = functional_weights[8]*(N[0]*source[0] + N[1]*source[1] + N[2]*source[2]);
const double crRHS5 = N[0]*t_adj[0] + N[1]*t_adj[1] + N[2]*t_adj[2];
const double crRHS6 = crRHS1*crRHS5;
const double crRHS7 = N[0]*t[0] + N[1]*t[1] + N[2]*t[2];
const double crRHS8 = 2.0*crRHS1*crRHS7*functional_weights[7];
const double crRHS9 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2];
const double crRHS10 = DN(0,0)*t_adj[0] + DN(1,0)*t_adj[1] + DN(2,0)*t_adj[2];
const double crRHS11 = DN(0,1)*t_adj[0] + DN(1,1)*t_adj[1] + DN(2,1)*t_adj[2];
const double crRHS12 = DN(0,0)*t[0] + DN(1,0)*t[1] + DN(2,0)*t[2];
const double crRHS13 = DN(0,1)*t[0] + DN(1,1)*t[1] + DN(2,1)*t[2];
const double crRHS14 = 2.0*crRHS9*functional_weights[5];
const double crRHS15 = N[0]*vconv_adj(0,0) + N[1]*vconv_adj(1,0) + N[2]*vconv_adj(2,0);
const double crRHS16 = DN(0,0)*crRHS15;
const double crRHS17 = N[0]*vconv_adj(0,1) + N[1]*vconv_adj(1,1) + N[2]*vconv_adj(2,1);
const double crRHS18 = DN(0,1)*crRHS17;
const double crRHS19 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2];
const double crRHS20 = crRHS19*crRHS5;
const double crRHS21 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crRHS22 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crRHS23 = crRHS12*crRHS21 + crRHS13*crRHS22;
const double crRHS24 = crRHS19*functional_weights[6];
const double crRHS25 = 1.0*(crRHS0 + crRHS10*crRHS15*crRHS19 + crRHS11*crRHS17*crRHS19 - crRHS2 - crRHS3 + crRHS4 - crRHS6 + crRHS7*functional_weights[6]*(crRHS21*(DN(0,0)*c[0] + DN(1,0)*c[1] + DN(2,0)*c[2]) + crRHS22*(DN(0,1)*c[0] + DN(1,1)*c[1] + DN(2,1)*c[2])) - crRHS8)*1.0/(crRHS1*stab_c3 + crRHS9*stab_c1*1.0/(h*h) + stab_c2*fabs(crRHS19)*1.0/h*sqrt(crRHS15*crRHS15 + crRHS17*crRHS17));
const double crRHS26 = crRHS1*crRHS25;
const double crRHS27 = crRHS19*crRHS25;
const double crRHS28 = DN(1,0)*crRHS15;
const double crRHS29 = DN(1,1)*crRHS17;
const double crRHS30 = DN(2,0)*crRHS15;
const double crRHS31 = DN(2,1)*crRHS17;
rRHS[0]+=-gauss_weight*(-N[0]*crRHS0 + N[0]*crRHS2 + N[0]*crRHS26 + N[0]*crRHS3 - N[0]*crRHS4 + N[0]*crRHS6 + N[0]*crRHS8 + crRHS14*(DN(0,0)*crRHS12 + DN(0,1)*crRHS13) + crRHS16*crRHS27 + crRHS18*crRHS27 + crRHS20*(crRHS16 + crRHS18) + crRHS24*(N[0]*crRHS23 + crRHS7*(DN(0,0)*crRHS21 + DN(0,1)*crRHS22)) + crRHS9*(DN(0,0)*crRHS10 + DN(0,1)*crRHS11));
rRHS[1]+=-gauss_weight*(-N[1]*crRHS0 + N[1]*crRHS2 + N[1]*crRHS26 + N[1]*crRHS3 - N[1]*crRHS4 + N[1]*crRHS6 + N[1]*crRHS8 + crRHS14*(DN(1,0)*crRHS12 + DN(1,1)*crRHS13) + crRHS20*(crRHS28 + crRHS29) + crRHS24*(N[1]*crRHS23 + crRHS7*(DN(1,0)*crRHS21 + DN(1,1)*crRHS22)) + crRHS27*crRHS28 + crRHS27*crRHS29 + crRHS9*(DN(1,0)*crRHS10 + DN(1,1)*crRHS11));
rRHS[2]+=-gauss_weight*(-N[2]*crRHS0 + N[2]*crRHS2 + N[2]*crRHS26 + N[2]*crRHS3 - N[2]*crRHS4 + N[2]*crRHS6 + N[2]*crRHS8 + crRHS14*(DN(2,0)*crRHS12 + DN(2,1)*crRHS13) + crRHS20*(crRHS30 + crRHS31) + crRHS24*(N[2]*crRHS23 + crRHS7*(DN(2,0)*crRHS21 + DN(2,1)*crRHS22)) + crRHS27*crRHS30 + crRHS27*crRHS31 + crRHS9*(DN(2,0)*crRHS10 + DN(2,1)*crRHS11));

}

template <>
void TransportTopologyOptimizationElement<TransportTopologyOptimizationElementData<2,4,true>>::ComputeGaussPointRHSContributionAdjoint(
    TransportTopologyOptimizationElementData<2,4,true>& rData,
    VectorType& rRHS)
{
    const array_1d<double,4> D = rData.Conductivity;
    const array_1d<double,4> k = rData.Decay;
    const array_1d<double,4> c = rData.ConvectionCoefficient;

    const double h = rData.ElementSize;
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,2>& DN = rData.DN_DX;

    // const double dyn_tau = rData.DynamicTau;
    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    Vector functional_weights = rData.Functional_Weights; //  functional terms weights

    const array_1d<double,4>& t_adj = rData.Temperature_adj;
    const array_1d<double,4>& tn_adj = rData.Temperature_adj_OldStep1;
    const array_1d<double,4>& tnn_adj = rData.Temperature_adj_OldStep2;
    const BoundedMatrix<double,2,4>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,2,4> vconv_adj = rData.ConvectiveVelocity - rData.MeshVelocity;
    const array_1d<double,4>& source_adj = rData.SourceTerm_adj;
    const array_1d<double,4>& opt_t = rData.Optimization_Temperature;
    const array_1d<double,4>& t = rData.Temperature;
    const BoundedMatrix<double,2,4> vconv = rData.ConvectiveVelocity - rData.MeshVelocity;
    const array_1d<double,4>& source = rData.SourceTerm;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT TRANSPORT ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*source_adj[0] + N[1]*source_adj[1] + N[2]*source_adj[2] + N[3]*source_adj[3];
const double crRHS1 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2] + N[3]*k[3];
const double crRHS2 = crRHS1*functional_weights[9];
const double crRHS3 = 2.0*functional_weights[4]*(N[0]*opt_t[0] + N[1]*opt_t[1] + N[2]*opt_t[2] + N[3]*opt_t[3]);
const double crRHS4 = functional_weights[8]*(N[0]*source[0] + N[1]*source[1] + N[2]*source[2] + N[3]*source[3]);
const double crRHS5 = N[0]*t_adj[0] + N[1]*t_adj[1] + N[2]*t_adj[2] + N[3]*t_adj[3];
const double crRHS6 = crRHS1*crRHS5;
const double crRHS7 = N[0]*t[0] + N[1]*t[1] + N[2]*t[2] + N[3]*t[3];
const double crRHS8 = 2.0*crRHS1*crRHS7*functional_weights[7];
const double crRHS9 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2] + D[3]*N[3];
const double crRHS10 = DN(0,0)*t_adj[0] + DN(1,0)*t_adj[1] + DN(2,0)*t_adj[2] + DN(3,0)*t_adj[3];
const double crRHS11 = DN(0,1)*t_adj[0] + DN(1,1)*t_adj[1] + DN(2,1)*t_adj[2] + DN(3,1)*t_adj[3];
const double crRHS12 = DN(0,0)*t[0] + DN(1,0)*t[1] + DN(2,0)*t[2] + DN(3,0)*t[3];
const double crRHS13 = DN(0,1)*t[0] + DN(1,1)*t[1] + DN(2,1)*t[2] + DN(3,1)*t[3];
const double crRHS14 = 2.0*crRHS9*functional_weights[5];
const double crRHS15 = N[0]*vconv_adj(0,0) + N[1]*vconv_adj(1,0) + N[2]*vconv_adj(2,0) + N[3]*vconv_adj(3,0);
const double crRHS16 = DN(0,0)*crRHS15;
const double crRHS17 = N[0]*vconv_adj(0,1) + N[1]*vconv_adj(1,1) + N[2]*vconv_adj(2,1) + N[3]*vconv_adj(3,1);
const double crRHS18 = DN(0,1)*crRHS17;
const double crRHS19 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3];
const double crRHS20 = crRHS19*crRHS5;
const double crRHS21 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crRHS22 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crRHS23 = crRHS12*crRHS21 + crRHS13*crRHS22;
const double crRHS24 = crRHS19*functional_weights[6];
const double crRHS25 = 1.0*(crRHS0 + crRHS10*crRHS15*crRHS19 + crRHS11*crRHS17*crRHS19 - crRHS2 - crRHS3 + crRHS4 - crRHS6 + crRHS7*functional_weights[6]*(crRHS21*(DN(0,0)*c[0] + DN(1,0)*c[1] + DN(2,0)*c[2] + DN(3,0)*c[3]) + crRHS22*(DN(0,1)*c[0] + DN(1,1)*c[1] + DN(2,1)*c[2] + DN(3,1)*c[3])) - crRHS8)*1.0/(crRHS1*stab_c3 + crRHS9*stab_c1*1.0/(h*h) + stab_c2*fabs(crRHS19)*1.0/h*sqrt(crRHS15*crRHS15 + crRHS17*crRHS17));
const double crRHS26 = crRHS1*crRHS25;
const double crRHS27 = crRHS19*crRHS25;
const double crRHS28 = DN(1,0)*crRHS15;
const double crRHS29 = DN(1,1)*crRHS17;
const double crRHS30 = DN(2,0)*crRHS15;
const double crRHS31 = DN(2,1)*crRHS17;
const double crRHS32 = DN(3,0)*crRHS15;
const double crRHS33 = DN(3,1)*crRHS17;
rRHS[0]+=-gauss_weight*(-N[0]*crRHS0 + N[0]*crRHS2 + N[0]*crRHS26 + N[0]*crRHS3 - N[0]*crRHS4 + N[0]*crRHS6 + N[0]*crRHS8 + crRHS14*(DN(0,0)*crRHS12 + DN(0,1)*crRHS13) + crRHS16*crRHS27 + crRHS18*crRHS27 + crRHS20*(crRHS16 + crRHS18) + crRHS24*(N[0]*crRHS23 + crRHS7*(DN(0,0)*crRHS21 + DN(0,1)*crRHS22)) + crRHS9*(DN(0,0)*crRHS10 + DN(0,1)*crRHS11));
rRHS[1]+=-gauss_weight*(-N[1]*crRHS0 + N[1]*crRHS2 + N[1]*crRHS26 + N[1]*crRHS3 - N[1]*crRHS4 + N[1]*crRHS6 + N[1]*crRHS8 + crRHS14*(DN(1,0)*crRHS12 + DN(1,1)*crRHS13) + crRHS20*(crRHS28 + crRHS29) + crRHS24*(N[1]*crRHS23 + crRHS7*(DN(1,0)*crRHS21 + DN(1,1)*crRHS22)) + crRHS27*crRHS28 + crRHS27*crRHS29 + crRHS9*(DN(1,0)*crRHS10 + DN(1,1)*crRHS11));
rRHS[2]+=-gauss_weight*(-N[2]*crRHS0 + N[2]*crRHS2 + N[2]*crRHS26 + N[2]*crRHS3 - N[2]*crRHS4 + N[2]*crRHS6 + N[2]*crRHS8 + crRHS14*(DN(2,0)*crRHS12 + DN(2,1)*crRHS13) + crRHS20*(crRHS30 + crRHS31) + crRHS24*(N[2]*crRHS23 + crRHS7*(DN(2,0)*crRHS21 + DN(2,1)*crRHS22)) + crRHS27*crRHS30 + crRHS27*crRHS31 + crRHS9*(DN(2,0)*crRHS10 + DN(2,1)*crRHS11));
rRHS[3]+=-gauss_weight*(-N[3]*crRHS0 + N[3]*crRHS2 + N[3]*crRHS26 + N[3]*crRHS3 - N[3]*crRHS4 + N[3]*crRHS6 + N[3]*crRHS8 + crRHS14*(DN(3,0)*crRHS12 + DN(3,1)*crRHS13) + crRHS20*(crRHS32 + crRHS33) + crRHS24*(N[3]*crRHS23 + crRHS7*(DN(3,0)*crRHS21 + DN(3,1)*crRHS22)) + crRHS27*crRHS32 + crRHS27*crRHS33 + crRHS9*(DN(3,0)*crRHS10 + DN(3,1)*crRHS11));

}

template <>
void TransportTopologyOptimizationElement<TransportTopologyOptimizationElementData<3,4,true>>::ComputeGaussPointRHSContributionAdjoint(
    TransportTopologyOptimizationElementData<3,4,true>& rData,
    VectorType& rRHS)
{
    const array_1d<double,4> D = rData.Conductivity;
    const array_1d<double,4> k = rData.Decay;
    const array_1d<double,4> c = rData.ConvectionCoefficient;

    const double h = rData.ElementSize;
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // const double dyn_tau = rData.DynamicTau;
    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    Vector functional_weights = rData.Functional_Weights; //  functional terms weights

    const array_1d<double,4>& t_adj = rData.Temperature_adj;
    const array_1d<double,4>& tn_adj = rData.Temperature_adj_OldStep1;
    const array_1d<double,4>& tnn_adj = rData.Temperature_adj_OldStep2;
    const BoundedMatrix<double,3,4>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,3,4> vconv_adj = rData.ConvectiveVelocity - rData.MeshVelocity;
    const array_1d<double,4>& source_adj = rData.SourceTerm_adj;
    const array_1d<double,4>& opt_t = rData.Optimization_Temperature;
    const array_1d<double,4>& t = rData.Temperature;
    const BoundedMatrix<double,3,4> vconv = rData.ConvectiveVelocity - rData.MeshVelocity;
    const array_1d<double,4>& source = rData.SourceTerm;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT TRANSPORT ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*source_adj[0] + N[1]*source_adj[1] + N[2]*source_adj[2] + N[3]*source_adj[3];
const double crRHS1 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2] + N[3]*k[3];
const double crRHS2 = crRHS1*functional_weights[9];
const double crRHS3 = 2.0*functional_weights[4]*(N[0]*opt_t[0] + N[1]*opt_t[1] + N[2]*opt_t[2] + N[3]*opt_t[3]);
const double crRHS4 = functional_weights[8]*(N[0]*source[0] + N[1]*source[1] + N[2]*source[2] + N[3]*source[3]);
const double crRHS5 = N[0]*t_adj[0] + N[1]*t_adj[1] + N[2]*t_adj[2] + N[3]*t_adj[3];
const double crRHS6 = crRHS1*crRHS5;
const double crRHS7 = N[0]*t[0] + N[1]*t[1] + N[2]*t[2] + N[3]*t[3];
const double crRHS8 = 2.0*crRHS1*crRHS7*functional_weights[7];
const double crRHS9 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2] + D[3]*N[3];
const double crRHS10 = DN(0,0)*t_adj[0] + DN(1,0)*t_adj[1] + DN(2,0)*t_adj[2] + DN(3,0)*t_adj[3];
const double crRHS11 = DN(0,1)*t_adj[0] + DN(1,1)*t_adj[1] + DN(2,1)*t_adj[2] + DN(3,1)*t_adj[3];
const double crRHS12 = DN(0,2)*t_adj[0] + DN(1,2)*t_adj[1] + DN(2,2)*t_adj[2] + DN(3,2)*t_adj[3];
const double crRHS13 = DN(0,0)*t[0] + DN(1,0)*t[1] + DN(2,0)*t[2] + DN(3,0)*t[3];
const double crRHS14 = DN(0,1)*t[0] + DN(1,1)*t[1] + DN(2,1)*t[2] + DN(3,1)*t[3];
const double crRHS15 = DN(0,2)*t[0] + DN(1,2)*t[1] + DN(2,2)*t[2] + DN(3,2)*t[3];
const double crRHS16 = 2.0*crRHS9*functional_weights[5];
const double crRHS17 = N[0]*vconv_adj(0,0) + N[1]*vconv_adj(1,0) + N[2]*vconv_adj(2,0) + N[3]*vconv_adj(3,0);
const double crRHS18 = DN(0,0)*crRHS17;
const double crRHS19 = N[0]*vconv_adj(0,1) + N[1]*vconv_adj(1,1) + N[2]*vconv_adj(2,1) + N[3]*vconv_adj(3,1);
const double crRHS20 = DN(0,1)*crRHS19;
const double crRHS21 = N[0]*vconv_adj(0,2) + N[1]*vconv_adj(1,2) + N[2]*vconv_adj(2,2) + N[3]*vconv_adj(3,2);
const double crRHS22 = DN(0,2)*crRHS21;
const double crRHS23 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3];
const double crRHS24 = crRHS23*crRHS5;
const double crRHS25 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crRHS26 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crRHS27 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crRHS28 = crRHS13*crRHS25 + crRHS14*crRHS26 + crRHS15*crRHS27;
const double crRHS29 = crRHS23*functional_weights[6];
const double crRHS30 = 1.0*(crRHS0 + crRHS10*crRHS17*crRHS23 + crRHS11*crRHS19*crRHS23 + crRHS12*crRHS21*crRHS23 - crRHS2 - crRHS3 + crRHS4 - crRHS6 + crRHS7*functional_weights[6]*(crRHS25*(DN(0,0)*c[0] + DN(1,0)*c[1] + DN(2,0)*c[2] + DN(3,0)*c[3]) + crRHS26*(DN(0,1)*c[0] + DN(1,1)*c[1] + DN(2,1)*c[2] + DN(3,1)*c[3]) + crRHS27*(DN(0,2)*c[0] + DN(1,2)*c[1] + DN(2,2)*c[2] + DN(3,2)*c[3])) - crRHS8)*1.0/(crRHS1*stab_c3 + crRHS9*stab_c1*1.0/(h*h) + stab_c2*fabs(crRHS23)*1.0/h*sqrt(crRHS17*crRHS17 + crRHS19*crRHS19 + crRHS21*crRHS21));
const double crRHS31 = crRHS1*crRHS30;
const double crRHS32 = crRHS23*crRHS30;
const double crRHS33 = DN(1,0)*crRHS17;
const double crRHS34 = DN(1,1)*crRHS19;
const double crRHS35 = DN(1,2)*crRHS21;
const double crRHS36 = DN(2,0)*crRHS17;
const double crRHS37 = DN(2,1)*crRHS19;
const double crRHS38 = DN(2,2)*crRHS21;
const double crRHS39 = DN(3,0)*crRHS17;
const double crRHS40 = DN(3,1)*crRHS19;
const double crRHS41 = DN(3,2)*crRHS21;
rRHS[0]+=-gauss_weight*(-N[0]*crRHS0 + N[0]*crRHS2 + N[0]*crRHS3 + N[0]*crRHS31 - N[0]*crRHS4 + N[0]*crRHS6 + N[0]*crRHS8 + crRHS16*(DN(0,0)*crRHS13 + DN(0,1)*crRHS14 + DN(0,2)*crRHS15) + crRHS18*crRHS32 + crRHS20*crRHS32 + crRHS22*crRHS32 + crRHS24*(crRHS18 + crRHS20 + crRHS22) + crRHS29*(N[0]*crRHS28 + crRHS7*(DN(0,0)*crRHS25 + DN(0,1)*crRHS26 + DN(0,2)*crRHS27)) + crRHS9*(DN(0,0)*crRHS10 + DN(0,1)*crRHS11 + DN(0,2)*crRHS12));
rRHS[1]+=-gauss_weight*(-N[1]*crRHS0 + N[1]*crRHS2 + N[1]*crRHS3 + N[1]*crRHS31 - N[1]*crRHS4 + N[1]*crRHS6 + N[1]*crRHS8 + crRHS16*(DN(1,0)*crRHS13 + DN(1,1)*crRHS14 + DN(1,2)*crRHS15) + crRHS24*(crRHS33 + crRHS34 + crRHS35) + crRHS29*(N[1]*crRHS28 + crRHS7*(DN(1,0)*crRHS25 + DN(1,1)*crRHS26 + DN(1,2)*crRHS27)) + crRHS32*crRHS33 + crRHS32*crRHS34 + crRHS32*crRHS35 + crRHS9*(DN(1,0)*crRHS10 + DN(1,1)*crRHS11 + DN(1,2)*crRHS12));
rRHS[2]+=-gauss_weight*(-N[2]*crRHS0 + N[2]*crRHS2 + N[2]*crRHS3 + N[2]*crRHS31 - N[2]*crRHS4 + N[2]*crRHS6 + N[2]*crRHS8 + crRHS16*(DN(2,0)*crRHS13 + DN(2,1)*crRHS14 + DN(2,2)*crRHS15) + crRHS24*(crRHS36 + crRHS37 + crRHS38) + crRHS29*(N[2]*crRHS28 + crRHS7*(DN(2,0)*crRHS25 + DN(2,1)*crRHS26 + DN(2,2)*crRHS27)) + crRHS32*crRHS36 + crRHS32*crRHS37 + crRHS32*crRHS38 + crRHS9*(DN(2,0)*crRHS10 + DN(2,1)*crRHS11 + DN(2,2)*crRHS12));
rRHS[3]+=-gauss_weight*(-N[3]*crRHS0 + N[3]*crRHS2 + N[3]*crRHS3 + N[3]*crRHS31 - N[3]*crRHS4 + N[3]*crRHS6 + N[3]*crRHS8 + crRHS16*(DN(3,0)*crRHS13 + DN(3,1)*crRHS14 + DN(3,2)*crRHS15) + crRHS24*(crRHS39 + crRHS40 + crRHS41) + crRHS29*(N[3]*crRHS28 + crRHS7*(DN(3,0)*crRHS25 + DN(3,1)*crRHS26 + DN(3,2)*crRHS27)) + crRHS32*crRHS39 + crRHS32*crRHS40 + crRHS32*crRHS41 + crRHS9*(DN(3,0)*crRHS10 + DN(3,1)*crRHS11 + DN(3,2)*crRHS12));
 
}

template <class TElementData>
void TransportTopologyOptimizationElement<TElementData>::AddVelocitySystem(
    TElementData& rData,
    MatrixType& rLocalLHS,
    VectorType& rLocalRHS) 
    {
    KRATOS_TRY;

    KRATOS_ERROR << "Calling base TransportTopologyOptimizationElement::AddVelocitySystem "
                    "implementation. This method is not supported by your "
                    "element."
                 << std::endl;

    KRATOS_CATCH("");
}

template <class TElementData>
void TransportTopologyOptimizationElement<TElementData>::AddMassLHS(
    TElementData& rData, MatrixType& rMassMatrix) {
    KRATOS_TRY;

    KRATOS_ERROR << "Calling base TransportTopologyOptimizationElement::AddMassLHS "
                    "implementation. This method is not supported by your "
                    "element."
                 << std::endl;

    KRATOS_CATCH("");
}

template <class TElementData>
void TransportTopologyOptimizationElement<TElementData>::GetCurrentValuesVector(
    const TElementData& rData,
    array_1d<double,LocalSize>& rValues) const {

    int local_index = 0;

    const auto& r_temperatures = rData.Temperature;

    for (unsigned int i = 0; i < NumNodes; ++i) {
        rValues[local_index++] = r_temperatures[i];  // Temperature Dof
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template< class TElementData >
void TransportTopologyOptimizationElement<TElementData>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    rSerializer.save("mpConstitutiveLaw",this->mpConstitutiveLaw);
}


template< class TElementData >
void TransportTopologyOptimizationElement<TElementData>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    rSerializer.load("mpConstitutiveLaw",this->mpConstitutiveLaw);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class TransportTopologyOptimizationElement< TransportTopologyOptimizationElementData<2,3,true> >;
template class TransportTopologyOptimizationElement< TransportTopologyOptimizationElementData<2,4,true> >;
template class TransportTopologyOptimizationElement< TransportTopologyOptimizationElementData<3,4,true> >;

///////////////////////////////////////////////////////////////////////////////////////////////////
}