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
    const double D = rData.Conductivity;
    const double k = rData.Decay;

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
    const double crLHS0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crLHS1 = DN(0,0)*crLHS0;
const double crLHS2 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crLHS3 = DN(0,1)*crLHS2;
const double crLHS4 = crLHS1 + crLHS3;
const double crLHS5 = N[0]*k;
const double crLHS6 = 1.0/(D*stab_c1/pow(h, 2) + k*stab_c3 + stab_c2*sqrt(pow(crLHS0, 2) + pow(crLHS2, 2))/h);
const double crLHS7 = crLHS6*(crLHS4 + crLHS5);
const double crLHS8 = DN(1,0)*crLHS0;
const double crLHS9 = DN(1,1)*crLHS2;
const double crLHS10 = crLHS8 + crLHS9;
const double crLHS11 = N[1]*k;
const double crLHS12 = crLHS6*(crLHS10 + crLHS11);
const double crLHS13 = D*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1)) + N[1]*crLHS5;
const double crLHS14 = DN(2,0)*crLHS0;
const double crLHS15 = DN(2,1)*crLHS2;
const double crLHS16 = crLHS14 + crLHS15;
const double crLHS17 = N[2]*k;
const double crLHS18 = crLHS6*(crLHS16 + crLHS17);
const double crLHS19 = D*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1)) + N[2]*crLHS5;
const double crLHS20 = D*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1)) + N[2]*crLHS11;
rLHS(0,0)+=gauss_weight*(D*(pow(DN(0,0), 2) + pow(DN(0,1), 2)) + pow(N[0], 2)*k + N[0]*crLHS4 + crLHS1*crLHS7 + crLHS3*crLHS7 - crLHS5*crLHS7);
rLHS(0,1)+=gauss_weight*(N[0]*crLHS10 + crLHS1*crLHS12 + crLHS12*crLHS3 - crLHS12*crLHS5 + crLHS13);
rLHS(0,2)+=gauss_weight*(N[0]*crLHS16 + crLHS1*crLHS18 + crLHS18*crLHS3 - crLHS18*crLHS5 + crLHS19);
rLHS(1,0)+=gauss_weight*(N[1]*crLHS4 - crLHS11*crLHS7 + crLHS13 + crLHS7*crLHS8 + crLHS7*crLHS9);
rLHS(1,1)+=gauss_weight*(D*(pow(DN(1,0), 2) + pow(DN(1,1), 2)) + pow(N[1], 2)*k + N[1]*crLHS10 - crLHS11*crLHS12 + crLHS12*crLHS8 + crLHS12*crLHS9);
rLHS(1,2)+=gauss_weight*(N[1]*crLHS16 - crLHS11*crLHS18 + crLHS18*crLHS8 + crLHS18*crLHS9 + crLHS20);
rLHS(2,0)+=gauss_weight*(N[2]*crLHS4 + crLHS14*crLHS7 + crLHS15*crLHS7 - crLHS17*crLHS7 + crLHS19);
rLHS(2,1)+=gauss_weight*(N[2]*crLHS10 + crLHS12*crLHS14 + crLHS12*crLHS15 - crLHS12*crLHS17 + crLHS20);
rLHS(2,2)+=gauss_weight*(D*(pow(DN(2,0), 2) + pow(DN(2,1), 2)) + pow(N[2], 2)*k + N[2]*crLHS16 + crLHS14*crLHS18 + crLHS15*crLHS18 - crLHS17*crLHS18);
 
    
}

template <>
void TransportTopologyOptimizationElement<TransportTopologyOptimizationElementData<3,4,true>>::ComputeGaussPointLHSContribution(
    TransportTopologyOptimizationElementData<3,4,true>& rData,
    MatrixType& rLHS)
{
    const double D = rData.Conductivity;
    const double k = rData.Decay;

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
    const double crLHS0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crLHS1 = DN(0,0)*crLHS0;
const double crLHS2 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crLHS3 = DN(0,1)*crLHS2;
const double crLHS4 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crLHS5 = DN(0,2)*crLHS4;
const double crLHS6 = crLHS1 + crLHS3 + crLHS5;
const double crLHS7 = N[0]*k;
const double crLHS8 = 1.0/(D*stab_c1/pow(h, 2) + k*stab_c3 + stab_c2*sqrt(pow(crLHS0, 2) + pow(crLHS2, 2) + pow(crLHS4, 2))/h);
const double crLHS9 = crLHS8*(crLHS6 + crLHS7);
const double crLHS10 = DN(1,0)*crLHS0;
const double crLHS11 = DN(1,1)*crLHS2;
const double crLHS12 = DN(1,2)*crLHS4;
const double crLHS13 = crLHS10 + crLHS11 + crLHS12;
const double crLHS14 = N[1]*k;
const double crLHS15 = crLHS8*(crLHS13 + crLHS14);
const double crLHS16 = D*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1) + DN(0,2)*DN(1,2)) + N[1]*crLHS7;
const double crLHS17 = DN(2,0)*crLHS0;
const double crLHS18 = DN(2,1)*crLHS2;
const double crLHS19 = DN(2,2)*crLHS4;
const double crLHS20 = crLHS17 + crLHS18 + crLHS19;
const double crLHS21 = N[2]*k;
const double crLHS22 = crLHS8*(crLHS20 + crLHS21);
const double crLHS23 = D*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1) + DN(0,2)*DN(2,2)) + N[2]*crLHS7;
const double crLHS24 = DN(3,0)*crLHS0;
const double crLHS25 = DN(3,1)*crLHS2;
const double crLHS26 = DN(3,2)*crLHS4;
const double crLHS27 = crLHS24 + crLHS25 + crLHS26;
const double crLHS28 = N[3]*k;
const double crLHS29 = crLHS8*(crLHS27 + crLHS28);
const double crLHS30 = D*(DN(0,0)*DN(3,0) + DN(0,1)*DN(3,1) + DN(0,2)*DN(3,2)) + N[3]*crLHS7;
const double crLHS31 = D*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1) + DN(1,2)*DN(2,2)) + N[2]*crLHS14;
const double crLHS32 = D*(DN(1,0)*DN(3,0) + DN(1,1)*DN(3,1) + DN(1,2)*DN(3,2)) + N[3]*crLHS14;
const double crLHS33 = D*(DN(2,0)*DN(3,0) + DN(2,1)*DN(3,1) + DN(2,2)*DN(3,2)) + N[3]*crLHS21;
rLHS(0,0)+=gauss_weight*(D*(pow(DN(0,0), 2) + pow(DN(0,1), 2) + pow(DN(0,2), 2)) + pow(N[0], 2)*k + N[0]*crLHS6 + crLHS1*crLHS9 + crLHS3*crLHS9 + crLHS5*crLHS9 - crLHS7*crLHS9);
rLHS(0,1)+=gauss_weight*(N[0]*crLHS13 + crLHS1*crLHS15 + crLHS15*crLHS3 + crLHS15*crLHS5 - crLHS15*crLHS7 + crLHS16);
rLHS(0,2)+=gauss_weight*(N[0]*crLHS20 + crLHS1*crLHS22 + crLHS22*crLHS3 + crLHS22*crLHS5 - crLHS22*crLHS7 + crLHS23);
rLHS(0,3)+=gauss_weight*(N[0]*crLHS27 + crLHS1*crLHS29 + crLHS29*crLHS3 + crLHS29*crLHS5 - crLHS29*crLHS7 + crLHS30);
rLHS(1,0)+=gauss_weight*(N[1]*crLHS6 + crLHS10*crLHS9 + crLHS11*crLHS9 + crLHS12*crLHS9 - crLHS14*crLHS9 + crLHS16);
rLHS(1,1)+=gauss_weight*(D*(pow(DN(1,0), 2) + pow(DN(1,1), 2) + pow(DN(1,2), 2)) + pow(N[1], 2)*k + N[1]*crLHS13 + crLHS10*crLHS15 + crLHS11*crLHS15 + crLHS12*crLHS15 - crLHS14*crLHS15);
rLHS(1,2)+=gauss_weight*(N[1]*crLHS20 + crLHS10*crLHS22 + crLHS11*crLHS22 + crLHS12*crLHS22 - crLHS14*crLHS22 + crLHS31);
rLHS(1,3)+=gauss_weight*(N[1]*crLHS27 + crLHS10*crLHS29 + crLHS11*crLHS29 + crLHS12*crLHS29 - crLHS14*crLHS29 + crLHS32);
rLHS(2,0)+=gauss_weight*(N[2]*crLHS6 + crLHS17*crLHS9 + crLHS18*crLHS9 + crLHS19*crLHS9 - crLHS21*crLHS9 + crLHS23);
rLHS(2,1)+=gauss_weight*(N[2]*crLHS13 + crLHS15*crLHS17 + crLHS15*crLHS18 + crLHS15*crLHS19 - crLHS15*crLHS21 + crLHS31);
rLHS(2,2)+=gauss_weight*(D*(pow(DN(2,0), 2) + pow(DN(2,1), 2) + pow(DN(2,2), 2)) + pow(N[2], 2)*k + N[2]*crLHS20 + crLHS17*crLHS22 + crLHS18*crLHS22 + crLHS19*crLHS22 - crLHS21*crLHS22);
rLHS(2,3)+=gauss_weight*(N[2]*crLHS27 + crLHS17*crLHS29 + crLHS18*crLHS29 + crLHS19*crLHS29 - crLHS21*crLHS29 + crLHS33);
rLHS(3,0)+=gauss_weight*(N[3]*crLHS6 + crLHS24*crLHS9 + crLHS25*crLHS9 + crLHS26*crLHS9 - crLHS28*crLHS9 + crLHS30);
rLHS(3,1)+=gauss_weight*(N[3]*crLHS13 + crLHS15*crLHS24 + crLHS15*crLHS25 + crLHS15*crLHS26 - crLHS15*crLHS28 + crLHS32);
rLHS(3,2)+=gauss_weight*(N[3]*crLHS20 + crLHS22*crLHS24 + crLHS22*crLHS25 + crLHS22*crLHS26 - crLHS22*crLHS28 + crLHS33);
rLHS(3,3)+=gauss_weight*(D*(pow(DN(3,0), 2) + pow(DN(3,1), 2) + pow(DN(3,2), 2)) + pow(N[3], 2)*k + N[3]*crLHS27 + crLHS24*crLHS29 + crLHS25*crLHS29 + crLHS26*crLHS29 - crLHS28*crLHS29);

}

template <>
void TransportTopologyOptimizationElement<TransportTopologyOptimizationElementData<2,3,true>>::ComputeGaussPointRHSContribution(
    TransportTopologyOptimizationElementData<2,3,true>& rData,
    VectorType& rRHS)
{
    const double D = rData.Conductivity;
    const double k = rData.Decay;

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
const double crRHS1 = k*(N[0]*t[0] + N[1]*t[1] + N[2]*t[2]);
const double crRHS2 = DN(0,0)*t[0] + DN(1,0)*t[1] + DN(2,0)*t[2];
const double crRHS3 = DN(0,1)*t[0] + DN(1,1)*t[1] + DN(2,1)*t[2];
const double crRHS4 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crRHS5 = crRHS2*crRHS4;
const double crRHS6 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crRHS7 = crRHS3*crRHS6;
const double crRHS8 = crRHS5 + crRHS7;
const double crRHS9 = 1.0*(crRHS0 - crRHS1 - crRHS5 - crRHS7)/(D*stab_c1/pow(h, 2) + k*stab_c3 + stab_c2*sqrt(pow(crRHS4, 2) + pow(crRHS6, 2))/h);
const double crRHS10 = crRHS9*k;
const double crRHS11 = crRHS4*crRHS9;
const double crRHS12 = crRHS6*crRHS9;
rRHS[0]+=-gauss_weight*(D*(DN(0,0)*crRHS2 + DN(0,1)*crRHS3) - DN(0,0)*crRHS11 - DN(0,1)*crRHS12 - N[0]*crRHS0 + N[0]*crRHS1 + N[0]*crRHS10 + N[0]*crRHS8);
rRHS[1]+=-gauss_weight*(D*(DN(1,0)*crRHS2 + DN(1,1)*crRHS3) - DN(1,0)*crRHS11 - DN(1,1)*crRHS12 - N[1]*crRHS0 + N[1]*crRHS1 + N[1]*crRHS10 + N[1]*crRHS8);
rRHS[2]+=-gauss_weight*(D*(DN(2,0)*crRHS2 + DN(2,1)*crRHS3) - DN(2,0)*crRHS11 - DN(2,1)*crRHS12 - N[2]*crRHS0 + N[2]*crRHS1 + N[2]*crRHS10 + N[2]*crRHS8);

}

template <>
void TransportTopologyOptimizationElement<TransportTopologyOptimizationElementData<3,4,true>>::ComputeGaussPointRHSContribution(
    TransportTopologyOptimizationElementData<3,4,true>& rData,
    VectorType& rRHS)
{
    const double D = rData.Conductivity;
    const double k = rData.Decay;

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
    const double crRHS0 = N[0]*source[0] + N[1]*source[1] + N[2]*source[2] + N[3]*source[3];
const double crRHS1 = k*(N[0]*t[0] + N[1]*t[1] + N[2]*t[2] + N[3]*t[3]);
const double crRHS2 = DN(0,0)*t[0] + DN(1,0)*t[1] + DN(2,0)*t[2] + DN(3,0)*t[3];
const double crRHS3 = DN(0,1)*t[0] + DN(1,1)*t[1] + DN(2,1)*t[2] + DN(3,1)*t[3];
const double crRHS4 = DN(0,2)*t[0] + DN(1,2)*t[1] + DN(2,2)*t[2] + DN(3,2)*t[3];
const double crRHS5 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crRHS6 = crRHS2*crRHS5;
const double crRHS7 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crRHS8 = crRHS3*crRHS7;
const double crRHS9 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crRHS10 = crRHS4*crRHS9;
const double crRHS11 = crRHS10 + crRHS6 + crRHS8;
const double crRHS12 = 1.0*(crRHS0 - crRHS1 - crRHS10 - crRHS6 - crRHS8)/(D*stab_c1/pow(h, 2) + k*stab_c3 + stab_c2*sqrt(pow(crRHS5, 2) + pow(crRHS7, 2) + pow(crRHS9, 2))/h);
const double crRHS13 = crRHS12*k;
const double crRHS14 = crRHS12*crRHS5;
const double crRHS15 = crRHS12*crRHS7;
const double crRHS16 = crRHS12*crRHS9;
rRHS[0]+=-gauss_weight*(D*(DN(0,0)*crRHS2 + DN(0,1)*crRHS3 + DN(0,2)*crRHS4) - DN(0,0)*crRHS14 - DN(0,1)*crRHS15 - DN(0,2)*crRHS16 - N[0]*crRHS0 + N[0]*crRHS1 + N[0]*crRHS11 + N[0]*crRHS13);
rRHS[1]+=-gauss_weight*(D*(DN(1,0)*crRHS2 + DN(1,1)*crRHS3 + DN(1,2)*crRHS4) - DN(1,0)*crRHS14 - DN(1,1)*crRHS15 - DN(1,2)*crRHS16 - N[1]*crRHS0 + N[1]*crRHS1 + N[1]*crRHS11 + N[1]*crRHS13);
rRHS[2]+=-gauss_weight*(D*(DN(2,0)*crRHS2 + DN(2,1)*crRHS3 + DN(2,2)*crRHS4) - DN(2,0)*crRHS14 - DN(2,1)*crRHS15 - DN(2,2)*crRHS16 - N[2]*crRHS0 + N[2]*crRHS1 + N[2]*crRHS11 + N[2]*crRHS13);
rRHS[3]+=-gauss_weight*(D*(DN(3,0)*crRHS2 + DN(3,1)*crRHS3 + DN(3,2)*crRHS4) - DN(3,0)*crRHS14 - DN(3,1)*crRHS15 - DN(3,2)*crRHS16 - N[3]*crRHS0 + N[3]*crRHS1 + N[3]*crRHS11 + N[3]*crRHS13);

}

// BUILD ADJOINT NS SYSTEM
template <>
void TransportTopologyOptimizationElement< TransportTopologyOptimizationElementData<2,3,true>>::ComputeGaussPointLHSContributionAdjoint(
    TransportTopologyOptimizationElementData<2,3,true> & rData,
    MatrixType& rLHS)
{
    const double D = rData.Conductivity;
    const double k = rData.Decay;

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
    const double crLHS0 = N[0]*vconv_adj(0,0) + N[1]*vconv_adj(1,0) + N[2]*vconv_adj(2,0);
const double crLHS1 = DN(0,0)*crLHS0;
const double crLHS2 = N[0]*vconv_adj(0,1) + N[1]*vconv_adj(1,1) + N[2]*vconv_adj(2,1);
const double crLHS3 = DN(0,1)*crLHS2;
const double crLHS4 = crLHS1 + crLHS3;
const double crLHS5 = N[0]*k;
const double crLHS6 = 1.0/(D*stab_c1/pow(h, 2) + k*stab_c3 + stab_c2*sqrt(pow(crLHS0, 2) + pow(crLHS2, 2))/h);
const double crLHS7 = crLHS6*(crLHS4 - crLHS5);
const double crLHS8 = DN(1,0)*crLHS0;
const double crLHS9 = DN(1,1)*crLHS2;
const double crLHS10 = crLHS8 + crLHS9;
const double crLHS11 = N[1]*k;
const double crLHS12 = crLHS6*(crLHS10 - crLHS11);
const double crLHS13 = D*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1)) + N[1]*crLHS5;
const double crLHS14 = DN(2,0)*crLHS0;
const double crLHS15 = DN(2,1)*crLHS2;
const double crLHS16 = crLHS14 + crLHS15;
const double crLHS17 = N[2]*k;
const double crLHS18 = crLHS6*(crLHS16 - crLHS17);
const double crLHS19 = D*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1)) + N[2]*crLHS5;
const double crLHS20 = D*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1)) + N[2]*crLHS11;
rLHS(0,0)+=gauss_weight*(D*(pow(DN(0,0), 2) + pow(DN(0,1), 2)) + pow(N[0], 2)*k - N[0]*crLHS4 - crLHS1*crLHS7 - crLHS3*crLHS7 + crLHS5*crLHS7);
rLHS(0,1)+=gauss_weight*(-N[0]*crLHS10 - crLHS1*crLHS12 - crLHS12*crLHS3 + crLHS12*crLHS5 + crLHS13);
rLHS(0,2)+=gauss_weight*(-N[0]*crLHS16 - crLHS1*crLHS18 - crLHS18*crLHS3 + crLHS18*crLHS5 + crLHS19);
rLHS(1,0)+=gauss_weight*(-N[1]*crLHS4 + crLHS11*crLHS7 + crLHS13 - crLHS7*crLHS8 - crLHS7*crLHS9);
rLHS(1,1)+=gauss_weight*(D*(pow(DN(1,0), 2) + pow(DN(1,1), 2)) + pow(N[1], 2)*k - N[1]*crLHS10 + crLHS11*crLHS12 - crLHS12*crLHS8 - crLHS12*crLHS9);
rLHS(1,2)+=gauss_weight*(-N[1]*crLHS16 + crLHS11*crLHS18 - crLHS18*crLHS8 - crLHS18*crLHS9 + crLHS20);
rLHS(2,0)+=gauss_weight*(-N[2]*crLHS4 - crLHS14*crLHS7 - crLHS15*crLHS7 + crLHS17*crLHS7 + crLHS19);
rLHS(2,1)+=gauss_weight*(-N[2]*crLHS10 - crLHS12*crLHS14 - crLHS12*crLHS15 + crLHS12*crLHS17 + crLHS20);
rLHS(2,2)+=gauss_weight*(D*(pow(DN(2,0), 2) + pow(DN(2,1), 2)) + pow(N[2], 2)*k - N[2]*crLHS16 - crLHS14*crLHS18 - crLHS15*crLHS18 + crLHS17*crLHS18);

}

template <>
void TransportTopologyOptimizationElement<TransportTopologyOptimizationElementData<3,4,true>>::ComputeGaussPointLHSContributionAdjoint(
    TransportTopologyOptimizationElementData<3,4,true>& rData,
    MatrixType& rLHS)
{
    const double D = rData.Conductivity;
    const double k = rData.Decay;

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
    const double crLHS0 = N[0]*vconv_adj(0,0) + N[1]*vconv_adj(1,0) + N[2]*vconv_adj(2,0) + N[3]*vconv_adj(3,0);
const double crLHS1 = DN(0,0)*crLHS0;
const double crLHS2 = N[0]*vconv_adj(0,1) + N[1]*vconv_adj(1,1) + N[2]*vconv_adj(2,1) + N[3]*vconv_adj(3,1);
const double crLHS3 = DN(0,1)*crLHS2;
const double crLHS4 = N[0]*vconv_adj(0,2) + N[1]*vconv_adj(1,2) + N[2]*vconv_adj(2,2) + N[3]*vconv_adj(3,2);
const double crLHS5 = DN(0,2)*crLHS4;
const double crLHS6 = crLHS1 + crLHS3 + crLHS5;
const double crLHS7 = N[0]*k;
const double crLHS8 = 1.0/(D*stab_c1/pow(h, 2) + k*stab_c3 + stab_c2*sqrt(pow(crLHS0, 2) + pow(crLHS2, 2) + pow(crLHS4, 2))/h);
const double crLHS9 = crLHS8*(crLHS6 - crLHS7);
const double crLHS10 = DN(1,0)*crLHS0;
const double crLHS11 = DN(1,1)*crLHS2;
const double crLHS12 = DN(1,2)*crLHS4;
const double crLHS13 = crLHS10 + crLHS11 + crLHS12;
const double crLHS14 = N[1]*k;
const double crLHS15 = crLHS8*(crLHS13 - crLHS14);
const double crLHS16 = -D*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1) + DN(0,2)*DN(1,2)) - N[1]*crLHS7;
const double crLHS17 = DN(2,0)*crLHS0;
const double crLHS18 = DN(2,1)*crLHS2;
const double crLHS19 = DN(2,2)*crLHS4;
const double crLHS20 = crLHS17 + crLHS18 + crLHS19;
const double crLHS21 = N[2]*k;
const double crLHS22 = crLHS8*(crLHS20 - crLHS21);
const double crLHS23 = -D*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1) + DN(0,2)*DN(2,2)) - N[2]*crLHS7;
const double crLHS24 = DN(3,0)*crLHS0;
const double crLHS25 = DN(3,1)*crLHS2;
const double crLHS26 = DN(3,2)*crLHS4;
const double crLHS27 = crLHS24 + crLHS25 + crLHS26;
const double crLHS28 = N[3]*k;
const double crLHS29 = crLHS8*(crLHS27 - crLHS28);
const double crLHS30 = -D*(DN(0,0)*DN(3,0) + DN(0,1)*DN(3,1) + DN(0,2)*DN(3,2)) - N[3]*crLHS7;
const double crLHS31 = -D*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1) + DN(1,2)*DN(2,2)) - N[2]*crLHS14;
const double crLHS32 = -D*(DN(1,0)*DN(3,0) + DN(1,1)*DN(3,1) + DN(1,2)*DN(3,2)) - N[3]*crLHS14;
const double crLHS33 = -D*(DN(2,0)*DN(3,0) + DN(2,1)*DN(3,1) + DN(2,2)*DN(3,2)) - N[3]*crLHS21;
rLHS(0,0)+=-gauss_weight*(-D*(pow(DN(0,0), 2) + pow(DN(0,1), 2) + pow(DN(0,2), 2)) - pow(N[0], 2)*k + N[0]*crLHS6 + crLHS1*crLHS9 + crLHS3*crLHS9 + crLHS5*crLHS9 - crLHS7*crLHS9);
rLHS(0,1)+=-gauss_weight*(N[0]*crLHS13 + crLHS1*crLHS15 + crLHS15*crLHS3 + crLHS15*crLHS5 - crLHS15*crLHS7 + crLHS16);
rLHS(0,2)+=-gauss_weight*(N[0]*crLHS20 + crLHS1*crLHS22 + crLHS22*crLHS3 + crLHS22*crLHS5 - crLHS22*crLHS7 + crLHS23);
rLHS(0,3)+=-gauss_weight*(N[0]*crLHS27 + crLHS1*crLHS29 + crLHS29*crLHS3 + crLHS29*crLHS5 - crLHS29*crLHS7 + crLHS30);
rLHS(1,0)+=-gauss_weight*(N[1]*crLHS6 + crLHS10*crLHS9 + crLHS11*crLHS9 + crLHS12*crLHS9 - crLHS14*crLHS9 + crLHS16);
rLHS(1,1)+=-gauss_weight*(-D*(pow(DN(1,0), 2) + pow(DN(1,1), 2) + pow(DN(1,2), 2)) - pow(N[1], 2)*k + N[1]*crLHS13 + crLHS10*crLHS15 + crLHS11*crLHS15 + crLHS12*crLHS15 - crLHS14*crLHS15);
rLHS(1,2)+=-gauss_weight*(N[1]*crLHS20 + crLHS10*crLHS22 + crLHS11*crLHS22 + crLHS12*crLHS22 - crLHS14*crLHS22 + crLHS31);
rLHS(1,3)+=-gauss_weight*(N[1]*crLHS27 + crLHS10*crLHS29 + crLHS11*crLHS29 + crLHS12*crLHS29 - crLHS14*crLHS29 + crLHS32);
rLHS(2,0)+=-gauss_weight*(N[2]*crLHS6 + crLHS17*crLHS9 + crLHS18*crLHS9 + crLHS19*crLHS9 - crLHS21*crLHS9 + crLHS23);
rLHS(2,1)+=-gauss_weight*(N[2]*crLHS13 + crLHS15*crLHS17 + crLHS15*crLHS18 + crLHS15*crLHS19 - crLHS15*crLHS21 + crLHS31);
rLHS(2,2)+=-gauss_weight*(-D*(pow(DN(2,0), 2) + pow(DN(2,1), 2) + pow(DN(2,2), 2)) - pow(N[2], 2)*k + N[2]*crLHS20 + crLHS17*crLHS22 + crLHS18*crLHS22 + crLHS19*crLHS22 - crLHS21*crLHS22);
rLHS(2,3)+=-gauss_weight*(N[2]*crLHS27 + crLHS17*crLHS29 + crLHS18*crLHS29 + crLHS19*crLHS29 - crLHS21*crLHS29 + crLHS33);
rLHS(3,0)+=-gauss_weight*(N[3]*crLHS6 + crLHS24*crLHS9 + crLHS25*crLHS9 + crLHS26*crLHS9 - crLHS28*crLHS9 + crLHS30);
rLHS(3,1)+=-gauss_weight*(N[3]*crLHS13 + crLHS15*crLHS24 + crLHS15*crLHS25 + crLHS15*crLHS26 - crLHS15*crLHS28 + crLHS32);
rLHS(3,2)+=-gauss_weight*(N[3]*crLHS20 + crLHS22*crLHS24 + crLHS22*crLHS25 + crLHS22*crLHS26 - crLHS22*crLHS28 + crLHS33);
rLHS(3,3)+=-gauss_weight*(-D*(pow(DN(3,0), 2) + pow(DN(3,1), 2) + pow(DN(3,2), 2)) - pow(N[3], 2)*k + N[3]*crLHS27 + crLHS24*crLHS29 + crLHS25*crLHS29 + crLHS26*crLHS29 - crLHS28*crLHS29);

}

template <>
void TransportTopologyOptimizationElement<TransportTopologyOptimizationElementData<2,3,true>>::ComputeGaussPointRHSContributionAdjoint(
    TransportTopologyOptimizationElementData<2,3,true>& rData,
    VectorType& rRHS)
{
    const double D = rData.Conductivity;
    const double k = rData.Decay;

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

    const array_1d<double,3>& t_adj = rData.Temperature_adj;
    const array_1d<double,3>& tn_adj = rData.Temperature_adj_OldStep1;
    const array_1d<double,3>& tnn_adj = rData.Temperature_adj_OldStep2;
    const BoundedMatrix<double,2,3>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,2,3> vconv_adj = rData.ConvectiveVelocity - rData.MeshVelocity;
    const array_1d<double,3>& source_adj = rData.SourceTerm_adj;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT TRANSPORT ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*source_adj[0] + N[1]*source_adj[1] + N[2]*source_adj[2];
const double crRHS1 = k*(N[0]*t_adj[0] + N[1]*t_adj[1] + N[2]*t_adj[2]);
const double crRHS2 = DN(0,0)*t_adj[0] + DN(1,0)*t_adj[1] + DN(2,0)*t_adj[2];
const double crRHS3 = DN(0,1)*t_adj[0] + DN(1,1)*t_adj[1] + DN(2,1)*t_adj[2];
const double crRHS4 = N[0]*vconv_adj(0,0) + N[1]*vconv_adj(1,0) + N[2]*vconv_adj(2,0);
const double crRHS5 = N[0]*vconv_adj(0,1) + N[1]*vconv_adj(1,1) + N[2]*vconv_adj(2,1);
const double crRHS6 = crRHS2*crRHS4 + crRHS3*crRHS5;
const double crRHS7 = 1.0*(crRHS0 - crRHS1 + crRHS6)/(D*stab_c1/pow(h, 2) + k*stab_c3 + stab_c2*sqrt(pow(crRHS4, 2) + pow(crRHS5, 2))/h);
const double crRHS8 = crRHS7*k;
const double crRHS9 = crRHS4*crRHS7;
const double crRHS10 = crRHS5*crRHS7;
rRHS[0]+=gauss_weight*(-D*(DN(0,0)*crRHS2 + DN(0,1)*crRHS3) + DN(0,0)*crRHS9 + DN(0,1)*crRHS10 + N[0]*crRHS0 - N[0]*crRHS1 + N[0]*crRHS6 - N[0]*crRHS8);
rRHS[1]+=gauss_weight*(-D*(DN(1,0)*crRHS2 + DN(1,1)*crRHS3) + DN(1,0)*crRHS9 + DN(1,1)*crRHS10 + N[1]*crRHS0 - N[1]*crRHS1 + N[1]*crRHS6 - N[1]*crRHS8);
rRHS[2]+=gauss_weight*(-D*(DN(2,0)*crRHS2 + DN(2,1)*crRHS3) + DN(2,0)*crRHS9 + DN(2,1)*crRHS10 + N[2]*crRHS0 - N[2]*crRHS1 + N[2]*crRHS6 - N[2]*crRHS8);

}

template <>
void TransportTopologyOptimizationElement<TransportTopologyOptimizationElementData<3,4,true>>::ComputeGaussPointRHSContributionAdjoint(
    TransportTopologyOptimizationElementData<3,4,true>& rData,
    VectorType& rRHS)
{
    const double D = rData.Conductivity;
    const double k = rData.Decay;

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

    const array_1d<double,3>& t_adj = rData.Temperature_adj;
    const array_1d<double,3>& tn_adj = rData.Temperature_adj_OldStep1;
    const array_1d<double,3>& tnn_adj = rData.Temperature_adj_OldStep2;
    const BoundedMatrix<double,2,3>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,2,3> vconv_adj = rData.ConvectiveVelocity - rData.MeshVelocity;
    const array_1d<double,3>& source_adj = rData.SourceTerm_adj;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT TRANSPORT ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*source_adj[0] + N[1]*source_adj[1] + N[2]*source_adj[2] + N[3]*source_adj[3];
const double crRHS1 = k*(N[0]*t_adj[0] + N[1]*t_adj[1] + N[2]*t_adj[2] + N[3]*t_adj[3]);
const double crRHS2 = DN(0,0)*t_adj[0] + DN(1,0)*t_adj[1] + DN(2,0)*t_adj[2] + DN(3,0)*t_adj[3];
const double crRHS3 = DN(0,1)*t_adj[0] + DN(1,1)*t_adj[1] + DN(2,1)*t_adj[2] + DN(3,1)*t_adj[3];
const double crRHS4 = DN(0,2)*t_adj[0] + DN(1,2)*t_adj[1] + DN(2,2)*t_adj[2] + DN(3,2)*t_adj[3];
const double crRHS5 = N[0]*vconv_adj(0,0) + N[1]*vconv_adj(1,0) + N[2]*vconv_adj(2,0) + N[3]*vconv_adj(3,0);
const double crRHS6 = N[0]*vconv_adj(0,1) + N[1]*vconv_adj(1,1) + N[2]*vconv_adj(2,1) + N[3]*vconv_adj(3,1);
const double crRHS7 = N[0]*vconv_adj(0,2) + N[1]*vconv_adj(1,2) + N[2]*vconv_adj(2,2) + N[3]*vconv_adj(3,2);
const double crRHS8 = crRHS2*crRHS5 + crRHS3*crRHS6 + crRHS4*crRHS7;
const double crRHS9 = 1.0*(crRHS0 - crRHS1 + crRHS8)/(D*stab_c1/pow(h, 2) + k*stab_c3 + stab_c2*sqrt(pow(crRHS5, 2) + pow(crRHS6, 2) + pow(crRHS7, 2))/h);
const double crRHS10 = crRHS9*k;
const double crRHS11 = crRHS5*crRHS9;
const double crRHS12 = crRHS6*crRHS9;
const double crRHS13 = crRHS7*crRHS9;
rRHS[0]+=gauss_weight*(-D*(DN(0,0)*crRHS2 + DN(0,1)*crRHS3 + DN(0,2)*crRHS4) + DN(0,0)*crRHS11 + DN(0,1)*crRHS12 + DN(0,2)*crRHS13 + N[0]*crRHS0 - N[0]*crRHS1 - N[0]*crRHS10 + N[0]*crRHS8);
rRHS[1]+=gauss_weight*(-D*(DN(1,0)*crRHS2 + DN(1,1)*crRHS3 + DN(1,2)*crRHS4) + DN(1,0)*crRHS11 + DN(1,1)*crRHS12 + DN(1,2)*crRHS13 + N[1]*crRHS0 - N[1]*crRHS1 - N[1]*crRHS10 + N[1]*crRHS8);
rRHS[2]+=gauss_weight*(-D*(DN(2,0)*crRHS2 + DN(2,1)*crRHS3 + DN(2,2)*crRHS4) + DN(2,0)*crRHS11 + DN(2,1)*crRHS12 + DN(2,2)*crRHS13 + N[2]*crRHS0 - N[2]*crRHS1 - N[2]*crRHS10 + N[2]*crRHS8);
rRHS[3]+=gauss_weight*(-D*(DN(3,0)*crRHS2 + DN(3,1)*crRHS3 + DN(3,2)*crRHS4) + DN(3,0)*crRHS11 + DN(3,1)*crRHS12 + DN(3,2)*crRHS13 + N[3]*crRHS0 - N[3]*crRHS1 - N[3]*crRHS10 + N[3]*crRHS8);
 
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
template class TransportTopologyOptimizationElement< TransportTopologyOptimizationElementData<3,4,true> >;

///////////////////////////////////////////////////////////////////////////////////////////////////
}