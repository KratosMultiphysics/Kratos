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
const double crLHS13 = 1.0/(crLHS0*stab_c3 + crLHS1*stab_c1/pow(h, 2) + crLHS7*stab_c2*sqrt(pow(crLHS2, 2) + pow(crLHS4, 2))/h);
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
rLHS(0,0)+=-gauss_weight*(-pow(N[0], 2)*crLHS0 + 1.0*N[0]*crLHS0*crLHS12*crLHS13 - crLHS1*(pow(DN(0,0), 2) + pow(DN(0,1), 2)) - crLHS10*crLHS15 - crLHS11*crLHS15 - crLHS6*crLHS8);
rLHS(0,1)+=-gauss_weight*(1.0*N[0]*crLHS0*crLHS13*crLHS22 - crLHS10*crLHS23 - crLHS11*crLHS23 - crLHS18*crLHS8 - crLHS24);
rLHS(0,2)+=-gauss_weight*(1.0*N[0]*crLHS0*crLHS13*crLHS30 - crLHS10*crLHS31 - crLHS11*crLHS31 - crLHS27*crLHS8 - crLHS32);
rLHS(1,0)+=-gauss_weight*(1.0*N[1]*crLHS0*crLHS12*crLHS13 - crLHS15*crLHS20 - crLHS15*crLHS21 - crLHS24 - crLHS33*crLHS6);
rLHS(1,1)+=-gauss_weight*(-pow(N[1], 2)*crLHS0 + 1.0*N[1]*crLHS0*crLHS13*crLHS22 - crLHS1*(pow(DN(1,0), 2) + pow(DN(1,1), 2)) - crLHS18*crLHS33 - crLHS20*crLHS23 - crLHS21*crLHS23);
rLHS(1,2)+=-gauss_weight*(1.0*N[1]*crLHS0*crLHS13*crLHS30 - crLHS20*crLHS31 - crLHS21*crLHS31 - crLHS27*crLHS33 - crLHS34);
rLHS(2,0)+=-gauss_weight*(1.0*N[2]*crLHS0*crLHS12*crLHS13 - crLHS15*crLHS28 - crLHS15*crLHS29 - crLHS32 - crLHS35*crLHS6);
rLHS(2,1)+=-gauss_weight*(1.0*N[2]*crLHS0*crLHS13*crLHS22 - crLHS18*crLHS35 - crLHS23*crLHS28 - crLHS23*crLHS29 - crLHS34);
rLHS(2,2)+=-gauss_weight*(-pow(N[2], 2)*crLHS0 + 1.0*N[2]*crLHS0*crLHS13*crLHS30 - crLHS1*(pow(DN(2,0), 2) + pow(DN(2,1), 2)) - crLHS27*crLHS35 - crLHS28*crLHS31 - crLHS29*crLHS31);
 
    
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
const double crLHS16 = 1.0/(crLHS0*stab_c3 + crLHS1*stab_c1/pow(h, 2) + crLHS9*stab_c2*sqrt(pow(crLHS2, 2) + pow(crLHS4, 2) + pow(crLHS6, 2))/h);
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
rLHS(0,0)+=-gauss_weight*(-pow(N[0], 2)*crLHS0 + 1.0*N[0]*crLHS0*crLHS15*crLHS16 - crLHS1*(pow(DN(0,0), 2) + pow(DN(0,1), 2) + pow(DN(0,2), 2)) - crLHS10*crLHS8 - crLHS12*crLHS18 - crLHS13*crLHS18 - crLHS14*crLHS18);
rLHS(0,1)+=-gauss_weight*(1.0*N[0]*crLHS0*crLHS16*crLHS27 - crLHS10*crLHS22 - crLHS12*crLHS28 - crLHS13*crLHS28 - crLHS14*crLHS28 - crLHS29);
rLHS(0,2)+=-gauss_weight*(1.0*N[0]*crLHS0*crLHS16*crLHS38 - crLHS10*crLHS33 - crLHS12*crLHS39 - crLHS13*crLHS39 - crLHS14*crLHS39 - crLHS40);
rLHS(0,3)+=-gauss_weight*(1.0*N[0]*crLHS0*crLHS16*crLHS48 - crLHS10*crLHS44 - crLHS12*crLHS49 - crLHS13*crLHS49 - crLHS14*crLHS49 - crLHS50);
rLHS(1,0)+=-gauss_weight*(1.0*N[1]*crLHS0*crLHS15*crLHS16 - crLHS18*crLHS24 - crLHS18*crLHS25 - crLHS18*crLHS26 - crLHS29 - crLHS51*crLHS8);
rLHS(1,1)+=-gauss_weight*(-pow(N[1], 2)*crLHS0 + 1.0*N[1]*crLHS0*crLHS16*crLHS27 - crLHS1*(pow(DN(1,0), 2) + pow(DN(1,1), 2) + pow(DN(1,2), 2)) - crLHS22*crLHS51 - crLHS24*crLHS28 - crLHS25*crLHS28 - crLHS26*crLHS28);
rLHS(1,2)+=-gauss_weight*(1.0*N[1]*crLHS0*crLHS16*crLHS38 - crLHS24*crLHS39 - crLHS25*crLHS39 - crLHS26*crLHS39 - crLHS33*crLHS51 - crLHS52);
rLHS(1,3)+=-gauss_weight*(1.0*N[1]*crLHS0*crLHS16*crLHS48 - crLHS24*crLHS49 - crLHS25*crLHS49 - crLHS26*crLHS49 - crLHS44*crLHS51 - crLHS53);
rLHS(2,0)+=-gauss_weight*(1.0*N[2]*crLHS0*crLHS15*crLHS16 - crLHS18*crLHS35 - crLHS18*crLHS36 - crLHS18*crLHS37 - crLHS40 - crLHS54*crLHS8);
rLHS(2,1)+=-gauss_weight*(1.0*N[2]*crLHS0*crLHS16*crLHS27 - crLHS22*crLHS54 - crLHS28*crLHS35 - crLHS28*crLHS36 - crLHS28*crLHS37 - crLHS52);
rLHS(2,2)+=-gauss_weight*(-pow(N[2], 2)*crLHS0 + 1.0*N[2]*crLHS0*crLHS16*crLHS38 - crLHS1*(pow(DN(2,0), 2) + pow(DN(2,1), 2) + pow(DN(2,2), 2)) - crLHS33*crLHS54 - crLHS35*crLHS39 - crLHS36*crLHS39 - crLHS37*crLHS39);
rLHS(2,3)+=-gauss_weight*(1.0*N[2]*crLHS0*crLHS16*crLHS48 - crLHS35*crLHS49 - crLHS36*crLHS49 - crLHS37*crLHS49 - crLHS44*crLHS54 - crLHS55);
rLHS(3,0)+=-gauss_weight*(1.0*N[3]*crLHS0*crLHS15*crLHS16 - crLHS18*crLHS45 - crLHS18*crLHS46 - crLHS18*crLHS47 - crLHS50 - crLHS56*crLHS8);
rLHS(3,1)+=-gauss_weight*(1.0*N[3]*crLHS0*crLHS16*crLHS27 - crLHS22*crLHS56 - crLHS28*crLHS45 - crLHS28*crLHS46 - crLHS28*crLHS47 - crLHS53);
rLHS(3,2)+=-gauss_weight*(1.0*N[3]*crLHS0*crLHS16*crLHS38 - crLHS33*crLHS56 - crLHS39*crLHS45 - crLHS39*crLHS46 - crLHS39*crLHS47 - crLHS55);
rLHS(3,3)+=-gauss_weight*(-pow(N[3], 2)*crLHS0 + 1.0*N[3]*crLHS0*crLHS16*crLHS48 - crLHS1*(pow(DN(3,0), 2) + pow(DN(3,1), 2) + pow(DN(3,2), 2)) - crLHS44*crLHS56 - crLHS45*crLHS49 - crLHS46*crLHS49 - crLHS47*crLHS49);

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
const double crRHS12 = 1.0*(crRHS0 - crRHS10*crRHS7 - crRHS10*crRHS9 - crRHS2)/(crRHS1*stab_c3 + crRHS10*stab_c2*sqrt(pow(crRHS6, 2) + pow(crRHS8, 2))/h + crRHS3*stab_c1/pow(h, 2));
const double crRHS13 = crRHS1*crRHS12;
const double crRHS14 = crRHS10*crRHS12;
const double crRHS15 = crRHS14*crRHS6;
const double crRHS16 = crRHS14*crRHS8;
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS15 - DN(0,1)*crRHS16 - N[0]*crRHS0 + N[0]*crRHS11 + N[0]*crRHS13 + N[0]*crRHS2 + crRHS3*(DN(0,0)*crRHS4 + DN(0,1)*crRHS5));
rRHS[1]+=-gauss_weight*(-DN(1,0)*crRHS15 - DN(1,1)*crRHS16 - N[1]*crRHS0 + N[1]*crRHS11 + N[1]*crRHS13 + N[1]*crRHS2 + crRHS3*(DN(1,0)*crRHS4 + DN(1,1)*crRHS5));
rRHS[2]+=-gauss_weight*(-DN(2,0)*crRHS15 - DN(2,1)*crRHS16 - N[2]*crRHS0 + N[2]*crRHS11 + N[2]*crRHS13 + N[2]*crRHS2 + crRHS3*(DN(2,0)*crRHS4 + DN(2,1)*crRHS5));

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
const double crRHS15 = 1.0*(crRHS0 - crRHS10*crRHS13 - crRHS12*crRHS13 - crRHS13*crRHS8 - crRHS2)/(crRHS1*stab_c3 + crRHS13*stab_c2*sqrt(pow(crRHS11, 2) + pow(crRHS7, 2) + pow(crRHS9, 2))/h + crRHS3*stab_c1/pow(h, 2));
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
const double crLHS2 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2];
const double crLHS3 = N[0]*vconv_adj(0,0) + N[1]*vconv_adj(1,0) + N[2]*vconv_adj(2,0);
const double crLHS4 = DN(0,0)*crLHS3;
const double crLHS5 = N[0]*vconv_adj(0,1) + N[1]*vconv_adj(1,1) + N[2]*vconv_adj(2,1);
const double crLHS6 = DN(0,1)*crLHS5;
const double crLHS7 = crLHS4 + crLHS6;
const double crLHS8 = N[0]*crLHS0;
const double crLHS9 = crLHS2*crLHS4;
const double crLHS10 = crLHS2*crLHS6;
const double crLHS11 = crLHS10 - crLHS8 + crLHS9;
const double crLHS12 = 1.0/(crLHS0*stab_c3 + crLHS1*stab_c1/pow(h, 2) + crLHS2*stab_c2*sqrt(pow(crLHS3, 2) + pow(crLHS5, 2))/h);
const double crLHS13 = crLHS12*crLHS8;
const double crLHS14 = crLHS11*crLHS12;
const double crLHS15 = DN(1,0)*crLHS3;
const double crLHS16 = DN(1,1)*crLHS5;
const double crLHS17 = crLHS15 + crLHS16;
const double crLHS18 = N[1]*crLHS0;
const double crLHS19 = crLHS15*crLHS2;
const double crLHS20 = crLHS16*crLHS2;
const double crLHS21 = -crLHS18 + crLHS19 + crLHS20;
const double crLHS22 = crLHS12*crLHS21;
const double crLHS23 = N[1]*crLHS8 + crLHS1*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1));
const double crLHS24 = DN(2,0)*crLHS3;
const double crLHS25 = DN(2,1)*crLHS5;
const double crLHS26 = crLHS24 + crLHS25;
const double crLHS27 = N[2]*crLHS0;
const double crLHS28 = crLHS2*crLHS24;
const double crLHS29 = crLHS2*crLHS25;
const double crLHS30 = -crLHS27 + crLHS28 + crLHS29;
const double crLHS31 = crLHS12*crLHS30;
const double crLHS32 = N[2]*crLHS8 + crLHS1*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1));
const double crLHS33 = crLHS12*crLHS18;
const double crLHS34 = N[2]*crLHS18 + crLHS1*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1));
const double crLHS35 = crLHS12*crLHS27;
rLHS(0,0)+=-gauss_weight*(-pow(N[0], 2)*crLHS0 + N[0]*crLHS2*crLHS7 - crLHS1*(pow(DN(0,0), 2) + pow(DN(0,1), 2)) - crLHS10*crLHS14 - crLHS11*crLHS13 - crLHS14*crLHS9);
rLHS(0,1)+=-gauss_weight*(N[0]*crLHS17*crLHS2 - crLHS10*crLHS22 - crLHS13*crLHS21 - crLHS22*crLHS9 - crLHS23);
rLHS(0,2)+=-gauss_weight*(N[0]*crLHS2*crLHS26 - crLHS10*crLHS31 - crLHS13*crLHS30 - crLHS31*crLHS9 - crLHS32);
rLHS(1,0)+=-gauss_weight*(N[1]*crLHS2*crLHS7 - crLHS11*crLHS33 - crLHS14*crLHS19 - crLHS14*crLHS20 - crLHS23);
rLHS(1,1)+=-gauss_weight*(-pow(N[1], 2)*crLHS0 + N[1]*crLHS17*crLHS2 - crLHS1*(pow(DN(1,0), 2) + pow(DN(1,1), 2)) - crLHS19*crLHS22 - crLHS20*crLHS22 - crLHS21*crLHS33);
rLHS(1,2)+=-gauss_weight*(N[1]*crLHS2*crLHS26 - crLHS19*crLHS31 - crLHS20*crLHS31 - crLHS30*crLHS33 - crLHS34);
rLHS(2,0)+=-gauss_weight*(N[2]*crLHS2*crLHS7 - crLHS11*crLHS35 - crLHS14*crLHS28 - crLHS14*crLHS29 - crLHS32);
rLHS(2,1)+=-gauss_weight*(N[2]*crLHS17*crLHS2 - crLHS21*crLHS35 - crLHS22*crLHS28 - crLHS22*crLHS29 - crLHS34);
rLHS(2,2)+=-gauss_weight*(-pow(N[2], 2)*crLHS0 + N[2]*crLHS2*crLHS26 - crLHS1*(pow(DN(2,0), 2) + pow(DN(2,1), 2)) - crLHS28*crLHS31 - crLHS29*crLHS31 - crLHS30*crLHS35);

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
const double crLHS2 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3];
const double crLHS3 = N[0]*vconv_adj(0,0) + N[1]*vconv_adj(1,0) + N[2]*vconv_adj(2,0) + N[3]*vconv_adj(3,0);
const double crLHS4 = DN(0,0)*crLHS3;
const double crLHS5 = N[0]*vconv_adj(0,1) + N[1]*vconv_adj(1,1) + N[2]*vconv_adj(2,1) + N[3]*vconv_adj(3,1);
const double crLHS6 = DN(0,1)*crLHS5;
const double crLHS7 = N[0]*vconv_adj(0,2) + N[1]*vconv_adj(1,2) + N[2]*vconv_adj(2,2) + N[3]*vconv_adj(3,2);
const double crLHS8 = DN(0,2)*crLHS7;
const double crLHS9 = crLHS4 + crLHS6 + crLHS8;
const double crLHS10 = N[0]*crLHS0;
const double crLHS11 = crLHS2*crLHS4;
const double crLHS12 = crLHS2*crLHS6;
const double crLHS13 = crLHS2*crLHS8;
const double crLHS14 = -crLHS10 + crLHS11 + crLHS12 + crLHS13;
const double crLHS15 = 1.0/(crLHS0*stab_c3 + crLHS1*stab_c1/pow(h, 2) + crLHS2*stab_c2*sqrt(pow(crLHS3, 2) + pow(crLHS5, 2) + pow(crLHS7, 2))/h);
const double crLHS16 = crLHS10*crLHS15;
const double crLHS17 = crLHS14*crLHS15;
const double crLHS18 = DN(1,0)*crLHS3;
const double crLHS19 = DN(1,1)*crLHS5;
const double crLHS20 = DN(1,2)*crLHS7;
const double crLHS21 = crLHS18 + crLHS19 + crLHS20;
const double crLHS22 = N[1]*crLHS0;
const double crLHS23 = crLHS18*crLHS2;
const double crLHS24 = crLHS19*crLHS2;
const double crLHS25 = crLHS2*crLHS20;
const double crLHS26 = -crLHS22 + crLHS23 + crLHS24 + crLHS25;
const double crLHS27 = crLHS15*crLHS26;
const double crLHS28 = N[1]*crLHS10 + crLHS1*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1) + DN(0,2)*DN(1,2));
const double crLHS29 = DN(2,0)*crLHS3;
const double crLHS30 = DN(2,1)*crLHS5;
const double crLHS31 = DN(2,2)*crLHS7;
const double crLHS32 = crLHS29 + crLHS30 + crLHS31;
const double crLHS33 = N[2]*crLHS0;
const double crLHS34 = crLHS2*crLHS29;
const double crLHS35 = crLHS2*crLHS30;
const double crLHS36 = crLHS2*crLHS31;
const double crLHS37 = -crLHS33 + crLHS34 + crLHS35 + crLHS36;
const double crLHS38 = crLHS15*crLHS37;
const double crLHS39 = N[2]*crLHS10 + crLHS1*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1) + DN(0,2)*DN(2,2));
const double crLHS40 = DN(3,0)*crLHS3;
const double crLHS41 = DN(3,1)*crLHS5;
const double crLHS42 = DN(3,2)*crLHS7;
const double crLHS43 = crLHS40 + crLHS41 + crLHS42;
const double crLHS44 = N[3]*crLHS0;
const double crLHS45 = crLHS2*crLHS40;
const double crLHS46 = crLHS2*crLHS41;
const double crLHS47 = crLHS2*crLHS42;
const double crLHS48 = -crLHS44 + crLHS45 + crLHS46 + crLHS47;
const double crLHS49 = crLHS15*crLHS48;
const double crLHS50 = N[3]*crLHS10 + crLHS1*(DN(0,0)*DN(3,0) + DN(0,1)*DN(3,1) + DN(0,2)*DN(3,2));
const double crLHS51 = crLHS15*crLHS22;
const double crLHS52 = N[2]*crLHS22 + crLHS1*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1) + DN(1,2)*DN(2,2));
const double crLHS53 = N[3]*crLHS22 + crLHS1*(DN(1,0)*DN(3,0) + DN(1,1)*DN(3,1) + DN(1,2)*DN(3,2));
const double crLHS54 = crLHS15*crLHS33;
const double crLHS55 = N[3]*crLHS33 + crLHS1*(DN(2,0)*DN(3,0) + DN(2,1)*DN(3,1) + DN(2,2)*DN(3,2));
const double crLHS56 = crLHS15*crLHS44;
rLHS(0,0)+=-gauss_weight*(-pow(N[0], 2)*crLHS0 + N[0]*crLHS2*crLHS9 - crLHS1*(pow(DN(0,0), 2) + pow(DN(0,1), 2) + pow(DN(0,2), 2)) - crLHS11*crLHS17 - crLHS12*crLHS17 - crLHS13*crLHS17 - crLHS14*crLHS16);
rLHS(0,1)+=-gauss_weight*(N[0]*crLHS2*crLHS21 - crLHS11*crLHS27 - crLHS12*crLHS27 - crLHS13*crLHS27 - crLHS16*crLHS26 - crLHS28);
rLHS(0,2)+=-gauss_weight*(N[0]*crLHS2*crLHS32 - crLHS11*crLHS38 - crLHS12*crLHS38 - crLHS13*crLHS38 - crLHS16*crLHS37 - crLHS39);
rLHS(0,3)+=-gauss_weight*(N[0]*crLHS2*crLHS43 - crLHS11*crLHS49 - crLHS12*crLHS49 - crLHS13*crLHS49 - crLHS16*crLHS48 - crLHS50);
rLHS(1,0)+=-gauss_weight*(N[1]*crLHS2*crLHS9 - crLHS14*crLHS51 - crLHS17*crLHS23 - crLHS17*crLHS24 - crLHS17*crLHS25 - crLHS28);
rLHS(1,1)+=-gauss_weight*(-pow(N[1], 2)*crLHS0 + N[1]*crLHS2*crLHS21 - crLHS1*(pow(DN(1,0), 2) + pow(DN(1,1), 2) + pow(DN(1,2), 2)) - crLHS23*crLHS27 - crLHS24*crLHS27 - crLHS25*crLHS27 - crLHS26*crLHS51);
rLHS(1,2)+=-gauss_weight*(N[1]*crLHS2*crLHS32 - crLHS23*crLHS38 - crLHS24*crLHS38 - crLHS25*crLHS38 - crLHS37*crLHS51 - crLHS52);
rLHS(1,3)+=-gauss_weight*(N[1]*crLHS2*crLHS43 - crLHS23*crLHS49 - crLHS24*crLHS49 - crLHS25*crLHS49 - crLHS48*crLHS51 - crLHS53);
rLHS(2,0)+=-gauss_weight*(N[2]*crLHS2*crLHS9 - crLHS14*crLHS54 - crLHS17*crLHS34 - crLHS17*crLHS35 - crLHS17*crLHS36 - crLHS39);
rLHS(2,1)+=-gauss_weight*(N[2]*crLHS2*crLHS21 - crLHS26*crLHS54 - crLHS27*crLHS34 - crLHS27*crLHS35 - crLHS27*crLHS36 - crLHS52);
rLHS(2,2)+=-gauss_weight*(-pow(N[2], 2)*crLHS0 + N[2]*crLHS2*crLHS32 - crLHS1*(pow(DN(2,0), 2) + pow(DN(2,1), 2) + pow(DN(2,2), 2)) - crLHS34*crLHS38 - crLHS35*crLHS38 - crLHS36*crLHS38 - crLHS37*crLHS54);
rLHS(2,3)+=-gauss_weight*(N[2]*crLHS2*crLHS43 - crLHS34*crLHS49 - crLHS35*crLHS49 - crLHS36*crLHS49 - crLHS48*crLHS54 - crLHS55);
rLHS(3,0)+=-gauss_weight*(N[3]*crLHS2*crLHS9 - crLHS14*crLHS56 - crLHS17*crLHS45 - crLHS17*crLHS46 - crLHS17*crLHS47 - crLHS50);
rLHS(3,1)+=-gauss_weight*(N[3]*crLHS2*crLHS21 - crLHS26*crLHS56 - crLHS27*crLHS45 - crLHS27*crLHS46 - crLHS27*crLHS47 - crLHS53);
rLHS(3,2)+=-gauss_weight*(N[3]*crLHS2*crLHS32 - crLHS37*crLHS56 - crLHS38*crLHS45 - crLHS38*crLHS46 - crLHS38*crLHS47 - crLHS55);
rLHS(3,3)+=-gauss_weight*(-pow(N[3], 2)*crLHS0 + N[3]*crLHS2*crLHS43 - crLHS1*(pow(DN(3,0), 2) + pow(DN(3,1), 2) + pow(DN(3,2), 2)) - crLHS45*crLHS49 - crLHS46*crLHS49 - crLHS47*crLHS49 - crLHS48*crLHS56);

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

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT TRANSPORT ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*source_adj[0] + N[1]*source_adj[1] + N[2]*source_adj[2];
const double crRHS1 = 2.0*functional_weights[4]*(N[0]*opt_t[0] + N[1]*opt_t[1] + N[2]*opt_t[2]);
const double crRHS2 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2];
const double crRHS3 = crRHS2*(N[0]*t_adj[0] + N[1]*t_adj[1] + N[2]*t_adj[2]);
const double crRHS4 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2];
const double crRHS5 = DN(0,0)*t_adj[0] + DN(1,0)*t_adj[1] + DN(2,0)*t_adj[2];
const double crRHS6 = DN(0,1)*t_adj[0] + DN(1,1)*t_adj[1] + DN(2,1)*t_adj[2];
const double crRHS7 = N[0]*vconv_adj(0,0) + N[1]*vconv_adj(1,0) + N[2]*vconv_adj(2,0);
const double crRHS8 = crRHS5*crRHS7;
const double crRHS9 = N[0]*vconv_adj(0,1) + N[1]*vconv_adj(1,1) + N[2]*vconv_adj(2,1);
const double crRHS10 = crRHS6*crRHS9;
const double crRHS11 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2];
const double crRHS12 = crRHS11*(crRHS10 + crRHS8);
const double crRHS13 = 1.0*(crRHS0 - crRHS1 + crRHS10*crRHS11 + crRHS11*crRHS8 - crRHS3)/(crRHS11*stab_c2*sqrt(pow(crRHS7, 2) + pow(crRHS9, 2))/h + crRHS2*stab_c3 + crRHS4*stab_c1/pow(h, 2));
const double crRHS14 = crRHS13*crRHS2;
const double crRHS15 = crRHS11*crRHS13;
const double crRHS16 = crRHS15*crRHS7;
const double crRHS17 = crRHS15*crRHS9;
rRHS[0]+=-gauss_weight*(DN(0,0)*crRHS16 + DN(0,1)*crRHS17 - N[0]*crRHS0 + N[0]*crRHS1 - N[0]*crRHS12 + N[0]*crRHS14 + N[0]*crRHS3 + crRHS4*(DN(0,0)*crRHS5 + DN(0,1)*crRHS6));
rRHS[1]+=-gauss_weight*(DN(1,0)*crRHS16 + DN(1,1)*crRHS17 - N[1]*crRHS0 + N[1]*crRHS1 - N[1]*crRHS12 + N[1]*crRHS14 + N[1]*crRHS3 + crRHS4*(DN(1,0)*crRHS5 + DN(1,1)*crRHS6));
rRHS[2]+=-gauss_weight*(DN(2,0)*crRHS16 + DN(2,1)*crRHS17 - N[2]*crRHS0 + N[2]*crRHS1 - N[2]*crRHS12 + N[2]*crRHS14 + N[2]*crRHS3 + crRHS4*(DN(2,0)*crRHS5 + DN(2,1)*crRHS6));

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

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT TRANSPORT ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*source_adj[0] + N[1]*source_adj[1] + N[2]*source_adj[2] + N[3]*source_adj[3];
const double crRHS1 = 2.0*functional_weights[4]*(N[0]*opt_t[0] + N[1]*opt_t[1] + N[2]*opt_t[2] + N[3]*opt_t[3]);
const double crRHS2 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2] + N[3]*k[3];
const double crRHS3 = crRHS2*(N[0]*t_adj[0] + N[1]*t_adj[1] + N[2]*t_adj[2] + N[3]*t_adj[3]);
const double crRHS4 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2] + D[3]*N[3];
const double crRHS5 = DN(0,0)*t_adj[0] + DN(1,0)*t_adj[1] + DN(2,0)*t_adj[2] + DN(3,0)*t_adj[3];
const double crRHS6 = DN(0,1)*t_adj[0] + DN(1,1)*t_adj[1] + DN(2,1)*t_adj[2] + DN(3,1)*t_adj[3];
const double crRHS7 = DN(0,2)*t_adj[0] + DN(1,2)*t_adj[1] + DN(2,2)*t_adj[2] + DN(3,2)*t_adj[3];
const double crRHS8 = N[0]*vconv_adj(0,0) + N[1]*vconv_adj(1,0) + N[2]*vconv_adj(2,0) + N[3]*vconv_adj(3,0);
const double crRHS9 = crRHS5*crRHS8;
const double crRHS10 = N[0]*vconv_adj(0,1) + N[1]*vconv_adj(1,1) + N[2]*vconv_adj(2,1) + N[3]*vconv_adj(3,1);
const double crRHS11 = crRHS10*crRHS6;
const double crRHS12 = N[0]*vconv_adj(0,2) + N[1]*vconv_adj(1,2) + N[2]*vconv_adj(2,2) + N[3]*vconv_adj(3,2);
const double crRHS13 = crRHS12*crRHS7;
const double crRHS14 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3];
const double crRHS15 = crRHS14*(crRHS11 + crRHS13 + crRHS9);
const double crRHS16 = 1.0*(crRHS0 - crRHS1 + crRHS11*crRHS14 + crRHS13*crRHS14 + crRHS14*crRHS9 - crRHS3)/(crRHS14*stab_c2*sqrt(pow(crRHS10, 2) + pow(crRHS12, 2) + pow(crRHS8, 2))/h + crRHS2*stab_c3 + crRHS4*stab_c1/pow(h, 2));
const double crRHS17 = crRHS16*crRHS2;
const double crRHS18 = crRHS14*crRHS16;
const double crRHS19 = crRHS18*crRHS8;
const double crRHS20 = crRHS10*crRHS18;
const double crRHS21 = crRHS12*crRHS18;
rRHS[0]+=-gauss_weight*(DN(0,0)*crRHS19 + DN(0,1)*crRHS20 + DN(0,2)*crRHS21 - N[0]*crRHS0 + N[0]*crRHS1 - N[0]*crRHS15 + N[0]*crRHS17 + N[0]*crRHS3 + crRHS4*(DN(0,0)*crRHS5 + DN(0,1)*crRHS6 + DN(0,2)*crRHS7));
rRHS[1]+=-gauss_weight*(DN(1,0)*crRHS19 + DN(1,1)*crRHS20 + DN(1,2)*crRHS21 - N[1]*crRHS0 + N[1]*crRHS1 - N[1]*crRHS15 + N[1]*crRHS17 + N[1]*crRHS3 + crRHS4*(DN(1,0)*crRHS5 + DN(1,1)*crRHS6 + DN(1,2)*crRHS7));
rRHS[2]+=-gauss_weight*(DN(2,0)*crRHS19 + DN(2,1)*crRHS20 + DN(2,2)*crRHS21 - N[2]*crRHS0 + N[2]*crRHS1 - N[2]*crRHS15 + N[2]*crRHS17 + N[2]*crRHS3 + crRHS4*(DN(2,0)*crRHS5 + DN(2,1)*crRHS6 + DN(2,2)*crRHS7));
rRHS[3]+=-gauss_weight*(DN(3,0)*crRHS19 + DN(3,1)*crRHS20 + DN(3,2)*crRHS21 - N[3]*crRHS0 + N[3]*crRHS1 - N[3]*crRHS15 + N[3]*crRHS17 + N[3]*crRHS3 + crRHS4*(DN(3,0)*crRHS5 + DN(3,1)*crRHS6 + DN(3,2)*crRHS7));
 
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