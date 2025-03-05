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
    //substitute_lhs_2D3N_T 
    
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
    //substitute_lhs_2D4N_T 
    
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
    //substitute_lhs_3D4N_T
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
    //substitute_rhs_2D3N_T
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
    //substitute_rhs_2D4N_T
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
    //substitute_rhs_3D4N_T
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
    //substitute_lhs_2D3N_ADJ_T
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
    //substitute_lhs_2D4N_ADJ_T
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
    //substitute_lhs_3D4N_ADJ_T
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
    //substitute_rhs_2D3N_ADJ_T
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
    //substitute_rhs_2D4N_ADJ_T
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
    //substitute_rhs_3D4N_ADJ_T 
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