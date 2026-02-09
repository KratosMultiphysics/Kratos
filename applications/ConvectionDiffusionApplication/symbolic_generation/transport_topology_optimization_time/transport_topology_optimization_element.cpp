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
    
    const double dt = rData.DeltaTime;
    const double theta = rData.TimeIntegrationTheta;
    const double time_coeff = rData.TopOptTimeCoefficient;

    const double dyn_tau = rData.DynamicTau;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const BoundedMatrix<double,2,3> vconv = rData.ConvectionVelocity - rData.MeshVelocity;

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;

    // TRANSPORT ELEMENTAL LHS MATRIX
    const double crLHS0 = pow(N[0], 2);
const double crLHS1 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2];
const double crLHS2 = crLHS1*theta;
const double crLHS3 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2];
const double crLHS4 = crLHS3*time_coeff/dt;
const double crLHS5 = 1.0*crLHS4;
const double crLHS6 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2];
const double crLHS7 = crLHS6*theta;
const double crLHS8 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crLHS9 = DN(0,0)*crLHS8;
const double crLHS10 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crLHS11 = DN(0,1)*crLHS10;
const double crLHS12 = crLHS11 + crLHS9;
const double crLHS13 = crLHS3*theta;
const double crLHS14 = N[0]*crLHS13;
const double crLHS15 = N[0]*crLHS2;
const double crLHS16 = 1.0*N[0]*crLHS4;
const double crLHS17 = crLHS11*crLHS13 + crLHS13*crLHS9 + crLHS15 + crLHS16;
const double crLHS18 = 1.0/(crLHS1*stab_c3 + crLHS4*dyn_tau + crLHS6*stab_c1/pow(h, 2) + stab_c2*sqrt(pow(crLHS10, 2) + pow(crLHS8, 2))*fabs(crLHS3)/h);
const double crLHS19 = 1.0*crLHS18*crLHS3;
const double crLHS20 = crLHS17*crLHS19;
const double crLHS21 = DN(1,0)*crLHS8;
const double crLHS22 = DN(1,1)*crLHS10;
const double crLHS23 = crLHS21 + crLHS22;
const double crLHS24 = N[1]*crLHS2;
const double crLHS25 = N[1]*crLHS5;
const double crLHS26 = crLHS13*crLHS21 + crLHS13*crLHS22 + crLHS24 + crLHS25;
const double crLHS27 = crLHS19*crLHS26;
const double crLHS28 = N[1]*crLHS15 + N[1]*crLHS16 + crLHS7*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1));
const double crLHS29 = DN(2,0)*crLHS8;
const double crLHS30 = DN(2,1)*crLHS10;
const double crLHS31 = crLHS29 + crLHS30;
const double crLHS32 = N[2]*crLHS2 + N[2]*crLHS5 + crLHS13*crLHS29 + crLHS13*crLHS30;
const double crLHS33 = crLHS19*crLHS32;
const double crLHS34 = N[2]*crLHS15 + N[2]*crLHS16 + crLHS7*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1));
const double crLHS35 = N[1]*crLHS13;
const double crLHS36 = pow(N[1], 2);
const double crLHS37 = N[2]*crLHS24 + N[2]*crLHS25 + crLHS7*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1));
const double crLHS38 = N[2]*crLHS13;
const double crLHS39 = pow(N[2], 2);
rLHS(0,0)+=-gauss_weight*(1.0*N[0]*crLHS1*crLHS17*crLHS18 - crLHS0*crLHS2 - crLHS0*crLHS5 - crLHS11*crLHS20 - crLHS12*crLHS14 - crLHS20*crLHS9 - crLHS7*(pow(DN(0,0), 2) + pow(DN(0,1), 2)));
rLHS(0,1)+=-gauss_weight*(1.0*N[0]*crLHS1*crLHS18*crLHS26 - crLHS11*crLHS27 - crLHS14*crLHS23 - crLHS27*crLHS9 - crLHS28);
rLHS(0,2)+=-gauss_weight*(1.0*N[0]*crLHS1*crLHS18*crLHS32 - crLHS11*crLHS33 - crLHS14*crLHS31 - crLHS33*crLHS9 - crLHS34);
rLHS(1,0)+=-gauss_weight*(1.0*N[1]*crLHS1*crLHS17*crLHS18 - crLHS12*crLHS35 - crLHS20*crLHS21 - crLHS20*crLHS22 - crLHS28);
rLHS(1,1)+=-gauss_weight*(1.0*N[1]*crLHS1*crLHS18*crLHS26 - crLHS2*crLHS36 - crLHS21*crLHS27 - crLHS22*crLHS27 - crLHS23*crLHS35 - crLHS36*crLHS5 - crLHS7*(pow(DN(1,0), 2) + pow(DN(1,1), 2)));
rLHS(1,2)+=-gauss_weight*(1.0*N[1]*crLHS1*crLHS18*crLHS32 - crLHS21*crLHS33 - crLHS22*crLHS33 - crLHS31*crLHS35 - crLHS37);
rLHS(2,0)+=-gauss_weight*(1.0*N[2]*crLHS1*crLHS17*crLHS18 - crLHS12*crLHS38 - crLHS20*crLHS29 - crLHS20*crLHS30 - crLHS34);
rLHS(2,1)+=-gauss_weight*(1.0*N[2]*crLHS1*crLHS18*crLHS26 - crLHS23*crLHS38 - crLHS27*crLHS29 - crLHS27*crLHS30 - crLHS37);
rLHS(2,2)+=-gauss_weight*(1.0*N[2]*crLHS1*crLHS18*crLHS32 - crLHS2*crLHS39 - crLHS29*crLHS33 - crLHS30*crLHS33 - crLHS31*crLHS38 - crLHS39*crLHS5 - crLHS7*(pow(DN(2,0), 2) + pow(DN(2,1), 2)));
 
    
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
    
    const double dt = rData.DeltaTime;
    const double theta = rData.TimeIntegrationTheta;
    const double time_coeff = rData.TopOptTimeCoefficient;

    const double dyn_tau = rData.DynamicTau;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,2>& DN = rData.DN_DX;

    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const BoundedMatrix<double,2,4> vconv = rData.ConvectionVelocity - rData.MeshVelocity;

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;

    // TRANSPORT ELEMENTAL LHS MATRIX
    const double crLHS0 = pow(N[0], 2);
const double crLHS1 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2] + N[3]*k[3];
const double crLHS2 = crLHS1*theta;
const double crLHS3 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3];
const double crLHS4 = crLHS3*time_coeff/dt;
const double crLHS5 = 1.0*crLHS4;
const double crLHS6 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2] + D[3]*N[3];
const double crLHS7 = crLHS6*theta;
const double crLHS8 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crLHS9 = DN(0,0)*crLHS8;
const double crLHS10 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crLHS11 = DN(0,1)*crLHS10;
const double crLHS12 = crLHS11 + crLHS9;
const double crLHS13 = crLHS3*theta;
const double crLHS14 = N[0]*crLHS13;
const double crLHS15 = N[0]*crLHS2;
const double crLHS16 = 1.0*N[0]*crLHS4;
const double crLHS17 = crLHS11*crLHS13 + crLHS13*crLHS9 + crLHS15 + crLHS16;
const double crLHS18 = 1.0/(crLHS1*stab_c3 + crLHS4*dyn_tau + crLHS6*stab_c1/pow(h, 2) + stab_c2*sqrt(pow(crLHS10, 2) + pow(crLHS8, 2))*fabs(crLHS3)/h);
const double crLHS19 = 1.0*crLHS18*crLHS3;
const double crLHS20 = crLHS17*crLHS19;
const double crLHS21 = DN(1,0)*crLHS8;
const double crLHS22 = DN(1,1)*crLHS10;
const double crLHS23 = crLHS21 + crLHS22;
const double crLHS24 = N[1]*crLHS2;
const double crLHS25 = N[1]*crLHS5;
const double crLHS26 = crLHS13*crLHS21 + crLHS13*crLHS22 + crLHS24 + crLHS25;
const double crLHS27 = crLHS19*crLHS26;
const double crLHS28 = N[1]*crLHS15 + N[1]*crLHS16 + crLHS7*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1));
const double crLHS29 = DN(2,0)*crLHS8;
const double crLHS30 = DN(2,1)*crLHS10;
const double crLHS31 = crLHS29 + crLHS30;
const double crLHS32 = N[2]*crLHS2;
const double crLHS33 = N[2]*crLHS5;
const double crLHS34 = crLHS13*crLHS29 + crLHS13*crLHS30 + crLHS32 + crLHS33;
const double crLHS35 = crLHS19*crLHS34;
const double crLHS36 = N[2]*crLHS15 + N[2]*crLHS16 + crLHS7*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1));
const double crLHS37 = DN(3,0)*crLHS8;
const double crLHS38 = DN(3,1)*crLHS10;
const double crLHS39 = crLHS37 + crLHS38;
const double crLHS40 = N[3]*crLHS2 + N[3]*crLHS5 + crLHS13*crLHS37 + crLHS13*crLHS38;
const double crLHS41 = crLHS19*crLHS40;
const double crLHS42 = N[3]*crLHS15 + N[3]*crLHS16 + crLHS7*(DN(0,0)*DN(3,0) + DN(0,1)*DN(3,1));
const double crLHS43 = N[1]*crLHS13;
const double crLHS44 = pow(N[1], 2);
const double crLHS45 = N[2]*crLHS24 + N[2]*crLHS25 + crLHS7*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1));
const double crLHS46 = N[3]*crLHS24 + N[3]*crLHS25 + crLHS7*(DN(1,0)*DN(3,0) + DN(1,1)*DN(3,1));
const double crLHS47 = N[2]*crLHS13;
const double crLHS48 = pow(N[2], 2);
const double crLHS49 = N[3]*crLHS32 + N[3]*crLHS33 + crLHS7*(DN(2,0)*DN(3,0) + DN(2,1)*DN(3,1));
const double crLHS50 = N[3]*crLHS13;
const double crLHS51 = pow(N[3], 2);
rLHS(0,0)+=-gauss_weight*(1.0*N[0]*crLHS1*crLHS17*crLHS18 - crLHS0*crLHS2 - crLHS0*crLHS5 - crLHS11*crLHS20 - crLHS12*crLHS14 - crLHS20*crLHS9 - crLHS7*(pow(DN(0,0), 2) + pow(DN(0,1), 2)));
rLHS(0,1)+=-gauss_weight*(1.0*N[0]*crLHS1*crLHS18*crLHS26 - crLHS11*crLHS27 - crLHS14*crLHS23 - crLHS27*crLHS9 - crLHS28);
rLHS(0,2)+=-gauss_weight*(1.0*N[0]*crLHS1*crLHS18*crLHS34 - crLHS11*crLHS35 - crLHS14*crLHS31 - crLHS35*crLHS9 - crLHS36);
rLHS(0,3)+=-gauss_weight*(1.0*N[0]*crLHS1*crLHS18*crLHS40 - crLHS11*crLHS41 - crLHS14*crLHS39 - crLHS41*crLHS9 - crLHS42);
rLHS(1,0)+=-gauss_weight*(1.0*N[1]*crLHS1*crLHS17*crLHS18 - crLHS12*crLHS43 - crLHS20*crLHS21 - crLHS20*crLHS22 - crLHS28);
rLHS(1,1)+=-gauss_weight*(1.0*N[1]*crLHS1*crLHS18*crLHS26 - crLHS2*crLHS44 - crLHS21*crLHS27 - crLHS22*crLHS27 - crLHS23*crLHS43 - crLHS44*crLHS5 - crLHS7*(pow(DN(1,0), 2) + pow(DN(1,1), 2)));
rLHS(1,2)+=-gauss_weight*(1.0*N[1]*crLHS1*crLHS18*crLHS34 - crLHS21*crLHS35 - crLHS22*crLHS35 - crLHS31*crLHS43 - crLHS45);
rLHS(1,3)+=-gauss_weight*(1.0*N[1]*crLHS1*crLHS18*crLHS40 - crLHS21*crLHS41 - crLHS22*crLHS41 - crLHS39*crLHS43 - crLHS46);
rLHS(2,0)+=-gauss_weight*(1.0*N[2]*crLHS1*crLHS17*crLHS18 - crLHS12*crLHS47 - crLHS20*crLHS29 - crLHS20*crLHS30 - crLHS36);
rLHS(2,1)+=-gauss_weight*(1.0*N[2]*crLHS1*crLHS18*crLHS26 - crLHS23*crLHS47 - crLHS27*crLHS29 - crLHS27*crLHS30 - crLHS45);
rLHS(2,2)+=-gauss_weight*(1.0*N[2]*crLHS1*crLHS18*crLHS34 - crLHS2*crLHS48 - crLHS29*crLHS35 - crLHS30*crLHS35 - crLHS31*crLHS47 - crLHS48*crLHS5 - crLHS7*(pow(DN(2,0), 2) + pow(DN(2,1), 2)));
rLHS(2,3)+=-gauss_weight*(1.0*N[2]*crLHS1*crLHS18*crLHS40 - crLHS29*crLHS41 - crLHS30*crLHS41 - crLHS39*crLHS47 - crLHS49);
rLHS(3,0)+=-gauss_weight*(1.0*N[3]*crLHS1*crLHS17*crLHS18 - crLHS12*crLHS50 - crLHS20*crLHS37 - crLHS20*crLHS38 - crLHS42);
rLHS(3,1)+=-gauss_weight*(1.0*N[3]*crLHS1*crLHS18*crLHS26 - crLHS23*crLHS50 - crLHS27*crLHS37 - crLHS27*crLHS38 - crLHS46);
rLHS(3,2)+=-gauss_weight*(1.0*N[3]*crLHS1*crLHS18*crLHS34 - crLHS31*crLHS50 - crLHS35*crLHS37 - crLHS35*crLHS38 - crLHS49);
rLHS(3,3)+=-gauss_weight*(1.0*N[3]*crLHS1*crLHS18*crLHS40 - crLHS2*crLHS51 - crLHS37*crLHS41 - crLHS38*crLHS41 - crLHS39*crLHS50 - crLHS5*crLHS51 - crLHS7*(pow(DN(3,0), 2) + pow(DN(3,1), 2)));
 
    
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
    
    const double dt = rData.DeltaTime;
    const double theta = rData.TimeIntegrationTheta;
    const double time_coeff = rData.TopOptTimeCoefficient;

    const double dyn_tau = rData.DynamicTau;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,3>& DN = rData.DN_DX;

    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const BoundedMatrix<double,3,4> vconv = rData.ConvectionVelocity - rData.MeshVelocity;

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;

    // TRANSPORT ELEMENTAL LHS MATRIX
    const double crLHS0 = pow(N[0], 2);
const double crLHS1 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2] + N[3]*k[3];
const double crLHS2 = crLHS1*theta;
const double crLHS3 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3];
const double crLHS4 = crLHS3*time_coeff/dt;
const double crLHS5 = 1.0*crLHS4;
const double crLHS6 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2] + D[3]*N[3];
const double crLHS7 = crLHS6*theta;
const double crLHS8 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crLHS9 = DN(0,0)*crLHS8;
const double crLHS10 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crLHS11 = DN(0,1)*crLHS10;
const double crLHS12 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crLHS13 = DN(0,2)*crLHS12;
const double crLHS14 = crLHS11 + crLHS13 + crLHS9;
const double crLHS15 = crLHS3*theta;
const double crLHS16 = N[0]*crLHS15;
const double crLHS17 = 1.0/(crLHS1*stab_c3 + crLHS4*dyn_tau + crLHS6*stab_c1/pow(h, 2) + stab_c2*sqrt(pow(crLHS10, 2) + pow(crLHS12, 2) + pow(crLHS8, 2))*fabs(crLHS3)/h);
const double crLHS18 = N[0]*crLHS2;
const double crLHS19 = 1.0*N[0]*crLHS4;
const double crLHS20 = crLHS11*crLHS15 + crLHS13*crLHS15 + crLHS15*crLHS9 + crLHS18 + crLHS19;
const double crLHS21 = 1.0*crLHS17*crLHS3;
const double crLHS22 = crLHS20*crLHS21;
const double crLHS23 = DN(1,0)*crLHS8;
const double crLHS24 = DN(1,1)*crLHS10;
const double crLHS25 = DN(1,2)*crLHS12;
const double crLHS26 = crLHS23 + crLHS24 + crLHS25;
const double crLHS27 = N[1]*crLHS2;
const double crLHS28 = N[1]*crLHS5;
const double crLHS29 = crLHS15*crLHS23 + crLHS15*crLHS24 + crLHS15*crLHS25 + crLHS27 + crLHS28;
const double crLHS30 = crLHS21*crLHS29;
const double crLHS31 = N[1]*crLHS18 + N[1]*crLHS19 + crLHS7*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1) + DN(0,2)*DN(1,2));
const double crLHS32 = DN(2,0)*crLHS8;
const double crLHS33 = DN(2,1)*crLHS10;
const double crLHS34 = DN(2,2)*crLHS12;
const double crLHS35 = crLHS32 + crLHS33 + crLHS34;
const double crLHS36 = N[2]*crLHS2;
const double crLHS37 = N[2]*crLHS5;
const double crLHS38 = crLHS15*crLHS32 + crLHS15*crLHS33 + crLHS15*crLHS34 + crLHS36 + crLHS37;
const double crLHS39 = crLHS21*crLHS38;
const double crLHS40 = N[2]*crLHS18 + N[2]*crLHS19 + crLHS7*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1) + DN(0,2)*DN(2,2));
const double crLHS41 = DN(3,0)*crLHS8;
const double crLHS42 = DN(3,1)*crLHS10;
const double crLHS43 = DN(3,2)*crLHS12;
const double crLHS44 = crLHS41 + crLHS42 + crLHS43;
const double crLHS45 = N[3]*crLHS2 + N[3]*crLHS5 + crLHS15*crLHS41 + crLHS15*crLHS42 + crLHS15*crLHS43;
const double crLHS46 = crLHS21*crLHS45;
const double crLHS47 = N[3]*crLHS18 + N[3]*crLHS19 + crLHS7*(DN(0,0)*DN(3,0) + DN(0,1)*DN(3,1) + DN(0,2)*DN(3,2));
const double crLHS48 = N[1]*crLHS15;
const double crLHS49 = pow(N[1], 2);
const double crLHS50 = N[2]*crLHS27 + N[2]*crLHS28 + crLHS7*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1) + DN(1,2)*DN(2,2));
const double crLHS51 = N[3]*crLHS27 + N[3]*crLHS28 + crLHS7*(DN(1,0)*DN(3,0) + DN(1,1)*DN(3,1) + DN(1,2)*DN(3,2));
const double crLHS52 = N[2]*crLHS15;
const double crLHS53 = pow(N[2], 2);
const double crLHS54 = N[3]*crLHS36 + N[3]*crLHS37 + crLHS7*(DN(2,0)*DN(3,0) + DN(2,1)*DN(3,1) + DN(2,2)*DN(3,2));
const double crLHS55 = N[3]*crLHS15;
const double crLHS56 = pow(N[3], 2);
rLHS(0,0)+=-gauss_weight*(1.0*N[0]*crLHS1*crLHS17*crLHS20 - crLHS0*crLHS2 - crLHS0*crLHS5 - crLHS11*crLHS22 - crLHS13*crLHS22 - crLHS14*crLHS16 - crLHS22*crLHS9 - crLHS7*(pow(DN(0,0), 2) + pow(DN(0,1), 2) + pow(DN(0,2), 2)));
rLHS(0,1)+=-gauss_weight*(1.0*N[0]*crLHS1*crLHS17*crLHS29 - crLHS11*crLHS30 - crLHS13*crLHS30 - crLHS16*crLHS26 - crLHS30*crLHS9 - crLHS31);
rLHS(0,2)+=-gauss_weight*(1.0*N[0]*crLHS1*crLHS17*crLHS38 - crLHS11*crLHS39 - crLHS13*crLHS39 - crLHS16*crLHS35 - crLHS39*crLHS9 - crLHS40);
rLHS(0,3)+=-gauss_weight*(1.0*N[0]*crLHS1*crLHS17*crLHS45 - crLHS11*crLHS46 - crLHS13*crLHS46 - crLHS16*crLHS44 - crLHS46*crLHS9 - crLHS47);
rLHS(1,0)+=-gauss_weight*(1.0*N[1]*crLHS1*crLHS17*crLHS20 - crLHS14*crLHS48 - crLHS22*crLHS23 - crLHS22*crLHS24 - crLHS22*crLHS25 - crLHS31);
rLHS(1,1)+=-gauss_weight*(1.0*N[1]*crLHS1*crLHS17*crLHS29 - crLHS2*crLHS49 - crLHS23*crLHS30 - crLHS24*crLHS30 - crLHS25*crLHS30 - crLHS26*crLHS48 - crLHS49*crLHS5 - crLHS7*(pow(DN(1,0), 2) + pow(DN(1,1), 2) + pow(DN(1,2), 2)));
rLHS(1,2)+=-gauss_weight*(1.0*N[1]*crLHS1*crLHS17*crLHS38 - crLHS23*crLHS39 - crLHS24*crLHS39 - crLHS25*crLHS39 - crLHS35*crLHS48 - crLHS50);
rLHS(1,3)+=-gauss_weight*(1.0*N[1]*crLHS1*crLHS17*crLHS45 - crLHS23*crLHS46 - crLHS24*crLHS46 - crLHS25*crLHS46 - crLHS44*crLHS48 - crLHS51);
rLHS(2,0)+=-gauss_weight*(1.0*N[2]*crLHS1*crLHS17*crLHS20 - crLHS14*crLHS52 - crLHS22*crLHS32 - crLHS22*crLHS33 - crLHS22*crLHS34 - crLHS40);
rLHS(2,1)+=-gauss_weight*(1.0*N[2]*crLHS1*crLHS17*crLHS29 - crLHS26*crLHS52 - crLHS30*crLHS32 - crLHS30*crLHS33 - crLHS30*crLHS34 - crLHS50);
rLHS(2,2)+=-gauss_weight*(1.0*N[2]*crLHS1*crLHS17*crLHS38 - crLHS2*crLHS53 - crLHS32*crLHS39 - crLHS33*crLHS39 - crLHS34*crLHS39 - crLHS35*crLHS52 - crLHS5*crLHS53 - crLHS7*(pow(DN(2,0), 2) + pow(DN(2,1), 2) + pow(DN(2,2), 2)));
rLHS(2,3)+=-gauss_weight*(1.0*N[2]*crLHS1*crLHS17*crLHS45 - crLHS32*crLHS46 - crLHS33*crLHS46 - crLHS34*crLHS46 - crLHS44*crLHS52 - crLHS54);
rLHS(3,0)+=-gauss_weight*(1.0*N[3]*crLHS1*crLHS17*crLHS20 - crLHS14*crLHS55 - crLHS22*crLHS41 - crLHS22*crLHS42 - crLHS22*crLHS43 - crLHS47);
rLHS(3,1)+=-gauss_weight*(1.0*N[3]*crLHS1*crLHS17*crLHS29 - crLHS26*crLHS55 - crLHS30*crLHS41 - crLHS30*crLHS42 - crLHS30*crLHS43 - crLHS51);
rLHS(3,2)+=-gauss_weight*(1.0*N[3]*crLHS1*crLHS17*crLHS38 - crLHS35*crLHS55 - crLHS39*crLHS41 - crLHS39*crLHS42 - crLHS39*crLHS43 - crLHS54);
rLHS(3,3)+=-gauss_weight*(1.0*N[3]*crLHS1*crLHS17*crLHS45 - crLHS2*crLHS56 - crLHS41*crLHS46 - crLHS42*crLHS46 - crLHS43*crLHS46 - crLHS44*crLHS55 - crLHS5*crLHS56 - crLHS7*(pow(DN(3,0), 2) + pow(DN(3,1), 2) + pow(DN(3,2), 2)));

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
    
    const double dt = rData.DeltaTime;
    const double theta = rData.TimeIntegrationTheta;
    const double time_coeff = rData.TopOptTimeCoefficient;

    const double dyn_tau = rData.DynamicTau;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const array_1d<double,3>& t = rData.Temperature;
    const array_1d<double,3>& tn = rData.Temperature_OldStep1;
    const BoundedMatrix<double,2,3>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,2,3> vconv = rData.ConvectionVelocity - rData.MeshVelocity;
    const array_1d<double,3>& source = rData.SourceTerm;
    const array_1d<double,3>& sourcen = rData.SourceTerm_OldStep1;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // TRANSPORT ELEMENTAL RHS VECTOR
    const double crRHS0 = theta - 1.0;
const double crRHS1 = -crRHS0*sourcen[0] + source[0]*theta;
const double crRHS2 = -crRHS0*sourcen[1] + source[1]*theta;
const double crRHS3 = -crRHS0*sourcen[2] + source[2]*theta;
const double crRHS4 = N[0]*crRHS1 + N[1]*crRHS2 + N[2]*crRHS3;
const double crRHS5 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2];
const double crRHS6 = crRHS5*time_coeff/dt;
const double crRHS7 = 1.0*crRHS6*(N[0]*(t[0] - tn[0]) + N[1]*(t[1] - tn[1]) + N[2]*(t[2] - tn[2]));
const double crRHS8 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2];
const double crRHS9 = -crRHS0*tn[0] + t[0]*theta;
const double crRHS10 = -crRHS0*tn[1] + t[1]*theta;
const double crRHS11 = -crRHS0*tn[2] + t[2]*theta;
const double crRHS12 = crRHS8*(N[0]*crRHS9 + N[1]*crRHS10 + N[2]*crRHS11);
const double crRHS13 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2];
const double crRHS14 = DN(0,0)*crRHS9 + DN(1,0)*crRHS10 + DN(2,0)*crRHS11;
const double crRHS15 = DN(0,1)*crRHS9 + DN(1,1)*crRHS10 + DN(2,1)*crRHS11;
const double crRHS16 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crRHS17 = crRHS14*crRHS16;
const double crRHS18 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crRHS19 = crRHS15*crRHS18;
const double crRHS20 = crRHS5*(crRHS17 + crRHS19);
const double crRHS21 = 1.0*(N[0]*crRHS1 + N[1]*crRHS2 + N[2]*crRHS3 - crRHS12 - crRHS17*crRHS5 - crRHS19*crRHS5 - crRHS7)/(crRHS13*stab_c1/pow(h, 2) + crRHS6*dyn_tau + crRHS8*stab_c3 + stab_c2*sqrt(pow(crRHS16, 2) + pow(crRHS18, 2))*fabs(crRHS5)/h);
const double crRHS22 = crRHS21*crRHS8;
const double crRHS23 = crRHS21*crRHS5;
const double crRHS24 = crRHS16*crRHS23;
const double crRHS25 = crRHS18*crRHS23;
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS24 - DN(0,1)*crRHS25 + N[0]*crRHS12 + N[0]*crRHS20 + N[0]*crRHS22 - N[0]*crRHS4 + N[0]*crRHS7 + crRHS13*(DN(0,0)*crRHS14 + DN(0,1)*crRHS15));
rRHS[1]+=-gauss_weight*(-DN(1,0)*crRHS24 - DN(1,1)*crRHS25 + N[1]*crRHS12 + N[1]*crRHS20 + N[1]*crRHS22 - N[1]*crRHS4 + N[1]*crRHS7 + crRHS13*(DN(1,0)*crRHS14 + DN(1,1)*crRHS15));
rRHS[2]+=-gauss_weight*(-DN(2,0)*crRHS24 - DN(2,1)*crRHS25 + N[2]*crRHS12 + N[2]*crRHS20 + N[2]*crRHS22 - N[2]*crRHS4 + N[2]*crRHS7 + crRHS13*(DN(2,0)*crRHS14 + DN(2,1)*crRHS15));

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
    
    const double dt = rData.DeltaTime;
    const double theta = rData.TimeIntegrationTheta;
    const double time_coeff = rData.TopOptTimeCoefficient;

    const double dyn_tau = rData.DynamicTau;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,2>& DN = rData.DN_DX;

    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const array_1d<double,4>& t = rData.Temperature;
    const array_1d<double,4>& tn = rData.Temperature_OldStep1;
    const BoundedMatrix<double,2,4>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,2,4> vconv = rData.ConvectionVelocity - rData.MeshVelocity;
    const array_1d<double,4>& source = rData.SourceTerm;
    const array_1d<double,4>& sourcen = rData.SourceTerm_OldStep1;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // TRANSPORT ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3];
const double crRHS1 = crRHS0*time_coeff/dt;
const double crRHS2 = 1.0*crRHS1*(N[0]*(t[0] - tn[0]) + N[1]*(t[1] - tn[1]) + N[2]*(t[2] - tn[2]) + N[3]*(t[3] - tn[3]));
const double crRHS3 = theta - 1.0;
const double crRHS4 = N[0]*(-crRHS3*sourcen[0] + source[0]*theta) + N[1]*(-crRHS3*sourcen[1] + source[1]*theta) + N[2]*(-crRHS3*sourcen[2] + source[2]*theta) + N[3]*(-crRHS3*sourcen[3] + source[3]*theta);
const double crRHS5 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2] + N[3]*k[3];
const double crRHS6 = -crRHS3*tn[0] + t[0]*theta;
const double crRHS7 = -crRHS3*tn[1] + t[1]*theta;
const double crRHS8 = -crRHS3*tn[2] + t[2]*theta;
const double crRHS9 = -crRHS3*tn[3] + t[3]*theta;
const double crRHS10 = crRHS5*(N[0]*crRHS6 + N[1]*crRHS7 + N[2]*crRHS8 + N[3]*crRHS9);
const double crRHS11 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2] + D[3]*N[3];
const double crRHS12 = DN(0,0)*crRHS6 + DN(1,0)*crRHS7 + DN(2,0)*crRHS8 + DN(3,0)*crRHS9;
const double crRHS13 = DN(0,1)*crRHS6 + DN(1,1)*crRHS7 + DN(2,1)*crRHS8 + DN(3,1)*crRHS9;
const double crRHS14 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crRHS15 = crRHS12*crRHS14;
const double crRHS16 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crRHS17 = crRHS13*crRHS16;
const double crRHS18 = crRHS0*(crRHS15 + crRHS17);
const double crRHS19 = 1.0*(-crRHS0*crRHS15 - crRHS0*crRHS17 - crRHS10 - crRHS2 + crRHS4)/(crRHS1*dyn_tau + crRHS11*stab_c1/pow(h, 2) + crRHS5*stab_c3 + stab_c2*sqrt(pow(crRHS14, 2) + pow(crRHS16, 2))*fabs(crRHS0)/h);
const double crRHS20 = crRHS19*crRHS5;
const double crRHS21 = crRHS0*crRHS19;
const double crRHS22 = crRHS14*crRHS21;
const double crRHS23 = crRHS16*crRHS21;
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS22 - DN(0,1)*crRHS23 + N[0]*crRHS10 + N[0]*crRHS18 + N[0]*crRHS2 + N[0]*crRHS20 - N[0]*crRHS4 + crRHS11*(DN(0,0)*crRHS12 + DN(0,1)*crRHS13));
rRHS[1]+=-gauss_weight*(-DN(1,0)*crRHS22 - DN(1,1)*crRHS23 + N[1]*crRHS10 + N[1]*crRHS18 + N[1]*crRHS2 + N[1]*crRHS20 - N[1]*crRHS4 + crRHS11*(DN(1,0)*crRHS12 + DN(1,1)*crRHS13));
rRHS[2]+=-gauss_weight*(-DN(2,0)*crRHS22 - DN(2,1)*crRHS23 + N[2]*crRHS10 + N[2]*crRHS18 + N[2]*crRHS2 + N[2]*crRHS20 - N[2]*crRHS4 + crRHS11*(DN(2,0)*crRHS12 + DN(2,1)*crRHS13));
rRHS[3]+=-gauss_weight*(-DN(3,0)*crRHS22 - DN(3,1)*crRHS23 + N[3]*crRHS10 + N[3]*crRHS18 + N[3]*crRHS2 + N[3]*crRHS20 - N[3]*crRHS4 + crRHS11*(DN(3,0)*crRHS12 + DN(3,1)*crRHS13));

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
    
    const double dt = rData.DeltaTime;
    const double theta = rData.TimeIntegrationTheta;
    const double time_coeff = rData.TopOptTimeCoefficient;

    const double dyn_tau = rData.DynamicTau;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,3>& DN = rData.DN_DX;

    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const array_1d<double,4>& t = rData.Temperature;
    const array_1d<double,4>& tn = rData.Temperature_OldStep1;
    const BoundedMatrix<double,3,4>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,3,4> vconv = rData.ConvectionVelocity - rData.MeshVelocity;
    const array_1d<double,4>& source = rData.SourceTerm;
    const array_1d<double,4>& sourcen = rData.SourceTerm_OldStep1;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // TRANSPORT ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3];
const double crRHS1 = crRHS0*time_coeff/dt;
const double crRHS2 = 1.0*crRHS1*(N[0]*(t[0] - tn[0]) + N[1]*(t[1] - tn[1]) + N[2]*(t[2] - tn[2]) + N[3]*(t[3] - tn[3]));
const double crRHS3 = theta - 1.0;
const double crRHS4 = -crRHS3*sourcen[0] + source[0]*theta;
const double crRHS5 = -crRHS3*sourcen[1] + source[1]*theta;
const double crRHS6 = -crRHS3*sourcen[2] + source[2]*theta;
const double crRHS7 = -crRHS3*sourcen[3] + source[3]*theta;
const double crRHS8 = N[0]*crRHS4 + N[1]*crRHS5 + N[2]*crRHS6 + N[3]*crRHS7;
const double crRHS9 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2] + N[3]*k[3];
const double crRHS10 = -crRHS3*tn[0] + t[0]*theta;
const double crRHS11 = -crRHS3*tn[1] + t[1]*theta;
const double crRHS12 = -crRHS3*tn[2] + t[2]*theta;
const double crRHS13 = -crRHS3*tn[3] + t[3]*theta;
const double crRHS14 = crRHS9*(N[0]*crRHS10 + N[1]*crRHS11 + N[2]*crRHS12 + N[3]*crRHS13);
const double crRHS15 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2] + D[3]*N[3];
const double crRHS16 = DN(0,0)*crRHS10 + DN(1,0)*crRHS11 + DN(2,0)*crRHS12 + DN(3,0)*crRHS13;
const double crRHS17 = DN(0,1)*crRHS10 + DN(1,1)*crRHS11 + DN(2,1)*crRHS12 + DN(3,1)*crRHS13;
const double crRHS18 = DN(0,2)*crRHS10 + DN(1,2)*crRHS11 + DN(2,2)*crRHS12 + DN(3,2)*crRHS13;
const double crRHS19 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crRHS20 = crRHS16*crRHS19;
const double crRHS21 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crRHS22 = crRHS17*crRHS21;
const double crRHS23 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crRHS24 = crRHS18*crRHS23;
const double crRHS25 = crRHS0*(crRHS20 + crRHS22 + crRHS24);
const double crRHS26 = 1.0*(N[0]*crRHS4 + N[1]*crRHS5 + N[2]*crRHS6 + N[3]*crRHS7 - crRHS0*crRHS20 - crRHS0*crRHS22 - crRHS0*crRHS24 - crRHS14 - crRHS2)/(crRHS1*dyn_tau + crRHS15*stab_c1/pow(h, 2) + crRHS9*stab_c3 + stab_c2*sqrt(pow(crRHS19, 2) + pow(crRHS21, 2) + pow(crRHS23, 2))*fabs(crRHS0)/h);
const double crRHS27 = crRHS26*crRHS9;
const double crRHS28 = crRHS0*crRHS26;
const double crRHS29 = crRHS19*crRHS28;
const double crRHS30 = crRHS21*crRHS28;
const double crRHS31 = crRHS23*crRHS28;
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS29 - DN(0,1)*crRHS30 - DN(0,2)*crRHS31 + N[0]*crRHS14 + N[0]*crRHS2 + N[0]*crRHS25 + N[0]*crRHS27 - N[0]*crRHS8 + crRHS15*(DN(0,0)*crRHS16 + DN(0,1)*crRHS17 + DN(0,2)*crRHS18));
rRHS[1]+=-gauss_weight*(-DN(1,0)*crRHS29 - DN(1,1)*crRHS30 - DN(1,2)*crRHS31 + N[1]*crRHS14 + N[1]*crRHS2 + N[1]*crRHS25 + N[1]*crRHS27 - N[1]*crRHS8 + crRHS15*(DN(1,0)*crRHS16 + DN(1,1)*crRHS17 + DN(1,2)*crRHS18));
rRHS[2]+=-gauss_weight*(-DN(2,0)*crRHS29 - DN(2,1)*crRHS30 - DN(2,2)*crRHS31 + N[2]*crRHS14 + N[2]*crRHS2 + N[2]*crRHS25 + N[2]*crRHS27 - N[2]*crRHS8 + crRHS15*(DN(2,0)*crRHS16 + DN(2,1)*crRHS17 + DN(2,2)*crRHS18));
rRHS[3]+=-gauss_weight*(-DN(3,0)*crRHS29 - DN(3,1)*crRHS30 - DN(3,2)*crRHS31 + N[3]*crRHS14 + N[3]*crRHS2 + N[3]*crRHS25 + N[3]*crRHS27 - N[3]*crRHS8 + crRHS15*(DN(3,0)*crRHS16 + DN(3,1)*crRHS17 + DN(3,2)*crRHS18));

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
    
    const double dt = rData.DeltaTime;
    const double theta = rData.TimeIntegrationTheta;
    const double time_coeff = rData.TopOptTimeCoefficient;

    const double dyn_tau = rData.DynamicTau;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const BoundedMatrix<double,2,3> vconv_adj = rData.ConvectionVelocity - rData.MeshVelocity;

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT TRANSPORT ELEMENTAL LHS MATRIX
    const double crLHS0 = pow(N[0], 2);
const double crLHS1 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2];
const double crLHS2 = crLHS1*theta;
const double crLHS3 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2];
const double crLHS4 = crLHS3*time_coeff/dt;
const double crLHS5 = 1.0*crLHS4;
const double crLHS6 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2];
const double crLHS7 = crLHS6*theta;
const double crLHS8 = N[0]*vconv_adj(0,0) + N[1]*vconv_adj(1,0) + N[2]*vconv_adj(2,0);
const double crLHS9 = DN(0,0)*crLHS8;
const double crLHS10 = N[0]*vconv_adj(0,1) + N[1]*vconv_adj(1,1) + N[2]*vconv_adj(2,1);
const double crLHS11 = DN(0,1)*crLHS10;
const double crLHS12 = crLHS3*theta;
const double crLHS13 = crLHS12*(crLHS11 + crLHS9);
const double crLHS14 = N[0]*crLHS1;
const double crLHS15 = crLHS3*crLHS9;
const double crLHS16 = crLHS11*crLHS3;
const double crLHS17 = -crLHS14 + crLHS15 + crLHS16;
const double crLHS18 = crLHS14*theta;
const double crLHS19 = 1.0/(crLHS1*stab_c3 + crLHS4*dyn_tau + crLHS6*stab_c1/pow(h, 2) + stab_c2*sqrt(pow(crLHS10, 2) + pow(crLHS8, 2))*fabs(crLHS3)/h);
const double crLHS20 = crLHS18*crLHS19;
const double crLHS21 = crLHS19*theta;
const double crLHS22 = crLHS17*crLHS21;
const double crLHS23 = N[1]*crLHS1;
const double crLHS24 = DN(1,0)*crLHS8;
const double crLHS25 = crLHS24*crLHS3;
const double crLHS26 = DN(1,1)*crLHS10;
const double crLHS27 = crLHS26*crLHS3;
const double crLHS28 = -crLHS23 + crLHS25 + crLHS27;
const double crLHS29 = crLHS21*crLHS28;
const double crLHS30 = N[0]*crLHS5;
const double crLHS31 = N[1]*crLHS18 + N[1]*crLHS30 + crLHS7*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1));
const double crLHS32 = N[2]*crLHS1;
const double crLHS33 = DN(2,0)*crLHS8;
const double crLHS34 = crLHS3*crLHS33;
const double crLHS35 = DN(2,1)*crLHS10;
const double crLHS36 = crLHS3*crLHS35;
const double crLHS37 = -crLHS32 + crLHS34 + crLHS36;
const double crLHS38 = crLHS21*crLHS37;
const double crLHS39 = N[2]*crLHS18 + N[2]*crLHS30 + crLHS7*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1));
const double crLHS40 = crLHS12*(crLHS24 + crLHS26);
const double crLHS41 = crLHS23*theta;
const double crLHS42 = crLHS19*crLHS41;
const double crLHS43 = pow(N[1], 2);
const double crLHS44 = N[1]*N[2]*crLHS5 + N[2]*crLHS41 + crLHS7*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1));
const double crLHS45 = crLHS12*(crLHS33 + crLHS35);
const double crLHS46 = crLHS21*crLHS32;
const double crLHS47 = pow(N[2], 2);
rLHS(0,0)+=gauss_weight*(N[0]*crLHS13 + crLHS0*crLHS2 + crLHS0*crLHS5 + crLHS15*crLHS22 + crLHS16*crLHS22 + crLHS17*crLHS20 + crLHS7*(pow(DN(0,0), 2) + pow(DN(0,1), 2)));
rLHS(0,1)+=gauss_weight*(N[1]*crLHS13 + crLHS15*crLHS29 + crLHS16*crLHS29 + crLHS20*crLHS28 + crLHS31);
rLHS(0,2)+=gauss_weight*(N[2]*crLHS13 + crLHS15*crLHS38 + crLHS16*crLHS38 + crLHS20*crLHS37 + crLHS39);
rLHS(1,0)+=gauss_weight*(N[0]*crLHS40 + crLHS17*crLHS42 + crLHS22*crLHS25 + crLHS22*crLHS27 + crLHS31);
rLHS(1,1)+=gauss_weight*(N[1]*crLHS40 + crLHS2*crLHS43 + crLHS25*crLHS29 + crLHS27*crLHS29 + crLHS28*crLHS42 + crLHS43*crLHS5 + crLHS7*(pow(DN(1,0), 2) + pow(DN(1,1), 2)));
rLHS(1,2)+=gauss_weight*(N[2]*crLHS40 + crLHS25*crLHS38 + crLHS27*crLHS38 + crLHS37*crLHS42 + crLHS44);
rLHS(2,0)+=gauss_weight*(N[0]*crLHS45 + crLHS17*crLHS46 + crLHS22*crLHS34 + crLHS22*crLHS36 + crLHS39);
rLHS(2,1)+=gauss_weight*(N[1]*crLHS45 + crLHS28*crLHS46 + crLHS29*crLHS34 + crLHS29*crLHS36 + crLHS44);
rLHS(2,2)+=gauss_weight*(N[2]*crLHS45 + crLHS2*crLHS47 + crLHS34*crLHS38 + crLHS36*crLHS38 + crLHS37*crLHS46 + crLHS47*crLHS5 + crLHS7*(pow(DN(2,0), 2) + pow(DN(2,1), 2)));

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
    
    const double dt = rData.DeltaTime;
    const double theta = rData.TimeIntegrationTheta;
    const double time_coeff = rData.TopOptTimeCoefficient;

    const double dyn_tau = rData.DynamicTau;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,2>& DN = rData.DN_DX;

    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const BoundedMatrix<double,2,4> vconv_adj = rData.ConvectionVelocity - rData.MeshVelocity;

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT TRANSPORT ELEMENTAL LHS MATRIX
    const double crLHS0 = pow(N[0], 2);
const double crLHS1 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2] + N[3]*k[3];
const double crLHS2 = crLHS1*theta;
const double crLHS3 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3];
const double crLHS4 = crLHS3*time_coeff/dt;
const double crLHS5 = 1.0*crLHS4;
const double crLHS6 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2] + D[3]*N[3];
const double crLHS7 = crLHS6*theta;
const double crLHS8 = N[0]*vconv_adj(0,0) + N[1]*vconv_adj(1,0) + N[2]*vconv_adj(2,0) + N[3]*vconv_adj(3,0);
const double crLHS9 = DN(0,0)*crLHS8;
const double crLHS10 = N[0]*vconv_adj(0,1) + N[1]*vconv_adj(1,1) + N[2]*vconv_adj(2,1) + N[3]*vconv_adj(3,1);
const double crLHS11 = DN(0,1)*crLHS10;
const double crLHS12 = crLHS3*theta;
const double crLHS13 = crLHS12*(crLHS11 + crLHS9);
const double crLHS14 = N[0]*crLHS1;
const double crLHS15 = crLHS3*crLHS9;
const double crLHS16 = crLHS11*crLHS3;
const double crLHS17 = -crLHS14 + crLHS15 + crLHS16;
const double crLHS18 = crLHS14*theta;
const double crLHS19 = 1.0/(crLHS1*stab_c3 + crLHS4*dyn_tau + crLHS6*stab_c1/pow(h, 2) + stab_c2*sqrt(pow(crLHS10, 2) + pow(crLHS8, 2))*fabs(crLHS3)/h);
const double crLHS20 = crLHS18*crLHS19;
const double crLHS21 = crLHS19*theta;
const double crLHS22 = crLHS17*crLHS21;
const double crLHS23 = N[1]*crLHS1;
const double crLHS24 = DN(1,0)*crLHS8;
const double crLHS25 = crLHS24*crLHS3;
const double crLHS26 = DN(1,1)*crLHS10;
const double crLHS27 = crLHS26*crLHS3;
const double crLHS28 = -crLHS23 + crLHS25 + crLHS27;
const double crLHS29 = crLHS21*crLHS28;
const double crLHS30 = N[0]*crLHS5;
const double crLHS31 = N[1]*crLHS18 + N[1]*crLHS30 + crLHS7*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1));
const double crLHS32 = N[2]*crLHS1;
const double crLHS33 = DN(2,0)*crLHS8;
const double crLHS34 = crLHS3*crLHS33;
const double crLHS35 = DN(2,1)*crLHS10;
const double crLHS36 = crLHS3*crLHS35;
const double crLHS37 = -crLHS32 + crLHS34 + crLHS36;
const double crLHS38 = crLHS21*crLHS37;
const double crLHS39 = N[2]*crLHS18 + N[2]*crLHS30 + crLHS7*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1));
const double crLHS40 = N[3]*crLHS1;
const double crLHS41 = DN(3,0)*crLHS8;
const double crLHS42 = crLHS3*crLHS41;
const double crLHS43 = DN(3,1)*crLHS10;
const double crLHS44 = crLHS3*crLHS43;
const double crLHS45 = -crLHS40 + crLHS42 + crLHS44;
const double crLHS46 = crLHS21*crLHS45;
const double crLHS47 = N[3]*crLHS18 + N[3]*crLHS30 + crLHS7*(DN(0,0)*DN(3,0) + DN(0,1)*DN(3,1));
const double crLHS48 = crLHS12*(crLHS24 + crLHS26);
const double crLHS49 = crLHS23*theta;
const double crLHS50 = crLHS19*crLHS49;
const double crLHS51 = pow(N[1], 2);
const double crLHS52 = N[1]*crLHS5;
const double crLHS53 = N[2]*crLHS49 + N[2]*crLHS52 + crLHS7*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1));
const double crLHS54 = N[3]*crLHS49 + N[3]*crLHS52 + crLHS7*(DN(1,0)*DN(3,0) + DN(1,1)*DN(3,1));
const double crLHS55 = crLHS12*(crLHS33 + crLHS35);
const double crLHS56 = crLHS32*theta;
const double crLHS57 = crLHS19*crLHS56;
const double crLHS58 = pow(N[2], 2);
const double crLHS59 = N[2]*N[3]*crLHS5 + N[3]*crLHS56 + crLHS7*(DN(2,0)*DN(3,0) + DN(2,1)*DN(3,1));
const double crLHS60 = crLHS12*(crLHS41 + crLHS43);
const double crLHS61 = crLHS21*crLHS40;
const double crLHS62 = pow(N[3], 2);
rLHS(0,0)+=gauss_weight*(N[0]*crLHS13 + crLHS0*crLHS2 + crLHS0*crLHS5 + crLHS15*crLHS22 + crLHS16*crLHS22 + crLHS17*crLHS20 + crLHS7*(pow(DN(0,0), 2) + pow(DN(0,1), 2)));
rLHS(0,1)+=gauss_weight*(N[1]*crLHS13 + crLHS15*crLHS29 + crLHS16*crLHS29 + crLHS20*crLHS28 + crLHS31);
rLHS(0,2)+=gauss_weight*(N[2]*crLHS13 + crLHS15*crLHS38 + crLHS16*crLHS38 + crLHS20*crLHS37 + crLHS39);
rLHS(0,3)+=gauss_weight*(N[3]*crLHS13 + crLHS15*crLHS46 + crLHS16*crLHS46 + crLHS20*crLHS45 + crLHS47);
rLHS(1,0)+=gauss_weight*(N[0]*crLHS48 + crLHS17*crLHS50 + crLHS22*crLHS25 + crLHS22*crLHS27 + crLHS31);
rLHS(1,1)+=gauss_weight*(N[1]*crLHS48 + crLHS2*crLHS51 + crLHS25*crLHS29 + crLHS27*crLHS29 + crLHS28*crLHS50 + crLHS5*crLHS51 + crLHS7*(pow(DN(1,0), 2) + pow(DN(1,1), 2)));
rLHS(1,2)+=gauss_weight*(N[2]*crLHS48 + crLHS25*crLHS38 + crLHS27*crLHS38 + crLHS37*crLHS50 + crLHS53);
rLHS(1,3)+=gauss_weight*(N[3]*crLHS48 + crLHS25*crLHS46 + crLHS27*crLHS46 + crLHS45*crLHS50 + crLHS54);
rLHS(2,0)+=gauss_weight*(N[0]*crLHS55 + crLHS17*crLHS57 + crLHS22*crLHS34 + crLHS22*crLHS36 + crLHS39);
rLHS(2,1)+=gauss_weight*(N[1]*crLHS55 + crLHS28*crLHS57 + crLHS29*crLHS34 + crLHS29*crLHS36 + crLHS53);
rLHS(2,2)+=gauss_weight*(N[2]*crLHS55 + crLHS2*crLHS58 + crLHS34*crLHS38 + crLHS36*crLHS38 + crLHS37*crLHS57 + crLHS5*crLHS58 + crLHS7*(pow(DN(2,0), 2) + pow(DN(2,1), 2)));
rLHS(2,3)+=gauss_weight*(N[3]*crLHS55 + crLHS34*crLHS46 + crLHS36*crLHS46 + crLHS45*crLHS57 + crLHS59);
rLHS(3,0)+=gauss_weight*(N[0]*crLHS60 + crLHS17*crLHS61 + crLHS22*crLHS42 + crLHS22*crLHS44 + crLHS47);
rLHS(3,1)+=gauss_weight*(N[1]*crLHS60 + crLHS28*crLHS61 + crLHS29*crLHS42 + crLHS29*crLHS44 + crLHS54);
rLHS(3,2)+=gauss_weight*(N[2]*crLHS60 + crLHS37*crLHS61 + crLHS38*crLHS42 + crLHS38*crLHS44 + crLHS59);
rLHS(3,3)+=gauss_weight*(N[3]*crLHS60 + crLHS2*crLHS62 + crLHS42*crLHS46 + crLHS44*crLHS46 + crLHS45*crLHS61 + crLHS5*crLHS62 + crLHS7*(pow(DN(3,0), 2) + pow(DN(3,1), 2)));

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
    
    const double dt = rData.DeltaTime;
    const double theta = rData.TimeIntegrationTheta;
    const double time_coeff = rData.TopOptTimeCoefficient;

    const double dyn_tau = rData.DynamicTau;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,3>& DN = rData.DN_DX;

    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const BoundedMatrix<double,3,4> vconv_adj = rData.ConvectionVelocity - rData.MeshVelocity;

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT TRANSPORT ELEMENTAL LHS MATRIX
    const double crLHS0 = pow(N[0], 2);
const double crLHS1 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2] + N[3]*k[3];
const double crLHS2 = crLHS1*theta;
const double crLHS3 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3];
const double crLHS4 = crLHS3*time_coeff/dt;
const double crLHS5 = 1.0*crLHS4;
const double crLHS6 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2] + D[3]*N[3];
const double crLHS7 = crLHS6*theta;
const double crLHS8 = N[0]*vconv_adj(0,0) + N[1]*vconv_adj(1,0) + N[2]*vconv_adj(2,0) + N[3]*vconv_adj(3,0);
const double crLHS9 = DN(0,0)*crLHS8;
const double crLHS10 = N[0]*vconv_adj(0,1) + N[1]*vconv_adj(1,1) + N[2]*vconv_adj(2,1) + N[3]*vconv_adj(3,1);
const double crLHS11 = DN(0,1)*crLHS10;
const double crLHS12 = N[0]*vconv_adj(0,2) + N[1]*vconv_adj(1,2) + N[2]*vconv_adj(2,2) + N[3]*vconv_adj(3,2);
const double crLHS13 = DN(0,2)*crLHS12;
const double crLHS14 = crLHS3*theta;
const double crLHS15 = crLHS14*(crLHS11 + crLHS13 + crLHS9);
const double crLHS16 = N[0]*crLHS1;
const double crLHS17 = crLHS3*crLHS9;
const double crLHS18 = crLHS11*crLHS3;
const double crLHS19 = crLHS13*crLHS3;
const double crLHS20 = -crLHS16 + crLHS17 + crLHS18 + crLHS19;
const double crLHS21 = crLHS16*theta;
const double crLHS22 = 1.0/(crLHS1*stab_c3 + crLHS4*dyn_tau + crLHS6*stab_c1/pow(h, 2) + stab_c2*sqrt(pow(crLHS10, 2) + pow(crLHS12, 2) + pow(crLHS8, 2))*fabs(crLHS3)/h);
const double crLHS23 = crLHS21*crLHS22;
const double crLHS24 = crLHS22*theta;
const double crLHS25 = crLHS20*crLHS24;
const double crLHS26 = N[1]*crLHS1;
const double crLHS27 = DN(1,0)*crLHS8;
const double crLHS28 = crLHS27*crLHS3;
const double crLHS29 = DN(1,1)*crLHS10;
const double crLHS30 = crLHS29*crLHS3;
const double crLHS31 = DN(1,2)*crLHS12;
const double crLHS32 = crLHS3*crLHS31;
const double crLHS33 = -crLHS26 + crLHS28 + crLHS30 + crLHS32;
const double crLHS34 = crLHS24*crLHS33;
const double crLHS35 = N[0]*crLHS5;
const double crLHS36 = N[1]*crLHS21 + N[1]*crLHS35 + crLHS7*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1) + DN(0,2)*DN(1,2));
const double crLHS37 = N[2]*crLHS1;
const double crLHS38 = DN(2,0)*crLHS8;
const double crLHS39 = crLHS3*crLHS38;
const double crLHS40 = DN(2,1)*crLHS10;
const double crLHS41 = crLHS3*crLHS40;
const double crLHS42 = DN(2,2)*crLHS12;
const double crLHS43 = crLHS3*crLHS42;
const double crLHS44 = -crLHS37 + crLHS39 + crLHS41 + crLHS43;
const double crLHS45 = crLHS24*crLHS44;
const double crLHS46 = N[2]*crLHS21 + N[2]*crLHS35 + crLHS7*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1) + DN(0,2)*DN(2,2));
const double crLHS47 = N[3]*crLHS1;
const double crLHS48 = DN(3,0)*crLHS8;
const double crLHS49 = crLHS3*crLHS48;
const double crLHS50 = DN(3,1)*crLHS10;
const double crLHS51 = crLHS3*crLHS50;
const double crLHS52 = DN(3,2)*crLHS12;
const double crLHS53 = crLHS3*crLHS52;
const double crLHS54 = -crLHS47 + crLHS49 + crLHS51 + crLHS53;
const double crLHS55 = crLHS24*crLHS54;
const double crLHS56 = N[3]*crLHS21 + N[3]*crLHS35 + crLHS7*(DN(0,0)*DN(3,0) + DN(0,1)*DN(3,1) + DN(0,2)*DN(3,2));
const double crLHS57 = crLHS14*(crLHS27 + crLHS29 + crLHS31);
const double crLHS58 = crLHS26*theta;
const double crLHS59 = crLHS22*crLHS58;
const double crLHS60 = pow(N[1], 2);
const double crLHS61 = N[1]*crLHS5;
const double crLHS62 = N[2]*crLHS58 + N[2]*crLHS61 + crLHS7*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1) + DN(1,2)*DN(2,2));
const double crLHS63 = N[3]*crLHS58 + N[3]*crLHS61 + crLHS7*(DN(1,0)*DN(3,0) + DN(1,1)*DN(3,1) + DN(1,2)*DN(3,2));
const double crLHS64 = crLHS14*(crLHS38 + crLHS40 + crLHS42);
const double crLHS65 = crLHS37*theta;
const double crLHS66 = crLHS22*crLHS65;
const double crLHS67 = pow(N[2], 2);
const double crLHS68 = N[2]*N[3]*crLHS5 + N[3]*crLHS65 + crLHS7*(DN(2,0)*DN(3,0) + DN(2,1)*DN(3,1) + DN(2,2)*DN(3,2));
const double crLHS69 = crLHS14*(crLHS48 + crLHS50 + crLHS52);
const double crLHS70 = crLHS24*crLHS47;
const double crLHS71 = pow(N[3], 2);
rLHS(0,0)+=gauss_weight*(N[0]*crLHS15 + crLHS0*crLHS2 + crLHS0*crLHS5 + crLHS17*crLHS25 + crLHS18*crLHS25 + crLHS19*crLHS25 + crLHS20*crLHS23 + crLHS7*(pow(DN(0,0), 2) + pow(DN(0,1), 2) + pow(DN(0,2), 2)));
rLHS(0,1)+=gauss_weight*(N[1]*crLHS15 + crLHS17*crLHS34 + crLHS18*crLHS34 + crLHS19*crLHS34 + crLHS23*crLHS33 + crLHS36);
rLHS(0,2)+=gauss_weight*(N[2]*crLHS15 + crLHS17*crLHS45 + crLHS18*crLHS45 + crLHS19*crLHS45 + crLHS23*crLHS44 + crLHS46);
rLHS(0,3)+=gauss_weight*(N[3]*crLHS15 + crLHS17*crLHS55 + crLHS18*crLHS55 + crLHS19*crLHS55 + crLHS23*crLHS54 + crLHS56);
rLHS(1,0)+=gauss_weight*(N[0]*crLHS57 + crLHS20*crLHS59 + crLHS25*crLHS28 + crLHS25*crLHS30 + crLHS25*crLHS32 + crLHS36);
rLHS(1,1)+=gauss_weight*(N[1]*crLHS57 + crLHS2*crLHS60 + crLHS28*crLHS34 + crLHS30*crLHS34 + crLHS32*crLHS34 + crLHS33*crLHS59 + crLHS5*crLHS60 + crLHS7*(pow(DN(1,0), 2) + pow(DN(1,1), 2) + pow(DN(1,2), 2)));
rLHS(1,2)+=gauss_weight*(N[2]*crLHS57 + crLHS28*crLHS45 + crLHS30*crLHS45 + crLHS32*crLHS45 + crLHS44*crLHS59 + crLHS62);
rLHS(1,3)+=gauss_weight*(N[3]*crLHS57 + crLHS28*crLHS55 + crLHS30*crLHS55 + crLHS32*crLHS55 + crLHS54*crLHS59 + crLHS63);
rLHS(2,0)+=gauss_weight*(N[0]*crLHS64 + crLHS20*crLHS66 + crLHS25*crLHS39 + crLHS25*crLHS41 + crLHS25*crLHS43 + crLHS46);
rLHS(2,1)+=gauss_weight*(N[1]*crLHS64 + crLHS33*crLHS66 + crLHS34*crLHS39 + crLHS34*crLHS41 + crLHS34*crLHS43 + crLHS62);
rLHS(2,2)+=gauss_weight*(N[2]*crLHS64 + crLHS2*crLHS67 + crLHS39*crLHS45 + crLHS41*crLHS45 + crLHS43*crLHS45 + crLHS44*crLHS66 + crLHS5*crLHS67 + crLHS7*(pow(DN(2,0), 2) + pow(DN(2,1), 2) + pow(DN(2,2), 2)));
rLHS(2,3)+=gauss_weight*(N[3]*crLHS64 + crLHS39*crLHS55 + crLHS41*crLHS55 + crLHS43*crLHS55 + crLHS54*crLHS66 + crLHS68);
rLHS(3,0)+=gauss_weight*(N[0]*crLHS69 + crLHS20*crLHS70 + crLHS25*crLHS49 + crLHS25*crLHS51 + crLHS25*crLHS53 + crLHS56);
rLHS(3,1)+=gauss_weight*(N[1]*crLHS69 + crLHS33*crLHS70 + crLHS34*crLHS49 + crLHS34*crLHS51 + crLHS34*crLHS53 + crLHS63);
rLHS(3,2)+=gauss_weight*(N[2]*crLHS69 + crLHS44*crLHS70 + crLHS45*crLHS49 + crLHS45*crLHS51 + crLHS45*crLHS53 + crLHS68);
rLHS(3,3)+=gauss_weight*(N[3]*crLHS69 + crLHS2*crLHS71 + crLHS49*crLHS55 + crLHS5*crLHS71 + crLHS51*crLHS55 + crLHS53*crLHS55 + crLHS54*crLHS70 + crLHS7*(pow(DN(3,0), 2) + pow(DN(3,1), 2) + pow(DN(3,2), 2)));

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
    
    const double dt = rData.DeltaTime;
    const double theta = rData.TimeIntegrationTheta;
    const double time_coeff = rData.TopOptTimeCoefficient;

    const double dyn_tau = rData.DynamicTau;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    Vector functional_weights = rData.Functional_Weights; //  functional terms weights

    const array_1d<double,3>& t_adj = rData.Temperature_adj;
    const array_1d<double,3>& tn_adj = rData.Temperature_adj_OldStep1;
    const BoundedMatrix<double,2,3>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,2,3> vconv_adj = rData.ConvectionVelocity - rData.MeshVelocity;
    const array_1d<double,3>& source_adj = rData.SourceTerm_adj;
    const array_1d<double,3>& sourcen_adj = rData.SourceTerm_adj_OldStep1;
    const array_1d<double,3>& opt_t = rData.Optimization_Temperature;
    const BoundedMatrix<double,2,3> functional_v = rData.Functional_derivative_velocity;
    const array_1d<double,3>& functional_t = rData.Functional_derivative_transport_scalar;
    const array_1d<double,3>& source = rData.SourceTerm;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT TRANSPORT ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2];
const double crRHS1 = crRHS0*functional_weights[9];
const double crRHS2 = 2.0*functional_weights[4]*(N[0]*opt_t[0] + N[1]*opt_t[1] + N[2]*opt_t[2]);
const double crRHS3 = functional_weights[8]*(N[0]*source[0] + N[1]*source[1] + N[2]*source[2]);
const double crRHS4 = N[0]*functional_t[0] + N[1]*functional_t[1] + N[2]*functional_t[2];
const double crRHS5 = 2.0*crRHS0*crRHS4*functional_weights[7];
const double crRHS6 = DN(0,0)*functional_t[0] + DN(1,0)*functional_t[1] + DN(2,0)*functional_t[2];
const double crRHS7 = DN(0,1)*functional_t[0] + DN(1,1)*functional_t[1] + DN(2,1)*functional_t[2];
const double crRHS8 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2];
const double crRHS9 = 2.0*crRHS8*functional_weights[5];
const double crRHS10 = theta - 1.0;
const double crRHS11 = N[0]*(-crRHS10*sourcen_adj[0] + source_adj[0]*theta) + N[1]*(-crRHS10*sourcen_adj[1] + source_adj[1]*theta) + N[2]*(-crRHS10*sourcen_adj[2] + source_adj[2]*theta);
const double crRHS12 = 1.0*N[0];
const double crRHS13 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2];
const double crRHS14 = crRHS13*time_coeff/dt;
const double crRHS15 = crRHS14*(N[0]*(t_adj[0] - tn_adj[0]) + N[1]*(t_adj[1] - tn_adj[1]) + N[2]*(t_adj[2] - tn_adj[2]));
const double crRHS16 = -crRHS10*tn_adj[0] + t_adj[0]*theta;
const double crRHS17 = -crRHS10*tn_adj[1] + t_adj[1]*theta;
const double crRHS18 = -crRHS10*tn_adj[2] + t_adj[2]*theta;
const double crRHS19 = N[0]*crRHS16 + N[1]*crRHS17 + N[2]*crRHS18;
const double crRHS20 = crRHS0*crRHS19;
const double crRHS21 = N[0]*vconv_adj(0,0) + N[1]*vconv_adj(1,0) + N[2]*vconv_adj(2,0);
const double crRHS22 = DN(0,0)*crRHS21;
const double crRHS23 = N[0]*vconv_adj(0,1) + N[1]*vconv_adj(1,1) + N[2]*vconv_adj(2,1);
const double crRHS24 = DN(0,1)*crRHS23;
const double crRHS25 = crRHS13*crRHS19;
const double crRHS26 = DN(0,0)*crRHS16 + DN(1,0)*crRHS17 + DN(2,0)*crRHS18;
const double crRHS27 = DN(0,1)*crRHS16 + DN(1,1)*crRHS17 + DN(2,1)*crRHS18;
const double crRHS28 = N[0]*functional_v(0,0) + N[1]*functional_v(1,0) + N[2]*functional_v(2,0);
const double crRHS29 = N[0]*functional_v(0,1) + N[1]*functional_v(1,1) + N[2]*functional_v(2,1);
const double crRHS30 = crRHS28*crRHS6 + crRHS29*crRHS7;
const double crRHS31 = crRHS13*functional_weights[6];
const double crRHS32 = (-crRHS1 + crRHS11 + crRHS13*crRHS21*crRHS26 + crRHS13*crRHS23*crRHS27 - crRHS2 - crRHS20 + crRHS3 + crRHS4*functional_weights[6]*(crRHS28*(DN(0,0)*c[0] + DN(1,0)*c[1] + DN(2,0)*c[2]) + crRHS29*(DN(0,1)*c[0] + DN(1,1)*c[1] + DN(2,1)*c[2])) - crRHS5)/(crRHS0*stab_c3 + crRHS14*dyn_tau + crRHS8*stab_c1/pow(h, 2) + stab_c2*sqrt(pow(crRHS21, 2) + pow(crRHS23, 2))*fabs(crRHS13)/h);
const double crRHS33 = crRHS0*crRHS32;
const double crRHS34 = 1.0*crRHS13*crRHS32;
const double crRHS35 = 1.0*crRHS15;
const double crRHS36 = DN(1,0)*crRHS21;
const double crRHS37 = DN(1,1)*crRHS23;
const double crRHS38 = 1.0*crRHS33;
const double crRHS39 = DN(2,0)*crRHS21;
const double crRHS40 = DN(2,1)*crRHS23;
rRHS[0]+=-gauss_weight*(N[0]*crRHS1 - N[0]*crRHS11 + N[0]*crRHS2 + N[0]*crRHS20 - N[0]*crRHS3 + N[0]*crRHS5 + crRHS12*crRHS15 + crRHS12*crRHS33 + crRHS22*crRHS34 + crRHS24*crRHS34 + crRHS25*(crRHS22 + crRHS24) + crRHS31*(N[0]*crRHS30 + crRHS4*(DN(0,0)*crRHS28 + DN(0,1)*crRHS29)) + crRHS8*(DN(0,0)*crRHS26 + DN(0,1)*crRHS27) + crRHS9*(DN(0,0)*crRHS6 + DN(0,1)*crRHS7));
rRHS[1]+=-gauss_weight*(N[1]*crRHS1 - N[1]*crRHS11 + N[1]*crRHS2 + N[1]*crRHS20 - N[1]*crRHS3 + N[1]*crRHS35 + N[1]*crRHS38 + N[1]*crRHS5 + crRHS25*(crRHS36 + crRHS37) + crRHS31*(N[1]*crRHS30 + crRHS4*(DN(1,0)*crRHS28 + DN(1,1)*crRHS29)) + crRHS34*crRHS36 + crRHS34*crRHS37 + crRHS8*(DN(1,0)*crRHS26 + DN(1,1)*crRHS27) + crRHS9*(DN(1,0)*crRHS6 + DN(1,1)*crRHS7));
rRHS[2]+=-gauss_weight*(N[2]*crRHS1 - N[2]*crRHS11 + N[2]*crRHS2 + N[2]*crRHS20 - N[2]*crRHS3 + N[2]*crRHS35 + N[2]*crRHS38 + N[2]*crRHS5 + crRHS25*(crRHS39 + crRHS40) + crRHS31*(N[2]*crRHS30 + crRHS4*(DN(2,0)*crRHS28 + DN(2,1)*crRHS29)) + crRHS34*crRHS39 + crRHS34*crRHS40 + crRHS8*(DN(2,0)*crRHS26 + DN(2,1)*crRHS27) + crRHS9*(DN(2,0)*crRHS6 + DN(2,1)*crRHS7));

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
    
    const double dt = rData.DeltaTime;
    const double theta = rData.TimeIntegrationTheta;
    const double time_coeff = rData.TopOptTimeCoefficient;

    const double dyn_tau = rData.DynamicTau;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,2>& DN = rData.DN_DX;

    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    Vector functional_weights = rData.Functional_Weights; //  functional terms weights

    const array_1d<double,4>& t_adj = rData.Temperature_adj;
    const array_1d<double,4>& tn_adj = rData.Temperature_adj_OldStep1;
    const BoundedMatrix<double,2,4>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,2,4> vconv_adj = rData.ConvectionVelocity - rData.MeshVelocity;
    const array_1d<double,4>& source_adj = rData.SourceTerm_adj;
    const array_1d<double,4>& sourcen_adj = rData.SourceTerm_adj_OldStep1;
    const array_1d<double,4>& opt_t = rData.Optimization_Temperature;
    const BoundedMatrix<double,2,4> functional_v = rData.Functional_derivative_velocity;
    const array_1d<double,4>& functional_t = rData.Functional_derivative_transport_scalar;
    const array_1d<double,4>& source = rData.SourceTerm;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT TRANSPORT ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2] + N[3]*k[3];
const double crRHS1 = crRHS0*functional_weights[9];
const double crRHS2 = 2.0*functional_weights[4]*(N[0]*opt_t[0] + N[1]*opt_t[1] + N[2]*opt_t[2] + N[3]*opt_t[3]);
const double crRHS3 = functional_weights[8]*(N[0]*source[0] + N[1]*source[1] + N[2]*source[2] + N[3]*source[3]);
const double crRHS4 = N[0]*functional_t[0] + N[1]*functional_t[1] + N[2]*functional_t[2] + N[3]*functional_t[3];
const double crRHS5 = 2.0*crRHS0*crRHS4*functional_weights[7];
const double crRHS6 = DN(0,0)*functional_t[0] + DN(1,0)*functional_t[1] + DN(2,0)*functional_t[2] + DN(3,0)*functional_t[3];
const double crRHS7 = DN(0,1)*functional_t[0] + DN(1,1)*functional_t[1] + DN(2,1)*functional_t[2] + DN(3,1)*functional_t[3];
const double crRHS8 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2] + D[3]*N[3];
const double crRHS9 = 2.0*crRHS8*functional_weights[5];
const double crRHS10 = 1.0*N[0];
const double crRHS11 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3];
const double crRHS12 = crRHS11*time_coeff/dt;
const double crRHS13 = crRHS12*(N[0]*(t_adj[0] - tn_adj[0]) + N[1]*(t_adj[1] - tn_adj[1]) + N[2]*(t_adj[2] - tn_adj[2]) + N[3]*(t_adj[3] - tn_adj[3]));
const double crRHS14 = theta - 1.0;
const double crRHS15 = N[0]*(-crRHS14*sourcen_adj[0] + source_adj[0]*theta) + N[1]*(-crRHS14*sourcen_adj[1] + source_adj[1]*theta) + N[2]*(-crRHS14*sourcen_adj[2] + source_adj[2]*theta) + N[3]*(-crRHS14*sourcen_adj[3] + source_adj[3]*theta);
const double crRHS16 = -crRHS14*tn_adj[0] + t_adj[0]*theta;
const double crRHS17 = -crRHS14*tn_adj[1] + t_adj[1]*theta;
const double crRHS18 = -crRHS14*tn_adj[2] + t_adj[2]*theta;
const double crRHS19 = -crRHS14*tn_adj[3] + t_adj[3]*theta;
const double crRHS20 = N[0]*crRHS16 + N[1]*crRHS17 + N[2]*crRHS18 + N[3]*crRHS19;
const double crRHS21 = crRHS0*crRHS20;
const double crRHS22 = N[0]*vconv_adj(0,0) + N[1]*vconv_adj(1,0) + N[2]*vconv_adj(2,0) + N[3]*vconv_adj(3,0);
const double crRHS23 = DN(0,0)*crRHS22;
const double crRHS24 = N[0]*vconv_adj(0,1) + N[1]*vconv_adj(1,1) + N[2]*vconv_adj(2,1) + N[3]*vconv_adj(3,1);
const double crRHS25 = DN(0,1)*crRHS24;
const double crRHS26 = crRHS11*crRHS20;
const double crRHS27 = DN(0,0)*crRHS16 + DN(1,0)*crRHS17 + DN(2,0)*crRHS18 + DN(3,0)*crRHS19;
const double crRHS28 = DN(0,1)*crRHS16 + DN(1,1)*crRHS17 + DN(2,1)*crRHS18 + DN(3,1)*crRHS19;
const double crRHS29 = N[0]*functional_v(0,0) + N[1]*functional_v(1,0) + N[2]*functional_v(2,0) + N[3]*functional_v(3,0);
const double crRHS30 = N[0]*functional_v(0,1) + N[1]*functional_v(1,1) + N[2]*functional_v(2,1) + N[3]*functional_v(3,1);
const double crRHS31 = crRHS29*crRHS6 + crRHS30*crRHS7;
const double crRHS32 = crRHS11*functional_weights[6];
const double crRHS33 = (-crRHS1 + crRHS11*crRHS22*crRHS27 + crRHS11*crRHS24*crRHS28 + crRHS15 - crRHS2 - crRHS21 + crRHS3 + crRHS4*functional_weights[6]*(crRHS29*(DN(0,0)*c[0] + DN(1,0)*c[1] + DN(2,0)*c[2] + DN(3,0)*c[3]) + crRHS30*(DN(0,1)*c[0] + DN(1,1)*c[1] + DN(2,1)*c[2] + DN(3,1)*c[3])) - crRHS5)/(crRHS0*stab_c3 + crRHS12*dyn_tau + crRHS8*stab_c1/pow(h, 2) + stab_c2*sqrt(pow(crRHS22, 2) + pow(crRHS24, 2))*fabs(crRHS11)/h);
const double crRHS34 = crRHS0*crRHS33;
const double crRHS35 = 1.0*crRHS11*crRHS33;
const double crRHS36 = 1.0*crRHS13;
const double crRHS37 = DN(1,0)*crRHS22;
const double crRHS38 = DN(1,1)*crRHS24;
const double crRHS39 = 1.0*crRHS34;
const double crRHS40 = DN(2,0)*crRHS22;
const double crRHS41 = DN(2,1)*crRHS24;
const double crRHS42 = DN(3,0)*crRHS22;
const double crRHS43 = DN(3,1)*crRHS24;
rRHS[0]+=-gauss_weight*(N[0]*crRHS1 - N[0]*crRHS15 + N[0]*crRHS2 + N[0]*crRHS21 - N[0]*crRHS3 + N[0]*crRHS5 + crRHS10*crRHS13 + crRHS10*crRHS34 + crRHS23*crRHS35 + crRHS25*crRHS35 + crRHS26*(crRHS23 + crRHS25) + crRHS32*(N[0]*crRHS31 + crRHS4*(DN(0,0)*crRHS29 + DN(0,1)*crRHS30)) + crRHS8*(DN(0,0)*crRHS27 + DN(0,1)*crRHS28) + crRHS9*(DN(0,0)*crRHS6 + DN(0,1)*crRHS7));
rRHS[1]+=-gauss_weight*(N[1]*crRHS1 - N[1]*crRHS15 + N[1]*crRHS2 + N[1]*crRHS21 - N[1]*crRHS3 + N[1]*crRHS36 + N[1]*crRHS39 + N[1]*crRHS5 + crRHS26*(crRHS37 + crRHS38) + crRHS32*(N[1]*crRHS31 + crRHS4*(DN(1,0)*crRHS29 + DN(1,1)*crRHS30)) + crRHS35*crRHS37 + crRHS35*crRHS38 + crRHS8*(DN(1,0)*crRHS27 + DN(1,1)*crRHS28) + crRHS9*(DN(1,0)*crRHS6 + DN(1,1)*crRHS7));
rRHS[2]+=-gauss_weight*(N[2]*crRHS1 - N[2]*crRHS15 + N[2]*crRHS2 + N[2]*crRHS21 - N[2]*crRHS3 + N[2]*crRHS36 + N[2]*crRHS39 + N[2]*crRHS5 + crRHS26*(crRHS40 + crRHS41) + crRHS32*(N[2]*crRHS31 + crRHS4*(DN(2,0)*crRHS29 + DN(2,1)*crRHS30)) + crRHS35*crRHS40 + crRHS35*crRHS41 + crRHS8*(DN(2,0)*crRHS27 + DN(2,1)*crRHS28) + crRHS9*(DN(2,0)*crRHS6 + DN(2,1)*crRHS7));
rRHS[3]+=-gauss_weight*(N[3]*crRHS1 - N[3]*crRHS15 + N[3]*crRHS2 + N[3]*crRHS21 - N[3]*crRHS3 + N[3]*crRHS36 + N[3]*crRHS39 + N[3]*crRHS5 + crRHS26*(crRHS42 + crRHS43) + crRHS32*(N[3]*crRHS31 + crRHS4*(DN(3,0)*crRHS29 + DN(3,1)*crRHS30)) + crRHS35*crRHS42 + crRHS35*crRHS43 + crRHS8*(DN(3,0)*crRHS27 + DN(3,1)*crRHS28) + crRHS9*(DN(3,0)*crRHS6 + DN(3,1)*crRHS7));

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
    
    const double dt = rData.DeltaTime;
    const double theta = rData.TimeIntegrationTheta;
    const double time_coeff = rData.TopOptTimeCoefficient;

    const double dyn_tau = rData.DynamicTau;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,3>& DN = rData.DN_DX;

    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    Vector functional_weights = rData.Functional_Weights; //  functional terms weights

    const array_1d<double,4>& t_adj = rData.Temperature_adj;
    const array_1d<double,4>& tn_adj = rData.Temperature_adj_OldStep1;
    const BoundedMatrix<double,3,4>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,3,4> vconv_adj = rData.ConvectionVelocity - rData.MeshVelocity;
    const array_1d<double,4>& source_adj = rData.SourceTerm_adj;
    const array_1d<double,4>& sourcen_adj = rData.SourceTerm_adj_OldStep1;
    const array_1d<double,4>& opt_t = rData.Optimization_Temperature;
    const BoundedMatrix<double,3,4> functional_v = rData.Functional_derivative_velocity;
    const array_1d<double,4>& functional_t = rData.Functional_derivative_transport_scalar;
    const array_1d<double,4>& source = rData.SourceTerm;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT TRANSPORT ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*k[0] + N[1]*k[1] + N[2]*k[2] + N[3]*k[3];
const double crRHS1 = crRHS0*functional_weights[9];
const double crRHS2 = 2.0*functional_weights[4]*(N[0]*opt_t[0] + N[1]*opt_t[1] + N[2]*opt_t[2] + N[3]*opt_t[3]);
const double crRHS3 = functional_weights[8]*(N[0]*source[0] + N[1]*source[1] + N[2]*source[2] + N[3]*source[3]);
const double crRHS4 = N[0]*functional_t[0] + N[1]*functional_t[1] + N[2]*functional_t[2] + N[3]*functional_t[3];
const double crRHS5 = 2.0*crRHS0*crRHS4*functional_weights[7];
const double crRHS6 = 1.0*N[0];
const double crRHS7 = N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3];
const double crRHS8 = crRHS7*time_coeff/dt;
const double crRHS9 = crRHS8*(N[0]*(t_adj[0] - tn_adj[0]) + N[1]*(t_adj[1] - tn_adj[1]) + N[2]*(t_adj[2] - tn_adj[2]) + N[3]*(t_adj[3] - tn_adj[3]));
const double crRHS10 = theta - 1.0;
const double crRHS11 = N[0]*(-crRHS10*sourcen_adj[0] + source_adj[0]*theta) + N[1]*(-crRHS10*sourcen_adj[1] + source_adj[1]*theta) + N[2]*(-crRHS10*sourcen_adj[2] + source_adj[2]*theta) + N[3]*(-crRHS10*sourcen_adj[3] + source_adj[3]*theta);
const double crRHS12 = DN(0,0)*functional_t[0] + DN(1,0)*functional_t[1] + DN(2,0)*functional_t[2] + DN(3,0)*functional_t[3];
const double crRHS13 = DN(0,1)*functional_t[0] + DN(1,1)*functional_t[1] + DN(2,1)*functional_t[2] + DN(3,1)*functional_t[3];
const double crRHS14 = DN(0,2)*functional_t[0] + DN(1,2)*functional_t[1] + DN(2,2)*functional_t[2] + DN(3,2)*functional_t[3];
const double crRHS15 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2] + D[3]*N[3];
const double crRHS16 = 2.0*crRHS15*functional_weights[5];
const double crRHS17 = -crRHS10*tn_adj[0] + t_adj[0]*theta;
const double crRHS18 = -crRHS10*tn_adj[1] + t_adj[1]*theta;
const double crRHS19 = -crRHS10*tn_adj[2] + t_adj[2]*theta;
const double crRHS20 = -crRHS10*tn_adj[3] + t_adj[3]*theta;
const double crRHS21 = N[0]*crRHS17 + N[1]*crRHS18 + N[2]*crRHS19 + N[3]*crRHS20;
const double crRHS22 = crRHS0*crRHS21;
const double crRHS23 = N[0]*vconv_adj(0,0) + N[1]*vconv_adj(1,0) + N[2]*vconv_adj(2,0) + N[3]*vconv_adj(3,0);
const double crRHS24 = DN(0,0)*crRHS23;
const double crRHS25 = N[0]*vconv_adj(0,1) + N[1]*vconv_adj(1,1) + N[2]*vconv_adj(2,1) + N[3]*vconv_adj(3,1);
const double crRHS26 = DN(0,1)*crRHS25;
const double crRHS27 = N[0]*vconv_adj(0,2) + N[1]*vconv_adj(1,2) + N[2]*vconv_adj(2,2) + N[3]*vconv_adj(3,2);
const double crRHS28 = DN(0,2)*crRHS27;
const double crRHS29 = crRHS21*crRHS7;
const double crRHS30 = N[0]*functional_v(0,0) + N[1]*functional_v(1,0) + N[2]*functional_v(2,0) + N[3]*functional_v(3,0);
const double crRHS31 = N[0]*functional_v(0,1) + N[1]*functional_v(1,1) + N[2]*functional_v(2,1) + N[3]*functional_v(3,1);
const double crRHS32 = N[0]*functional_v(0,2) + N[1]*functional_v(1,2) + N[2]*functional_v(2,2) + N[3]*functional_v(3,2);
const double crRHS33 = crRHS12*crRHS30 + crRHS13*crRHS31 + crRHS14*crRHS32;
const double crRHS34 = crRHS7*functional_weights[6];
const double crRHS35 = DN(0,0)*crRHS17 + DN(1,0)*crRHS18 + DN(2,0)*crRHS19 + DN(3,0)*crRHS20;
const double crRHS36 = DN(0,1)*crRHS17 + DN(1,1)*crRHS18 + DN(2,1)*crRHS19 + DN(3,1)*crRHS20;
const double crRHS37 = DN(0,2)*crRHS17 + DN(1,2)*crRHS18 + DN(2,2)*crRHS19 + DN(3,2)*crRHS20;
const double crRHS38 = (-crRHS1 + crRHS11 - crRHS2 - crRHS22 + crRHS23*crRHS35*crRHS7 + crRHS25*crRHS36*crRHS7 + crRHS27*crRHS37*crRHS7 + crRHS3 + crRHS4*functional_weights[6]*(crRHS30*(DN(0,0)*c[0] + DN(1,0)*c[1] + DN(2,0)*c[2] + DN(3,0)*c[3]) + crRHS31*(DN(0,1)*c[0] + DN(1,1)*c[1] + DN(2,1)*c[2] + DN(3,1)*c[3]) + crRHS32*(DN(0,2)*c[0] + DN(1,2)*c[1] + DN(2,2)*c[2] + DN(3,2)*c[3])) - crRHS5)/(crRHS0*stab_c3 + crRHS15*stab_c1/pow(h, 2) + crRHS8*dyn_tau + stab_c2*sqrt(pow(crRHS23, 2) + pow(crRHS25, 2) + pow(crRHS27, 2))*fabs(crRHS7)/h);
const double crRHS39 = crRHS0*crRHS38;
const double crRHS40 = 1.0*crRHS38*crRHS7;
const double crRHS41 = 1.0*crRHS9;
const double crRHS42 = DN(1,0)*crRHS23;
const double crRHS43 = DN(1,1)*crRHS25;
const double crRHS44 = DN(1,2)*crRHS27;
const double crRHS45 = 1.0*crRHS39;
const double crRHS46 = DN(2,0)*crRHS23;
const double crRHS47 = DN(2,1)*crRHS25;
const double crRHS48 = DN(2,2)*crRHS27;
const double crRHS49 = DN(3,0)*crRHS23;
const double crRHS50 = DN(3,1)*crRHS25;
const double crRHS51 = DN(3,2)*crRHS27;
rRHS[0]+=-gauss_weight*(N[0]*crRHS1 - N[0]*crRHS11 + N[0]*crRHS2 + N[0]*crRHS22 - N[0]*crRHS3 + N[0]*crRHS5 + crRHS15*(DN(0,0)*crRHS35 + DN(0,1)*crRHS36 + DN(0,2)*crRHS37) + crRHS16*(DN(0,0)*crRHS12 + DN(0,1)*crRHS13 + DN(0,2)*crRHS14) + crRHS24*crRHS40 + crRHS26*crRHS40 + crRHS28*crRHS40 + crRHS29*(crRHS24 + crRHS26 + crRHS28) + crRHS34*(N[0]*crRHS33 + crRHS4*(DN(0,0)*crRHS30 + DN(0,1)*crRHS31 + DN(0,2)*crRHS32)) + crRHS39*crRHS6 + crRHS6*crRHS9);
rRHS[1]+=-gauss_weight*(N[1]*crRHS1 - N[1]*crRHS11 + N[1]*crRHS2 + N[1]*crRHS22 - N[1]*crRHS3 + N[1]*crRHS41 + N[1]*crRHS45 + N[1]*crRHS5 + crRHS15*(DN(1,0)*crRHS35 + DN(1,1)*crRHS36 + DN(1,2)*crRHS37) + crRHS16*(DN(1,0)*crRHS12 + DN(1,1)*crRHS13 + DN(1,2)*crRHS14) + crRHS29*(crRHS42 + crRHS43 + crRHS44) + crRHS34*(N[1]*crRHS33 + crRHS4*(DN(1,0)*crRHS30 + DN(1,1)*crRHS31 + DN(1,2)*crRHS32)) + crRHS40*crRHS42 + crRHS40*crRHS43 + crRHS40*crRHS44);
rRHS[2]+=-gauss_weight*(N[2]*crRHS1 - N[2]*crRHS11 + N[2]*crRHS2 + N[2]*crRHS22 - N[2]*crRHS3 + N[2]*crRHS41 + N[2]*crRHS45 + N[2]*crRHS5 + crRHS15*(DN(2,0)*crRHS35 + DN(2,1)*crRHS36 + DN(2,2)*crRHS37) + crRHS16*(DN(2,0)*crRHS12 + DN(2,1)*crRHS13 + DN(2,2)*crRHS14) + crRHS29*(crRHS46 + crRHS47 + crRHS48) + crRHS34*(N[2]*crRHS33 + crRHS4*(DN(2,0)*crRHS30 + DN(2,1)*crRHS31 + DN(2,2)*crRHS32)) + crRHS40*crRHS46 + crRHS40*crRHS47 + crRHS40*crRHS48);
rRHS[3]+=-gauss_weight*(N[3]*crRHS1 - N[3]*crRHS11 + N[3]*crRHS2 + N[3]*crRHS22 - N[3]*crRHS3 + N[3]*crRHS41 + N[3]*crRHS45 + N[3]*crRHS5 + crRHS15*(DN(3,0)*crRHS35 + DN(3,1)*crRHS36 + DN(3,2)*crRHS37) + crRHS16*(DN(3,0)*crRHS12 + DN(3,1)*crRHS13 + DN(3,2)*crRHS14) + crRHS29*(crRHS49 + crRHS50 + crRHS51) + crRHS34*(N[3]*crRHS33 + crRHS4*(DN(3,0)*crRHS30 + DN(3,1)*crRHS31 + DN(3,2)*crRHS32)) + crRHS40*crRHS49 + crRHS40*crRHS50 + crRHS40*crRHS51);
 
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