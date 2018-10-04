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

#include "fluid_element.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"

#include "custom_utilities/qsvms_data.h"
#include "custom_utilities/time_integrated_qsvms_data.h"
#include "custom_utilities/fic_data.h"
#include "custom_utilities/time_integrated_fic_data.h"
#include "custom_utilities/symbolic_navier_stokes_data.h"
#include "custom_utilities/two_fluid_navier_stokes_data.h"
#include "custom_utilities/element_size_calculator.h"
#include "custom_utilities/vorticity_utilities.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TElementData >
FluidElement<TElementData>::FluidElement(IndexType NewId):
    Element(NewId)
{}

template< class TElementData >
FluidElement<TElementData>::FluidElement(IndexType NewId, const NodesArrayType& ThisNodes):
    Element(NewId,ThisNodes)
{}


template< class TElementData >
FluidElement<TElementData>::FluidElement(IndexType NewId, GeometryType::Pointer pGeometry):
    Element(NewId,pGeometry)
{}


template< class TElementData >
FluidElement<TElementData>::FluidElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
    Element(NewId,pGeometry,pProperties)
{}


template< class TElementData >
FluidElement<TElementData>::~FluidElement()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TElementData >
Element::Pointer FluidElement<TElementData>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY;
    KRATOS_ERROR << "Attempting to Create base FluidElement instances." << std::endl;
    KRATOS_CATCH("");
}

template< class TElementData >
Element::Pointer FluidElement<TElementData>::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY;
    KRATOS_ERROR << "Attempting to Create base FluidElement instances." << std::endl;
    KRATOS_CATCH("");
}

template< class TElementData >
void FluidElement<TElementData>::Initialize() {
    KRATOS_TRY;

    // If we are restarting, the constitutive law will be already defined
    if (mpConstitutiveLaw == nullptr) {
        const Properties& r_properties = this->GetProperties();
        KRATOS_ERROR_IF_NOT(r_properties.Has(CONSTITUTIVE_LAW))
            << "In initialization of Element " << this->Info()
            << ": No CONSTITUTIVE_LAW defined for property "
            << r_properties.Id() << "." << std::endl;

        mpConstitutiveLaw = r_properties[CONSTITUTIVE_LAW]->Clone();

        const GeometryType& r_geometry = this->GetGeometry();
        const auto& r_shape_functions = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
        mpConstitutiveLaw->InitializeMaterial(r_properties,r_geometry,row(r_shape_functions,0));
    }

    KRATOS_CATCH("");
}

template <class TElementData>
void FluidElement<TElementData>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                      VectorType& rRightHandSideVector,
                                                      ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    if (TElementData::ElementManagesTimeIntegration) {
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

            data.UpdateGeometryValues(g, gauss_weights[g], row(shape_functions, g),
                shape_derivatives[g]);

            this->CalculateMaterialResponse(data);

            this->AddTimeIntegratedSystem(
                data, rLeftHandSideMatrix, rRightHandSideVector);
        }
    }
}

template <class TElementData>
void FluidElement<TElementData>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                                       ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);

    if (TElementData::ElementManagesTimeIntegration) {
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
            data.UpdateGeometryValues(g, gauss_weights[g], row(shape_functions, g),
                shape_derivatives[g]);

            this->CalculateMaterialResponse(data);

            this->AddTimeIntegratedLHS(data, rLeftHandSideMatrix);
        }
    }
}

template <class TElementData>
void FluidElement<TElementData>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                        ProcessInfo& rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);

    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    if (TElementData::ElementManagesTimeIntegration) {
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
            data.UpdateGeometryValues(g, gauss_weights[g], row(shape_functions, g),
                shape_derivatives[g]);

            this->CalculateMaterialResponse(data);

            this->AddTimeIntegratedRHS(data, rRightHandSideVector);
        }
    }
}

template <class TElementData>
void FluidElement<TElementData>::CalculateLocalVelocityContribution(
    MatrixType& rDampMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
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
            data.UpdateGeometryValues(
                g, gauss_weights[g], row(shape_functions, g), r_dndx);

            this->CalculateMaterialResponse(data);

            this->AddVelocitySystem(data, rDampMatrix, rRightHandSideVector);
        }
    }
}

template <class TElementData>
void FluidElement<TElementData>::CalculateMassMatrix(MatrixType& rMassMatrix,
                                                     ProcessInfo& rCurrentProcessInfo)
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
            data.UpdateGeometryValues(g, gauss_weights[g], row(shape_functions, g),
                shape_derivatives[g]);

            this->CalculateMaterialResponse(data);

            this->AddMassLHS(data, rMassMatrix);
        }
    }
}

template< class TElementData >
void FluidElement< TElementData >::EquationIdVector(EquationIdVectorType &rResult, ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& r_geometry = this->GetGeometry();

    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    const unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        rResult[LocalIndex++] = r_geometry[i].GetDof(VELOCITY_X,xpos).EquationId();
        rResult[LocalIndex++] = r_geometry[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
        if (Dim == 3) rResult[LocalIndex++] = r_geometry[i].GetDof(VELOCITY_Z,xpos+2).EquationId();
        rResult[LocalIndex++] = r_geometry[i].GetDof(PRESSURE,ppos).EquationId();
    }
}


template< class TElementData >
void FluidElement< TElementData >::GetDofList(DofsVectorType &rElementalDofList, ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& r_geometry = this->GetGeometry();

     if (rElementalDofList.size() != LocalSize)
         rElementalDofList.resize(LocalSize);

     unsigned int LocalIndex = 0;

     for (unsigned int i = 0; i < NumNodes; ++i)
     {
         rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(VELOCITY_X);
         rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(VELOCITY_Y);
         if (Dim == 3) rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(VELOCITY_Z);
         rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(PRESSURE);
     }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FluidElement<TElementData>::GetFirstDerivativesVector(Vector &rValues, int Step)
{
    GeometryType& r_geometry = this->GetGeometry();

    if (rValues.size() != LocalSize)
        rValues.resize(LocalSize,false);

    unsigned int Index = 0;

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        const array_1d<double,3>& rVel = r_geometry[i].FastGetSolutionStepValue(VELOCITY,Step);
        for (unsigned int d = 0; d < Dim; d++)
            rValues[Index++] = rVel[d];
        rValues[Index++] = r_geometry[i].FastGetSolutionStepValue(PRESSURE,Step);
    }
}


template< class TElementData >
void FluidElement<TElementData>::GetSecondDerivativesVector(Vector &rValues, int Step)
{
    GeometryType& r_geometry = this->GetGeometry();

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
GeometryData::IntegrationMethod FluidElement<TElementData>::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_2;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template< class TElementData >
int FluidElement<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo)
{
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

    // Extra variables used in computing projections
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);

    const GeometryType& r_geometry = this->GetGeometry();

    for(unsigned int i=0; i<NumNodes; ++i)
    {
        const Node<3>& rNode = r_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rNode);

        // Check that required dofs exist
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X,rNode);
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y,rNode);
        if (Dim == 3) KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Z,rNode);
        KRATOS_CHECK_DOF_IN_NODE(PRESSURE,rNode);
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
    KRATOS_ERROR_IF(mpConstitutiveLaw == nullptr) << "Constitutive Law not initialized for Element " << this->Info() << std::endl;

    constexpr auto dimension = Dim;  // I need to set this here otherwise it gives me a linking error when attempting to '<<' Dim.

    KRATOS_ERROR_IF(mpConstitutiveLaw->WorkingSpaceDimension() != Dim)
        << "Wrong dimension: The " << mpConstitutiveLaw->WorkingSpaceDimension()
        << "D constitutive law " << mpConstitutiveLaw->Info()
        << " is not compatible with " << dimension << "D element " << this->Info()
        << "." << std::endl;

    out = mpConstitutiveLaw->Check(this->GetProperties(),r_geometry,rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT( out == 0) << "The Constitutive Law provided for Element " << this->Info() << " is not correct." << std::endl;

    return out;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FluidElement<TElementData>::GetValueOnIntegrationPoints(
    Variable<array_1d<double, 3 > > const& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    if (rVariable == VORTICITY)
    {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_function_gradients;
        this->CalculateGeometryData(gauss_weights,shape_functions,shape_function_gradients);

        VorticityUtilities<Dim>::CalculateVorticityVector(this->GetGeometry(),shape_function_gradients,rValues);
    }
}


template< class TElementData >
void FluidElement<TElementData>::GetValueOnIntegrationPoints(
    Variable<double> const& rVariable,
    std::vector<double>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    if (rVariable == Q_VALUE)
    {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_function_gradients;
        this->CalculateGeometryData(gauss_weights,shape_functions,shape_function_gradients);

        VorticityUtilities<Dim>::CalculateQValue(this->GetGeometry(),shape_function_gradients,rValues);
    }
    else if (rVariable == VORTICITY_MAGNITUDE)
    {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_function_gradients;
        this->CalculateGeometryData(gauss_weights,shape_functions,shape_function_gradients);

        VorticityUtilities<Dim>::CalculateVorticityMagnitude(this->GetGeometry(),shape_function_gradients,rValues);
    }
    else if (rVariable == UPDATE_STATISTICS)
    {
        KRATOS_DEBUG_ERROR_IF(!rCurrentProcessInfo.Has(STATISTICS_CONTAINER)) << "Trying to compute turbulent statistics, but ProcessInfo does not have STATISTICS_CONTAINER defined." << std::endl;
        rCurrentProcessInfo.GetValue(STATISTICS_CONTAINER)->UpdateStatistics(this);
    }
}

template <class TElementData>
void FluidElement<TElementData>::GetValueOnIntegrationPoints(
    Variable<array_1d<double, 6>> const& rVariable,
    std::vector<array_1d<double, 6>>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{}

template <class TElementData>
void FluidElement<TElementData>::GetValueOnIntegrationPoints(
    Variable<Vector> const& rVariable,
    std::vector<Vector>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{}

template <class TElementData>
void FluidElement<TElementData>::GetValueOnIntegrationPoints(
    Variable<Matrix> const& rVariable,
    std::vector<Matrix>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output


template< class TElementData >
std::string FluidElement<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "FluidElement #" << Id();
    return buffer.str();
}


template< class TElementData >
void FluidElement<TElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "FluidElement" << Dim << "D";
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template <class TElementData>
double FluidElement<TElementData>::GetAtCoordinate(
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
array_1d<double, 3> FluidElement<TElementData>::GetAtCoordinate(
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
double FluidElement<TElementData>::GetAtCoordinate(
    const double Value,
    const typename TElementData::ShapeFunctionsType& rN) const
{
    return Value;
}

template <class TElementData>
void FluidElement<TElementData>::CalculateMaterialResponse(TElementData& rData) const {

    Internals::StrainRateSpecialization<TElementData,Dim>::Calculate(rData.StrainRate,rData.Velocity,rData.DN_DX);

    auto& Values = rData.ConstitutiveLawValues;

    const Vector shape_functions_vector(rData.N);
    Values.SetShapeFunctionsValues(shape_functions_vector);

    //ATTENTION: here we assume that only one constitutive law is employed for all of the gauss points in the element.
    //this is ok under the hypothesis that no history dependent behavior is employed
    mpConstitutiveLaw->CalculateMaterialResponseCauchy(Values);

    mpConstitutiveLaw->CalculateValue(Values,EFFECTIVE_VISCOSITY,rData.EffectiveViscosity);
}

template< class TElementData >
void FluidElement<TElementData>::CalculateGeometryData(Vector &rGaussWeights,
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

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FluidElement<TElementData>::ConvectionOperator(Vector &rResult,
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
void FluidElement<TElementData>::AddTimeIntegratedSystem(
    TElementData& rData, MatrixType& rLHS, VectorType& rRHS) {
    KRATOS_TRY;

    KRATOS_ERROR << "Calling base FluidElement::AddTimeIntegratedSystem "
                    "implementation. This method is not supported by your "
                    "element."
                 << std::endl;

    KRATOS_CATCH("");
}

template <class TElementData>
void FluidElement<TElementData>::AddTimeIntegratedLHS(
    TElementData& rData, MatrixType& rLHS) {
    KRATOS_TRY;

    KRATOS_ERROR << "Calling base FluidElement::AddTimeIntegratedLHS "
                    "implementation. This method is not supported by your "
                    "element."
                 << std::endl;

    KRATOS_CATCH("");
}

template <class TElementData>
void FluidElement<TElementData>::AddTimeIntegratedRHS(
    TElementData& rData, VectorType& rRHS) {
    KRATOS_TRY;

    KRATOS_ERROR << "Calling base FluidElement::AddTimeIntegratedRHS "
                    "implementation. This method is not supported by your "
                    "element."
                 << std::endl;

    KRATOS_CATCH("");
}

template <class TElementData>
void FluidElement<TElementData>::AddVelocitySystem(
    TElementData& rData,
    MatrixType& rLocalLHS,
    VectorType& rLocalRHS) {
    KRATOS_TRY;

    KRATOS_ERROR << "Calling base FluidElement::AddVelocitySystem "
                    "implementation. This method is not supported by your "
                    "element."
                 << std::endl;

    KRATOS_CATCH("");
}

template <class TElementData>
void FluidElement<TElementData>::AddMassLHS(
    TElementData& rData, MatrixType& rMassMatrix) {
    KRATOS_TRY;

    KRATOS_ERROR << "Calling base FluidElement::AddMassLHS "
                    "implementation. This method is not supported by your "
                    "element."
                 << std::endl;

    KRATOS_CATCH("");
}

template <class TElementData>
void FluidElement<TElementData>::AddBoundaryIntegral(TElementData& rData,
    const Vector& rUnitNormal, MatrixType& rLHS, VectorType& rRHS) {

    KRATOS_TRY;

    KRATOS_ERROR << "Calling base FluidElement::AddBoundaryIntegral "
                    "implementation. This method is not supported by your "
                    "element."
                 << std::endl;

    KRATOS_CATCH("");
}

template <class TElementData>
void FluidElement<TElementData>::GetCurrentValuesVector(
    const TElementData& rData,
    array_1d<double,LocalSize>& rValues) const {

    int local_index = 0;

    const auto& r_velocities = rData.Velocity;
    const auto& r_pressures = rData.Pressure;

    for (unsigned int i = 0; i < NumNodes; ++i) {
        for (unsigned int d = 0; d < Dim; ++d)  // Velocity Dofs
            rValues[local_index++] = r_velocities(i, d);
        rValues[local_index++] = r_pressures[i];  // Pressure Dof
    }
}

template< class TElementData >
ConstitutiveLaw::Pointer FluidElement<TElementData>::GetConstitutiveLaw() {
    return this->mpConstitutiveLaw;
}

template< class TElementData >
const ConstitutiveLaw::Pointer FluidElement<TElementData>::GetConstitutiveLaw() const {
    return this->mpConstitutiveLaw;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template< class TElementData >
void FluidElement<TElementData>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    rSerializer.save("mpConstitutiveLaw",this->mpConstitutiveLaw);
}


template< class TElementData >
void FluidElement<TElementData>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    rSerializer.load("mpConstitutiveLaw",this->mpConstitutiveLaw);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

namespace Internals {

template< class TElementData >
void StrainRateSpecialization<TElementData,2>::Calculate(
    Vector& rStrainRate,
    const typename TElementData::NodalVectorData& rVelocities,
    const typename TElementData::ShapeDerivativesType& rDNDX) {

    noalias(rStrainRate) = ZeroVector(3);
    for (unsigned int i = 0; i < TElementData::NumNodes; i++) {
        rStrainRate[0] += rDNDX(i,0)*rVelocities(i,0);
        rStrainRate[1] += rDNDX(i,1)*rVelocities(i,1);
        rStrainRate[2] += rDNDX(i,0)*rVelocities(i,1) + rDNDX(i,1)*rVelocities(i,0);
    }
}

template< class TElementData >
void StrainRateSpecialization<TElementData,3>::Calculate(
    Vector& rStrainRate,
    const typename TElementData::NodalVectorData& rVelocities,
    const typename TElementData::ShapeDerivativesType& rDNDX) {

    noalias(rStrainRate) = ZeroVector(6);
    for (unsigned int i = 0; i < TElementData::NumNodes; i++) {
        rStrainRate[0] += rDNDX(i,0)*rVelocities(i,0);
        rStrainRate[1] += rDNDX(i,1)*rVelocities(i,1);
        rStrainRate[2] += rDNDX(i,2)*rVelocities(i,2);
        rStrainRate[3] += rDNDX(i,0)*rVelocities(i,1) + rDNDX(i,1)*rVelocities(i,0);
        rStrainRate[4] += rDNDX(i,1)*rVelocities(i,2) + rDNDX(i,2)*rVelocities(i,1);
        rStrainRate[5] += rDNDX(i,0)*rVelocities(i,2) + rDNDX(i,2)*rVelocities(i,0);
    }
}

} // namespace Internals

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class FluidElement< SymbolicNavierStokesData<2,3> >;
template class FluidElement< SymbolicNavierStokesData<3,4> >;

template class FluidElement< QSVMSData<2,3> >;
template class FluidElement< QSVMSData<3,4> >;

template class FluidElement< QSVMSData<2,4> >;
template class FluidElement< QSVMSData<3,8> >;

template class FluidElement< TimeIntegratedQSVMSData<2,3> >;
template class FluidElement< TimeIntegratedQSVMSData<3,4> >;

template class FluidElement< FICData<2,3> >;
template class FluidElement< FICData<3,4> >;

template class FluidElement< FICData<2,4> >;
template class FluidElement< FICData<3,8> >;

template class FluidElement< TimeIntegratedFICData<2,3> >;
template class FluidElement< TimeIntegratedFICData<3,4> >;

template class FluidElement< TwoFluidNavierStokesData<2, 3> >;
template class FluidElement< TwoFluidNavierStokesData<3, 4> >;

///////////////////////////////////////////////////////////////////////////////////////////////////

}
