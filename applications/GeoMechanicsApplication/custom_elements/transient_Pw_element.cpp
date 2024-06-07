// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

// Application includes
#include "custom_elements/transient_Pw_element.hpp"
#include "custom_utilities/dof_utilities.h"
#include "custom_utilities/transport_equation_utilities.hpp"

namespace Kratos
{

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer TransientPwElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                             NodesArrayType const& ThisNodes,
                                                             PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new TransientPwElement(NewId, this->GetGeometry().Create(ThisNodes),
                                                   pProperties, this->GetStressStatePolicy().Clone()));
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer TransientPwElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                             GeometryType::Pointer pGeom,
                                                             PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(
        new TransientPwElement(NewId, pGeom, pProperties, this->GetStressStatePolicy().Clone()));
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const
{
    rElementalDofList = GetDofs();
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix,
                                                              const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int N_DOF = this->GetNumberOfDOF();

    // Resizing mass matrix
    if (rMassMatrix.size1() != N_DOF) rMassMatrix.resize(N_DOF, N_DOF, false);
    noalias(rMassMatrix) = ZeroMatrix(N_DOF, N_DOF);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateDampingMatrix(MatrixType& rDampingMatrix,
                                                                 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int N_DOF = this->GetNumberOfDOF();

    // Compute Damping Matrix
    if (rDampingMatrix.size1() != N_DOF) rDampingMatrix.resize(N_DOF, N_DOF, false);
    noalias(rDampingMatrix) = ZeroMatrix(N_DOF, N_DOF);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::GetValuesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY

    const unsigned int N_DOF = this->GetNumberOfDOF();

    if (rValues.size() != N_DOF) rValues.resize(N_DOF, false);

    // Why are we constructing a zero vector here?
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rValues[i] = 0.0;
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::GetFirstDerivativesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY

    const unsigned int N_DOF = this->GetNumberOfDOF();

    if (rValues.size() != N_DOF) rValues.resize(N_DOF, false);

    // Why are we constructing a zero vector here?
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rValues[i] = 0.0;
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::GetSecondDerivativesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY

    const unsigned int N_DOF = this->GetNumberOfDOF();

    if (rValues.size() != N_DOF) rValues.resize(N_DOF, false);

    // Why are we constructing a zero vector here?
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rValues[i] = 0.0;
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const PropertiesType& Prop       = this->GetProperties();
    const GeometryType&   Geom       = this->GetGeometry();
    const unsigned int    NumGPoints = Geom.IntegrationPointsNumber(this->GetIntegrationMethod());

    // pointer to constitutive laws
    if (mConstitutiveLawVector.size() != NumGPoints) mConstitutiveLawVector.resize(NumGPoints);

    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i) {
        mConstitutiveLawVector[i] = nullptr;
    }

    if (mRetentionLawVector.size() != NumGPoints) mRetentionLawVector.resize(NumGPoints);
    for (unsigned int i = 0; i < mRetentionLawVector.size(); ++i) {
        mRetentionLawVector[i] = RetentionLawFactory::Clone(Prop);
        mRetentionLawVector[i]->InitializeMaterial(
            Prop, Geom, row(Geom.ShapeFunctionsValues(this->GetIntegrationMethod()), i));
    }

    mIsInitialised = true;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
int TransientPwElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const PropertiesType& Prop = this->GetProperties();
    const GeometryType&   Geom = this->GetGeometry();

    if (Geom.DomainSize() < 1.0e-15)
        KRATOS_ERROR << "DomainSize < 1.0e-15 for the element " << this->Id() << std::endl;

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        if (Geom[i].SolutionStepsDataHas(WATER_PRESSURE) == false)
            KRATOS_ERROR << "missing variable WATER_PRESSURE on node " << Geom[i].Id() << std::endl;

        if (Geom[i].SolutionStepsDataHas(DT_WATER_PRESSURE) == false)
            KRATOS_ERROR << "missing variable DT_WATER_PRESSURE on node " << Geom[i].Id() << std::endl;

        if (Geom[i].SolutionStepsDataHas(VOLUME_ACCELERATION) == false)
            KRATOS_ERROR << "missing variable VOLUME_ACCELERATION on node " << Geom[i].Id() << std::endl;

        if (Geom[i].HasDofFor(WATER_PRESSURE) == false)
            KRATOS_ERROR << "missing variable WATER_PRESSURE on node " << Geom[i].Id() << std::endl;
    }

    // Verify ProcessInfo variables

    // Verify properties
    if (Prop.Has(DENSITY_WATER) == false || Prop[DENSITY_WATER] < 0.0)
        KRATOS_ERROR << "DENSITY_WATER does not exist in the material "
                        "properties or has an invalid value at element"
                     << this->Id() << std::endl;

    if (Prop.Has(BULK_MODULUS_SOLID) == false || Prop[BULK_MODULUS_SOLID] < 0.0)
        KRATOS_ERROR << "BULK_MODULUS_SOLID does not exist in the material "
                        "properties or has an invalid value at element"
                     << this->Id() << std::endl;

    if (Prop.Has(POROSITY) == false || Prop[POROSITY] < 0.0 || Prop[POROSITY] > 1.0)
        KRATOS_ERROR << "POROSITY does not exist in the material properties or "
                        "has an invalid value at element"
                     << this->Id() << std::endl;

    if (TDim == 2) {
        // If this is a 2D problem, nodes must be in XY plane
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            if (Geom[i].Z() != 0.0)
                KRATOS_ERROR << " Node with non-zero Z coordinate found. Id: " << Geom[i].Id() << std::endl;
        }
    }

    // Verify specific properties
    if (Prop.Has(BULK_MODULUS_FLUID) == false || Prop[BULK_MODULUS_FLUID] < 0.0)
        KRATOS_ERROR << "BULK_MODULUS_FLUID does not exist in the material "
                        "properties or has an invalid value at element"
                     << this->Id() << std::endl;

    if (Prop.Has(DYNAMIC_VISCOSITY) == false || Prop[DYNAMIC_VISCOSITY] < 0.0)
        KRATOS_ERROR << "DYNAMIC_VISCOSITY does not exist in the material "
                        "properties or has an invalid value at element"
                     << this->Id() << std::endl;

    if (Prop.Has(PERMEABILITY_XX) == false || Prop[PERMEABILITY_XX] < 0.0)
        KRATOS_ERROR << "PERMEABILITY_XX does not exist in the material "
                        "properties or has an invalid value at element"
                     << this->Id() << std::endl;

    if (Prop.Has(PERMEABILITY_YY) == false || Prop[PERMEABILITY_YY] < 0.0)
        KRATOS_ERROR << "PERMEABILITY_YY does not exist in the material "
                        "properties or has an invalid value at element"
                     << this->Id() << std::endl;

    if (Prop.Has(PERMEABILITY_XY) == false || Prop[PERMEABILITY_XY] < 0.0)
        KRATOS_ERROR << "PERMEABILITY_XY does not exist in the material "
                        "properties or has an invalid value at element"
                     << this->Id() << std::endl;

    if (!Prop.Has(BIOT_COEFFICIENT))
        KRATOS_ERROR << "BIOT_COEFFICIENT does not exist in the material "
                        "properties in element"
                     << this->Id() << std::endl;

    if constexpr (TDim > 2) {
        if (Prop.Has(PERMEABILITY_ZZ) == false || Prop[PERMEABILITY_ZZ] < 0.0)
            KRATOS_ERROR << "PERMEABILITY_ZZ does not exist in the material "
                            "properties or has an invalid value at element"
                         << this->Id() << std::endl;

        if (Prop.Has(PERMEABILITY_YZ) == false || Prop[PERMEABILITY_YZ] < 0.0)
            KRATOS_ERROR << "PERMEABILITY_YZ does not exist in the material "
                            "properties or has an invalid value at element"
                         << this->Id() << std::endl;

        if (Prop.Has(PERMEABILITY_ZX) == false || Prop[PERMEABILITY_ZX] < 0.0)
            KRATOS_ERROR << "PERMEABILITY_ZX does not exist in the material "
                            "properties or has an invalid value at element"
                         << this->Id() << std::endl;
    }

    // Verify that the constitutive law has the correct dimension

    // Check constitutive law
    if (mRetentionLawVector.size() > 0) {
        return mRetentionLawVector[0]->Check(Prop, rCurrentProcessInfo);
    }

    return 0;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (!mIsInitialised) this->Initialize(rCurrentProcessInfo);

    // Defining necessary variables
    const GeometryType& Geom       = this->GetGeometry();
    const unsigned int  NumGPoints = Geom.IntegrationPointsNumber(this->GetIntegrationMethod());

    RetentionLaw::Parameters RetentionParameters(this->GetProperties());

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        // Initialize retention law
        mRetentionLawVector[GPoint]->InitializeSolutionStep(RetentionParameters);
    }

    // reset hydraulic discharge
    this->ResetHydraulicDischarge();

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::InitializeNonLinearIteration(const ProcessInfo&)
{
    // nothing
}

//----------------------------------------------------------------------------------------------------

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::FinalizeNonLinearIteration(const ProcessInfo&)
{
    // nothing
}

//----------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    this->CalculateHydraulicDischarge(rCurrentProcessInfo);

    // Defining necessary variables
    const GeometryType& Geom       = this->GetGeometry();
    const unsigned int  NumGPoints = Geom.IntegrationPointsNumber(this->GetIntegrationMethod());

    RetentionLaw::Parameters RetentionParameters(this->GetProperties());

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        // retention law
        mRetentionLawVector[GPoint]->FinalizeSolutionStep(RetentionParameters);
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                                                       std::vector<double>& rOutput,
                                                                       const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == DEGREE_OF_SATURATION || rVariable == EFFECTIVE_SATURATION ||
        rVariable == BISHOP_COEFFICIENT || rVariable == DERIVATIVE_OF_SATURATION ||
        rVariable == RELATIVE_PERMEABILITY || rVariable == HYDRAULIC_HEAD) {
        BaseType::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    } else {
        if (rOutput.size() != mRetentionLawVector.size())
            rOutput.resize(mRetentionLawVector.size());

        std::fill(rOutput.begin(), rOutput.end(), 0.0);
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                                                       std::vector<array_1d<double, 3>>& rOutput,
                                                                       const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& rGeom      = this->GetGeometry();
    const IndexType     NumGPoints = rGeom.IntegrationPointsNumber(this->GetIntegrationMethod());
    if (rOutput.size() != NumGPoints) rOutput.resize(NumGPoints);

    if (rVariable == FLUID_FLUX_VECTOR) {
        std::vector<double> permeability_update_factors(NumGPoints, 1.0);
        const auto fluid_fluxes = this->CalculateFluidFluxes(permeability_update_factors, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            GeoElementUtilities::FillArray1dOutput(rOutput[GPoint], fluid_fluxes[GPoint]);
        }
    } else {
        if (rOutput.size() != mRetentionLawVector.size())
            rOutput.resize(mRetentionLawVector.size());

        for (unsigned int i = 0; i < mRetentionLawVector.size(); ++i) {
            noalias(rOutput[i]) = ZeroVector(3);
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                                       std::vector<Matrix>& rOutput,
                                                                       const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == PERMEABILITY_MATRIX) {
        BaseType::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    } else {
        if (rOutput.size() != mRetentionLawVector.size())
            rOutput.resize(mRetentionLawVector.size());

        for (unsigned int i = 0; i < mRetentionLawVector.size(); ++i) {
            rOutput[i].resize(TDim, TDim, false);
            noalias(rOutput[i]) = ZeroMatrix(TDim, TDim);
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateAll(MatrixType&        rLeftHandSideMatrix,
                                                       VectorType&        rRightHandSideVector,
                                                       const ProcessInfo& rCurrentProcessInfo,
                                                       bool CalculateStiffnessMatrixFlag,
                                                       bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    // Previous definitions
    const PropertiesType&                           Prop = this->GetProperties();
    const GeometryType&                             Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        Geom.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int NumGPoints = IntegrationPoints.size();

    // Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    RetentionLaw::Parameters RetentionParameters(this->GetProperties());

    const auto fluid_pressures = GeoTransportEquationUtilities::CalculateFluidPressures(
        Variables.NContainer, Variables.PressureVector);
    const auto relative_permeability_values = this->CalculateRelativePermeabilityValues(fluid_pressures);
    const auto bishop_coefficients = this->CalculateBishopCoefficients(fluid_pressures);
    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(IntegrationPoints, Variables.detJContainer);
    std::vector<double> biot_coefficients(NumGPoints, Prop[BIOT_COEFFICIENT]);
    const auto          degrees_of_saturation = this->CalculateDegreesOfSaturation(fluid_pressures);
    const auto derivatives_of_saturation = this->CalculateDerivativesOfSaturation(fluid_pressures);
    const auto biot_moduli_inverse = GeoTransportEquationUtilities::CalculateInverseBiotModuli(
        biot_coefficients, degrees_of_saturation, derivatives_of_saturation, Prop);

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        // Compute GradNpT, B and StrainVector
        this->CalculateKinematics(Variables, GPoint);

        // Compute Nu and BodyAcceleration
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Variables.Nu, Variables.NContainer, GPoint);
        GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
            Variables.BodyAcceleration, Variables.NContainer, Variables.VolumeAcceleration, GPoint);

        Variables.RelativePermeability = relative_permeability_values[GPoint];
        Variables.BishopCoefficient    = bishop_coefficients[GPoint];

        Variables.BiotCoefficient    = biot_coefficients[GPoint];
        Variables.BiotModulusInverse = biot_moduli_inverse[GPoint];
        Variables.DegreeOfSaturation = degrees_of_saturation[GPoint];

        Variables.IntegrationCoefficient = integration_coefficients[GPoint];

        // Contributions to the left hand side
        if (CalculateStiffnessMatrixFlag) this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        // Contributions to the right hand side
        if (CalculateResidualVectorFlag)
            this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::InitializeElementVariables(ElementVariables& rVariables,
                                                                     const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Properties variables
    this->InitializeProperties(rVariables);

    // ProcessInfo variables
    rVariables.DtPressureCoefficient = rCurrentProcessInfo[DT_PRESSURE_COEFFICIENT];

    // Nodal Variables
    this->InitializeNodalPorePressureVariables(rVariables);
    this->InitializeNodalVolumeAccelerationVariables(rVariables);

    // Variables computed at each GP
    noalias(rVariables.Nu) = ZeroMatrix(TDim, TNumNodes * TDim);
    rVariables.Np.resize(TNumNodes, false);
    rVariables.GradNpT.resize(TNumNodes, TDim, false);

    // General Variables
    const GeometryType& Geom       = this->GetGeometry();
    const unsigned int  NumGPoints = Geom.IntegrationPointsNumber(this->GetIntegrationMethod());

    // shape functions
    (rVariables.NContainer).resize(NumGPoints, TNumNodes, false);
    rVariables.NContainer = Geom.ShapeFunctionsValues(this->GetIntegrationMethod());

    // gradient of shape functions and determinant of Jacobian
    (rVariables.detJContainer).resize(NumGPoints, false);

    Geom.ShapeFunctionsIntegrationPointsGradients(
        rVariables.DN_DXContainer, rVariables.detJContainer, this->GetIntegrationMethod());

    // Retention law
    rVariables.DegreeOfSaturation   = 1.0;
    rVariables.RelativePermeability = 1.0;
    rVariables.BishopCoefficient    = 1.0;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateAndAddLHS(MatrixType&       rLeftHandSideMatrix,
                                                             ElementVariables& rVariables)
{
    KRATOS_TRY

    this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix, rVariables);

    const auto permeability_matrix = GeoTransportEquationUtilities::CalculatePermeabilityMatrix<TDim, TNumNodes>(
        rVariables.GradNpT, rVariables.DynamicViscosityInverse, rVariables.PermeabilityMatrix,
        rVariables.RelativePermeability, rVariables.IntegrationCoefficient);
    rLeftHandSideMatrix += permeability_matrix;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix,
                                                                               ElementVariables& rVariables)
{
    KRATOS_TRY

    BoundedMatrix<double, TNumNodes, TNumNodes> compressibility_matrix =
        GeoTransportEquationUtilities::CalculateCompressibilityMatrix(
            rVariables.Np, rVariables.BiotModulusInverse, rVariables.IntegrationCoefficient);

    // Distribute compressibility block matrix into the elemental matrix
    rLeftHandSideMatrix += compressibility_matrix * rVariables.DtPressureCoefficient;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateAndAddRHS(VectorType&       rRightHandSideVector,
                                                             ElementVariables& rVariables,
                                                             unsigned int      GPoint)
{
    KRATOS_TRY

    this->CalculateAndAddCompressibilityFlow(rRightHandSideVector, rVariables);
    this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);
    this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector,
                                                                          ElementVariables& rVariables)
{
    KRATOS_TRY

    auto permeability_flow = this->CalculatePermeabilityFlow(rVariables);

    // Distribute permeability block vector into elemental vector
    rRightHandSideVector += permeability_flow;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector,
                                                                       ElementVariables& rVariables)
{
    KRATOS_TRY

    auto fluid_body_flow = this->CalculateFluidBodyFlow(rVariables);

    // Distribute fluid body flow block vector into elemental vector
    rRightHandSideVector += fluid_body_flow;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector,
                                                                             ElementVariables& rVariables)
{
    KRATOS_TRY

    auto compressibility_flow = this->CalculateCompressibilityFlow(rVariables);

    // Distribute compressibility block vector into elemental vector
    rRightHandSideVector += compressibility_flow;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateKinematics(ElementVariables& rVariables, unsigned int PointNumber)

{
    KRATOS_TRY

    // Setting the vector of shape functions and the matrix of the shape functions global gradients
    noalias(rVariables.Np)      = row(rVariables.NContainer, PointNumber);
    noalias(rVariables.GradNpT) = rVariables.DN_DXContainer[PointNumber];

    rVariables.detJ = rVariables.detJContainer[PointNumber];

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
unsigned int TransientPwElement<TDim, TNumNodes>::GetNumberOfDOF() const
{
    return TNumNodes;
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::DofsVectorType TransientPwElement<TDim, TNumNodes>::GetDofs() const
{
    return Geo::DofUtilities::ExtractDofsFromNodes(this->GetGeometry(), WATER_PRESSURE);
}

//----------------------------------------------------------------------------------------------------

template class TransientPwElement<2, 3>;
template class TransientPwElement<2, 4>;
template class TransientPwElement<3, 4>;
template class TransientPwElement<3, 8>;

template class TransientPwElement<2, 6>;
template class TransientPwElement<2, 8>;
template class TransientPwElement<2, 9>;
template class TransientPwElement<2, 10>;
template class TransientPwElement<2, 15>;
template class TransientPwElement<3, 10>;
template class TransientPwElement<3, 20>;
template class TransientPwElement<3, 27>;

} // Namespace Kratos
