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
#include "custom_utilities/check_utilities.h"
#include "custom_utilities/dof_utilities.h"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "includes/cfd_variables.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer TransientPwElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                             NodesArrayType const& ThisNodes,
                                                             PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new TransientPwElement(NewId, this->GetGeometry().Create(ThisNodes),
                                                   pProperties, this->GetStressStatePolicy().Clone()));
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer TransientPwElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                             GeometryType::Pointer pGeom,
                                                             PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(
        new TransientPwElement(NewId, pGeom, pProperties, this->GetStressStatePolicy().Clone()));
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix,
                                                              const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int n_DoF = this->GetNumberOfDOF();

    // Resizing mass matrix
    if (rMassMatrix.size1() != n_DoF) rMassMatrix.resize(n_DoF, n_DoF, false);
    noalias(rMassMatrix) = ZeroMatrix(n_DoF, n_DoF);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateDampingMatrix(MatrixType& rDampingMatrix,
                                                                 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int n_DoF = this->GetNumberOfDOF();

    // Compute Damping Matrix
    if (rDampingMatrix.size1() != n_DoF) rDampingMatrix.resize(n_DoF, n_DoF, false);
    noalias(rDampingMatrix) = ZeroMatrix(n_DoF, n_DoF);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::GetValuesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY

    const unsigned int n_DoF = this->GetNumberOfDOF();

    if (rValues.size() != n_DoF) rValues.resize(n_DoF, false);

    // Why are we constructing a zero vector here?
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rValues[i] = 0.0;
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::GetFirstDerivativesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY

    const unsigned int n_DoF = this->GetNumberOfDOF();

    if (rValues.size() != n_DoF) rValues.resize(n_DoF, false);

    // Why are we constructing a zero vector here?
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rValues[i] = 0.0;
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::GetSecondDerivativesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY

    const unsigned int n_DoF = this->GetNumberOfDOF();

    if (rValues.size() != n_DoF) rValues.resize(n_DoF, false);

    // Why are we constructing a zero vector here?
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rValues[i] = 0.0;
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const PropertiesType& r_properties = this->GetProperties();
    const GeometryType&   r_geom       = this->GetGeometry();
    const unsigned int    number_of_integration_points =
        r_geom.IntegrationPointsNumber(this->GetIntegrationMethod());

    if (mConstitutiveLawVector.size() != number_of_integration_points)
        mConstitutiveLawVector.resize(number_of_integration_points);
    for (auto& constitutive_law : mConstitutiveLawVector) {
        constitutive_law = nullptr;
    }

    if (mRetentionLawVector.size() != number_of_integration_points)
        mRetentionLawVector.resize(number_of_integration_points);
    for (auto& r_retention_law : mRetentionLawVector) {
        r_retention_law = RetentionLawFactory::Clone(r_properties);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
int TransientPwElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const PropertiesType& r_properties = this->GetProperties();
    const GeometryType&   r_geom       = this->GetGeometry();

    CheckUtilities::CheckDomainSize(r_geom.DomainSize(), this->Id());

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        if (r_geom[i].SolutionStepsDataHas(WATER_PRESSURE) == false)
            KRATOS_ERROR << "Missing variable WATER_PRESSURE on node " << r_geom[i].Id() << std::endl;

        if (r_geom[i].SolutionStepsDataHas(DT_WATER_PRESSURE) == false)
            KRATOS_ERROR << "Missing variable DT_WATER_PRESSURE on node " << r_geom[i].Id() << std::endl;

        if (r_geom[i].SolutionStepsDataHas(VOLUME_ACCELERATION) == false)
            KRATOS_ERROR << "Missing variable VOLUME_ACCELERATION on node " << r_geom[i].Id() << std::endl;

        if (r_geom[i].HasDofFor(WATER_PRESSURE) == false)
            KRATOS_ERROR << "Missing variable WATER_PRESSURE on node " << r_geom[i].Id() << std::endl;
    }

    // Verify ProcessInfo variables

    // Verify properties
    if (r_properties.Has(DENSITY_WATER) == false || r_properties[DENSITY_WATER] < 0.0)
        KRATOS_ERROR << "DENSITY_WATER does not exist in the material "
                        "properties or has an invalid value at element "
                     << this->Id() << std::endl;

    if (r_properties.Has(BULK_MODULUS_SOLID) == false || r_properties[BULK_MODULUS_SOLID] < 0.0)
        KRATOS_ERROR << "BULK_MODULUS_SOLID does not exist in the material "
                        "properties or has an invalid value at element "
                     << this->Id() << std::endl;

    if (r_properties.Has(POROSITY) == false || r_properties[POROSITY] < 0.0 || r_properties[POROSITY] > 1.0)
        KRATOS_ERROR << "POROSITY does not exist in the material properties or "
                        "has an invalid value at element "
                     << this->Id() << std::endl;

    if (TDim == 2) {
        // If this is a 2D problem, nodes must be in XY plane
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            if (r_geom[i].Z() != 0.0)
                KRATOS_ERROR << " Node with non-zero Z coordinate found. Id: " << r_geom[i].Id() << std::endl;
        }
    }

    // Verify specific properties
    if (r_properties.Has(BULK_MODULUS_FLUID) == false || r_properties[BULK_MODULUS_FLUID] < 0.0)
        KRATOS_ERROR << "BULK_MODULUS_FLUID does not exist in the material "
                        "properties or has an invalid value at element "
                     << this->Id() << std::endl;

    if (r_properties.Has(DYNAMIC_VISCOSITY) == false || r_properties[DYNAMIC_VISCOSITY] < 0.0)
        KRATOS_ERROR << "DYNAMIC_VISCOSITY does not exist in the material "
                        "properties or has an invalid value at element "
                     << this->Id() << std::endl;

    if (r_properties.Has(PERMEABILITY_XX) == false || r_properties[PERMEABILITY_XX] < 0.0)
        KRATOS_ERROR << "PERMEABILITY_XX does not exist in the material "
                        "properties or has an invalid value at element "
                     << this->Id() << std::endl;

    if (r_properties.Has(PERMEABILITY_YY) == false || r_properties[PERMEABILITY_YY] < 0.0)
        KRATOS_ERROR << "PERMEABILITY_YY does not exist in the material "
                        "properties or has an invalid value at element "
                     << this->Id() << std::endl;

    if (r_properties.Has(PERMEABILITY_XY) == false || r_properties[PERMEABILITY_XY] < 0.0)
        KRATOS_ERROR << "PERMEABILITY_XY does not exist in the material "
                        "properties or has an invalid value at element "
                     << this->Id() << std::endl;

    if (!r_properties.Has(BIOT_COEFFICIENT))
        KRATOS_ERROR << "BIOT_COEFFICIENT does not exist in the material "
                        "properties in element "
                     << this->Id() << std::endl;

    if constexpr (TDim > 2) {
        if (r_properties.Has(PERMEABILITY_ZZ) == false || r_properties[PERMEABILITY_ZZ] < 0.0)
            KRATOS_ERROR << "PERMEABILITY_ZZ does not exist in the material "
                            "properties or has an invalid value at element "
                         << this->Id() << std::endl;

        if (r_properties.Has(PERMEABILITY_YZ) == false || r_properties[PERMEABILITY_YZ] < 0.0)
            KRATOS_ERROR << "PERMEABILITY_YZ does not exist in the material "
                            "properties or has an invalid value at element "
                         << this->Id() << std::endl;

        if (r_properties.Has(PERMEABILITY_ZX) == false || r_properties[PERMEABILITY_ZX] < 0.0)
            KRATOS_ERROR << "PERMEABILITY_ZX does not exist in the material "
                            "properties or has an invalid value at element "
                         << this->Id() << std::endl;
    }

    if (!mRetentionLawVector.empty()) {
        return mRetentionLawVector[0]->Check(r_properties, rCurrentProcessInfo);
    }

    return 0;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    this->ResetHydraulicDischarge();

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::InitializeNonLinearIteration(const ProcessInfo&)
{
    // nothing
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::FinalizeNonLinearIteration(const ProcessInfo&)
{
    // nothing
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    this->CalculateHydraulicDischarge(rCurrentProcessInfo);

    KRATOS_CATCH("")
}

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

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                                                       std::vector<array_1d<double, 3>>& rOutput,
                                                                       const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& r_geom = this->GetGeometry();
    const IndexType     number_of_integration_points =
        r_geom.IntegrationPointsNumber(this->GetIntegrationMethod());
    if (rOutput.size() != number_of_integration_points)
        rOutput.resize(number_of_integration_points);

    if (rVariable == FLUID_FLUX_VECTOR) {
        std::vector<double> permeability_update_factors(number_of_integration_points, 1.0);
        const auto fluid_fluxes = this->CalculateFluidFluxes(permeability_update_factors, rCurrentProcessInfo);

        for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
             ++integration_point) {
            GeoElementUtilities::FillArray1dOutput(rOutput[integration_point], fluid_fluxes[integration_point]);
        }
    } else {
        for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
             ++integration_point) {
            noalias(rOutput[integration_point]) = ZeroVector(3);
        }
    }

    KRATOS_CATCH("")
}

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

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateAll(MatrixType&        rLeftHandSideMatrix,
                                                       VectorType&        rRightHandSideVector,
                                                       const ProcessInfo& rCurrentProcessInfo,
                                                       bool CalculateStiffnessMatrixFlag,
                                                       bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    // Previous definitions
    const PropertiesType&                           r_properties = this->GetProperties();
    const GeometryType&                             r_geom       = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        r_geom.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int number_of_integration_points = IntegrationPoints.size();

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
    std::vector<double> biot_coefficients(number_of_integration_points, r_properties[BIOT_COEFFICIENT]);
    const auto degrees_of_saturation     = this->CalculateDegreesOfSaturation(fluid_pressures);
    const auto derivatives_of_saturation = this->CalculateDerivativesOfSaturation(fluid_pressures);
    const auto biot_moduli_inverse = GeoTransportEquationUtilities::CalculateInverseBiotModuli(
        biot_coefficients, degrees_of_saturation, derivatives_of_saturation, r_properties);

    // Loop over integration points
    for (unsigned int integration_point = 0; integration_point < number_of_integration_points; ++integration_point) {
        // Compute GradNpT, B and StrainVector
        this->CalculateKinematics(Variables, integration_point);

        // Compute Nu and BodyAcceleration
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Variables.Nu, Variables.NContainer, integration_point);
        GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
            Variables.BodyAcceleration, Variables.NContainer, Variables.VolumeAcceleration, integration_point);

        Variables.RelativePermeability = relative_permeability_values[integration_point];
        Variables.BishopCoefficient    = bishop_coefficients[integration_point];

        Variables.BiotCoefficient    = biot_coefficients[integration_point];
        Variables.BiotModulusInverse = biot_moduli_inverse[integration_point];
        Variables.DegreeOfSaturation = degrees_of_saturation[integration_point];

        Variables.IntegrationCoefficient = integration_coefficients[integration_point];

        // Contributions to the left hand side
        if (CalculateStiffnessMatrixFlag) this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        // Contributions to the right hand side
        if (CalculateResidualVectorFlag)
            this->CalculateAndAddRHS(rRightHandSideVector, Variables, integration_point);
    }

    KRATOS_CATCH("")
}

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
    const GeometryType& r_geom = this->GetGeometry();
    const unsigned int  number_of_integration_points =
        r_geom.IntegrationPointsNumber(this->GetIntegrationMethod());

    // shape functions
    (rVariables.NContainer).resize(number_of_integration_points, TNumNodes, false);
    rVariables.NContainer = r_geom.ShapeFunctionsValues(this->GetIntegrationMethod());

    // gradient of shape functions and determinant of Jacobian
    (rVariables.detJContainer).resize(number_of_integration_points, false);

    r_geom.ShapeFunctionsIntegrationPointsGradients(
        rVariables.DN_DXContainer, rVariables.detJContainer, this->GetIntegrationMethod());

    // Retention law
    rVariables.DegreeOfSaturation   = 1.0;
    rVariables.RelativePermeability = 1.0;
    rVariables.BishopCoefficient    = 1.0;

    KRATOS_CATCH("")
}

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

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix,
                                                                               const ElementVariables& rVariables)
{
    KRATOS_TRY

    const auto compressibility_matrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix(
        rVariables.Np, rVariables.BiotModulusInverse, rVariables.IntegrationCoefficient);

    rLeftHandSideMatrix += compressibility_matrix * rVariables.DtPressureCoefficient;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateAndAddRHS(VectorType&       rRightHandSideVector,
                                                             ElementVariables& rVariables,
                                                             unsigned int)
{
    KRATOS_TRY

    this->CalculateAndAddCompressibilityFlow(rRightHandSideVector, rVariables);
    this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);
    this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector,
                                                                          const ElementVariables& rVariables)
{
    KRATOS_TRY

    rRightHandSideVector += this->CalculatePermeabilityFlow(rVariables);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector,
                                                                       const ElementVariables& rVariables)
{
    KRATOS_TRY

    noalias(rRightHandSideVector) += this->CalculateFluidBodyFlow(rVariables);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector,
                                                                             const ElementVariables& rVariables)
{
    KRATOS_TRY

    rRightHandSideVector += this->CalculateCompressibilityFlow(rVariables);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientPwElement<TDim, TNumNodes>::CalculateKinematics(ElementVariables& rVariables,
                                                              unsigned int IntegrationPointIndex)

{
    KRATOS_TRY

    // Setting the vector of shape functions and the matrix of the shape functions global gradients
    noalias(rVariables.Np)      = row(rVariables.NContainer, IntegrationPointIndex);
    noalias(rVariables.GradNpT) = rVariables.DN_DXContainer[IntegrationPointIndex];

    rVariables.detJ = rVariables.detJContainer[IntegrationPointIndex];

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
std::size_t TransientPwElement<TDim, TNumNodes>::GetNumberOfDOF() const
{
    return TNumNodes;
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::DofsVectorType TransientPwElement<TDim, TNumNodes>::GetDofs() const
{
    return Geo::DofUtilities::ExtractDofsFromNodes(this->GetGeometry(), WATER_PRESSURE);
}

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
