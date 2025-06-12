// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

// Application includes
#include "custom_elements/U_Pw_base_element.hpp"
#include "custom_utilities/dof_utilities.h"
#include "custom_utilities/equation_of_motion_utilities.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{

int UPwBaseElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Base class checks for positive area and Id > 0
    if (int ierr = Element::Check(rCurrentProcessInfo); ierr != 0) return ierr;

    const PropertiesType& rProp = this->GetProperties();
    const GeometryType&   rGeom = this->GetGeometry();

    // verify nodal variables and dofs
    for (unsigned int i = 0; i < this->GetGeometry().PointsNumber(); ++i) {
        if (!rGeom[i].SolutionStepsDataHas(DISPLACEMENT))
            KRATOS_ERROR << "missing variable DISPLACEMENT on node " << rGeom[i].Id() << std::endl;

        if (!rGeom[i].SolutionStepsDataHas(VELOCITY))
            KRATOS_ERROR << "missing variable VELOCITY on node " << rGeom[i].Id() << std::endl;

        if (!rGeom[i].SolutionStepsDataHas(ACCELERATION))
            KRATOS_ERROR << "missing variable ACCELERATION on node " << rGeom[i].Id() << std::endl;

        if (!rGeom[i].SolutionStepsDataHas(WATER_PRESSURE))
            KRATOS_ERROR << "missing variable WATER_PRESSURE on node " << rGeom[i].Id() << std::endl;

        if (!rGeom[i].SolutionStepsDataHas(DT_WATER_PRESSURE))
            KRATOS_ERROR << "missing variable DT_WATER_PRESSURE on node " << rGeom[i].Id() << std::endl;

        if (!rGeom[i].SolutionStepsDataHas(VOLUME_ACCELERATION))
            KRATOS_ERROR << "missing variable VOLUME_ACCELERATION on node " << rGeom[i].Id() << std::endl;

        if (!rGeom[i].HasDofFor(DISPLACEMENT_X) || !rGeom[i].HasDofFor(DISPLACEMENT_Y) ||
            (this->GetGeometry().WorkingSpaceDimension() > 2 && !rGeom[i].HasDofFor(DISPLACEMENT_Z)))
            KRATOS_ERROR << "missing one of the dofs for the variable DISPLACEMENT on node "
                         << rGeom[i].Id() << std::endl;

        if (!rGeom[i].HasDofFor(WATER_PRESSURE))
            KRATOS_ERROR << "missing the dof for the variable WATER_PRESSURE on node "
                         << rGeom[i].Id() << std::endl;
    }

    // Verify properties
    if (!rProp.Has(DENSITY_SOLID) || rProp[DENSITY_SOLID] < 0.0)
        KRATOS_ERROR << "DENSITY_SOLID has Key zero, is not defined or has an "
                        "invalid value at element"
                     << this->Id() << std::endl;

    if (!rProp.Has(DENSITY_WATER) || rProp[DENSITY_WATER] < 0.0)
        KRATOS_ERROR << "DENSITY_WATER has Key zero, is not defined or has an "
                        "invalid value at element"
                     << this->Id() << std::endl;

    if (!rProp.Has(BULK_MODULUS_SOLID) || rProp[BULK_MODULUS_SOLID] < 0.0)
        KRATOS_ERROR
            << "BULK_MODULUS_SOLID has Key zero, is not defined or has an invalid value at element"
            << this->Id() << std::endl;

    if (!rProp.Has(POROSITY) || rProp[POROSITY] < 0.0 || rProp[POROSITY] > 1.0)
        KRATOS_ERROR << "POROSITY has Key zero, is not defined or has an invalid value at element"
                     << this->Id() << std::endl;

    if (this->GetGeometry().WorkingSpaceDimension() == 2) {
        // If this is a 2D problem, nodes must be in XY plane
        for (unsigned int i = 0; i < this->GetGeometry().PointsNumber(); ++i) {
            if (rGeom[i].Z() != 0.0)
                KRATOS_ERROR << " Node with non-zero Z coordinate found. Id: " << rGeom[i].Id() << std::endl;
        }
    }

    return 0;

    KRATOS_CATCH("")
}

void UPwBaseElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_properties = this->GetProperties();
    const auto& r_geometry   = this->GetGeometry();
    const auto number_of_integration_points = r_geometry.IntegrationPointsNumber(mThisIntegrationMethod);

    mConstitutiveLawVector.resize(number_of_integration_points);
    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i) {
        mConstitutiveLawVector[i] = r_properties[CONSTITUTIVE_LAW]->Clone();
        mConstitutiveLawVector[i]->InitializeMaterial(
            r_properties, r_geometry, row(r_geometry.ShapeFunctionsValues(mThisIntegrationMethod), i));
    }

    mRetentionLawVector.resize(number_of_integration_points);
    for (auto& r_retention_law : mRetentionLawVector) {
        r_retention_law = RetentionLawFactory::Clone(r_properties);
    }

    if (mStressVector.size() != number_of_integration_points) {
        mStressVector.resize(number_of_integration_points);
        for (auto& r_stress_vector : mStressVector) {
            r_stress_vector.resize(GetStressStatePolicy().GetVoigtSize());
            std::fill(r_stress_vector.begin(), r_stress_vector.end(), 0.0);
        }
    }
    std::vector<Vector> strain_vectors(number_of_integration_points,
                                       ZeroVector(GetStressStatePolicy().GetVoigtSize()));

    mStateVariablesFinalized.resize(number_of_integration_points);
    ConstitutiveLaw::Parameters cl_values;
    cl_values.SetProcessInfo(rCurrentProcessInfo);
    cl_values.SetMaterialProperties(r_properties);
    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i) {
        cl_values.SetStrainVector(strain_vectors[i]);
        cl_values.SetStressVector(mStressVector[i]);
        if (r_properties[CONSTITUTIVE_LAW]->Has(STATE_VARIABLES))
            mConstitutiveLawVector[i]->SetValue(STATE_VARIABLES, mStateVariablesFinalized[i], rCurrentProcessInfo);
        mConstitutiveLawVector[i]->InitializeMaterialResponseCauchy(cl_values);
    }

    KRATOS_CATCH("")
}

void UPwBaseElement::ResetConstitutiveLaw()
{
    KRATOS_TRY

    for (auto& r_stress_vector : mStressVector) {
        r_stress_vector.clear();
    }
    mStressVector.clear();

    for (auto& r_state_variables : mStateVariablesFinalized) {
        r_state_variables.clear();
    }
    mStateVariablesFinalized.clear();

    KRATOS_CATCH("")
}

void UPwBaseElement::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const
{
    rElementalDofList = GetDofs();
}

GeometryData::IntegrationMethod UPwBaseElement::GetIntegrationMethod() const
{
    switch (this->GetGeometry().GetGeometryOrderType()) {
    case GeometryData::Kratos_Cubic_Order:
        return GeometryData::IntegrationMethod::GI_GAUSS_4;
    case GeometryData::Kratos_Quartic_Order:
        return GeometryData::IntegrationMethod::GI_GAUSS_5;
    default:
        return GeometryData::IntegrationMethod::GI_GAUSS_2;
    }
}

void UPwBaseElement::CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                                          VectorType&        rRightHandSideVector,
                                          const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rLeftHandSideMatrix  = ZeroMatrix{this->GetNumberOfDOF(), this->GetNumberOfDOF()};
    rRightHandSideVector = ZeroVector{this->GetNumberOfDOF()};
    const auto CalculateStiffnessMatrixFlag = true;
    const auto CalculateResidualVectorFlag  = true;
    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                 CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);

    KRATOS_CATCH("")
}

void UPwBaseElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rLeftHandSideMatrix = ZeroMatrix{this->GetNumberOfDOF(), this->GetNumberOfDOF()};
    auto       dummy_right_hand_side_vector = Vector{};
    const auto CalculateStiffnessMatrixFlag = true;
    const auto CalculateResidualVectorFlag  = false;
    CalculateAll(rLeftHandSideMatrix, dummy_right_hand_side_vector, rCurrentProcessInfo,
                 CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);

    KRATOS_CATCH("")
}

void UPwBaseElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    auto dummy_left_hand_side_matrix        = Matrix{};
    rRightHandSideVector                    = ZeroVector{this->GetNumberOfDOF()};
    const auto CalculateStiffnessMatrixFlag = false;
    const auto CalculateResidualVectorFlag  = true;
    CalculateAll(dummy_left_hand_side_matrix, rRightHandSideVector, rCurrentProcessInfo,
                 CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);

    KRATOS_CATCH("")
}

void UPwBaseElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

void UPwBaseElement::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateMassMatrix method for a "
                    "particular element ... illegal operation!!"
                 << this->Id() << std::endl;

    KRATOS_CATCH("")
}

void UPwBaseElement::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    MatrixType mass_matrix = ZeroMatrix{this->GetNumberOfDOF(), this->GetNumberOfDOF()};
    this->CalculateMassMatrix(mass_matrix, rCurrentProcessInfo);

    MatrixType stiffness_matrix = ZeroMatrix{this->GetNumberOfDOF(), this->GetNumberOfDOF()};
    this->CalculateMaterialStiffnessMatrix(stiffness_matrix, rCurrentProcessInfo);

    const auto& r_prop = this->GetProperties();
    rDampingMatrix     = GeoEquationOfMotionUtilities::CalculateDampingMatrix(
        r_prop.Has(RAYLEIGH_ALPHA) ? r_prop[RAYLEIGH_ALPHA] : rCurrentProcessInfo[RAYLEIGH_ALPHA],
        r_prop.Has(RAYLEIGH_BETA) ? r_prop[RAYLEIGH_BETA] : rCurrentProcessInfo[RAYLEIGH_BETA],
        mass_matrix, stiffness_matrix);

    KRATOS_CATCH("")
}

void UPwBaseElement::GetValuesVector(Vector& rValues, int Step) const
{
    rValues = Geo::DofUtilities::ExtractSolutionStepValuesOfUPwDofs(GetDofs(), Step);
}

void UPwBaseElement::GetFirstDerivativesVector(Vector& rValues, int Step) const
{
    rValues = Geo::DofUtilities::ExtractFirstTimeDerivativesOfUPwDofs(GetDofs(), Step);
}

void UPwBaseElement::GetSecondDerivativesVector(Vector& rValues, int Step) const
{
    rValues = Geo::DofUtilities::ExtractSecondTimeDerivativesOfUPwDofs(GetDofs(), Step);
}

void UPwBaseElement::SetValuesOnIntegrationPoints(const Variable<Vector>&    rVariable,
                                                  const std::vector<Vector>& rValues,
                                                  const ProcessInfo&         rCurrentProcessInfo)
{
    KRATOS_TRY

    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
        mConstitutiveLawVector[i]->SetValue(rVariable, rValues[i], rCurrentProcessInfo);

    KRATOS_CATCH("")
}

void UPwBaseElement::SetValuesOnIntegrationPoints(const Variable<Matrix>&    rVariable,
                                                  const std::vector<Matrix>& rValues,
                                                  const ProcessInfo&         rCurrentProcessInfo)
{
    KRATOS_TRY

    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i) {
        mConstitutiveLawVector[i]->SetValue(rVariable, rValues[i], rCurrentProcessInfo);
    }

    KRATOS_CATCH("")
}

void UPwBaseElement::SetValuesOnIntegrationPoints(const Variable<double>&    rVariable,
                                                  const std::vector<double>& rValues,
                                                  const ProcessInfo&         rCurrentProcessInfo)
{
    KRATOS_TRY

    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i) {
        mConstitutiveLawVector[i]->SetValue(rVariable, rValues[i], rCurrentProcessInfo);
    }

    KRATOS_CATCH("")
}

void UPwBaseElement::CalculateOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
                                                  std::vector<ConstitutiveLaw::Pointer>& rValues,
                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == CONSTITUTIVE_LAW) {
        rValues.resize(mConstitutiveLawVector.size());
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i) {
            rValues[i] = mConstitutiveLawVector[i];
        }
    }

    KRATOS_CATCH("")
}

void UPwBaseElement::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                                  std::vector<array_1d<double, 3>>&    rValues,
                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateOnIntegrationPoints (array_1d<double, "
                    "3>) method for a particular element ... illegal operation!!"
                 << this->Id() << std::endl;

    KRATOS_CATCH("")
}

void UPwBaseElement::CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                  std::vector<Matrix>&    rValues,
                                                  const ProcessInfo&      rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateOnIntegrationPoints (Matrix) "
                    "method for a particular element ... illegal operation!!"
                 << this->Id() << std::endl;

    KRATOS_CATCH("")
}

void UPwBaseElement::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                  std::vector<Vector>&    rValues,
                                                  const ProcessInfo&      rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateOnIntegrationPoints (Vector) "
                    "method for a particular element ... illegal operation!!"
                 << this->Id() << std::endl;

    KRATOS_CATCH("")
}

void UPwBaseElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                                  std::vector<double>&    rValues,
                                                  const ProcessInfo&      rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateOnIntegrationPoints (double) "
                    "method for a particular element ... illegal operation!!"
                 << this->Id() << std::endl;

    KRATOS_CATCH("")
}

void UPwBaseElement::CalculateMaterialStiffnessMatrix(MatrixType& rStiffnessMatrix, const ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateMaterialStiffnessMatrix "
                    "method for a particular element ... illegal operation!!"
                 << this->Id() << std::endl;

    KRATOS_CATCH("")
}

void UPwBaseElement::CalculateAll(MatrixType&        rLeftHandSideMatrix,
                                  VectorType&        rRightHandSideVector,
                                  const ProcessInfo& CurrentProcessInfo,
                                  const bool         CalculateStiffnessMatrixFlag,
                                  const bool         CalculateResidualVectorFlag)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateAll method for a particular "
                    "element ... illegal operation!!"
                 << this->Id() << std::endl;

    KRATOS_CATCH("")
}

std::vector<double> UPwBaseElement::CalculateIntegrationCoefficients(
    const GeometryType::IntegrationPointsArrayType& rIntegrationPoints, const Vector& rDetJs) const
{
    return mIntegrationCoefficientsCalculator.Run<>(rIntegrationPoints, rDetJs, this);
}

void UPwBaseElement::CalculateDerivativesOnInitialConfiguration(
    double& rDetJ, Matrix& rJ0, Matrix& rInvJ0, Matrix& rDNu_DX0, unsigned int IntegrationPointIndex) const
{
    KRATOS_TRY

    const auto& r_geometry           = this->GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints(mThisIntegrationMethod);

    GeometryUtils::JacobianOnInitialConfiguration(r_geometry, r_integration_points[IntegrationPointIndex], rJ0);
    const auto& r_dn_de = r_geometry.ShapeFunctionsLocalGradients(mThisIntegrationMethod)[IntegrationPointIndex];
    MathUtils<>::InvertMatrix(rJ0, rInvJ0, rDetJ);
    GeometryUtils::ShapeFunctionsGradients(r_dn_de, rInvJ0, rDNu_DX0);

    KRATOS_CATCH("")
}

void UPwBaseElement::CalculateJacobianOnCurrentConfiguration(double& detJ, Matrix& rJ, Matrix& rInvJ, unsigned int GPoint) const
{
    KRATOS_TRY

    rJ = this->GetGeometry().Jacobian(rJ, GPoint, mThisIntegrationMethod);
    MathUtils<double>::InvertMatrix(rJ, rInvJ, detJ);

    KRATOS_CATCH("")
}

std::size_t UPwBaseElement::GetNumberOfDOF() const
{
    return this->GetGeometry().PointsNumber() * (this->GetGeometry().WorkingSpaceDimension() + 1);
}

Element::DofsVectorType UPwBaseElement::GetDofs() const
{
    return Geo::DofUtilities::ExtractUPwDofsFromNodes(this->GetGeometry(),
                                                      this->GetGeometry().WorkingSpaceDimension());
}

StressStatePolicy& UPwBaseElement::GetStressStatePolicy() const { return *mpStressStatePolicy; }

std::unique_ptr<IntegrationCoefficientModifier> UPwBaseElement::CloneIntegrationCoefficientModifier() const
{
    return mIntegrationCoefficientsCalculator.CloneModifier();
}

} // Namespace Kratos
