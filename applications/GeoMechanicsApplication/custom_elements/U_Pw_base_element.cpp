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
#include "custom_elements/U_Pw_base_element.h"
#include "custom_utilities/check_utilities.h"
#include "custom_utilities/dof_utilities.h"
#include "custom_utilities/equation_of_motion_utilities.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{

int UPwBaseElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Base class checks for positive area and Id > 0
    if (int ierr = Element::Check(rCurrentProcessInfo); ierr != 0) return ierr;

    const auto& r_geometry = this->GetGeometry();

    CheckUtilities::CheckHasNodalSolutionStepData(
        r_geometry, {std::cref(DISPLACEMENT), std::cref(VELOCITY), std::cref(ACCELERATION),
                     std::cref(WATER_PRESSURE), std::cref(DT_WATER_PRESSURE), std::cref(VOLUME_ACCELERATION)});
    CheckUtilities::CheckHasDofs(
        r_geometry, {std::cref(DISPLACEMENT_X), std::cref(DISPLACEMENT_Y), std::cref(WATER_PRESSURE)});
    if (this->GetGeometry().WorkingSpaceDimension() > 2)
        CheckUtilities::CheckHasDofs(r_geometry, {std::cref(DISPLACEMENT_Z)});

    const CheckProperties check_properties(this->GetProperties(), "material properties", this->Id(),
                                           CheckProperties::Bounds::AllInclusive);
    check_properties.Check(DENSITY_SOLID);
    check_properties.Check(DENSITY_WATER);
    check_properties.Check(BULK_MODULUS_SOLID);
    constexpr auto max_value_porosity = 1.0;
    check_properties.Check(POROSITY, max_value_porosity);

    if (this->GetGeometry().WorkingSpaceDimension() == 2)
        CheckUtilities::CheckForNonZeroZCoordinateIn2D(this->GetGeometry());

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
        using enum GeometryData::KratosGeometryOrderType;
        using enum GeometryData::IntegrationMethod;
    case Kratos_Cubic_Order:
        return GI_GAUSS_4;
    case Kratos_Quartic_Order:
        return GI_GAUSS_5;
    default:
        return GI_GAUSS_2;
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

void UPwBaseElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
    rSerializer.save("StressStatePolicy", mpStressStatePolicy);
    rSerializer.save("RetentionLawVector", mRetentionLawVector);
    rSerializer.save("StateVariablesFinalized", mStateVariablesFinalized);
    rSerializer.save("StressVector", mStressVector);
    rSerializer.save("ThisIntegrationMethod", static_cast<int>(mThisIntegrationMethod));
    rSerializer.save("IntegrationCoefficientsCalculator", mIntegrationCoefficientsCalculator);
}

void UPwBaseElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
    rSerializer.load("StressStatePolicy", mpStressStatePolicy);
    rSerializer.load("RetentionLawVector", mRetentionLawVector);
    rSerializer.load("StateVariablesFinalized", mStateVariablesFinalized);
    rSerializer.load("StressVector", mStressVector);
    int integration_method;
    rSerializer.load("ThisIntegrationMethod", integration_method);
    mThisIntegrationMethod = static_cast<IntegrationMethod>(integration_method);
    rSerializer.load("IntegrationCoefficientsCalculator", mIntegrationCoefficientsCalculator);
}

StressStatePolicy& UPwBaseElement::GetStressStatePolicy() const { return *mpStressStatePolicy; }

std::unique_ptr<IntegrationCoefficientModifier> UPwBaseElement::CloneIntegrationCoefficientModifier() const
{
    return mIntegrationCoefficientsCalculator.CloneModifier();
}

std::string UPwBaseElement::Info() const
{
    const std::string constitutive_info =
        !mConstitutiveLawVector.empty() ? mConstitutiveLawVector[0]->Info() : "not defined";

    std::ostringstream oss;
    oss << "U-Pw Base class Element #" << Id() << "\nConstitutive law: " << constitutive_info;

    return oss.str();
}

void UPwBaseElement::PrintInfo(std::ostream& rOStream) const { rOStream << Info(); }

} // Namespace Kratos
