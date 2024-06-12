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

namespace Kratos
{

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
int UPwBaseElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Base class checks for positive area and Id > 0
    int ierr = Element::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;

    const PropertiesType& rProp = this->GetProperties();
    const GeometryType&   rGeom = this->GetGeometry();

    // verify nodal variables and dofs
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        if (rGeom[i].SolutionStepsDataHas(DISPLACEMENT) == false)
            KRATOS_ERROR << "missing variable DISPLACEMENT on node " << rGeom[i].Id() << std::endl;

        if (rGeom[i].SolutionStepsDataHas(VELOCITY) == false)
            KRATOS_ERROR << "missing variable VELOCITY on node " << rGeom[i].Id() << std::endl;

        if (rGeom[i].SolutionStepsDataHas(ACCELERATION) == false)
            KRATOS_ERROR << "missing variable ACCELERATION on node " << rGeom[i].Id() << std::endl;

        if (rGeom[i].SolutionStepsDataHas(WATER_PRESSURE) == false)
            KRATOS_ERROR << "missing variable WATER_PRESSURE on node " << rGeom[i].Id() << std::endl;

        if (rGeom[i].SolutionStepsDataHas(DT_WATER_PRESSURE) == false)
            KRATOS_ERROR << "missing variable DT_WATER_PRESSURE on node " << rGeom[i].Id() << std::endl;

        if (rGeom[i].SolutionStepsDataHas(VOLUME_ACCELERATION) == false)
            KRATOS_ERROR << "missing variable VOLUME_ACCELERATION on node " << rGeom[i].Id() << std::endl;

        if (rGeom[i].HasDofFor(DISPLACEMENT_X) == false ||
            rGeom[i].HasDofFor(DISPLACEMENT_Y) == false || rGeom[i].HasDofFor(DISPLACEMENT_Z) == false)
            KRATOS_ERROR << "missing one of the dofs for the variable "
                            "DISPLACEMENT on node "
                         << rGeom[i].Id() << std::endl;

        if (rGeom[i].HasDofFor(WATER_PRESSURE) == false)
            KRATOS_ERROR << "missing the dof for the variable WATER_PRESSURE "
                            "on node "
                         << rGeom[i].Id() << std::endl;
    }

    // Verify ProcessInfo variables

    // Verify properties
    if (rProp.Has(DENSITY_SOLID) == false || rProp[DENSITY_SOLID] < 0.0)
        KRATOS_ERROR << "DENSITY_SOLID has Key zero, is not defined or has an "
                        "invalid value at element"
                     << this->Id() << std::endl;

    if (rProp.Has(DENSITY_WATER) == false || rProp[DENSITY_WATER] < 0.0)
        KRATOS_ERROR << "DENSITY_WATER has Key zero, is not defined or has an "
                        "invalid value at element"
                     << this->Id() << std::endl;

    if (rProp.Has(YOUNG_MODULUS) == false) {
        if (rProp.Has(UDSM_NAME) == false)
            KRATOS_ERROR << "YOUNG_MODULUS has Key zero or is not defined at "
                            "element"
                         << this->Id() << std::endl;
    } else {
        if (rProp[YOUNG_MODULUS] <= 0.0)
            KRATOS_ERROR << "YOUNG_MODULUS has an invalid value at element" << this->Id() << std::endl;
    }

    if (rProp.Has(POISSON_RATIO) == false) {
        if (rProp.Has(UDSM_NAME) == false)
            KRATOS_ERROR << "POISSON_RATIO has Key zero or is not defined at "
                            "element"
                         << this->Id() << std::endl;
    } else {
        const double& PoissonRatio = rProp[POISSON_RATIO];
        if (PoissonRatio < 0.0 || PoissonRatio >= 0.5)
            KRATOS_ERROR << "POISSON_RATIO has an invalid value at element" << this->Id() << std::endl;
    }

    if (rProp.Has(BULK_MODULUS_SOLID) == false || rProp[BULK_MODULUS_SOLID] < 0.0)
        KRATOS_ERROR << "BULK_MODULUS_SOLID has Key zero, is not defined or "
                        "has an invalid value at element"
                     << this->Id() << std::endl;

    if (rProp.Has(POROSITY) == false || rProp[POROSITY] < 0.0 || rProp[POROSITY] > 1.0)
        KRATOS_ERROR << "POROSITY has Key zero, is not defined or has an "
                        "invalid value at element"
                     << this->Id() << std::endl;

    if (TDim == 2) {
        // If this is a 2D problem, nodes must be in XY plane
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            if (rGeom[i].Z() != 0.0)
                KRATOS_ERROR << " Node with non-zero Z coordinate found. Id: " << rGeom[i].Id() << std::endl;
        }
    }

    return 0;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
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
    for (unsigned int i = 0; i < mRetentionLawVector.size(); ++i) {
        mRetentionLawVector[i] = RetentionLawFactory::Clone(r_properties);
        mRetentionLawVector[i]->InitializeMaterial(
            r_properties, r_geometry, row(r_geometry.ShapeFunctionsValues(mThisIntegrationMethod), i));
    }

    if (mStressVector.size() != number_of_integration_points) {
        mStressVector.resize(number_of_integration_points);
        for (unsigned int i = 0; i < mStressVector.size(); ++i) {
            mStressVector[i].resize(GetStressStatePolicy().GetVoigtSize());
            std::fill(mStressVector[i].begin(), mStressVector[i].end(), 0.0);
        }
    }

    mStateVariablesFinalized.resize(number_of_integration_points);
    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i) {
        int nStateVariables = 0;
        nStateVariables = mConstitutiveLawVector[i]->GetValue(NUMBER_OF_UMAT_STATE_VARIABLES, nStateVariables);
        if (nStateVariables > 0) {
            mConstitutiveLawVector[i]->SetValue(STATE_VARIABLES, mStateVariablesFinalized[i], rCurrentProcessInfo);
        }
    }

    mIsInitialised = true;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::ResetConstitutiveLaw()
{
    KRATOS_TRY

    // erasing stress vectors
    for (unsigned int i = 0; i < mStressVector.size(); ++i) {
        mStressVector[i].clear();
    }
    mStressVector.clear();

    for (unsigned int i = 0; i < mStateVariablesFinalized.size(); ++i) {
        mStateVariablesFinalized[i].clear();
    }
    mStateVariablesFinalized.clear();

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const
{
    rElementalDofList = GetDofs();
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod UPwBaseElement<TDim, TNumNodes>::GetIntegrationMethod() const
{
    GeometryData::IntegrationMethod GI_GAUSS;

    switch (TNumNodes) {
    case 3:
        GI_GAUSS = GeometryData::IntegrationMethod::GI_GAUSS_2;
        break;
    case 6:
        GI_GAUSS = GeometryData::IntegrationMethod::GI_GAUSS_2;
        break;
    case 10:
        GI_GAUSS = GeometryData::IntegrationMethod::GI_GAUSS_4;
        break;
    case 15:
        GI_GAUSS = GeometryData::IntegrationMethod::GI_GAUSS_5;
        break;
    default:
        GI_GAUSS = GeometryData::IntegrationMethod::GI_GAUSS_2;
        break;
    }

    return GI_GAUSS;
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
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

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::CalculateLeftHandSide(MatrixType&        rLeftHandSideMatrix,
                                                            const ProcessInfo& rCurrentProcessInfo)
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

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                             const ProcessInfo& rCurrentProcessInfo)
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

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateMassMatrix method for a "
                    "particular element ... illegal operation!!"
                 << this->Id() << std::endl;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::CalculateDampingMatrix(MatrixType&        rDampingMatrix,
                                                             const ProcessInfo& rCurrentProcessInfo)
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

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::GetValuesVector(Vector& rValues, int Step) const
{
    rValues = Geo::DofUtilities::ExtractSolutionStepValuesOfUPwDofs(GetDofs(), Step);
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::GetFirstDerivativesVector(Vector& rValues, int Step) const
{
    rValues = Geo::DofUtilities::ExtractFirstTimeDerivativesOfUPwDofs(GetDofs(), Step);
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::GetSecondDerivativesVector(Vector& rValues, int Step) const
{
    rValues = Geo::DofUtilities::ExtractSecondTimeDerivativesOfUPwDofs(GetDofs(), Step);
}

//-------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::SetValuesOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                                   const std::vector<Vector>& rValues,
                                                                   const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
        mConstitutiveLawVector[i]->SetValue(rVariable, rValues[i], rCurrentProcessInfo);

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::SetValuesOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                                   const std::vector<Matrix>& rValues,
                                                                   const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i) {
        mConstitutiveLawVector[i]->SetValue(rVariable, rValues[i], rCurrentProcessInfo);
    }

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::SetValuesOnIntegrationPoints(const Variable<double>& rVariable,
                                                                   const std::vector<double>& rValues,
                                                                   const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i) {
        mConstitutiveLawVector[i]->SetValue(rVariable, rValues[i], rCurrentProcessInfo);
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
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

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                                                   std::vector<array_1d<double, 3>>& rValues,
                                                                   const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateOnIntegrationPoints (array_1d<double, "
                    "3>) method for a particular element ... illegal operation!!"
                 << this->Id() << std::endl;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                                   std::vector<Matrix>& rValues,
                                                                   const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateOnIntegrationPoints (Matrix) "
                    "method for a particular element ... illegal operation!!"
                 << this->Id() << std::endl;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                                   std::vector<Vector>& rValues,
                                                                   const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateOnIntegrationPoints (Vector) "
                    "method for a particular element ... illegal operation!!"
                 << this->Id() << std::endl;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                                                   std::vector<double>& rValues,
                                                                   const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateOnIntegrationPoints (double) "
                    "method for a particular element ... illegal operation!!"
                 << this->Id() << std::endl;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::CalculateMaterialStiffnessMatrix(MatrixType& rStiffnessMatrix,
                                                                       const ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateMaterialStiffnessMatrix "
                    "method for a particular element ... illegal operation!!"
                 << this->Id() << std::endl;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------

template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::CalculateAll(MatrixType&        rLeftHandSideMatrix,
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

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
double UPwBaseElement<TDim, TNumNodes>::CalculateIntegrationCoefficient(const GeometryType::IntegrationPointType& rIntegrationPoint,
                                                                        double detJ) const

{
    return mpStressStatePolicy->CalculateIntegrationCoefficient(rIntegrationPoint, detJ, GetGeometry());
}

template <unsigned int TDim, unsigned int TNumNodes>
std::vector<double> UPwBaseElement<TDim, TNumNodes>::CalculateIntegrationCoefficients(
    const GeometryType::IntegrationPointsArrayType& rIntegrationPoints, const Vector& rDetJs) const
{
    auto result = std::vector<double>{};
    std::transform(rIntegrationPoints.begin(), rIntegrationPoints.end(), rDetJs.begin(),
                   std::back_inserter(result), [this](const auto& rIntegrationPoint, const auto& rDetJ) {
        return this->CalculateIntegrationCoefficient(rIntegrationPoint, rDetJ);
    });
    return result;
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::CalculateDerivativesOnInitialConfiguration(
    double& detJ, Matrix& J0, Matrix& InvJ0, Matrix& DNu_DX0, unsigned int GPoint) const
{
    KRATOS_TRY

    const GeometryType&                             rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        rGeom.IntegrationPoints(mThisIntegrationMethod);

    GeometryUtils::JacobianOnInitialConfiguration(rGeom, IntegrationPoints[GPoint], J0);
    const Matrix& DN_De = rGeom.ShapeFunctionsLocalGradients(mThisIntegrationMethod)[GPoint];
    MathUtils<double>::InvertMatrix(J0, InvJ0, detJ);
    GeometryUtils::ShapeFunctionsGradients(DN_De, InvJ0, DNu_DX0);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwBaseElement<TDim, TNumNodes>::CalculateJacobianOnCurrentConfiguration(double& detJ,
                                                                              Matrix& rJ,
                                                                              Matrix& rInvJ,
                                                                              unsigned int GPoint) const
{
    KRATOS_TRY

    rJ = this->GetGeometry().Jacobian(rJ, GPoint, mThisIntegrationMethod);
    MathUtils<double>::InvertMatrix(rJ, rInvJ, detJ);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
unsigned int UPwBaseElement<TDim, TNumNodes>::GetNumberOfDOF() const
{
    return TNumNodes * (TDim + 1);
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::DofsVectorType UPwBaseElement<TDim, TNumNodes>::GetDofs() const
{
    return Geo::DofUtilities::ExtractUPwDofsFromNodes(this->GetGeometry(),
                                                      this->GetGeometry().WorkingSpaceDimension());
}

template <unsigned int TDim, unsigned int TNumNodes>
StressStatePolicy& UPwBaseElement<TDim, TNumNodes>::GetStressStatePolicy() const
{
    return *mpStressStatePolicy;
}

//----------------------------------------------------------------------------------------
template class UPwBaseElement<2, 3>;
template class UPwBaseElement<2, 4>;
template class UPwBaseElement<3, 4>;
template class UPwBaseElement<3, 6>;
template class UPwBaseElement<3, 8>;

template class UPwBaseElement<2, 6>;
template class UPwBaseElement<2, 8>;
template class UPwBaseElement<2, 9>;
template class UPwBaseElement<2, 10>;
template class UPwBaseElement<2, 15>;
template class UPwBaseElement<3, 10>;
template class UPwBaseElement<3, 20>;
template class UPwBaseElement<3, 27>;

} // Namespace Kratos
