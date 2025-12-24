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
#include "custom_elements/U_Pw_small_strain_interface_element.h"
#include "custom_utilities/check_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/interface_element_utilities.h"
#include "custom_utilities/stress_strain_utilities.h"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

void SetConstitutiveParameters(ConstitutiveLaw::Parameters& rConstitutiveParameters,
                               Matrix&                      rConstitutiveMatrix,
                               Vector&                      rStrainVector,
                               const Vector&                rNp,
                               const Matrix&                rGradNpT,
                               const double&                rDetF,
                               const Matrix&                rF)
{
    rConstitutiveParameters.SetConstitutiveMatrix(rConstitutiveMatrix);
    rConstitutiveParameters.SetStrainVector(rStrainVector);
    rConstitutiveParameters.SetShapeFunctionsValues(rNp);
    rConstitutiveParameters.SetShapeFunctionsDerivatives(rGradNpT);
    rConstitutiveParameters.SetDeterminantF(rDetF);
    rConstitutiveParameters.SetDeformationGradientF(rF);
    rConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
}

template <unsigned int TDim, unsigned int TNumNodes>
UPwSmallStrainInterfaceElement<TDim, TNumNodes>::UPwSmallStrainInterfaceElement(
    IndexType                                       NewId,
    const NodesArrayType&                           ThisNodes,
    std::unique_ptr<StressStatePolicy>              pStressStatePolicy,
    std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier)
    : UPwBaseElement(NewId, ThisNodes, std::move(pStressStatePolicy), std::move(pCoefficientModifier))
{
}

template <unsigned int TDim, unsigned int TNumNodes>
UPwSmallStrainInterfaceElement<TDim, TNumNodes>::UPwSmallStrainInterfaceElement(
    IndexType                                       NewId,
    GeometryType::Pointer                           pGeometry,
    std::unique_ptr<StressStatePolicy>              pStressStatePolicy,
    std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier)
    : UPwBaseElement(NewId, pGeometry, std::move(pStressStatePolicy), std::move(pCoefficientModifier))
{
}

template <unsigned int TDim, unsigned int TNumNodes>
UPwSmallStrainInterfaceElement<TDim, TNumNodes>::UPwSmallStrainInterfaceElement(
    IndexType                                       NewId,
    GeometryType::Pointer                           pGeometry,
    PropertiesType::Pointer                         pProperties,
    std::unique_ptr<StressStatePolicy>              pStressStatePolicy,
    std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier)
    : UPwBaseElement(NewId, pGeometry, pProperties, std::move(pStressStatePolicy), std::move(pCoefficientModifier))
{
    /// Lobatto integration method with the integration points located at the "mid plane nodes" of the interface
    mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwSmallStrainInterfaceElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                                         NodesArrayType const& ThisNodes,
                                                                         PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwSmallStrainInterfaceElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties,
        this->GetStressStatePolicy().Clone(), this->CloneIntegrationCoefficientModifier()));
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwSmallStrainInterfaceElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                                         GeometryType::Pointer pGeom,
                                                                         PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwSmallStrainInterfaceElement(
        NewId, pGeom, pProperties, this->GetStressStatePolicy().Clone(),
        this->CloneIntegrationCoefficientModifier()));
}

template <unsigned int TDim, unsigned int TNumNodes>
int UPwSmallStrainInterfaceElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const auto element_Id = this->Id();
    KRATOS_ERROR_IF(element_Id < 1)
        << "Element found with Id 0 or negative, element: " << element_Id << std::endl;

    // Verify generic variables
    int ierr = UPwBaseElement::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;

    const auto&           r_properties = this->GetProperties();
    const CheckProperties check_properties(r_properties, "property", element_Id,
                                           CheckProperties::Bounds::AllInclusive);
    check_properties.SingleUseBounds(CheckProperties::Bounds::AllExclusive).Check(MINIMUM_JOINT_WIDTH);
    check_properties.CheckAvailability(IGNORE_UNDRAINED);
    if (!r_properties[IGNORE_UNDRAINED]) {
        check_properties.Check(TRANSVERSAL_PERMEABILITY);
        check_properties.SingleUseBounds(CheckProperties::Bounds::AllExclusive).Check(BULK_MODULUS_FLUID);
        check_properties.SingleUseBounds(CheckProperties::Bounds::AllExclusive).Check(DYNAMIC_VISCOSITY);
    }

    check_properties.CheckAvailabilityAndSpecified(CONSTITUTIVE_LAW);
    ierr = r_properties[CONSTITUTIVE_LAW]->Check(r_properties, this->GetGeometry(), rCurrentProcessInfo);

    ConstitutiveLawUtilities::CheckHasStrainMeasure_Infinitesimal(r_properties, element_Id);

    constexpr auto expected_size = TDim == 2 ? VOIGT_SIZE_2D_INTERFACE : VOIGT_SIZE_3D_INTERFACE;
    ConstitutiveLawUtilities::CheckStrainSize(r_properties, expected_size, element_Id);

    return ierr;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    UPwBaseElement::Initialize(rCurrentProcessInfo);

    // Compute initial gap of the joint
    this->CalculateInitialGap(this->GetGeometry());

    // resize mStressVector
    const GeometryType&                             r_geometry = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        r_geometry.IntegrationPoints(mThisIntegrationMethod);
    const unsigned int NumGPoints = IntegrationPoints.size();
    if (const auto VoigtSize = TDim;
        (mStressVector.size() != NumGPoints) || (mStressVector[0].size() != VoigtSize)) {
        mStressVector.resize(NumGPoints);
        for (unsigned int i = 0; i < mStressVector.size(); ++i) {
            mStressVector[i].resize(VoigtSize);
            std::ranges::fill(mStressVector[i], 0.0);
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix,
                                                                          const ProcessInfo& rCurrentProcessInfo)
{
    // Resizing mass matrix, no mass in these interface elements.
    rMassMatrix = ZeroMatrix(this->GetNumberOfDOF(), this->GetNumberOfDOF());
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& r_geometry = this->GetGeometry();
    const Matrix&       NContainer = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
    array_1d<double, TNumNodes * TDim> DisplacementVector;
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(DisplacementVector, r_geometry, DISPLACEMENT);
    BoundedMatrix<double, TDim, TDim> RotationMatrix;
    this->CalculateRotationMatrix(RotationMatrix, r_geometry);
    BoundedMatrix<double, TDim, TNumNodes * TDim> Nu = ZeroMatrix(TDim, TNumNodes * TDim);
    array_1d<double, TDim>                        RelDispVector;

    // Create constitutive law parameters:
    Vector StrainVector(TDim);
    Vector StressVector(TDim);
    Matrix ConstitutiveMatrix(TDim, TDim);
    Vector Np(TNumNodes);
    Matrix GradNpT(TNumNodes, TDim);
    Matrix F    = identity_matrix<double>(TDim);
    double detF = 1.0;
    ConstitutiveLaw::Parameters ConstitutiveParameters(r_geometry, this->GetProperties(), rCurrentProcessInfo);
    SetConstitutiveParameters(ConstitutiveParameters, ConstitutiveMatrix, StrainVector, Np, GradNpT, detF, F);

    // Auxiliary output variables
    unsigned int NumGPoints = mConstitutiveLawVector.size();

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        InterfaceElementUtilities::CalculateNuMatrix(Nu, NContainer, GPoint);

        noalias(RelDispVector) = prod(Nu, DisplacementVector);

        noalias(StrainVector) = prod(RotationMatrix, RelDispVector);

        // Initialize constitutive law
        noalias(StressVector) = mStressVector[GPoint];
        ConstitutiveParameters.SetStressVector(StressVector);
        mConstitutiveLawVector[GPoint]->InitializeMaterialResponseCauchy(ConstitutiveParameters);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto&   r_properties = this->GetProperties();
    const auto&   r_geometry   = this->GetGeometry();
    const Matrix& NContainer   = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
    array_1d<double, TNumNodes * TDim> DisplacementVector;
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(DisplacementVector, r_geometry, DISPLACEMENT);
    BoundedMatrix<double, TDim, TDim> RotationMatrix;
    this->CalculateRotationMatrix(RotationMatrix, r_geometry);
    BoundedMatrix<double, TDim, TNumNodes * TDim> Nu = ZeroMatrix(TDim, TNumNodes * TDim);
    array_1d<double, TDim>                        RelDispVector;
    const double& MinimumJointWidth = r_properties[MINIMUM_JOINT_WIDTH];
    double        JointWidth;

    // Create constitutive law parameters:
    Vector StrainVector(TDim);
    Vector StressVector(TDim);
    Matrix ConstitutiveMatrix(TDim, TDim);
    Vector Np(TNumNodes);
    Matrix GradNpT(TNumNodes, TDim);
    Matrix F    = identity_matrix<double>(TDim);
    double detF = 1.0;
    ConstitutiveLaw::Parameters ConstitutiveParameters(r_geometry, r_properties, rCurrentProcessInfo);
    SetConstitutiveParameters(ConstitutiveParameters, ConstitutiveMatrix, StrainVector, Np, GradNpT, detF, F);

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
        InterfaceElementUtilities::CalculateNuMatrix(Nu, NContainer, GPoint);

        noalias(RelDispVector) = prod(Nu, DisplacementVector);

        noalias(StrainVector) = prod(RotationMatrix, RelDispVector);

        this->CheckAndCalculateJointWidth(JointWidth, ConstitutiveParameters,
                                          StrainVector[TDim - 1], MinimumJointWidth, GPoint);

        noalias(Np) = row(NContainer, GPoint);

        // compute constitutive tensor and/or stresses
        noalias(StressVector) = mStressVector[GPoint];
        ConstitutiveParameters.SetStressVector(StressVector);
        mConstitutiveLawVector[GPoint]->FinalizeMaterialResponseCauchy(ConstitutiveParameters);

        ModifyInactiveElementStress(JointWidth, mStressVector[GPoint]);

        mStateVariablesFinalized[GPoint] =
            mConstitutiveLawVector[GPoint]->GetValue(STATE_VARIABLES, mStateVariablesFinalized[GPoint]);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::ModifyInactiveElementStress(const double& JointWidth,
                                                                                  Vector& StressVector)
{
    KRATOS_TRY

    const auto&   r_properties      = this->GetProperties();
    const double& MinimumJointWidth = r_properties[MINIMUM_JOINT_WIDTH];

    if (JointWidth > MinimumJointWidth) {
        bool ConsiderGapClosure =
            r_properties.Has(CONSIDER_GAP_CLOSURE) ? r_properties[CONSIDER_GAP_CLOSURE] : false;
        if (ConsiderGapClosure) {
            const double decayFactor = 1.0;
            const double x           = (JointWidth / MinimumJointWidth) - 1.0;
            StressVector *= std::max(0.01, exp(-x * decayFactor));
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
Vector UPwSmallStrainInterfaceElement<TDim, TNumNodes>::SetFullStressVector(const Vector& rStressVector)
{
    Vector full_stress_vector(6, 0);

    if constexpr (TDim == 2) {
        full_stress_vector[INDEX_3D_ZZ] = rStressVector[INDEX_2D_INTERFACE_ZZ];
        full_stress_vector[INDEX_3D_XZ] = rStressVector[INDEX_2D_INTERFACE_XZ];
    } else if constexpr (TDim == 3) {
        full_stress_vector[INDEX_3D_ZZ] = rStressVector[INDEX_3D_INTERFACE_ZZ];
        full_stress_vector[INDEX_3D_YZ] = rStressVector[INDEX_3D_INTERFACE_YZ];
        full_stress_vector[INDEX_3D_XZ] = rStressVector[INDEX_3D_INTERFACE_XZ];
    }
    return full_stress_vector;
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& rGeom      = this->GetGeometry();
    const unsigned int  NumGPoints = rGeom.IntegrationPointsNumber(mThisIntegrationMethod);
    std::vector<double> GPValues(NumGPoints);

    // Printed on standard GiD Gauss points
    const unsigned int OutputGPoints = rGeom.IntegrationPointsNumber(this->GetIntegrationMethod());
    if (rValues.size() != OutputGPoints) rValues.resize(OutputGPoints);

    if (rVariable == VON_MISES_STRESS) {
        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            StressStrainUtilities EquivalentStress;
            Vector full_stress_vector = this->SetFullStressVector(mStressVector[GPoint]);
            GPValues[GPoint] = EquivalentStress.CalculateVonMisesStress(full_stress_vector);
        }

        this->InterpolateOutputDoubles(rValues, GPValues);
    } else if (rVariable == MEAN_EFFECTIVE_STRESS) {
        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            StressStrainUtilities EquivalentStress;
            Vector full_stress_vector = this->SetFullStressVector(mStressVector[GPoint]);
            GPValues[GPoint]          = EquivalentStress.CalculateMeanStress(full_stress_vector);
        }
        this->InterpolateOutputDoubles(rValues, GPValues);
    } else if (rVariable == MEAN_STRESS) {
        std::vector<Vector> StressVector;
        CalculateOnLobattoIntegrationPoints(TOTAL_STRESS_VECTOR, StressVector, rCurrentProcessInfo);

        // loop integration points
        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            StressStrainUtilities EquivalentStress;
            Vector full_stress_vector = this->SetFullStressVector(StressVector[GPoint]);
            GPValues[GPoint]          = EquivalentStress.CalculateMeanStress(full_stress_vector);
        }
        this->InterpolateOutputDoubles(rValues, GPValues);
    } else if (rVariable == ENGINEERING_VON_MISES_STRAIN || rVariable == ENGINEERING_VOLUMETRIC_STRAIN ||
               rVariable == GREEN_LAGRANGE_VON_MISES_STRAIN || rVariable == GREEN_LAGRANGE_VOLUMETRIC_STRAIN) {
        // Current variable is not calculated in U_Pw interface elements the output is set to 0
        for (unsigned int i = 0; i < OutputGPoints; ++i) {
            rValues[i] = 0;
        }
    } else if (rVariable == DAMAGE_VARIABLE) {
        // Variables computed on Lobatto points

        for (unsigned int i = 0; i < NumGPoints; ++i)
            GPValues[i] = mConstitutiveLawVector[i]->GetValue(rVariable, GPValues[i]);

        this->InterpolateOutputDoubles(rValues, GPValues);
    } else if (rVariable == STATE_VARIABLE) {
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
            rValues[i] = mConstitutiveLawVector[i]->GetValue(rVariable, rValues[i]);
    } else if (rVariable == JOINT_WIDTH) {
        // Variables computed on Lobatto points
        std::vector<array_1d<double, 3>> GPAuxValues(NumGPoints);
        this->CalculateOnIntegrationPoints(LOCAL_RELATIVE_DISPLACEMENT_VECTOR, GPAuxValues, rCurrentProcessInfo);

        for (unsigned int i = 0; i < NumGPoints; ++i) {
            GPValues[i] = mInitialGap[i] + GPAuxValues[i][TDim - 1];
        }

        // Printed on standard GiD Gauss points
        this->InterpolateOutputDoubles(rValues, GPValues);
    } else if (rVariable == DEGREE_OF_SATURATION || rVariable == EFFECTIVE_SATURATION ||
               rVariable == BISHOP_COEFFICIENT || rVariable == DERIVATIVE_OF_SATURATION ||
               rVariable == RELATIVE_PERMEABILITY) {
        // Defining necessary variables

        // Element variables
        InterfaceElementVariables Variables;
        this->InitializeElementVariables(Variables, rGeom, this->GetProperties(), rCurrentProcessInfo);

        RetentionLaw::Parameters RetentionParameters(this->GetProperties());

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            Variables.FluidPressure = GeoTransportEquationUtilities::CalculateFluidPressure(
                Variables.Np, Variables.PressureVector);
            RetentionParameters.SetFluidPressure(Variables.FluidPressure);

            if (rVariable == DEGREE_OF_SATURATION)
                GPValues[GPoint] = mRetentionLawVector[GPoint]->CalculateSaturation(RetentionParameters);
            if (rVariable == EFFECTIVE_SATURATION)
                GPValues[GPoint] =
                    mRetentionLawVector[GPoint]->CalculateEffectiveSaturation(RetentionParameters);
            if (rVariable == BISHOP_COEFFICIENT)
                GPValues[GPoint] = mRetentionLawVector[GPoint]->CalculateBishopCoefficient(RetentionParameters);
            if (rVariable == DERIVATIVE_OF_SATURATION)
                GPValues[GPoint] =
                    mRetentionLawVector[GPoint]->CalculateDerivativeOfSaturation(RetentionParameters);
            if (rVariable == RELATIVE_PERMEABILITY)
                GPValues[GPoint] =
                    mRetentionLawVector[GPoint]->CalculateRelativePermeability(RetentionParameters);
        }

        this->InterpolateOutputDoubles(rValues, GPValues);
    } else if (rVariable == CONFINED_STIFFNESS || rVariable == SHEAR_STIFFNESS) {
        // set the correct index of the variable in the constitutive matrix
        size_t variable_index;
        if (rVariable == CONFINED_STIFFNESS) {
            if (TDim == 2) {
                variable_index = INDEX_2D_INTERFACE_ZZ;
            } else if (TDim == 3) {
                variable_index = INDEX_3D_INTERFACE_ZZ;
            } else {
                KRATOS_ERROR << "CONFINED_STIFFNESS can not be retrieved for dim " << TDim
                             << " in element: " << this->Id() << std::endl;
            }
        } else if (rVariable == SHEAR_STIFFNESS) {
            if (TDim == 2) {
                variable_index = INDEX_2D_INTERFACE_XZ;
            } else if (TDim == 3) {
                variable_index = INDEX_3D_INTERFACE_XZ;
            } else {
                KRATOS_ERROR << "SHEAR_STIFFNESS can not be retrieved for dim " << TDim
                             << " in element: " << this->Id() << std::endl;
            }
        }

        InterfaceElementVariables Variables;
        const PropertiesType&     rProp = this->GetProperties();

        this->InitializeElementVariables(Variables, rGeom, rProp, rCurrentProcessInfo);

        // Containers of variables at all integration points
        const Matrix& NContainer = rGeom.ShapeFunctionsValues(mThisIntegrationMethod);
        const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer =
            rGeom.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
        GeometryType::JacobiansType JContainer(NumGPoints);
        rGeom.Jacobian(JContainer, mThisIntegrationMethod);

        ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, rProp, rCurrentProcessInfo);
        ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        ConstitutiveLawUtilities::SetConstitutiveParameters(
            ConstitutiveParameters, Variables.StrainVector, Variables.ConstitutiveMatrix,
            Variables.Np, Variables.GradNpT, Variables.F, Variables.detF);

        // Auxiliary variables
        const double&          MinimumJointWidth = rProp[MINIMUM_JOINT_WIDTH];
        array_1d<double, TDim> RelDispVector;
        SFGradAuxVariables     SFGradAuxVars;

        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            // Compute Np, StrainVector, JointWidth, GradNpT
            noalias(Variables.Np) = row(NContainer, GPoint);
            InterfaceElementUtilities::CalculateNuMatrix(Variables.Nu, NContainer, GPoint);
            noalias(RelDispVector)          = prod(Variables.Nu, Variables.DisplacementVector);
            noalias(Variables.StrainVector) = prod(Variables.RotationMatrix, RelDispVector);

            this->CheckAndCalculateJointWidth(Variables.JointWidth, ConstitutiveParameters,
                                              Variables.StrainVector[TDim - 1], MinimumJointWidth, GPoint);

            this->CalculateShapeFunctionsGradients<Matrix>(
                Variables.GradNpT, SFGradAuxVars, JContainer[GPoint], Variables.RotationMatrix,
                DN_DeContainer[GPoint], NContainer, Variables.JointWidth, GPoint);

            // Compute constitutive tensor
            noalias(Variables.StressVector) = mStressVector[GPoint];
            ConstitutiveParameters.SetStressVector(Variables.StressVector);

            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            // get variable from constitutive matrix
            GPValues[GPoint] = Variables.ConstitutiveMatrix(variable_index, variable_index);
        }

        // Printed on standard GiD Gauss points
        this->InterpolateOutputDoubles(rValues, GPValues);
    } else {
        // Variables computed on Lobatto points

        for (unsigned int i = 0; i < NumGPoints; ++i)
            GPValues[i] = mConstitutiveLawVector[i]->GetValue(rVariable, GPValues[i]);

        // Printed on standard GiD Gauss points
        this->InterpolateOutputDoubles(rValues, GPValues);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable, std::vector<array_1d<double, 3>>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == FLUID_FLUX_VECTOR || rVariable == LOCAL_STRESS_VECTOR ||
        rVariable == LOCAL_RELATIVE_DISPLACEMENT_VECTOR || rVariable == LOCAL_FLUID_FLUX_VECTOR) {
        // Variables computed on Lobatto points
        const GeometryType& r_geometry = this->GetGeometry();

        const unsigned int nLobottoGPoints = r_geometry.IntegrationPointsNumber(mThisIntegrationMethod);
        std::vector<array_1d<double, 3>> GPValues(nLobottoGPoints);
        this->CalculateOnLobattoIntegrationPoints(rVariable, GPValues, rCurrentProcessInfo);

        // Printed on standard GiD Gauss points
        const unsigned int nOutputGPoints = r_geometry.IntegrationPointsNumber(this->GetIntegrationMethod());
        if (rValues.size() != nOutputGPoints) rValues.resize(nOutputGPoints);

        this->InterpolateOutputValues<array_1d<double, 3>>(rValues, GPValues);
    } else {
        // Variables computed on Lobatto points
        const GeometryType& r_geometry = this->GetGeometry();
        const unsigned int  NumGPoints = r_geometry.IntegrationPointsNumber(mThisIntegrationMethod);
        std::vector<array_1d<double, 3>> GPValues(NumGPoints);

        for (unsigned int i = 0; i < NumGPoints; ++i)
            GPValues[i] = mConstitutiveLawVector[i]->GetValue(rVariable, GPValues[i]);

        // Printed on standard GiD Gauss points
        const unsigned int OutputGPoints = r_geometry.IntegrationPointsNumber(this->GetIntegrationMethod());
        if (rValues.size() != OutputGPoints) rValues.resize(OutputGPoints);

        this->InterpolateOutputValues<array_1d<double, 3>>(rValues, GPValues);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == PERMEABILITY_MATRIX || rVariable == LOCAL_PERMEABILITY_MATRIX) {
        // Variables computed on Lobatto points
        const GeometryType& r_geometry = this->GetGeometry();
        std::vector<Matrix> GPValues(r_geometry.IntegrationPointsNumber(mThisIntegrationMethod));

        this->CalculateOnLobattoIntegrationPoints(rVariable, GPValues, rCurrentProcessInfo);

        // Printed on standard GiD Gauss points
        const unsigned int OutputGPoints = r_geometry.IntegrationPointsNumber(this->GetIntegrationMethod());
        if (rValues.size() != OutputGPoints) rValues.resize(OutputGPoints);

        for (unsigned int GPoint = 0; GPoint < OutputGPoints; ++GPoint)
            rValues[GPoint].resize(TDim, TDim, false);

        this->InterpolateOutputValues<Matrix>(rValues, GPValues);
    } else {
        // Printed on standard GiD Gauss points
        const unsigned int OutputGPoints =
            this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
        if (rValues.size() != OutputGPoints) rValues.resize(OutputGPoints);

        for (unsigned int i = 0; i < OutputGPoints; ++i) {
            rValues[i].resize(TDim, TDim, false);
            noalias(rValues[i]) = ZeroMatrix(TDim, TDim);
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateOnLobattoIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable, std::vector<array_1d<double, 3>>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == FLUID_FLUX_VECTOR) {
        const auto&        r_properties = this->GetProperties();
        const auto&        r_geometry   = this->GetGeometry();
        const unsigned int NumGPoints = r_geometry.IntegrationPointsNumber(mThisIntegrationMethod);

        // Defining the shape functions, the Jacobian and the shape functions local gradients Containers
        const Matrix& NContainer = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
        const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer =
            r_geometry.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
        GeometryType::JacobiansType JContainer(NumGPoints);
        r_geometry.Jacobian(JContainer, mThisIntegrationMethod);

        // Defining necessary variables
        BoundedMatrix<double, TDim, TDim> RotationMatrix;
        this->CalculateRotationMatrix(RotationMatrix, r_geometry);

        array_1d<double, TDim> LocalRelDispVector;
        array_1d<double, TDim> RelDispVector;
        const double&          MinimumJointWidth = r_properties[MINIMUM_JOINT_WIDTH];

        BoundedMatrix<double, TNumNodes, TDim> GradNpT;

        const double& TransversalPermeability = r_properties[TRANSVERSAL_PERMEABILITY];

        array_1d<double, TDim> LocalFluidFlux;
        array_1d<double, TDim> GradPressureTerm;
        array_1d<double, TDim> FluidFlux;
        SFGradAuxVariables     SFGradAuxVars;

        // Element variables
        InterfaceElementVariables Variables;
        this->InitializeElementVariables(Variables, r_geometry, r_properties, rCurrentProcessInfo);

        RetentionLaw::Parameters RetentionParameters(r_properties);

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            InterfaceElementUtilities::CalculateNuMatrix(Variables.Nu, NContainer, GPoint);

            noalias(RelDispVector) = prod(Variables.Nu, Variables.DisplacementVector);

            noalias(LocalRelDispVector) = prod(RotationMatrix, RelDispVector);

            double JointWidth;
            this->CalculateJointWidth(JointWidth, LocalRelDispVector[TDim - 1], MinimumJointWidth, GPoint);

            this->CalculateShapeFunctionsGradients<BoundedMatrix<double, TNumNodes, TDim>>(
                GradNpT, SFGradAuxVars, JContainer[GPoint], RotationMatrix, DN_DeContainer[GPoint],
                NContainer, JointWidth, GPoint);

            GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
                Variables.BodyAcceleration, NContainer, Variables.VolumeAcceleration, GPoint);

            InterfaceElementUtilities::FillPermeabilityMatrix(Variables.LocalPermeabilityMatrix,
                                                              JointWidth, TransversalPermeability);

            noalias(GradPressureTerm) = prod(trans(GradNpT), Variables.PressureVector);
            noalias(GradPressureTerm) +=
                PORE_PRESSURE_SIGN_FACTOR * r_properties[DENSITY_WATER] * Variables.BodyAcceleration;

            noalias(LocalFluidFlux) = PORE_PRESSURE_SIGN_FACTOR * Variables.DynamicViscosityInverse *
                                      Variables.RelativePermeability *
                                      prod(Variables.LocalPermeabilityMatrix, GradPressureTerm);

            noalias(FluidFlux) = prod(trans(RotationMatrix), LocalFluidFlux);

            GeoElementUtilities::FillArray1dOutput(rOutput[GPoint], FluidFlux);
        }
    } else if (rVariable == LOCAL_STRESS_VECTOR) {
        array_1d<double, TDim> LocalStressVector;

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            noalias(LocalStressVector) = mStressVector[GPoint];
            GeoElementUtilities::FillArray1dOutput(rOutput[GPoint], LocalStressVector);
        }
    } else if (rVariable == LOCAL_RELATIVE_DISPLACEMENT_VECTOR) {
        // Defining necessary variables
        const GeometryType& r_geometry = this->GetGeometry();
        const Matrix&       NContainer = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
        array_1d<double, TNumNodes * TDim> DisplacementVector;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(DisplacementVector, r_geometry, DISPLACEMENT);
        BoundedMatrix<double, TDim, TDim> RotationMatrix;
        this->CalculateRotationMatrix(RotationMatrix, r_geometry);
        BoundedMatrix<double, TDim, TNumNodes * TDim> Nu = ZeroMatrix(TDim, TNumNodes * TDim);
        array_1d<double, TDim>                        LocalRelDispVector;
        array_1d<double, TDim>                        RelDispVector;

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            InterfaceElementUtilities::CalculateNuMatrix(Nu, NContainer, GPoint);

            noalias(RelDispVector) = prod(Nu, DisplacementVector);

            noalias(LocalRelDispVector) = prod(RotationMatrix, RelDispVector);

            GeoElementUtilities::FillArray1dOutput(rOutput[GPoint], LocalRelDispVector);
        }
    } else if (rVariable == LOCAL_FLUID_FLUX_VECTOR) {
        const auto&        r_properties = this->GetProperties();
        const auto&        r_geometry   = this->GetGeometry();
        const unsigned int NumGPoints = r_geometry.IntegrationPointsNumber(mThisIntegrationMethod);

        // Defining the shape functions, the Jacobian and the shape functions local gradients Containers
        const Matrix& NContainer = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
        const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer =
            r_geometry.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
        GeometryType::JacobiansType JContainer(NumGPoints);
        r_geometry.Jacobian(JContainer, mThisIntegrationMethod);

        // Defining necessary variables
        array_1d<double, TNumNodes> PressureVector;
        for (unsigned int i = 0; i < TNumNodes; ++i)
            PressureVector[i] = r_geometry[i].FastGetSolutionStepValue(WATER_PRESSURE);
        array_1d<double, TNumNodes * TDim> VolumeAcceleration;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(VolumeAcceleration, r_geometry,
                                                                     VOLUME_ACCELERATION);
        array_1d<double, TDim>             BodyAcceleration;
        array_1d<double, TNumNodes * TDim> DisplacementVector;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(DisplacementVector, r_geometry, DISPLACEMENT);
        BoundedMatrix<double, TDim, TDim> RotationMatrix;
        this->CalculateRotationMatrix(RotationMatrix, r_geometry);
        BoundedMatrix<double, TDim, TNumNodes * TDim> Nu = ZeroMatrix(TDim, TNumNodes * TDim);
        array_1d<double, TDim>                        LocalRelDispVector;
        array_1d<double, TDim>                        RelDispVector;
        const double& MinimumJointWidth = r_properties[MINIMUM_JOINT_WIDTH];
        double        JointWidth;
        BoundedMatrix<double, TNumNodes, TDim> GradNpT;
        const double& TransversalPermeability = r_properties[TRANSVERSAL_PERMEABILITY];
        BoundedMatrix<double, TDim, TDim> LocalPermeabilityMatrix = ZeroMatrix(TDim, TDim);
        const double           DynamicViscosityInverse = 1.0 / r_properties[DYNAMIC_VISCOSITY];
        const double&          FluidDensity            = r_properties[DENSITY_WATER];
        array_1d<double, TDim> LocalFluidFlux;
        array_1d<double, TDim> GradPressureTerm;
        SFGradAuxVariables     SFGradAuxVars;

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            InterfaceElementUtilities::CalculateNuMatrix(Nu, NContainer, GPoint);

            noalias(RelDispVector) = prod(Nu, DisplacementVector);

            noalias(LocalRelDispVector) = prod(RotationMatrix, RelDispVector);

            this->CalculateJointWidth(JointWidth, LocalRelDispVector[TDim - 1], MinimumJointWidth, GPoint);

            this->CalculateShapeFunctionsGradients<BoundedMatrix<double, TNumNodes, TDim>>(
                GradNpT, SFGradAuxVars, JContainer[GPoint], RotationMatrix, DN_DeContainer[GPoint],
                NContainer, JointWidth, GPoint);

            GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
                BodyAcceleration, NContainer, VolumeAcceleration, GPoint);

            InterfaceElementUtilities::FillPermeabilityMatrix(LocalPermeabilityMatrix, JointWidth,
                                                              TransversalPermeability);

            noalias(GradPressureTerm) = prod(trans(GradNpT), PressureVector);
            noalias(GradPressureTerm) += PORE_PRESSURE_SIGN_FACTOR * FluidDensity * BodyAcceleration;

            noalias(LocalFluidFlux) = -DynamicViscosityInverse * prod(LocalPermeabilityMatrix, GradPressureTerm);

            GeoElementUtilities::FillArray1dOutput(rOutput[GPoint], LocalFluidFlux);
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateOnLobattoIntegrationPoints(
    const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Defining necessary variables
    const GeometryType& rGeom      = this->GetGeometry();
    const IndexType     NumGPoints = rGeom.IntegrationPointsNumber(mThisIntegrationMethod);

    // calculated on Lobatto points
    rOutput.resize(NumGPoints);

    if (rVariable == TOTAL_STRESS_VECTOR) {
        // Defining necessary variables
        const PropertiesType& rProp = this->GetProperties();

        // Containers of variables at all integration points
        const Matrix& NContainer = rGeom.ShapeFunctionsValues(mThisIntegrationMethod);
        const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer =
            rGeom.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
        GeometryType::JacobiansType JContainer(NumGPoints);
        rGeom.Jacobian(JContainer, mThisIntegrationMethod);

        ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, rProp, rCurrentProcessInfo);
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

        // Element variables
        InterfaceElementVariables Variables;
        this->InitializeElementVariables(Variables, rGeom, rProp, rCurrentProcessInfo);

        RetentionLaw::Parameters RetentionParameters(this->GetProperties());

        Vector VoigtVector(mStressVector[0].size());
        noalias(VoigtVector) = ZeroVector(VoigtVector.size());

        for (unsigned int i = 0; i < mStressVector[0].size(); ++i)
            VoigtVector[i] = 1.0;

        const bool hasBiotCoefficient = rProp.Has(BIOT_COEFFICIENT);

        Vector TotalStressVector(mStressVector[0].size());

        ConstitutiveLawUtilities::SetConstitutiveParameters(
            ConstitutiveParameters, Variables.StrainVector, Variables.ConstitutiveMatrix,
            Variables.Np, Variables.GradNpT, Variables.F, Variables.detF);

        // Auxiliary variables
        const double&          MinimumJointWidth = rProp[MINIMUM_JOINT_WIDTH];
        array_1d<double, TDim> RelDispVector;
        SFGradAuxVariables     SFGradAuxVars;

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            // Compute Np, StrainVector, JointWidth, GradNpT
            noalias(Variables.Np) = row(NContainer, GPoint);
            InterfaceElementUtilities::CalculateNuMatrix(Variables.Nu, NContainer, GPoint);
            noalias(RelDispVector)          = prod(Variables.Nu, Variables.DisplacementVector);
            noalias(Variables.StrainVector) = prod(Variables.RotationMatrix, RelDispVector);

            this->CheckAndCalculateJointWidth(Variables.JointWidth, ConstitutiveParameters,
                                              Variables.StrainVector[TDim - 1], MinimumJointWidth, GPoint);

            this->CalculateShapeFunctionsGradients<Matrix>(
                Variables.GradNpT, SFGradAuxVars, JContainer[GPoint], Variables.RotationMatrix,
                DN_DeContainer[GPoint], NContainer, Variables.JointWidth, GPoint);

            // compute constitutive tensor and/or stresses
            noalias(Variables.StressVector) = mStressVector[GPoint];
            ConstitutiveParameters.SetStressVector(Variables.StressVector);
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            InitializeBiotCoefficients(Variables, hasBiotCoefficient);

            this->CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

            noalias(TotalStressVector) = mStressVector[GPoint];
            noalias(TotalStressVector) += PORE_PRESSURE_SIGN_FACTOR * Variables.BiotCoefficient *
                                          Variables.BishopCoefficient * Variables.FluidPressure * VoigtVector;

            rOutput[GPoint] = TotalStressVector;
        }
    }
    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateOnLobattoIntegrationPoints(
    const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == PERMEABILITY_MATRIX) {
        const auto& r_geometry   = this->GetGeometry();
        const auto& r_properties = this->GetProperties();

        // Defining the shape functions container
        const Matrix& NContainer = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);

        // Defining necessary variables
        array_1d<double, TNumNodes * TDim> DisplacementVector;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(DisplacementVector, r_geometry, DISPLACEMENT);
        BoundedMatrix<double, TDim, TDim> RotationMatrix;
        this->CalculateRotationMatrix(RotationMatrix, r_geometry);
        BoundedMatrix<double, TDim, TNumNodes * TDim> Nu = ZeroMatrix(TDim, TNumNodes * TDim);
        array_1d<double, TDim>                        LocalRelDispVector;
        array_1d<double, TDim>                        RelDispVector;
        const double& MinimumJointWidth = r_properties[MINIMUM_JOINT_WIDTH];
        double        JointWidth;
        const double& TransversalPermeability = r_properties[TRANSVERSAL_PERMEABILITY];
        BoundedMatrix<double, TDim, TDim> LocalPermeabilityMatrix = ZeroMatrix(TDim, TDim);
        BoundedMatrix<double, TDim, TDim> PermeabilityMatrix;

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            InterfaceElementUtilities::CalculateNuMatrix(Nu, NContainer, GPoint);

            noalias(RelDispVector) = prod(Nu, DisplacementVector);

            noalias(LocalRelDispVector) = prod(RotationMatrix, RelDispVector);

            this->CalculateJointWidth(JointWidth, LocalRelDispVector[TDim - 1], MinimumJointWidth, GPoint);

            InterfaceElementUtilities::FillPermeabilityMatrix(LocalPermeabilityMatrix, JointWidth,
                                                              TransversalPermeability);

            noalias(PermeabilityMatrix) =
                prod(trans(RotationMatrix),
                     BoundedMatrix<double, TDim, TDim>(prod(LocalPermeabilityMatrix, RotationMatrix)));

            rOutput[GPoint].resize(TDim, TDim, false);
            noalias(rOutput[GPoint]) = PermeabilityMatrix;
        }
    } else if (rVariable == LOCAL_PERMEABILITY_MATRIX) {
        const auto& r_geometry   = this->GetGeometry();
        const auto& r_properties = this->GetProperties();

        // Defining the shape functions container
        const Matrix& NContainer = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);

        // Defining necessary variables
        array_1d<double, TNumNodes * TDim> DisplacementVector;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(DisplacementVector, r_geometry, DISPLACEMENT);
        BoundedMatrix<double, TDim, TDim> RotationMatrix;
        this->CalculateRotationMatrix(RotationMatrix, r_geometry);
        BoundedMatrix<double, TDim, TNumNodes * TDim> Nu = ZeroMatrix(TDim, TNumNodes * TDim);
        array_1d<double, TDim>                        LocalRelDispVector;
        array_1d<double, TDim>                        RelDispVector;
        const double& MinimumJointWidth = r_properties[MINIMUM_JOINT_WIDTH];
        double        JointWidth;
        const double& TransversalPermeability = r_properties[TRANSVERSAL_PERMEABILITY];
        BoundedMatrix<double, TDim, TDim> LocalPermeabilityMatrix = ZeroMatrix(TDim, TDim);

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            InterfaceElementUtilities::CalculateNuMatrix(Nu, NContainer, GPoint);

            noalias(RelDispVector) = prod(Nu, DisplacementVector);

            noalias(LocalRelDispVector) = prod(RotationMatrix, RelDispVector);

            this->CalculateJointWidth(JointWidth, LocalRelDispVector[TDim - 1], MinimumJointWidth, GPoint);

            InterfaceElementUtilities::FillPermeabilityMatrix(LocalPermeabilityMatrix, JointWidth,
                                                              TransversalPermeability);

            rOutput[GPoint].resize(TDim, TDim, false);
            noalias(rOutput[GPoint]) = LocalPermeabilityMatrix;
        }
    }

    KRATOS_CATCH("")
}

template <>
void UPwSmallStrainInterfaceElement<2, 4>::CalculateInitialGap(const GeometryType& rGeometry)
{
    KRATOS_TRY

    const double& MinimumJointWidth = this->GetProperties()[MINIMUM_JOINT_WIDTH];

    mInitialGap.resize(2);
    mIsOpen.resize(2);

    array_1d<double, 3> Vx;
    noalias(Vx)    = rGeometry.GetPoint(3) - rGeometry.GetPoint(0);
    mInitialGap[0] = norm_2(Vx);

    noalias(Vx)    = rGeometry.GetPoint(2) - rGeometry.GetPoint(1);
    mInitialGap[1] = norm_2(Vx);

    for (unsigned i = 0; i < mIsOpen.size(); ++i) {
        mIsOpen[i] = mInitialGap[i] >= MinimumJointWidth;
    }

    KRATOS_CATCH("")
}

template <>
void UPwSmallStrainInterfaceElement<3, 6>::CalculateInitialGap(const GeometryType& rGeometry)
{
    KRATOS_TRY

    const double& MinimumJointWidth = this->GetProperties()[MINIMUM_JOINT_WIDTH];

    mInitialGap.resize(3);
    mIsOpen.resize(3);

    array_1d<double, 3> Vx;
    noalias(Vx)    = rGeometry.GetPoint(3) - rGeometry.GetPoint(0);
    mInitialGap[0] = norm_2(Vx);

    noalias(Vx)    = rGeometry.GetPoint(4) - rGeometry.GetPoint(1);
    mInitialGap[1] = norm_2(Vx);

    noalias(Vx)    = rGeometry.GetPoint(5) - rGeometry.GetPoint(2);
    mInitialGap[2] = norm_2(Vx);

    for (unsigned i = 0; i < mIsOpen.size(); ++i) {
        mIsOpen[i] = mInitialGap[i] >= MinimumJointWidth;
    }

    KRATOS_CATCH("")
}

template <>
void UPwSmallStrainInterfaceElement<3, 8>::CalculateInitialGap(const GeometryType& rGeometry)
{
    KRATOS_TRY

    const double& MinimumJointWidth = this->GetProperties()[MINIMUM_JOINT_WIDTH];

    mInitialGap.resize(4);
    mIsOpen.resize(4);

    array_1d<double, 3> Vx;
    noalias(Vx)    = rGeometry.GetPoint(4) - rGeometry.GetPoint(0);
    mInitialGap[0] = norm_2(Vx);

    noalias(Vx)    = rGeometry.GetPoint(5) - rGeometry.GetPoint(1);
    mInitialGap[1] = norm_2(Vx);

    noalias(Vx)    = rGeometry.GetPoint(6) - rGeometry.GetPoint(2);
    mInitialGap[2] = norm_2(Vx);

    noalias(Vx)    = rGeometry.GetPoint(7) - rGeometry.GetPoint(3);
    mInitialGap[3] = norm_2(Vx);

    for (unsigned i = 0; i < mIsOpen.size(); ++i) {
        mIsOpen[i] = !(mInitialGap[i] < MinimumJointWidth);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateMaterialStiffnessMatrix(MatrixType& rStiffnessMatrix,
                                                                                       const ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int N_DOF = this->GetNumberOfDOF();

    // Resizing mass matrix
    if (rStiffnessMatrix.size1() != N_DOF) rStiffnessMatrix.resize(N_DOF, N_DOF, false);
    noalias(rStiffnessMatrix) = ZeroMatrix(N_DOF, N_DOF);

    // Previous definitions
    const PropertiesType&                           r_properties = this->GetProperties();
    const GeometryType&                             r_geometry   = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        r_geometry.IntegrationPoints(mThisIntegrationMethod);
    const unsigned int NumGPoints = IntegrationPoints.size();

    // Containers of variables at all integration points
    const Matrix& NContainer = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer =
        r_geometry.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
    GeometryType::JacobiansType JContainer(NumGPoints);
    r_geometry.Jacobian(JContainer, mThisIntegrationMethod);
    Vector detJContainer(NumGPoints);
    r_geometry.DeterminantOfJacobian(detJContainer, mThisIntegrationMethod);
    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(IntegrationPoints, detJContainer);

    // Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(r_geometry, r_properties, CurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    // Element variables
    InterfaceElementVariables Variables;
    this->InitializeElementVariables(Variables, r_geometry, r_properties, CurrentProcessInfo);

    ConstitutiveLawUtilities::SetConstitutiveParameters(ConstitutiveParameters, Variables.StrainVector,
                                                        Variables.ConstitutiveMatrix, Variables.Np,
                                                        Variables.GradNpT, Variables.F, Variables.detF);

    // Auxiliary variables
    const double&          MinimumJointWidth = r_properties[MINIMUM_JOINT_WIDTH];
    array_1d<double, TDim> RelDispVector;
    SFGradAuxVariables     SFGradAuxVars;

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        // Compute Np, StrainVector, JointWidth, GradNpT
        noalias(Variables.Np) = row(NContainer, GPoint);
        InterfaceElementUtilities::CalculateNuMatrix(Variables.Nu, NContainer, GPoint);
        noalias(RelDispVector)          = prod(Variables.Nu, Variables.DisplacementVector);
        noalias(Variables.StrainVector) = prod(Variables.RotationMatrix, RelDispVector);

        this->CheckAndCalculateJointWidth(Variables.JointWidth, ConstitutiveParameters,
                                          Variables.StrainVector[TDim - 1], MinimumJointWidth, GPoint);

        this->CalculateShapeFunctionsGradients<Matrix>(
            Variables.GradNpT, SFGradAuxVars, JContainer[GPoint], Variables.RotationMatrix,
            DN_DeContainer[GPoint], NContainer, Variables.JointWidth, GPoint);

        // Compute constitutive tensor
        noalias(Variables.StressVector) = mStressVector[GPoint];
        ConstitutiveParameters.SetStressVector(Variables.StressVector);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        Variables.IntegrationCoefficient = integration_coefficients[GPoint];

        // Compute stiffness matrix
        this->CalculateAndAddStiffnessMatrix(rStiffnessMatrix, Variables);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                                                   VectorType& rRightHandSideVector,
                                                                   const ProcessInfo& CurrentProcessInfo,
                                                                   bool CalculateStiffnessMatrixFlag,
                                                                   bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    // Previous definitions
    const PropertiesType&                           r_properties = this->GetProperties();
    const GeometryType&                             r_geometry   = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        r_geometry.IntegrationPoints(mThisIntegrationMethod);
    const unsigned int NumGPoints = IntegrationPoints.size();

    // Containers of variables at all integration points
    const Matrix& NContainer = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer =
        r_geometry.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
    GeometryType::JacobiansType JContainer(NumGPoints);
    r_geometry.Jacobian(JContainer, mThisIntegrationMethod);
    Vector detJContainer(NumGPoints);
    r_geometry.DeterminantOfJacobian(detJContainer, mThisIntegrationMethod);

    // Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(r_geometry, r_properties, CurrentProcessInfo);

    // stiffness matrix is needed to calculate Biot coefficient
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if (CalculateResidualVectorFlag) ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    // Element variables
    InterfaceElementVariables Variables;
    this->InitializeElementVariables(Variables, r_geometry, r_properties, CurrentProcessInfo);

    ConstitutiveLawUtilities::SetConstitutiveParameters(ConstitutiveParameters, Variables.StrainVector,
                                                        Variables.ConstitutiveMatrix, Variables.Np,
                                                        Variables.GradNpT, Variables.F, Variables.detF);

    // Auxiliary variables
    const double&          MinimumJointWidth       = r_properties[MINIMUM_JOINT_WIDTH];
    const double&          TransversalPermeability = r_properties[TRANSVERSAL_PERMEABILITY];
    array_1d<double, TDim> RelDispVector;
    SFGradAuxVariables     SFGradAuxVars;

    RetentionLaw::Parameters RetentionParameters(this->GetProperties());

    const bool hasBiotCoefficient = r_properties.Has(BIOT_COEFFICIENT);

    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(IntegrationPoints, detJContainer);

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        // Compute Np, StrainVector, JointWidth, GradNpT
        noalias(Variables.Np) = row(NContainer, GPoint);

        InterfaceElementUtilities::CalculateNuMatrix(Variables.Nu, NContainer, GPoint);

        noalias(RelDispVector) = prod(Variables.Nu, Variables.DisplacementVector);

        noalias(Variables.StrainVector) = prod(Variables.RotationMatrix, RelDispVector);

        this->CheckAndCalculateJointWidth(Variables.JointWidth, ConstitutiveParameters,
                                          Variables.StrainVector[TDim - 1], MinimumJointWidth, GPoint);

        this->CalculateShapeFunctionsGradients<Matrix>(
            Variables.GradNpT, SFGradAuxVars, JContainer[GPoint], Variables.RotationMatrix,
            DN_DeContainer[GPoint], NContainer, Variables.JointWidth, GPoint);

        // Compute BodyAcceleration and Permeability Matrix
        GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
            Variables.BodyAcceleration, NContainer, Variables.VolumeAcceleration, GPoint);

        InterfaceElementUtilities::FillPermeabilityMatrix(
            Variables.LocalPermeabilityMatrix, Variables.JointWidth, TransversalPermeability);

        // Compute constitutive tensor and stresses
        ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        ModifyInactiveElementStress(Variables.JointWidth, mStressVector[GPoint]);

        this->CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

        this->InitializeBiotCoefficients(Variables, hasBiotCoefficient);

        Variables.IntegrationCoefficient = integration_coefficients[GPoint];

        // Contributions to the left hand side
        if (CalculateStiffnessMatrixFlag) this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        // Contributions to the right hand side
        if (CalculateResidualVectorFlag)
            this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::InitializeElementVariables(InterfaceElementVariables& rVariables,
                                                                                 const GeometryType& rGeometry,
                                                                                 const PropertiesType& rProperties,
                                                                                 const ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    // Properties variables
    rVariables.IgnoreUndrained = rProperties[IGNORE_UNDRAINED];

    rVariables.DynamicViscosityInverse = 1.0 / rProperties[DYNAMIC_VISCOSITY];

    // ProcessInfo variables
    rVariables.VelocityCoefficient   = CurrentProcessInfo[VELOCITY_COEFFICIENT];
    rVariables.DtPressureCoefficient = CurrentProcessInfo[DT_PRESSURE_COEFFICIENT];

    // Nodal Variables
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rVariables.PressureVector[i]   = rGeometry[i].FastGetSolutionStepValue(WATER_PRESSURE);
        rVariables.DtPressureVector[i] = rGeometry[i].FastGetSolutionStepValue(DT_WATER_PRESSURE);
    }
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.DisplacementVector,
                                                                 rGeometry, DISPLACEMENT);
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.VelocityVector, rGeometry, VELOCITY);
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.VolumeAcceleration,
                                                                 rGeometry, VOLUME_ACCELERATION);

    // General Variables
    this->CalculateRotationMatrix(rVariables.RotationMatrix, rGeometry);
    InterfaceElementUtilities::CalculateVoigtVector(rVariables.VoigtVector);

    // Variables computed at each GP

    // Constitutive Law parameters
    rVariables.StressVector.resize(TDim, false);
    rVariables.StrainVector.resize(TDim, false);
    rVariables.ConstitutiveMatrix.resize(TDim, TDim, false);
    rVariables.Np.resize(TNumNodes, false);
    rVariables.GradNpT.resize(TNumNodes, TDim, false);
    rVariables.F.resize(TDim, TDim, false);
    rVariables.detF = 1.0;

    // Auxiliary variables
    noalias(rVariables.Nu)                      = ZeroMatrix(TDim, TNumNodes * TDim);
    noalias(rVariables.LocalPermeabilityMatrix) = ZeroMatrix(TDim, TDim);

    // Retention law
    rVariables.FluidPressure          = 0.0;
    rVariables.DegreeOfSaturation     = 1.0;
    rVariables.DerivativeOfSaturation = 0.0;
    rVariables.RelativePermeability   = 1.0;
    rVariables.BishopCoefficient      = 1.0;

    KRATOS_CATCH("")
}

template <>
void UPwSmallStrainInterfaceElement<2, 4>::CalculateRotationMatrix(BoundedMatrix<double, 2, 2>& rRotationMatrix,
                                                                   const GeometryType& rGeometry)
{
    KRATOS_TRY

    // Define mid-plane points for quadrilateral_interface_2d_4
    array_1d<double, 3> pmid0;
    array_1d<double, 3> pmid1;
    noalias(pmid0) = 0.5 * (rGeometry.GetPoint(0) + rGeometry.GetPoint(3));
    noalias(pmid1) = 0.5 * (rGeometry.GetPoint(1) + rGeometry.GetPoint(2));

    // Unitary vector in local x direction
    array_1d<double, 3> Vx;
    noalias(Vx)       = pmid1 - pmid0;
    double inv_norm_x = 1.0 / norm_2(Vx);
    Vx[0] *= inv_norm_x;
    Vx[1] *= inv_norm_x;

    // Rotation Matrix
    rRotationMatrix(0, 0) = Vx[0];
    rRotationMatrix(0, 1) = Vx[1];

    // We need to determine the unitary vector in local y direction pointing towards the TOP face of the joint

    // Unitary vector in local x direction (3D)
    array_1d<double, 3> Vx3D;
    Vx3D[0] = Vx[0];
    Vx3D[1] = Vx[1];
    Vx3D[2] = 0.0;

    // Unitary vector in local y direction (first option)
    array_1d<double, 3> Vy3D;
    Vy3D[0] = -Vx[1];
    Vy3D[1] = Vx[0];
    Vy3D[2] = 0.0;

    // Vector in global z direction (first option)
    array_1d<double, 3> Vz;
    MathUtils<double>::CrossProduct(Vz, Vx3D, Vy3D);

    // Vz must have the same sign as vector (0,0,1)
    if (Vz[2] > 0.0) {
        rRotationMatrix(1, 0) = -Vx[1];
        rRotationMatrix(1, 1) = Vx[0];
    } else {
        rRotationMatrix(1, 0) = Vx[1];
        rRotationMatrix(1, 1) = -Vx[0];
    }

    KRATOS_CATCH("")
}

template <>
void UPwSmallStrainInterfaceElement<3, 6>::CalculateRotationMatrix(BoundedMatrix<double, 3, 3>& rRotationMatrix,
                                                                   const GeometryType& rGeometry)
{
    KRATOS_TRY

    // Define mid-plane points for prism_interface_3d_6
    array_1d<double, 3> pmid0;
    array_1d<double, 3> pmid1;
    array_1d<double, 3> pmid2;
    noalias(pmid0) = 0.5 * (rGeometry.GetPoint(0) + rGeometry.GetPoint(3));
    noalias(pmid1) = 0.5 * (rGeometry.GetPoint(1) + rGeometry.GetPoint(4));
    noalias(pmid2) = 0.5 * (rGeometry.GetPoint(2) + rGeometry.GetPoint(5));

    // Unitary vector in local x direction
    array_1d<double, 3> Vx;
    noalias(Vx)     = pmid1 - pmid0;
    double inv_norm = 1.0 / norm_2(Vx);
    Vx[0] *= inv_norm;
    Vx[1] *= inv_norm;
    Vx[2] *= inv_norm;

    // Unitary vector in local z direction
    array_1d<double, 3> Vy;
    noalias(Vy) = pmid2 - pmid0;
    array_1d<double, 3> Vz;
    MathUtils<double>::CrossProduct(Vz, Vx, Vy);
    inv_norm = 1.0 / norm_2(Vz);
    Vz[0] *= inv_norm;
    Vz[1] *= inv_norm;
    Vz[2] *= inv_norm;

    // Unitary vector in local y direction
    MathUtils<double>::CrossProduct(Vy, Vz, Vx);

    // Rotation Matrix
    rRotationMatrix(0, 0) = Vx[0];
    rRotationMatrix(0, 1) = Vx[1];
    rRotationMatrix(0, 2) = Vx[2];

    rRotationMatrix(1, 0) = Vy[0];
    rRotationMatrix(1, 1) = Vy[1];
    rRotationMatrix(1, 2) = Vy[2];

    rRotationMatrix(2, 0) = Vz[0];
    rRotationMatrix(2, 1) = Vz[1];
    rRotationMatrix(2, 2) = Vz[2];

    KRATOS_CATCH("")
}

template <>
void UPwSmallStrainInterfaceElement<3, 8>::CalculateRotationMatrix(BoundedMatrix<double, 3, 3>& rRotationMatrix,
                                                                   const GeometryType& rGeometry)
{
    KRATOS_TRY

    // Define mid-plane points for hexahedra_interface_3d_8
    array_1d<double, 3> pmid0;
    array_1d<double, 3> pmid1;
    array_1d<double, 3> pmid2;
    noalias(pmid0) = 0.5 * (rGeometry.GetPoint(0) + rGeometry.GetPoint(4));
    noalias(pmid1) = 0.5 * (rGeometry.GetPoint(1) + rGeometry.GetPoint(5));
    noalias(pmid2) = 0.5 * (rGeometry.GetPoint(2) + rGeometry.GetPoint(6));

    // Unitary vector in local x direction
    array_1d<double, 3> Vx;
    noalias(Vx)     = pmid1 - pmid0;
    double inv_norm = 1.0 / norm_2(Vx);
    Vx[0] *= inv_norm;
    Vx[1] *= inv_norm;
    Vx[2] *= inv_norm;

    // Unitary vector in local z direction
    array_1d<double, 3> Vy;
    noalias(Vy) = pmid2 - pmid0;
    array_1d<double, 3> Vz;
    MathUtils<double>::CrossProduct(Vz, Vx, Vy);
    inv_norm = 1.0 / norm_2(Vz);
    Vz[0] *= inv_norm;
    Vz[1] *= inv_norm;
    Vz[2] *= inv_norm;

    // Unitary vector in local y direction
    MathUtils<double>::CrossProduct(Vy, Vz, Vx);

    // Rotation Matrix
    rRotationMatrix(0, 0) = Vx[0];
    rRotationMatrix(0, 1) = Vx[1];
    rRotationMatrix(0, 2) = Vx[2];

    rRotationMatrix(1, 0) = Vy[0];
    rRotationMatrix(1, 1) = Vy[1];
    rRotationMatrix(1, 2) = Vy[2];

    rRotationMatrix(2, 0) = Vz[0];
    rRotationMatrix(2, 1) = Vz[1];
    rRotationMatrix(2, 2) = Vz[2];

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateJointWidth(double& rJointWidth,
                                                                          const double& NormalRelDisp,
                                                                          const double& MinimumJointWidth,
                                                                          const unsigned int& GPoint)
{
    KRATOS_TRY

    rJointWidth = std::max(mInitialGap[GPoint] + NormalRelDisp, MinimumJointWidth);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CheckAndCalculateJointWidth(double& rJointWidth,
                                                                                  ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                                                                  double& rNormalRelDisp,
                                                                                  double MinimumJointWidth,
                                                                                  unsigned int GPoint)
{
    KRATOS_TRY

    rJointWidth = mInitialGap[GPoint] + rNormalRelDisp;

    // Ignore contact between interfaces
    rConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY);

    if (mIsOpen[GPoint]) {
        // Initially open joint
        if (rJointWidth < MinimumJointWidth) {
            // consider contact between interfaces
            rConstitutiveParameters.Reset(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY);
            rNormalRelDisp = rJointWidth - MinimumJointWidth;
            rJointWidth    = MinimumJointWidth;
        }
    } else {
        // Initially closed joint
        if (rJointWidth < 0.0) {
            // consider contact between interfaces
            rConstitutiveParameters.Reset(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY);
            rNormalRelDisp = rJointWidth;
            rJointWidth    = MinimumJointWidth;
        } else if (rJointWidth < MinimumJointWidth) {
            rJointWidth = MinimumJointWidth;
        }
    }

    KRATOS_CATCH("")
}

template <>
template <class TMatrixType>
void UPwSmallStrainInterfaceElement<2, 4>::CalculateShapeFunctionsGradients(TMatrixType& rGradNpT,
                                                                            SFGradAuxVariables& rAuxVariables,
                                                                            const Matrix& Jacobian,
                                                                            const BoundedMatrix<double, 2, 2>& RotationMatrix,
                                                                            const Matrix& DN_De,
                                                                            const Matrix& Ncontainer,
                                                                            const double& JointWidth,
                                                                            const unsigned int& GPoint)
{
    KRATOS_TRY

    // Quadrilateral_interface_2d_4
    rAuxVariables.GlobalCoordinatesGradients[0] = Jacobian(0, 0);
    rAuxVariables.GlobalCoordinatesGradients[1] = Jacobian(1, 0);
    noalias(rAuxVariables.LocalCoordinatesGradients) =
        prod(RotationMatrix, rAuxVariables.GlobalCoordinatesGradients);

    rGradNpT(0, 0) = DN_De(0, 0) / rAuxVariables.LocalCoordinatesGradients[0];
    rGradNpT(0, 1) = -Ncontainer(GPoint, 0) / JointWidth;
    rGradNpT(1, 0) = DN_De(1, 0) / rAuxVariables.LocalCoordinatesGradients[0];
    rGradNpT(1, 1) = -Ncontainer(GPoint, 1) / JointWidth;
    rGradNpT(2, 0) = DN_De(2, 0) / rAuxVariables.LocalCoordinatesGradients[0];
    rGradNpT(2, 1) = Ncontainer(GPoint, 2) / JointWidth;
    rGradNpT(3, 0) = DN_De(3, 0) / rAuxVariables.LocalCoordinatesGradients[0];
    rGradNpT(3, 1) = Ncontainer(GPoint, 3) / JointWidth;

    KRATOS_CATCH("")
}

template <>
template <class TMatrixType>
void UPwSmallStrainInterfaceElement<3, 6>::CalculateShapeFunctionsGradients(TMatrixType& rGradNpT,
                                                                            SFGradAuxVariables& rAuxVariables,
                                                                            const Matrix& Jacobian,
                                                                            const BoundedMatrix<double, 3, 3>& RotationMatrix,
                                                                            const Matrix& DN_De,
                                                                            const Matrix& Ncontainer,
                                                                            const double& JointWidth,
                                                                            const unsigned int& GPoint)
{
    KRATOS_TRY

    // Prism_interface_3d_6
    for (unsigned int i = 0; i < 6; ++i) {
        rAuxVariables.ShapeFunctionsNaturalGradientsMatrix(i, 0) = DN_De(i, 0);
        rAuxVariables.ShapeFunctionsNaturalGradientsMatrix(i, 1) = DN_De(i, 1);
    }

    rAuxVariables.GlobalCoordinatesGradients[0] = Jacobian(0, 0);
    rAuxVariables.GlobalCoordinatesGradients[1] = Jacobian(1, 0);
    rAuxVariables.GlobalCoordinatesGradients[2] = Jacobian(2, 0);
    noalias(rAuxVariables.LocalCoordinatesGradients) =
        prod(RotationMatrix, rAuxVariables.GlobalCoordinatesGradients);

    rAuxVariables.LocalCoordinatesGradientsMatrix(0, 0) = rAuxVariables.LocalCoordinatesGradients[0];
    rAuxVariables.LocalCoordinatesGradientsMatrix(1, 0) = rAuxVariables.LocalCoordinatesGradients[1];

    rAuxVariables.GlobalCoordinatesGradients[0] = Jacobian(0, 1);
    rAuxVariables.GlobalCoordinatesGradients[1] = Jacobian(1, 1);
    rAuxVariables.GlobalCoordinatesGradients[2] = Jacobian(2, 1);
    noalias(rAuxVariables.LocalCoordinatesGradients) =
        prod(RotationMatrix, rAuxVariables.GlobalCoordinatesGradients);

    rAuxVariables.LocalCoordinatesGradientsMatrix(0, 1) = rAuxVariables.LocalCoordinatesGradients[0];
    rAuxVariables.LocalCoordinatesGradientsMatrix(1, 1) = rAuxVariables.LocalCoordinatesGradients[1];

    double determinant;
    MathUtils<>::InvertMatrix2(rAuxVariables.LocalCoordinatesGradientsInvMatrix,
                               rAuxVariables.LocalCoordinatesGradientsMatrix, determinant);

    noalias(rAuxVariables.ShapeFunctionsGradientsMatrix) =
        prod(rAuxVariables.ShapeFunctionsNaturalGradientsMatrix, rAuxVariables.LocalCoordinatesGradientsInvMatrix);

    rGradNpT(0, 0) = rAuxVariables.ShapeFunctionsGradientsMatrix(0, 0);
    rGradNpT(0, 1) = rAuxVariables.ShapeFunctionsGradientsMatrix(0, 1);
    rGradNpT(0, 2) = -Ncontainer(GPoint, 0) / JointWidth;
    rGradNpT(1, 0) = rAuxVariables.ShapeFunctionsGradientsMatrix(1, 0);
    rGradNpT(1, 1) = rAuxVariables.ShapeFunctionsGradientsMatrix(1, 1);
    rGradNpT(1, 2) = -Ncontainer(GPoint, 1) / JointWidth;
    rGradNpT(2, 0) = rAuxVariables.ShapeFunctionsGradientsMatrix(2, 0);
    rGradNpT(2, 1) = rAuxVariables.ShapeFunctionsGradientsMatrix(2, 1);
    rGradNpT(2, 2) = -Ncontainer(GPoint, 2) / JointWidth;
    rGradNpT(3, 0) = rAuxVariables.ShapeFunctionsGradientsMatrix(3, 0);
    rGradNpT(3, 1) = rAuxVariables.ShapeFunctionsGradientsMatrix(3, 1);
    rGradNpT(3, 2) = Ncontainer(GPoint, 3) / JointWidth;
    rGradNpT(4, 0) = rAuxVariables.ShapeFunctionsGradientsMatrix(4, 0);
    rGradNpT(4, 1) = rAuxVariables.ShapeFunctionsGradientsMatrix(4, 1);
    rGradNpT(4, 2) = Ncontainer(GPoint, 4) / JointWidth;
    rGradNpT(5, 0) = rAuxVariables.ShapeFunctionsGradientsMatrix(5, 0);
    rGradNpT(5, 1) = rAuxVariables.ShapeFunctionsGradientsMatrix(5, 1);
    rGradNpT(5, 2) = Ncontainer(GPoint, 5) / JointWidth;

    KRATOS_CATCH("")
}

template <>
template <class TMatrixType>
void UPwSmallStrainInterfaceElement<3, 8>::CalculateShapeFunctionsGradients(TMatrixType& rGradNpT,
                                                                            SFGradAuxVariables& rAuxVariables,
                                                                            const Matrix& Jacobian,
                                                                            const BoundedMatrix<double, 3, 3>& RotationMatrix,
                                                                            const Matrix& DN_De,
                                                                            const Matrix& Ncontainer,
                                                                            const double& JointWidth,
                                                                            const unsigned int& GPoint)
{
    KRATOS_TRY

    // Hexahedral_interface_3d_8
    for (unsigned int i = 0; i < 8; ++i) {
        rAuxVariables.ShapeFunctionsNaturalGradientsMatrix(i, 0) = DN_De(i, 0);
        rAuxVariables.ShapeFunctionsNaturalGradientsMatrix(i, 1) = DN_De(i, 1);
    }

    rAuxVariables.GlobalCoordinatesGradients[0] = Jacobian(0, 0);
    rAuxVariables.GlobalCoordinatesGradients[1] = Jacobian(1, 0);
    rAuxVariables.GlobalCoordinatesGradients[2] = Jacobian(2, 0);
    noalias(rAuxVariables.LocalCoordinatesGradients) =
        prod(RotationMatrix, rAuxVariables.GlobalCoordinatesGradients);

    rAuxVariables.LocalCoordinatesGradientsMatrix(0, 0) = rAuxVariables.LocalCoordinatesGradients[0];
    rAuxVariables.LocalCoordinatesGradientsMatrix(1, 0) = rAuxVariables.LocalCoordinatesGradients[1];

    rAuxVariables.GlobalCoordinatesGradients[0] = Jacobian(0, 1);
    rAuxVariables.GlobalCoordinatesGradients[1] = Jacobian(1, 1);
    rAuxVariables.GlobalCoordinatesGradients[2] = Jacobian(2, 1);
    noalias(rAuxVariables.LocalCoordinatesGradients) =
        prod(RotationMatrix, rAuxVariables.GlobalCoordinatesGradients);

    rAuxVariables.LocalCoordinatesGradientsMatrix(0, 1) = rAuxVariables.LocalCoordinatesGradients[0];
    rAuxVariables.LocalCoordinatesGradientsMatrix(1, 1) = rAuxVariables.LocalCoordinatesGradients[1];

    double determinant;
    MathUtils<>::InvertMatrix2(rAuxVariables.LocalCoordinatesGradientsInvMatrix,
                               rAuxVariables.LocalCoordinatesGradientsMatrix, determinant);

    noalias(rAuxVariables.ShapeFunctionsGradientsMatrix) =
        prod(rAuxVariables.ShapeFunctionsNaturalGradientsMatrix, rAuxVariables.LocalCoordinatesGradientsInvMatrix);

    rGradNpT(0, 0) = rAuxVariables.ShapeFunctionsGradientsMatrix(0, 0);
    rGradNpT(0, 1) = rAuxVariables.ShapeFunctionsGradientsMatrix(0, 1);
    rGradNpT(0, 2) = -Ncontainer(GPoint, 0) / JointWidth;
    rGradNpT(1, 0) = rAuxVariables.ShapeFunctionsGradientsMatrix(1, 0);
    rGradNpT(1, 1) = rAuxVariables.ShapeFunctionsGradientsMatrix(1, 1);
    rGradNpT(1, 2) = -Ncontainer(GPoint, 1) / JointWidth;
    rGradNpT(2, 0) = rAuxVariables.ShapeFunctionsGradientsMatrix(2, 0);
    rGradNpT(2, 1) = rAuxVariables.ShapeFunctionsGradientsMatrix(2, 1);
    rGradNpT(2, 2) = -Ncontainer(GPoint, 2) / JointWidth;
    rGradNpT(3, 0) = rAuxVariables.ShapeFunctionsGradientsMatrix(3, 0);
    rGradNpT(3, 1) = rAuxVariables.ShapeFunctionsGradientsMatrix(3, 1);
    rGradNpT(3, 2) = -Ncontainer(GPoint, 3) / JointWidth;
    rGradNpT(4, 0) = rAuxVariables.ShapeFunctionsGradientsMatrix(4, 0);
    rGradNpT(4, 1) = rAuxVariables.ShapeFunctionsGradientsMatrix(4, 1);
    rGradNpT(4, 2) = Ncontainer(GPoint, 4) / JointWidth;
    rGradNpT(5, 0) = rAuxVariables.ShapeFunctionsGradientsMatrix(5, 0);
    rGradNpT(5, 1) = rAuxVariables.ShapeFunctionsGradientsMatrix(5, 1);
    rGradNpT(5, 2) = Ncontainer(GPoint, 5) / JointWidth;
    rGradNpT(6, 0) = rAuxVariables.ShapeFunctionsGradientsMatrix(6, 0);
    rGradNpT(6, 1) = rAuxVariables.ShapeFunctionsGradientsMatrix(6, 1);
    rGradNpT(6, 2) = Ncontainer(GPoint, 6) / JointWidth;
    rGradNpT(7, 0) = rAuxVariables.ShapeFunctionsGradientsMatrix(7, 0);
    rGradNpT(7, 1) = rAuxVariables.ShapeFunctionsGradientsMatrix(7, 1);
    rGradNpT(7, 2) = Ncontainer(GPoint, 7) / JointWidth;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
                                                                         InterfaceElementVariables& rVariables)
{
    KRATOS_TRY

    this->CalculateAndAddStiffnessMatrix(rLeftHandSideMatrix, rVariables);

    if (!rVariables.IgnoreUndrained) {
        this->CalculateAndAddCouplingMatrix(rLeftHandSideMatrix, rVariables);

        this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix, rVariables);

        this->CalculateAndAddPermeabilityMatrix(rLeftHandSideMatrix, rVariables);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateAndAddStiffnessMatrix(
    MatrixType& rLeftHandSideMatrix, const InterfaceElementVariables& rVariables)
{
    KRATOS_TRY

    const BoundedMatrix<double, TDim, TDim> dim_matrix = prod(
        trans(rVariables.RotationMatrix),
        BoundedMatrix<double, TDim, TDim>(prod(rVariables.ConstitutiveMatrix, rVariables.RotationMatrix)));
    const BoundedMatrix<double, TNumNodes * TDim, TDim> temp_matrix = prod(trans(rVariables.Nu), dim_matrix);
    const BoundedMatrix<double, TNumNodes * TDim, TNumNodes * TDim> stiffness_matrix =
        prod(temp_matrix, rVariables.Nu) * rVariables.IntegrationCoefficient;

    // Distribute stiffness block matrix into the elemental matrix
    GeoElementUtilities::AssembleUUBlockMatrix(rLeftHandSideMatrix, stiffness_matrix);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateAndAddCouplingMatrix(
    MatrixType& rLeftHandSideMatrix, const InterfaceElementVariables& rVariables)
{
    KRATOS_TRY

    const Matrix b_matrix = prod(rVariables.RotationMatrix, rVariables.Nu);

    const auto coupling_matrix = GeoTransportEquationUtilities::CalculateCouplingMatrix<TDim, TNumNodes>(
        b_matrix, rVariables.VoigtVector, rVariables.Np, rVariables.BiotCoefficient,
        rVariables.BishopCoefficient, rVariables.IntegrationCoefficient);

    GeoElementUtilities::AssembleUPBlockMatrix(rLeftHandSideMatrix, coupling_matrix);

    if (!rVariables.IgnoreUndrained) {
        const double SaturationCoefficient = rVariables.DegreeOfSaturation / rVariables.BishopCoefficient;
        const BoundedMatrix<double, TNumNodes, TNumNodes * TDim> transposed_coupling_matrix =
            PORE_PRESSURE_SIGN_FACTOR * SaturationCoefficient * rVariables.VelocityCoefficient *
            trans(coupling_matrix);

        GeoElementUtilities::AssemblePUBlockMatrix(rLeftHandSideMatrix, transposed_coupling_matrix);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateAndAddCompressibilityMatrix(
    MatrixType& rLeftHandSideMatrix, const InterfaceElementVariables& rVariables)
{
    KRATOS_TRY

    const auto compressibility_matrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix(
        rVariables.Np, rVariables.BiotModulusInverse, rVariables.IntegrationCoefficient);

    GeoElementUtilities::AssemblePPBlockMatrix(
        rLeftHandSideMatrix, compressibility_matrix * rVariables.DtPressureCoefficient * rVariables.JointWidth);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateAndAddPermeabilityMatrix(
    MatrixType& rLeftHandSideMatrix, const InterfaceElementVariables& rVariables)
{
    KRATOS_TRY

    const auto permeability_matrix = GeoTransportEquationUtilities::CalculatePermeabilityMatrix<TDim, TNumNodes>(
        rVariables.GradNpT, rVariables.DynamicViscosityInverse, rVariables.LocalPermeabilityMatrix,
        rVariables.RelativePermeability * rVariables.JointWidth, rVariables.IntegrationCoefficient);

    GeoElementUtilities::AssemblePPBlockMatrix(rLeftHandSideMatrix, permeability_matrix);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                                                         InterfaceElementVariables& rVariables,
                                                                         unsigned int GPoint)
{
    KRATOS_TRY

    this->CalculateAndAddStiffnessForce(rRightHandSideVector, rVariables, GPoint);

    this->CalculateAndAddMixBodyForce(rRightHandSideVector, rVariables);

    this->CalculateAndAddCouplingTerms(rRightHandSideVector, rVariables);

    if (!rVariables.IgnoreUndrained) {
        this->CalculateAndAddCompressibilityFlow(rRightHandSideVector, rVariables);

        this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);

        this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateAndAddStiffnessForce(
    VectorType& rRightHandSideVector, const InterfaceElementVariables& rVariables, unsigned int GPoint)
{
    KRATOS_TRY
    const BoundedMatrix<double, TNumNodes * TDim, TDim> u_dim_matrix =
        prod(trans(rVariables.Nu), trans(rVariables.RotationMatrix));

    const array_1d<double, TNumNodes * TDim> stiffness_force =
        -1.0 * prod(u_dim_matrix, mStressVector[GPoint]) * rVariables.IntegrationCoefficient;

    GeoElementUtilities::AssembleUBlockVector(rRightHandSideVector, stiffness_force);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector,
                                                                                  InterfaceElementVariables& rVariables)
{
    KRATOS_TRY

    this->CalculateSoilGamma(rVariables);

    const array_1d<double, TNumNodes * TDim> mix_body_force =
        prod(trans(rVariables.Nu), rVariables.SoilGamma) * rVariables.JointWidth * rVariables.IntegrationCoefficient;

    // Distribute body force block vector into elemental vector
    GeoElementUtilities::AssembleUBlockVector(rRightHandSideVector, mix_body_force);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateSoilGamma(InterfaceElementVariables& rVariables)
{
    KRATOS_TRY

    noalias(rVariables.SoilGamma) = GeoTransportEquationUtilities::CalculateSoilDensity(
                                        rVariables.DegreeOfSaturation, this->GetProperties()) *
                                    rVariables.BodyAcceleration;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector,
                                                                                   InterfaceElementVariables& rVariables)
{
    KRATOS_TRY

    const BoundedMatrix<double, TNumNodes * TDim, TDim> temp_matrix =
        prod(trans(rVariables.Nu), trans(rVariables.RotationMatrix));

    array_1d<double, TNumNodes * TDim> u_vector = prod(temp_matrix, rVariables.VoigtVector);

    const BoundedMatrix<double, TNumNodes * TDim, TNumNodes> coupling_matrix =
        -PORE_PRESSURE_SIGN_FACTOR * rVariables.BiotCoefficient * rVariables.BishopCoefficient *
        outer_prod(u_vector, rVariables.Np) * rVariables.IntegrationCoefficient;

    noalias(u_vector) = prod(coupling_matrix, rVariables.PressureVector); // finally it is the coupling force vector

    GeoElementUtilities::AssembleUBlockVector(rRightHandSideVector, u_vector);

    if (!rVariables.IgnoreUndrained) {
        const double SaturationCoefficient = rVariables.DegreeOfSaturation / rVariables.BishopCoefficient;
        const array_1d<double, TNumNodes> coupling_flow_vector =
            PORE_PRESSURE_SIGN_FACTOR * SaturationCoefficient *
            prod(trans(coupling_matrix), rVariables.VelocityVector);

        GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, coupling_flow_vector);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateAndAddCompressibilityFlow(
    VectorType& rRightHandSideVector, const InterfaceElementVariables& rVariables)
{
    KRATOS_TRY

    const auto compressibility_matrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix(
        rVariables.Np, rVariables.BiotModulusInverse, rVariables.IntegrationCoefficient);

    const array_1d<double, TNumNodes> compressibility_flow =
        -1.0 * prod(compressibility_matrix * rVariables.JointWidth, rVariables.DtPressureVector);

    GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, compressibility_flow);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateAndAddPermeabilityFlow(
    VectorType& rRightHandSideVector, const InterfaceElementVariables& rVariables)
{
    KRATOS_TRY

    const BoundedMatrix<double, TNumNodes, TDim> temp_matrix =
        prod(rVariables.GradNpT, rVariables.LocalPermeabilityMatrix);

    const BoundedMatrix<double, TNumNodes, TNumNodes> permeability_matrix =
        -PORE_PRESSURE_SIGN_FACTOR * rVariables.DynamicViscosityInverse *
        rVariables.RelativePermeability * prod(temp_matrix, trans(rVariables.GradNpT)) *
        rVariables.JointWidth * rVariables.IntegrationCoefficient;

    const array_1d<double, TNumNodes> permeability_flow =
        -1.0 * prod(permeability_matrix, rVariables.PressureVector);

    // Distribute permeability block vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, permeability_flow);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateAndAddFluidBodyFlow(
    VectorType& rRightHandSideVector, const InterfaceElementVariables& rVariables)
{
    KRATOS_TRY

    const BoundedMatrix<double, TNumNodes, TDim> temp_matrix =
        -PORE_PRESSURE_SIGN_FACTOR * prod(rVariables.GradNpT, rVariables.LocalPermeabilityMatrix) *
        rVariables.JointWidth * rVariables.IntegrationCoefficient;

    const array_1d<double, TNumNodes> fluid_body_flow =
        rVariables.DynamicViscosityInverse * this->GetProperties()[DENSITY_WATER] *
        rVariables.RelativePermeability * prod(temp_matrix, rVariables.BodyAcceleration);

    GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, fluid_body_flow);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateRetentionResponse(
    InterfaceElementVariables& rVariables, RetentionLaw::Parameters& rRetentionParameters, unsigned int GPoint)
{
    KRATOS_TRY

    rVariables.FluidPressure =
        GeoTransportEquationUtilities::CalculateFluidPressure(rVariables.Np, rVariables.PressureVector);
    rRetentionParameters.SetFluidPressure(rVariables.FluidPressure);

    rVariables.DegreeOfSaturation = mRetentionLawVector[GPoint]->CalculateSaturation(rRetentionParameters);
    rVariables.DerivativeOfSaturation =
        mRetentionLawVector[GPoint]->CalculateDerivativeOfSaturation(rRetentionParameters);
    rVariables.RelativePermeability =
        mRetentionLawVector[GPoint]->CalculateRelativePermeability(rRetentionParameters);
    rVariables.BishopCoefficient = mRetentionLawVector[GPoint]->CalculateBishopCoefficient(rRetentionParameters);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
double UPwSmallStrainInterfaceElement<TDim, TNumNodes>::CalculateBulkModulus(const Matrix& ConstitutiveMatrix)
{
    KRATOS_TRY

    const int    IndexM = ConstitutiveMatrix.size1() - 1;
    const double M      = ConstitutiveMatrix(IndexM, IndexM);
    const double G      = ConstitutiveMatrix(0, 0);

    return M - (4.0 / 3.0) * G;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::InitializeBiotCoefficients(InterfaceElementVariables& rVariables,
                                                                                 const bool& hasBiotCoefficient)
{
    KRATOS_TRY

    const auto& r_properties = this->GetProperties();

    // Properties variables
    if (hasBiotCoefficient) {
        rVariables.BiotCoefficient = r_properties[BIOT_COEFFICIENT];
    } else {
        // calculate Bulk modulus from stiffness matrix
        const double BulkModulus   = CalculateBulkModulus(rVariables.ConstitutiveMatrix);
        rVariables.BiotCoefficient = 1.0 - BulkModulus / r_properties[BULK_MODULUS_SOLID];
    }

    rVariables.BiotModulusInverse =
        (rVariables.BiotCoefficient - r_properties[POROSITY]) / r_properties[BULK_MODULUS_SOLID] +
        r_properties[POROSITY] / r_properties[BULK_MODULUS_FLUID];

    rVariables.BiotModulusInverse *= rVariables.DegreeOfSaturation;
    rVariables.BiotModulusInverse -= rVariables.DerivativeOfSaturation * r_properties[POROSITY];

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainInterfaceElement<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
}

template <>
void UPwSmallStrainInterfaceElement<2, 4>::InterpolateOutputDoubles(std::vector<double>& rOutput,
                                                                    const std::vector<double>& GPValues)
{
    KRATOS_TRY
#ifdef KRATOS_DEBUG
    KRATOS_ERROR_IF_NOT(rOutput.size() == 4) << "size of rOutput must be 4 " << this->Id() << std::endl;
#endif

    // Interpolation of computed values at Lobatto GP to the standard GiD gauss points
    rOutput[0] = 0.6220084679281462 * GPValues[0] + 0.16666666666666663 * GPValues[1] +
                 0.044658198738520435 * GPValues[1] + 0.16666666666666663 * GPValues[0];
    rOutput[1] = 0.16666666666666663 * GPValues[0] + 0.6220084679281462 * GPValues[1] +
                 0.16666666666666663 * GPValues[1] + 0.044658198738520435 * GPValues[0];
    rOutput[2] = 0.044658198738520435 * GPValues[0] + 0.16666666666666663 * GPValues[1] +
                 0.6220084679281462 * GPValues[1] + 0.16666666666666663 * GPValues[0];
    rOutput[3] = 0.16666666666666663 * GPValues[0] + 0.044658198738520435 * GPValues[1] +
                 0.16666666666666663 * GPValues[1] + 0.6220084679281462 * GPValues[0];

    KRATOS_CATCH("")
}

template <>
void UPwSmallStrainInterfaceElement<3, 6>::InterpolateOutputDoubles(std::vector<double>& rOutput,
                                                                    const std::vector<double>& GPValues)
{
    KRATOS_TRY

#ifdef KRATOS_DEBUG
    KRATOS_ERROR_IF_NOT(rOutput.size() == 6) << "size of rOutput must be 6 " << this->Id() << std::endl;
#endif

    // Interpolation of computed values at Lobatto GP to the standard GiD gauss points
    rOutput[0] = 0.5257834230632086 * GPValues[0] + 0.13144585576580214 * GPValues[1] +
                 0.13144585576580214 * GPValues[2] + 0.14088324360345805 * GPValues[0] +
                 0.03522081090086451 * GPValues[1] + 0.03522081090086451 * GPValues[2];

    rOutput[1] = 0.13144585576580214 * GPValues[0] + 0.5257834230632086 * GPValues[1] +
                 0.13144585576580214 * GPValues[2] + 0.03522081090086451 * GPValues[0] +
                 0.14088324360345805 * GPValues[1] + 0.03522081090086451 * GPValues[2];

    rOutput[2] = 0.13144585576580214 * GPValues[0] + 0.13144585576580214 * GPValues[1] +
                 0.5257834230632086 * GPValues[2] + 0.03522081090086451 * GPValues[0] +
                 0.03522081090086451 * GPValues[1] + 0.14088324360345805 * GPValues[2];

    rOutput[3] = 0.14088324360345805 * GPValues[0] + 0.03522081090086451 * GPValues[1] +
                 0.03522081090086451 * GPValues[2] + 0.5257834230632086 * GPValues[0] +
                 0.13144585576580214 * GPValues[1] + 0.13144585576580214 * GPValues[2];

    rOutput[4] = 0.03522081090086451 * GPValues[0] + 0.14088324360345805 * GPValues[1] +
                 0.03522081090086451 * GPValues[2] + 0.13144585576580214 * GPValues[0] +
                 0.5257834230632086 * GPValues[1] + 0.13144585576580214 * GPValues[2];

    rOutput[5] = 0.03522081090086451 * GPValues[0] + 0.03522081090086451 * GPValues[1] +
                 0.14088324360345805 * GPValues[2] + 0.13144585576580214 * GPValues[0] +
                 0.13144585576580214 * GPValues[1] + 0.5257834230632086 * GPValues[2];

    KRATOS_CATCH("")
}

template <>
void UPwSmallStrainInterfaceElement<3, 8>::InterpolateOutputDoubles(std::vector<double>& rOutput,
                                                                    const std::vector<double>& GPValues)
{
    KRATOS_TRY
#ifdef KRATOS_DEBUG
    KRATOS_ERROR_IF_NOT(rOutput.size() == 8) << "size of rOutput must be 8 " << this->Id() << std::endl;
#endif

    // Interpolation of computed values at Lobatto GP to the standard GiD gauss points
    rOutput[0] = 0.4905626121623441 * GPValues[0] + 0.13144585576580212 * GPValues[1] +
                 0.035220810900864506 * GPValues[2] + 0.13144585576580212 * GPValues[3] +
                 0.13144585576580212 * GPValues[0] + 0.035220810900864506 * GPValues[1] +
                 0.009437387837655926 * GPValues[2] + 0.035220810900864506 * GPValues[3];

    rOutput[1] = 0.13144585576580212 * GPValues[0] + 0.4905626121623441 * GPValues[1] +
                 0.13144585576580212 * GPValues[2] + 0.035220810900864506 * GPValues[3] +
                 0.035220810900864506 * GPValues[0] + 0.13144585576580212 * GPValues[1] +
                 0.035220810900864506 * GPValues[2] + 0.009437387837655926 * GPValues[3];

    rOutput[2] = 0.035220810900864506 * GPValues[0] + 0.13144585576580212 * GPValues[1] +
                 0.4905626121623441 * GPValues[2] + 0.13144585576580212 * GPValues[3] +
                 0.009437387837655926 * GPValues[0] + 0.035220810900864506 * GPValues[1] +
                 0.13144585576580212 * GPValues[2] + 0.035220810900864506 * GPValues[3];

    rOutput[3] = 0.13144585576580212 * GPValues[0] + 0.035220810900864506 * GPValues[1] +
                 0.13144585576580212 * GPValues[2] + 0.4905626121623441 * GPValues[3] +
                 0.035220810900864506 * GPValues[0] + 0.009437387837655926 * GPValues[1] +
                 0.035220810900864506 * GPValues[2] + 0.13144585576580212 * GPValues[3];

    rOutput[4] = 0.13144585576580212 * GPValues[0] + 0.035220810900864506 * GPValues[1] +
                 0.009437387837655926 * GPValues[2] + 0.035220810900864506 * GPValues[3] +
                 0.4905626121623441 * GPValues[0] + 0.13144585576580212 * GPValues[1] +
                 0.035220810900864506 * GPValues[2] + 0.13144585576580212 * GPValues[3];

    rOutput[5] = 0.035220810900864506 * GPValues[0] + 0.13144585576580212 * GPValues[1] +
                 0.035220810900864506 * GPValues[2] + 0.009437387837655926 * GPValues[3] +
                 0.13144585576580212 * GPValues[0] + 0.4905626121623441 * GPValues[1] +
                 0.13144585576580212 * GPValues[2] + 0.035220810900864506 * GPValues[3];

    rOutput[6] = 0.009437387837655926 * GPValues[0] + 0.035220810900864506 * GPValues[1] +
                 0.13144585576580212 * GPValues[2] + 0.035220810900864506 * GPValues[3] +
                 0.035220810900864506 * GPValues[0] + 0.13144585576580212 * GPValues[1] +
                 0.4905626121623441 * GPValues[2] + 0.13144585576580212 * GPValues[3];

    rOutput[7] = 0.035220810900864506 * GPValues[0] + 0.009437387837655926 * GPValues[1] +
                 0.035220810900864506 * GPValues[2] + 0.13144585576580212 * GPValues[3] +
                 0.13144585576580212 * GPValues[0] + 0.035220810900864506 * GPValues[1] +
                 0.13144585576580212 * GPValues[2] + 0.4905626121623441 * GPValues[3];

    KRATOS_CATCH("")
}

template <>
template <class TValueType>
void UPwSmallStrainInterfaceElement<2, 4>::InterpolateOutputValues(std::vector<TValueType>& rOutput,
                                                                   const std::vector<TValueType>& GPValues)
{
    KRATOS_TRY

#ifdef KRATOS_DEBUG
    KRATOS_ERROR_IF_NOT(rOutput.size() == 4) << "size of rOutput must be 4 " << this->Id() << std::endl;
#endif

    // Interpolation of computed values at Lobatto GP to the standard GiD gauss points
    noalias(rOutput[0]) = 0.6220084679281462 * GPValues[0] + 0.16666666666666663 * GPValues[1] +
                          0.044658198738520435 * GPValues[1] + 0.16666666666666663 * GPValues[0];

    noalias(rOutput[1]) = 0.16666666666666663 * GPValues[0] + 0.6220084679281462 * GPValues[1] +
                          0.16666666666666663 * GPValues[1] + 0.044658198738520435 * GPValues[0];

    noalias(rOutput[2]) = 0.044658198738520435 * GPValues[0] + 0.16666666666666663 * GPValues[1] +
                          0.6220084679281462 * GPValues[1] + 0.16666666666666663 * GPValues[0];

    noalias(rOutput[3]) = 0.16666666666666663 * GPValues[0] + 0.044658198738520435 * GPValues[1] +
                          0.16666666666666663 * GPValues[1] + 0.6220084679281462 * GPValues[0];

    KRATOS_CATCH("")
}

template <>
template <class TValueType>
void UPwSmallStrainInterfaceElement<3, 6>::InterpolateOutputValues(std::vector<TValueType>& rOutput,
                                                                   const std::vector<TValueType>& GPValues)
{
    KRATOS_TRY

#ifdef KRATOS_DEBUG
    KRATOS_ERROR_IF_NOT(rOutput.size() == 6) << "size of rOutput must be 6 " << this->Id() << std::endl;
#endif

    // Interpolation of computed values at Lobatto GP to the standard GiD gauss points
    noalias(rOutput[0]) = 0.5257834230632086 * GPValues[0] + 0.13144585576580214 * GPValues[1] +
                          0.13144585576580214 * GPValues[2] + 0.14088324360345805 * GPValues[0] +
                          0.03522081090086451 * GPValues[1] + 0.03522081090086451 * GPValues[2];

    noalias(rOutput[1]) = 0.13144585576580214 * GPValues[0] + 0.5257834230632086 * GPValues[1] +
                          0.13144585576580214 * GPValues[2] + 0.03522081090086451 * GPValues[0] +
                          0.14088324360345805 * GPValues[1] + 0.03522081090086451 * GPValues[2];

    noalias(rOutput[2]) = 0.13144585576580214 * GPValues[0] + 0.13144585576580214 * GPValues[1] +
                          0.5257834230632086 * GPValues[2] + 0.03522081090086451 * GPValues[0] +
                          0.03522081090086451 * GPValues[1] + 0.14088324360345805 * GPValues[2];

    noalias(rOutput[3]) = 0.14088324360345805 * GPValues[0] + 0.03522081090086451 * GPValues[1] +
                          0.03522081090086451 * GPValues[2] + 0.5257834230632086 * GPValues[0] +
                          0.13144585576580214 * GPValues[1] + 0.13144585576580214 * GPValues[2];

    noalias(rOutput[4]) = 0.03522081090086451 * GPValues[0] + 0.14088324360345805 * GPValues[1] +
                          0.03522081090086451 * GPValues[2] + 0.13144585576580214 * GPValues[0] +
                          0.5257834230632086 * GPValues[1] + 0.13144585576580214 * GPValues[2];

    noalias(rOutput[5]) = 0.03522081090086451 * GPValues[0] + 0.03522081090086451 * GPValues[1] +
                          0.14088324360345805 * GPValues[2] + 0.13144585576580214 * GPValues[0] +
                          0.13144585576580214 * GPValues[1] + 0.5257834230632086 * GPValues[2];

    KRATOS_CATCH("")
}

template <>
template <class TValueType>
void UPwSmallStrainInterfaceElement<3, 8>::InterpolateOutputValues(std::vector<TValueType>& rOutput,
                                                                   const std::vector<TValueType>& GPValues)
{
    KRATOS_TRY

#ifdef KRATOS_DEBUG
    KRATOS_ERROR_IF_NOT(rOutput.size() == 8) << "size of rOutput must be 8 " << this->Id() << std::endl;
#endif

    // Interpolation of computed values at Lobatto GP to the standard GiD gauss points
    noalias(rOutput[0]) = 0.4905626121623441 * GPValues[0] + 0.13144585576580212 * GPValues[1] +
                          0.035220810900864506 * GPValues[2] + 0.13144585576580212 * GPValues[3] +
                          0.13144585576580212 * GPValues[0] + 0.035220810900864506 * GPValues[1] +
                          0.009437387837655926 * GPValues[2] + 0.035220810900864506 * GPValues[3];

    noalias(rOutput[1]) = 0.13144585576580212 * GPValues[0] + 0.4905626121623441 * GPValues[1] +
                          0.13144585576580212 * GPValues[2] + 0.035220810900864506 * GPValues[3] +
                          0.035220810900864506 * GPValues[0] + 0.13144585576580212 * GPValues[1] +
                          0.035220810900864506 * GPValues[2] + 0.009437387837655926 * GPValues[3];

    noalias(rOutput[2]) = 0.035220810900864506 * GPValues[0] + 0.13144585576580212 * GPValues[1] +
                          0.4905626121623441 * GPValues[2] + 0.13144585576580212 * GPValues[3] +
                          0.009437387837655926 * GPValues[0] + 0.035220810900864506 * GPValues[1] +
                          0.13144585576580212 * GPValues[2] + 0.035220810900864506 * GPValues[3];

    noalias(rOutput[3]) = 0.13144585576580212 * GPValues[0] + 0.035220810900864506 * GPValues[1] +
                          0.13144585576580212 * GPValues[2] + 0.4905626121623441 * GPValues[3] +
                          0.035220810900864506 * GPValues[0] + 0.009437387837655926 * GPValues[1] +
                          0.035220810900864506 * GPValues[2] + 0.13144585576580212 * GPValues[3];

    noalias(rOutput[4]) = 0.13144585576580212 * GPValues[0] + 0.035220810900864506 * GPValues[1] +
                          0.009437387837655926 * GPValues[2] + 0.035220810900864506 * GPValues[3] +
                          0.4905626121623441 * GPValues[0] + 0.13144585576580212 * GPValues[1] +
                          0.035220810900864506 * GPValues[2] + 0.13144585576580212 * GPValues[3];

    noalias(rOutput[5]) = 0.035220810900864506 * GPValues[0] + 0.13144585576580212 * GPValues[1] +
                          0.035220810900864506 * GPValues[2] + 0.009437387837655926 * GPValues[3] +
                          0.13144585576580212 * GPValues[0] + 0.4905626121623441 * GPValues[1] +
                          0.13144585576580212 * GPValues[2] + 0.035220810900864506 * GPValues[3];

    noalias(rOutput[6]) = 0.009437387837655926 * GPValues[0] + 0.035220810900864506 * GPValues[1] +
                          0.13144585576580212 * GPValues[2] + 0.035220810900864506 * GPValues[3] +
                          0.035220810900864506 * GPValues[0] + 0.13144585576580212 * GPValues[1] +
                          0.4905626121623441 * GPValues[2] + 0.13144585576580212 * GPValues[3];

    noalias(rOutput[7]) = 0.035220810900864506 * GPValues[0] + 0.009437387837655926 * GPValues[1] +
                          0.035220810900864506 * GPValues[2] + 0.13144585576580212 * GPValues[3] +
                          0.13144585576580212 * GPValues[0] + 0.035220810900864506 * GPValues[1] +
                          0.13144585576580212 * GPValues[2] + 0.4905626121623441 * GPValues[3];

    KRATOS_CATCH("")
}

template void UPwSmallStrainInterfaceElement<2, 4>::CalculateShapeFunctionsGradients<BoundedMatrix<double, 4, 2>>(
    BoundedMatrix<double, 4, 2>&       rGradNpT,
    SFGradAuxVariables&                rAuxVariables,
    const Matrix&                      Jacobian,
    const BoundedMatrix<double, 2, 2>& RotationMatrix,
    const Matrix&                      DN_De,
    const Matrix&                      Ncontainer,
    const double&                      JointWidth,
    const unsigned int&                GPoint);

template void UPwSmallStrainInterfaceElement<2, 4>::CalculateShapeFunctionsGradients<Matrix>(
    Matrix&                            rGradNpT,
    SFGradAuxVariables&                rAuxVariables,
    const Matrix&                      Jacobian,
    const BoundedMatrix<double, 2, 2>& RotationMatrix,
    const Matrix&                      DN_De,
    const Matrix&                      Ncontainer,
    const double&                      JointWidth,
    const unsigned int&                GPoint);

template void UPwSmallStrainInterfaceElement<3, 6>::CalculateShapeFunctionsGradients<BoundedMatrix<double, 6, 3>>(
    BoundedMatrix<double, 6, 3>&       rGradNpT,
    SFGradAuxVariables&                rAuxVariables,
    const Matrix&                      Jacobian,
    const BoundedMatrix<double, 3, 3>& RotationMatrix,
    const Matrix&                      DN_De,
    const Matrix&                      Ncontainer,
    const double&                      JointWidth,
    const unsigned int&                GPoint);

template void UPwSmallStrainInterfaceElement<3, 6>::CalculateShapeFunctionsGradients<Matrix>(
    Matrix&                            rGradNpT,
    SFGradAuxVariables&                rAuxVariables,
    const Matrix&                      Jacobian,
    const BoundedMatrix<double, 3, 3>& RotationMatrix,
    const Matrix&                      DN_De,
    const Matrix&                      Ncontainer,
    const double&                      JointWidth,
    const unsigned int&                GPoint);

template void UPwSmallStrainInterfaceElement<3, 8>::CalculateShapeFunctionsGradients<BoundedMatrix<double, 8, 3>>(
    BoundedMatrix<double, 8, 3>&       rGradNpT,
    SFGradAuxVariables&                rAuxVariables,
    const Matrix&                      Jacobian,
    const BoundedMatrix<double, 3, 3>& RotationMatrix,
    const Matrix&                      DN_De,
    const Matrix&                      Ncontainer,
    const double&                      JointWidth,
    const unsigned int&                GPoint);

template void UPwSmallStrainInterfaceElement<3, 8>::CalculateShapeFunctionsGradients<Matrix>(
    Matrix&                            rGradNpT,
    SFGradAuxVariables&                rAuxVariables,
    const Matrix&                      Jacobian,
    const BoundedMatrix<double, 3, 3>& RotationMatrix,
    const Matrix&                      DN_De,
    const Matrix&                      Ncontainer,
    const double&                      JointWidth,
    const unsigned int&                GPoint);

template class UPwSmallStrainInterfaceElement<2, 4>;
template class UPwSmallStrainInterfaceElement<3, 6>;
template class UPwSmallStrainInterfaceElement<3, 8>;

template void UPwSmallStrainInterfaceElement<2, 4>::InterpolateOutputValues<array_1d<double, 3>>(
    std::vector<array_1d<double, 3>>& rOutput, const std::vector<array_1d<double, 3>>& GPValues);
template void UPwSmallStrainInterfaceElement<3, 6>::InterpolateOutputValues<array_1d<double, 3>>(
    std::vector<array_1d<double, 3>>& rOutput, const std::vector<array_1d<double, 3>>& GPValues);
template void UPwSmallStrainInterfaceElement<3, 8>::InterpolateOutputValues<array_1d<double, 3>>(
    std::vector<array_1d<double, 3>>& rOutput, const std::vector<array_1d<double, 3>>& GPValues);

template void UPwSmallStrainInterfaceElement<2, 4>::InterpolateOutputValues<Matrix>(
    std::vector<Matrix>& rOutput, const std::vector<Matrix>& GPValues);
template void UPwSmallStrainInterfaceElement<3, 6>::InterpolateOutputValues<Matrix>(
    std::vector<Matrix>& rOutput, const std::vector<Matrix>& GPValues);
template void UPwSmallStrainInterfaceElement<3, 8>::InterpolateOutputValues<Matrix>(
    std::vector<Matrix>& rOutput, const std::vector<Matrix>& GPValues);

} // namespace Kratos
