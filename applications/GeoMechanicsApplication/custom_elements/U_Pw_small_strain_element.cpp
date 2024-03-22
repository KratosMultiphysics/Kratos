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
#include "custom_elements/U_Pw_small_strain_element.hpp"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwSmallStrainElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                NodesArrayType const& ThisNodes,
                                                                PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwSmallStrainElement(NewId, this->GetGeometry().Create(ThisNodes),
                                                      pProperties, this->GetStressStatePolicy().Clone()));
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwSmallStrainElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                GeometryType::Pointer pGeom,
                                                                PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(
        new UPwSmallStrainElement(NewId, pGeom, pProperties, this->GetStressStatePolicy().Clone()));
}

template <unsigned int TDim, unsigned int TNumNodes>
int UPwSmallStrainElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Base class checks for positive area and Id > 0
    // Verify generic variables
    int ierr = UPwBaseElement<TDim, TNumNodes>::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;

    const PropertiesType& rProp = this->GetProperties();
    const GeometryType&   rGeom = this->GetGeometry();

    KRATOS_ERROR_IF(rGeom.DomainSize() < 1.0e-15)
        << "DomainSize < 1.0e-15 for the element " << this->Id() << std::endl;

    // Verify specific properties
    KRATOS_ERROR_IF_NOT(rProp.Has(IGNORE_UNDRAINED))
        << "IGNORE_UNDRAINED does not exist in the parameter list" << this->Id() << std::endl;

    if (!rProp[IGNORE_UNDRAINED]) {
        KRATOS_ERROR_IF(!rProp.Has(BULK_MODULUS_FLUID) || rProp[BULK_MODULUS_FLUID] < 0.0)
            << "BULK_MODULUS_FLUID has Key zero, is not defined or has an invalid value at element"
            << this->Id() << std::endl;

        KRATOS_ERROR_IF(!rProp.Has(DYNAMIC_VISCOSITY) || rProp[DYNAMIC_VISCOSITY] < 0.0)
            << "DYNAMIC_VISCOSITY has Key zero, is not defined or has an invalid value at element"
            << this->Id() << std::endl;

        KRATOS_ERROR_IF(!rProp.Has(PERMEABILITY_XX) || rProp[PERMEABILITY_XX] < 0.0)
            << "PERMEABILITY_XX has Key zero, is not defined or has an invalid value at element"
            << this->Id() << std::endl;

        KRATOS_ERROR_IF(!rProp.Has(PERMEABILITY_YY) || rProp[PERMEABILITY_YY] < 0.0)
            << "PERMEABILITY_YY has Key zero, is not defined or has an invalid value at element"
            << this->Id() << std::endl;

        KRATOS_ERROR_IF(!rProp.Has(PERMEABILITY_XY) || rProp[PERMEABILITY_XY] < 0.0)
            << "PERMEABILITY_XY has Key zero, is not defined or has an invalid value at element"
            << this->Id() << std::endl;

        if constexpr (TDim > 2) {
            KRATOS_ERROR_IF(!rProp.Has(PERMEABILITY_ZZ) || rProp[PERMEABILITY_ZZ] < 0.0)
                << "PERMEABILITY_ZZ has Key zero, is not defined or has an invalid value at element"
                << this->Id() << std::endl;

            KRATOS_ERROR_IF(!rProp.Has(PERMEABILITY_YZ) || rProp[PERMEABILITY_YZ] < 0.0)
                << "PERMEABILITY_YZ has Key zero, is not defined or has an invalid value at element"
                << this->Id() << std::endl;

            KRATOS_ERROR_IF(!rProp.Has(PERMEABILITY_ZX) || rProp[PERMEABILITY_ZX] < 0.0)
                << "PERMEABILITY_ZX has Key zero, is not defined or has an invalid value at element"
                << this->Id() << std::endl;
        }
    }

    // Verify that the constitutive law exists
    KRATOS_ERROR_IF_NOT(this->GetProperties().Has(CONSTITUTIVE_LAW))
        << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;

    // Verify that the constitutive law has the correct dimension
    const SizeType strainSize = this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    if (TDim == 2) {
        KRATOS_ERROR_IF_NOT(strainSize == VOIGT_SIZE_2D_PLANE_STRAIN)
            << "Wrong constitutive law used. This is a 2D element! expected "
               "strain size is "
            << VOIGT_SIZE_2D_PLANE_STRAIN << " But received: " << strainSize
            << " in element id: " << this->Id() << std::endl;
    } else {
        KRATOS_ERROR_IF_NOT(strainSize == VOIGT_SIZE_3D)
            << "Wrong constitutive law used. This is a 3D element! expected "
               "strain size is "
            << VOIGT_SIZE_3D << " But received: " << strainSize << " in element id: " << this->Id()
            << std::endl;
    }

    // Check constitutive law
    if (mConstitutiveLawVector.size() > 0) {
        return mConstitutiveLawVector[0]->Check(rProp, rGeom, rCurrentProcessInfo);
    }

    // Check retention law
    if (mRetentionLawVector.size() > 0) {
        return mRetentionLawVector[0]->Check(rProp, rCurrentProcessInfo);
    }

    return ierr;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (!mIsInitialised) this->Initialize(rCurrentProcessInfo);

    // Defining necessary variables
    const GeometryType& rGeom      = this->GetGeometry();
    const IndexType     NumGPoints = rGeom.IntegrationPointsNumber(mThisIntegrationMethod);

    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, this->GetProperties(), rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveParameters.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE); // Note: this is for nonlocal damage

    // Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    // Create general parameters of retention law
    RetentionLaw::Parameters RetentionParameters(this->GetProperties(), rCurrentProcessInfo);

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        // Compute Np, GradNpT, B and StrainVector
        this->CalculateKinematics(Variables, GPoint);

        // Compute infinitesimal strain
        this->CalculateStrain(Variables, GPoint);

        // Set Gauss points variables to constitutive law parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        // Initialize constitutive law
        noalias(Variables.StressVector) = mStressVector[GPoint];
        ConstitutiveParameters.SetStressVector(Variables.StressVector);
        mConstitutiveLawVector[GPoint]->InitializeMaterialResponseCauchy(ConstitutiveParameters);

        // Initialize retention law
        mRetentionLawVector[GPoint]->InitializeSolutionStep(RetentionParameters);
    }

    // Reset hydraulic discharge
    this->ResetHydraulicDischarge();

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::ResetHydraulicDischarge()
{
    KRATOS_TRY

    // Reset hydraulic discharge
    GeometryType& rGeom = this->GetGeometry();
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        ThreadSafeNodeWrite(rGeom[i], HYDRAULIC_DISCHARGE, 0.0);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateHydraulicDischarge(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    std::vector<array_1d<double, 3>> FluidFlux;
    this->CalculateOnIntegrationPoints(FLUID_FLUX_VECTOR, FluidFlux, rCurrentProcessInfo);

    const GeometryType& rGeom      = this->GetGeometry();
    const IndexType     NumGPoints = rGeom.IntegrationPointsNumber(mThisIntegrationMethod);
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        rGeom.IntegrationPoints(mThisIntegrationMethod);

    ElementVariables Variables;
    // Gradient of shape functions and determinant of Jacobian
    Variables.GradNpTInitialConfiguration.resize(TNumNodes, TDim, false);
    Variables.GradNpT.resize(TNumNodes, TDim, false);
    Variables.detJContainer.resize(NumGPoints, false);
    rGeom.ShapeFunctionsIntegrationPointsGradients(Variables.DN_DXContainer,
                                                   Variables.detJContainer, mThisIntegrationMethod);

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        noalias(Variables.GradNpT) = Variables.DN_DXContainer[GPoint];
        Variables.detJ             = Variables.detJContainer[GPoint];

        // Compute weighting coefficient for integration
        Variables.IntegrationCoefficient =
            this->CalculateIntegrationCoefficient(IntegrationPoints, GPoint, Variables.detJ);

        for (unsigned int node = 0; node < TNumNodes; ++node) {
            double HydraulicDischarge = 0;
            for (unsigned int iDir = 0; iDir < TDim; ++iDir) {
                HydraulicDischarge += Variables.GradNpT(node, iDir) * FluidFlux[GPoint][iDir];
            }

            HydraulicDischarge *= Variables.IntegrationCoefficient;
            HydraulicDischarge += rGeom[node].FastGetSolutionStepValue(HYDRAULIC_DISCHARGE);
            ThreadSafeNodeWrite(this->GetGeometry()[node], HYDRAULIC_DISCHARGE, HydraulicDischarge);
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Defining necessary variables
    const GeometryType& rGeom      = this->GetGeometry();
    const IndexType     NumGPoints = rGeom.IntegrationPointsNumber(mThisIntegrationMethod);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, this->GetProperties(), rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveParameters.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE); // Note: this is for nonlocal damage

    // Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        this->CalculateKinematics(Variables, GPoint);

        // Compute infinitesimal strain
        this->CalculateStrain(Variables, GPoint);

        // Set gauss points variables to constitutive law parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        // Compute Stress
        ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    this->InitializeNonLinearIteration(rCurrentProcessInfo);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    this->CalculateHydraulicDischarge(rCurrentProcessInfo);

    // Defining necessary variables
    const GeometryType& rGeom      = this->GetGeometry();
    const IndexType     NumGPoints = rGeom.IntegrationPointsNumber(mThisIntegrationMethod);

    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, this->GetProperties(), rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    // Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    // Create general parameters of retention law
    RetentionLaw::Parameters RetentionParameters(this->GetProperties(), rCurrentProcessInfo);

    Matrix StressContainer(NumGPoints, mStressVector[0].size());
    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        this->CalculateKinematics(Variables, GPoint);

        // Compute infinitesimal strain
        this->CalculateStrain(Variables, GPoint);

        // Set Gauss points variables to constitutive law parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        // Compute constitutive tensor and/or stresses
        noalias(Variables.StressVector) = mStressVector[GPoint];
        ConstitutiveParameters.SetStressVector(Variables.StressVector);
        mConstitutiveLawVector[GPoint]->FinalizeMaterialResponseCauchy(ConstitutiveParameters);
        mStateVariablesFinalized[GPoint] =
            mConstitutiveLawVector[GPoint]->GetValue(STATE_VARIABLES, mStateVariablesFinalized[GPoint]);

        // retention law
        mRetentionLawVector[GPoint]->FinalizeSolutionStep(RetentionParameters);

        if (rCurrentProcessInfo[NODAL_SMOOTHING])
            this->SaveGPStress(StressContainer, mStressVector[GPoint], GPoint);
    }
    if (rCurrentProcessInfo[NODAL_SMOOTHING]) this->ExtrapolateGPValues(StressContainer);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::SaveGPStress(Matrix&       rStressContainer,
                                                          const Vector& StressVector,
                                                          unsigned int  GPoint)
{
    KRATOS_TRY

    for (unsigned int i = 0; i < StressVector.size(); ++i) {
        rStressContainer(GPoint, i) = StressVector[i];
    }

    /* INFO: (Quadrilateral_2D_4 with GI_GAUSS_2)
     *
     *                      |S0-0 S1-0 S2-0|
     * rStressContainer =   |S0-1 S1-1 S2-1|
     *                      |S0-2 S1-2 S2-2|
     *                      |S0-3 S1-3 S2-3|
     *
     * S1-0 = S[1] at GP 0
     */

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::ExtrapolateGPValues(const Matrix& StressContainer)
{
    KRATOS_TRY

    array_1d<double, TNumNodes> DamageContainer;

    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode) {
        DamageContainer[iNode] = 0.0;
        DamageContainer[iNode] =
            mConstitutiveLawVector[iNode]->GetValue(DAMAGE_VARIABLE, DamageContainer[iNode]);
    }

    GeometryType&               rGeom = this->GetGeometry();
    const double&               Area  = rGeom.Area(); // In 3D this is volume
    array_1d<Vector, TNumNodes> NodalStressVector;    // List with stresses at each node
    array_1d<Matrix, TNumNodes> NodalStressTensor;

    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode) {
        NodalStressVector[iNode].resize(VoigtSize);
        NodalStressTensor[iNode].resize(StressTensorSize, StressTensorSize);
    }

    BoundedMatrix<double, TNumNodes, TNumNodes> ExtrapolationMatrix;
    this->CalculateExtrapolationMatrix(ExtrapolationMatrix);

    Matrix AuxNodalStress;
    AuxNodalStress.resize(TNumNodes, VoigtSize);
    noalias(AuxNodalStress) = prod(ExtrapolationMatrix, StressContainer);

    /* INFO:
     *
     *                  |S0-0 S1-0 S2-0|
     * AuxNodalStress = |S0-1 S1-1 S2-1|
     *                  |S0-2 S1-2 S2-2|
     *
     * S1-0 = S[1] at node 0
     */

    array_1d<double, TNumNodes> NodalDamage;
    noalias(NodalDamage) = prod(ExtrapolationMatrix, DamageContainer);

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        noalias(NodalStressVector[i]) = row(AuxNodalStress, i) * Area;
        noalias(NodalStressTensor[i]) = MathUtils<double>::StressVectorToTensor(NodalStressVector[i]);

        rGeom[i].SetLock();
        noalias(rGeom[i].FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)) += NodalStressTensor[i];
        rGeom[i].FastGetSolutionStepValue(NODAL_DAMAGE_VARIABLE) += NodalDamage[i] * Area;
        rGeom[i].FastGetSolutionStepValue(NODAL_AREA) += Area;
        rGeom[i].UnSetLock();
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::SetValuesOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                                          const std::vector<Vector>& rValues,
                                                                          const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == CAUCHY_STRESS_VECTOR) {
        KRATOS_ERROR_IF(rValues.size() != mStressVector.size())
            << "Unexpected number of values for "
               "UPwSmallStrainElement::SetValuesOnIntegrationPoints"
            << std::endl;
        std::copy(rValues.begin(), rValues.end(), mStressVector.begin());
    } else {
        KRATOS_ERROR_IF(rValues.size() < mConstitutiveLawVector.size())
            << "Insufficient number of values for "
               "UPwSmallStrainElement::SetValuesOnIntegrationPoints"
            << std::endl;
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            mConstitutiveLawVector[GPoint]->SetValue(rVariable, rValues[GPoint], rCurrentProcessInfo);
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                                                          std::vector<double>& rOutput,
                                                                          const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    // Defining necessary variables
    const GeometryType& rGeom      = this->GetGeometry();
    const IndexType     NumGPoints = rGeom.IntegrationPointsNumber(mThisIntegrationMethod);

    auto& r_prop = this->GetProperties();
    if (rOutput.size() != NumGPoints) rOutput.resize(NumGPoints);

    if (rVariable == VON_MISES_STRESS) {
        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            StressStrainUtilities EquivalentStress;
            rOutput[GPoint] = EquivalentStress.CalculateVonMisesStress(mStressVector[GPoint]);
        }
    } else if (rVariable == MEAN_EFFECTIVE_STRESS) {
        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            StressStrainUtilities EquivalentStress;
            rOutput[GPoint] = EquivalentStress.CalculateMeanStress(mStressVector[GPoint]);
        }
    } else if (rVariable == MEAN_STRESS) {
        std::vector<Vector> StressVector;
        CalculateOnIntegrationPoints(TOTAL_STRESS_VECTOR, StressVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            StressStrainUtilities EquivalentStress;
            rOutput[GPoint] = EquivalentStress.CalculateMeanStress(StressVector[GPoint]);
        }
    } else if (rVariable == ENGINEERING_VON_MISES_STRAIN) {
        std::vector<Vector> StrainVector;
        CalculateOnIntegrationPoints(ENGINEERING_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            StressStrainUtilities EquivalentStress;
            rOutput[GPoint] = EquivalentStress.CalculateVonMisesStrain(StrainVector[GPoint]);
        }
    } else if (rVariable == ENGINEERING_VOLUMETRIC_STRAIN) {
        std::vector<Vector> StrainVector;
        CalculateOnIntegrationPoints(ENGINEERING_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            StressStrainUtilities EquivalentStress;
            rOutput[GPoint] = EquivalentStress.CalculateTrace(StrainVector[GPoint]);
        }
    } else if (rVariable == GREEN_LAGRANGE_VON_MISES_STRAIN) {
        std::vector<Vector> StrainVector;
        CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            StressStrainUtilities EquivalentStress;
            rOutput[GPoint] = EquivalentStress.CalculateVonMisesStrain(StrainVector[GPoint]);
        }
    } else if (rVariable == GREEN_LAGRANGE_VOLUMETRIC_STRAIN) {
        std::vector<Vector> StrainVector;
        CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            StressStrainUtilities EquivalentStress;
            rOutput[GPoint] = EquivalentStress.CalculateTrace(StrainVector[GPoint]);
        }
    } else if (rVariable == DEGREE_OF_SATURATION || rVariable == EFFECTIVE_SATURATION ||
               rVariable == BISHOP_COEFFICIENT || rVariable == DERIVATIVE_OF_SATURATION ||
               rVariable == RELATIVE_PERMEABILITY) {
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        // create general parameters of retention law
        RetentionLaw::Parameters RetentionParameters(this->GetProperties(), rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            // Compute Np, GradNpT, B and StrainVector
            this->CalculateKinematics(Variables, GPoint);

            this->CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

            if (rVariable == DEGREE_OF_SATURATION) rOutput[GPoint] = Variables.DegreeOfSaturation;
            if (rVariable == EFFECTIVE_SATURATION) rOutput[GPoint] = Variables.EffectiveSaturation;
            if (rVariable == BISHOP_COEFFICIENT) rOutput[GPoint] = Variables.BishopCoefficient;
            if (rVariable == DERIVATIVE_OF_SATURATION)
                rOutput[GPoint] = Variables.DerivativeOfSaturation;
            if (rVariable == RELATIVE_PERMEABILITY)
                rOutput[GPoint] = Variables.RelativePermeability;
        }
    } else if (rVariable == HYDRAULIC_HEAD) {
        const PropertiesType& rProp = this->GetProperties();

        // Defining the shape functions, the Jacobian and the shape functions local gradients containers
        const Matrix& NContainer = rGeom.ShapeFunctionsValues(mThisIntegrationMethod);

        const auto NodalHydraulicHead =
            GeoElementUtilities::CalculateNodalHydraulicHeadFromWaterPressures(rGeom, rProp);

        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            double HydraulicHead = 0.0;
            for (unsigned int node = 0; node < TNumNodes; ++node)
                HydraulicHead += NContainer(GPoint, node) * NodalHydraulicHead[node];

            rOutput[GPoint] = HydraulicHead;
        }
    } else if (rVariable == CONFINED_STIFFNESS || rVariable == SHEAR_STIFFNESS) {
        size_t variable_index;
        if (rVariable == CONFINED_STIFFNESS) {
            if (TDim == 2) {
                variable_index = INDEX_2D_PLANE_STRAIN_XX;
            } else if (TDim == 3) {
                variable_index = INDEX_3D_XX;
            } else {
                KRATOS_ERROR << "CONFINED_STIFFNESS can not be retrieved for dim " << TDim
                             << " in element: " << this->Id() << std::endl;
            }
        } else if (rVariable == SHEAR_STIFFNESS) {
            if (TDim == 2) {
                variable_index = INDEX_2D_PLANE_STRAIN_XY;
            } else if (TDim == 3) {
                variable_index = INDEX_3D_XZ;
            } else {
                KRATOS_ERROR << "SHEAR_STIFFNESS can not be retrieved for dim " << TDim
                             << " in element: " << this->Id() << std::endl;
            }
        }

        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        if (rOutput.size() != mConstitutiveLawVector.size())
            rOutput.resize(mConstitutiveLawVector.size());

        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            const PropertiesType& rProp = this->GetProperties();

            ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, rProp, rCurrentProcessInfo);
            ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
            ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

            // Compute Np, GradNpT, B and StrainVector
            this->CalculateKinematics(Variables, GPoint);

            // Compute infinitesimal strain
            this->CalculateStrain(Variables, GPoint);

            // set Gauss points variables to constitutive law parameters
            this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

            // compute constitutive tensor and/or stresses
            noalias(Variables.StressVector) = mStressVector[GPoint];
            ConstitutiveParameters.SetStressVector(Variables.StressVector);

            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            rOutput[GPoint] = Variables.ConstitutiveMatrix(variable_index, variable_index);
        }
    } else if (r_prop.Has(rVariable)) {
        // Map initial material property to gauss points, as required for the output
        rOutput.clear();
        std::fill_n(std::back_inserter(rOutput), mConstitutiveLawVector.size(), r_prop.GetValue(rVariable));
    } else {
        if (rOutput.size() != mConstitutiveLawVector.size())
            rOutput.resize(mConstitutiveLawVector.size());

        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i) {
            rOutput[i] = mConstitutiveLawVector[i]->GetValue(rVariable, rOutput[i]);
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable, std::vector<array_1d<double, 3>>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& rGeom      = this->GetGeometry();
    const IndexType     NumGPoints = rGeom.IntegrationPointsNumber(mThisIntegrationMethod);
    if (rOutput.size() != NumGPoints) rOutput.resize(NumGPoints);

    if (rVariable == FLUID_FLUX_VECTOR) {
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        const PropertiesType& rProp = this->GetProperties();

        array_1d<double, TDim> GradPressureTerm;
        array_1d<double, TDim> FluidFlux;

        // Create general parameters of retention law
        RetentionLaw::Parameters RetentionParameters(this->GetProperties(), rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            // Compute Np, GradNpT, B and StrainVector
            this->CalculateKinematics(Variables, GPoint);

            // Compute strain, needed for update permeability
            this->CalculateStrain(Variables, GPoint);
            this->CalculatePermeabilityUpdateFactor(Variables);

            GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
                Variables.BodyAcceleration, Variables.NContainer, Variables.VolumeAcceleration, GPoint);
            Variables.FluidPressure = CalculateFluidPressure(Variables);

            RetentionParameters.SetFluidPressure(Variables.FluidPressure);

            Variables.RelativePermeability =
                mRetentionLawVector[GPoint]->CalculateRelativePermeability(RetentionParameters);

            noalias(GradPressureTerm) = prod(trans(Variables.GradNpT), Variables.PressureVector);

            noalias(GradPressureTerm) +=
                PORE_PRESSURE_SIGN_FACTOR * rProp[DENSITY_WATER] * Variables.BodyAcceleration;

            noalias(FluidFlux) = PORE_PRESSURE_SIGN_FACTOR * Variables.DynamicViscosityInverse *
                                 Variables.RelativePermeability * Variables.PermeabilityUpdateFactor *
                                 prod(Variables.PermeabilityMatrix, GradPressureTerm);

            GeoElementUtilities::FillArray1dOutput(rOutput[GPoint], FluidFlux);
        }
    } else {
        if (rOutput.size() != mConstitutiveLawVector.size())
            rOutput.resize(mConstitutiveLawVector.size());

        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i) {
            noalias(rOutput[i]) = ZeroVector(3);
            rOutput[i]          = mConstitutiveLawVector[i]->GetValue(rVariable, rOutput[i]);
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                                          std::vector<Vector>& rOutput,
                                                                          const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Defining necessary variables
    const GeometryType&   rGeom      = this->GetGeometry();
    const IndexType       NumGPoints = rGeom.IntegrationPointsNumber(mThisIntegrationMethod);
    const PropertiesType& rProp      = this->GetProperties();

    if (rOutput.size() != NumGPoints) rOutput.resize(NumGPoints);

    if (rVariable == CAUCHY_STRESS_VECTOR) {
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            if (rOutput[GPoint].size() != mStressVector[GPoint].size())
                rOutput[GPoint].resize(mStressVector[GPoint].size(), false);

            rOutput[GPoint] = mStressVector[GPoint];
        }
    } else if (rVariable == TOTAL_STRESS_VECTOR) {
        // Defining necessary variables
        ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, rProp, rCurrentProcessInfo);
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

        // Element variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        // Create general parameters of retention law
        RetentionLaw::Parameters RetentionParameters(this->GetProperties(), rCurrentProcessInfo);

        Vector VoigtVector(mStressVector[0].size());
        noalias(VoigtVector) = ZeroVector(VoigtVector.size());

        for (unsigned int i = 0; i < StressTensorSize; ++i)
            VoigtVector[i] = 1.0;

        const bool hasBiotCoefficient = rProp.Has(BIOT_COEFFICIENT);

        Vector TotalStressVector(mStressVector[0].size());

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            // Compute Np, GradNpT, B and StrainVector
            this->CalculateKinematics(Variables, GPoint);

            // Compute infinitesimal strain
            this->CalculateStrain(Variables, GPoint);

            // Set Gauss points variables to constitutive law parameters
            this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

            // Compute constitutive tensor and/or stresses
            noalias(Variables.StressVector) = mStressVector[GPoint];
            ConstitutiveParameters.SetStressVector(Variables.StressVector);
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            Variables.BiotCoefficient = CalculateBiotCoefficient(Variables, hasBiotCoefficient);

            this->CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

            noalias(TotalStressVector) = mStressVector[GPoint];
            noalias(TotalStressVector) += PORE_PRESSURE_SIGN_FACTOR * Variables.BiotCoefficient *
                                          Variables.BishopCoefficient * Variables.FluidPressure * VoigtVector;

            if (rOutput[GPoint].size() != TotalStressVector.size())
                rOutput[GPoint].resize(TotalStressVector.size(), false);

            rOutput[GPoint] = TotalStressVector;
        }
    } else if (rVariable == ENGINEERING_STRAIN_VECTOR) {
        // Definition of variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            noalias(Variables.Np) = row(Variables.NContainer, GPoint);

            Matrix J0, InvJ0;
            this->CalculateDerivativesOnInitialConfiguration(
                Variables.detJInitialConfiguration, J0, InvJ0, Variables.GradNpTInitialConfiguration, GPoint);

            // Calculating operator B
            this->CalculateBMatrix(Variables.B, Variables.GradNpTInitialConfiguration, Variables.Np);

            // Compute infinitesimal strain
            this->CalculateCauchyStrain(Variables);

            if (rOutput[GPoint].size() != Variables.StrainVector.size())
                rOutput[GPoint].resize(Variables.StrainVector.size(), false);

            rOutput[GPoint] = Variables.StrainVector;
        }
    } else if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR) {
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            // Compute element kinematics (Np, gradNpT, |J|, B, strains)
            this->CalculateKinematics(Variables, GPoint);

            this->CalculateStrain(Variables, GPoint);

            if (rOutput[GPoint].size() != Variables.StrainVector.size())
                rOutput[GPoint].resize(Variables.StrainVector.size(), false);

            rOutput[GPoint] = Variables.StrainVector;
        }
    } else if (rProp.Has(rVariable)) {
        // Map initial material property to Gauss points, as required for the output
        rOutput.clear();
        std::fill_n(std::back_inserter(rOutput), mConstitutiveLawVector.size(), rProp.GetValue(rVariable));
    } else {
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
            rOutput[i] = mConstitutiveLawVector[i]->GetValue(rVariable, rOutput[i]);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                                          std::vector<Matrix>& rOutput,
                                                                          const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& rGeom      = this->GetGeometry();
    const IndexType     NumGPoints = rGeom.IntegrationPointsNumber(mThisIntegrationMethod);

    if (rOutput.size() != NumGPoints) rOutput.resize(NumGPoints);

    if (rVariable == CAUCHY_STRESS_TENSOR) {
        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            rOutput[GPoint].resize(StressTensorSize, StressTensorSize, false);
            rOutput[GPoint] = MathUtils<double>::StressVectorToTensor(mStressVector[GPoint]);
        }
    } else if (rVariable == TOTAL_STRESS_TENSOR) {
        std::vector<Vector> StressVector;

        this->CalculateOnIntegrationPoints(TOTAL_STRESS_VECTOR, StressVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            if (rOutput[GPoint].size2() != TDim) rOutput[GPoint].resize(TDim, TDim, false);

            rOutput[GPoint] = MathUtils<double>::StressVectorToTensor(StressVector[GPoint]);
        }
    } else if (rVariable == ENGINEERING_STRAIN_TENSOR) {
        std::vector<Vector> StrainVector;

        CalculateOnIntegrationPoints(ENGINEERING_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            if (rOutput[GPoint].size2() != TDim) rOutput[GPoint].resize(TDim, TDim, false);

            rOutput[GPoint] = MathUtils<double>::StrainVectorToTensor(StrainVector[GPoint]);
        }
    } else if (rVariable == GREEN_LAGRANGE_STRAIN_TENSOR) {
        std::vector<Vector> StrainVector;

        CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            if (rOutput[GPoint].size2() != TDim) rOutput[GPoint].resize(TDim, TDim, false);

            rOutput[GPoint] = MathUtils<double>::StrainVectorToTensor(StrainVector[GPoint]);
        }
    } else if (rVariable == PERMEABILITY_MATRIX) {
        // If the permeability of the element is a given property
        BoundedMatrix<double, TDim, TDim> PermeabilityMatrix;
        GeoElementUtilities::FillPermeabilityMatrix(PermeabilityMatrix, this->GetProperties());

        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            rOutput[GPoint].resize(TDim, TDim, false);
            noalias(rOutput[GPoint]) = PermeabilityMatrix;
        }
    } else {
        if (rOutput.size() != mConstitutiveLawVector.size())
            rOutput.resize(mConstitutiveLawVector.size());

        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i) {
            rOutput[i].resize(TDim, TDim, false);
            noalias(rOutput[i]) = ZeroMatrix(TDim, TDim);
            rOutput[i]          = mConstitutiveLawVector[i]->GetValue(rVariable, rOutput[i]);
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateMaterialStiffnessMatrix(MatrixType& rStiffnessMatrix,
                                                                              const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const IndexType N_DOF = TNumNodes * (TDim + 1);

    if (rStiffnessMatrix.size1() != N_DOF) rStiffnessMatrix.resize(N_DOF, N_DOF, false);
    noalias(rStiffnessMatrix) = ZeroMatrix(N_DOF, N_DOF);

    // Previous definitions
    const PropertiesType&                           rProp = this->GetProperties();
    const GeometryType&                             rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        rGeom.IntegrationPoints(mThisIntegrationMethod);
    const IndexType NumGPoints = IntegrationPoints.size();

    // Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, rProp, rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        // Compute Np, GradNpT, B and StrainVector
        this->CalculateKinematics(Variables, GPoint);

        // Compute infinitesimal strain
        this->CalculateStrain(Variables, GPoint);

        // Set Gauss points variables to constitutive law parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        // Compute constitutive tensor
        noalias(Variables.StressVector) = mStressVector[GPoint];
        ConstitutiveParameters.SetStressVector(Variables.StressVector);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        // Compute weighting coefficient for integration
        Variables.IntegrationCoefficient =
            this->CalculateIntegrationCoefficient(IntegrationPoints, GPoint, Variables.detJ);

        // Compute stiffness matrix
        this->CalculateAndAddStiffnessMatrix(rStiffnessMatrix, Variables);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddGeometricStiffnessMatrix(
    MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables, unsigned int GPoint)
{
    KRATOS_TRY

    Matrix StressTensor = MathUtils<double>::StressVectorToTensor(mStressVector[GPoint]);
    Matrix ReducedKgMatrix =
        prod(rVariables.GradNpT,
             rVariables.IntegrationCoefficient * Matrix(prod(StressTensor, trans(rVariables.GradNpT)))); // to be optimized

    Matrix UUMatrix(TNumNodes * TDim, TNumNodes * TDim);
    noalias(UUMatrix) = ZeroMatrix(TNumNodes * TDim, TNumNodes * TDim);
    MathUtils<double>::ExpandAndAddReducedMatrix(UUMatrix, ReducedKgMatrix, TDim);

    // Distribute stiffness block matrix into the elemental matrix
    GeoElementUtilities::AssembleUUBlockMatrix(rLeftHandSideMatrix, UUMatrix);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix,
                                                                 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const IndexType N_DOF = this->GetNumberOfDOF();

    if (rMassMatrix.size1() != N_DOF) rMassMatrix.resize(N_DOF, N_DOF, false);
    noalias(rMassMatrix) = ZeroMatrix(N_DOF, N_DOF);

    const GeometryType&                             rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        rGeom.IntegrationPoints(this->mThisIntegrationMethod);
    const IndexType NumGPoints = IntegrationPoints.size();

    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    // Create general parameters of retention law
    RetentionLaw::Parameters RetentionParameters(this->GetProperties(), rCurrentProcessInfo);

    // Defining shape functions at all integration points
    // Defining necessary variables
    BoundedMatrix<double, TDim, TNumNodes * TDim> Nut = ZeroMatrix(TDim, TNumNodes * TDim);
    BoundedMatrix<double, TDim, TNumNodes * TDim> AuxDensityMatrix =
        ZeroMatrix(TDim, TNumNodes * TDim);
    BoundedMatrix<double, TDim, TDim> DensityMatrix = ZeroMatrix(TDim, TDim);

    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Nut, Variables.NContainer, GPoint);

        Matrix J0, InvJ0, DNu_DX0;
        this->CalculateDerivativesOnInitialConfiguration(Variables.detJInitialConfiguration, J0,
                                                         InvJ0, DNu_DX0, GPoint);

        // Calculating weighting coefficient for integration
        Variables.IntegrationCoefficientInitialConfiguration = this->CalculateIntegrationCoefficient(
            IntegrationPoints, GPoint, Variables.detJInitialConfiguration);

        CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

        this->CalculateSoilDensity(Variables);

        GeoElementUtilities::AssembleDensityMatrix(DensityMatrix, Variables.Density);

        noalias(AuxDensityMatrix) = prod(DensityMatrix, Nut);

        // Adding contribution to Mass matrix
        GeoElementUtilities::AssembleUUBlockMatrix(
            rMassMatrix, prod(trans(Nut), AuxDensityMatrix) * Variables.IntegrationCoefficientInitialConfiguration);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAll(MatrixType&        rLeftHandSideMatrix,
                                                          VectorType&        rRightHandSideVector,
                                                          const ProcessInfo& rCurrentProcessInfo,
                                                          bool CalculateStiffnessMatrixFlag,
                                                          bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    // Previous definitions
    const PropertiesType&                           rProp = this->GetProperties();
    const GeometryType&                             rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        rGeom.IntegrationPoints(this->mThisIntegrationMethod);
    const IndexType NumGPoints = IntegrationPoints.size();

    // Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, rProp, rCurrentProcessInfo);

    // Stiffness matrix is needed to calculate Biot coefficient
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if (CalculateResidualVectorFlag) ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    // Create general parameters of retention law
    RetentionLaw::Parameters RetentionParameters(this->GetProperties(), rCurrentProcessInfo);

    const bool hasBiotCoefficient = rProp.Has(BIOT_COEFFICIENT);

    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        // Compute GradNpT, B and StrainVector
        this->CalculateKinematics(Variables, GPoint);

        // Compute infinitesimal strain
        this->CalculateStrain(Variables, GPoint);

        // Set Gauss points variables to constitutive law parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        // Compute Nu and BodyAcceleration
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Variables.Nu, Variables.NContainer, GPoint);
        GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
            Variables.BodyAcceleration, Variables.NContainer, Variables.VolumeAcceleration, GPoint);

        // Compute constitutive tensor and stresses
        ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

        this->InitializeBiotCoefficients(Variables, hasBiotCoefficient);
        this->CalculatePermeabilityUpdateFactor(Variables);

        // Compute weighting coefficient for integration
        Variables.IntegrationCoefficient =
            this->CalculateIntegrationCoefficient(IntegrationPoints, GPoint, Variables.detJ);

        Variables.IntegrationCoefficientInitialConfiguration = this->CalculateIntegrationCoefficient(
            IntegrationPoints, GPoint, Variables.detJInitialConfiguration);

        // Contributions to the left hand side
        if (CalculateStiffnessMatrixFlag) this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        // Contributions to the right hand side
        if (CalculateResidualVectorFlag)
            this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
double UPwSmallStrainElement<TDim, TNumNodes>::CalculateBiotCoefficient(const ElementVariables& rVariables,
                                                                        bool hasBiotCoefficient) const
{
    KRATOS_TRY

    const PropertiesType& rProp = this->GetProperties();

    // Properties variables
    if (hasBiotCoefficient) {
        return rProp[BIOT_COEFFICIENT];
    } else {
        // Calculate Bulk modulus from stiffness matrix
        const double BulkModulus = CalculateBulkModulus(rVariables.ConstitutiveMatrix);
        return 1.0 - BulkModulus / rProp[BULK_MODULUS_SOLID];
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::InitializeBiotCoefficients(ElementVariables& rVariables,
                                                                        bool hasBiotCoefficient)
{
    KRATOS_TRY

    const PropertiesType& rProp = this->GetProperties();

    // Properties variables
    rVariables.BiotCoefficient = CalculateBiotCoefficient(rVariables, hasBiotCoefficient);

    if (!rProp[IGNORE_UNDRAINED]) {
        rVariables.BiotModulusInverse =
            (rVariables.BiotCoefficient - rProp[POROSITY]) / rProp[BULK_MODULUS_SOLID] +
            rProp[POROSITY] / rProp[BULK_MODULUS_FLUID];
    } else {
        rVariables.BiotModulusInverse =
            (rVariables.BiotCoefficient - rProp[POROSITY]) / rProp[BULK_MODULUS_SOLID] + rProp[POROSITY] / TINY;
    }

    rVariables.BiotModulusInverse *= rVariables.DegreeOfSaturation;
    rVariables.BiotModulusInverse -= rVariables.DerivativeOfSaturation * rProp[POROSITY];

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculatePermeabilityUpdateFactor(ElementVariables& rVariables)
{
    KRATOS_TRY

    const PropertiesType& rProp = this->GetProperties();

    if (rProp[PERMEABILITY_CHANGE_INVERSE_FACTOR] > 0.0) {
        const double          InverseCK = rProp[PERMEABILITY_CHANGE_INVERSE_FACTOR];
        StressStrainUtilities EquivalentStress;
        const double          epsV      = EquivalentStress.CalculateTrace(rVariables.StrainVector);
        const double          ePrevious = rProp[POROSITY] / (1.0 - rProp[POROSITY]);
        const double          eCurrent  = (1.0 + ePrevious) * std::exp(epsV) - 1.0;
        const double          permLog10 = (eCurrent - ePrevious) * InverseCK;
        rVariables.PermeabilityUpdateFactor = pow(10.0, permLog10);
    } else {
        rVariables.PermeabilityUpdateFactor = 1.0;
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
double UPwSmallStrainElement<TDim, TNumNodes>::CalculateBulkModulus(const Matrix& ConstitutiveMatrix) const
{
    KRATOS_TRY

    const int IndexG = ConstitutiveMatrix.size1() - 1;
    return ConstitutiveMatrix(0, 0) - (4.0 / 3.0) * ConstitutiveMatrix(IndexG, IndexG);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::InitializeElementVariables(ElementVariables& rVariables,
                                                                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Properties variables
    this->InitializeProperties(rVariables);

    // ProcessInfo variables
    rVariables.VelocityCoefficient   = rCurrentProcessInfo[VELOCITY_COEFFICIENT];
    rVariables.DtPressureCoefficient = rCurrentProcessInfo[DT_PRESSURE_COEFFICIENT];

    // Nodal Variables
    this->InitializeNodalDisplacementVariables(rVariables);
    this->InitializeNodalPorePressureVariables(rVariables);
    this->InitializeNodalVolumeAccelerationVariables(rVariables);

    // Variables computed at each GP
    rVariables.Nu = ZeroMatrix(TDim, TNumNodes * TDim);
    rVariables.Np.resize(TNumNodes, false);
    rVariables.GradNpT.resize(TNumNodes, TDim, false);
    rVariables.F = identity_matrix<double>(TDim);

    rVariables.detF = 1.0;

    // General Variables
    rVariables.VoigtVector = ZeroVector(VoigtSize);
    std::fill_n(rVariables.VoigtVector.begin(), StressTensorSize, 1.0);

    rVariables.B = ZeroMatrix(VoigtSize, TNumNodes * TDim);

    const GeometryType& rGeom      = this->GetGeometry();
    const IndexType     NumGPoints = rGeom.IntegrationPointsNumber(mThisIntegrationMethod);

    // Shape functions
    rVariables.NContainer = rGeom.ShapeFunctionsValues(mThisIntegrationMethod);

    // Gradient of shape functions and determinant of Jacobian
    rVariables.detJContainer.resize(NumGPoints, false);

    rGeom.ShapeFunctionsIntegrationPointsGradients(
        rVariables.DN_DXContainer, rVariables.detJContainer, mThisIntegrationMethod);

    // Constitutive Law parameters
    rVariables.StressVector.resize(VoigtSize, false);
    rVariables.StrainVector.resize(VoigtSize, false);
    rVariables.ConstitutiveMatrix.resize(VoigtSize, VoigtSize, false);

    // Auxiliary variables
    rVariables.UVoigtMatrix.resize(TNumNodes * TDim, VoigtSize, false);

    // Retention law
    rVariables.FluidPressure          = 0.0;
    rVariables.DegreeOfSaturation     = 1.0;
    rVariables.DerivativeOfSaturation = 0.0;
    rVariables.RelativePermeability   = 1.0;
    rVariables.BishopCoefficient      = 1.0;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateBMatrix(Matrix& rB, const Matrix& GradNpT, const Vector& Np)
{
    rB = this->GetStressStatePolicy().CalculateBMatrix(GradNpT, Np, this->GetGeometry());
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
                                                                ElementVariables& rVariables)
{
    KRATOS_TRY

    this->CalculateAndAddStiffnessMatrix(rLeftHandSideMatrix, rVariables);
    this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix, rVariables);

    if (!rVariables.IgnoreUndrained) {
        this->CalculateAndAddCouplingMatrix(rLeftHandSideMatrix, rVariables);
        this->CalculateAndAddPermeabilityMatrix(rLeftHandSideMatrix, rVariables);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix,
                                                                            ElementVariables& rVariables)
{
    KRATOS_TRY

    noalias(rVariables.UVoigtMatrix) = prod(trans(rVariables.B), rVariables.ConstitutiveMatrix);
    noalias(rVariables.UUMatrix) = prod(rVariables.UVoigtMatrix, rVariables.B) * rVariables.IntegrationCoefficient;

    // Distribute stiffness block matrix into the elemental matrix
    GeoElementUtilities::AssembleUUBlockMatrix(rLeftHandSideMatrix, rVariables.UUMatrix);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix,
                                                                           ElementVariables& rVariables)
{
    KRATOS_TRY

    noalias(rVariables.UVector) = prod(trans(rVariables.B), rVariables.VoigtVector);

    noalias(rVariables.UPMatrix) =
        PORE_PRESSURE_SIGN_FACTOR * rVariables.BiotCoefficient * rVariables.BishopCoefficient *
        outer_prod(rVariables.UVector, rVariables.Np) * rVariables.IntegrationCoefficient;

    // Distribute coupling block matrix into the elemental matrix
    GeoElementUtilities::AssembleUPBlockMatrix(rLeftHandSideMatrix, rVariables.UPMatrix);

    if (!rVariables.IgnoreUndrained) {
        const double SaturationCoefficient = rVariables.DegreeOfSaturation / rVariables.BishopCoefficient;
        noalias(rVariables.PUMatrix) = PORE_PRESSURE_SIGN_FACTOR * SaturationCoefficient *
                                       rVariables.VelocityCoefficient * trans(rVariables.UPMatrix);

        // Distribute transposed coupling block matrix into the elemental matrix
        GeoElementUtilities::AssemblePUBlockMatrix(rLeftHandSideMatrix, rVariables.PUMatrix);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateCompressibilityMatrix(
    BoundedMatrix<double, TNumNodes, TNumNodes>& rPMatrix, const ElementVariables& rVariables) const
{
    KRATOS_TRY

    noalias(rPMatrix) = -PORE_PRESSURE_SIGN_FACTOR * rVariables.DtPressureCoefficient *
                        rVariables.BiotModulusInverse * outer_prod(rVariables.Np, rVariables.Np) *
                        rVariables.IntegrationCoefficient;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix,
                                                                                  ElementVariables& rVariables)
{
    KRATOS_TRY

    this->CalculateCompressibilityMatrix(rVariables.PPMatrix, rVariables);

    // Distribute compressibility block matrix into the elemental matrix
    GeoElementUtilities::AssemblePPBlockMatrix(rLeftHandSideMatrix, rVariables.PPMatrix);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculatePermeabilityMatrix(
    BoundedMatrix<double, TNumNodes, TDim>&      rPDimMatrix,
    BoundedMatrix<double, TNumNodes, TNumNodes>& rPMatrix,
    const ElementVariables&                      rVariables) const
{
    KRATOS_TRY

    noalias(rPDimMatrix) = -PORE_PRESSURE_SIGN_FACTOR * prod(rVariables.GradNpT, rVariables.PermeabilityMatrix);

    noalias(rPMatrix) = rVariables.DynamicViscosityInverse * rVariables.RelativePermeability *
                        rVariables.PermeabilityUpdateFactor *
                        prod(rPDimMatrix, trans(rVariables.GradNpT)) * rVariables.IntegrationCoefficient;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix,
                                                                               ElementVariables& rVariables)
{
    KRATOS_TRY

    this->CalculatePermeabilityMatrix(rVariables.PDimMatrix, rVariables.PPMatrix, rVariables);

    // Distribute permeability block matrix into the elemental matrix
    GeoElementUtilities::AssemblePPBlockMatrix(rLeftHandSideMatrix, rVariables.PPMatrix);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                                                ElementVariables& rVariables,
                                                                unsigned int      GPoint)
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
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddStiffnessForce(VectorType& rRightHandSideVector,
                                                                           ElementVariables& rVariables,
                                                                           unsigned int GPoint)
{
    KRATOS_TRY

    noalias(rVariables.UVector) =
        -1.0 * prod(trans(rVariables.B), mStressVector[GPoint]) * rVariables.IntegrationCoefficient;

    // Distribute stiffness block vector into elemental vector
    GeoElementUtilities::AssembleUBlockVector(rRightHandSideVector, rVariables.UVector);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector,
                                                                         ElementVariables& rVariables)
{
    KRATOS_TRY

    this->CalculateSoilGamma(rVariables);

    noalias(rVariables.UVector) = prod(trans(rVariables.Nu), rVariables.SoilGamma) *
                                  rVariables.IntegrationCoefficientInitialConfiguration;

    // Distribute body force block vector into elemental vector
    GeoElementUtilities::AssembleUBlockVector(rRightHandSideVector, rVariables.UVector);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateSoilDensity(ElementVariables& rVariables)
{
    KRATOS_TRY

    rVariables.Density = rVariables.DegreeOfSaturation * rVariables.Porosity * rVariables.FluidDensity +
                         (1.0 - rVariables.Porosity) * rVariables.SolidDensity;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateSoilGamma(ElementVariables& rVariables)
{
    KRATOS_TRY

    this->CalculateSoilDensity(rVariables);

    noalias(rVariables.SoilGamma) = rVariables.Density * rVariables.BodyAcceleration;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector,
                                                                          ElementVariables& rVariables)
{
    KRATOS_TRY

    noalias(rVariables.UVector) = prod(trans(rVariables.B), rVariables.VoigtVector);

    noalias(rVariables.UPMatrix) =
        -PORE_PRESSURE_SIGN_FACTOR * rVariables.BiotCoefficient * rVariables.BishopCoefficient *
        outer_prod(rVariables.UVector, rVariables.Np) * rVariables.IntegrationCoefficient;

    noalias(rVariables.UVector) = prod(rVariables.UPMatrix, rVariables.PressureVector);

    // Distribute coupling block vector 1 into elemental vector
    GeoElementUtilities::AssembleUBlockVector(rRightHandSideVector, rVariables.UVector);

    if (!rVariables.IgnoreUndrained) {
        const double SaturationCoefficient = rVariables.DegreeOfSaturation / rVariables.BishopCoefficient;
        noalias(rVariables.PVector) = PORE_PRESSURE_SIGN_FACTOR * SaturationCoefficient *
                                      prod(trans(rVariables.UPMatrix), rVariables.VelocityVector);

        // Distribute coupling block vector 2 into elemental vector
        GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, rVariables.PVector);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateCompressibilityFlow(
    BoundedMatrix<double, TNumNodes, TNumNodes>& rPMatrix,
    array_1d<double, TNumNodes>&                 rPVector,
    const ElementVariables&                      rVariables) const
{
    KRATOS_TRY

    noalias(rPMatrix) = -PORE_PRESSURE_SIGN_FACTOR * rVariables.BiotModulusInverse *
                        outer_prod(rVariables.Np, rVariables.Np) * rVariables.IntegrationCoefficient;

    noalias(rPVector) = -prod(rPMatrix, rVariables.DtPressureVector);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector,
                                                                                ElementVariables& rVariables)
{
    KRATOS_TRY

    this->CalculateCompressibilityFlow(rVariables.PPMatrix, rVariables.PVector, rVariables);

    // Distribute compressibility block vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, rVariables.PVector);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculatePermeabilityFlow(
    BoundedMatrix<double, TNumNodes, TDim>&      rPDimMatrix,
    BoundedMatrix<double, TNumNodes, TNumNodes>& rPMatrix,
    array_1d<double, TNumNodes>&                 rPVector,
    const ElementVariables&                      rVariables) const
{
    KRATOS_TRY

    noalias(rPDimMatrix) = prod(rVariables.GradNpT, rVariables.PermeabilityMatrix);

    noalias(rPMatrix) = -PORE_PRESSURE_SIGN_FACTOR * rVariables.DynamicViscosityInverse *
                        rVariables.RelativePermeability * rVariables.PermeabilityUpdateFactor *
                        prod(rPDimMatrix, trans(rVariables.GradNpT)) * rVariables.IntegrationCoefficient;

    noalias(rPVector) = -prod(rPMatrix, rVariables.PressureVector);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector,
                                                                             ElementVariables& rVariables)
{
    KRATOS_TRY

    this->CalculatePermeabilityFlow(rVariables.PDimMatrix, rVariables.PPMatrix, rVariables.PVector, rVariables);

    // Distribute permeability block vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, rVariables.PVector);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateFluidBodyFlow(BoundedMatrix<double, TNumNodes, TDim>& rPDimMatrix,
                                                                    array_1d<double, TNumNodes>& rPVector,
                                                                    const ElementVariables& rVariables) const
{
    KRATOS_TRY

    noalias(rPDimMatrix) = prod(rVariables.GradNpT, rVariables.PermeabilityMatrix) * rVariables.IntegrationCoefficient;

    noalias(rPVector) = rVariables.DynamicViscosityInverse * rVariables.FluidDensity *
                        rVariables.RelativePermeability * rVariables.PermeabilityUpdateFactor *
                        prod(rPDimMatrix, rVariables.BodyAcceleration);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector,
                                                                          ElementVariables& rVariables)
{
    KRATOS_TRY

    this->CalculateFluidBodyFlow(rVariables.PDimMatrix, rVariables.PVector, rVariables);

    // Distribute fluid body flow block vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, rVariables.PVector);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateStrain(ElementVariables& rVariables, unsigned int GPoint)
{
    if (rVariables.UseHenckyStrain) {
        this->CalculateDeformationGradient(rVariables, GPoint);
        noalias(rVariables.StrainVector) = StressStrainUtilities::CalculateHenckyStrain(rVariables.F, VoigtSize);
    } else {
        this->CalculateCauchyStrain(rVariables);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateCauchyStrain(ElementVariables& rVariables)
{
    noalias(rVariables.StrainVector) = prod(rVariables.B, rVariables.DisplacementVector);
}

template <unsigned int TDim, unsigned int TNumNodes>
Vector UPwSmallStrainElement<TDim, TNumNodes>::CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient)
{
    return this->GetStressStatePolicy().CalculateGreenLagrangeStrain(rDeformationGradient);
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateDeformationGradient(ElementVariables& rVariables,
                                                                          unsigned int GPoint)
{
    KRATOS_TRY

    // Calculation of derivative of shape function with respect to reference
    // configuration derivative of shape function (displacement)
    Matrix J0, InvJ0, DNu_DX0;
    double detJ0;
    this->CalculateDerivativesOnInitialConfiguration(detJ0, J0, InvJ0, DNu_DX0, GPoint);

    // Calculating current Jacobian in order to find deformation gradient
    Matrix J, InvJ, DNu_DX;
    double detJ;
    this->CalculateJacobianOnCurrentConfiguration(detJ, J, InvJ, GPoint);

    KRATOS_ERROR_IF(detJ < 0.0) << "ERROR:: ELEMENT ID: " << this->Id() << " INVERTED. DETJ: " << detJ
                                << " nodes:" << this->GetGeometry() << std::endl;

    // Deformation gradient
    noalias(rVariables.F) = prod(J, InvJ0);
    rVariables.detF       = MathUtils<double>::Det(rVariables.F);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::InitializeNodalPorePressureVariables(ElementVariables& rVariables)
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();

    // Nodal variables
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rVariables.PressureVector[i]   = rGeom[i].FastGetSolutionStepValue(WATER_PRESSURE);
        rVariables.DtPressureVector[i] = rGeom[i].FastGetSolutionStepValue(DT_WATER_PRESSURE);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::InitializeNodalDisplacementVariables(ElementVariables& rVariables)
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();

    // Nodal variables
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.DisplacementVector, rGeom, DISPLACEMENT);
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.VelocityVector, rGeom, VELOCITY);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::InitializeNodalVolumeAccelerationVariables(ElementVariables& rVariables)
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();

    // Nodal variables
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.VolumeAcceleration,
                                                                 rGeom, VOLUME_ACCELERATION);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::InitializeProperties(ElementVariables& rVariables)
{
    KRATOS_TRY

    const PropertiesType& rProp = this->GetProperties();

    rVariables.IgnoreUndrained = rProp[IGNORE_UNDRAINED];
    rVariables.UseHenckyStrain = false;
    if (rProp.Has(USE_HENCKY_STRAIN)) rVariables.UseHenckyStrain = rProp[USE_HENCKY_STRAIN];

    rVariables.ConsiderGeometricStiffness = false;
    if (rProp.Has(CONSIDER_GEOMETRIC_STIFFNESS))
        rVariables.ConsiderGeometricStiffness = rProp[CONSIDER_GEOMETRIC_STIFFNESS];

    rVariables.DynamicViscosityInverse = 1.0 / rProp[DYNAMIC_VISCOSITY];
    rVariables.FluidDensity            = rProp[DENSITY_WATER];
    rVariables.SolidDensity            = rProp[DENSITY_SOLID];
    rVariables.Porosity                = rProp[POROSITY];
    GeoElementUtilities::FillPermeabilityMatrix(rVariables.PermeabilityMatrix, rProp);

    rVariables.PermeabilityUpdateFactor = 1.0;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateKinematics(ElementVariables& rVariables, unsigned int GPoint)
{
    KRATOS_TRY

    // Setting the vector of shape functions and the matrix of the shape functions global gradients
    noalias(rVariables.Np)      = row(rVariables.NContainer, GPoint);
    noalias(rVariables.GradNpT) = rVariables.DN_DXContainer[GPoint];

    rVariables.detJ = rVariables.detJContainer[GPoint];

    // Compute the deformation matrix B
    this->CalculateBMatrix(rVariables.B, rVariables.GradNpT, rVariables.Np);

    Matrix J0, InvJ0;
    this->CalculateDerivativesOnInitialConfiguration(rVariables.detJInitialConfiguration, J0, InvJ0,
                                                     rVariables.GradNpTInitialConfiguration, GPoint);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::SetConstitutiveParameters(ElementVariables& rVariables,
                                                                       ConstitutiveLaw::Parameters& rConstitutiveParameters)
{
    KRATOS_TRY

    rConstitutiveParameters.SetStrainVector(rVariables.StrainVector);
    rConstitutiveParameters.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
    rConstitutiveParameters.SetShapeFunctionsValues(rVariables.Np);
    rConstitutiveParameters.SetShapeFunctionsDerivatives(rVariables.GradNpT);
    rConstitutiveParameters.SetDeformationGradientF(rVariables.F);
    rConstitutiveParameters.SetDeterminantF(rVariables.detF);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::SetRetentionParameters(const ElementVariables& rVariables,
                                                                    RetentionLaw::Parameters& rRetentionParameters)
{
    KRATOS_TRY

    rRetentionParameters.SetFluidPressure(rVariables.FluidPressure);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
double UPwSmallStrainElement<TDim, TNumNodes>::CalculateFluidPressure(const ElementVariables& rVariables)
{
    KRATOS_TRY

    return inner_prod(rVariables.Np, rVariables.PressureVector);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateRetentionResponse(ElementVariables& rVariables,
                                                                        RetentionLaw::Parameters& rRetentionParameters,
                                                                        unsigned int GPoint)
{
    KRATOS_TRY

    rVariables.FluidPressure = CalculateFluidPressure(rVariables);
    SetRetentionParameters(rVariables, rRetentionParameters);

    rVariables.DegreeOfSaturation = mRetentionLawVector[GPoint]->CalculateSaturation(rRetentionParameters);
    rVariables.DerivativeOfSaturation =
        mRetentionLawVector[GPoint]->CalculateDerivativeOfSaturation(rRetentionParameters);
    rVariables.RelativePermeability =
        mRetentionLawVector[GPoint]->CalculateRelativePermeability(rRetentionParameters);
    rVariables.BishopCoefficient = mRetentionLawVector[GPoint]->CalculateBishopCoefficient(rRetentionParameters);
    rVariables.EffectiveSaturation =
        mRetentionLawVector[GPoint]->CalculateEffectiveSaturation(rRetentionParameters);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateExtrapolationMatrix(BoundedMatrix<double, TNumNodes, TNumNodes>& rExtrapolationMatrix)
{
    KRATOS_TRY

    KRATOS_ERROR << "undefined number of nodes in CalculateExtrapolationMatrix "
                    "... TNumNodes:"
                 << TNumNodes << " element: " << this->Id() << std::endl;

    KRATOS_CATCH("")
}

template <>
void UPwSmallStrainElement<2, 3>::CalculateExtrapolationMatrix(BoundedMatrix<double, 3, 3>& rExtrapolationMatrix)
{
    // The matrix contains the shape functions at each GP evaluated at each
    // node. Rows: nodes Columns: GP

    // Triangle_2d_3
    // GI_GAUSS_2
    rExtrapolationMatrix(0, 0) = 1.6666666666666666666;
    rExtrapolationMatrix(0, 1) = -0.33333333333333333333;
    rExtrapolationMatrix(0, 2) = -0.33333333333333333333;
    rExtrapolationMatrix(1, 0) = -0.33333333333333333333;
    rExtrapolationMatrix(1, 1) = 1.6666666666666666666;
    rExtrapolationMatrix(1, 2) = -0.33333333333333333333;
    rExtrapolationMatrix(2, 0) = -0.33333333333333333333;
    rExtrapolationMatrix(2, 1) = -0.33333333333333333333;
    rExtrapolationMatrix(2, 2) = 1.6666666666666666666;
}

template <>
void UPwSmallStrainElement<2, 4>::CalculateExtrapolationMatrix(BoundedMatrix<double, 4, 4>& rExtrapolationMatrix)
{
    // Quadrilateral_2d_4
    // GI_GAUSS_2
    rExtrapolationMatrix(0, 0) = 1.8660254037844386;
    rExtrapolationMatrix(0, 1) = -0.5;
    rExtrapolationMatrix(0, 2) = 0.13397459621556132;
    rExtrapolationMatrix(0, 3) = -0.5;
    rExtrapolationMatrix(1, 0) = -0.5;
    rExtrapolationMatrix(1, 1) = 1.8660254037844386;
    rExtrapolationMatrix(1, 2) = -0.5;
    rExtrapolationMatrix(1, 3) = 0.13397459621556132;
    rExtrapolationMatrix(2, 0) = 0.13397459621556132;
    rExtrapolationMatrix(2, 1) = -0.5;
    rExtrapolationMatrix(2, 2) = 1.8660254037844386;
    rExtrapolationMatrix(2, 3) = -0.5;
    rExtrapolationMatrix(3, 0) = -0.5;
    rExtrapolationMatrix(3, 1) = 0.13397459621556132;
    rExtrapolationMatrix(3, 2) = -0.5;
    rExtrapolationMatrix(3, 3) = 1.8660254037844386;
}

template <>
void UPwSmallStrainElement<3, 4>::CalculateExtrapolationMatrix(BoundedMatrix<double, 4, 4>& rExtrapolationMatrix)
{
    // Tetrahedra_3d_4
    // GI_GAUSS_2
    rExtrapolationMatrix(0, 0) = -0.309016988749894905;
    rExtrapolationMatrix(0, 1) = -0.3090169887498949046;
    rExtrapolationMatrix(0, 2) = -0.309016988749894905;
    rExtrapolationMatrix(0, 3) = 1.9270509662496847144;
    rExtrapolationMatrix(1, 0) = 1.9270509662496847144;
    rExtrapolationMatrix(1, 1) = -0.30901698874989490481;
    rExtrapolationMatrix(1, 2) = -0.3090169887498949049;
    rExtrapolationMatrix(1, 3) = -0.30901698874989490481;
    rExtrapolationMatrix(2, 0) = -0.30901698874989490473;
    rExtrapolationMatrix(2, 1) = 1.9270509662496847143;
    rExtrapolationMatrix(2, 2) = -0.3090169887498949049;
    rExtrapolationMatrix(2, 3) = -0.30901698874989490481;
    rExtrapolationMatrix(3, 0) = -0.3090169887498949048;
    rExtrapolationMatrix(3, 1) = -0.30901698874989490471;
    rExtrapolationMatrix(3, 2) = 1.9270509662496847143;
    rExtrapolationMatrix(3, 3) = -0.30901698874989490481;
}

template <>
void UPwSmallStrainElement<3, 8>::CalculateExtrapolationMatrix(BoundedMatrix<double, 8, 8>& rExtrapolationMatrix)
{
    // Hexahedra_3d_8
    // GI_GAUSS_2
    rExtrapolationMatrix(0, 0) = 2.549038105676658;
    rExtrapolationMatrix(0, 1) = -0.6830127018922192;
    rExtrapolationMatrix(0, 2) = 0.18301270189221927;
    rExtrapolationMatrix(0, 3) = -0.6830127018922192;
    rExtrapolationMatrix(0, 4) = -0.6830127018922192;
    rExtrapolationMatrix(0, 5) = 0.18301270189221927;
    rExtrapolationMatrix(0, 6) = -0.04903810567665795;
    rExtrapolationMatrix(0, 7) = 0.18301270189221927;

    rExtrapolationMatrix(1, 0) = -0.6830127018922192;
    rExtrapolationMatrix(1, 1) = 2.549038105676658;
    rExtrapolationMatrix(1, 2) = -0.6830127018922192;
    rExtrapolationMatrix(1, 3) = 0.18301270189221927;
    rExtrapolationMatrix(1, 4) = 0.18301270189221927;
    rExtrapolationMatrix(1, 5) = -0.6830127018922192;
    rExtrapolationMatrix(1, 6) = 0.18301270189221927;
    rExtrapolationMatrix(1, 7) = -0.04903810567665795;

    rExtrapolationMatrix(2, 0) = 0.18301270189221927;
    rExtrapolationMatrix(2, 1) = -0.6830127018922192;
    rExtrapolationMatrix(2, 2) = 2.549038105676658;
    rExtrapolationMatrix(2, 3) = -0.6830127018922192;
    rExtrapolationMatrix(2, 4) = -0.04903810567665795;
    rExtrapolationMatrix(2, 5) = 0.18301270189221927;
    rExtrapolationMatrix(2, 6) = -0.6830127018922192;
    rExtrapolationMatrix(2, 7) = 0.18301270189221927;

    rExtrapolationMatrix(3, 0) = -0.6830127018922192;
    rExtrapolationMatrix(3, 1) = 0.18301270189221927;
    rExtrapolationMatrix(3, 2) = -0.6830127018922192;
    rExtrapolationMatrix(3, 3) = 2.549038105676658;
    rExtrapolationMatrix(3, 4) = 0.18301270189221927;
    rExtrapolationMatrix(3, 5) = -0.04903810567665795;
    rExtrapolationMatrix(3, 6) = 0.18301270189221927;
    rExtrapolationMatrix(3, 7) = -0.6830127018922192;

    rExtrapolationMatrix(4, 0) = -0.6830127018922192;
    rExtrapolationMatrix(4, 1) = 0.18301270189221927;
    rExtrapolationMatrix(4, 2) = -0.04903810567665795;
    rExtrapolationMatrix(4, 3) = 0.18301270189221927;
    rExtrapolationMatrix(4, 4) = 2.549038105676658;
    rExtrapolationMatrix(4, 5) = -0.6830127018922192;
    rExtrapolationMatrix(4, 6) = 0.18301270189221927;
    rExtrapolationMatrix(4, 7) = -0.6830127018922192;

    rExtrapolationMatrix(5, 0) = 0.18301270189221927;
    rExtrapolationMatrix(5, 1) = -0.6830127018922192;
    rExtrapolationMatrix(5, 2) = 0.18301270189221927;
    rExtrapolationMatrix(5, 3) = -0.04903810567665795;
    rExtrapolationMatrix(5, 4) = -0.6830127018922192;
    rExtrapolationMatrix(5, 5) = 2.549038105676658;
    rExtrapolationMatrix(5, 6) = -0.6830127018922192;
    rExtrapolationMatrix(5, 7) = 0.18301270189221927;

    rExtrapolationMatrix(6, 0) = -0.04903810567665795;
    rExtrapolationMatrix(6, 1) = 0.18301270189221927;
    rExtrapolationMatrix(6, 2) = -0.6830127018922192;
    rExtrapolationMatrix(6, 3) = 0.18301270189221927;
    rExtrapolationMatrix(6, 4) = 0.18301270189221927;
    rExtrapolationMatrix(6, 5) = -0.6830127018922192;
    rExtrapolationMatrix(6, 6) = 2.549038105676658;
    rExtrapolationMatrix(6, 7) = -0.6830127018922192;

    rExtrapolationMatrix(7, 0) = 0.18301270189221927;
    rExtrapolationMatrix(7, 1) = -0.04903810567665795;
    rExtrapolationMatrix(7, 2) = 0.18301270189221927;
    rExtrapolationMatrix(7, 3) = -0.6830127018922192;
    rExtrapolationMatrix(7, 4) = -0.6830127018922192;
    rExtrapolationMatrix(7, 5) = 0.18301270189221927;
    rExtrapolationMatrix(7, 6) = -0.6830127018922192;
    rExtrapolationMatrix(7, 7) = 2.549038105676658;
}

template class UPwSmallStrainElement<2, 3>;
template class UPwSmallStrainElement<2, 4>;
template class UPwSmallStrainElement<3, 4>;
template class UPwSmallStrainElement<3, 8>;

template class UPwSmallStrainElement<2, 6>;
template class UPwSmallStrainElement<2, 8>;
template class UPwSmallStrainElement<2, 9>;
template class UPwSmallStrainElement<2, 10>;
template class UPwSmallStrainElement<2, 15>;
template class UPwSmallStrainElement<3, 10>;
template class UPwSmallStrainElement<3, 20>;
template class UPwSmallStrainElement<3, 27>;

} // Namespace Kratos
