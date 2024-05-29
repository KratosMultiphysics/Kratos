// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Athira Vadakkekkara
//  Collaborators:
//
//ToDo: Transform the calculated dHdEpseq_tension and dHdEpseq_compression in to global coordinates.
// Project includes
#include "elastic_anisotropic_damage_3d_nonlocal_equivalent_strains_tension_compression.h"
#include "constitutive_laws_application_variables.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "includes/checks.h"
#include "utilities/math_utils.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC()
    : ElasticIsotropic3D()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC(const ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC &rOther)
    : ElasticIsotropic3D(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::Clone() const
{
    return Kratos::make_shared<ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC>(ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC(*this));
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::~ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC()
{
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

int ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(POISSON_RATIO));
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS_TENSION));
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS_COMPRESSION));
    KRATOS_CHECK(rMaterialProperties.Has(DAMAGE_MODEL_PARAMETER_BETA1_TENSION));
    KRATOS_CHECK(rMaterialProperties.Has(DAMAGE_MODEL_PARAMETER_BETA2_TENSION));
    KRATOS_CHECK(rMaterialProperties.Has(DAMAGE_MODEL_PARAMETER_BETA1_COMPRESSION));
    KRATOS_CHECK(rMaterialProperties.Has(DAMAGE_MODEL_PARAMETER_BETA2_COMPRESSION));
    return 0;
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    mMaterialProperties = ZeroVector(6);
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double ft = rMaterialProperties[YIELD_STRESS_TENSION];
    const double nu = rMaterialProperties[POISSON_RATIO];
    mMaterialProperties[0] = (rMaterialProperties.Has(DAMAGE_THRESHOLD_TENSION)==true) ? rMaterialProperties[DAMAGE_THRESHOLD_TENSION] : ft/E;
    mMaterialProperties[1] = (rMaterialProperties.Has(DAMAGE_THRESHOLD_COMPRESSION)==true) ? rMaterialProperties[DAMAGE_THRESHOLD_COMPRESSION] : (1 + nu) * (20*ft)/(9*E);
    mMaterialProperties[2] = rMaterialProperties[DAMAGE_MODEL_PARAMETER_BETA1_TENSION];
    mMaterialProperties[3] = rMaterialProperties[DAMAGE_MODEL_PARAMETER_BETA2_TENSION];
    mMaterialProperties[4] = rMaterialProperties[DAMAGE_MODEL_PARAMETER_BETA1_COMPRESSION];
    mMaterialProperties[5] = rMaterialProperties[DAMAGE_MODEL_PARAMETER_BETA2_COMPRESSION];
}

//************************************************************************************
//************************************************************************************

Vector& ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValues
    )
{
    KRATOS_TRY

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValues,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::FinalizeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rParametersValues)
{
    KRATOS_TRY

    InternalDamageVariables rInternalDamageParameters;
    KRATOS_WATCH("finalize")
    this->CalculateStressResponse(rParametersValues, rInternalDamageParameters);
    mInternalDamageVariables.DamageMatrix = rInternalDamageParameters.DamageMatrix;
    mInternalDamageVariables.EquivalentStrainsTC = rInternalDamageParameters.EquivalentStrainsTC;
    mInternalDamageVariables.kappa = rInternalDamageParameters.kappa;
    mInternalDamageVariables.StrainHistoryParameter = rInternalDamageParameters.StrainHistoryParameter;
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rParametersValues)
{
    KRATOS_TRY

    InternalDamageVariables rInternalDamageVariables;
    this->CalculateStressResponse(rParametersValues, rInternalDamageVariables);
    //check whether internal parameters are set properly
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::CalculateStressResponse(
    ConstitutiveLaw::Parameters& rParametersValues,
    InternalDamageVariables& rInternalDamageVariables)
{
    KRATOS_TRY
    Flags& r_constitutive_law_options = rParametersValues.GetOptions();
    Vector& r_strain_vector = rParametersValues.GetStrainVector();
    CalculateValue(rParametersValues, STRAIN, r_strain_vector);

    // If we compute the tangent moduli or the stress
    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ||
        r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ))
        {
        Vector& r_stress_vector       = rParametersValues.GetStressVector();
        const Vector& r_strain_vector = rParametersValues.GetStrainVector();
        KRATOS_WATCH(r_strain_vector)
        Matrix& r_constitutive_matrix = rParametersValues.GetConstitutiveMatrix();
        Matrix r_elasticity_matrix;
        r_elasticity_matrix.resize(6, 6, false);
        CalculateElasticMatrix(r_elasticity_matrix, rParametersValues);
        noalias(r_stress_vector) = prod(r_elasticity_matrix, r_strain_vector);
        ModelParameters rModelParameters(Dimension);
        BoundedMatrixVoigtType EffStiffnessMatrix = ZeroMatrix(6, 6);
        BoundedMatrixVoigtType damage_effect_tensor = ZeroMatrix(6, 6);
        BoundedMatrixVoigtType transformed_damage_effect_tensor = ZeroMatrix(6, 6);
        BoundedMatrixVoigtType Inv_M  = ZeroMatrix(6,6);
        BoundedVectorVoigtType H_uNL1 = ZeroVector(6);
        BoundedVectorVoigtType H_uNL2 = ZeroVector(6);
        BoundedVectorVoigtType H_NL1u = ZeroVector(6);
        BoundedVectorVoigtType H_NL2u = ZeroVector(6);
        const double H_NLNL = 0.0;
        double det_M = 0.0;

        BoundedVectorType principal_strains = ZeroVector(3);
        BoundedVectorType principal_damage_vector;
        GetEigenValues(principal_strains, STRAIN, r_strain_vector);
        BoundedVectorType principal_damage_increment_vector;
        BoundedMatrixType damage_increment_tensor;
        BoundedMatrixType total_damage_tensor = mInternalDamageVariables.DamageMatrix;
        rModelParameters.kappa = mInternalDamageVariables.kappa;
        rModelParameters.strain_history_parameter = mInternalDamageVariables.StrainHistoryParameter;
        Vector local_equivalent_strains = ZeroVector(2);
        ComputePrincipalDamageIncrement(principal_damage_increment_vector, rParametersValues, rModelParameters, principal_strains);
        TransformPrincipalDamageIncrementToGlobal(damage_increment_tensor, principal_damage_increment_vector, r_strain_vector);
        total_damage_tensor += damage_increment_tensor;
        GetTotalPrincipalDamageVector(principal_damage_vector,total_damage_tensor);
        if(principal_damage_increment_vector[0] > 0.0 || principal_damage_increment_vector[1] > 0.0 || principal_damage_increment_vector[2] > 0.0){
            GetDamageEffectTensor(damage_effect_tensor, principal_damage_vector);
            GetTransformedDamageEffectTensor(transformed_damage_effect_tensor, damage_effect_tensor, r_strain_vector);
            MathUtils<double>::InvertMatrix(transformed_damage_effect_tensor, Inv_M, det_M);
            const BoundedMatrixVoigtType temp = prod(r_elasticity_matrix,trans(Inv_M));
            noalias(r_elasticity_matrix) = prod(Inv_M, temp);
            noalias(r_stress_vector)  = prod(r_elasticity_matrix, r_strain_vector);
            Calculate_tangent_HuNL(H_uNL1, H_uNL2, rParametersValues, rModelParameters, principal_damage_vector, principal_strains);
            Calculate_tangent_HNLu(H_NL1u, H_NL2u, rParametersValues, principal_strains);
            KRATOS_WATCH("there is damage")
        }
        AssembleConstitutiveMatrix(r_constitutive_matrix, r_elasticity_matrix, H_NL1u, H_NL2u, H_uNL1, H_uNL2, H_NLNL);
        rInternalDamageVariables.DamageMatrix = total_damage_tensor;
        CalculateLocalEquivalentStrains(local_equivalent_strains, rParametersValues, principal_strains);
        rInternalDamageVariables.EquivalentStrainsTC = local_equivalent_strains;
        rInternalDamageVariables.kappa = rModelParameters.kappa ;
        rInternalDamageVariables.StrainHistoryParameter = rModelParameters.strain_history_parameter;
        KRATOS_WATCH(r_stress_vector)
        }
    KRATOS_WATCH(rParametersValues.GetProcessInfo()[TIME])
    KRATOS_WATCH("----------------------------------------------------------------------------------------")
    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::GetEigenValues(
    BoundedVectorType& Pri_Values,
    const Variable<Vector>& rThisVariable,
    const Vector& VectorForm)
{
    KRATOS_TRY
    BoundedMatrixType MatrixForm = ZeroMatrix(3,3);
    BoundedMatrixType EigenVectors;
    BoundedMatrixType EigenValues = ZeroMatrix(3,3);
    VectorToTensor(MatrixForm, VectorForm, rThisVariable);
    //prinicpal values, max and min
    MathUtils<double>::GaussSeidelEigenSystem(MatrixForm, EigenVectors, EigenValues);
    Pri_Values[0] = EigenValues(0,0);
    Pri_Values[1] = EigenValues(1,1);
    Pri_Values[2] = EigenValues(2,2);
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::ComputePrincipalDamageIncrement(
    BoundedVectorType& PrincipalDamageIncrement,
    ConstitutiveLaw::Parameters& rParametersValues,
    ModelParameters& rModelParameters,
    const BoundedVectorType& PrincipalStrains
    )
{
    KRATOS_TRY

    Vector nonlocal_equivalent_strains = ZeroVector(2);
    BoundedVectorType principal_nonlocal_strains;
    BoundedVectorType principal_damage_n = ZeroVector(3);
    BoundedMatrixType MatrixForm = ZeroMatrix(3,3);
    BoundedMatrixType EigenVectors;
    BoundedMatrixType EigenValues = ZeroMatrix(3,3);
    Vector local_equivalent_strains_tc = ZeroVector(2);
    MathUtils<double>::GaussSeidelEigenSystem(mInternalDamageVariables.DamageMatrix, EigenVectors, EigenValues);
    principal_damage_n[0] = EigenValues(0,0);
    principal_damage_n[1] = EigenValues(1,1);
    principal_damage_n[2] = EigenValues(2,2);
    CalculateLocalEquivalentStrains(local_equivalent_strains_tc, rParametersValues, PrincipalStrains);
    const auto& N = rParametersValues.GetShapeFunctionsValues();
    for (SizeType i = 0; i < N.size(); ++i) {
        nonlocal_equivalent_strains[0] += N[i] * rParametersValues.GetElementGeometry()[i].FastGetSolutionStepValue(NON_LOCAL_EQUIVALENT_STRAIN_TENSION, 0);
    }
    for (SizeType i = 0; i < N.size(); ++i) {
        nonlocal_equivalent_strains[1] += N[i] * rParametersValues.GetElementGeometry()[i].FastGetSolutionStepValue(NON_LOCAL_EQUIVALENT_STRAIN_COMPRESSION, 0);
    }
    KRATOS_WATCH(nonlocal_equivalent_strains)
    ScaleNonlocalEquivalentStrain(principal_nonlocal_strains, rParametersValues, PrincipalStrains, nonlocal_equivalent_strains, local_equivalent_strains_tc);
    KRATOS_WATCH(local_equivalent_strains_tc)
    SetModelParameters(rModelParameters, rParametersValues, PrincipalStrains, principal_nonlocal_strains);
    KRATOS_WATCH(rModelParameters.kappa0)
    KRATOS_WATCH(rModelParameters.strain_history_parameter)
    KRATOS_WATCH(principal_nonlocal_strains)
    CheckDamageCriteria(PrincipalDamageIncrement, rModelParameters, principal_damage_n);
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::CheckDamageCriteria(
    BoundedVectorType& PrincipalDamageIncrement,
    const ModelParameters& rModelParameters,
    const BoundedVectorType& PrincipalDamageVector_n
    )
{
    KRATOS_TRY
    for (SizeType i = 0; i < Dimension; ++i) {
        if (rModelParameters.kappa[i]>= 0 && rModelParameters.kappa[i]<=rModelParameters.kappa0[i]){
            PrincipalDamageIncrement[i]=0.0;
            KRATOS_WATCH("no damage")
        }else if (rModelParameters.kappa[i]> rModelParameters.kappa0[i]){
            PrincipalDamageIncrement[i] = (1-PrincipalDamageVector_n[i]) * ((rModelParameters.Beta1[i]/rModelParameters.kappa[i]) + (rModelParameters.Beta2[i]/rModelParameters.kappa0[i])) * rModelParameters.del_kappa[i];
            KRATOS_WATCH("damage")
            KRATOS_WATCH(rModelParameters.kappa)
            KRATOS_WATCH(rModelParameters.del_kappa)
            KRATOS_WATCH(PrincipalDamageVector_n)
        }
        if(PrincipalDamageIncrement[i] < 0.0){
            PrincipalDamageIncrement[i] = 0.0;
        }

    }
    KRATOS_WATCH(PrincipalDamageIncrement)
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::SetModelParameters(
    ModelParameters& rModelParameters,
    ConstitutiveLaw::Parameters& rParametersValues,
    const BoundedVectorType& PrincipalStrains,
    const BoundedVectorType& PrincipalNonlocalStrains
    )
{
    KRATOS_TRY

    BoundedVectorType current_max_strain_history_parameter = ZeroVector(3);
    BoundedVectorType previous_kappa = rModelParameters.kappa;
    for(SizeType i = 0; i < Dimension; ++i){
        rModelParameters.Beta1[i] = (PrincipalStrains[i] > eps) ? mMaterialProperties[2] : mMaterialProperties[4];
        rModelParameters.Beta2[i] = (PrincipalStrains[i] > eps) ? mMaterialProperties[3] : mMaterialProperties[5];
        rModelParameters.kappa0[i] = (PrincipalStrains[i] > eps) ? mMaterialProperties[0] : mMaterialProperties[1];
        current_max_strain_history_parameter[i] = std::max(rModelParameters.strain_history_parameter[i], PrincipalNonlocalStrains[i]);
        rModelParameters.strain_history_parameter[i] = current_max_strain_history_parameter[i];
        rModelParameters.kappa[i] = std::max(current_max_strain_history_parameter[i], rModelParameters.kappa0[i]);
        rModelParameters.del_kappa[i] = rModelParameters.kappa[i] - previous_kappa[i];
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::CalculateLocalEquivalentStrains(
    Vector& LocalEquivalentStrains,
    ConstitutiveLaw::Parameters& rParametersValues,
    const BoundedVectorType& PrincipalStrains
    )
{
    KRATOS_TRY
    Vector& r_strain_vector = rParametersValues.GetStrainVector();
    BoundedVectorType PrincipalDeviatoricStrains = ZeroVector(3);
    GetPrincipalDeviatoricStrains(PrincipalDeviatoricStrains, r_strain_vector, PrincipalStrains);
    KRATOS_WATCH(PrincipalDeviatoricStrains)
    LocalEquivalentStrains[0] = sqrt(pow(MacaulayBrackets(PrincipalStrains[0]),2) + pow(MacaulayBrackets(PrincipalStrains[1]),2) + pow(MacaulayBrackets(PrincipalStrains[2]),2));
    LocalEquivalentStrains[1] = sqrt(pow(MacaulayBrackets(PrincipalDeviatoricStrains[0]),2) + pow(MacaulayBrackets(PrincipalDeviatoricStrains[1]),2) + pow(MacaulayBrackets(PrincipalDeviatoricStrains[2]),2));
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::GetPrincipalDeviatoricStrains(
    BoundedVectorType& PrincipalDeviatoricStrains,
    const Vector& StrainVector,
    const BoundedVectorType& PrincipalStrains
    )
{
    KRATOS_TRY
    const double trace = StrainVector[0] + StrainVector[1] + StrainVector[2];
    const double hydrostatic_strain = trace/3;
    for (SizeType i = 0; i < Dimension; ++i){
        PrincipalDeviatoricStrains[i] = PrincipalStrains[i] - hydrostatic_strain;
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//***********************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::TransformPrincipalDamageIncrementToGlobal(
    BoundedMatrixType& DamageIncrementTensor,
    const BoundedVectorType& PrinicipalDamageIncrementVector,
    const Vector& StrainVector
    )
{
    KRATOS_TRY
    BoundedVectorType l = ZeroVector(3);
    BoundedVectorType m = ZeroVector(3);
    BoundedVectorType n = ZeroVector(3);
    BoundedMatrixType MatrixForm = ZeroMatrix(3,3);
    BoundedMatrixType EigenVectors;
    BoundedMatrixType EigenValues = ZeroMatrix(3,3);
    VectorToTensor(MatrixForm, StrainVector, STRAIN);
    MathUtils<double>::GaussSeidelEigenSystem(MatrixForm, EigenVectors, EigenValues);
        for(SizeType i = 0; i < Dimension; ++i){
            l[i] = EigenVectors(i,0);
            m[i] = EigenVectors(i,1);
            n[i] = EigenVectors(i,2);
        }
    for(SizeType i = 0; i < Dimension; ++i){
        DamageIncrementTensor(i,i) = PrinicipalDamageIncrementVector[0] * l[i] * l[i] + PrinicipalDamageIncrementVector[1] * m[i] * m[i] + PrinicipalDamageIncrementVector[2] * n[i] * n[i];
    }
    DamageIncrementTensor(0,1) = PrinicipalDamageIncrementVector[0] * l[0] * l[1] + PrinicipalDamageIncrementVector[1] * m[0] * m[1] + PrinicipalDamageIncrementVector[2] * n[0] * n[1];
    DamageIncrementTensor(1,0) = PrinicipalDamageIncrementVector[0] * l[0] * l[1] + PrinicipalDamageIncrementVector[1] * m[0] * m[1] + PrinicipalDamageIncrementVector[2] * n[0] * n[1];
    DamageIncrementTensor(0,2) = PrinicipalDamageIncrementVector[0] * l[0] * l[2] + PrinicipalDamageIncrementVector[1] * m[0] * m[2] + PrinicipalDamageIncrementVector[2] * n[0] * n[2];
    DamageIncrementTensor(2,0) = PrinicipalDamageIncrementVector[0] * l[0] * l[2] + PrinicipalDamageIncrementVector[1] * m[0] * m[2] + PrinicipalDamageIncrementVector[2] * n[0] * n[2];
    DamageIncrementTensor(1,2) = PrinicipalDamageIncrementVector[0] * l[1] * l[2] + PrinicipalDamageIncrementVector[1] * m[1] * m[2] + PrinicipalDamageIncrementVector[2] * n[1] * n[2];
    DamageIncrementTensor(2,1) = PrinicipalDamageIncrementVector[0] * l[1] * l[2] + PrinicipalDamageIncrementVector[1] * m[1] * m[2] + PrinicipalDamageIncrementVector[2] * n[1] * n[2];
    KRATOS_CATCH("")
}

//************************************************************************************
//***********************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::GetTotalPrincipalDamageVector(
    BoundedVectorType& PrincipaldamageVector,
    const BoundedMatrixType& TotalDamageTensor
    )
{
    KRATOS_TRY
    BoundedMatrixType EigenVectors;
    BoundedMatrixType EigenValues = ZeroMatrix(3,3);
    MathUtils<double>::GaussSeidelEigenSystem(TotalDamageTensor, EigenVectors, EigenValues);
    PrincipaldamageVector[0] = EigenValues(0,0);
    PrincipaldamageVector[1] = EigenValues(1,1);
    PrincipaldamageVector[2] = EigenValues(2,2);
    KRATOS_CATCH("")
}

//************************************************************************************
//***********************************************************************************

Vector& ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::CalculateValue(
    ConstitutiveLaw::Parameters& rParametersValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    KRATOS_TRY
    //Explicitly having STRAIN and INITAL_STRAIN_VECTOR calculated in base class
    if (rThisVariable == EQUIVALENT_STRAINS_TC) {
        InternalDamageVariables rInternalDamageVariables;
        this->CalculateStressResponse( rParametersValues, rInternalDamageVariables);
        rValue = rInternalDamageVariables.EquivalentStrainsTC;
    }else{
        ElasticIsotropic3D::CalculateValue(rParametersValues, rThisVariable, rValue);
    }
    return(rValue);
    KRATOS_CATCH("")
}

//************************************************************************************
//***********************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::GetDamageEffectTensor(
    BoundedMatrixVoigtType& DamageEffectTensor,
    const BoundedVectorType& DamageVector
)
{
    KRATOS_TRY
    const double D1 =  DamageVector[0];
    const double D2 =  DamageVector[1];
    const double D3 =  DamageVector[2];
    DamageEffectTensor(0,0) = pow((1-D1),(-1));
    DamageEffectTensor(1,1) = pow((1-D2),(-1));
    DamageEffectTensor(2,2) = pow((1-D3),(-1));
    DamageEffectTensor(3,3) = pow((1-((D2+D3)*0.5)),(-1));
    DamageEffectTensor(4,4) = pow((1-((D3+D1)*0.5)),(-1));
    DamageEffectTensor(5,5) = pow((1-((D1+D2)*0.5)),(-1));
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::GetTransformedDamageEffectTensor(
    BoundedMatrixVoigtType& TransformedDamageEffectTensor,
    const BoundedMatrixVoigtType& DamageEffectTensor,
    const Vector& StrainVector
)
{
    KRATOS_TRY

    BoundedVectorVoigtType W = ZeroVector(6);
    BoundedVectorType l = ZeroVector(3);
    BoundedVectorType m = ZeroVector(3);
    BoundedVectorType n = ZeroVector(3);
    BoundedMatrixType MatrixForm = ZeroMatrix(3,3);
    BoundedMatrixType EigenVectors;
    BoundedMatrixType EigenValues = ZeroMatrix(3,3);
    for(SizeType i = 0; i < VoigtSize; ++i){
        W[i] = 1.0/DamageEffectTensor(i,i);
    }
    VectorToTensor(MatrixForm, StrainVector, STRAIN);
    MathUtils<double>::GaussSeidelEigenSystem(MatrixForm, EigenVectors, EigenValues);
    for(SizeType i = 0; i < Dimension; ++i){
        l[i] = EigenVectors(i,0);
        m[i] = EigenVectors(i,1);
        n[i] = EigenVectors(i,2);
    }
    for(SizeType i = 0; i < Dimension; ++i) {
        for(SizeType j = 0; j < Dimension; ++j) {
            TransformedDamageEffectTensor(i,j)  = ( (l[i] * l[i]) * (l[j] * l[j]) ) / W[0];
            TransformedDamageEffectTensor(i,j) += ( (m[i] * m[i]) * (m[j] * m[j]) ) / W[1] ;
            TransformedDamageEffectTensor(i,j) += ( (n[i] * n[i]) * (n[j] * n[j]) ) / W[2] ;
            TransformedDamageEffectTensor(i,j) += ( 2 * m[i] * n[i] * m[j] * n[j] ) / W[3] ;
            TransformedDamageEffectTensor(i,j) += ( 2 * n[i] * l[i] * n[j] * l[j] ) / W[4] ;
            TransformedDamageEffectTensor(i,j) += ( 2 * l[i] * m[i] * l[j] * m[j] ) / W[5] ;
        }
    }

    for(SizeType i = 3; i < 5; ++i ){
        TransformedDamageEffectTensor(i,i)  = ( 2.0 * l[2] * l[2] * l[4-i] * l[4-i] ) / W[0] ;
        TransformedDamageEffectTensor(i,i) += ( 2.0 * m[2] * m[2] * m[4-i] * m[4-i] ) / W[1] ;
        TransformedDamageEffectTensor(i,i) += ( 2.0 * n[2] * n[2] * n[4-i] * n[4-i] ) / W[2] ;
        TransformedDamageEffectTensor(i,i) += std::pow(( m[2] * n[4-i] + n[2] * m[4-i] ), 2.0) / W[3] ;
        TransformedDamageEffectTensor(i,i) += std::pow(( n[2] * l[4-i] + l[2] * n[4-i] ), 2.0) / W[4] ;
        TransformedDamageEffectTensor(i,i) += std::pow(( l[2] * m[4-i] + m[2] * l[4-i] ), 2.0) / W[5] ;
    }

    TransformedDamageEffectTensor(5,5)  =  ( 2.0 * l[0] * l[0] * l[1] * l[1] ) / W[0] ;
    TransformedDamageEffectTensor(5,5) +=  ( 2.0 * m[0] * m[0] * m[1] * m[1] ) / W[1] ;
    TransformedDamageEffectTensor(5,5) +=  ( 2.0 * n[0] * n[0] * n[1] * n[1] ) / W[2] ;
    TransformedDamageEffectTensor(5,5) +=  std::pow(( m[0] * n[1] + n[0] * m[1] ), 2.0) / W[3] ;
    TransformedDamageEffectTensor(5,5) +=  std::pow(( n[0] * l[1] + l[0] * n[1] ), 2.0) / W[4] ;
    TransformedDamageEffectTensor(5,5) +=  std::pow(( l[0] * m[1] + m[0] * l[1] ), 2.0) / W[5] ;

    for(SizeType i = 3; i < 5; ++i) {
        for(SizeType j = 0; j < Dimension; ++j) {
            TransformedDamageEffectTensor(i,j)  = ( l[2] * l[4-i] * l[j] * l[j] ) / W[0] ;
            TransformedDamageEffectTensor(i,j) += ( m[2] * m[4-i] * m[j] * m[j] ) / W[1] ;
            TransformedDamageEffectTensor(i,j) += ( n[2] * n[4-i] * n[j] * n[j] ) / W[2] ;
            TransformedDamageEffectTensor(i,j) += ( m[j] * n[j] * ( m[2] * n[4-i] +  n[2] * m[4-i] ) ) / W[3] ;
            TransformedDamageEffectTensor(i,j) += ( n[j] * l[j] * ( n[2] * l[4-i] +  l[2] * n[4-i] ) ) / W[4] ;
            TransformedDamageEffectTensor(i,j) += ( l[j] * m[j] * ( l[2] * m[4-i] +  m[2] * l[4-i] ) ) / W[5] ;
        }
    }

    TransformedDamageEffectTensor(4,3)  = ( 2 * l[0] * l[1] * l[2] * l[2] ) / W[0] ;
    TransformedDamageEffectTensor(4,3) += ( 2 * m[0] * m[1] * m[2] * m[2] ) / W[1] ;
    TransformedDamageEffectTensor(4,3) += ( 2 * n[0] * n[1] * n[2] * n[2] ) / W[2] ;
    TransformedDamageEffectTensor(4,3) += ( ( m[1] * n[2] + m[2] * n[1] ) * ( m[0] * n[2] + m[2] * n[0]) ) / W[3] ;
    TransformedDamageEffectTensor(4,3) += ( ( n[1] * l[2] + n[2] * l[1] ) * ( n[2] * l[0] + l[2] * n[0]) ) / W[4] ;
    TransformedDamageEffectTensor(4,3) += ( ( l[1] * m[2] + l[2] * m[1] ) * ( l[2] * m[0] + m[2] * l[0]) ) / W[5] ;

    for(SizeType j = 0; j < Dimension; ++j){
       TransformedDamageEffectTensor(5,j) =  ( l[0] * l[1] * l[j] * l[j] ) / W[0] ;
       TransformedDamageEffectTensor(5,j) += ( m[0] * m[1] * m[j] * m[j] ) / W[1] ;
       TransformedDamageEffectTensor(5,j) += ( n[0] * n[1] * n[j] * n[j] ) / W[2] ;
       TransformedDamageEffectTensor(5,j) += ( m[j] * n[j] * ( m[0] * n[1] + m[1] * n[0] ) ) / W[3] ;
       TransformedDamageEffectTensor(5,j) += ( n[j] * l[j] * ( n[0] * l[1] + l[0] * n[1] ) ) / W[4] ;
       TransformedDamageEffectTensor(5,j) += ( l[j] * m[j] * ( l[0] * m[1] + m[0] * l[1] ) ) / W[5] ;
    }

    for(SizeType j = 3; j < 5; ++j){
       TransformedDamageEffectTensor(5,j) =  ( 2 * l[0] * l[1] * l[2] * l[4-j] ) / W[0] ;
       TransformedDamageEffectTensor(5,j) += ( 2 * m[0] * m[1] * m[2] * m[4-j] ) / W[1] ;
       TransformedDamageEffectTensor(5,j) += ( 2 * n[0] * n[1] * n[2] * n[4-j] ) / W[2] ;
       TransformedDamageEffectTensor(5,j) += ( ( m[2] * n[4-j] + m[4-j] * n[2] ) * ( m[0] * n[1] + n[1] * m[0] ) ) / W[3] ;
       TransformedDamageEffectTensor(5,j) += ( ( n[2] * l[4-j] + n[4-j] * l[2])  * ( n[0] * l[1] + l[0] * n[1] ) ) / W[4] ;
       TransformedDamageEffectTensor(5,j) += ( ( l[2] * m[4-j] + l[4-j] * m[2] ) * ( l[0] * m[1] + m[0] * l[1] ) ) / W[5] ;
    }

    for(SizeType i = 3; i < VoigtSize; ++i){
        for(SizeType j = 0; i < Dimension; ++j){
            TransformedDamageEffectTensor(j,i) = 2 * TransformedDamageEffectTensor(i,j);
        }
    }

    for(SizeType i = 4; i < VoigtSize; ++i){
        for(SizeType j = 3; j < 5; ++j){
            TransformedDamageEffectTensor(j,i) = TransformedDamageEffectTensor(i,j);
        }
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::VectorToTensor(
    BoundedMatrixType& TensorForm,
    const Vector& VectorForm,
    const Variable<Vector>& rThisVariable
)
{
    KRATOS_TRY
    if(rThisVariable == STRESSES){
        TensorForm(0,1)= TensorForm(1,0)= VectorForm[5];
        TensorForm(0,2)= TensorForm(2,0)= VectorForm[4];
        TensorForm(2,1)= TensorForm(1,2)= VectorForm[3];
    }else if(rThisVariable == STRAIN){
        TensorForm(0,1)= TensorForm(1,0)= 0.5 * VectorForm[5];
        TensorForm(0,2)= TensorForm(2,0)= 0.5 * VectorForm[4];
        TensorForm(2,1)= TensorForm(1,2)= 0.5 * VectorForm[3];
    }
    for(SizeType i = 0; i < Dimension; ++i){
        TensorForm(i,i) = VectorForm[i];
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::ScaleNonlocalEquivalentStrain(
    BoundedVectorType& Principal_Nonlocal_Strains,
    ConstitutiveLaw::Parameters& rParametersValues,
    const BoundedVectorType& Principal_Strains,
    const Vector& Nonlocal_Equivalent_Strains,
    const Vector& Local_Equivalent_Strains
)
{
    KRATOS_TRY

    const Vector& r_strain_vector = rParametersValues.GetStrainVector();
    BoundedVectorType PrincipalDeviatoricStrains = ZeroVector(3);
    GetPrincipalDeviatoricStrains(PrincipalDeviatoricStrains, r_strain_vector, Principal_Strains);
    for (SizeType i = 0; i < Dimension; ++i){
        if (Principal_Strains[i]>0.0){
            Principal_Nonlocal_Strains[i] = Principal_Strains[i] * Nonlocal_Equivalent_Strains[0]/(Local_Equivalent_Strains[0]+eps);
        }else{
            Principal_Nonlocal_Strains[i] = fabs(PrincipalDeviatoricStrains[i]) * Nonlocal_Equivalent_Strains[1]/(Local_Equivalent_Strains[1]+eps);
        }
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::Calculate_tangent_HuNL(
    BoundedVectorVoigtType& H_uNL1,
    BoundedVectorVoigtType& H_uNL2,
    ConstitutiveLaw::Parameters& rParametersValues,
    const ModelParameters& rModelParameters,
    const Vector& Damage_Vector,
    const BoundedVectorType& Principal_Strains
)
{
    KRATOS_TRY
    const Vector& r_strain_vector = rParametersValues.GetStrainVector();
    const Properties& r_material_properties = rParametersValues.GetMaterialProperties();
    const double E   = r_material_properties[YOUNG_MODULUS];
    const double nu  = r_material_properties[POISSON_RATIO];
    const double E_factor = E/((1 + nu ) * (1- 2 * nu));
    BoundedVectorType Kappa = rModelParameters.kappa;
    BoundedVectorType beta1 = rModelParameters.Beta1;
    BoundedVectorType beta2 = rModelParameters.Beta2;
    BoundedVectorType k0 = rModelParameters.kappa0;
    BoundedVectorType dDdkappa = ZeroVector(3);
    BoundedVectorType dkappadNLEpseq_tension = ZeroVector(3);
    BoundedVectorType dkappadNLEpseq_compression = ZeroVector(3);
    BoundedVectorType dDdNLEpseq_tension = ZeroVector(3);
    BoundedVectorType dDdNLEpseq_compression = ZeroVector(3);
    BoundedMatrix6x3Type dHdD = ZeroMatrix(6,3);
    BoundedMatrixVoigtType dHdEpseq_tension = ZeroMatrix(6,6);
    BoundedMatrixVoigtType dHdEpseq_compression = ZeroMatrix(6,6);
    Vector local_equivalent_strains = ZeroVector(2);
    CalculateLocalEquivalentStrains(local_equivalent_strains, rParametersValues, Principal_Strains);
    double Local_Equivalent_Strain_tension = local_equivalent_strains[0];
    double Local_Equivalent_Strain_compression = local_equivalent_strains[1];
    for(SizeType i =0; i < Dimension; ++i){
        dDdkappa[i]   = (1- Damage_Vector[i]) * (beta1[i]/Kappa[i] + beta2[i]/k0[i]);
        dkappadNLEpseq_tension[i] = MacaulayBrackets(Principal_Strains[i]) /(Local_Equivalent_Strain_tension+eps);
        dkappadNLEpseq_compression[i] = MacaulayBrackets(Principal_Strains[i]) /(Local_Equivalent_Strain_compression+eps);
        dDdNLEpseq_tension[i] = dDdkappa[i] * dkappadNLEpseq_tension[i];
        dDdNLEpseq_compression[i] = dDdkappa[i] * dkappadNLEpseq_compression[i];
        dHdD(i,i) = E_factor * (-2 * (1-Damage_Vector[i]) * (1-nu));
    }
    //KRATOS_WATCH(dDdkappa)
    //KRATOS_WATCH(dkappadNLEpseq_tension)
    //KRATOS_WATCH(dkappadNLEpseq_compression)
    //KRATOS_WATCH(dDdNLEpseq_tension)
    //KRATOS_WATCH(dDdNLEpseq_tension)
    dHdD(3,1) = E_factor * 0.5 * (1-2*nu) * ((0.5 * (Damage_Vector[1]+Damage_Vector[2]))-1) ;
    dHdD(3,2) = E_factor * 0.5 * (1-2*nu) * ((0.5 * (Damage_Vector[1]+Damage_Vector[2]))-1) ;
    dHdD(4,0) = E_factor * 0.5 * (1-2*nu) * ((0.5 * (Damage_Vector[0]+Damage_Vector[2]))-1) ;
    dHdD(4,2) = E_factor * 0.5 * (1-2*nu) * ((0.5 * (Damage_Vector[0]+Damage_Vector[2]))-1) ;
    dHdD(5,0) = E_factor * 0.5 * (1-2*nu) * ((0.5 * (Damage_Vector[0]+Damage_Vector[1]))-1) ;
    dHdD(5,1) = E_factor * 0.5 * (1-2*nu) * ((0.5 * (Damage_Vector[0]+Damage_Vector[1]))-1) ;
    for(SizeType i =0; i < VoigtSize; ++i){
        BoundedVectorVoigtType dHdEpseq_tensoin_diagonal = prod(dHdD, dDdNLEpseq_tension);
        BoundedVectorVoigtType dHdEpseq_compression_diagonal = prod(dHdD, dDdNLEpseq_compression);
        dHdEpseq_tension(i,i) = dHdEpseq_tensoin_diagonal[i];
        dHdEpseq_compression(i,i) = dHdEpseq_compression_diagonal[i];
    }
    H_uNL1 = prod(dHdEpseq_tension, r_strain_vector);
    H_uNL2 = prod(dHdEpseq_compression, r_strain_vector);
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::Calculate_tangent_HNLu(
    BoundedVectorVoigtType& H_NL1u,
    BoundedVectorVoigtType& H_NL2u,
    ConstitutiveLaw::Parameters& rParametersValues,
    const BoundedVectorType& Principal_Strains
)
{
    KRATOS_TRY
    BoundedMatrix3x6Type derivatives_of_eigen_values;
    const Vector& r_strain_vector = rParametersValues.GetStrainVector();
    Vector local_equivalent_strains = ZeroVector(2);
    CalculateLocalEquivalentStrains(local_equivalent_strains, rParametersValues, Principal_Strains);
    double Local_Equivalent_Strain_tension = local_equivalent_strains[0];
    double Local_Equivalent_Strain_compression = local_equivalent_strains[1];
    CalculateDerivativesofEigenvalues(derivatives_of_eigen_values, Principal_Strains, r_strain_vector, STRAIN);
    for(SizeType i = 0; i < Dimension; ++i){
        if (Principal_Strains[i]>0.0){
            for (SizeType j = 0; j < VoigtSize; ++j){
                H_NL1u[j]   +=  Principal_Strains[i] * derivatives_of_eigen_values(i,j);
                H_NL2u[j]   +=  Principal_Strains[i] * derivatives_of_eigen_values(i,j);
            }
        }
    }
    H_NL1u /= (Local_Equivalent_Strain_tension +eps) ;
    H_NL1u /= (Local_Equivalent_Strain_compression +eps) ;
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::CalculateDerivativesofEigenvalues(
     BoundedMatrix3x6Type &DerivativesofEigenvalues,
     const BoundedVectorType &EigenvaluesVector,
     const BoundedVectorVoigtType &Voigtform,
     const Variable<Vector>& rThisVariable
    )
{
    KRATOS_TRY
    BoundedMatrixType Matrixform;
    VectorToTensor(Matrixform, Voigtform, rThisVariable);
    BoundedMatrixType DerivtivesMatrix;
    for(SizeType i = 0; i < Dimension; ++i){
        for(SizeType j = 0; j < Dimension; ++j){
            if(i != j && Matrixform(i,j)< eps){
                Matrixform(i,j) = eps;
            }
        }
    }
    for(SizeType i = 0; i < Dimension; ++i){
        BoundedMatrixType AminusLambdaMatrix = Matrixform - EigenvaluesVector[i] * IdentityMatrix(Dimension, Dimension);
        BoundedMatrixType cofactor_matrix = MathUtils<double>::CofactorMatrix(AminusLambdaMatrix);
        const double trace = cofactor_matrix(0,0) + cofactor_matrix(1,1) + cofactor_matrix(2,2);
        DerivtivesMatrix= (1/trace) * cofactor_matrix;
        DerivativesofEigenvalues(i,0) = DerivtivesMatrix(0,0);
        DerivativesofEigenvalues(i,1) = DerivtivesMatrix(1,1);
        DerivativesofEigenvalues(i,2) = DerivtivesMatrix(2,2);
        DerivativesofEigenvalues(i,3) = DerivtivesMatrix(1,2);
        DerivativesofEigenvalues(i,4) = DerivtivesMatrix(0,2);
        DerivativesofEigenvalues(i,5) = DerivtivesMatrix(0,1);
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::AssembleConstitutiveMatrix(
    Matrix& ConstitutiveMatrix,
    const Matrix& H_uu,
    const Vector& H_NL1u,
    const Vector& H_NL2u,
    const Vector& H_uNL1,
    const Vector& H_uNL2,
    const double& H_NLNL
    )
{
    KRATOS_TRY
    const SizeType num_rows_C = H_uu.size1();
    const SizeType num_cols_C = H_uu.size2();
    const SizeType size_CuNL  = H_uNL1.size();
    const SizeType size_CNLu  = H_NL1u.size();

    for (SizeType i = 0; i < num_rows_C ; ++i) {
        for (SizeType j = 0; j < num_cols_C ; ++j) {
            ConstitutiveMatrix(i, j) = H_uu(i,j);
        }
    }
    for (SizeType i = 0; i < size_CuNL; ++i) {
        ConstitutiveMatrix(i, num_cols_C) = H_uNL1[i] ;
        ConstitutiveMatrix(i, num_cols_C) = H_uNL2[i] ;
    }
    for (SizeType i = 0; i < size_CNLu; ++i) {
        ConstitutiveMatrix(num_rows_C, i) = H_NL1u[i] ;
        ConstitutiveMatrix(num_rows_C +1, i) = H_NL2u[i] ;
    }
    ConstitutiveMatrix(num_rows_C, num_cols_C) = H_NLNL;
    ConstitutiveMatrix(num_rows_C+1, num_cols_C) = H_NLNL;
    ConstitutiveMatrix(num_rows_C, num_cols_C+1) = H_NLNL;
    ConstitutiveMatrix(num_rows_C+1, num_cols_C+1) = H_NLNL;
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainSize = this->GetStrainSize();
    rFeatures.mSpaceDimension = this->WorkingSpaceDimension();
}

//************************************************************************************
//************************************************************************************



void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("DamageMatrix", mInternalDamageVariables.DamageMatrix);
    rSerializer.save("EquivalentStrainsTC", mInternalDamageVariables.EquivalentStrainsTC);
    rSerializer.save("kappa", mInternalDamageVariables.kappa);
    rSerializer.save("StrainHistoryParameter", mInternalDamageVariables.StrainHistoryParameter);
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.load("DamageMatrix", mInternalDamageVariables.DamageMatrix);
    rSerializer.load("EquivalentStrainsTC", mInternalDamageVariables.EquivalentStrainsTC);
    rSerializer.load("kappa", mInternalDamageVariables.kappa);
    rSerializer.load("StrainHistoryParameter", mInternalDamageVariables.StrainHistoryParameter);
}



} /* namespace Kratos.*/
