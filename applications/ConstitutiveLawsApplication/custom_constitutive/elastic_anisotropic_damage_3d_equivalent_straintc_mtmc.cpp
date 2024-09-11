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
//  Main authors:    athira vadakkekkara
//  Collaborators:
//

// Project includes
#include "elastic_anisotropic_damage_3d_equivalent_straintc_mtmc.h"
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

ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::ElasticAnisotropicDamage3DEquivalentStrainTCMTMC()
    : ElasticIsotropic3D()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::ElasticAnisotropicDamage3DEquivalentStrainTCMTMC(const ElasticAnisotropicDamage3DEquivalentStrainTCMTMC &rOther)
    : ElasticIsotropic3D(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::Clone() const
{
    return Kratos::make_shared<ElasticAnisotropicDamage3DEquivalentStrainTCMTMC>(ElasticAnisotropicDamage3DEquivalentStrainTCMTMC(*this));
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::~ElasticAnisotropicDamage3DEquivalentStrainTCMTMC()
{
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

int ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::Check(
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

bool ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::Has(const Variable<Vector>& rThisVariable)
{

    if(rThisVariable == STRAIN){
        // explicitly returning "false", so the element calls CalculateValue(...)
        return false;
    }else if(rThisVariable == DAMAGE_VECTOR){
        // explicitly returning "false", so the element calls CalculateValue(...)
        return true;
    } else if(rThisVariable == EQUIVALENT_STRAINS_TC){
        return true;
    }

    return false;
}

//************************************************************************************
//************************************************************************************

Vector& ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValues
    )
{
    KRATOS_TRY
    if(rThisVariable == DAMAGE_VECTOR){
        rValues = mDamageVector;
    }else if(rThisVariable == EQUIVALENT_STRAINS_TC) {
        rValues = mEquivalentStrains;
    }
    return rValues;
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValues,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY
    if(rThisVariable == DAMAGE_VECTOR){
        mDamageVector = rValues;
    }else if(rThisVariable == EQUIVALENT_STRAINS_TC) {
        mEquivalentStrains = rValues;
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::FinalizeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rParametersValues)
{
    KRATOS_TRY
    Vector equivalent_strains;
    Vector damage_vector;
    this->CalculateStressResponse(rParametersValues, damage_vector, equivalent_strains);
    mEquivalentStrains = equivalent_strains;
    mDamageVector = damage_vector;
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rParametersValues)
{
    KRATOS_TRY
    Vector equivalent_strains;
    Vector damage_vector;
    CalculateStressResponse(rParametersValues, damage_vector, equivalent_strains);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::CalculateStressResponse(
    ConstitutiveLaw::Parameters& rParametersValues,
    Vector& rDamageVector,
    Vector& rEquivalentStrains)
{
    KRATOS_TRY
    const Properties& r_material_properties = rParametersValues.GetMaterialProperties();
    Flags& r_constitutive_law_options = rParametersValues.GetOptions();
    Vector& r_strain_vector = rParametersValues.GetStrainVector();
    CalculateValue(rParametersValues, STRAIN, r_strain_vector);

    // If we compute the tangent moduli or the stress
    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS )||
        r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ))
        {
        Vector& r_stress_vector       = rParametersValues.GetStressVector();
        const Vector& r_strain_vector = rParametersValues.GetStrainVector();
        Matrix& r_constitutive_matrix = rParametersValues.GetConstitutiveMatrix();
        Matrix r_elastic_tensor;
        r_elastic_tensor.resize(6, 6, false);
        CalculateElasticMatrix(r_elastic_tensor, rParametersValues);
        noalias(r_stress_vector)      = prod(r_elastic_tensor, r_strain_vector);

        BoundedMatrixVoigtType EffStiffnessMatrix = ZeroMatrix(6, 6);
        BoundedMatrixVoigtType damage_effect_tensor = ZeroMatrix(6, 6);
        BoundedMatrixVoigtType transformed_damage_effect_tensor = ZeroMatrix(6, 6);
        BoundedMatrixVoigtType Inv_M = ZeroMatrix(6, 6);
        BoundedVectorVoigtType H_uNL1 = ZeroVector(6);
        BoundedVectorVoigtType H_NL1u = ZeroVector(6);
        BoundedVectorVoigtType H_uNL2 = ZeroVector(6);
        BoundedVectorVoigtType H_NL2u = ZeroVector(6);
        BoundedVectorType Spr = ZeroVector(3);
        BoundedVectorType principal_nonlocal_equivalent_strains = ZeroVector(3);
        BoundedVectorType kappa = ZeroVector(3);
        BoundedVectorType k0 = ZeroVector(3);
        BoundedVectorType beta1 = ZeroVector(3);
        BoundedVectorType beta2 = ZeroVector(3);
        BoundedVectorType F = ZeroVector(3);
        BoundedVectorType principal_strains = ZeroVector(3);
        Vector damage_vector = ZeroVector(3);
        Vector local_equivalent_strains = ZeroVector(2);
        Vector nonlocal_equivalent_strains = ZeroVector(2);
        const auto& N = rParametersValues.GetShapeFunctionsValues();
        const double beta1t = r_material_properties[DAMAGE_MODEL_PARAMETER_BETA1_TENSION];
        const double beta2t = r_material_properties[DAMAGE_MODEL_PARAMETER_BETA2_TENSION];
        const double beta1c = r_material_properties[DAMAGE_MODEL_PARAMETER_BETA1_COMPRESSION];
        const double beta2c = r_material_properties[DAMAGE_MODEL_PARAMETER_BETA2_COMPRESSION];
        const double E = r_material_properties[YOUNG_MODULUS];
        const double nu = r_material_properties[POISSON_RATIO];
        const double ft = r_material_properties[YIELD_STRESS_TENSION];
        const double H_NLNL = 0.0;
        double det_M = 0.0;

        GetEigenValues(principal_strains, STRAIN, r_strain_vector);
        for (SizeType i = 0; i < N.size(); ++i) {
            nonlocal_equivalent_strains[0] += N[i] * rParametersValues.GetElementGeometry()[i].FastGetSolutionStepValue(NON_LOCAL_EQUIVALENT_STRAIN_TENSION, 0);
        }
        for (SizeType i = 0; i < N.size(); ++i) {
            nonlocal_equivalent_strains[1] += N[i] * rParametersValues.GetElementGeometry()[i].FastGetSolutionStepValue(NON_LOCAL_EQUIVALENT_STRAIN_COMPRESSION, 0);
        }
        CalculateLocalEquivalentStrains(local_equivalent_strains, rParametersValues, principal_strains);
        ScaleNonlocalEquivalentStrain(principal_nonlocal_equivalent_strains, rParametersValues, principal_strains, nonlocal_equivalent_strains, local_equivalent_strains);
        const double k0t = (r_material_properties.Has(DAMAGE_THRESHOLD_TENSION)==true) ? r_material_properties[DAMAGE_THRESHOLD_TENSION] : ft/E;
        const double k0c = (r_material_properties.Has(DAMAGE_THRESHOLD_COMPRESSION)==true) ? r_material_properties[DAMAGE_THRESHOLD_COMPRESSION] : (1 + nu) * (2.15*ft)/(E);

        for(SizeType i = 0; i < Dimension; ++i) {
            // k0[i] = (principal_strains[i] > eps) ? k0t : k0c;
            // beta1[i] = (principal_strains[i] > eps) ? beta1t : beta1c;
            // beta2[i] = (principal_strains[i] > eps) ? beta2t : beta2c;
            //changed to see if the notched specimen test is working
            k0[i] = k0t;
            beta1[i] =beta1t;
            beta2[i] =beta2t;
            kappa[i] = std::max(principal_nonlocal_equivalent_strains[i],k0[i]);
            F[i]     = principal_nonlocal_equivalent_strains[i]-kappa[i];
        }
        //Compute damage in principal directions
        for (SizeType i = 0; i < Dimension; ++i) {
            if (kappa[i]>= 0 && kappa[i]<=k0[i]){
                damage_vector[i]=0.0;
            }else if (kappa[i]> k0[i]){
                const double var1      = pow((k0[i]/kappa[i]),beta1[i]);
                const double var2      = exp(-beta2[i]*((kappa[i]-k0[i])/(k0[i])));
                damage_vector[i] = 1.0 - var1 * var2;
            }
            if(damage_vector[i] < 0.0){
                damage_vector[i] = 0.0;
            }
        }
        if(damage_vector[0] > 0.0 || damage_vector[1] > 0.0 || damage_vector[2] > 0.0){
            GetDamageEffectTensor(damage_effect_tensor, damage_vector);
            GetTransformedDamageEffectTensor(transformed_damage_effect_tensor, damage_effect_tensor, r_strain_vector);
            MathUtils<double>::InvertMatrix(transformed_damage_effect_tensor, Inv_M, det_M);
            const BoundedMatrixVoigtType temp = prod(r_elastic_tensor,trans(Inv_M));
            noalias(r_elastic_tensor) = prod(Inv_M, temp);
            noalias(r_stress_vector)  = prod(r_elastic_tensor, r_strain_vector);
            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                Calculate_tangent_Huu(r_elastic_tensor, rParametersValues);
                Calculate_tangent_HuNL(H_uNL1, H_uNL2, rParametersValues, kappa, beta1, beta2, k0, damage_vector, principal_strains);
                Calculate_tangent_HNLu(H_NL1u, H_NL2u, rParametersValues, principal_strains);
            }
        }
        AssembleConstitutiveMatrix(r_constitutive_matrix, r_elastic_tensor, H_NL1u, H_NL2u, H_uNL1, H_uNL2, H_NLNL);
        rDamageVector = damage_vector;
        rEquivalentStrains = local_equivalent_strains;
        }
    //KRATOS_WATCH(rParametersValues.GetProcessInfo()[TIME])
    //KRATOS_WATCH("----------------------------------------------------------------------------------------")
    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::GetEigenValues(
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

Vector& ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::CalculateValue(
    ConstitutiveLaw::Parameters& rParametersValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    KRATOS_TRY
    //Explicitly having STRAIN and INITAL_STRAIN_VECTOR calculated in base class
    if (rThisVariable == DAMAGE_VECTOR) {
        Vector damage_vector;
        Vector equivalent_strains;
        this->CalculateStressResponse( rParametersValues, damage_vector, equivalent_strains);
        rValue = damage_vector;
    }else if (rThisVariable == EQUIVALENT_STRAINS_TC) {
        Vector equivalent_strains;
        Vector damage_vector;
        this->CalculateStressResponse( rParametersValues, damage_vector, equivalent_strains);
        rValue = equivalent_strains;
    }else{
        ElasticIsotropic3D::CalculateValue(rParametersValues, rThisVariable, rValue);
    }
    return(rValue);
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::CalculateLocalEquivalentStrains(
    Vector& LocalEquivalentStrains,
    ConstitutiveLaw::Parameters& rParametersValues,
    const BoundedVectorType& PrincipalStrains
    )
{
    KRATOS_TRY
    Vector& r_strain_vector = rParametersValues.GetStrainVector();
    BoundedVectorType PrincipalDeviatoricStrains = ZeroVector(3);
    GetPrincipalDeviatoricStrains(PrincipalDeviatoricStrains, r_strain_vector, PrincipalStrains);
    LocalEquivalentStrains[0] = sqrt(pow(MacaulayBrackets(PrincipalStrains[0]),2) + pow(MacaulayBrackets(PrincipalStrains[1]),2) + pow(MacaulayBrackets(PrincipalStrains[2]),2));
    LocalEquivalentStrains[1] = sqrt(pow(MacaulayBrackets(PrincipalDeviatoricStrains[0]),2) + pow(MacaulayBrackets(PrincipalDeviatoricStrains[1]),2) + pow(MacaulayBrackets(PrincipalDeviatoricStrains[2]),2));
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::GetPrincipalDeviatoricStrains(
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
//************************************************************************************

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::GetDamageEffectTensor(
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
    // DamageEffectTensor(3,3) = pow((1-((D2+D3)*0.5)),(-1));
    // DamageEffectTensor(4,4) = pow((1-((D3+D1)*0.5)),(-1));
    // DamageEffectTensor(5,5) = pow((1-((D1+D2)*0.5)),(-1));
    DamageEffectTensor(3,3) = pow((sqrt ((1-D1)*(1-D2))),(-1));
    DamageEffectTensor(4,4) = pow((sqrt ((1-D2)*(1-D3))),(-1));
    DamageEffectTensor(5,5) = pow((sqrt ((1-D1)*(1-D3))),(-1));
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::GetTransformedDamageEffectTensor(
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
            TransformedDamageEffectTensor(i,j) += ( 2 * m[i] * n[i] * m[j] * n[j] ) / W[4] ;
            TransformedDamageEffectTensor(i,j) += ( 2 * n[i] * l[i] * n[j] * l[j] ) / W[5] ;
            TransformedDamageEffectTensor(i,j) += ( 2 * l[i] * m[i] * l[j] * m[j] ) / W[3] ;
        }
    }

    for(SizeType i = 3; i < 5; ++i ){
        TransformedDamageEffectTensor(i,i)  = ( 2.0 * l[2] * l[2] * l[4-i] * l[4-i] ) / W[0] ;
        TransformedDamageEffectTensor(i,i) += ( 2.0 * m[2] * m[2] * m[4-i] * m[4-i] ) / W[1] ;
        TransformedDamageEffectTensor(i,i) += ( 2.0 * n[2] * n[2] * n[4-i] * n[4-i] ) / W[2] ;
        TransformedDamageEffectTensor(i,i) += std::pow(( m[2] * n[4-i] + n[2] * m[4-i] ), 2.0) / W[4] ;
        TransformedDamageEffectTensor(i,i) += std::pow(( n[2] * l[4-i] + l[2] * n[4-i] ), 2.0) / W[5] ;
        TransformedDamageEffectTensor(i,i) += std::pow(( l[2] * m[4-i] + m[2] * l[4-i] ), 2.0) / W[3] ;
    }

    TransformedDamageEffectTensor(5,5)  =  ( 2.0 * l[0] * l[0] * l[1] * l[1] ) / W[0] ;
    TransformedDamageEffectTensor(5,5) +=  ( 2.0 * m[0] * m[0] * m[1] * m[1] ) / W[1] ;
    TransformedDamageEffectTensor(5,5) +=  ( 2.0 * n[0] * n[0] * n[1] * n[1] ) / W[2] ;
    TransformedDamageEffectTensor(5,5) +=  std::pow(( m[0] * n[1] + n[0] * m[1] ), 2.0) / W[4] ;
    TransformedDamageEffectTensor(5,5) +=  std::pow(( n[0] * l[1] + l[0] * n[1] ), 2.0) / W[5] ;
    TransformedDamageEffectTensor(5,5) +=  std::pow(( l[0] * m[1] + m[0] * l[1] ), 2.0) / W[3] ;

    for(SizeType i = 3; i < 5; ++i) {
        for(SizeType j = 0; j < Dimension; ++j) {
            TransformedDamageEffectTensor(i,j)  = ( l[2] * l[4-i] * l[j] * l[j] ) / W[0] ;
            TransformedDamageEffectTensor(i,j) += ( m[2] * m[4-i] * m[j] * m[j] ) / W[1] ;
            TransformedDamageEffectTensor(i,j) += ( n[2] * n[4-i] * n[j] * n[j] ) / W[2] ;
            TransformedDamageEffectTensor(i,j) += ( m[j] * n[j] * ( m[2] * n[4-i] +  n[2] * m[4-i] ) ) / W[4] ;
            TransformedDamageEffectTensor(i,j) += ( n[j] * l[j] * ( n[2] * l[4-i] +  l[2] * n[4-i] ) ) / W[5] ;
            TransformedDamageEffectTensor(i,j) += ( l[j] * m[j] * ( l[2] * m[4-i] +  m[2] * l[4-i] ) ) / W[3] ;
        }
    }

    TransformedDamageEffectTensor(4,3)  = ( 2 * l[0] * l[1] * l[2] * l[2] ) / W[0] ;
    TransformedDamageEffectTensor(4,3) += ( 2 * m[0] * m[1] * m[2] * m[2] ) / W[1] ;
    TransformedDamageEffectTensor(4,3) += ( 2 * n[0] * n[1] * n[2] * n[2] ) / W[2] ;
    TransformedDamageEffectTensor(4,3) += ( ( m[1] * n[2] + m[2] * n[1] ) * ( m[0] * n[2] + m[2] * n[0]) ) / W[4] ;
    TransformedDamageEffectTensor(4,3) += ( ( n[1] * l[2] + n[2] * l[1] ) * ( n[2] * l[0] + l[2] * n[0]) ) / W[5] ;
    TransformedDamageEffectTensor(4,3) += ( ( l[1] * m[2] + l[2] * m[1] ) * ( l[2] * m[0] + m[2] * l[0]) ) / W[3] ;

    for(SizeType j = 0; j < Dimension; ++j){
       TransformedDamageEffectTensor(5,j) =  ( l[0] * l[1] * l[j] * l[j] ) / W[0] ;
       TransformedDamageEffectTensor(5,j) += ( m[0] * m[1] * m[j] * m[j] ) / W[1] ;
       TransformedDamageEffectTensor(5,j) += ( n[0] * n[1] * n[j] * n[j] ) / W[2] ;
       TransformedDamageEffectTensor(5,j) += ( m[j] * n[j] * ( m[0] * n[1] + m[1] * n[0] ) ) / W[4] ;
       TransformedDamageEffectTensor(5,j) += ( n[j] * l[j] * ( n[0] * l[1] + l[0] * n[1] ) ) / W[5] ;
       TransformedDamageEffectTensor(5,j) += ( l[j] * m[j] * ( l[0] * m[1] + m[0] * l[1] ) ) / W[3] ;
    }

    for(SizeType j = 3; j < 5; ++j){
       TransformedDamageEffectTensor(5,j) =  ( 2 * l[0] * l[1] * l[2] * l[4-j] ) / W[0] ;
       TransformedDamageEffectTensor(5,j) += ( 2 * m[0] * m[1] * m[2] * m[4-j] ) / W[1] ;
       TransformedDamageEffectTensor(5,j) += ( 2 * n[0] * n[1] * n[2] * n[4-j] ) / W[2] ;
       TransformedDamageEffectTensor(5,j) += ( ( m[2] * n[4-j] + m[4-j] * n[2] ) * ( m[0] * n[1] + n[1] * m[0] ) ) / W[4] ;
       TransformedDamageEffectTensor(5,j) += ( ( n[2] * l[4-j] + n[4-j] * l[2])  * ( n[0] * l[1] + l[0] * n[1] ) ) / W[5] ;
       TransformedDamageEffectTensor(5,j) += ( ( l[2] * m[4-j] + l[4-j] * m[2] ) * ( l[0] * m[1] + m[0] * l[1] ) ) / W[3] ;
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

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::VectorToTensor(
    BoundedMatrixType& TensorForm,
    const Vector& VectorForm,
    const Variable<Vector>& rThisVariable
)
{
    KRATOS_TRY
    if(rThisVariable == STRESSES){
        TensorForm(0,1)= TensorForm(1,0)= VectorForm[3];
        TensorForm(0,2)= TensorForm(2,0)= VectorForm[5];
        TensorForm(2,1)= TensorForm(1,2)= VectorForm[4];
    }else if(rThisVariable == STRAIN){
        TensorForm(0,1)= TensorForm(1,0)= 0.5 * VectorForm[3];
        TensorForm(0,2)= TensorForm(2,0)= 0.5 * VectorForm[5];
        TensorForm(2,1)= TensorForm(1,2)= 0.5 * VectorForm[4];
    }
    for(SizeType i = 0; i < Dimension; ++i){
        TensorForm(i,i) = VectorForm[i];
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::ScaleNonlocalEquivalentStrain(
    BoundedVectorType& Principal_Nonlocal_Strains,
    ConstitutiveLaw::Parameters& rParametersValues,
    const BoundedVectorType& Principal_Strains,
    const Vector& Nonlocal_Equivalent_Strains,
    Vector& Local_Equivalent_Strains
)
{
    KRATOS_TRY

    const Vector& r_strain_vector = rParametersValues.GetStrainVector();
    BoundedVectorType PrincipalDeviatoricStrains = ZeroVector(3);
    GetPrincipalDeviatoricStrains(PrincipalDeviatoricStrains, r_strain_vector, Principal_Strains);
    if(Local_Equivalent_Strains[0] == 0){
        Local_Equivalent_Strains[0] = 2 * eps;
    }
    if(Local_Equivalent_Strains[1] == 0){
        Local_Equivalent_Strains[1] = 2 * eps;
    }
    for (SizeType i = 0; i < Dimension; ++i){
        if (Principal_Strains[i]>0.0){
            Principal_Nonlocal_Strains[i] = Principal_Strains[i] * Nonlocal_Equivalent_Strains[0]/(Local_Equivalent_Strains[0]);
        }else if(Principal_Strains[i]<0.0){
            Principal_Nonlocal_Strains[i] = fabs(PrincipalDeviatoricStrains[i]) * Nonlocal_Equivalent_Strains[1]/(Local_Equivalent_Strains[1]);
        }else{
            Principal_Nonlocal_Strains[i] = 0.0;
        }
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::Calculate_tangent_Huu(
    Matrix& r_elasticity_matrix,
    ConstitutiveLaw::Parameters& rParametersValues)
{
    // Converged values to be storaged
    const Vector unperturbed_strain_vector_gp = Vector(rParametersValues.GetStrainVector());
    const Vector unperturbed_stress_vector_gp = Vector(rParametersValues.GetStressVector());

    // The number of components
    const SizeType num_components = unperturbed_stress_vector_gp.size();

    // The tangent tensor_Huu
    r_elasticity_matrix.clear();
    Matrix auxiliar_tensor = ZeroMatrix(num_components,num_components);

    // Calculate the perturbation
    double perturbation = PerturbationThreshold;
    if (rParametersValues.GetMaterialProperties().Has(PERTURBATION_SIZE)) {
        perturbation = rParametersValues.GetMaterialProperties()[PERTURBATION_SIZE];
    } else {
        for (IndexType i_component = 0; i_component < num_components; ++i_component) {
            double component_perturbation;
            CalculatePerturbation(unperturbed_strain_vector_gp, i_component, component_perturbation);
            perturbation = std::max(component_perturbation, perturbation);
        }
        // We check that the perturbation has a threshold value of PerturbationThreshold
        if (perturbation < PerturbationThreshold) perturbation = PerturbationThreshold;
    }

    // Loop over components of the strain
    Vector& r_perturbed_strain = rParametersValues.GetStrainVector();
    Vector& r_perturbed_integrated_stress = rParametersValues.GetStressVector();

    for (IndexType i_component = 0; i_component < num_components; ++i_component) {
        // Apply the perturbation (positive)
        PerturbStrainVector(r_perturbed_strain, unperturbed_strain_vector_gp, perturbation, i_component);

        // We continue with the calculations
        IntegratePerturbedStrain(rParametersValues, this, StressMeasure_PK2);

        // Compute stress (plus)
        const Vector strain_plus = r_perturbed_strain;
        const Vector stress_plus = r_perturbed_integrated_stress;

        // Reset the values to the initial ones
        noalias(r_perturbed_strain) = unperturbed_strain_vector_gp;
        noalias(r_perturbed_integrated_stress) = unperturbed_stress_vector_gp;

        // Apply the perturbation (negative)
        PerturbStrainVector(r_perturbed_strain, unperturbed_strain_vector_gp, - perturbation, i_component);

        // We continue with the calculations
        IntegratePerturbedStrain(rParametersValues, this,StressMeasure_PK2);

        // Compute stress (minus)
        const Vector strain_minus = r_perturbed_strain;
        const Vector stress_minus = r_perturbed_integrated_stress;

        // Finally we compute the components
        CalculateComponentsToTangentTensorSecondOrder(auxiliar_tensor, strain_plus, strain_minus, stress_plus, stress_minus, i_component);

        // Reset the values to the initial ones
        noalias(r_perturbed_strain) = unperturbed_strain_vector_gp;
        noalias(r_perturbed_integrated_stress) = unperturbed_stress_vector_gp;
    }

    noalias(r_elasticity_matrix) = 0.5*(auxiliar_tensor + trans(auxiliar_tensor));
}
//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::CalculatePerturbation(
    const Vector& rStrainVector,
    const IndexType Component,
    double& rPerturbation)
{
    double perturbation_1, perturbation_2;
    if (std::abs(rStrainVector[Component]) > tolerance) {
        perturbation_1 = PerturbationCoefficient1 * rStrainVector[Component];
    } else {
        double min_strain_component;
        GetMinAbsValue(rStrainVector, min_strain_component);
        perturbation_1 = PerturbationCoefficient1 * min_strain_component;
    }
    double max_strain_component;
    GetMaxAbsValue(rStrainVector, max_strain_component);
    perturbation_2 = PerturbationCoefficient2 * max_strain_component;
    rPerturbation = std::max(perturbation_1, perturbation_2);
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::GetMinAbsValue(
    const Vector& rArrayValues,
    double& rMinValue)
{
    const SizeType dimension = rArrayValues.size();

    IndexType counter = 0;
    double aux = std::numeric_limits<double>::max();
    for (IndexType i = 0; i < dimension; ++i) {
        if (std::abs(rArrayValues[i]) < aux) {
            aux = std::abs(rArrayValues[i]);
            ++counter;
        }
    }

    rMinValue = aux;
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::GetMaxAbsValue(
    const Vector& rArrayValues,
    double& rMaxValue)
{
    const SizeType dimension = rArrayValues.size();

    IndexType counter = 0;
    double aux = 0.0;
    for (IndexType i = 0; i < dimension; ++i) {
        if (std::abs(rArrayValues[i]) > aux) {
            aux = std::abs(rArrayValues[i]);
            ++counter;
        }
    }

    rMaxValue = aux;
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::PerturbStrainVector(
    Vector& rPerturbedStrainVector,
    const Vector& rStrainVectorGP,
    const double Perturbation,
    const IndexType Component
    )
{
    noalias(rPerturbedStrainVector) = rStrainVectorGP;
    rPerturbedStrainVector[Component] += Perturbation;
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::IntegratePerturbedStrain(
    ConstitutiveLaw::Parameters& rValues,
    ConstitutiveLaw* pConstitutiveLaw,
    const ConstitutiveLaw::StressMeasure& rStressMeasure)
{
    Flags& cl_options = rValues.GetOptions();

    // In order to avoid recursivity...
    const bool flag_back_up_1 = cl_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_back_up_2 = cl_options.Is(ConstitutiveLaw::COMPUTE_STRESS);

    cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    pConstitutiveLaw->CalculateMaterialResponse(rValues, rStressMeasure);

    cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_back_up_1);
    cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_back_up_2);
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::CalculateComponentsToTangentTensorSecondOrder(
    Matrix& rTangentTensor,
    const Vector& rVectorStrainPlus,
    const Vector& rVectorStrainMinus,
    const Vector& rStressPlus,
    const Vector& rStressMinus,
    const IndexType Component
    )
{
    const double perturbation = (rVectorStrainPlus[Component] - rVectorStrainMinus[Component]);
    const SizeType voigt_size = rStressPlus.size();
    for (IndexType row = 0; row < voigt_size; ++row) {
        rTangentTensor(row, Component) = (rStressPlus[row] - rStressMinus[row]) / perturbation;
    }
}

//************************************************************************************
//************************************************************************************
void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::Calculate_tangent_HuNL(
    BoundedVectorVoigtType& H_uNL1,
    BoundedVectorVoigtType& H_uNL2,
    ConstitutiveLaw::Parameters& rParametersValues,
    const BoundedVectorType& Kappa,
    const BoundedVectorType& beta1,
    const BoundedVectorType& beta2,
    const BoundedVectorType& k0,
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
    BoundedMatrixType dDdkappa = ZeroMatrix(3,3);
    BoundedVectorType dkappadNLEpseq_tension = ZeroVector(3);
    BoundedVectorType dkappadNLEpseq_compression = ZeroVector(3);
    BoundedVectorType dDdNLEpseq_tension = ZeroVector(3);
    BoundedVectorType dDdNLEpseq_compression = ZeroVector(3);
    array_1d<BoundedMatrix<double, 6, 6>, 3> dHdD;
    for (SizeType i = 0; i < Dimension; ++i){
        dHdD[i] = ZeroMatrix(6,6);
    }
    BoundedMatrixVoigtType dHdEpseq_tension = ZeroMatrix(6,6);
    BoundedMatrixVoigtType dHdEpseq_compression = ZeroMatrix(6,6);
    Vector local_equivalent_strains = ZeroVector(2);
    CalculateLocalEquivalentStrains(local_equivalent_strains, rParametersValues, Principal_Strains);
    double Local_Equivalent_Strain_tension = local_equivalent_strains[0];
    double Local_Equivalent_Strain_compression = local_equivalent_strains[1];
    for(SizeType i =0; i < Dimension; ++i){
        dDdkappa(i,i)    = (1- Damage_Vector[i]) * (beta1[i]/Kappa[i] + beta2[i]/k0[i]);
        dkappadNLEpseq_tension[i] = MacaulayBrackets(Principal_Strains[i]) /(Local_Equivalent_Strain_tension+eps);
        dkappadNLEpseq_compression[i] = MacaulayBrackets(Principal_Strains[i]) /(Local_Equivalent_Strain_compression+eps);
        dDdNLEpseq_tension[i] = dDdkappa(i,i)  * dkappadNLEpseq_tension[i];
        dDdNLEpseq_compression[i] = dDdkappa(i,i)  * dkappadNLEpseq_compression[i];
        dHdD[i](i,i) = E_factor * (-2 * (1-Damage_Vector[i]) * (1-nu));
    }
    for(SizeType i =0; i < Dimension; ++i){
        for(SizeType j = 0; j < Dimension; ++j){
            for(SizeType k = 0; k < Dimension; ++k){
                if(j != k && (i==j ||i ==k)){
                    SizeType m = (i==j) ? k : j;
                    dHdD[i](j,k) = -E_factor * nu * (1 - Damage_Vector[m]);
                }
            }
        }
    }
    // dHdD[0](4,4) = E_factor * 0.5 * (1-2*nu) * ((0.5 * (Damage_Vector[0]+Damage_Vector[2]))-1) ;
    // dHdD[0](5,5) = E_factor * 0.5 * (1-2*nu) * ((0.5 * (Damage_Vector[0]+Damage_Vector[1]))-1) ;
    // dHdD[1](3,3) = E_factor * 0.5 * (1-2*nu) * ((0.5 * (Damage_Vector[1]+Damage_Vector[2]))-1) ;
    // dHdD[1](5,5) = E_factor * 0.5 * (1-2*nu) * ((0.5 * (Damage_Vector[1]+Damage_Vector[0]))-1) ;
    // dHdD[2](3,3) = E_factor * 0.5 * (1-2*nu) * ((0.5 * (Damage_Vector[2]+Damage_Vector[1]))-1) ;
    // dHdD[2](4,4) = E_factor * 0.5 * (1-2*nu) * ((0.5 * (Damage_Vector[2]+Damage_Vector[0]))-1) ;
    dHdD[0](3,3) = E_factor * 0.5 * (1-2*nu) * (Damage_Vector[1]-1) ;
    dHdD[0](5,5) = E_factor * 0.5 * (1-2*nu) * (Damage_Vector[2]-1) ;
    dHdD[1](3,3) = E_factor * 0.5 * (1-2*nu) * (Damage_Vector[0]-1) ;
    dHdD[1](4,4) = E_factor * 0.5 * (1-2*nu) * (Damage_Vector[2]-1) ;
    dHdD[2](4,4) = E_factor * 0.5 * (1-2*nu) * (Damage_Vector[1]-1) ;
    dHdD[2](5,5) = E_factor * 0.5 * (1-2*nu) * (Damage_Vector[0]-1) ;
    TensorProduct(dHdEpseq_tension, dHdEpseq_compression, dHdD, dDdNLEpseq_tension, dDdNLEpseq_compression);
    H_uNL1 = prod(dHdEpseq_tension, r_strain_vector);
    H_uNL2 = prod(dHdEpseq_compression, r_strain_vector);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::TensorProduct(
    BoundedMatrixVoigtType& dHdNL1,
    BoundedMatrixVoigtType& dHdNL2,
    const array_1d<BoundedMatrix<double, 6, 6>, 3>& dHdD,
    const BoundedVectorType& Vector1,
    const BoundedVectorType& Vector2
    )
{
    KRATOS_TRY
    for (SizeType i = 0; i < VoigtSize; ++i) {
        for (SizeType j = 0; j < VoigtSize; ++j) {
            for (SizeType k = 0; k < Dimension; ++k) {
                dHdNL1(i,j) += dHdD[k](i,j)* Vector1[k];
                dHdNL2(i,j) += dHdD[k](i,j)* Vector2[k];
            }
        }
    }
    KRATOS_CATCH("")
}
//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::Calculate_tangent_HNLu(
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
    BoundedVectorType PrincipalDeviatoricStrains = ZeroVector(3);
    BoundedVectorVoigtType dHydrostaticdEpsilon = ZeroVector(6);
    dHydrostaticdEpsilon[0] = 1/3;
    dHydrostaticdEpsilon[1] = 1/3;
    dHydrostaticdEpsilon[2] = 1/3;
    GetPrincipalDeviatoricStrains(PrincipalDeviatoricStrains, r_strain_vector, Principal_Strains);
    for(SizeType i = 0; i < Dimension; ++i){
        for (SizeType j = 0; j < VoigtSize; ++j){
            H_NL1u[j]   +=  MacaulayBrackets(Principal_Strains[i]) * derivatives_of_eigen_values(i,j);
            H_NL2u[j]   +=  MacaulayBrackets(PrincipalDeviatoricStrains[i]) * (derivatives_of_eigen_values(i,j) - dHydrostaticdEpsilon[j]);
        }
    }
    H_NL1u /= (Local_Equivalent_Strain_tension +eps) ;
    H_NL2u /= (Local_Equivalent_Strain_compression +eps) ;
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::CalculateDerivativesofEigenvalues(
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
    // for(SizeType i = 0; i < Dimension; ++i){
    //     for(SizeType j = 0; j < Dimension; ++j){
    //         if(i != j && Matrixform(i,j)< eps){
    //             Matrixform(i,j) = eps;
    //         }
    //     }
    // }
    // for(SizeType i = 0; i < Dimension; ++i){
    //     BoundedMatrixType AminusLambdaMatrix = Matrixform - EigenvaluesVector[i] * IdentityMatrix(Dimension, Dimension);
    //     BoundedMatrixType cofactor_matrix = MathUtils<double>::CofactorMatrix(AminusLambdaMatrix);
    //     const double trace = cofactor_matrix(0,0) + cofactor_matrix(1,1) + cofactor_matrix(2,2);
    //     DerivtivesMatrix= (1/trace) * cofactor_matrix;
    //     DerivativesofEigenvalues(i,0) = DerivtivesMatrix(0,0);
    //     DerivativesofEigenvalues(i,1) = DerivtivesMatrix(1,1);
    //     DerivativesofEigenvalues(i,2) = DerivtivesMatrix(2,2);
    //     DerivativesofEigenvalues(i,3) = DerivtivesMatrix(0,1);
    //     DerivativesofEigenvalues(i,4) = DerivtivesMatrix(1,2);
    //     DerivativesofEigenvalues(i,5) = DerivtivesMatrix(0,2);
    // }
    BoundedMatrixType EigenVectors;
    BoundedMatrixType EigenValues = ZeroMatrix(3,3);
    MathUtils<double>::GaussSeidelEigenSystem(Matrixform, EigenVectors, EigenValues);
    BoundedVectorType Eigen_vector_column = ZeroVector (3) ;
    BoundedMatrixType DerivativesMatrix = ZeroMatrix(3,3);
    for(SizeType i = 0; i < Dimension; ++i){
        Eigen_vector_column = column(EigenVectors,i);
        DerivativesMatrix = outer_prod(Eigen_vector_column, Eigen_vector_column);
        DerivativesofEigenvalues(i,0) = DerivativesMatrix(0,0);
        DerivativesofEigenvalues(i,1) = DerivativesMatrix(1,1);
        DerivativesofEigenvalues(i,2) = DerivativesMatrix(2,2);
        DerivativesofEigenvalues(i,3) = DerivativesMatrix(0,1);
        DerivativesofEigenvalues(i,4) = DerivativesMatrix(1,2);
        DerivativesofEigenvalues(i,5) = DerivativesMatrix(2,0);
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::AssembleConstitutiveMatrix(
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
        ConstitutiveMatrix(i, num_cols_C + 1) = H_uNL2[i] ;
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

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::GetLawFeatures(Features& rFeatures)
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



void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mEquivalentStrains", mEquivalentStrains);
    rSerializer.save("mDamageVector", mDamageVector);
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DEquivalentStrainTCMTMC::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.load("mEquivalentStrains", mEquivalentStrains);
    rSerializer.load("mDamageVector", mDamageVector);
}



} /* namespace Kratos.*/
