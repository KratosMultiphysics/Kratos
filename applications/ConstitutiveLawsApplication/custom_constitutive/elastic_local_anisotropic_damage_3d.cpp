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
#include "elastic_local_anisotropic_damage_3d.h"
#include "constitutive_laws_application_variables.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "includes/checks.h"
#include "utilities/math_utils.h"
#include "math.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ElasticAnisotropicDamage::ElasticAnisotropicDamage()
    : ElasticIsotropic3D()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

ElasticAnisotropicDamage::ElasticAnisotropicDamage(const ElasticAnisotropicDamage &rOther)
    : ElasticIsotropic3D(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer ElasticAnisotropicDamage::Clone() const
{
    return Kratos::make_shared<ElasticAnisotropicDamage>(ElasticAnisotropicDamage(*this));
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

ElasticAnisotropicDamage::~ElasticAnisotropicDamage()
{
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

int ElasticAnisotropicDamage::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    // const double tolerance = std::numeric_limits<double>::epsilon();

    KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(POISSON_RATIO));
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS_TENSION));
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS_COMPRESSION));
    return 0;
}

//************************************************************************************
//************************************************************************************

bool ElasticAnisotropicDamage::Has(const Variable<double>& rThisVariable)
{
    if(rThisVariable == STRAIN_ENERGY){
        // explicitly returning "false", so the element calls CalculateValue(...)
        return false;
    } else if(rThisVariable == DAMAGE_VARIABLE){
        return false;
    }

    return false;
}

//************************************************************************************
//************************************************************************************

bool ElasticAnisotropicDamage::Has(const Variable<Vector>& rThisVariable)
{
    if(rThisVariable == INTERNAL_VARIABLES){
        // explicitly returning "false", so the element calls CalculateValue(...)s
        return false;
    } else if(rThisVariable == STRAIN){
        // explicitly returning "false", so the element calls CalculateValue(...)
        return false;
    } else if(rThisVariable == DAMAGE_VECTOR){
        return true;
    }

    return false;
}

//************************************************************************************
//************************************************************************************

Vector& ElasticAnisotropicDamage::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValues
    )
{
    KRATOS_TRY
    if(rThisVariable == DAMAGE_VECTOR){
        rValues = mDamageVector;
    }
    return rValues;
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValues,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY
    if(rThisVariable == DAMAGE_VECTOR){
        mDamageVector = rValues;
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // const double yield_stress = rMaterialProperties[STRESS_LIMITS](0);
    // const double young_modulus = rMaterialProperties[YOUNG_MODULUS];
    // mStrainVariable = yield_stress / std::sqrt(young_modulus);
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::FinalizeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rParametersValues)
{
    KRATOS_TRY
    Vector damage_vector;
    this->CalculateStressResponse(rParametersValues, damage_vector);
    mDamageVector = damage_vector;

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rParametersValues)
{
    KRATOS_TRY
    Vector damage_vector;
    CalculateStressResponse(rParametersValues, damage_vector);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::CalculateStressResponse(
    ConstitutiveLaw::Parameters& rParametersValues,
    Vector& rDamageVector)
{
    KRATOS_TRY
    const Properties& r_material_properties = rParametersValues.GetMaterialProperties();
    Flags& r_constitutive_law_options = rParametersValues.GetOptions();
    Vector& r_strain_vector = rParametersValues.GetStrainVector();
    KRATOS_WATCH(r_strain_vector);
    CalculateValue(rParametersValues, STRAIN, r_strain_vector);

    // If we compute the tangent moduli or the stress
    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ||
        r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR )) 
        {
        Vector& r_stress_vector       = rParametersValues.GetStressVector();
        const Vector& r_strain_vector = rParametersValues.GetStrainVector();   
        Matrix& r_constitutive_matrix = rParametersValues.GetConstitutiveMatrix();
        CalculateElasticMatrix(r_constitutive_matrix, rParametersValues);
        noalias(r_stress_vector)      = prod(r_constitutive_matrix, r_strain_vector);
        KRATOS_WATCH(r_stress_vector);
        const double eps = 1e-8;
        const double beta1t = r_material_properties[DAMAGE_MODEL_PARAMETER_BETA1_TENSION];
        const double beta2t = r_material_properties[DAMAGE_MODEL_PARAMETER_BETA2_TENSION];
        const double beta1c = r_material_properties[DAMAGE_MODEL_PARAMETER_BETA1_COMPRESSION];
        const double beta2c = r_material_properties[DAMAGE_MODEL_PARAMETER_BETA2_COMPRESSION];
        const double E   = r_material_properties[YOUNG_MODULUS];
        const double fck = r_material_properties[YIELD_STRESS_COMPRESSION];
        const double ft  = r_material_properties[YIELD_STRESS_TENSION];
        const double k0t0 = r_material_properties[DAMAGE_THRESHOLD_TENSION];
        const double k0c0 = r_material_properties[DAMAGE_THRESHOLD_COMPRESSION];
        double k0t, k0c ;
        double SprMax = 0.0;
        double Max_principal_strain = 0.0;
        BoundedVectorType damage_vector= ZeroVector(3);
        BoundedVectorType Spr = ZeroVector(3);
        BoundedVectorType kappa = ZeroVector(3);
        BoundedVectorType F     = ZeroVector(3);
        BoundedMatrixVoigtType EffStiffnessMatrix = ZeroMatrix(6, 6); 
        BoundedMatrix3x6Type dEprdE = ZeroMatrix(3,6);
        BoundedMatrixType dkdEpr = ZeroMatrix(3,3);
        BoundedMatrix6x3Type dHdk   = ZeroMatrix(6,3);

        GetEigenValues(Spr, SprMax, STRESSES, r_stress_vector);
        ManipulationOfZeroEntries(Spr, eps);         //manipulation of zero entries in stress
        const double H = (SprMax < eps) ? 0.0 : 1.0;
        double r=0.0; 
        GetStressWeightFactor(r,Spr);
        const double del_r = (r > 0 && r < 0.1) ? 1.0 : 0.0;
        if (k0t0==0 || k0c0==0){
            k0t = fck/(E*7.93103);//to be modified later
            k0c = (10./3.) * ft/E;
        }else{
            k0t = std::min(k0t0,fck/E);
            k0c = std::min(k0c0,(10./3.) * ft/E);
        }
        const double k0 = k0t  * H * (1-del_r) + (1.0- H + del_r) * k0c;
        const double beta1 = beta1t  * H * (1-del_r) + (1.-H + del_r) * beta1c;
        const double beta2 = beta2t  * H * (1-del_r) + (1.-H + del_r) * beta2c;
        BoundedVectorType principal_strains = ZeroVector(3);

        GetEigenValues(principal_strains, Max_principal_strain, STRAIN, r_strain_vector);  //calculate prinicpal strains
        for(SizeType i = 0; i < Dimension; ++i) {
            kappa[i]       = std::max(abs(principal_strains[i]),k0);
            F[i]           = abs(principal_strains[i])-kappa[i];
        }
        //Compute damage in principal directions
        for (SizeType i = 0; i < Dimension; ++i) {
            if (F[i] < 0) { 
                damage_vector[i]=0.0;
            }else{
                const double var1      = pow((k0/kappa[i]),beta1);
                const double var2      = exp(-beta2*((kappa[i]-k0)/(k0)));
                damage_vector[i] = 1.0 - var1 * var2;              
            }
            if(damage_vector[i] < 0.0){
                damage_vector[i] = 0.0;
            }
        } 
        CalculateParameters(EffStiffnessMatrix, dEprdE, dkdEpr, rParametersValues, damage_vector);
        CalculatePartialDerivatives(dHdk, r_material_properties, damage_vector, k0, beta1, beta2, kappa);
        const BoundedMatrix3x6Type b = prod(dkdEpr, dEprdE);
        noalias(r_constitutive_matrix) = EffStiffnessMatrix + prod(dHdk, b);
        noalias(r_stress_vector)       = prod(r_constitutive_matrix,r_strain_vector);
        KRATOS_WATCH(r_stress_vector);
        KRATOS_WATCH(damage_vector);
        rDamageVector = damage_vector;
        }
        KRATOS_WATCH(rParametersValues.GetProcessInfo()[TIME]);
        KRATOS_WATCH("-------------------------------------------");
    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::GetEigenValues(
    BoundedVectorType& Pri_Values,
    double& MaxValue,
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
    MaxValue    = std::max({Pri_Values[0],Pri_Values[1],Pri_Values[2]}, [](const double a, const double b){return std::abs(b) > std::abs(a);});
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::ComputedSprdS(
    BoundedMatrix3x6Type& dSprdS,
    const Vector& VectorForm,
    const BoundedVectorType& PrincipalVector,
    const Variable<Vector>& rThisVariable
    )
{   
    KRATOS_TRY
    BoundedMatrixType indx = ZeroMatrix(3, 3);
    BoundedVectorVoigtType stress2;
    BoundedVectorVoigtType dSprdS_entries;
    BoundedVectorVoigtType dI1dS = ZeroVector(6);
    dI1dS[0] = dI1dS[1] = dI1dS[2]= 1.0; 
    indx(0,1)= indx(1,0)= indx(2,2)=1.0;
    indx(0,2)= indx(1,1)= indx(2,0)=2.0;
    BoundedMatrixType MatrixForm = ZeroMatrix(3, 3);

    VectorToTensor(MatrixForm, VectorForm, rThisVariable);
    BoundedMatrixType stress(MatrixForm);
    stress = prod(stress, stress);
    stress2[0]= stress(0,0);
    stress2[1]= stress(1,1); 
    stress2[2]= stress(2,2);
    stress2[3]= stress(1,2);
    stress2[4]= stress(0,2); 
    stress2[5]= stress(0,1);
    for(SizeType i = 0; i < Dimension; ++i){  
        dSprdS_entries  = stress2 - ( PrincipalVector(indx(1,i)) + PrincipalVector(indx(2,i)) )*(VectorForm) + PrincipalVector(indx(1,i)) * PrincipalVector(indx(2,i)) * dI1dS;
        dSprdS_entries /=  (PrincipalVector(indx(0,i)) - PrincipalVector(indx(1,i)) ) * ( PrincipalVector(indx(0,i)) - PrincipalVector(indx(2,i)) );
        for(SizeType j = 0; j < VoigtSize; ++j){
            dSprdS(i,j)= dSprdS_entries[j];
        }
        for(SizeType k = Dimension; k < VoigtSize; ++k){
            dSprdS(i,k) *= 2.;
        }   
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::GetStressWeightFactor(
    double &w, 
    const BoundedVectorType &s_pr) const
{
    Vector N1(3);
    const double eps = 1.e-8;
    SizeType kk;
    for(kk=0; kk < Dimension; ++kk ){
        N1(kk) = 0.5 * ( abs(s_pr(kk)) + s_pr(kk) );
    }
    double N11 = N1(0) + N1(1) + N1(2);
    double D11 = eps + abs(s_pr(0)) + abs(s_pr(1)) + abs(s_pr(2));
    w = N11 / D11 ;
}

//************************************************************************************
//************************************************************************************

Vector& ElasticAnisotropicDamage::CalculateValue(
    ConstitutiveLaw::Parameters& rParametersValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValues
    )
{
    KRATOS_TRY
    if (rThisVariable == DAMAGE_VECTOR) {
        Vector damage_vector;
        this->CalculateStressResponse( rParametersValues, damage_vector);
        rValues = damage_vector;
    }else{
        ElasticIsotropic3D::CalculateValue(rParametersValues, rThisVariable, rValues);
    }
    return( rValues );
    KRATOS_CATCH("")
}

//************************************************************************************
//***********************************************************************************

void ElasticAnisotropicDamage::GetDamageEffectTensor(
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
    DamageEffectTensor(3,3) = 0.5*(pow((1-D2),(-1))+pow((1-D3),(-1)));
    DamageEffectTensor(4,4) = 0.5*(pow((1-D3),(-1))+pow((1-D1),(-1)));
    DamageEffectTensor(5,5) = 0.5*(pow((1-D1),(-1))+pow((1-D2),(-1)));

    KRATOS_CATCH("")
}
//************************************************************************************
//***********************************************************************************

void ElasticAnisotropicDamage::CalculateParameters(
    BoundedMatrixVoigtType& EffStiffnessMatrix,
    BoundedMatrix3x6Type& dEprdE,
    BoundedMatrixType& dkdEpr,
    ConstitutiveLaw::Parameters& rParametersValues, 
    const BoundedVectorType& DamageVector
)
{
    KRATOS_TRY

    Matrix& r_constitutive_matrix = rParametersValues.GetConstitutiveMatrix();
    const Vector& r_strain_vector = rParametersValues.GetStrainVector();   
    const double eps = 1e-8;     
    BoundedMatrixVoigtType M     = ZeroMatrix(6,6);
    BoundedMatrixVoigtType Inv_M = ZeroMatrix(6,6);
    BoundedVectorType principal_strains = ZeroVector(3);
    double det_M;
    double Max_principal_strain = 0.0;

    CalculateElasticMatrix(r_constitutive_matrix, rParametersValues);
    GetDamageEffectTensor(M, DamageVector);
    MathUtils<double>::InvertMatrix(M, Inv_M, det_M);
    const BoundedMatrixVoigtType a = prod(r_constitutive_matrix,trans(Inv_M));
    EffStiffnessMatrix = prod(Inv_M,a);    
    GetEigenValues(principal_strains, Max_principal_strain, STRAIN, r_strain_vector);
    ManipulationOfZeroEntries(principal_strains, eps);
    ComputedSprdS(dEprdE, r_strain_vector, principal_strains, STRAIN); 
    for(SizeType i = 0; i < Dimension; ++i){  
        if(DamageVector[i] > 0){
            dkdEpr(i,i) = 1; 
        } 
    } 

    KRATOS_CATCH("") 
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::CalculatePartialDerivatives(
    BoundedMatrix6x3Type& dHdk,
    const Properties& rMaterialProperties,
    const BoundedVectorType& DamageVector,
    const double Kappa0, 
    const double Beta1, 
    const double Beta2, 
    const BoundedVectorType& Kappa
)
{
    const double E   = rMaterialProperties[YOUNG_MODULUS];
    const double nu  = rMaterialProperties[POISSON_RATIO];
    const double E_factor = ((1 + nu ) * (1- nu))/E;
    double dDdk;
    for(SizeType i =0; i < VoigtSize; ++i){
        for(SizeType j =0; j < Dimension; ++j){
            if (i==j) {
                dDdk      = (1- DamageVector[j]) * (Beta1/Kappa[j] + Beta2/Kappa0);
                dHdk(j,j) = E_factor * (-2 * (1-DamageVector[j]) * (1-nu) * dDdk);
            } else if (i ==0 || i == 1 || i == 2){
                dHdk(i,j) = 0;
            }else{
                if(i-j == 3){
                    dHdk(i,j) = 0;
                }else{
                    dDdk      = (1- DamageVector[j]) * (Beta1/Kappa[j] + Beta2/Kappa0);
                    dHdk(i,j) = E_factor * 0.5 * (1-2*nu) * ((0.5 * (DamageVector[(j+1)%3]+DamageVector[(j+2)%3]))-1) * dDdk;
                }
            }
        }
    }
}

//************************************************************************************
//************************************************************************************
void ElasticAnisotropicDamage::ManipulationOfZeroEntries(
    BoundedVectorType& PrincipalVector,
    const double eps
)
{
    if( (abs(PrincipalVector(0)-PrincipalVector(1)) < 2.0*eps ) && ( abs(PrincipalVector(0)-PrincipalVector(2)) < 2.0*eps ) && ( abs(PrincipalVector(1)-PrincipalVector(2)) < 2.0*eps ) ){
            PrincipalVector(0) = PrincipalVector(0) + 1.1*eps;
            PrincipalVector(1) = PrincipalVector(1) + 0.75*eps;
            PrincipalVector(2) = PrincipalVector(2) + 0.5*eps;
    }else{
            if( abs(PrincipalVector(0)-PrincipalVector(1)) < eps ) PrincipalVector(1) = PrincipalVector(1) - eps;
            if( abs(PrincipalVector(1)-PrincipalVector(2)) < eps ) PrincipalVector(2) = PrincipalVector(2) - eps;
    }
}
//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::VectorToTensor(
    BoundedMatrixType& TensorForm,
    const Vector& VectorForm, 
    const Variable<Vector>& rThisVariable
)
{
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
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::GetLawFeatures(Features& rFeatures)
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



void ElasticAnisotropicDamage::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mDamageVector", mDamageVector);
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mDamageVector", mDamageVector);
}



} /* namespace Kratos.*/
