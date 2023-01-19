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
//  Main authors:    
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

        const double eps = 1e-8;
        const double beta1t = 0.85;
        const double beta2t = 0.18;
        const double beta1c = 0.0;
        const double beta2c = 0.095;
        const double E   = r_material_properties[YOUNG_MODULUS];
        const double fck = r_material_properties[YIELD_STRESS_COMPRESSION];
        const double ft  = r_material_properties[YIELD_STRESS_TENSION];
        const double k0t0 = r_material_properties[DAMAGE_THRESHOLD_TENSION];
        const double k0c0 = r_material_properties[DAMAGE_THRESHOLD_COMPRESSION];
        double SprMax, H, Max_principal_strain, k0t, k0c ;
        Vector dI1dS  = ZeroVector(6);
        Vector dJ2dS  = ZeroVector(6); 
        Vector damage_vector= ZeroVector(3);
        Vector Spr = ZeroVector(3);
        Vector dSmaxdSig = ZeroVector(6);;
        Matrix dJ2ddS = ZeroMatrix(6, 6); 
        Matrix dSprdS = ZeroMatrix(3,6);
        Matrix EffStiffnessMatrix = ZeroMatrix(6, 6); 
        Matrix dEprdE = ZeroMatrix(3,6);
        Matrix dkdEpr = ZeroMatrix(3,3);
        Matrix dHdk   = ZeroMatrix(6,3);

        GetDerivatives( r_stress_vector, dI1dS, dJ2ddS,dJ2dS);
        GetEigenValues(r_stress_vector, Spr, SprMax);
        //manipulation of zero entries in stress
        if( (fabs(Spr(0)-Spr(1)) < 2.0*eps ) && ( fabs(Spr(0)-Spr(2)) < 2.0*eps ) && ( fabs(Spr(1)-Spr(2)) < 2.0*eps ) ){
            Spr(0) = Spr(0) + 1.1*eps;
            Spr(1) = Spr(1) + 0.75*eps;
            Spr(2) = Spr(2) + 0.5*eps;
        }else{
            if( fabs(Spr(0)-Spr(1)) < eps ) Spr(1) = Spr(1) - eps;
            if( fabs(Spr(1)-Spr(2)) < eps ) Spr(2) = Spr(2) - eps;
        }
        if(SprMax < eps){
            H = 0.0;
        }else{
            H = 1.0;  
            ComputedSprdS(r_stress_vector, Spr, dSprdS); 
            for(int i=0; i<6; ++i){
                dSmaxdSig[i]= dSprdS(0,i); 
            }
        }

       //calculate prinicpal strains
        Vector principal_strains = ZeroVector(3);
        GetEigenValues(r_strain_vector, principal_strains, Max_principal_strain);
        double r=0.0; GetStressWeightFactor(r,Spr);
        KRATOS_WATCH(r);
        double del_r = 0.0;
        if (r>0 && r < 0.1){del_r = 1.0;}

        if (k0t0==0 || k0c0==0){
            k0t = fck/E;
            k0c = (10./3.) * ft/E;
        }else{
            k0t = std::min(k0t0,fck/E);
            k0c = std::min(k0c0,(10./3.) * ft/E);
        }
        double k0 = k0t  * H * (1-del_r) + (1.0- H + del_r) * k0c;
        double beta1 = beta1t  * H * (1-del_r) + (1.-H + del_r) * beta1c;
        double beta2 = beta2t  * H * (1-del_r) + (1.-H + del_r) * beta2c;

        Vector kappa = ZeroVector(3);
        Vector F     = ZeroVector(3);
        for(unsigned int i = 0; i < 3; i++) {
            kappa[i]       = std::max(fabs(principal_strains[i]),k0);
            F[i]           = fabs(principal_strains[i])-kappa[i];
        }
        //Compute damage in principal directions
        for (unsigned int i = 0; i < 3; i++) {
            if (F[i] < 0) { 
                damage_vector[i]=0.0;
            }else{
                double var1      = pow((k0/kappa[i]),beta1);
                double var2      = exp(-beta2*((kappa[i]-k0)/(k0)));
                damage_vector[i] = 1.0 - var1 * var2;              
            }
            if(damage_vector[i] < 0.0){damage_vector[i] = 0.0;}
        } 
        CalculateParameters(rParametersValues, damage_vector, EffStiffnessMatrix, dEprdE, dkdEpr);
        CalculatePartialDerivatives(r_material_properties, damage_vector, k0, beta1, beta2, kappa, dHdk);
        Matrix b = prod(dkdEpr, dEprdE);
        r_constitutive_matrix = EffStiffnessMatrix + prod(dHdk, b);
        r_stress_vector       = prod(r_constitutive_matrix,r_strain_vector);
        rDamageVector = damage_vector;
        }

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::GetDerivatives( 
    const Vector StressVector,
    Vector& dI1dS,
    Matrix& dJ2ddS,
    Vector& dJ2dS)
{
    KRATOS_TRY
    dI1dS[0] = dI1dS[1] = dI1dS[2]= 1.0; 
    dJ2ddS(0,0) = dJ2ddS(1,1) = dJ2ddS(2,2) =  2.0/3.0;
    dJ2ddS(0,1) = dJ2ddS(0,2) = dJ2ddS(1,0) = dJ2ddS(1,2) = dJ2ddS(2,0) = dJ2ddS(2,1)= -1.0/3.0;
    dJ2ddS(3,3) = dJ2ddS(4,4) = dJ2ddS(5,5) =  2.0;
    dJ2dS = prod(dJ2ddS, StressVector);
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::GetEigenValues(
    const Vector& StressVector,
    Vector& Pri_Values,
    double& MaxValue)
{
    KRATOS_TRY
    Matrix stress_matrix = ZeroMatrix(3,3);
    stress_matrix(0,1)= stress_matrix(1,0)= StressVector[5];
    stress_matrix(0,2)= stress_matrix(2,0)= StressVector[4];
    stress_matrix(2,1)= stress_matrix(1,2)= StressVector[3];
    for(int i=0; i<3; ++i){
        stress_matrix(i,i) = StressVector[i];
    }
    Matrix EigenVectors;
    Matrix EigenValues = ZeroMatrix(3,3);
    //prinicpal values, max and min
    MathUtils<double>::GaussSeidelEigenSystem(stress_matrix, EigenVectors, EigenValues, 1.0e-18, 20);
    Pri_Values[0] = EigenValues(0,0);
    Pri_Values[1] = EigenValues(1,1);
    Pri_Values[2] = EigenValues(2,2);
    MaxValue    = std::max({Pri_Values[0],Pri_Values[1],Pri_Values[2]}, [](const double &a, const double &b){return std::abs(b) > std::abs(a);});
    std::sort(Pri_Values.begin(), Pri_Values.end(), std::greater<double>());
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::ComputedSprdS(
    const Vector StressVector,
    const Vector Spr,
    Matrix& dSprdS)
{   
    KRATOS_TRY
    Matrix indx = ZeroMatrix(3, 3);
    Vector stress2(6);
    Vector dSprdS_entries(6);
    Vector dI1dS = ZeroVector(6);
    dI1dS[0] = dI1dS[1] = dI1dS[2]= 1.0; 

    indx(0,1)= indx(1,0)= indx(2,2)=1.0;
    indx(0,2)= indx(1,1)= indx(2,0)=2.0;
   
    Matrix stress_matrix = ZeroMatrix(3, 3);
    stress_matrix(0,1)= stress_matrix(1,0)= StressVector[5];
    stress_matrix(0,2)= stress_matrix(2,0)= StressVector[4];
    stress_matrix(2,1)= stress_matrix(1,2)= StressVector[3];
    for(int i=0; i<3; ++i){
        stress_matrix(i,i) = StressVector[i];
    }

    Matrix stress(stress_matrix);
    stress = prod(stress, stress);
    stress2[0]= stress(0,0);
    stress2[1]= stress(1,1); 
    stress2[2]= stress(2,2);
    stress2[3]= stress(1,2);
    stress2[4]= stress(0,2); 
    stress2[5]= stress(0,1);
    
    for(int i=0; i<3; ++i){  
        dSprdS_entries  = stress2 - ( Spr(indx(1,i)) + Spr(indx(2,i)) )*(StressVector) + Spr(indx(1,i)) * Spr(indx(2,i)) * dI1dS;
        dSprdS_entries /=  (Spr(indx(0,i)) - Spr(indx(1,i)) ) * ( Spr(indx(0,i)) - Spr(indx(2,i)) );
        for(int j=0; j<6; ++j){
            dSprdS(i,j)= dSprdS_entries[j];
        }
        for(int k=3; k<6; ++k){
            dSprdS(i,k) *= 2.;
        }   
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::GetInvariants(
    const Vector& StressVector,
    double& I1,
    double& J2)
{
    KRATOS_TRY
    I1 = StressVector(0)+StressVector(1)+StressVector(2);
    double I2 = (StressVector(0)*StressVector(1) + StressVector(1)*StressVector(2) + StressVector(2)* StressVector(0)) - ( StressVector(3)*StressVector(3) + StressVector(4)*StressVector(4) + StressVector(5)*StressVector(5) );
    J2 = (1.0/3.0)*(std::pow(I1,2.0))-I2;   
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::GetStressWeightFactor( double &w, const Vector &s_pr) const
{
    Vector N1(3);
    double eps = 1.e-8;
    int kk;
    for(kk=0; kk < 3; ++kk ){
        N1(kk) = 0.5 * ( fabs(s_pr(kk)) + s_pr(kk) );
    }
    double N11 = N1(0) + N1(1) + N1(2);
    double D11 = eps + fabs(s_pr(0)) + fabs(s_pr(1)) + fabs(s_pr(2));
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
    const Vector& DamageVector, 
    Matrix& DamageEffectTensor 
)
{
    KRATOS_TRY

    double D1 =  DamageVector[0];
    double D2 =  DamageVector[1];
    double D3 =  DamageVector[2];
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
    ConstitutiveLaw::Parameters& rParametersValues, 
    const Vector& DamageVector,
    Matrix& EffStiffnessMatrix,
    Matrix& dEprdE,
    Matrix& dkdEpr
)
{
    KRATOS_TRY

    Matrix& r_constitutive_matrix = rParametersValues.GetConstitutiveMatrix();
    const Vector& r_strain_vector = rParametersValues.GetStrainVector();        

    Matrix M     = ZeroMatrix(6,6);
    Matrix Inv_M = ZeroMatrix(6,6);
    Vector principal_strains = ZeroVector(3);
    double det_M, Max_principal_strain;


    CalculateElasticMatrix(r_constitutive_matrix, rParametersValues);
    GetDamageEffectTensor(DamageVector, M);
    MathUtils<double>::InvertMatrix(M, Inv_M, det_M);
    Matrix a = prod(r_constitutive_matrix,trans(Inv_M));
    EffStiffnessMatrix = prod(Inv_M,a);    

    GetEigenValues(r_strain_vector, principal_strains, Max_principal_strain);
    ComputedSprdS(r_strain_vector, principal_strains, dEprdE); 
    for(int i = 0; i < 3; ++i){  
        if(DamageVector[i] > 0){
            dkdEpr(i,i) = 1; 
        } 
    } 

    KRATOS_CATCH("") 
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::CalculatePartialDerivatives(
    const Properties& rMaterialProperties,
    const Vector& DamageVector,
    const double& Kappa0, 
    const double& Beta1, 
    const double& Beta2, 
    const Vector& Kappa,
    Matrix& dHdk
)
{
    const double E   = rMaterialProperties[YOUNG_MODULUS];
    const double nu  = rMaterialProperties[POISSON_RATIO];
    const double E_factor = ((1 + nu ) * (1- nu))/E;
    double dDdk;
    for(size_t i =0; i < dHdk.size1(); ++i){
        for(size_t j =0; j < DamageVector.size(); ++j){
            if(i==j){
                dDdk      = (1- DamageVector[j]) * (Beta1/Kappa[j] + Beta2/Kappa0);
                dHdk(j,j) = E_factor * (-2 * (1-DamageVector[j]) * (1-nu) * dDdk);
            }else if(i==0 || i==1 || i==2){
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
