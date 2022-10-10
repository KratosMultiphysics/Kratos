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
//  Main authors:    Marcelo Raschi
//  Collaborators:
//

// Project includes
#include "non_local_elastic_isotropic_damage.h"
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

NonLocalElasticIsotropicDamage::NonLocalElasticIsotropicDamage()
    : ElasticIsotropic3D()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

NonLocalElasticIsotropicDamage::NonLocalElasticIsotropicDamage(const NonLocalElasticIsotropicDamage &rOther)
    : ElasticIsotropic3D(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer NonLocalElasticIsotropicDamage::Clone() const
{
    return Kratos::make_shared<NonLocalElasticIsotropicDamage>(NonLocalElasticIsotropicDamage(*this));
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

NonLocalElasticIsotropicDamage::~NonLocalElasticIsotropicDamage()
{
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

int NonLocalElasticIsotropicDamage::Check(
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

bool NonLocalElasticIsotropicDamage::Has(const Variable<double>& rThisVariable)
{
    if(rThisVariable == STRAIN_ENERGY){
        // explicitly returning "false", so the element calls CalculateValue(...)
        return false;
    } else if(rThisVariable == DAMAGE_VARIABLE){
        return true;
    }

    return false;
}

//************************************************************************************
//************************************************************************************

bool NonLocalElasticIsotropicDamage::Has(const Variable<Vector>& rThisVariable)
{
    if(rThisVariable == INTERNAL_VARIABLES){
        // explicitly returning "false", so the element calls CalculateValue(...)s
        return false;
    } else if(rThisVariable == STRAIN){
        // explicitly returning "false", so the element calls CalculateValue(...)
        return false;
    }

    return false;
}

//************************************************************************************
//************************************************************************************

double& NonLocalElasticIsotropicDamage::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    KRATOS_TRY
    if(rThisVariable == DAMAGE_VARIABLE){
        rValue = mDamageVariable;
    }
    return rValue;
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void NonLocalElasticIsotropicDamage::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY
    if(rThisVariable == DAMAGE_VARIABLE){
        mDamageVariable = rValue;
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void NonLocalElasticIsotropicDamage::InitializeMaterial(
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

void NonLocalElasticIsotropicDamage::FinalizeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rParametersValues)
{
    KRATOS_TRY
    double damage_variable;
    this->CalculateStressResponse(rParametersValues, damage_variable);
    mDamageVariable = damage_variable;

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void NonLocalElasticIsotropicDamage::CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rParametersValues)
{
    KRATOS_TRY
    double damage_variable;
    CalculateStressResponse(rParametersValues, damage_variable);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void NonLocalElasticIsotropicDamage::CalculateStressResponse(
    ConstitutiveLaw::Parameters& rParametersValues,
    double& rDamageVariable)
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
        Matrix r_elastic_tensor;
        r_elastic_tensor.resize(6, 6, false);
        CalculateElasticMatrix(r_elastic_tensor, rParametersValues);
        noalias(r_stress_vector)      = prod(r_elastic_tensor, r_strain_vector);
        Vector HuDNL        = ZeroVector(6);
        Vector HDNLu        = ZeroVector(6);
        double HDNLDNL      = 0.0;
        const double eps = 1e-8;
        const double E   = r_material_properties[YOUNG_MODULUS];
        const double fck = r_material_properties[YIELD_STRESS_COMPRESSION];
        const double ft  = r_material_properties[YIELD_STRESS_TENSION];
        const double fcb = 1.16 * fck;  
        const double alphaL = (((fcb/fck)-1.0)/(2.0*(fcb/fck)-1.0));
        const double betaL  = ((1.0 - alphaL) * (fck/(ft))) - (1. + alphaL); 
        Vector dI1dS  = ZeroVector(6);
        Vector dJ2dS  = ZeroVector(6); 
        Matrix dJ2ddS = ZeroMatrix(6, 6); 
        Matrix dSprdS = ZeroMatrix(3,6);
        GetDerivatives( r_stress_vector, dI1dS, dJ2ddS,dJ2dS);
        Vector Spr = ZeroVector(3);;
        double SprMax, SprMin;
        GetEigenValues(r_stress_vector, Spr, SprMax, SprMin);
        
        //manipulation of zero entries in stress
        if( (fabs(Spr(0)-Spr(1)) < 2.0*eps ) && ( fabs(Spr(0)-Spr(2)) < 2.0*eps ) && ( fabs(Spr(1)-Spr(2)) < 2.0*eps ) )
        {
            Spr(0) = Spr(0) + 1.1*eps;
            Spr(1) = Spr(1) + 0.75*eps;
            Spr(2) = Spr(2) + 0.5*eps;
            dSprdS(0,0) = dSprdS(1,1) = dSprdS(2,2) = 1.0;
        }
        else{
            if( fabs(Spr(0)-Spr(1)) < eps ) Spr(1) = Spr(1) - eps;
            if( fabs(Spr(1)-Spr(2)) < eps ) Spr(2) = Spr(2) - eps;
            ComputedSprdS(r_stress_vector, Spr, dSprdS);
        }
        Vector dSmaxdSig(6);
        for(int i=0; i<6; ++i){
          dSmaxdSig[i]= dSprdS(0,i); 
        }
        //computing invariants of stress
        double I1, J2, D, DN, H, Eps_eq;
        GetInvariants(r_stress_vector, I1, J2);
        
        if(SprMax < eps){
            H = 0.0;
        }else{
            H = 1.0;   
        }
        //local damage equivalent strain
        Eps_eq = (std::sqrt( 3. * J2 ) + alphaL * I1 + betaL * H * SprMax) /(E * ( 1. - alphaL)) ;
        const double beta1t = 0.85;
        const double beta2t = 0.18;
        const double beta1c = 0.0;
        const double beta2c = 0.095;
        double k0t,k0c;

        const double k0t0 = r_material_properties[DAMAGE_THRESHOLD_TENSION];
        const double k0c0 = r_material_properties[DAMAGE_THRESHOLD_COMPRESSION];
        if (k0t0==0 || k0c0==0){
            k0t = fck/E;
            k0c = (10./3.) * ft/E;
        }else{
            k0t = std::min(k0t0,fck/E);
            k0c = std::min(k0c0,(10./3.) * ft/E);
        }
        double k0 = k0t  * H  + (1.-H) * k0c;
        double beta1 = beta1t  * H  + (1.-H) * beta1c;
        double beta2 = beta2t  * H  + (1.-H) * beta2c;

        double Loc_e = std::fabs(Eps_eq);
        double kappa = std::max(Loc_e,k0);
        double f_d   = Loc_e - kappa;

        if (f_d < 0.0){
            D = 0.0;
            //AssembleConstitutiveMatrix(r_constitutive_matrix, r_elastic_tensor, HDNLu, HuDNL, HDNLDNL);
        }else if (f_d == 0.0) {       
            double var1      = pow((k0/kappa),beta1);
            double var2      = exp(-beta2*((kappa-k0)/(k0)));
            D                = 1.0 - var1 * var2;
        }
        
        const auto& N    = rParametersValues.GetShapeFunctionsValues();
        double NL_Damage_GP = 0.0;
        for (size_t i = 0; i < N.size(); ++i) {
            NL_Damage_GP += N[i] * rParametersValues.GetElementGeometry()[i].FastGetSolutionStepValue(NON_LOCAL_DAMAGE, 0);
        }
        if (NL_Damage_GP == 0){

            AssembleConstitutiveMatrix(r_constitutive_matrix, r_elastic_tensor, HDNLu, HuDNL, HDNLDNL);
        
        }else{       
            DN                  = std::pow((1.0 - NL_Damage_GP),2);
            r_stress_vector    *= DN;
            //Huu = dSigmadEps; HuDNL = dSigmadDNL; HDNLu = dDlocaldEps; HDNLDNL = dDlocal/dDNL;        
            double dKappadSmax  = betaL * H;
            Vector dKappadSig   = 1./(E * (1.-alphaL)) * (alphaL * dI1dS + 1.5/std::sqrt(3.*J2) * dJ2dS);
            dKappadSig         += 1./(E * (1.-alphaL)) * (dKappadSmax * dSmaxdSig);
            double dDdKappa     = (1.- D) * (beta1/kappa + beta2/k0);

            r_elastic_tensor   *= DN ;
            HuDNL               = 2.0*(NL_Damage_GP-1.0) * prod(r_elastic_tensor, r_strain_vector) ;  
            HDNLu               = dDdKappa * prod(dKappadSig, r_elastic_tensor );
            HDNLDNL             = 0.0;
            AssembleConstitutiveMatrix(r_constitutive_matrix, r_elastic_tensor, HDNLu, HuDNL, HDNLDNL);
        }
        if(D < 0.0) D = 0.0;
        rDamageVariable = D;
        KRATOS_WATCH(r_stress_vector);
        KRATOS_WATCH(rParametersValues.GetProcessInfo()[TIME]);
        KRATOS_WATCH("-------------------------------------------");
        //std::exit(-1);
    }
    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void NonLocalElasticIsotropicDamage::GetDerivatives( 
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

void NonLocalElasticIsotropicDamage::GetEigenValues(
    const Vector& StressVector,
    Vector& Pri_Values,
    double& MaxValue,
    double& MinValue)
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
    std::sort(Pri_Values.begin(), Pri_Values.end(), std::greater<double>());
    MaxValue = std::max(std::max(EigenValues(0,0),EigenValues(1,1)),EigenValues(2,2));
    MinValue = std::min(std::min(EigenValues(0,0),EigenValues(1,1)),EigenValues(2,2));
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void NonLocalElasticIsotropicDamage::ComputedSprdS(
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

void NonLocalElasticIsotropicDamage::GetInvariants(
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

void NonLocalElasticIsotropicDamage::AssembleConstitutiveMatrix(
    Matrix& ConstitutiveMatrix, 
    const Matrix& Cuu, 
    const Vector& CDNLu, 
    const Vector& CuDNL, 
    const double& CDNLDNL
    )const
{
    KRATOS_TRY
    const std::size_t num_rows_C = Cuu.size1();
    const std::size_t num_cols_C = Cuu.size2();
    const std::size_t size_CuNL  = CuDNL.size();
    const std::size_t size_CNLu  = CDNLu.size();

    for (std::size_t i = 0; i < num_rows_C ; ++i) {
        for (std::size_t j = 0; j < num_cols_C ; ++j) {
            ConstitutiveMatrix(i, j) = Cuu(i,j);
        }
    }
    for (std::size_t i = 0; i < size_CuNL; ++i) {
        ConstitutiveMatrix(i, num_cols_C) = CuDNL[i] ;
    }
    for (std::size_t i = 0; i < size_CNLu; ++i) {
        ConstitutiveMatrix(num_rows_C, i) = CDNLu[i] ;
    }
    ConstitutiveMatrix(num_rows_C, num_cols_C) = CDNLDNL;
    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

double& NonLocalElasticIsotropicDamage::CalculateValue(
    ConstitutiveLaw::Parameters& rParametersValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    KRATOS_TRY
    if (rThisVariable == DAMAGE_VARIABLE) {
        double damage_variable;
        this->CalculateStressResponse( rParametersValues, damage_variable);
        rValue = damage_variable;
    }
    return( rValue );
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

Vector& NonLocalElasticIsotropicDamage::CalculateValue(
    ConstitutiveLaw::Parameters& rParametersValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    KRATOS_TRY
    //Explicitly having STRAIN and INITAL_STRAIN_VECTOR calculated in base class
    ElasticIsotropic3D::CalculateValue(rParametersValues, rThisVariable, rValue);
    return(rValue);
    KRATOS_CATCH("")
}

//************************************************************************************
//***********************************************************************************

void NonLocalElasticIsotropicDamage::GetLawFeatures(Features& rFeatures)
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



void NonLocalElasticIsotropicDamage::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mDamageVariable", mDamageVariable);
}

//************************************************************************************
//************************************************************************************

void NonLocalElasticIsotropicDamage::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mDamageVariable", mDamageVariable);
}



} /* namespace Kratos.*/
