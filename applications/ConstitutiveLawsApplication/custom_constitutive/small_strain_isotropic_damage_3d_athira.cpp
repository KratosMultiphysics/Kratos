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
#include "small_strain_isotropic_damage_3d_athira.h"
#include "constitutive_laws_application_variables.h"
#include "structural_mechanics_application_variables.h"
#include "constitutive_laws_application_variables.h"
#include "includes/checks.h"
#include "utilities/math_utils.h"
#include "math.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallStrainIsotropicDamageAthira3D::SmallStrainIsotropicDamageAthira3D()
    : ElasticIsotropic3D()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

SmallStrainIsotropicDamageAthira3D::SmallStrainIsotropicDamageAthira3D(const SmallStrainIsotropicDamageAthira3D &rOther)
    : ElasticIsotropic3D(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer SmallStrainIsotropicDamageAthira3D::Clone() const
{
    return Kratos::make_shared<SmallStrainIsotropicDamageAthira3D>(SmallStrainIsotropicDamageAthira3D(*this));
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

SmallStrainIsotropicDamageAthira3D::~SmallStrainIsotropicDamageAthira3D()
{
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

int SmallStrainIsotropicDamageAthira3D::Check(
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

bool SmallStrainIsotropicDamageAthira3D::Has(const Variable<double>& rThisVariable)
{
    if(rThisVariable == STRAIN_ENERGY){
        // explicitly returning "false", so the element calls CalculateValue(...)
        return false;
    } else if(rThisVariable == DAMAGE_VARIABLE){
        // explicitly returning "false", so the element calls CalculateValue(...)
        return false;

    }

    return false;
}

//************************************************************************************
//************************************************************************************

bool SmallStrainIsotropicDamageAthira3D::Has(const Variable<Vector>& rThisVariable)
{
    if(rThisVariable == INTERNAL_VARIABLES){
        return true;
    } else if(rThisVariable == STRAIN){
        // explicitly returning "false", so the element calls CalculateValue(...)
        return false;
    }

    return false;
}

//************************************************************************************
//************************************************************************************

Vector& SmallStrainIsotropicDamageAthira3D::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if(rThisVariable == INTERNAL_VARIABLES){
        rValue.resize(1);
        rValue[0] = mStrainVariable;
    }

    return rValue;
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamageAthira3D::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rProcessInfo
    )
{
    if(rThisVariable == INTERNAL_VARIABLES){
        mStrainVariable = rValue[0];
    }
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamageAthira3D::InitializeMaterial(
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

void SmallStrainIsotropicDamageAthira3D::FinalizeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rParametersValues)
{
    Vector internal_variables(1);
    this->CalculateStressResponse(rParametersValues, internal_variables);
    mStrainVariable = internal_variables[0];
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamageAthira3D::CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rParametersValues)
{
    Vector internal_variables(1);
    CalculateStressResponse(rParametersValues, internal_variables);
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamageAthira3D::CalculateStressResponse(
    ConstitutiveLaw::Parameters& rParametersValues,
    Vector& rInternalVariables)
{
    const Properties& r_material_properties = rParametersValues.GetMaterialProperties();
    Flags& r_constitutive_law_options = rParametersValues.GetOptions();
    Vector& r_strain_vector = rParametersValues.GetStrainVector();
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
        const double E   = r_material_properties[YOUNG_MODULUS];
        const double fck = r_material_properties[YIELD_STRESS_COMPRESSION];
        const double ft  = r_material_properties[YIELD_STRESS_TENSION];
        const double fcb = 1.16 * fck;  
        const double alphaL = (((fcb/fck)-1.0)/(2.0*(fcb/fck)-1.0));
        const double betaL  = ((1.0 - alphaL) * (fck/(ft))) - (1. + alphaL);          
        Matrix dJ2ddS = ZeroMatrix(6, 6);
        dJ2ddS(0,0) = dJ2ddS(1,1) = dJ2ddS(2,2) =  2.0/3.0;
        dJ2ddS(0,1) = dJ2ddS(0,2) = dJ2ddS(1,0) = dJ2ddS(1,2) = dJ2ddS(2,0) = dJ2ddS(2,1)= -1.0/3.0;
        dJ2ddS(3,3) = dJ2ddS(4,4) = dJ2ddS(5,5) =  2.0;

        Vector dI1dS = ZeroVector(6);
        dI1dS[0] = dI1dS[1] = dI1dS[2]= 1.0; 
    
        //converting stress to matrix form//
        Matrix r_stress_matrix = ZeroMatrix(3, 3);
        r_stress_matrix(0,1)= r_stress_matrix(1,0)= r_stress_vector[5];
        r_stress_matrix(0,2)= r_stress_matrix(2,0)= r_stress_vector[4];
        r_stress_matrix(2,1)= r_stress_matrix(1,2)= r_stress_vector[3];
        for(int i=0; i<3; ++i){
             r_stress_matrix(i,i) = r_stress_vector[i];
        }

        Matrix EigenVectorsMatrix, EigenValuesMatrix;
        //computing pricipal stress
        MathUtils<double>::GaussSeidelEigenSystem(r_stress_matrix, EigenVectorsMatrix, EigenValuesMatrix, 1.0e-18, 20);
        Vector Spr = ZeroVector(3);
        Matrix indx = ZeroMatrix(6, 6);
        indx(0,1)= indx(1,0)= indx(2,2)=1.0;
        indx(0,2)= indx(1,1)= indx(2,0)=2.0;
        Spr[0] = EigenValuesMatrix(0,0);
        Spr[1] = EigenValuesMatrix(1,1);
        Spr[2] = EigenValuesMatrix(2,2);
        
        double pr_str_max = std::max(std::max(EigenValuesMatrix(0,0),EigenValuesMatrix(1,1)),EigenValuesMatrix(2,2));
        Matrix dSprdS(3,6);
        Vector stress2(6);
        Vector dSprdS_entries(6);
        Matrix stress(r_stress_matrix);
        stress = prod(stress, stress);
        stress2[0]= stress(0,0);
        stress2[1]= stress(1,1); 
        stress2[2]= stress(2,2);
        stress2[3]= stress(1,2);
        stress2[4]= stress(0,2); 
        stress2[5]= stress(0,1);
        
        if( (fabs(Spr(0)-Spr(1)) < 2.0*eps ) && ( fabs(Spr(0)-Spr(2)) < 2.0*eps ) && ( fabs(Spr(1)-Spr(2)) < 2.0*eps ) )
        {
            Spr(0) = Spr(0) + 1.1*eps;
            Spr(1) = Spr(1) + 0.75*eps;
            Spr(2) = Spr(2) + 0.5*eps;
        }
        else{
            if( fabs(Spr(0)-Spr(1)) < eps ) Spr(1) = Spr(1) - eps;
            if( fabs(Spr(1)-Spr(2)) < eps ) Spr(2) = Spr(2) - eps;
        }

        for(int i=0; i<3; ++i){  
            dSprdS_entries  = stress2 - ( Spr(indx(1,i)) + Spr(indx(2,i)) )*(r_stress_vector) + Spr(indx(1,i)) * Spr(indx(2,i)) * dI1dS;
            dSprdS_entries /=  (Spr(indx(0,i)) - Spr(indx(1,i)) ) * ( Spr(indx(0,i)) - Spr(indx(2,i)) );
            for(int j=0; j<6; ++j){
                dSprdS(i,j)= dSprdS_entries[j];
            }
            for(int k=3; k<6; ++k){
                dSprdS(i,k) *= 2.;
            }   
        }
        Vector dSmaxdSig(6);
        for(int i=0; i<6; ++i){
          dSmaxdSig[i]= dSprdS(0,i); 
        }
        //computing invariants of stress
        double I1 = r_stress_matrix(0,0)+r_stress_matrix(1,1)+r_stress_matrix(2,2);
        double I2 = (r_stress_matrix(0,0)*r_stress_matrix(1,1) + r_stress_matrix(1,1)*r_stress_matrix(2,2) + r_stress_matrix(2,2)*r_stress_matrix(0,0)) - ( r_stress_matrix(1,2)*r_stress_matrix(1,2) + r_stress_matrix(0,2)*r_stress_matrix(0,2) + r_stress_matrix(0,1)*r_stress_matrix(0,1) );
        double J2 = (1.0/3.0)*(std::pow(I1,2.0))-I2;
        Vector dJ2dS  = prod(dJ2ddS, r_stress_vector);
        double D, DN, H, Eps_eq;
        
        if(pr_str_max < eps){
            H = 0.0;
        }else{
            H = 1.0;   
        }
        
        //local damage equivalent strain
        Eps_eq = (std::sqrt( 3. * J2 ) + alphaL * I1 + betaL * H * pr_str_max) /(E * ( 1. - alphaL)) ;
        
        const double beta1t = 0.85;
        const double beta2t = 0.18;
        const double beta1c = 0.0;
        const double beta2c = 0.095;

        const double k0t = fck/E;
        const double k0c = (10./3.) * ft/E;

        double k0 = k0t  * H  + (1.-H) * k0c;
        double beta1 = beta1t  * H  + (1.-H) * beta1c;
        double beta2 = beta2t  * H  + (1.-H) * beta2c;

        double Loc_e = std::fabs(Eps_eq);
        double kappa = std::max(Loc_e,k0);
        double f_d   = Loc_e - kappa;

        if (f_d < 0.0){
          D = 0.0;
         }else if (f_d == 0.0) {
       
            double var1      = pow((k0/kappa),beta1);
            double var2      = exp(-beta2*((kappa-k0)/(k0)));
            D                = 1.0 - var1 * var2;
            DN               = std::pow((1. -D),2);
            r_stress_vector *= DN;
 
            double factor       = 2.0*(D-1.0);         
            Vector dsdeii       = factor * prod(r_constitutive_matrix, r_strain_vector) ;
            double dKappadSmax  = betaL * H;
            Vector dKappadSig   = 1./(E * (1.-alphaL)) * (alphaL * dI1dS + 1.5/std::sqrt(3.*J2) * dJ2dS);
            dKappadSig         += 1./(E * (1.-alphaL)) * (dKappadSmax * dSmaxdSig);
            double dDdKappa     = (1.- D) * (beta1/kappa + beta2/k0);
            Matrix product(6,6);
            TensorProduct6(product,dsdeii,dKappadSig);
            r_constitutive_matrix = DN * r_constitutive_matrix + dDdKappa * prod(product,r_constitutive_matrix);
        }
        if(D < 0.0) D = 0.0;
        KRATOS_WATCH(H);
        KRATOS_WATCH(D);
        KRATOS_WATCH(r_stress_vector);
        KRATOS_WATCH(r_strain_vector);
        KRATOS_WATCH(r_constitutive_matrix);
        KRATOS_WATCH(rParametersValues.GetProcessInfo()[TIME]);
        KRATOS_WATCH("-------------------------------------------");
        //std::exit(-1);
    }
}

//************************************************************************************
//************************************************************************************


double& SmallStrainIsotropicDamageAthira3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParametersValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
        ElasticIsotropic3D::CalculateValue(rParametersValues, rThisVariable, rValue);
    return(rValue);
}

//************************************************************************************
//************************************************************************************

Vector& SmallStrainIsotropicDamageAthira3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParametersValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    //Explicitly having STRAIN and INITAL_STRAIN_VECTOR calculated in base class
    ElasticIsotropic3D::CalculateValue(rParametersValues, rThisVariable, rValue);
    return(rValue);
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamageAthira3D::GetLawFeatures(Features& rFeatures)
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



void SmallStrainIsotropicDamageAthira3D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mStrainVariable", mStrainVariable);
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamageAthira3D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mStrainVariable", mStrainVariable);
}

void SmallStrainIsotropicDamageAthira3D::AssembleSubMatrixToMatrix(
    Matrix& rOutput, 
    const Matrix& rInput, 
    const int StartRowIndex, 
    const int StartColIndex)
{
    KRATOS_TRY

    const std::size_t number_of_rows = rInput.size1();
    const std::size_t number_of_cols = rInput.size2();

    KRATOS_DEBUG_ERROR_IF(number_of_rows > rOutput.size1() || number_of_cols > rOutput.size2()) << "Size of matrices not matching  ";

    for (std::size_t i = 0; i < number_of_rows ; ++i) {
            for (std::size_t j = 0; j < number_of_cols ; ++j) {
                rOutput(StartRowIndex+i, StartColIndex+j) = rInput(i,j) ;
            }
        }

    KRATOS_CATCH("")
}
void SmallStrainIsotropicDamageAthira3D::TensorProduct6( 
    Matrix& rOutput,
    const Vector& rVector1,
    const Vector& rVector2)
{
    KRATOS_TRY
    KRATOS_DEBUG_ERROR_IF(rVector1.size()!=6 || rVector2.size()!=6) << "check the size of vectors  ";

    if (rOutput.size1() != 6 || rOutput.size1() != 6) {
        rOutput.resize(6, 6,false);
    }

    for(int i=0; i<6; ++i){
        for(int j=0; j<6; ++j){
            rOutput(i,j)= rVector1[i] * rVector2[j];
        }
    }

    KRATOS_CATCH("")    
}

} /* namespace Kratos.*/
