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
//                   Manuel Caicedo
//                   Alfredo Huespe
//  Collaborator:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "small_strain_plasticity_3d_ath.h"
#include "constitutive_laws_application_variables.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "includes/checks.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallStrainPlasticity3DAth::SmallStrainPlasticity3DAth()
    : ConstitutiveLaw()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

SmallStrainPlasticity3DAth::SmallStrainPlasticity3DAth(const SmallStrainPlasticity3DAth &rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer SmallStrainPlasticity3DAth::Clone() const
{
    return Kratos::make_shared<SmallStrainPlasticity3DAth>(SmallStrainPlasticity3DAth(*this));
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

SmallStrainPlasticity3DAth::~SmallStrainPlasticity3DAth()
{
}

//************************************************************************************
//************************************************************************************

bool SmallStrainPlasticity3DAth::Has(const Variable<double>& rThisVariable)
{
    if(rThisVariable == ACCUMULATED_PLASTIC_STRAIN){
        return true;
    }
    return false;
}

//************************************************************************************
//************************************************************************************

void SmallStrainPlasticity3DAth::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if(rThisVariable == ACCUMULATED_PLASTIC_STRAIN){
        mAccumulatedPlasticStrain = rValue;
    }
}

//************************************************************************************
//************************************************************************************

double& SmallStrainPlasticity3DAth::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if(rThisVariable == ACCUMULATED_PLASTIC_STRAIN){
        rValue = mAccumulatedPlasticStrain;
    }

    return rValue;
}

//************************************************************************************
//************************************************************************************

void SmallStrainPlasticity3DAth::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    mPlasticStrain = ZeroVector(this->GetStrainSize());
    mAccumulatedPlasticStrain = 0.0;
    mhmax = 0.0;
    mhmin = 0.0;
}

//************************************************************************************
//************************************************************************************

void SmallStrainPlasticity3DAth::InitializeMaterialResponseCauchy(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
}

//************************************************************************************
//************************************************************************************

void SmallStrainPlasticity3DAth::InitializeMaterialResponsePK2(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    // In small deformation is the same as compute Cauchy
    InitializeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainPlasticity3DAth::InitializeMaterialResponsePK1(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    // In small deformation is the same as compute Cauchy
    InitializeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainPlasticity3DAth::InitializeMaterialResponseKirchhoff(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    // In small deformation is the same as compute Cauchy
    InitializeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainPlasticity3DAth::FinalizeMaterialResponseCauchy(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    Vector plastic_strain;
    double accumulated_plastic_strain, hmax, hmin;
    this->CalculateStressResponse(rValues, plastic_strain, accumulated_plastic_strain, hmax, hmin);
    mPlasticStrain = plastic_strain;
    mAccumulatedPlasticStrain = accumulated_plastic_strain;
    mhmax = hmax;
    mhmin = hmin;
}

//************************************************************************************
//************************************************************************************

void SmallStrainPlasticity3DAth::FinalizeMaterialResponsePK2(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    // In small deformation is the same as compute Cauchy
    FinalizeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainPlasticity3DAth::FinalizeMaterialResponsePK1(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    // In small deformation is the same as compute Cauchy
    FinalizeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainPlasticity3DAth::FinalizeMaterialResponseKirchhoff(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    // In small deformation is the same as compute Cauchy
    FinalizeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainPlasticity3DAth::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainPlasticity3DAth::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // In small deformation is the same as compute Cauchy
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainPlasticity3DAth::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // In small deformation is the same as compute Cauchy
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainPlasticity3DAth::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    Vector plastic_strain;
    double accumulated_plastic_strain, hmax, hmin;
    this->CalculateStressResponse(rValues, plastic_strain, accumulated_plastic_strain, hmax, hmin);
}

//************************************************************************************
//************************************************************************************

void SmallStrainPlasticity3DAth::CalculateStressResponse(
    ConstitutiveLaw::Parameters& rValues,
    Vector& rPlasticStrain,
    double& rAccumulatedPlasticStrain,
    double& rhmax,
    double& rhmin)
{
    KRATOS_TRY

    rPlasticStrain.resize(6, false);
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    Flags& r_constitutive_law_options = rValues.GetOptions();
    Vector& r_strain_vector = rValues.GetStrainVector();
    if (rValues.GetProcessInfo().Has(INITIAL_STRAIN)) {
        noalias(r_strain_vector) += rValues.GetProcessInfo()[INITIAL_STRAIN];
    }

    // If we compute the tangent moduli or the stress
    if( r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS) ||
        r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        Vector& r_stress_vector = rValues.GetStressVector();
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        double trial_yield_function;
        rPlasticStrain = mPlasticStrain;
        rAccumulatedPlasticStrain = mAccumulatedPlasticStrain;
        rhmax = mhmax;
        rhmin = mhmin;

        Matrix elastic_tensor;
        elastic_tensor.resize(6, 6, false);
        CalculateElasticMatrix(r_material_properties, elastic_tensor);
        Vector trial_stress(6);
        noalias(trial_stress) = prod(elastic_tensor, (r_strain_vector - mPlasticStrain));
        trial_yield_function = this->YieldFunction(trial_stress, r_material_properties, mhmax, mhmin);
        if (trial_yield_function <= 0.) {
            // ELASTIC
            if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
                r_stress_vector = trial_stress;
            }
            if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
                r_constitutive_matrix = elastic_tensor;
            }
        } else {
            // INELASTIC
            Vector r = ZeroVector(7);
            Vector r1;
            Vector dv = ZeroVector(7);
            double dlamda = dv[6];
            Vector v = ZeroVector(7);
            Matrix drdv = ZeroMatrix(7,7);
            Matrix inverse_drdv = ZeroMatrix(7,7);
            v = AddSegmentToVector(v, trial_stress, 0);
            v[6] = mAccumulatedPlasticStrain;
            const int maxiter = 100;
            const double eps = 1e-8;
            double dFdlam;
            Matrix dSprdS(3,6);
            Vector delEpIn, dQdS, dFdS;
            Matrix dQ2ddS;
            double hmax = mhmax;
            double hmin = mhmin;
            Vector PlasticStrain = mPlasticStrain;
            Matrix rhsC = ZeroMatrix(7,6);
            for(int i=0; i<6; ++i) {rhsC(i,i) = -1.0;}
            //NR iterations
            for(int n = 0; n < maxiter; ++n){
                GetDerivatives(trial_stress, r_material_properties, PlasticStrain, dlamda, dQdS, dQ2ddS, dFdS, dFdlam, hmax, hmin);
                double det, det_drdv;
                Matrix inverse_C_tensor, drdv1;
                MathUtils<double>::InvertMatrix(elastic_tensor, inverse_C_tensor, det);
                r1 = r_strain_vector - prod(inverse_C_tensor, trial_stress) - PlasticStrain;
                r = AddSegmentToVector(r, r1, 0);
                r[6] = YieldFunction(trial_stress, r_material_properties, hmax, hmin);
                drdv1 = -inverse_C_tensor - (dlamda * dQ2ddS);
                AssembleSubMatrixToMatrix(drdv, drdv1, 0, 0);
                dQdS = -1.0 * dQdS;
                AddColumnToMatrix(drdv, dQdS, 0, 6);
                AddRowToMatrix(drdv, dFdS, 6, 0);
                drdv(6,6) = dFdlam;
                MathUtils<double>::InvertMatrix(drdv, inverse_drdv, det_drdv);
                dv = prod(inverse_drdv, r);
                dlamda = dv[6];
                v -= dv;
                for (size_t i = 0; i < trial_stress.size(); ++i) {trial_stress(i) = v(i);}
                mAccumulatedPlasticStrain = v[6];

                double norm1 = MathUtils<double>::Norm(dv);
                double norm2 = MathUtils<double>::Norm(v);
                if( ((norm1/norm2)< eps) ||(norm1< eps )) {
                    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
                        Matrix Sol;
                        Sol = prod(inverse_drdv, rhsC);
                        for (int i = 0; i < 6 ; ++i) {
                            for (int j = 0; j < 6 ; ++j) {
                                r_constitutive_matrix(i, j) = Sol(i,j);
                            }
                        }
                    }
                    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
                        for (size_t i = 0; i < r_stress_vector.size(); ++i) {r_stress_vector(i) = v(i);}
                    }
                }
                if( ((norm1/norm2)< eps) ||(norm1< eps ))  break;
            }
            // // We update the stress tensor
            rAccumulatedPlasticStrain = v[6];
            rPlasticStrain = PlasticStrain;
            rhmax = hmax;
            rhmin = hmin;
        }
    }

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************

double& SmallStrainPlasticity3DAth::CalculateValue(
    ConstitutiveLaw::Parameters& rValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    return(rValue);
}

//************************************************************************************
//************************************************************************************

Vector& SmallStrainPlasticity3DAth::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
        if (rThisVariable == STRAIN){
        const SizeType space_dimension = this->WorkingSpaceDimension();

        // Compute total deformation gradient
        const Matrix& r_F = rParameterValues.GetDeformationGradientF();
        KRATOS_DEBUG_ERROR_IF(r_F.size1()!= space_dimension || r_F.size2() != space_dimension)
            << "expected size of F " << space_dimension << "x" << space_dimension
            << ", got " << r_F.size1() << "x" << r_F.size2() << std::endl;

        const Matrix C_tensor = prod(trans(r_F),r_F);
        ConstitutiveLawUtilities<6>::CalculateGreenLagrangianStrain(C_tensor, rValue);

    } else if (rThisVariable == STRESSES ||
        rThisVariable == CAUCHY_STRESS_VECTOR ||
        rThisVariable == KIRCHHOFF_STRESS_VECTOR ||
        rThisVariable == PK2_STRESS_VECTOR) {

        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        // Set flags to only compute the stress
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, true );

        this->CalculateMaterialResponsePK2(rParameterValues);
        if (rValue.size() != GetStrainSize()) {
            rValue.resize(GetStrainSize());
        }
        noalias(rValue) = rParameterValues.GetStressVector();

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );

    } else if (rThisVariable == INITIAL_STRAIN_VECTOR) {
        if (this->HasInitialState()) {
	    if (rValue.size() != GetStrainSize()) {
	        rValue.resize(GetStrainSize());
	    }
	    noalias(rValue) = GetInitialState().GetInitialStrainVector();
        } else {
            noalias(rValue) = ZeroVector(0);
        }
    }
    return(rValue);
}

//************************************************************************************
//************************************************************************************

double SmallStrainPlasticity3DAth::YieldFunction(
    const Vector StressVector,
    const Properties& rMaterialProperties,
    const double hmax,
    const double hmin)
{
    const double eps = 1e-8;
    const double fck = rMaterialProperties[YIELD_STRESS_COMPRESSION];
    const double ft  = rMaterialProperties[YIELD_STRESS_TENSION];
    const double Hc = rMaterialProperties[HARDENING_MODULUS_COMPRESSION];
    const double Ht = rMaterialProperties[HARDENING_MODULUS_TENSION];
    const double fc0 =  1.0/3.0 * fck;
    const double fcb = 1.16 * fc0;
    const double alphaL = (((fcb/fc0)-1.0)/(2.0*(fcb/fc0)-1.0));
    Vector Spr = ZeroVector(3);
    double pr_str_max, pr_str_min;
    double H, I1, J2;
    GetEigenValues(StressVector,Spr, pr_str_max, pr_str_min);
    if(pr_str_max < eps){
        H = 0.0;
    }else{
        H = 1.0;
    }
    double fbarc  = fc0  + Hc * hmin;
    double fbart  = ft   + Ht * hmax;
    double betaL  = ((1.0 - alphaL) * (fbarc/(fbart))) - (1.0 + alphaL);
    GetInvariants(StressVector, I1, J2);
    return std::sqrt( 3.0 * J2 ) + alphaL * I1 + betaL * H * pr_str_max - fbarc * (1.0 - alphaL);
}

//************************************************************************************
//************************************************************************************

void SmallStrainPlasticity3DAth::CalculateElasticMatrix(
    const Properties &rMaterialProperties, Matrix &rElasticMatrix)
{
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double poisson_ratio = rMaterialProperties[POISSON_RATIO];
    const double lambda = E * poisson_ratio / ((1. + poisson_ratio) * (1. - 2. * poisson_ratio));
    const double mu = E / (2. + 2. * poisson_ratio);

    if (rElasticMatrix.size1() != 6 || rElasticMatrix.size2() != 6)
        rElasticMatrix.resize(6, 6, false);
    rElasticMatrix.clear();

    rElasticMatrix(0, 0) = lambda + 2. * mu;
    rElasticMatrix(0, 1) = lambda;
    rElasticMatrix(0, 2) = lambda;
    rElasticMatrix(1, 0) = lambda;
    rElasticMatrix(1, 1) = lambda + 2. * mu;
    rElasticMatrix(1, 2) = lambda;
    rElasticMatrix(2, 0) = lambda;
    rElasticMatrix(2, 1) = lambda;
    rElasticMatrix(2, 2) = lambda + 2. * mu;
    rElasticMatrix(3, 3) = mu;
    rElasticMatrix(4, 4) = mu;
    rElasticMatrix(5, 5) = mu;
}

//************************************************************************************
//************************************************************************************

void SmallStrainPlasticity3DAth::GetInvariants(
    const Vector& StressVector,
    double& I1,
    double& J2)
{
    I1 = StressVector(0)+StressVector(1)+StressVector(2);
    double I2 = (StressVector(0)*StressVector(1) + StressVector(1)*StressVector(2) + StressVector(2)* StressVector(0)) - ( StressVector(3)*StressVector(3) + StressVector(4)*StressVector(4) + StressVector(5)*StressVector(5) );
    J2 = (1.0/3.0)*(std::pow(I1,2.0))-I2;
}

//************************************************************************************
//************************************************************************************

void SmallStrainPlasticity3DAth::GetEigenValues(
    const Vector& StressVector,
    Vector& Pri_Values,
    double& MaxValue,
    double& MinValue)
{
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
}

//************************************************************************************
//************************************************************************************

void SmallStrainPlasticity3DAth::ComputedSprdS(
    const Vector StressVector,
    const Vector Spr,
    Matrix& dSprdS)
{
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
}

//************************************************************************************
//************************************************************************************
Vector& SmallStrainPlasticity3DAth::AddSegmentToVector(
    Vector& rOutput,
    const Vector& rInput,
    const size_t StartRowIndex)
{
    KRATOS_TRY
    for (std::size_t i = 0; i < rInput.size(); ++i) {
        rOutput(StartRowIndex+i) = rInput(i) ;
    }

    KRATOS_CATCH("")
    return rOutput;
}
//************************************************************************************
//************************************************************************************
void SmallStrainPlasticity3DAth::GetDerivatives(
    const Vector StressVector,
    const Properties &rMaterialProperties,
    Vector& PlasticStrain,
    const double& DelLamda,
    Vector& dQdS,
    Matrix& dQ2ddS,
    Vector& dFdS,
    double& dFdlam,
    double& hmax,
    double& hmin)
{
    const double eps  = 1e-8;
    const double fck = rMaterialProperties[YIELD_STRESS_COMPRESSION];
    const double ft  = rMaterialProperties[YIELD_STRESS_TENSION];
    const double Hc = rMaterialProperties[HARDENING_MODULUS_COMPRESSION];
    const double Ht = rMaterialProperties[HARDENING_MODULUS_TENSION];
    const double fc0 =  1.0/3.0 * fck;
    const double fcb = 1.16 * fc0;
    const double alphaL = (((fcb/fc0)-1.0)/(2.0*(fcb/fc0)-1.0));
    const double alphaQ = 0.15 * alphaL;

    Vector dI1dS = ZeroVector(6);
    Matrix dJ2ddS = ZeroMatrix(6, 6);
    Vector Spr = ZeroVector(3);
    Vector dEinpr = ZeroVector(3);
    Vector Einpr = ZeroVector(3);
    Vector dEinprmaxdEpsin = ZeroVector(6);
    Vector dEinprmindEpsin = ZeroVector(6);
    Vector dEinprmindS, dEinprmaxdS, dbetaLdS ;
    Matrix dEinprdEps = ZeroMatrix(3,6);
    Matrix dSprdS = ZeroMatrix(3,6);
    double H, I1, J2, SprMax, SprMin, w, EinprMax, EinprMin, dEinprMax, dEinprMin;

    dI1dS[0] = dI1dS[1] = dI1dS[2]= 1.0;
    dJ2ddS(0,0) = dJ2ddS(1,1) = dJ2ddS(2,2) =  2.0/3.0;
    dJ2ddS(0,1) = dJ2ddS(0,2) = dJ2ddS(1,0) = dJ2ddS(1,2) = dJ2ddS(2,0) = dJ2ddS(2,1)= -1.0/3.0;
    dJ2ddS(3,3) = dJ2ddS(4,4) = dJ2ddS(5,5) =  2.0;
    Vector dJ2dS = prod(dJ2ddS, StressVector);
    GetInvariants(StressVector, I1, J2);
    dQdS = alphaQ * dI1dS + 0.5  * (3.0/std::sqrt(3.0 * J2)) * dJ2dS;
    Matrix outerpro = ZeroMatrix(6,6);
    outerpro = outer_prod(dJ2dS, dJ2dS);
    dQ2ddS = 0.5  * (3.0/std::sqrt(3.0 * J2)) * dJ2ddS - 0.25 * (3.0*3.0/std::pow(3.0 * J2,(1.5))) * outer_prod(dJ2dS, dJ2dS);

    GetEigenValues(StressVector, Spr, SprMax, SprMin);
    if(SprMax < eps){
        H = 0.0;
    }else{
        H = 1.0;
    }

    //check and manipulate principal stresses
    if( ( std::abs(Spr(0)-Spr(1)) < 2.0*eps ) && ( std::abs(Spr(0)-Spr(2)) < 2.0*eps ) && ( std::abs(Spr(1)-Spr(2)) < 2.0*eps ) ){
        dSprdS(0,0) = dSprdS(1,1) = dSprdS(2,2) = 1.0;
    }else{
        if( std::abs(Spr(0)-Spr(1)) < eps ) Spr(1) = Spr(1) - eps;
        if( std::abs(Spr(1)-Spr(2)) < eps ) Spr(2) = Spr(2) - eps;
        ComputedSprdS(StressVector, Spr, dSprdS);
    }
    Vector dSmaxdSig(6);
    for(int i=0; i<6; ++i){
        dSmaxdSig[i]= dSprdS(0,i);
    }
    GetStressWeightFactor(w, Spr);
    Vector dPlasticStrain = (DelLamda * dQdS);
    GetEigenValues(dPlasticStrain, dEinpr, dEinprMax, dEinprMin);
    hmax += (w * dEinprMax);
    hmin -= ((1-w ) * dEinprMin);
    double fbarc  = fc0  + Hc * hmin;
    double fbart  = ft   + Ht * hmax;
    double betaL  = ((1.0 - alphaL) * (fbarc/(fbart))) - (1.0 + alphaL);

    PlasticStrain -= dPlasticStrain;
    GetEigenValues(PlasticStrain, Einpr, EinprMax, EinprMin);
    //check and manipulate principal strains
    if( ( std::abs(Einpr(0)-Einpr(1)) < 2.0*eps ) && ( std::abs(Einpr(0)-Einpr(2)) < 2.0*eps ) && ( std::abs(Einpr(1)-Einpr(2)) < 2.0*eps ) ){
        dEinprdEps(0,0) = dEinprdEps(1,1) = dEinprdEps(2,2) = 1.0;
    }else{
        if( std::abs(Einpr(0)-Einpr(1)) < eps ) Einpr(1) = Einpr(1) - eps;
        if( std::abs(Einpr(1)-Einpr(2)) < eps ) Einpr(2) = Einpr(2) - eps;
        ComputedSprdS(PlasticStrain, Einpr, dEinprdEps);
    }

    for(int i=0; i<6; ++i){
        dEinprmaxdEpsin[i]= dEinprdEps(0,i);
        dEinprmindEpsin[i]= dEinprdEps(2,i);
    }

    dEinprmaxdS = prod(dEinprmaxdEpsin, (DelLamda * dQ2ddS));
    dEinprmindS = prod(dEinprmindEpsin, (DelLamda * dQ2ddS));
    dbetaLdS = ((1.0-alphaL) * (Hc/fbart))*dEinprmindS;
    dbetaLdS -= (1.0-alphaL)*( (Ht*fbarc)/(fbart*fbart) )* dEinprmaxdS;
    dFdS = alphaL * dI1dS + 0.5  * (3.0/std::sqrt(3.0 * J2)) * dJ2dS + H * SprMax * dbetaLdS + betaL * H * dSmaxdSig  - (1.0- alphaL) * Hc * dEinprmindS;
    dFdlam  = H * SprMax * (1-alphaL) * Hc / fbart * inner_prod(dEinprmindEpsin, dQdS);
    dFdlam += H * SprMax * (1-alphaL) * (fbarc / (fbart*fbart)) * (-Ht) * inner_prod(dEinprmaxdEpsin,dQdS);
    dFdlam += - (1-alphaL) * Hc * inner_prod(dEinprmindEpsin, dQdS);
}

//************************************************************************************
//************************************************************************************
void SmallStrainPlasticity3DAth::GetStressWeightFactor( double &w, const Vector &s_pr) const
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

void SmallStrainPlasticity3DAth::AssembleSubMatrixToMatrix(
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
            rOutput(StartRowIndex+i, StartColIndex+j) = rInput(i,j);
        }
    }

    KRATOS_CATCH("")
}
//************************************************************************************
//************************************************************************************
void SmallStrainPlasticity3DAth::AddColumnToMatrix(
    Matrix& rOutput,
    const Vector& rInput,
    const int StartRowIndex,
    const int StartColIndex)
{
    KRATOS_TRY
    for (std::size_t i = 0; i < rInput.size(); ++i) {
        rOutput(StartRowIndex+i, StartColIndex) = rInput[i] ;
    }
    KRATOS_CATCH("")
}
//************************************************************************************
//************************************************************************************
void SmallStrainPlasticity3DAth::AddRowToMatrix(
    Matrix& rOutput,
    const Vector& rInput,
    const int StartRowIndex,
    const int StartColIndex)
{
    KRATOS_TRY
    for (std::size_t i = 0; i < rInput.size(); ++i) {
        rOutput(StartRowIndex, StartColIndex+i) = rInput[i] ;
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SmallStrainPlasticity3DAth::GetLawFeatures(Features& rFeatures)
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

int SmallStrainPlasticity3DAth::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(POISSON_RATIO));
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS_COMPRESSION));
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS_TENSION));
    KRATOS_CHECK(rMaterialProperties.Has(HARDENING_MODULUS_COMPRESSION));
    KRATOS_CHECK(rMaterialProperties.Has(HARDENING_MODULUS_TENSION));

    return 0;
}

//************************************************************************************
//************************************************************************************

void SmallStrainPlasticity3DAth::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mPlasticStrain", mPlasticStrain);
    rSerializer.save("mAccumulatedPlasticStrain", mAccumulatedPlasticStrain);
    rSerializer.save("mhmax", mhmax);
    rSerializer.save("mhmin", mhmin);
}

//************************************************************************************
//************************************************************************************

void SmallStrainPlasticity3DAth::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.load("mPlasticStrain", mPlasticStrain);
    rSerializer.load("mAccumulatedPlasticStrain", mAccumulatedPlasticStrain);
    rSerializer.load("mhmax", mhmax);
    rSerializer.load("mhmin", mhmin);
}

} /* namespace Kratos.*/
