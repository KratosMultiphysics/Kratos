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
//  Main authors:   Athira Vadakkekkara
//  Collaborators:
//

// System includes
#include <algorithm>

// Project includes
#include "small_strain_isotropic_damage_3d_mazars.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"
#include "constitutive_laws_application_variables.h"
#include "structural_mechanics_application_variables.h"
#include "constitutive_laws_application_variables.h"
#include "includes/checks.h"
#include "utilities/math_utils.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallStrainIsotropicDamage3DMazars::SmallStrainIsotropicDamage3DMazars()
    : ElasticIsotropic3D()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

SmallStrainIsotropicDamage3DMazars::SmallStrainIsotropicDamage3DMazars(const SmallStrainIsotropicDamage3DMazars &rOther)
    : ElasticIsotropic3D(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer SmallStrainIsotropicDamage3DMazars::Clone() const
{
    return Kratos::make_shared<SmallStrainIsotropicDamage3DMazars>(SmallStrainIsotropicDamage3DMazars(*this));
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

SmallStrainIsotropicDamage3DMazars::~SmallStrainIsotropicDamage3DMazars()
{
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

int SmallStrainIsotropicDamage3DMazars::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(POISSON_RATIO));
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS_TENSION));
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS_COMPRESSION));
    return 0;
}

//************************************************************************************
//************************************************************************************

bool SmallStrainIsotropicDamage3DMazars::Has(const Variable<double>& rThisVariable)
{
    if(rThisVariable == STRAIN_ENERGY){
        // explicitly returning "false", so the element calls CalculateValue(...)
        return false;
    } else if(rThisVariable == DAMAGE_VARIABLE){
        // explicitly returning "false", so the element calls CalculateValue(...)
        return true;

    }

    return false;
}

//************************************************************************************
//************************************************************************************

bool SmallStrainIsotropicDamage3DMazars::Has(const Variable<Vector>& rThisVariable)
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

double& SmallStrainIsotropicDamage3DMazars::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if(rThisVariable == INTERNAL_VARIABLES){
        rValue = mInternalVariables[0];
    }else if(rThisVariable == DAMAGE_VARIABLE){
        rValue = mInternalVariables[1];
    }

    return rValue;
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamage3DMazars::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rProcessInfo
    )
{
    if(rThisVariable == INTERNAL_VARIABLES){
        mInternalVariables[0] = rValue;
    }else if(rThisVariable == DAMAGE_VARIABLE){
        mInternalVariables[1] = rValue;
    }
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamage3DMazars::InitializeMaterial(
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

void SmallStrainIsotropicDamage3DMazars::FinalizeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rParametersValues)
{
    Vector internal_variables(2);
    this->CalculateStressResponse(rParametersValues, internal_variables);
    mInternalVariables = internal_variables;
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamage3DMazars::CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rParametersValues)
{
    Vector internal_variables(2);
    CalculateStressResponse(rParametersValues, internal_variables);
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamage3DMazars::CalculateStressResponse(
    ConstitutiveLaw::Parameters& rParametersValues,
    Vector& rInternalVariables)
{
    const Properties& r_material_properties = rParametersValues.GetMaterialProperties();
    Flags& r_constitutive_law_options = rParametersValues.GetOptions();
    Vector& r_strain_vector = rParametersValues.GetStrainVector();
    CalculateValue(rParametersValues, STRAIN, r_strain_vector);
    rInternalVariables = mInternalVariables;
    // If we compute the tangent moduli or the stress
    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ||
        r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ))
        {
        Vector& r_stress_vector       = rParametersValues.GetStressVector();
        const Vector& r_strain_vector = rParametersValues.GetStrainVector();
        Matrix& r_constitutive_matrix = rParametersValues.GetConstitutiveMatrix();
        CalculateElasticMatrix(r_constitutive_matrix, rParametersValues);
        noalias(r_stress_vector)      = prod(r_constitutive_matrix, r_strain_vector);
        KRATOS_WATCH(r_strain_vector)
        BoundedMatrix3x6Type derivatives_of_eigen_values;
        BoundedVectorType Spr = ZeroVector(3);
        BoundedVectorType principal_strains = ZeroVector(3);
        double r, D, local_equivalent_strain, k0t, k0c, SprMax, max_principal_strain;
        const double eps = 1e-8;
        const double E   = r_material_properties[YOUNG_MODULUS];
        const double ft  = r_material_properties[YIELD_STRESS_TENSION];
        const double beta1t = r_material_properties[DAMAGE_MODEL_PARAMETER_BETA1_TENSION];
        const double beta2t = r_material_properties[DAMAGE_MODEL_PARAMETER_BETA2_TENSION];
        const double beta1c = r_material_properties[DAMAGE_MODEL_PARAMETER_BETA1_COMPRESSION];
        const double beta2c = r_material_properties[DAMAGE_MODEL_PARAMETER_BETA2_COMPRESSION];
        GetEigenValues(Spr, SprMax, STRESSES, r_stress_vector);
        const double H = (SprMax < eps) ? 0.0 : 1.0;
        GetEigenValues(principal_strains, max_principal_strain, STRAIN, r_strain_vector);

        //local equivalent strain
        local_equivalent_strain = sqrt(pow(MacaulayBrackets(principal_strains[0]),2) + pow(MacaulayBrackets(principal_strains[1]),2) + pow(MacaulayBrackets(principal_strains[2]),2));
        GetStressWeightFactor(r,Spr);
        const double del_r = (r > 0 && r < 0.1) ? 1.0 : 0.0;

        if(r_material_properties.Has(DAMAGE_THRESHOLD_TENSION)==true){
            k0t = r_material_properties[DAMAGE_THRESHOLD_TENSION];
        }else{
            k0t = ft/E;
        }
        if(r_material_properties.Has(DAMAGE_THRESHOLD_COMPRESSION)==true){
            k0c = r_material_properties[DAMAGE_THRESHOLD_COMPRESSION];
        }else{
            k0c = ft/E;
        }

        const double k0 = k0t  * H * (1-del_r) + (1.0- H + del_r) * k0c;
        const double beta1 = beta1t  * H * (1-del_r) + (1.-H + del_r) * beta1c;
        const double beta2 = beta2t  * H * (1-del_r) + (1.-H + del_r) * beta2c;
        const double max_e = std::max(mInternalVariables[0], local_equivalent_strain);
        const double kappa = std::max(max_e, k0);

        if ((kappa >= 0.0) && (kappa <= k0)){
            D = 0.0;
        }else if (kappa > k0) {
            const double var1     = pow((k0/kappa),beta1);
            const double var2     = exp(-beta2*((kappa-k0)/(k0)));
            D                     = 1.0 - var1 * var2;
            const double DN       = pow((1.0 - D),2);
            r_stress_vector      *= DN;
            const double factor   = 2.0*(1.0 - D);
            BoundedVectorVoigtType dsde           = factor * prod(r_constitutive_matrix, r_strain_vector) ;
            const double dDdKappa                 = (1.0- D) * (beta1/kappa + beta2/k0);
            BoundedVectorVoigtType dKappadEpsilon;
            CalculateDerivativesofEigenvalues(derivatives_of_eigen_values, principal_strains, r_strain_vector, STRAIN);


            const double dKappadEpseq = 1.0;
            BoundedVectorType dEpseqdEpspr = ZeroVector(3);
            for(SizeType i = 0; i < Dimension; ++i){
             dEpseqdEpspr[i] = MacaulayBrackets(principal_strains[i]) /local_equivalent_strain;
            }
            BoundedVectorVoigtType dEpseqdEps = prod(dEpseqdEpspr,derivatives_of_eigen_values);
            dKappadEpsilon = dKappadEpseq * dEpseqdEps;
            BoundedVectorVoigtType dDdEpsilon = dDdKappa * dKappadEpsilon;
            BoundedMatrixVoigtType product = ZeroMatrix(6,6);
            TensorProduct6(product, dsde, dDdEpsilon);
            KRATOS_WATCH(DN * r_constitutive_matrix )
            r_constitutive_matrix = DN * r_constitutive_matrix - product;

            KRATOS_WATCH(product)
            KRATOS_WATCH(dsde)
            KRATOS_WATCH(dDdKappa)
            KRATOS_WATCH(dEpseqdEpspr)
            KRATOS_WATCH(derivatives_of_eigen_values)
            KRATOS_WATCH(dEpseqdEps)

            KRATOS_WATCH(r_constitutive_matrix)

            // for(SizeType i = 0; i < Dimension; ++i){
            //     if (MacaulayBrackets(principal_strains[i])!=0.0){
            //         for (SizeType j = 0; j < VoigtSize; ++j){
            //             dKappadEpsilon[j]   +=  principal_strains[i] * derivatives_of_eigen_values(i,j);
            //         }
            //     }
            // }
            // dKappadEpsilon /= local_equivalent_strain ;
            // BoundedMatrixVoigtType product(6,6);
            // TensorProduct6(product,dsde,dKappadEpsilon);
            // r_constitutive_matrix = DN * r_constitutive_matrix + dDdKappa * prod(product,r_constitutive_matrix);
            // if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            //     this->CalculateTangentTensor(rParametersValues);
            //     KRATOS_WATCH(rParametersValues.GetConstitutiveMatrix())
            // }
        }else{
            KRATOS_ERROR << "check the damage loading function" << std::endl;
        }

        KRATOS_WATCH(D)
        KRATOS_WATCH(r_stress_vector)
        if(D < 0.0) D = 0.0;
        rInternalVariables[0] = max_e;
        rInternalVariables[1] = D;
        KRATOS_WATCH(rParametersValues.GetProcessInfo()[TIME]);
        KRATOS_WATCH("-------------------------------------------");
        //std::exit(-1);
    }
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamage3DMazars::CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ? r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ? static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    if (tangent_operator_estimation == TangentOperatorEstimation::Analytic) {
        KRATOS_ERROR << "Analytic solution not available" << std::endl;
    } else if (tangent_operator_estimation == TangentOperatorEstimation::FirstOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (first order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 1);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 2);
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

double& SmallStrainIsotropicDamage3DMazars::CalculateValue(
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

Vector& SmallStrainIsotropicDamage3DMazars::CalculateValue(
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

void SmallStrainIsotropicDamage3DMazars::TensorProduct6(
    BoundedMatrixVoigtType& rOutput,
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

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamage3DMazars::GetDerivatives(
    const Vector StressVector,
    Vector& dI1dS,
    Matrix& dJ2ddS,
    Vector& dJ2dS)
{
    dI1dS[0] = dI1dS[1] = dI1dS[2]= 1.0;
    dJ2ddS(0,0) = dJ2ddS(1,1) = dJ2ddS(2,2) =  2.0/3.0;
    dJ2ddS(0,1) = dJ2ddS(0,2) = dJ2ddS(1,0) = dJ2ddS(1,2) = dJ2ddS(2,0) = dJ2ddS(2,1)= -1.0/3.0;
    dJ2ddS(3,3) = dJ2ddS(4,4) = dJ2ddS(5,5) =  2.0;
    dJ2dS = prod(dJ2ddS, StressVector);
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamage3DMazars::GetEigenValues(
    BoundedVectorType& Principal_Values_vector,
    double& MaxValue,
    const Variable<Vector>& rThisVariable,
    const Vector& VectorForm)
{
    KRATOS_TRY
    BoundedMatrixType MatrixForm = ZeroMatrix(3,3);
    BoundedMatrixType EigenVectors;
    BoundedMatrixType EigenValues = ZeroMatrix(3,3);
    VectorToTensor(MatrixForm, VectorForm, rThisVariable);
    MathUtils<double>::GaussSeidelEigenSystem(MatrixForm, EigenVectors, EigenValues);
    KRATOS_WATCH(EigenVectors)
    Principal_Values_vector[0] = EigenValues(0,0);
    Principal_Values_vector[1] = EigenValues(1,1);
    Principal_Values_vector[2] = EigenValues(2,2);
    MaxValue    = std::max({Principal_Values_vector[0],Principal_Values_vector[1],Principal_Values_vector[2]}, [](const double a, const double b){return std::abs(b) > std::abs(a);});
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamage3DMazars::ComputedSprdS(
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

void SmallStrainIsotropicDamage3DMazars::GetInvariants(
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

void SmallStrainIsotropicDamage3DMazars::VectorToTensor(
    BoundedMatrixType& TensorForm,
    const Vector& VectorForm,
    const Variable<Vector>& rThisVariable
)
{
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
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamage3DMazars::GetStressWeightFactor(
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

void SmallStrainIsotropicDamage3DMazars::CalculateDerivativesofEigenvalues(
     BoundedMatrix3x6Type &DerivativesofEigenvalues,
     BoundedVectorType &EigenvaluesVector,
     const BoundedVectorVoigtType &Voigtform,
     const Variable<Vector>& rThisVariable
    )
{
    BoundedMatrixType Matrixform;
    const double eps = 1e-8;
    VectorToTensor(Matrixform, Voigtform, rThisVariable);
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

    // Method 1
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
    //     DerivativesMatrix= (1/trace) * cofactor_matrix;
    //     DerivativesofEigenvalues(i,0) = DerivativesMatrix(0,0);
    //     DerivativesofEigenvalues(i,1) = DerivativesMatrix(1,1);
    //     DerivativesofEigenvalues(i,2) = DerivativesMatrix(2,2);
    //     DerivativesofEigenvalues(i,3) = DerivativesMatrix(1,2);
    //     DerivativesofEigenvalues(i,4) = DerivativesMatrix(0,2);
    //     DerivativesofEigenvalues(i,5) = DerivativesMatrix(0,1);
    // }


    // Method 2
    // if( (fabs(EigenvaluesVector(0)-EigenvaluesVector(1)) < 2.0*eps ) && ( fabs(EigenvaluesVector(0)-EigenvaluesVector(2)) < 2.0*eps ) && ( fabs(EigenvaluesVector(1)-EigenvaluesVector(2)) < 2.0*eps ) ){
    //     EigenvaluesVector(0) = EigenvaluesVector(0) + 1.1*eps;
    //     EigenvaluesVector(1) = EigenvaluesVector(1) + 0.75*eps;
    //     EigenvaluesVector(2) = EigenvaluesVector(2) + 0.5*eps;
    //     DerivativesofEigenvalues(0,0) = 1.0;
    //     DerivativesofEigenvalues(1,1) = 1.0;
    //     DerivativesofEigenvalues(2,2) = 1.0;
    // }else{
    //     if( fabs(EigenvaluesVector(0)-EigenvaluesVector(1)) < eps ) EigenvaluesVector(1) = EigenvaluesVector(1) - eps;
    //     if( fabs(EigenvaluesVector(1)-EigenvaluesVector(2)) < eps ) EigenvaluesVector(2) = EigenvaluesVector(2) - eps;

    // }
    // Matrix indx = ZeroMatrix(3, 3);
    // Vector Vectorform_squared(6);
    // Vector dSprdS_entries(6);
    // Vector dI1dS = ZeroVector(6);
    // dI1dS[0] = dI1dS[1] = dI1dS[2]= 1.0;

    // indx(0,1)= indx(1,0)= indx(2,2)=1.0;
    // indx(0,2)= indx(1,1)= indx(2,0)=2.0;
    // Matrix Matrixform_squared = ZeroMatrix(3,3);
    // Matrixform_squared = prod(Matrixform, Matrixform);
    // Vectorform_squared[0]= Matrixform_squared(0,0);
    // Vectorform_squared[1]= Matrixform_squared(1,1);
    // Vectorform_squared[2]= Matrixform_squared(2,2);
    // Vectorform_squared[3]= Matrixform_squared(1,2);
    // Vectorform_squared[4]= Matrixform_squared(0,2);
    // Vectorform_squared[5]= Matrixform_squared(0,1);

    // for(int i=0; i<3; ++i){
    //     dSprdS_entries  = Vectorform_squared - ( EigenvaluesVector(indx(1,i)) + EigenvaluesVector(indx(2,i)) )*(Voigtform) + EigenvaluesVector(indx(1,i)) * EigenvaluesVector(indx(2,i)) * dI1dS;
    //     dSprdS_entries /=  (EigenvaluesVector(indx(0,i)) - EigenvaluesVector(indx(1,i)) ) * ( EigenvaluesVector(indx(0,i)) - EigenvaluesVector(indx(2,i)) );
    //     for(int j=0; j<6; ++j){
    //         DerivativesofEigenvalues(i,j)= dSprdS_entries[j];
    //     }
    //     for(int k=3; k<6; ++k){
    //         DerivativesofEigenvalues(i,k) *= 2.;
    //     }
    // }

}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamage3DMazars::GetLawFeatures(Features& rFeatures)
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

void SmallStrainIsotropicDamage3DMazars::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mInternalVariables", mInternalVariables);
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamage3DMazars::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.load("mInternalVariables", mInternalVariables);
}


} /* namespace Kratos.*/
