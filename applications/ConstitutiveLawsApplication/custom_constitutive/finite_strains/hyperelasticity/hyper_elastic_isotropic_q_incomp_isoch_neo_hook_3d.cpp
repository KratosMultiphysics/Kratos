// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//
// System includes
#include <iostream>

// External includes

// Project includes
#include "hyper_elastic_isotropic_q_incomp_isoch_neo_hook_3d.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D()
    : HyperElasticIsotropicNeoHookean3D()
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D(
    const HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D& rOther
    ) : HyperElasticIsotropicNeoHookean3D(rOther)
{
}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::Clone() const
{
    return Kratos::make_shared<HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D>(*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::
~HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D()
{
};

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::CalculateMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues
    )
{
    CalculateMaterialResponsePK2(rValues);

    Vector& r_stress_vector                = rValues.GetStressVector();
    const Matrix& r_deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f             = rValues.GetDeterminantF();

    TransformStresses(r_stress_vector, r_deformation_gradient_f, determinant_f, StressMeasure_PK2, StressMeasure_PK1);
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    KRATOS_TRY;

    // Get Values to compute the constitutive law:
    Flags &r_flags = rValues.GetOptions();

    const SizeType dimension = WorkingSpaceDimension();

    const Properties& material_properties  = rValues.GetMaterialProperties();
    Vector& r_strain_vector                = rValues.GetStrainVector();

    // The material properties
    const double young_modulus = material_properties[YOUNG_MODULUS];
    const double poisson_coefficient = material_properties[POISSON_RATIO];

    // The deformation gradient
    const Matrix& r_deformation_gradient_f = rValues.GetDeformationGradientF();
    double determinant_f = rValues.GetDeterminantF();
    KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;

    const double mu = young_modulus / (2.0 * (1.0 + poisson_coefficient));
    const double C1 = 0.5 * mu;

    Matrix C_tensor(dimension, dimension);

    if(r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangianStrain(rValues, r_strain_vector);
        noalias(C_tensor) = prod(trans(r_deformation_gradient_f), r_deformation_gradient_f);
    } else {
        Matrix strain_tensor(dimension, dimension);
        noalias(strain_tensor) = MathUtils<double>::StrainVectorToTensor(r_strain_vector);
        noalias(C_tensor) = 2.0 * strain_tensor + IdentityMatrix(dimension);
        determinant_f = std::sqrt(MathUtils<double>::Det(C_tensor));
    }

    if (r_flags.Is(COMPUTE_CONSTITUTIVE_TENSOR) || r_flags.Is(COMPUTE_STRESS))
        CalculateStressAndConstitutiveMatrixPK2(C_tensor, rValues.GetElementGeometry().GetValue(PRESSURE), C1, rValues.GetStressVector(), rValues.GetConstitutiveMatrix(), r_flags);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::CalculateMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
    CalculateMaterialResponsePK2(rValues);

    const auto &r_flags = rValues.GetOptions();
    if (r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        Vector &r_integrated_stress_vector = rValues.GetStressVector();
        Matrix stress_matrix(3, 3);
        noalias(stress_matrix) = MathUtils<double>::StressVectorToTensor(r_integrated_stress_vector);
        ContraVariantPushForward(stress_matrix, rValues.GetDeformationGradientF()); // Kirchhoff
        noalias(r_integrated_stress_vector) = MathUtils<double>::StressTensorToVector(stress_matrix, r_integrated_stress_vector.size());
    }
    if (r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        PushForwardConstitutiveMatrix(rValues.GetConstitutiveMatrix(), rValues.GetDeformationGradientF());
    }
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::CalculateMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues
    )
{
    CalculateMaterialResponseKirchhoff(rValues);

    Vector& r_stress_vector       = rValues.GetStressVector();
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
    const double determinant_f    = rValues.GetDeterminantF();

    // Set to Cauchy Stress:
    r_stress_vector       /= determinant_f;
    r_constitutive_matrix /= determinant_f;
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::InitializeMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::InitializeMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::InitializeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::InitializeMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::FinalizeMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::FinalizeMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::FinalizeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::FinalizeMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::GetLawFeatures(
    Features& rFeatures
    )
{
    //Set the type of law
    rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
    rFeatures.mOptions.Set( FINITE_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_GreenLagrange);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = VoigtSize;

    //Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

/***********************************************************************************/
/***********************************************************************************/

int HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS] <= 0.0) << "YOUNG_MODULUS is null or negative." << std::endl;

    const double tolerance = 1.0e-12;
    const double nu_upper_bound = 0.5;
    const double nu_lower_bound = -1.0;
    const double nu = rMaterialProperties[POISSON_RATIO];
    KRATOS_ERROR_IF((nu_upper_bound - nu) < tolerance) << "POISSON_RATIO is above the upper bound 0.5." << std::endl;
    KRATOS_ERROR_IF((nu - nu_lower_bound) < tolerance) << "POISSON_RATIO is below the lower bound -1.0." << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[DENSITY] < 0.0) << "DENSITY is negative." << std::endl;
    KRATOS_ERROR_IF_NOT(rElementGeometry.Has(PRESSURE)) << "Elemental PRESSURE is not set by the element..." << std::endl;

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/
void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::CalculateStressAndConstitutiveMatrixPK2(
    const Matrix& rC,
    const double Pressure,
    const double C1,
    Vector& rStress,
    Matrix &rTangentTensor,
    const Flags& rFlags
)
{
    KRATOS_TRY
    double v[308];
    const double one_third = 1.0 / 3.0;
    v[1]=rC(0,0);
    v[2]=rC(0,1);
    v[3]=rC(0,2);
    v[4]=rC(1,0);
    v[5]=rC(1,1);
    v[58]=-(v[2]*v[4])+v[1]*v[5];
    v[6]=rC(1,2);
    v[167]=-(v[3]*v[5])+v[2]*v[6];
    v[154]=v[3]*v[4]-v[1]*v[6];
    v[7]=rC(2,0);
    v[8]=rC(2,1);
    v[82]=-(v[5]*v[7])+v[4]*v[8];
    v[74]=v[2]*v[7]-v[1]*v[8];
    v[9]=rC(2,2);
    v[135]=v[3]*v[8]-v[2]*v[9];
    v[66]=v[6]*v[7]-v[4]*v[9];
    v[46]=-(v[3]*v[7])+v[1]*v[9];
    v[38]=-(v[6]*v[8])+v[5]*v[9];
    v[14]=v[1]*v[38]+v[2]*v[66]+v[3]*v[82];
    v[200]=0.5/std::sqrt(v[14]);
    v[201]=v[200];
    v[84]=v[200]*v[82];
    v[76]=v[200]*v[74];
    v[68]=v[201]*v[66];
    v[60]=v[201]*v[58];
    v[50]=v[201]*v[46];
    v[47]=1.0/(v[14]*v[14]);
    v[83]=-(v[47]*v[82]);
    v[75]=-(v[47]*v[74]);
    v[152]=v[38]*v[75];
    v[67]=-(v[47]*v[66]);
    v[165]=v[58]*v[67];
    v[59]=-(v[47]*v[58]);
    v[178]=v[38]*v[59];
    v[176]=v[135]*v[59];
    v[163]=v[46]*v[59];
    v[48]=-(v[46]*v[47]);
    v[146]=v[38]*v[48];
    v[40]=v[201]*v[38];
    v[39]=-(v[38]*v[47]);
    v[174]=v[154]*v[39];
    v[99]=v[1]/v[14]+v[163];
    v[63]=v[178]+v[5]/v[14];
    v[54]=v[146]+v[9]/v[14];
    v[10]=C1;
    v[11]=Pressure;
    v[226]=v[11]*v[40];
    v[218]=v[11]*v[60];
    v[217]=-(v[11]*v[84]);
    v[216]=v[11]*v[76];
    v[214]=v[11]*v[68];
    v[212]=v[11]*v[50];
    v[12]=std::sqrt(v[14]);
    v[202]=-(v[11]*v[12]);
    v[203]=((-4*one_third)*v[10])/std::pow(v[12],5*one_third);
    v[86]=v[203]*v[84];
    v[78]=v[203]*v[76];
    v[70]=v[203]*v[68];
    v[62]=v[203]*v[60];
    v[53]=v[203]*v[50];
    v[42]=v[203]*v[40];
    v[31]=(2*v[10])/std::pow(v[12],2.0*one_third);
    v[138]=-one_third*v[31];
    v[13]=v[1]+v[5]+v[9];
    v[207]=v[13];
    v[204]=-one_third*v[13];
    v[205]=v[204];
    v[225]=v[217]+v[204]*v[86];
    v[148]=v[204]*v[78];
    v[229]=v[148]-v[216];
    v[145]=v[205]*v[70];
    v[224]=v[145]-v[214];
    v[142]=v[138]+v[205]*v[62];
    v[227]=v[142]+v[218];
    v[139]=v[138]+v[205]*v[53];
    v[223]=v[139]+v[212];
    v[134]=v[138]+v[205]*v[42];
    v[220]=v[134]+v[226];
    v[117]=v[205]*v[99];
    v[114]=v[205]*v[63];
    v[93]=v[205]*v[54];
    v[34]=v[205]*v[31];
    v[222]=v[34]*v[48];
    v[219]=v[34]*v[58];
    v[215]=v[34]*v[46];
    v[211]=v[34]*v[38];
    v[228]=v[34]*v[59];
    v[221]=v[34]*v[39];
    v[208]=v[202]+v[34];
    v[15]=v[38]/v[14];
    v[206]=-(v[11]*v[15]);
    v[210]=v[206]*v[60]+v[202]*v[63];
    v[209]=v[206]*v[50]+v[202]*v[54];
    v[55]=-one_third*v[15];
    v[44]=1e0+v[207]*v[55];
    v[17]=v[135]/v[14];
    v[18]=v[167]/v[14];
    v[20]=v[46]/v[14];
    v[213]=-(v[20]*v[218])+v[202]*v[99];
    v[96]=-one_third*v[20];
    v[91]=1e0+v[207]*v[96];
    v[21]=v[154]/v[14];
    v[24]=v[58]/v[14];
    v[116]=-one_third*v[24];
    v[112]=1e0+v[116]*v[207];
    if (rFlags.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        rStress[0]=v[12]*v[206]+v[31]*v[44];
        rStress[1]=v[20]*v[202]+v[31]*v[91];
        rStress[2]=v[202]*v[24]+v[112]*v[31];
        rStress[3]=v[17]*v[208];
        rStress[4]=v[208]*v[21];
        rStress[5]=v[18]*v[208];
    }
    if (rFlags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        rTangentTensor(0,0)=2*(-(v[206]*v[40])+v[42]*v[44]+v[31]*(v[205]*v[38]*v[39]+v[55]));
        rTangentTensor(0,1)=2*(v[209]+v[44]*v[53]+v[31]*(v[55]+v[93]));
        rTangentTensor(0,2)=2*(v[210]+v[31]*(v[114]+v[55])+v[44]*v[62]);
        rTangentTensor(0,3)=v[211]*v[67]-v[206]*v[68]+v[44]*v[70];
        rTangentTensor(0,4)=v[206]*v[76]+v[44]*v[78]+v[208]*(v[152]-v[8]/v[14]);
        rTangentTensor(0,5)=v[211]*v[83]-v[206]*v[84]+v[44]*v[86];
        rTangentTensor(1,0)=2*(v[209]+v[42]*v[91]+v[31]*(v[93]+v[96]));
        rTangentTensor(1,1)=2*(v[20]*v[212]+v[53]*v[91]+v[31]*(v[205]*v[46]*v[48]+v[96]));
        rTangentTensor(1,2)=2*(v[213]+v[62]*v[91]+v[31]*(v[117]+v[96]));
        rTangentTensor(1,3)=v[20]*v[214]+v[215]*v[67]+v[70]*v[91];
        rTangentTensor(1,4)=v[20]*v[216]+v[215]*v[75]+v[78]*v[91];
        rTangentTensor(1,5)=v[20]*v[217]+v[208]*(-(v[7]/v[14])+v[46]*v[83])+v[86]*v[91];
        rTangentTensor(2,0)=2*(v[210]+(v[114]+v[116])*v[31]+v[112]*v[42]);
        rTangentTensor(2,1)=2*(v[213]+(v[116]+v[117])*v[31]+v[112]*v[53]);
        rTangentTensor(2,2)=2*(v[218]*v[24]+v[31]*(v[116]+v[205]*v[58]*v[59])+v[112]*v[62]);
        rTangentTensor(2,3)=-(v[214]*v[24])+v[208]*(v[165]-v[4]/v[14])+v[112]*v[70];
        rTangentTensor(2,4)=v[216]*v[24]+v[219]*v[75]+v[112]*v[78];
        rTangentTensor(2,5)=-(v[217]*v[24])+v[219]*v[83]+v[112]*v[86];
        rTangentTensor(3,0)=2*(v[17]*v[220]+v[135]*v[221]);
        rTangentTensor(3,1)=2*(v[135]*v[222]+v[17]*v[223]);
        rTangentTensor(3,2)=2*((v[176]-v[2]/v[14])*v[208]+v[17]*(v[142]-v[218]));
        rTangentTensor(3,3)=v[146]*v[208]+v[17]*v[224];
        rTangentTensor(3,4)=v[17]*(4*v[148]+v[216]);
        rTangentTensor(3,5)=v[152]*v[208]+v[17]*v[225];
        rTangentTensor(4,0)=2*(v[21]*(v[134]-v[226])+v[208]*(v[174]-v[6]/v[14]));
        rTangentTensor(4,1)=2*(v[154]*v[222]+v[21]*v[223]);
        rTangentTensor(4,2)=2*(v[21]*v[227]+v[154]*v[228]);
        rTangentTensor(4,3)=v[21]*(4*v[145]+v[214]);
        rTangentTensor(4,4)=v[163]*v[208]+v[21]*v[229];
        rTangentTensor(4,5)=v[165]*v[208]+v[21]*v[225];
        rTangentTensor(5,0)=2*(v[18]*v[220]+v[167]*v[221]);
        rTangentTensor(5,1)=2*(v[18]*(v[139]-v[212])+v[208]*(-(v[3]/v[14])+v[167]*v[48]));
        rTangentTensor(5,2)=2*(v[18]*v[227]+v[167]*v[228]);
        rTangentTensor(5,3)=v[174]*v[208]+v[18]*v[224];
        rTangentTensor(5,4)=v[176]*v[208]+v[18]*v[229];
        rTangentTensor(5,5)=v[178]*v[208]+v[18]*v[225];
    }
    KRATOS_CATCH("");
}

} // Namespace Kratos
