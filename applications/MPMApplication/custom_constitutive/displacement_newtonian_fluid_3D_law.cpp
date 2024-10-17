//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Contri Alessandro
//
//  References:    This class is adapted from applications/MPMApplication/custom_constitutive/hyperelastic_3D_law.cpp


// System includes

// External includes

// Project includes
#include "custom_constitutive/displacement_newtonian_fluid_3D_law.hpp"

#include "mpm_application_variables.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

DispNewtonianFluid3DLaw::DispNewtonianFluid3DLaw()
    : ConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

DispNewtonianFluid3DLaw::DispNewtonianFluid3DLaw(const DispNewtonianFluid3DLaw& rOther)
    : ConstitutiveLaw(rOther)
    ,mInverseDeformationGradientF0(rOther.mInverseDeformationGradientF0)
    ,mDeterminantF0(rOther.mDeterminantF0)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer DispNewtonianFluid3DLaw::Clone() const
{
    return Kratos::make_shared<DispNewtonianFluid3DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

DispNewtonianFluid3DLaw::~DispNewtonianFluid3DLaw()
{
}

//*******************************OPERATIONS FROM BASE CLASS***************************
//************************************************************************************

//***********************HAS : DOUBLE - VECTOR - MATRIX*******************************
//************************************************************************************
bool DispNewtonianFluid3DLaw::Has(const Variable<double>& rThisVariable)
{
    if (rThisVariable==DYNAMIC_VISCOSITY)
    {
        return true;
    }
    else if (rThisVariable == PRESSURE_COEFFICIENT)
    {
        return true;
    }
    else if (rThisVariable == BULK_MODULUS)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool DispNewtonianFluid3DLaw::Has(const Variable<Vector>& rThisVariable)
{
    return false;
}

bool DispNewtonianFluid3DLaw::Has(const Variable<Matrix>& rThisVariable)
{
    return false;
}

//******************CALCULATE VALUE: DOUBLE - VECTOR - MATRIX*************************
//************************************************************************************

double& DispNewtonianFluid3DLaw::CalculateValue(Parameters& rParametersValues, const Variable<double>& rThisVariable, double& rValue)
{
    return (this->GetValue(rThisVariable, rValue));
}


//***********************GET VALUE: DOUBLE - VECTOR - MATRIX**************************
//************************************************************************************

double& DispNewtonianFluid3DLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
{
    rValue=0;
    return(rValue);
}

Vector& DispNewtonianFluid3DLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    return(rValue);
}

Matrix& DispNewtonianFluid3DLaw::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
{
    return(rValue);
}


//***********************SET VALUE: DOUBLE - VECTOR - MATRIX**************************
//************************************************************************************

void DispNewtonianFluid3DLaw::SetValue(const Variable<double>& rThisVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo)
{
    if (rThisVariable == DETERMINANT_F)
        mDeterminantF0 = rValue;
}

void DispNewtonianFluid3DLaw::SetValue(const Variable<Vector>& rThisVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
{

}

void DispNewtonianFluid3DLaw::SetValue(const Variable<Matrix>& rThisVariable, const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
{

}

//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************

void DispNewtonianFluid3DLaw::InitializeMaterial(const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues)
{
  mDeterminantF0                = 1;
  mInverseDeformationGradientF0 = IdentityMatrix(3);
}

//***********************MATERIAL RESPONSE********************************************
//************************************************************************************

void DispNewtonianFluid3DLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{

    //-----------------------------//

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

   //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const Properties& MaterialProperties=rValues.GetMaterialProperties();
    const Matrix& DeformationGradientF=rValues.GetDeformationGradientF();
    const double& DeterminantF=rValues.GetDeterminantF();

    Vector& StrainVector=rValues.GetStrainVector();
    Vector& StressVector=rValues.GetStressVector();
    Matrix& ConstitutiveMatrix=rValues.GetConstitutiveMatrix();

    const GeometryType& domain_geometry = rValues.GetElementGeometry();
    const Vector& shape_functions = rValues.GetShapeFunctionsValues();
    const ProcessInfo& current_process_info = rValues.GetProcessInfo();

    //0.- Initialize parameters
    MaterialResponseVariables ViscousVariables;
    ViscousVariables.Identity=IdentityMatrix(3);

    ViscousVariables.SetElementGeometry(domain_geometry);
    ViscousVariables.SetShapeFunctionsValues(shape_functions);

    //1.- Material constants
    ViscousVariables.Mu = MaterialProperties[DYNAMIC_VISCOSITY];
    ViscousVariables.BulkModulus = MaterialProperties[BULK_MODULUS];

    ViscousVariables.DeltaTime = current_process_info[DELTA_TIME];

    //3.-DeformationGradient Tensor 3D
    ViscousVariables.DeformationGradientF=DeformationGradientF;
    ViscousVariables.DeformationGradientF = Transform2DTo3D(ViscousVariables.DeformationGradientF);

    //4.-Determinant of the Total Deformation Gradient
    ViscousVariables.DeterminantF=DeterminantF;


    //5.-Left Cauchy Green tensor b: (stored in the CauchyGreenMatrix)
    ViscousVariables.CauchyGreenMatrix.resize(3,3,false);
    noalias(ViscousVariables.CauchyGreenMatrix) = prod(ViscousVariables.DeformationGradientF,trans(ViscousVariables.DeformationGradientF));


    //6.- Updating ViscousVariables.DeformationRate with an extimate of the deformation rate
    this->CalculateDeformationRate(ViscousVariables);

    //7.-Almansi Strain:
    if (Options.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN))
    {
        // e= 0.5*(1-invbT*invb)
        this->CalculateAlmansiStrain(ViscousVariables.CauchyGreenMatrix,StrainVector);
    }

    //8.-Calculate Total Cauchy stress
    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
    {
        this->CalculateStress( ViscousVariables, StressMeasure_Cauchy, StressVector );

    }

    //9.-Calculate Constitutive Matrix related to Total Cauchy stress
    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
        this->CalculateConstitutiveMatrix ( ViscousVariables, ConstitutiveMatrix );
    }

}

//************************************************************************************
//************************************************************************************

void DispNewtonianFluid3DLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);

    const double& DeterminantF=rValues.GetDeterminantF();
    Vector& StressVector=rValues.GetStressVector();
    Matrix& ConstitutiveMatrix=rValues.GetConstitutiveMatrix();

    //Set to Kirchhoff Stress:
    StressVector=StressVector*DeterminantF;
    ConstitutiveMatrix=ConstitutiveMatrix*DeterminantF;
}

//***********************************UPDATE*******************************************
//************************************************************************************


void DispNewtonianFluid3DLaw::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponseKirchhoff(rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    UpdateInternalVariables( rValues );
}

//************************************************************************************
//************************************************************************************

void DispNewtonianFluid3DLaw::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponseCauchy(rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    UpdateInternalVariables( rValues );
}


//************************************************************************************
//************************************************************************************

void DispNewtonianFluid3DLaw::UpdateInternalVariables(Parameters& rValues)
{
    const Matrix& DeformationGradientF    = rValues.GetDeformationGradientF();
    const double& DeterminantF            = rValues.GetDeterminantF();

    Matrix DeformationGradientF0          = DeformationGradientF;
    DeformationGradientF0 = Transform2DTo3D(DeformationGradientF0);
    MathUtils<double>::InvertMatrix( DeformationGradientF0, this->mInverseDeformationGradientF0, mDeterminantF0);
    mDeterminantF0 = DeterminantF; //special treatment of the determinant
}


//***********************COMPUTE TOTAL STRAIN VECTOR**********************************
//************************************************************************************

void DispNewtonianFluid3DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
						Vector& rStrainVector )
{
    // e = 0.5*(1-invFT*invF) or e = 0.5*(1-inv(b))
    //Calculating the inverse of the jacobian
    Matrix InverseLeftCauchyGreen = ZeroMatrix( 3, 3 );
    double det_b=0;
    MathUtils<double>::InvertMatrix( rLeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    rStrainVector[0] = 0.5 * (1.00 - InverseLeftCauchyGreen(0, 0));
    rStrainVector[1] = 0.5 * (1.00 - InverseLeftCauchyGreen(1, 1));
    rStrainVector[2] = 0.5 * (1.00 - InverseLeftCauchyGreen(2, 2));
    rStrainVector[3] = -InverseLeftCauchyGreen(0, 1); // xy
    rStrainVector[4] = -InverseLeftCauchyGreen(1, 2); // yz
    rStrainVector[5] = -InverseLeftCauchyGreen(0, 2); // xz

}


//******************************* COMPUTE DEFORMATION RATE  ************************
//************************************************************************************

void DispNewtonianFluid3DLaw::CalculateDeformationRate(MaterialResponseVariables& rViscousVariables)
{

    // Computation of a first order extimate of the deformation rate
    // as done in Simo&Hughes, Computational Inelasticity, 1998

    // Incremental deformation rate
    Matrix f = prod(rViscousVariables.DeformationGradientF, mInverseDeformationGradientF0);

    // incremental LeftCauchyGreen
    Matrix b = prod(f, trans(f));

    // computation of inverse of incremental LeftCauchyGreen
    Matrix inv_b;
    double aux_det;
    MathUtils<double>::InvertMatrix(b, inv_b, aux_det);

    rViscousVariables.DeformationRate.resize(3, 3, false);
    noalias(rViscousVariables.DeformationRate) = 0.5 / rViscousVariables.DeltaTime * (rViscousVariables.Identity - inv_b);

}


//***************************** COMPUTE VOLUMETRIC PRESSURE  *************************
//************************************************************************************


double& DispNewtonianFluid3DLaw::CalculateVolumetricPressure(const MaterialResponseVariables& rViscousVariables,
    double& rPressure)
{

    rPressure = -1.0 * rViscousVariables.BulkModulus * (pow(rViscousVariables.DeterminantF, -1.0) - 1.0);

    return rPressure;
}


//************************* COMPUTE VOLUMETRIC PRESSURE FACTORS***********************
//************************************************************************************

Vector& DispNewtonianFluid3DLaw::CalculateVolumetricPressureFactors(const MaterialResponseVariables& rViscousVariables,
    Vector& rFactors)

{

    if (rFactors.size() != 3) rFactors.resize(3, false);

    double Pressure = 0;
    this->CalculateVolumetricPressure(rViscousVariables, Pressure);

    double Pressure_tilde = rViscousVariables.BulkModulus;

    rFactors[0] = Pressure_tilde;
    rFactors[1] = (2.0 * Pressure);
    rFactors[2] = 1.0;

    return rFactors;
}

//***********************COMPUTE TOTAL STRESS*****************************************
//************************************************************************************

void DispNewtonianFluid3DLaw::CalculateStress(MaterialResponseVariables& rViscousVariables, StressMeasure rStressMeasure, Vector& rStressVector)
{

    rViscousVariables.StressMatrix.resize(3,3,false);

    if (rStressMeasure==StressMeasure_Cauchy)
    {
        double Pressure = 0;
        Pressure = this->CalculateVolumetricPressure(rViscousVariables, Pressure);

        rViscousVariables.StressMatrix = Pressure * IdentityMatrix(3);

        Matrix DeviatoricPart;
        this->CalculateDeviatoricPart(rViscousVariables.DeformationRate, DeviatoricPart);

        rViscousVariables.StressMatrix += 2.0 * rViscousVariables.Mu * DeviatoricPart;

    }

    rStressVector=MathUtils<double>::StressTensorToVector(rViscousVariables.StressMatrix, rStressVector.size());
}

//******************************* COMPUTE VOLUMETRIC STRESS  *************************
//************************************************************************************

void DispNewtonianFluid3DLaw::CalculateVolumetricStress(const MaterialResponseVariables& rViscousVariables,
    Vector& rVolStressVector)
{

    Matrix VolStressMatrix = ZeroMatrix(3, 3);

    double Pressure = 0;
    Pressure = this->CalculateVolumetricPressure(rViscousVariables, Pressure);

    VolStressMatrix = Pressure * rViscousVariables.Identity;
    rVolStressVector = MathUtils<double>::StressTensorToVector(VolStressMatrix, rVolStressVector.size());

}


//******************************* COMPUTE ISOCHORIC STRESS  **************************
//************************************************************************************
void DispNewtonianFluid3DLaw::CalculateIsochoricStress(const MaterialResponseVariables& rViscousVariables,
    StressMeasure rStressMeasure,
    Vector& rIsoStressVector)
{

    Matrix IsoStressMatrix(3, 3);

    if (rStressMeasure == StressMeasure_Cauchy)
    {
        IsoStressMatrix = 2.0 * rViscousVariables.Mu * rViscousVariables.DeformationRate;
    }

    rIsoStressVector = MathUtils<double>::StressTensorToVector(IsoStressMatrix, rIsoStressVector.size());

}

//***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
//************************************************************************************

void DispNewtonianFluid3DLaw::CalculateConstitutiveMatrix (const MaterialResponseVariables& rViscousVariables, Matrix& rConstitutiveMatrix)
{

    rConstitutiveMatrix.clear();

    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
            rConstitutiveMatrix( i, j ) = ConstitutiveComponent(rConstitutiveMatrix( i, j ), rViscousVariables,
                                          this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }
    }

}

//***********************COMPUTE TENSOR COMPONENTS************************************
//************************************************************************************

double& DispNewtonianFluid3DLaw::ConstitutiveComponent (double& rCabcd, const MaterialResponseVariables& rViscousVariables,
    const unsigned int& a, const unsigned int& b, const unsigned int& c, const unsigned int& d)
{

    double IdotI = rViscousVariables.Identity(a, b) * rViscousVariables.Identity(c, d);
    double Isym = (rViscousVariables.Identity(a, c) * rViscousVariables.Identity(b, d) +
        rViscousVariables.Identity(a, d) * rViscousVariables.Identity(b, c)) / 2.0;
    Matrix I = rViscousVariables.Identity;

    Matrix f = prod(rViscousVariables.DeformationGradientF, mInverseDeformationGradientF0);
    Matrix FingerMatrix = prod(f, trans(f));
    Matrix inv_b;
    double aux_det;
    MathUtils<double>::InvertMatrix(FingerMatrix, inv_b, aux_det);
    double trace_inv_b = inv_b(0, 0) + inv_b(1, 1) + inv_b(2, 2);

    // Volumetric part
    Vector Factors(3);
    noalias(Factors) = ZeroVector(3);
    Factors = this->CalculateVolumetricPressureFactors(rViscousVariables, Factors);

    rCabcd = Factors[0] * IdotI;
    rCabcd -= Factors[1] * Isym;
    rCabcd *= Factors[2];


    // Deviatoric part
    rCabcd += rViscousVariables.Mu / rViscousVariables.DeltaTime * (I(a,c)*inv_b(b,d)+
        I(a, d) * inv_b(b, c)+ I(b, d) * inv_b(a, c)+ I(b, c) * inv_b(a, d)-
        I(c, d) * inv_b(a, b));

    rCabcd -= rViscousVariables.Mu / rViscousVariables.DeltaTime * (2.0 / 3.0 * trace_inv_b * Isym
        - 1.0 / 3.0 * trace_inv_b * I(a, b) * I(c, d) + 2.0 / 3.0 * I(a, b) * inv_b(c, d));

    return rCabcd;

}


//***********************COMPUTE VOLUMETRIC CONSTITUTIVE MATRIX***********************
//************************************************************************************

void DispNewtonianFluid3DLaw::CalculateVolumetricConstitutiveMatrix(const MaterialResponseVariables& rViscousVariables,
    Matrix& rConstitutiveMatrix)
{
    rConstitutiveMatrix.clear();

    Vector Factors(3);
    noalias(Factors) = ZeroVector(3);
    Factors = this->CalculateVolumetricPressureFactors(rViscousVariables, Factors);

    for (unsigned int i = 0; i < 6; i++)
    {
        for (unsigned int j = 0; j < 6; j++)
        {
            rConstitutiveMatrix(i, j) = VolumetricConstitutiveComponent(rConstitutiveMatrix(i, j), rViscousVariables, Factors,
                this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }
    }

}

//********************CONSTITUTIVE MATRIX VOLUMETRIC COMPONENT************************
//************************************************************************************


double& DispNewtonianFluid3DLaw::VolumetricConstitutiveComponent(double& rCabcd,
    const MaterialResponseVariables& rViscousVariables,
    const Vector& rFactors,
    const unsigned int& a, const unsigned int& b,
    const unsigned int& c, const unsigned int& d)
{
    double IdotI = rViscousVariables.Identity(a, b) * rViscousVariables.Identity(c, d);
    double Isym = (rViscousVariables.Identity(a, c) * rViscousVariables.Identity(b, d) +
        rViscousVariables.Identity(a, d) * rViscousVariables.Identity(b, c)) / 2.0;

    rCabcd = rFactors[0] * IdotI;
    rCabcd -= rFactors[1] * Isym;
    rCabcd *= rFactors[2];

    return rCabcd;
}

//***********************COMPUTE ISOCHORIC CONSTITUTIVE MATRIX************************
//************************************************************************************

void DispNewtonianFluid3DLaw::CalculateIsochoricConstitutiveMatrix(const MaterialResponseVariables& rViscousVariables,
    const Matrix& rIsoStressMatrix,
    Matrix& rConstitutiveMatrix)
{
    rConstitutiveMatrix.clear();

    for (unsigned int i = 0; i < 6; i++)
    {
        for (unsigned int j = 0; j < 6; j++)
        {
            rConstitutiveMatrix(i, j) = IsochoricConstitutiveComponent(rConstitutiveMatrix(i, j), rViscousVariables, rIsoStressMatrix,
                this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }
    }
}


//********************CONSTITUTIVE MATRIX ISOCHORIC COMPONENT*************************
//************************************************************************************

double& DispNewtonianFluid3DLaw::IsochoricConstitutiveComponent(double& rCabcd,
    const MaterialResponseVariables& rViscousVariables,
    const Matrix& rIsoStressMatrix,
    const unsigned int& a, const unsigned int& b,
    const unsigned int& c, const unsigned int& d)
{
    double IdotI = rViscousVariables.Identity(a, b) * rViscousVariables.Identity(c, d);
    double Isym = (rViscousVariables.Identity(a, c) * rViscousVariables.Identity(b, d) +
        rViscousVariables.Identity(a, d) * rViscousVariables.Identity(b, c)) / 2.0;
    Matrix I = rViscousVariables.Identity;

    Matrix f = prod(rViscousVariables.DeformationGradientF, mInverseDeformationGradientF0);
    Matrix FingerMatrix = prod(f, trans(f));
    Matrix inv_b;
    double aux_det;
    MathUtils<double>::InvertMatrix(FingerMatrix, inv_b, aux_det);

    rCabcd = rViscousVariables.Mu / rViscousVariables.DeltaTime * (I(a, c) * inv_b(b, d) +
        I(a, d) * inv_b(b, c) + I(b, d) * inv_b(a, c) + I(b, c) * inv_b(a, d) -
        I(c, d) * inv_b(a, b));

    rCabcd += rViscousVariables.Mu / rViscousVariables.DeltaTime * IdotI
        - 2.0 * rViscousVariables.Mu / rViscousVariables.DeltaTime * Isym;

    return rCabcd;
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

/**
 * Takes a matrix 2x2 and transforms it to a 3x3 adding a 3rd row and a 3rd column with a 1 in the diagonal
 * if the matrix passed is 3D is does nothing
 * if the matrix passed is bigger or smaller throws an error
 * @param rMatrix : usually the DeformationGradientF
 */
Matrix& DispNewtonianFluid3DLaw::Transform2DTo3D(Matrix& rMatrix)
{
    if (rMatrix.size1() == 2 && rMatrix.size2() == 2)
    {
        const BoundedMatrix<double, 2, 2> temp_matrix = rMatrix;
        rMatrix.resize(3, 3, false);
        rMatrix = IdentityMatrix(3);

        rMatrix(0, 0) = temp_matrix(0, 0);
        rMatrix(1, 1) = temp_matrix(1, 1);

        rMatrix(0, 1) = temp_matrix(0, 1);
        rMatrix(1, 0) = temp_matrix(1, 0);
    }
    else if (rMatrix.size1() != 3 && rMatrix.size2() != 3)
    {
        KRATOS_ERROR << "Matrix Dimensions are not correct !" << std::endl;
    }

    return rMatrix;
}

/**
 * Compute the deviatoric part of a 3x3 Matrix
 * @param rMatrix : Matrix whose deviatoric part should be computed
 * @param rDevMatrix : Deviatoric part of rMatrix
 */
void DispNewtonianFluid3DLaw::CalculateDeviatoricPart(const Matrix& rMatrix,
    Matrix& rDevMatrix)
{
    rDevMatrix = rMatrix;

    double trace = (rMatrix(0, 0) + rMatrix(1, 1) + rMatrix(2, 2));

    rDevMatrix(0, 0) -= trace / 3.0;
    rDevMatrix(1, 1) -= trace / 3.0;
    rDevMatrix(2, 2) = -(rDevMatrix(0, 0) + rDevMatrix(1, 1));

}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void DispNewtonianFluid3DLaw::GetLawFeatures(Features& rFeatures)
{
    	//Set the type of law
	rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
	rFeatures.mOptions.Set( FINITE_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	//Set the strain size
	rFeatures.mStrainSize = GetStrainSize();

	//Set the spacedimension
	rFeatures.mSpaceDimension = WorkingSpaceDimension();

}

//*************************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW *****************
//************************************************************************************

bool DispNewtonianFluid3DLaw::CheckParameters(Parameters& rValues)
{
     return rValues.CheckAllParameters();
}

int DispNewtonianFluid3DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) const
{

    KRATOS_ERROR_IF(DENSITY.Key()==0 || rMaterialProperties[DENSITY]<=0.00)<<"DENSITY has Key zero or invalid value"<< std::endl;

    KRATOS_ERROR_IF(DYNAMIC_VISCOSITY.Key() == 0 || rMaterialProperties[DYNAMIC_VISCOSITY] < 0.00) << "DYNAMIC_VISCOSITY has Key zero or invalid value" << std::endl;

    KRATOS_ERROR_IF(BULK_MODULUS.Key() == 0 || rMaterialProperties[BULK_MODULUS] <= 0.00) << "BULK_MODULUS has Key zero or invalid value" << std::endl;

    return 0;
}

} // Namespace Kratos
