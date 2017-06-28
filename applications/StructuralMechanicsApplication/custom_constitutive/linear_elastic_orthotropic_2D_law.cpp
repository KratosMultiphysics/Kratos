// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Peter Wilson
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/linear_elastic_orthotropic_2D_law.hpp"

#include "structural_mechanics_application_variables.h"

namespace Kratos
{
	//******************************CONSTRUCTOR*********************************
	//**************************************************************************

	LinearElasticOrthotropic2DLaw::LinearElasticOrthotropic2DLaw()
		: ConstitutiveLaw()
	{
	}

	//******************************COPY CONSTRUCTOR****************************
	//**************************************************************************

	LinearElasticOrthotropic2DLaw::LinearElasticOrthotropic2DLaw
		(const LinearElasticOrthotropic2DLaw& rOther)
		: ConstitutiveLaw(rOther)
		, mInverseDeformationGradientF0(rOther.mInverseDeformationGradientF0)
		, mDeterminantF0(rOther.mDeterminantF0)
		, mStrainEnergy(rOther.mStrainEnergy)
	{
	}

	//********************************CLONE*************************************
	//**************************************************************************

	ConstitutiveLaw::Pointer LinearElasticOrthotropic2DLaw::Clone() const
	{
		LinearElasticOrthotropic2DLaw::Pointer p_clone
			(new LinearElasticOrthotropic2DLaw(*this));
		return p_clone;
	}

	//*******************************DESTRUCTOR*********************************
	//**************************************************************************

	LinearElasticOrthotropic2DLaw::~LinearElasticOrthotropic2DLaw()
	{
	}

	void LinearElasticOrthotropic2DLaw::testString()
	{
		std::cout << "Printing LinearElasticOrthotropic2DLaw test string" << std::endl;
	}

	//*******************************OPERATIONS FROM BASE CLASS***************************
	//************************************************************************************

	//***********************HAS : DOUBLE - VECTOR - MATRIX*******************************
	//************************************************************************************
	
	// pwdebug
	bool LinearElasticOrthotropic2DLaw::Has(const Variable<double>& rThisVariable)
	{
		return false;
	}

	bool LinearElasticOrthotropic2DLaw::Has(const Variable<Vector>& rThisVariable)
	{
		return false;
	}

	bool LinearElasticOrthotropic2DLaw::Has(const Variable<Matrix>& rThisVariable)
	{
		return false;
	}
	

	//***********************GET VALUE: DOUBLE - VECTOR - MATRIX**************************
	//************************************************************************************
	
	// pwdebug
	double& LinearElasticOrthotropic2DLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
	{
		if (rThisVariable == STRAIN_ENERGY)
		{
			rValue = mStrainEnergy;
		}
		else {
			rValue = 0;
		}


		return(rValue);
	}

	Vector& LinearElasticOrthotropic2DLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
	{
		return(rValue);
	}

	Matrix& LinearElasticOrthotropic2DLaw::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
	{
		return(rValue);
	}
	

	//***********************SET VALUE: DOUBLE - VECTOR - MATRIX**************************
	//************************************************************************************

	
	// pwdebug
	void LinearElasticOrthotropic2DLaw::SetValue(const Variable<double>& rThisVariable, const double& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{

		if (rThisVariable == DETERMINANT_F)
		{
			mDeterminantF0 = rValue;
		}
	}

	void LinearElasticOrthotropic2DLaw::SetValue(const Variable<Vector>& rThisVariable, const Vector& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{

	}

	void LinearElasticOrthotropic2DLaw::SetValue(const Variable<Matrix>& rThisVariable, const Matrix& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{

	}
	


	//************* STARTING - ENDING  METHODS
	//************************************************************************************
	//************************************************************************************

	
	// pwdebug
	void LinearElasticOrthotropic2DLaw::InitializeMaterial(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues)
	{
		mDeterminantF0 = 1;
		mInverseDeformationGradientF0 = identity_matrix<double>(3);
		mStrainEnergy = 0;

	}
	
	//************************************************************************************
	//************************************************************************************

	
	// pwdebug
	void LinearElasticOrthotropic2DLaw::InitializeSolutionStep(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry, //this is just to give the array of nodes
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{

	}
	
	//************************************************************************************
	//************************************************************************************

	
	// pwdebug
	void LinearElasticOrthotropic2DLaw::FinalizeSolutionStep(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry, //this is just to give the array of nodes
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{

	}
	


	//************* COMPUTING  METHODS
	//**************************************************************************
	//**************************************************************************








	//*****************************MATERIAL RESPONSES***************************
	//**************************************************************************

	void  LinearElasticOrthotropic2DLaw::CalculateMaterialResponsePK2(Parameters& rValues)
	{
		//1.- Lame constants
		//const double& YoungModulus = MaterialProperties[YOUNG_MODULUS];
		//const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];

		//a.-Check if the constitutive parameters are passed correctly to the law calculation
		//CheckParameters(rValues);

		//b.- Get Values to compute the constitutive law:
		Flags &Options = rValues.GetOptions();

		const Properties& MaterialProperties = rValues.GetMaterialProperties();

		Vector& StrainVector = rValues.GetStrainVector();
		Vector& StressVector = rValues.GetStressVector();

		//-----------------------------//

		if (Options.Is(ConstitutiveLaw::COMPUTE_STRAIN))
		{
			//only needed
			const Matrix& DeformationGradientF = rValues.GetDeformationGradientF();

			//4.-Right Cauchy Green
			Matrix RightCauchyGreen = prod(trans(DeformationGradientF), DeformationGradientF);

			//5.-Green-Lagrange Strain:

			//E= 0.5*(FT*F-1)
			this->CalculateGreenLagrangeStrain(RightCauchyGreen, StrainVector);
		}

		//7.-Calculate Total PK2 stress

		if (Options.Is(ConstitutiveLaw::COMPUTE_STRESS))
		{
			if (Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
			{
				Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
				this->CalculateLinearElasticMatrix(ConstitutiveMatrix, MaterialProperties);
				this->CalculateStress(StrainVector, ConstitutiveMatrix, StressVector);
			}
			else {
				Matrix ConstitutiveMatrix(StrainVector.size(), StrainVector.size());
				noalias(ConstitutiveMatrix) = ZeroMatrix(StrainVector.size(), StrainVector.size());

				this->CalculateLinearElasticMatrix(ConstitutiveMatrix, MaterialProperties);
				this->CalculateStress(StrainVector, ConstitutiveMatrix, StressVector);
			}
		}
		else if (Options.IsNot(ConstitutiveLaw::COMPUTE_STRESS) && Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
		{
			Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
			this->CalculateLinearElasticMatrix(ConstitutiveMatrix, MaterialProperties);
		}
	}

	

	//***********************************UPDATE*******************************************
	//************************************************************************************
	/*
	// pwdebug
	void LinearElasticOrthotropic2DLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
	{

		rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
		this->CalculateMaterialResponsePK2(rValues);
		rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

		UpdateInternalVariables(rValues);
	}
	*/

	//************************************************************************************
	//************************************************************************************
	/*
	// pwdebug
	void LinearElasticOrthotropic2DLaw::UpdateInternalVariables(Parameters& rValues)
	{
		const Matrix& DeformationGradientF = rValues.GetDeformationGradientF();
		const double& DeterminantF = rValues.GetDeterminantF();

		Matrix DeformationGradientF0 = DeformationGradientF;
		DeformationGradientF0 = Transform2DTo3D(DeformationGradientF0);
		MathUtils<double>::InvertMatrix(DeformationGradientF0, this->mInverseDeformationGradientF0, mDeterminantF0);
		mDeterminantF0 = DeterminantF; //special treatment of the determinant 
	}
	*/
	
	//***********************COMPUTE TOTAL STRAIN*****************************************
	//************************************************************************************

	void LinearElasticOrthotropic2DLaw::CalculateGreenLagrangeStrain(const Matrix & rRightCauchyGreen,
		Vector& rStrainVector)
	{
		// TAKEN FROM LinearElasticPlasticPlaneStrain2DLaw
		//E= 0.5*(FT*F-1)
		rStrainVector[0] = 0.5 * (rRightCauchyGreen(0, 0) - 1.00);
		rStrainVector[1] = 0.5 * (rRightCauchyGreen(1, 1) - 1.00);
		rStrainVector[2] = rRightCauchyGreen(0, 1);
	}

	//***********************COMPUTE TOTAL STRESS PK2*************************************
	//************************************************************************************

	void LinearElasticOrthotropic2DLaw::CalculateStress(const Vector & rStrainVector,
		const Matrix & rConstitutiveMatrix,
		Vector& rStressVector)
	{
		//1.-2nd Piola Kirchhoff StressVector increment
		if (rStressVector.size() != rStrainVector.size())
			rStressVector.resize(rStrainVector.size(), false);

		noalias(rStressVector) = prod(rConstitutiveMatrix, rStrainVector);
	}

	//***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
	//************************************************************************************

	void LinearElasticOrthotropic2DLaw::CalculateLinearElasticMatrix(Matrix& rConstitutiveMatrix,
		const Properties& rMaterialProperties)
	{
		double E1 = rMaterialProperties[YOUNG_MODULUS_X];
		double E2 = rMaterialProperties[YOUNG_MODULUS_Y];

		double G12 = rMaterialProperties[SHEAR_MODULUS_XY];

		//double G13 = G12;	// currently handled through "shell_cross_section.cpp"
		//double G23 = G12;	// currently handled through "shell_cross_section.cpp"

		double v12 = rMaterialProperties[POISSON_RATIO_XY];
		double v21 = v12*E2 / E1;

		double Q11 = E1 / (1.0 - v12*v21);
		double Q12 = v12*E2 / (1.0 - v12*v21);
		double Q22 = E2 / (1.0 - v12*v21);
		double Q66 = G12;
		//double Q44 = G23;
		//double Q55 = G13;

		double theta = 0.0;	// rotation currently handled through "shell_cross_section.cpp" variable iPlyAngle
		double c = cos(theta);
		double c2 = c*c;
		double c4 = c2 * c2;
		double s = sin(theta);
		double s2 = s*s;
		double s4 = s2*s2;

		rConstitutiveMatrix.clear();

		rConstitutiveMatrix(0, 0) = Q11*c4 + 2.0*(Q12 + 2.0*Q66)*s2*c2 + Q22*s4;				// Q11_hat
		rConstitutiveMatrix(0, 1) = (Q11 + Q22 - 4.0*Q66)*s2*c2 + Q12*(s4 + c4);				// Q12_hat
		rConstitutiveMatrix(0, 2) = (Q11 - Q12 - 2.0*Q66)*s*c2*c + (Q12 - Q22 + 2.0*Q66)*s*s2*c;// Q16_hat

		rConstitutiveMatrix(1, 0) = rConstitutiveMatrix(0, 1);
		rConstitutiveMatrix(1, 1) = Q11*s4 + 2.0 * (Q12 + 2.0*Q66)*s2*c2 + Q22*c4;				// Q22_hat
		rConstitutiveMatrix(1, 2) = (Q11 - Q12 - 2.0*Q66)*s2*s*c + (Q12 - Q22 + 2.0*Q66)*c2*c*s;// Q16_hat

		rConstitutiveMatrix(2, 0) = rConstitutiveMatrix(0, 2);
		rConstitutiveMatrix(2, 1) = rConstitutiveMatrix(1, 2);
		rConstitutiveMatrix(2, 2) = (Q11 + Q22 - 2.0*Q12 - 2.0*Q66)*s2*c2 + Q66*(s4 + c4);		//Q66_hat
	}

	//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
	//************************************************************************************

	
	void LinearElasticOrthotropic2DLaw::GetLawFeatures(Features& rFeatures)
	{
		//Set the type of law
		rFeatures.mOptions.Set(PLANE_STRESS_LAW);
		rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
		rFeatures.mOptions.Set(ANISOTROPIC);

		//Set strain measure required by the consitutive law
		rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
		rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

		//Set the strain size
		rFeatures.mStrainSize = GetStrainSize();

		//Set the spacedimension
		rFeatures.mSpaceDimension = WorkingSpaceDimension();
	}
	

	//******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW*************************
	//************************************************************************************

	bool LinearElasticOrthotropic2DLaw::CheckParameters(Parameters& rValues)
	{
		return rValues.CheckAllParameters();
	}

	int LinearElasticOrthotropic2DLaw::Check(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const ProcessInfo& rCurrentProcessInfo)
	{
		if (YOUNG_MODULUS_X.Key() == 0 || !rMaterialProperties.Has(YOUNG_MODULUS_X))
			KRATOS_THROW_ERROR(std::invalid_argument, "YOUNG_MODULUS_X has Key zero or invalid value ", "")

			if (YOUNG_MODULUS_Y.Key() == 0 || !rMaterialProperties.Has(YOUNG_MODULUS_Y))
				KRATOS_THROW_ERROR(std::invalid_argument, "YOUNG_MODULUS_Y has Key zero or invalid value ", "")

				if (POISSON_RATIO_XY.Key() == 0 || !rMaterialProperties.Has(POISSON_RATIO_XY))
					KRATOS_THROW_ERROR(std::invalid_argument, "POISSON_RATIO_XY has Key zero invalid value ", "")

					if (DENSITY.Key() == 0 || !rMaterialProperties.Has(DENSITY))
						KRATOS_THROW_ERROR(std::invalid_argument, "DENSITY has Key zero or invalid value ", "")

						return 0;
	}
} // Namespace Kratos