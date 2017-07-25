// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Peter Wilson
//       Contact:    A.Winterstein@tum.de
//

#if !defined (KRATOS_LINEAR_ELASTIC_ORTHOTROPIC_2D_LAW_H_INCLUDED)
#define  KRATOS_LINEAR_ELASTIC_ORTHOTROPIC_2D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos
{
	/**
	* Defines a linear orthotropic constitutive law in 2D (Plane Stress)
	* This material law is defined by the parameters:
	* 1) YOUNG MODULUS
	* 2) POISSON RATIO
	* As there are no further parameters the functionality is valid
	* for small and large displacements elasticity.
	*/

	//class KRATOS_API(SOLID_MECHANICS_APPLICATION) LinearElasticOrthotropic2DLaw : public HyperElastic3DLaw
	class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) LinearElasticOrthotropic2DLaw : public ConstitutiveLaw
	{
	public:
		/**
		* Type Definitions
		*/
		typedef ProcessInfo      ProcessInfoType;
		typedef ConstitutiveLaw         BaseType;
		typedef std::size_t             SizeType;
		/**
		* Counted pointer of LinearElasticOrthotropic3DLaw
		*/

		KRATOS_CLASS_POINTER_DEFINITION(LinearElasticOrthotropic2DLaw);

		/**
		* Life Cycle
		*/

		/**
		* Default constructor.
		*/
		LinearElasticOrthotropic2DLaw();

		/**
		* Clone function (has to be implemented by any derived class)
		* @return a pointer to a new instance of this constitutive law
		*/
		ConstitutiveLaw::Pointer Clone() const;

		/**
		* Copy constructor.
		*/
		LinearElasticOrthotropic2DLaw(const LinearElasticOrthotropic2DLaw& rOther);

		/**
		* Assignment operator.
		*/

		//LinearElasticOrthotropic2DLaw& operator=(const LinearElasticOrthotropic2DLaw& rOther);

		/**
		* Destructor.
		*/
		virtual ~LinearElasticOrthotropic2DLaw();

		/**
		* Operators
		*/

		/**
		* Operations needed by the base class:
		*/
		/**
		* Voigt tensor size:
		*/
		SizeType GetStrainSize()
		{
			return 3;
		};

		void testString();

		/**
		* This function is designed to be called once to check compatibility with element
		* @param rFeatures
		*/
		void GetLawFeatures(Features& rFeatures); //update this

		bool Has(const Variable<double>& rThisVariable);
		bool Has(const Variable<Vector>& rThisVariable);
		bool Has(const Variable<Matrix>& rThisVariable);

		double& GetValue(const Variable<double>& rThisVariable, double& rValue);
		Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue);
		Matrix& GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue);

		void SetValue(const Variable<double>& rVariable,
			const double& rValue,
			const ProcessInfo& rCurrentProcessInfo);
		void SetValue(const Variable<Vector>& rThisVariable,
			const Vector& rValue,
			const ProcessInfo& rCurrentProcessInfo);
		void SetValue(const Variable<Matrix>& rThisVariable,
			const Matrix& rValue,
			const ProcessInfo& rCurrentProcessInfo);
		/**
		* Material parameters are inizialized
		*/
		void InitializeMaterial(const Properties& rMaterialProperties,
			const GeometryType& rElementGeometry,
			const Vector& rShapeFunctionsValues);

		void InitializeSolutionStep(const Properties& rMaterialProperties,
			const GeometryType& rElementGeometry, //this is just to give the array of nodes
			const Vector& rShapeFunctionsValues,
			const ProcessInfo& rCurrentProcessInfo);

		void FinalizeSolutionStep(const Properties& rMaterialProperties,
			const GeometryType& rElementGeometry, //this is just to give the array of nodes
			const Vector& rShapeFunctionsValues,
			const ProcessInfo& rCurrentProcessInfo);

		/**
		* Computes the material response:
		* PK2 stresses and algorithmic ConstitutiveMatrix
		* @param rValues
		* @see   Parameters
		*/
		void CalculateMaterialResponsePK2(Parameters & rValues);

		/**
		* Computes the material response:
		* Kirchhoff stresses and algorithmic ConstitutiveMatrix
		* @param rValues
		* @see   Parameters
		*/
		//void CalculateMaterialResponseKirchhoff(Parameters & rValues);

		/**
		* This function is designed to be called once to perform all the checks needed
		* on the input provided. Checks can be "expensive" as the function is designed
		* to catch user's errors.
		* @param rMaterialProperties
		* @param rElementGeometry
		* @param rCurrentProcessInfo
		* @return
		*/
		int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo);

		/**
		* Input and output
		*/
		/**
		* Turn back information as a string.
		*/
		//virtual String Info() const;
		/**
		* Print information about this object.
		*/
		//virtual void PrintInfo(std::ostream& rOStream) const;
		/**
		* Print object's data.
		*/
		//virtual void PrintData(std::ostream& rOStream) const;

	protected:

		///@name Protected static Member Variables
		///@{
		///@}
		///@name Protected member Variables
		///@{
		Matrix mInverseDeformationGradientF0;

		double mDeterminantF0;

		double mStrainEnergy;
		///@}
		///@name Protected Operators
		///@{
		///@}
		///@name Protected Operations
		///@{
		///@}
		/**
		* Calculates the GreenLagrange strains
		* @param rRightCauchyGreen
		* @param rStrainVector
		*/
		void CalculateGreenLagrangeStrain(const Matrix & rRightCauchyGreen,
			Vector& rStrainVector);

		/**
		* Calculates the stresses for given strain state
		* @param rStrainVector
		* @param rConstitutiveMatrix
		* @param rStressVector the stress vector corresponding to the deformation
		*/
		virtual void CalculateStress(const Vector &rStrainVector,
			const Matrix &rConstitutiveMatrix,
			Vector& rStressVector);

		/**
		* calculates the linear elastic constitutive matrix in terms of Young's modulus and
		* Poisson ratio
		* @param E the Young's modulus
		* @param NU the Poisson ratio
		* @return the linear elastic constitutive matrix
		*/

		virtual void CalculateLinearElasticMatrix(Matrix& rConstitutiveMatrix,
			const Properties& rMaterialProperties);

		/**
		* This function is designed to be called when before the material response
		* to check if all needed parameters for the constitutive are initialized
		* @param Parameters
		* @return
		*/
		bool CheckParameters(Parameters& rValues);

	private:

		///@name Static Member Variables
		///@{
		///@}
		///@name Member Variables
		///@{
		///@}
		///@name Private Operators
		///@{
		///@}
		///@name Private Operations
		///@{
		///@}

		///@}
		///@name Private  Access
		///@{
		///@}

		///@}
		///@name Serialization
		///@{
		friend class Serializer;

		virtual void save(Serializer& rSerializer) const
		{
			KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
				rSerializer.save("mInverseDeformationGradientF0", mInverseDeformationGradientF0);
			rSerializer.save("mDeterminantF0", mDeterminantF0);
			rSerializer.save("mStrainEnergy", mStrainEnergy);
		}

		virtual void load(Serializer& rSerializer)
		{
			KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
				rSerializer.load("mInverseDeformationGradientF0", mInverseDeformationGradientF0);
			rSerializer.load("mDeterminantF0", mDeterminantF0);
			rSerializer.load("mStrainEnergy", mStrainEnergy);
		}

		///@}
	}; // Class LinearElasticOrthotropic2DLaw
}  // namespace Kratos.
#endif // KRATOS_LINEAR_ELASTIC_ORTHOTROPIC_2D_LAW_H_INCLUDED  defined