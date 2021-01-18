// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Peter Wilson
//       Contact:    A.Winterstein [at] tum.de
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
	* Defines a linear elastic orthotropic constitutive law in 2D (Plane Stress)
	* This material law is defined by the parameters:
	* 1) ELASTIC MODULI (E1, E2, G12)
	* 2) POISSON RATIOS (v12, v21)
	* Valid for small strains.
	*/

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
		ConstitutiveLaw::Pointer Clone() const override;

		/**
		* Copy constructor.
		*/
		LinearElasticOrthotropic2DLaw(const LinearElasticOrthotropic2DLaw& rOther);

		/**
		* Destructor.
		*/
		~LinearElasticOrthotropic2DLaw() override;

		/**
		* Operators
		*/

		/**
		* Operations needed by the base class:
		*/
		/**
		* Voigt tensor size:
		*/
		SizeType GetStrainSize() override
		{
			return 3;
		};

		/**
		* This function is designed to be called once to check compatibility with element
		* @param rFeatures
		*/
		void GetLawFeatures(Features& rFeatures) override;

		/**
		* Computes the material response:
		* PK2 stresses and algorithmic ConstitutiveMatrix
		* @param rValues
		* @see   Parameters
		*/
		void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters & rValues) override;

		/**
		 * returns the value of a specified variable
		 * @param rThisVariable the variable to be returned
		 * @param rValue a reference to the returned value
		 * @param rValue output: the value of the specified variable
		 */
		bool& GetValue(const Variable<bool>& rThisVariable, bool& rValue) override;

		/**
		* This function is designed to be called once to perform all the checks needed
		* on the input provided. Checks can be "expensive" as the function is designed
		* to catch user's errors.
		* @param rMaterialProperties
		* @param rElementGeometry
		* @param rCurrentProcessInfo
		* @return
		*/
		int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;

	protected:

		///@name Protected static Member Variables
		///@{
		///@}
		///@name Protected member Variables
		///@{
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
		bool CheckParameters(ConstitutiveLaw::Parameters& rValues);

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

		void save(Serializer& rSerializer) const override
		{
			KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
		}

		void load(Serializer& rSerializer) override
		{
			KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
		}

		///@}
	}; // Class LinearElasticOrthotropic2DLaw
}  // namespace Kratos.
#endif // KRATOS_LINEAR_ELASTIC_ORTHOTROPIC_2D_LAW_H_INCLUDED  defined