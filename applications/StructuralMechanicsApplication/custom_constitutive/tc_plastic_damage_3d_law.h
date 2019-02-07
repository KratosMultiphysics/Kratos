// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Michael Loibl
//					 Tobias Teschemacher
//					 Riccardo Rossi
//
//  Based on work of Tesser and Talledo (University of Padua)		 

// Written: Tesser, Talledo

//  Send strains in following format :
// 
//     strain_vec = {   eps_00
//                      eps_11
//                    2 eps_01   }   <--- note the 2

#if !defined (KRATOS_TC_PLASTIC_DAMAGE_3D_LAW_H_INCLUDED)
#define  KRATOS_TC_PLASTIC_DAMAGE_3D_LAW_H_INCLUDED

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos
{
	class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) TCPlasticDamage3DLaw
		: public ConstitutiveLaw
	{
	public:

		///@name Type Definitions
		///@{

		/// Definition of the base class
		typedef ConstitutiveLaw BaseType;

		/// The index definition
		typedef std::size_t IndexType;

		/// The size definition
		typedef std::size_t SizeType;

		/// Counted pointer of GenericYieldSurface
		KRATOS_CLASS_POINTER_DEFINITION(TCPlasticDamage3DLaw);

		///@}
		///@name Life Cycle
		///@{

		/**
		* Default constructor.
		*/
		TCPlasticDamage3DLaw()
		{
		}

		/**
		* Clone.
		*/
		ConstitutiveLaw::Pointer Clone() const override
		{
			return Kratos::make_shared<TCPlasticDamage3DLaw>(*this);
		}

		/**
		* Copy constructor.
		*/
		TCPlasticDamage3DLaw(const TCPlasticDamage3DLaw& rOther)
			: ConstitutiveLaw(rOther)
		{
		}

		/**
		 * Destructor.
		 */
		~TCPlasticDamage3DLaw() override
		{
		}

		///@}
		///@name Operators
		///@{

		///@}
		///@name Operations
		///@{

		/**
		 * @brief Dimension of the law:
		 */

		SizeType WorkingSpaceDimension() override
		{
			return 3;
		};

		/**
		 * @brief Voigt tensor size:
		 */
		SizeType GetStrainSize() override
		{
			return 6;
		};

		/**
			* returns the value of a specified variable
			* @param rThisVariable the variable to be returned
			* @param rValue a reference to the returned value
			* @param rValue output: the value of the specified variable
			*/
		virtual double& GetValue(
			const Variable<double>& rThisVariable,
			double& rValue) override;

		/**
			* @brief Returns the value of a specified variable (Vector)
			* @param rThisVariable the variable to be returned
			* @param rValue a reference to the returned value
			* @return rValue output: the value of the specified variable
			*/
		virtual Vector& GetValue(
			const Variable<Vector>& rThisVariable,
			Vector& rValue) override;

		/**
			* @brief Returns the value of a specified variable (Matrix)
			* @param rThisVariable the variable to be returned
			* @return rValue output: the value of the specified variable
		 */
		virtual Matrix& GetValue(
			const Variable<Matrix>& rThisVariable,
			Matrix& rValue) override;

		void InitializeMaterial(
			const Properties& rMaterialProperties,
			const GeometryType& rElementGeometry,
			const Vector& rShapeFunctionsValues);
		
		/**
		 * Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
		 * @see Parameters
		 */
		void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

		/**
		 * Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
		 * @see Parameters
		 */
		void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

		/**
		 * Computes the material response in terms of Kirchhoff stresses and constitutive tensor
		 * @see Parameters
		 */
		void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

		/**
		 * Computes the material response in terms of Cauchy stresses and constitutive tensor
		 * @see Parameters
		 */
		void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

		/**
		 * to be called at the end of each solution step
		 * (e.g. from Element::FinalizeSolutionStep)
		 * @param rMaterialProperties the Properties instance of the current element
		 * @param rElementGeometry the geometry of the current element
		 * @param rShapeFunctionsValues the shape functions values in the current integration point
		 * @param the current ProcessInfo instance
		 */
		//void FinalizeSolutionStep(
		//    const Properties& rMaterialProperties,
		//    const GeometryType& rElementGeometry,
		//    const Vector& rShapeFunctionsValues,
		//    const ProcessInfo& rCurrentProcessInfo) override;

		/**
		 * Finalize the material response in terms of 1st Piola-Kirchhoff stresses
		 * @see Parameters
		 */
		void FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

		/**
		 * Finalize the material response in terms of 2nd Piola-Kirchhoff stresses
		 * @see Parameters
		 */
		void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

		/**
		 * Finalize the material response in terms of Kirchhoff stresses
		 * @see Parameters
		 */
		void FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

		/**
		 * Finalize the material response in terms of Cauchy stresses
		 * @see Parameters
		 */
		void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

		/**
		 * @brief Calculates the value of a specified variable (Matrix)
		 * @param rParameterValues the needed parameters for the CL calculation
		 * @param rThisVariable the variable to be returned
		 * @param rValue a reference to the returned value
		 * @param rValue output: the value of the specified variable
		 */
		//Matrix& CalculateValue(
		//    ConstitutiveLaw::Parameters& rParameterValues,
		//    const Variable<Matrix>& rThisVariable,
		//    Matrix& rValue
		//    ) override;

		///@}
		///@name Access
		///@{

		///@}
		///@name Inquiry
		///@{

		///@}
		///@name Input and output
		///@{

		///@}
		///@name Friends
		///@{

		///@}


	protected:
	private:
		///@name Static Member Variables
		///@{

		///@}
		///@name Member Variables
		///@{
		double m_f_01cc, m_f_ct, m_f_02cc, m_E, m_nu, m_Gf;
		
		Matrix m_D0;
		
		Vector m_elastic_strain;
		Vector m_plastic_strain;
		
		double m_beta;

		// material property which accounts for the increase of compressive strength due to biaxial compression
		double m_K;
		
		/** @brief distinction between three proposals for the damage surface
		 * @detail 1=ORIGINAL, 2=COMPDYN/proposal A, 3=HOMOGENEOUS/proposal B
		 */
		int usedEquivalentTensionDefinition;

		/** @brief damage variables
		 * @ m_d_n...negative damage variable
		 * m_d_p...positive damage variable
		 * m_r_0n...initial negative damage threshold
		 * m_r_0p...initial positive damage threshold
		 * m_r_n...current negative damage threshold
		 * m_r_p...current positive damage threshold
		 * m_r_n1...updated negative damage threshold
		 * m_r_p1...updated positive damage threshold
		 * m_A_n...compression parameter A
		 * m_B_n...compression parameter B
		 * m_A_p...tension parameter A
		 */
		double m_d_n, m_d_p, m_r_0n, m_r_0p, m_r_n, m_r_p, m_r_n1, m_r_p1, m_A_n, m_B_n, m_A_p;

		// @brief shear retention factor
		double m_SRF12, m_SRF13, m_SRF23;
		// @brief reference value of the strain used for the evolution law of the shear retention factor
		double m_strain_ref;

		/** variables not used so far (ML)
		  Matrix m_Di;

		  double m_compressive_strength_plastic;
		  double m_compressive_strength_elastic;

		  Matrix m_eigen_vectors;
		  Vector m_eigen_values;

		  double m_gamma_C;

		  int model;

		//// Converged values
		//Vector mPrevStressVector = ZeroVector(6);
		//Vector mPrevInelasticStrainVector = ZeroVector(6);

		//// Non Converged values
		//Vector mNonConvPrevStressVector = ZeroVector(6);
		//Vector mNonConvPrevInelasticStrainVector = ZeroVector(6);

		*/

		///@}
		///@name Private Operators
		///@{

		///@}
		///@name Private Operations
		///@{

		/**
		 * @brief This method is the core of the constitutive law implementation which orders the single methods
		 */
		void CalculateMaterialResponseInternal(
			const Vector& rStrainVector,
			Vector& rStressVector,
			Matrix& rConstitutiveLaw);

		/**
			* @brief This method calculates the 3D Elasticity Matrix
			* @details voigt notation, 
			* note that the shear strains are considered as 2*strain_ij in the strain vector
		 */ 
		void TCPlasticDamage3DLaw::CalculateElasticityMatrix(
    		Matrix& rElasticityMatrix);

		/**
        	* @brief This method performs Spectral Decomposition of the Stress Vector/Tensor
			* @details see "An energy-Equivalent" d+/d- Damage model with Enhanced
			* Microcrack Closure/Reopening Capabilities for Cohesive-Frictional
			* Materials" - M. Cervera and C. Tesei.
			* @param rStressVector The Stress Vector
			* @param rStressVectorTension The Tension Part of the Stress Vector
			* @param rStressVectorCompression The Compression Part of the Stress Vector
			* @param PMatrixTension The Tensile P Matrix
			* @param PMatrixCompression The Compressive P Matrix
         */
        void SpectralDecomposition(
			const Vector& rStressVector,
			Vector& rStressVectorTension,
			Vector& rStressVectorCompression,
			Vector& rStressEigenvalues,
			Matrix& rPMatrixTension,
			Matrix& rPMatrixCompression);

		void ComputeTau(
			const Vector& rStressEigenvalues,
			double& rtau_n,
			double& rtau_p);

		void DamageCriterion(
			const double& rtau_n,
			const double& rtau_p,
			const double& rtolerance);

		/** @brief method computes the negative damage variable m_d_n
		 * @param rd_n...dummy for damage variable
		 */
		void ComputeDamageCompression();

		/** @brief method computes the positive damage variable m_d_p
		 * @param rd_p...dummy for damage variable
		 */
		void ComputeDamageTension();

		/** @brief method computes the three shear retention factors based on an evolution law
		 * @details see A scalar damage model with a shear retention factor for the
		 * analysis of reinforced concrete structures: theory and
		 * validation - Scotta(2001)
		 */
		void ComputeSRF(const Vector& rStrainVector);

		/** @brief method updates stiffness matrix
		 * @param rConstitutiveLaw...stiffness matrix
		 */
		void ComputeStiffnessMatrix(Matrix& rConstitutiveLaw,
			const Matrix& rPMatrixTension,
			const Matrix& rPMatrixCompression);
	};
}

#endif