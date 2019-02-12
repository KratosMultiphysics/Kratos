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

#if !defined (KRATOS_TC_PLASTIC_DAMAGE_3D_LAW_H_INCLUDED)
#define  KRATOS_TC_PLASTIC_DAMAGE_3D_LAW_H_INCLUDED

// System includes

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
		 * @brief Initialize the material response in terms of 1st Piola-Kirchhoff stresses
		 * @see Parameters
		 */
		void InitializeMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues) override;

		/**
		 * @brief Initialize the material response in terms of 2nd Piola-Kirchhoff stresses
		 * @see Parameters
		 */
		void InitializeMaterialResponsePK2 (ConstitutiveLaw::Parameters& rValues) override;

		/**
		 * @brief Initialize the material response in terms of Kirchhoff stresses
		 * @see Parameters
		 */
		void InitializeMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues) override;

		/**
		 * @brief Initialize the material response in terms of Cauchy stresses
		 * @see Parameters
		 */
		void InitializeMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues) override;

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

		/**
		 * @brief This function provides the place to perform checks on the completeness of the input.
		 * @details It is designed to be called only once (or anyway, not often) typically at the beginning of the calculations, so to verify that nothing is missing from the input or that no common error is found.
		 * @param rMaterialProperties The properties of the material
		 * @param rElementGeometry The geometry of the element
		 * @param rCurrentProcessInfo The current process info instance
		 * @return 0 if OK, 1 otherwise
		 */
		int Check(
			const Properties& rMaterialProperties,
			const GeometryType& rElementGeometry,
			const ProcessInfo& rCurrentProcessInfo
			) override;

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

		/** @brief general material parameters
		 * @ m_elastic_uniaxial_compressive_strength...elastic limit uniaxial compressive strength (unequal to compressive strength)
		 * @ m_elastic_biaxial_compressive_strength...elastic limit biaxial compressive strength (unequal to compressive strength)
		 * */	
		double m_elastic_uniaxial_compressive_strength, m_tensile_strength, m_elastic_biaxial_compressive_strength, m_E, m_nu, m_Gf;
		
		Matrix m_D0;
		
		Vector m_elastic_strain;
		Vector m_plastic_strain;
		
		// rate intensity of plastic strain
		double m_beta;

		// material property which accounts for the increase of compressive strength due to biaxial compression
		double m_K;
		
		/** @brief distinction between three proposals for the damage surface
		 * @detail 1=ORIGINAL, 2=COMPDYN/proposal A, 3=HOMOGENEOUS/proposal B
		 */
		int usedEquivalentEffectiveStressDefinition;

		/** @brief damage variables
		 * @ m_damage_compression...negative damage variable
		 * m_damage_tension...positive damage variable
		 * m_initial_damage_threshold_compression...initial negative damage threshold
		 * m_initial_damage_threshold_tension...initial positive damage threshold
		 * m_damage_threshold_compression...current negative damage threshold
		 * m_damage_threshold_tension...current positive damage threshold
		 * m_damage_threshold_compression1...updated negative damage threshold
		 * m_damage_threshold_tension1...updated positive damage threshold
		 * m_compression_parameter_A...compression parameter A
		 * m_compression_parameter_B...compression parameter B
		 * m_tension_parameter_A...tension parameter A
		 */
		double m_damage_compression, m_damage_tension, m_initial_damage_threshold_compression, m_initial_damage_threshold_tension, 
		m_damage_threshold_compression, m_damage_threshold_tension, m_damage_threshold_compression1, m_damage_threshold_tension1, 
		m_compression_parameter_A, m_compression_parameter_B, m_tension_parameter_A;

		// @brief shear retention factor
		double m_SRF12, m_SRF13, m_SRF23;
		// @brief reference value of the strain used for the evolution law of the shear retention factor
		double m_strain_ref;

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
		void CalculateElasticityMatrix(
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
			double& rEquEffStressCompression,
			double& rEquEffStressTension);

		void DamageCriterion(
			const double& rEquEffStressCompression,
			const double& rEquEffStressTension,
			const double& rtolerance);

		/** @brief method computes the negative damage variable m_damage_compression
		 * @param rd_n...dummy for damage variable
		 */
		void ComputeDamageCompression();

		/** @brief method computes the positive damage variable m_damage_tension
		 * @param rd_p...dummy for damage variable
		 */
		void ComputeDamageTension();

		/** @brief method computes the three shear retention factors based on an evolution law
		 * @details see A scalar damage model with a shear retention factor for the
		 * analysis of reinforced concrete structures: theory and
		 * validation - Scotta(2001)
		 */
		void ComputeSRF(const Vector& rStrainVector);
		///@}
		///@name Private  Access
		///@{

		///@}
		///@name Private Inquiry
		///@{

		///@}
		///@name Un accessible methods
		///@{

		// Serialization

		friend class Serializer;

		void save(Serializer& rSerializer) const override
		{
			KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
			//rSerializer.save("PrevStressVector", mPrevStressVector);
			//rSerializer.save("PrevInelasticStrainVector", mPrevInelasticStrainVector);
			//rSerializer.save("NonConvPrevStressVector", mNonConvPrevStressVector);
			//rSerializer.save("NonConvPrevInelasticStrainVector", mNonConvPrevInelasticStrainVector);
		}

		void load(Serializer& rSerializer) override
		{
			KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
			//rSerializer.load("PrevStressVector", mPrevStressVector);
			//rSerializer.load("PrevInelasticStrainVector", mPrevInelasticStrainVector);
			//rSerializer.load("NonConvPrevStressVector", mNonConvPrevStressVector);
			//rSerializer.load("NonConvPrevInelasticStrainVector", mNonConvPrevInelasticStrainVector);
		}

		///@}

	}; // Class TCPlasticDamage3DLaw
} // namespace Kratos
#endif