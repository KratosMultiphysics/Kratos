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
		 * Computes the material response in terms of Cauchy stresses and constitutive tensor
		 * @see Parameters
		 */
		void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;


	protected:
	private:
		///@name Static Member Variables
		///@{

		///@}
		///@name Member Variables
		///@{
		static constexpr SizeType VoigtSize = 6;

		double m_compressive_strength;
		double m_tensile_strength;

		Vector m_elastic_strain;
		Vector m_plastic_strain;

		Matrix m_D0;

		/** variables not used so far (ML)
		  Matrix m_Di;

		  double m_compression_parameter_A;
		  double m_compression_parameter_B;

		  double m_tension_parameter_A;

		  double m_compressive_strength_plastic;
		  double m_compressive_strength_elastic;

		  double m_rate_biaxial_uniaxial;

		  double m_beta;
		  double m_Gf_t;
		  double m_Gf_c;

		  double m_E;
		  double m_nu;

		  double m_K;

		  Matrix m_eigen_vectors;
		  Vector m_eigen_values;

		  double m_treshold_compression_initial;
		  double m_treshold_tension_initial;

		  double m_treshold_compression;
		  double m_treshold_tension;

		  double m_gamma_C;

		  double m_damage_t;
		  double m_damage_c;

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

		void CalculateMaterialResponseInternal(
			const Vector& rStrainVector,
			Vector& rStressVector,
			Matrix& rConstitutiveLaw);

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
            static void SpectralDecomposition(
			const Vector& rStressVector,
			Vector& rStressVectorTension,
			Vector& rStressVectorCompression,
			Matrix& PMatrixTension,
			Matrix& PMatrixCompression);
	};
}

#endif