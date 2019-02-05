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
		double m_f_01cc;
		double m_f_ct;
		double m_f_02cc;

		Vector m_elastic_strain;
		Vector m_plastic_strain;

		Matrix m_D0;

		double m_E;
		double m_nu;

		double m_K;
		std::string usedEquivalentTensionDefinition;
		std::string COMPDYN;
		std::string ORIGINAL;
		std::string HOMOGENEOUS;
		double d_n, d_p, r_0n, r_0p, r_n, r_p, r_n1, r_p1;

		/** variables not used so far (ML)
		  Matrix m_Di;

		  double m_compression_parameter_A;
		  double m_compression_parameter_B;

		  double m_tension_parameter_A;

		  double m_compressive_strength_plastic;
		  double m_compressive_strength_elastic;

		  double m_beta;
		  double m_Gf_t;
		  double m_Gf_c;

		  double m_K;

		  Matrix m_eigen_vectors;
		  Vector m_eigen_values;

		  double m_treshold_compression_initial;
		  double m_treshold_tension_initial;

		  double m_treshold_compression;
		  double m_treshold_tension;

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
			Matrix& PMatrixTension,
			Matrix& PMatrixCompression);

		void ComputeTau(
			const Vector& rStressEigenvalues,
			double& rtau_n,
			double& rtau_p);

		void DamageCriterion(
			const double& rtau_n,
			const double& rtau_p,
			const double& rtolerance);
	};
}

#endif