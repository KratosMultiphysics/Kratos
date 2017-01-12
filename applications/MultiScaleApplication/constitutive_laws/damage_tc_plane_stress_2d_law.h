/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

/* *********************************************************
*
*   Last Modified by:    $Author:   Massimo Petracca$
*   Date:                $Date:     19-09-2014$
*   Revision:            $Revision: 1.0$
*
* ***********************************************************/

#if !defined(KRATOS_DAMAGE_TC_PLANE_STRESS_2D_LAW )
#define  KRATOS_DAMAGE_TC_PLANE_STRESS_2D_LAW

/* System includes */

/* External includes */

/* Project includes */
#include "includes/constitutive_law.h"

//#define DAM_TC_2D_IMPLEX
//#define DAM_TC_2D_INCREMENTAL_REGULARIZATION_V2

namespace Kratos
{

	/**
	* Base class of constitutive laws.
	*/
	class DamageTCPlaneStress2DLaw : public ConstitutiveLaw
	{

	public:

		KRATOS_CLASS_POINTER_DEFINITION(DamageTCPlaneStress2DLaw);

	public:

		/**
		* Constructor.
		*/
		DamageTCPlaneStress2DLaw();

		/**
		* Destructor.
		*/
		virtual ~DamageTCPlaneStress2DLaw(){};

		/**
		* Clone function (has to be implemented by any derived class)
		* @return a pointer to a new instance of this constitutive law
		* NOTE: implementation scheme:
		*      ConstitutiveLaw::Pointer p_clone(new DamageTCPlaneStress2DLaw());
		*      return p_clone;
		*/
		virtual ConstitutiveLaw::Pointer Clone() const;

		/**
		* @return the working space dimension of the current constitutive law
		* NOTE: this function HAS TO BE IMPLEMENTED by any derived class
		*/
		virtual SizeType WorkingSpaceDimension();

		/**
		* returns the size of the strain vector of the current constitutive law
		* NOTE: this function HAS TO BE IMPLEMENTED by any derived class
		*/
		virtual SizeType GetStrainSize();

		/**
		* returns whether this constitutive Law has specified variable
		* @param rThisVariable the variable to be checked for
		* @return true if the variable is defined in the constitutive law
		*/
		virtual bool Has(const Variable<double>& rThisVariable);

		/**
		* returns whether this constitutive Law has specified variable
		* @param rThisVariable the variable to be checked for
		* @return true if the variable is defined in the constitutive law
		*/
		virtual bool Has(const Variable<Vector>& rThisVariable);

		/**
		* returns whether this constitutive Law has specified variable
		* @param rThisVariable the variable to be checked for
		* @return true if the variable is defined in the constitutive law
		*/
		virtual bool Has(const Variable<Matrix>& rThisVariable);

		/**
		* returns whether this constitutive Law has specified variable
		* @param rThisVariable the variable to be checked for
		* @return true if the variable is defined in the constitutive law
		* NOTE: fixed size array of 3 doubles (e.g. for 2D stresses, plastic strains, ...)
		*/
		virtual bool Has(const Variable<array_1d<double, 3 > >& rThisVariable);

		/**
		* returns whether this constitutive Law has specified variable
		* @param rThisVariable the variable to be checked for
		* @return true if the variable is defined in the constitutive law
		* NOTE: fixed size array of 6 doubles (e.g. for stresses, plastic strains, ...)
		*/
		virtual bool Has(const Variable<array_1d<double, 6 > >& rThisVariable);

		/**
		* returns the value of a specified variable
		* @param rThisVariable the variable to be returned
		* @param rValue a reference to the returned value
		* @param rValue output: the value of the specified variable
		*/
		virtual double& GetValue(
			const Variable<double>& rThisVariable, 
			double& rValue);

		/**
		* returns the value of a specified variable
		* @param rThisVariable the variable to be returned
		* @param rValue a reference to the returned value
		* @return the value of the specified variable
		*/
		virtual Vector& GetValue(
			const Variable<Vector>& rThisVariable, 
			Vector& rValue);

		/**
		* returns the value of a specified variable
		* @param rThisVariable the variable to be returned
		* @return the value of the specified variable
		*/
		virtual Matrix& GetValue(
			const Variable<Matrix>& rThisVariable, 
			Matrix& rValue);

		/**
		* returns the value of a specified variable
		* @param rThisVariable the variable to be returned
		* @param rValue a reference to the returned value
		* @return the value of the specified variable
		*/
		virtual array_1d<double, 3 > & GetValue(
			const Variable<array_1d<double, 3 > >& rVariable,
			array_1d<double, 3 > & rValue);

		/**
		* returns the value of a specified variable
		* @param rThisVariable the variable to be returned
		* @param rValue a reference to the returned value
		* @return the value of the specified variable
		*/
		virtual array_1d<double, 6 > & GetValue(
			const Variable<array_1d<double, 6 > >& rVariable,
			array_1d<double, 6 > & rValue);

		/**
		* sets the value of a specified variable
		* @param rVariable the variable to be returned
		* @param rValue new value of the specified variable
		* @param rCurrentProcessInfo the process info
		*/
		virtual void SetValue(
			const Variable<double>& rVariable,
			const double& rValue,
			const ProcessInfo& rCurrentProcessInfo);

		/**
		* sets the value of a specified variable
		* @param rVariable the variable to be returned
		* @param rValue new value of the specified variable
		* @param rCurrentProcessInfo the process info
		*/
		virtual void SetValue(
			const Variable<Vector >& rVariable,
			const Vector& rValue, 
			const ProcessInfo& rCurrentProcessInfo);

		/**
		* sets the value of a specified variable
		* @param rVariable the variable to be returned
		* @param rValue new value of the specified variable
		* @param rCurrentProcessInfo the process info
		*/
		virtual void SetValue(
			const Variable<Matrix >& rVariable,
			const Matrix& rValue, 
			const ProcessInfo& rCurrentProcessInfo);

		/**
		* sets the value of a specified variable
		* @param rVariable the variable to be returned
		* @param rValue new value of the specified variable
		* @param rCurrentProcessInfo the process info
		*/
		virtual void SetValue(
			const Variable<array_1d<double, 3 > >& rVariable,
			const array_1d<double, 3 > & rValue,
			const ProcessInfo& rCurrentProcessInfo);

		/**
		* sets the value of a specified variable
		* @param rVariable the variable to be returned
		* @param rValue new value of the specified variable
		* @param rCurrentProcessInfo the process info
		*/
		virtual void SetValue(
			const Variable<array_1d<double, 6 > >& rVariable,
			const array_1d<double, 6 > & rValue,
			const ProcessInfo& rCurrentProcessInfo);

		/**
		* Is called to check whether the provided material parameters in the Properties
		* match the requirements of current constitutive model.
		* @param rMaterialProperties the current Properties to be validated against.
		* @return true, if parameters are correct; false, if parameters are insufficient / faulty
		* NOTE: this has to be implemented by each constitutive model. Returns false in base class since
		* no valid implementation is contained here.
		*/
		virtual bool ValidateInput(const Properties& rMaterialProperties);

		/**
		* returns the expected strain measure of this constitutive law (by default linear strains)
		* @return the expected strain measure
		*/
		virtual StrainMeasure GetStrainMeasure();

		/**
		* returns the stress measure of this constitutive law (by default 1st Piola-Kirchhoff stress in voigt notation)
		* @return the expected stress measure
		*/
		virtual StressMeasure GetStressMeasure();

		/**
		* returns whether this constitutive model is formulated in incremental strains/stresses
		* NOTE: by default, all constitutive models should be formulated in total strains
		* @return true, if formulated in incremental strains/stresses, false otherwise
		*/
		virtual bool IsIncremental();

		/**
		* This is to be called at the very beginning of the calculation
		* (e.g. from InitializeElement) in order to initialize all relevant
		* attributes of the constitutive law
		* @param rMaterialProperties the Properties instance of the current element
		* @param rElementGeometry the geometry of the current element
		* @param rShapeFunctionsValues the shape functions values in the current integration point
		*/
		virtual void InitializeMaterial(
			const Properties& rMaterialProperties,
			const GeometryType& rElementGeometry,
			const Vector& rShapeFunctionsValues);

		/**
		* to be called at the beginning of each solution step
		* (e.g. from Element::InitializeSolutionStep)
		* @param rMaterialProperties the Properties instance of the current element
		* @param rElementGeometry the geometry of the current element
		* @param rShapeFunctionsValues the shape functions values in the current integration point
		* @param the current ProcessInfo instance
		*/
		virtual void InitializeSolutionStep(
			const Properties& rMaterialProperties,
			const GeometryType& rElementGeometry,
			const Vector& rShapeFunctionsValues,
			const ProcessInfo& rCurrentProcessInfo);

		/**
		* to be called at the end of each solution step
		* (e.g. from Element::FinalizeSolutionStep)
		* @param rMaterialProperties the Properties instance of the current element
		* @param rElementGeometry the geometry of the current element
		* @param rShapeFunctionsValues the shape functions values in the current integration point
		* @param the current ProcessInfo instance
		*/
		virtual void FinalizeSolutionStep(
			const Properties& rMaterialProperties,
			const GeometryType& rElementGeometry,
			const Vector& rShapeFunctionsValues,
			const ProcessInfo& rCurrentProcessInfo);

		/**
		* to be called at the beginning of each step iteration
		* (e.g. from Element::InitializeNonLinearIteration)
		* @param rMaterialProperties the Properties instance of the current element
		* @param rElementGeometry the geometry of the current element
		* @param rShapeFunctionsValues the shape functions values in the current integration point
		* @param the current ProcessInfo instance
		*/
		virtual void InitializeNonLinearIteration(
			const Properties& rMaterialProperties,
			const GeometryType& rElementGeometry,
			const Vector& rShapeFunctionsValues,
			const ProcessInfo& rCurrentProcessInfo);

		/**
		* to be called at the end of each step iteration
		* (e.g. from Element::FinalizeNonLinearIteration)
		* @param rMaterialProperties the Properties instance of the current element
		* @param rElementGeometry the geometry of the current element
		* @param rShapeFunctionsValues the shape functions values in the current integration point
		* @param the current ProcessInfo instance
		*/
		virtual void FinalizeNonLinearIteration(
			const Properties& rMaterialProperties,
			const GeometryType& rElementGeometry,
			const Vector& rShapeFunctionsValues,
			const ProcessInfo& rCurrentProcessInfo);

		/**
		* Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
		* @see Parameters
		*/
		virtual void CalculateMaterialResponsePK1 (Parameters& rValues);

		/**
		* Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
		* @see Parameters
		*/
		virtual void CalculateMaterialResponsePK2 (Parameters& rValues);

		/**
		* Computes the material response in terms of Kirchhoff stresses and constitutive tensor
		* @see Parameters
		*/
		virtual void CalculateMaterialResponseKirchhoff (Parameters& rValues);

		/**
		* Computes the material response in terms of Cauchy stresses and constitutive tensor
		* @see Parameters
		*/
		virtual void CalculateMaterialResponseCauchy (Parameters& rValues);

		/**
		* Updates the material response in terms of 1st Piola-Kirchhoff stresses
		* @see Parameters
		*/
		virtual void FinalizeMaterialResponsePK1 (Parameters& rValues);

		/**
		* Updates the material response in terms of 2nd Piola-Kirchhoff stresses
		* @see Parameters
		*/
		virtual void FinalizeMaterialResponsePK2 (Parameters& rValues);

		/**
		* Updates the material response in terms of Kirchhoff stresses
		* @see Parameters
		*/
		virtual void FinalizeMaterialResponseKirchhoff (Parameters& rValues);

		/**
		* Updates the material response in terms of Cauchy stresses
		* @see Parameters
		*/
		virtual void FinalizeMaterialResponseCauchy (Parameters& rValues);

		/**
		* This can be used in order to reset all internal variables of the
		* constitutive law (e.g. if a model should be reset to its reference state)
		* @param rMaterialProperties the Properties instance of the current element
		* @param rElementGeometry the geometry of the current element
		* @param rShapeFunctionsValues the shape functions values in the current integration point
		* @param the current ProcessInfo instance
		*/
		virtual void ResetMaterial(
			const Properties& rMaterialProperties,
			const GeometryType& rElementGeometry,
			const Vector& rShapeFunctionsValues);

		/**
		* This function is designed to be called once to check compatibility with element
		* @param rFeatures
		*/
		virtual void GetLawFeatures(Features& rFeatures);

		/**
		* This function is designed to be called once to perform all the checks needed
		* on the input provided. Checks can be "expensive" as the function is designed
		* to catch user's errors.
		* @param rMaterialProperties
		* @param rElementGeometry
		* @param rCurrentProcessInfo
		* @return
		*/
		virtual int Check(
			const Properties& rMaterialProperties,
			const GeometryType& rElementGeometry,
			const ProcessInfo& rCurrentProcessInfo);

		// OLD METHOD
		virtual void CalculateMaterialResponse(const Vector& StrainVector,
			const Matrix& DeformationGradient,
			Vector& StressVector,
			Matrix& AlgorithmicTangent,
			const ProcessInfo& rCurrentProcessInfo,
			const Properties& rMaterialProperties,
			const GeometryType& rElementGeometry,
			const Vector& rShapeFunctionsValues,
			bool CalculateStresses = true,
			int CalculateTangent = true,
			bool SaveInternalVariables = true);

	public:

		struct CalculationData
		{
			// elasticity
			double E;
			double nu;
			Matrix C0;
			// tension
			double ft;
			double Gt;
			double m2; // tensile surface s1
			// compression
			double fc0;
			double fcp;
			double fcr;
			double ep;
			double c1;
			double c2;
			double c3;
			double Gc;
			double bm;
			double m1; // compressive reduction
			// effective stress data
			Vector S;
			Vector Si;
			Vector ST;
			Vector SC;
			Matrix PT;
			Matrix PC;
			// misc
			double eta;
			double lch;
			double dTime;
			double rate_coeff_1;
			double rate_coeff_2;
			// different tensile models
			int tensile_damage_model;
			// generalized rankine model
			double grank_a;
			double grank_b;
			double grank_c;
		};

	protected:

		///@name Protected static Member Variables
		///@{
		///@}

		///@name Protected member Variables
		///@{

		bool   m_initialized;
		double m_rt;
		double m_rt_converged;
		double m_rc;
		double m_rc_converged;
		double m_damage_t;
		double m_damage_c;
		double m_lch;
		double m_lch_multiplier;
		Vector m_initial_strain;

#ifdef DAM_TC_2D_IMPLEX
		// testing: for IMPLEX
		double m_rt_converged_old;
		double m_rc_converged_old;
		Vector m_strain;
		Vector m_strain_converged;
		double m_dTime_n;
		double m_dTime_n_converged;
		double m_rt_impl_temp;
		double m_rc_impl_temp;
#endif // DAM_TC_2D_IMPLEX

		double m_suggested_time_step;
		double m_error_code;

		double m_localized;
		double m_localized_converged;
		double m_localization_angle;

#ifdef DAM_TC_2D_INCREMENTAL_REGULARIZATION_V2
		double m_E;
		bool   m_has_changed_reg;
		double m_change_reg_t_x;
		double m_change_reg_c_x;
#endif // DAM_TC_2D_INCREMENTAL_REGULARIZATION_V2


		///@}

		///@name Protected Operators
		///@{

		void InitializeCalculationData(const Properties& props, 
									   const GeometryType& geom, 
									   const Vector& N,
									   const ProcessInfo& pinfo,
									   CalculationData& data);

		void CalculateElasticityMatrix(CalculationData& data);

		void TensionCompressionSplit(CalculationData& data);

		void TrialEquivalentStressTension(CalculationData& data, double& rt_trial);

		void TrialEquivalentStressCompression(CalculationData& data, double& rc_trial);

		void CalculateDamageTension(CalculationData& data, double rt, double& dt);

		void CalculateDamageCompression(CalculationData& data, double rc, double& dc);

		void CalculateMaterialResponseInternal(const Vector& strain_vector, Vector& stress_vector, CalculationData& data);

		///@}

		///@name Protected Operations
		///@{
		///@}


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

		///@name Private  Access
		///@{
		///@}

		///@name Serialization
		///@{

		friend class Serializer;


		virtual void save(Serializer& rSerializer) const
		{
			KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw );
		}

		virtual void load(Serializer& rSerializer)
		{
			KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw );
		}

		///@}

	}; /* Class DamageTCPlaneStress2DLaw */

} /* namespace Kratos.*/
#endif /* KRATOS_DAMAGE_TC_PLANE_STRESS_2D_LAW  defined */
