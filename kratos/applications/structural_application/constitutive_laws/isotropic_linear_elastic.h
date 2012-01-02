/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

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
*   Last Modified by:    $Author: janosch $
*   Date:                $Date: 2007-09-20 10:52:55 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/

#if !defined(KRATOS_ISOTROPIC_LINEAR_ELASTIC_H_INCLUDED )
#define  KRATOS_ISOTROPIC_LINEAR_ELASTIC_H_INCLUDED

// System includes 

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"


namespace Kratos
{
	/**
	 * Defines a linear elastic isotropic constitutive law in 3D space.
	 * This material law is defined by the parameters E (Young's modulus) 
	 * and NU (Poisson ratio)
	 * As there are no further parameters the functionality is limited 
	 * to linear elasticity.
	 */
	class IsotropicLinearElastic : public ConstitutiveLaw
	{
		public:
			/**
			 * Type Definitions
			 */
			typedef ConstitutiveLaw BaseType;
			/**
			 * Counted pointer of IsotropicLinearElastic
			 */
			typedef boost::shared_ptr<IsotropicLinearElastic> Pointer;
			
			/**
			 * Life Cycle 
			 */
			/**
			 * Default constructor.
			 */
			IsotropicLinearElastic();

			/**
			 * Destructor.
			 */
			virtual ~IsotropicLinearElastic();
			
			/**
			 * Clone function
			 */
			virtual boost::shared_ptr<ConstitutiveLaw> Clone() const
			{
				boost::shared_ptr<ConstitutiveLaw> p_clone(new IsotropicLinearElastic());
				return p_clone;
			}
			
			/**
			 * Operators 
			 */
			bool Has( const Variable<double>& rThisVariable );
			bool Has( const Variable<Vector>& rThisVariable );
			bool Has( const Variable<Matrix>& rThisVariable );
			double GetValue( const Variable<double>& rThisVariable );
			Vector GetValue( const Variable<Vector>& rThisVariable );
			Matrix GetValue( const Variable<Matrix>& rThisVariable );
			void SetValue( const Variable<double>& rThisVariable, double rValue );
			void SetValue( const Variable<Vector>& rThisVariable, Vector rValue );
			void SetValue( const Variable<Matrix>& rThisVariable, Matrix rValue );
			/**
			 * Operations
			 */
			/**
			 * Material parameters are inizialized
			 */ 
			void InitializeMaterial( const Properties& props,
					const GeometryType& geom,
					const Vector& ShapeFunctionsValues );
			
			/**
			 * Calculates the constitutive matrix for a given strain vector
			 * @param StrainVector the current vector of strains the constitutive 
			 * matrix is to be generated for
			 * @param rResult Matrix the result will be stored in
			 */
			void CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult);
			
			/**
			 * Calculates the stresses for given strain state
			 * @param StrainVector the current vector of strains
			 * @param rResult the stress vector corresponding to the given strains
			 */
			void CalculateStress(const Vector& StrainVector, Vector& rResult);
			
			/**
			 * As this constitutive law describes only linear elastic material properties
			 * this function is rather useless and in fact does nothing
			 */ 
			void UpdateMaterial( const Vector& StrainVector,
					const Properties& props,
					const GeometryType& geom,
					const Vector& ShapeFunctionsValues ,
					const ProcessInfo& CurrentProcessInfo);
			
			void InitializeSolutionStep( const Properties& props,
					const GeometryType& geom,
					const Vector& ShapeFunctionsValues ,
					const ProcessInfo& CurrentProcessInfo);
			
			void FinalizeSolutionStep( const Properties& props,
					const GeometryType& geom,
					const Vector& ShapeFunctionsValues ,
					const ProcessInfo& CurrentProcessInfo);
			
			/**
			 * Calculates the cauchy stresses. For a given deformation and stress state
			 * the cauchy stress vector is calculated
			 * @param Cauchy_StressVector the result vector of cauchy stresses (will be overwritten)
			 * @param F the current deformation gradient
			 * @param PK2_StressVector the current second Piola-Kirchhoff-Stress vector
			 * @param GreenLagrangeStrainVector the current Green-Lagrangian strains
			 */
			void CalculateCauchyStresses( Vector& Cauchy_StressVector,
					const Matrix& F,
					const Vector& PK2_StressVector,
					const Vector& GreenLagrangeStrainVector);
			
			
			/**
			 * converts a strain vector styled variable into its form, which the
			 * deviatoric parts are no longer multiplied by 2
			 */
			void Calculate(const Variable<Matrix >& rVariable, Matrix& rResult);
			
			void Calculate(const Variable<Vector >& rVariable, Vector& rResult);

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
			/**
			 * there are no protected class members
			 */
		private:

		///@}
		///@name Serialization
		///@{	
		friend class Serializer;

		virtual void save(Serializer& rSerializer) const
		{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
		}

		virtual void load(Serializer& rSerializer)
		{
		KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
		}
			
			/**
			 * Static Member Variables 
			 */
			double mE;
			double mNU;
			Matrix mCtangent;
			Vector mEw;
			Vector mCurrentStress;
			Vector mInSituStress;
			/**
			 * Member Variables 
			 */
			/**
			 * calculates the eigenvectors and eigenvalues of the current stress state.
			 * Eigenvalues represent the principal stresses while the eigenvectors
			 * represent the directions of principal stresses
			 * @param StressVector the current stress state vector
			 * @param smin the minimal principal stress (will be overwritten by 
			 * the solution)
			 * @param mineigenvect the eigenvector corresponding to smin (will be 
			 * overwritten by the solution)
			 * @param smid the intermediate principal stress (will be overwritten 
			 * by the solution)
			 * @param mideigenvect the eigenvector corresponding to smid (will be 
			 * overwritten by the solution)
			 * @param smax the maximum principal stress (will be overwritten by 
			 * the solution)
			 * @param maxeigenvect the eigenvector corresponding to smax (will be 
			 * overwritten by the solution)
			 */
			 void CalculateStressEigNonNormalized( const Vector& StressVector,
					 double& smin, Vector& mineigenvect,
					 double& smid, Vector& mideigenvect,
					 double& smax, Vector& maxeigenvect);
			 
			 /**
			  * calculates the normalized eigenvectors and eigenvalues of the current 
			  * stress state.
			  * Eigenvalues represent the principal stresses while the eigenvectors
			  * represent the directions of principal stresses. All eigenvectors 
			  * are normalized to length '1'. 
			  * @param StressVector the current stress state vector
			  * @param smin the minimal principal stress (will be overwritten by 
			  * the solution)
			  * @param mineigenvect the eigenvector corresponding to smin (will  
			  * be overwritten by the solution)
			  * @param smid the intermediate principal stress (will be overwritten 
			  * by the solution)
			  * @param mideigenvect the eigenvector corresponding to smid (will 
			  * be overwritten by the solution)
			  * @param smax the maximum principal stress (will be overwritten by 
			  * the solution)
			  * @param maxeigenvect the eigenvector corresponding to smax (will 
			  * be overwritten by the solution)
			  * @see CalculateStressEigNonNormalized
			  */    
			 void CalculateStressEig( const Vector& StressVector,
					 double& smin, Vector& mineigenvect,
					 double& smid, Vector& mideigenvect,
					 double& smax, Vector& maxeigenvect);
			  /**
			   * calculates the principal strains for given strain state. the principal strains
			   * are given sorted by size 
			   * @param StrainVector the current strain state
			   * @param eps1 the minimal principal strain (will be overwritten)
			   * @param eps2 the intermediate principal strain (will be overwritten)
			   * @param eps3 the maximum principal strain (will be overwritten)
			   * @see PrincipSTRESS
			   */
			 void PrincipSTRAIN( const Vector& StrainVector, double& eps1, 
					     double& eps2, double& eps3);
			 /**
			  * calculates the principal stresses for given stress state. 
			  * The principal stresses are given sorted by size 
			  * @param StressVector the current stress state
			  * @param str1 the minimal principal stress (will be overwritten)
			  * @param str2 the intermediate principal stress (will be overwritten)
			  * @param str3 the maximum principal stress (will be overwritten) 
			  */
			 void PrincipSTRESS( const Vector& StressVector, double& str1, 
					     double& str2, double& str3);
			 /**
			  * as the tangent matrix equals the elastic constitutive matrix, 
			  * this function does nothing
			  */
			 void CalculateTangentMatrix(const Vector& StrainVector);
			 
			 /**
			  * calculates the linear elastic constitutive matrix in terms of Young's modulus and
			  * Poisson ratio
			  * @param E the Young's modulus
			  * @param NU the Poisson ratio
			  * @return the linear elastic constitutive matrix
			  */
// 			 Matrix CalculateElasticMatrix(const double E, const double NU);
// 			 Matrix CalculateDirectionVariationTerm( const Matrix& Celastic, 
// 					 const Vector& ElasticStrain);
	}; // Class IsotropicLinearElastic 
}  // namespace Kratos.
#endif // KRATOS_ISOTROPIC_LINEAR_ELASTIC_H_INCLUDED  defined 
