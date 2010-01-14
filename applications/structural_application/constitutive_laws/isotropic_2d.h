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
*   Last Modified by:    $Author: kazem $
*   Date:                $Date: 2009-01-16 10:50:17 $
*   Revision:            $Revision: 1.11 $
*
* ***********************************************************/

#if !defined(KRATOS_ISOTROPIC_2D_H_INCLUDED )
#define  KRATOS_ISOTROPIC_2D_H_INCLUDED

// System includes 

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
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
	class Isotropic2D : public ConstitutiveLaw<Node<3> >
	{
		public:
			/**
			 * Type Definitions
			 */
			typedef ConstitutiveLaw<Node<3> > BaseType;
			/**
			 * Counted pointer of Isotropic2D
			 */
			typedef boost::shared_ptr<Isotropic2D> Pointer;
			
			/**
			 * Life Cycle 
			 */
			/**
			 * Default constructor.
			 */
			Isotropic2D();
			
			virtual boost::shared_ptr<ConstitutiveLaw<Node<3> > > Clone() const
			{
				boost::shared_ptr<ConstitutiveLaw<Node<3> > > p_clone(new Isotropic2D());
				return p_clone;
			}

			/**
			 * Destructor.
			 */
			virtual ~Isotropic2D();
			
			/**
			 * Operators 
			 */
			/**
			 * Operations
			 */
			bool Has( const Variable<double>& rThisVariable );
			bool Has( const Variable<Vector>& rThisVariable );
			bool Has( const Variable<Matrix>& rThisVariable );
			
			double GetValue( const Variable<double>& rThisVariable );
			Vector GetValue( const Variable<Vector>& rThisVariable );
			Matrix GetValue( const Variable<Matrix>& rThisVariable );
			
            void SetValue( const Variable<double>& rVariable, 
                           const double& Value, 
                           const ProcessInfo& rCurrentProcessInfo );
            void SetValue( const Variable<Vector>& rThisVariable, 
                           const Vector& rValue, 
                           const ProcessInfo& rCurrentProcessInfo );
            void SetValue( const Variable<Matrix>& rThisVariable, 
                           const Matrix& rValue, 
                           const ProcessInfo& rCurrentProcessInfo );
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
 	//		void PlaneStrainConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult);
			
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
/*			void InitializeSolutionStep( const Properties& props,
					const GeometryType& geom, //this is just to give the array of nodes
					const Vector& ShapeFunctionsValues ,
					const ProcessInfo& CurrentProcessInfo);
			
			void FinalizeSolutionStep( const Properties& props,
					const GeometryType& geom, //this is just to give the array of nodes
					const Vector& ShapeFunctionsValues ,
					const ProcessInfo& CurrentProcessInfo);
*/			
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
            void Calculate(const Variable<Matrix >& rVariable, Matrix& rResult, 
                           const ProcessInfo& rCurrentProcessInfo);
            
            void Calculate( const Variable<double>& rVariable, 
                                    double& Output, 
                                    const ProcessInfo& rCurrentProcessInfo);

			void CalculateStressAndTangentMatrix(Vector& StressVector,
			      const Vector& StrainVector,
			      Matrix& algorithmicTangent);

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
			
			/**
			 * Static Member Variables 
			 */

			 /**
			  * calculates the linear elastic constitutive matrix in terms of Young's modulus and
			  * Poisson ratio
			  * @param E the Young's modulus
			  * @param NU the Poisson ratio
			  * @return the linear elastic constitutive matrix
			  */
			 void CalculateElasticMatrix(Matrix& C, const double E, const double NU);

		         double mE,mNU,mDE;

 	//		 void CalculatePlaneElasticMatrix(Matrix& C, const double E, const double NU);


			 /**
			  * Un accessible methods 
			  */
			 /**
			  * Assignment operator.
			  */
			 //Isotropic2D& operator=(const IsotropicPlaneStressWrinklingNew& rOther);
			 /**
			  * Copy constructor.
			  */
			 //Isotropic2D(const IsotropicPlaneStressWrinklingNew& rOther);
	}; // Class Isotropic2D 
}  // namespace Kratos.
#endif // KRATOS_ISOTROPIC_2D_H_INCLUDED  defined 
