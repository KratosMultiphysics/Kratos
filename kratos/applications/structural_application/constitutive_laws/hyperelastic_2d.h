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
*   Last Modified by:    $Author: virginia $
*   Date:                $Date: 2008-12-04 16:14:58 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/

#if !defined(KRATOS_HYPERELASTIC_2D_H_INCLUDED )
#define  KRATOS_HYPERELASTIC_2D_H_INCLUDED

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
	 * Defines a hyperelastic constitutive law in 2D space for Membrane and EBST elements.
	 * This material law, also called Green-elastic material
	 * postulates the existance of a Helmholtz free-energy
	 * function (psi), which is defined per unit reference volume.
	 * For convenience, the strain-energy function vanishes in the reference
	 * configuration: psi(I)=0
	 * "is defined by the parameters E (Young's modulus) 
	 * and NU (Poisson ratio), and also, by the evaluation of the 
	 * Helmholtz free-energy function.
	 * As there are no further parameters the functionality is applied 
	 * to non-linear elasticity.
	 */   

	// First attempt: ISOTROPIC Hyperelastic Material


	class Hyperelastic2D : public ConstitutiveLaw

	{
		public:

				/**
			 * Type Definitions
			 */

			typedef ConstitutiveLaw BaseType;

			/**
			 * Counted pointer of Hyperelastic3D
			 */
			typedef boost::shared_ptr<Hyperelastic2D> Pointer; // 
			
			/**
			 * Life Cycle 
			 */
			/**
			 * Default constructor.
			 */
			Hyperelastic2D();
			


			virtual boost::shared_ptr<ConstitutiveLaw > Clone() const

			{

				boost::shared_ptr<ConstitutiveLaw > p_clone(new Hyperelastic2D());

				return p_clone;
			}

			/**
			 * Destructor.
			 */
			virtual ~Hyperelastic2D();
			


			//*****************
			/**
			 * Operators 
			 */
			/**
			 * Operations
			 */
			bool Has( const Variable<double>& rThisVariable );
			bool Has( const Variable<Vector>& rThisVariable );
			bool Has( const Variable<Matrix>& rThisVariable );

// NEW
			double GetValue( const Variable<double>& rThisVariable );
			Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue );
			Matrix GetValue( const Variable<Matrix>& rThisVariable );
// NEW
			
// 			double GetValue( const Variable<double>& rThisVariable );
// 			Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue );
// 			Matrix GetValue( const Variable<Matrix>& rThisVariable  );
// 			
// 			void SetValue( const Variable<double>& rThisVariable, const double& rValue, 
// 							  const ProcessInfo& rCurrentProcessInfo );
// 			void SetValue( const Variable<array_1d<double, 3> >& rThisVariable, 
// 							  const array_1d<double, 3>& rValue, const ProcessInfo& rCurrentProcessInfo );
// 			void SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
// 							  const ProcessInfo& rCurrentProcessInfo );
// 			void SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
// 							  const ProcessInfo& rCurrentProcessInfo );
	
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
			
			/**
			 * Calculates the stresses for given strain state
			 * @param StrainVector the current vector of strains
			 * @param rResult the stress vector corresponding to the given strains
			 */


			void CalculateStress(const Vector& StrainVector, Vector& rResult);
			

			void CalculateStressVector(const Vector& StrainVector, Vector& rResult);
                        //void CalculateStressVector(const Vector& StrainVector, array_1d<double, 6>& rResult);
 

			void CalculateEnergy(const Vector& StrainVector, Vector& Energy);

			//******************************************

			/**
			 * As this constitutive law describes Hyperelastic material properties
			 * this function is necessary  /// MODIFFY BELOW!!!!!!!!!!!!!!
			 */ 		
			void InitializeSolutionStep( const Properties& props,
					const GeometryType& geom, //this is just to give the array of nodes
					const Vector& ShapeFunctionsValues ,
					const ProcessInfo& CurrentProcessInfo);
			
            void ResetMaterial(  const Properties& props,
                                 const GeometryType& geom, //this is just to give the array of nodes
                                 const Vector& ShapeFunctionsValues);
            
            void FinalizeSolutionStep( const Properties& props,
					const GeometryType& geom, //this is just to give the array of nodes
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
			

           void CalculateStressAndTangentMatrix( Vector& StressVector,
                                 const Vector& StrainVector,
				 Matrix& algorithmicTangent);

          
            void Calculate(const Variable<Matrix >& rVariable, Matrix& rResult, 
                           const ProcessInfo& rCurrentProcessInfo);
            
            void Calculate( const Variable<double>& rVariable, 
                                    double& Output, 
                                    const ProcessInfo& rCurrentProcessInfo);


//  void Calculate( const Variable<double>& rVariable,
//                                  double& Output,
//                                  const ProcessInfo& rCurrentProcessInfo);

           void CalculateMaterialResponse( const Vector& StrainVector,
                                 const Matrix& DeformationGradient,
				 Vector& StressVector,
				 Matrix& AlgorithmicTangent,
				 const ProcessInfo& CurrentProcessInfo,
				 const Properties& props,
                                 const GeometryType& geom, //this is just to give the array of nodes
                                 const Vector& ShapeFunctionsValues ,
				 bool CalculateStresses = true,
				 int CalculateTangent = true,
				 bool SaveInternalVariables = true);

	
	   int Check(const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo);
			
	                 /**
			 * converts a strain vector styled variable into its form, which the
			 * deviatoric parts are no longer multiplied by 2
			 */
//             void Calculate(const Variable<Matrix >& rVariable, Matrix& rResult, const ProcessInfo& rCurrentProcessInfo);

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
		  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw );
		}

		virtual void load(Serializer& rSerializer)
		{
	 	   KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw );
		}
		  
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
	/*		 void CalculateElasticMatrix(Matrix& C, const double E, const double NU);

		         double mE,mNU;
				 double mMU;
				 Vector mInSituStress;
				 Matrix mCtangent;
				 Vector mCurrentStress;*/

			double mMU;
			double mThickness; //shell thickness
			double mcurrentThickness; // updated shell thickenss
			double mA0; //initial area
			double mA; // current area
			double mAlpha; // angle used for retraction in E0 // RETRACTION
			double mRetractionTime; // time at which retraction starts to act  // RETRACTION
			double mRetraction;  // RETRACTION

			void CalculateElasticMatrix(Matrix& C, const double E, const double NU); //ISOTROPIC
			double mK; //ISOTROPIC
			double mE,mNU; //ISOTROPIC
			Matrix mCtangent;	 //ISOTROPIC
			Vector mMaterialParameters;	 
			Vector mInSituStress;
			Vector auxStressVector;
			//array_1d<double,6> auxStressVector; // member of the class  mauxStressVector; 
			Matrix Ctang;
			Matrix c;
			Matrix CC;
			
				 
			 /**
			  * Un accessible methods 
			  */
			 /**
			  * Assignment operator.
			  */
			 //Isotropic3D& operator=(const IsotropicPlaneStressWrinklingNew& rOther);
			 /**
			  * Copy constructor.
			  */
			 //Isotropic3D(const IsotropicPlaneStressWrinklingNew& rOther);
//******************************************************

	}; // Class Hyperelastic3D 
}  // namespace Kratos.
#endif // KRATOS_HYPERELASTIC_2D_H_INCLUDED  defined 
