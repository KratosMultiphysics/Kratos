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
*   Last Modified by:    $Author: nelson $
*   Date:                $Date: 2008-09-03 
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/

#if !defined(KRATOS_ISOTROPIC_DAMAGE_H_INCLUDED )
#define  KRATOS_ISOTROPIC_DAMAGE_H_INCLUDED

// System includes 

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "fluency_criteria/fluency_criteria.h"
#include "soft_hard_behavior/softening_hardening_criteria.h"


namespace Kratos
{
	/**
	 * Defines a linear elastic isotropic constitutive law in 3D space.
	 * This material law is defined by the parameters E (Young's modulus) 
	 * and NU (Poisson ratio)
	 * As there are no further parameters the functionality is limited 
	 * to linear elasticity.
	 */
	class Isotropic_Damage : public ConstitutiveLaw<Node<3> >
	{
		public:

		      
		        ///@name Type Definitions
		        typedef ConstitutiveLaw<Node<3> > BaseType;
			/**
			 * Counted pointer of Isotropic_Damage
			 */
			typedef boost::shared_ptr<Isotropic_Damage> Pointer;

			typedef typename FluencyCriteria::Pointer FluencyCriteriaPointer;  

                        typedef  typename SofteningHardeningCriteria::Pointer SofteningHardeningCriteriaPointer;   
                        
                        typedef typename Properties::Pointer PropertiesPointer;

			
			/**
			 * Life Cycle 
			 */
			/**
			 * Default constructor.
			 */
			Isotropic_Damage();
                        Isotropic_Damage(FluencyCriteriaPointer FluencyCriteria, SofteningHardeningCriteriaPointer SofteningBehavior, PropertiesPointer Property);
			
			virtual boost::shared_ptr<ConstitutiveLaw<Node<3> > > Clone() const
			{
				boost::shared_ptr<ConstitutiveLaw<Node<3> > > p_clone(new Isotropic_Damage(mFluencyCriteria,mSofteningBehavior, mProperties));
				return p_clone;
			}

			/**
			 * Destructor.
			 */
			virtual ~Isotropic_Damage();
			
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
			

			  void InitializeMaterial( const Properties& props,
					const GeometryType& geom,
					const Vector& ShapeFunctionsValues );

			  void CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult);

			  void CalculateStress( const Vector& StrainVector, 
                                          Vector& StressVector);
          	
			  void CalculateCauchyStresses( Vector& Cauchy_StressVector,
					const Matrix& F,
					const Vector& PK2_StressVector,
					const Vector& GreenLagrangeStrainVector);
			

    void Calculate(const Variable<Matrix >& rVariable, Matrix& rResult, 
                           const ProcessInfo& rCurrentProcessInfo);

			 void InitializeSolutionStep( const Properties& props,
                    const GeometryType& geom,
                    const Vector& ShapeFunctionsValues ,
                    const ProcessInfo& CurrentProcessInfo);

    void FinalizeSolutionStep( const Properties& props,
                                               const GeometryType& geom, 
                                               const Vector& ShapeFunctionsValues ,
                                               const ProcessInfo& CurrentProcessInfo);

	  void CalculateStressAndTangentMatrix( Vector& StressVector,
                    const Vector& StrainVector,
                    Matrix& algorithmicTangent);


			// Metodo para calcular los factores de dano
			

  // d = variable de dano a calcular. 
 
		
		private:

			/**
			 * there are no protected class members
			 */
		protected:
			
				// Atributos  Privados
    double mEc;
    double mEt; 
    double mFc;
    double mFt;
    double mGE;
    double mNU;
    double ml;
    double md;
    double mr_old;
    double mr_new; 
    double mr_o;
    double mArea;


/*    enum myState
         {PlaneStress=1, PlaneStrain, Tri_D }; 
     
    myState State; */ 
    FluencyCriteriaPointer mFluencyCriteria;
    SofteningHardeningCriteriaPointer mSofteningBehavior;
    PropertiesPointer mProperties;


			// Miembros Privados
   void CalculateNoDamageElasticMatrix(Matrix& C, const double E, const double NU);
   void CalculateNoDamageStress(const Vector& StrainVector, Vector& rResult);
   void CalculateDamage(const Matrix& ConstitutiveMatrix, const Vector& Stress,const Vector& StressVector);


	}; // Class Isotropic_Damage 
}  // namespace Kratos.
#endif // KRATOS_ISOTROPIC_2D_H_INCLUDED  defined 
