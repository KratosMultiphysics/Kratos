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
#include "structural_application.h"
#include "includes/define.h"
#include "includes/serializer.h"
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
	class Isotropic_Damage_2D: public ConstitutiveLaw
	{
		public:

		      
		        ///@name Type Definitions
		        typedef ConstitutiveLaw BaseType;
			/**
			 * Counted pointer of Isotropic_Damage_2D
			 */
			typedef boost::shared_ptr<Isotropic_Damage_2D> Pointer;

			typedef FluencyCriteria::Pointer FluencyCriteriaPointer;  

                        typedef SofteningHardeningCriteria::Pointer SofteningHardeningCriteriaPointer;   
                        
                        typedef Properties::Pointer PropertiesPointer;

			
			/**
			 * Life Cycle 
			 */
			/**
			 * Default constructor.
			 */
			Isotropic_Damage_2D();
                        Isotropic_Damage_2D(FluencyCriteriaPointer FluencyCriteria, SofteningHardeningCriteriaPointer SofteningBehavior, PropertiesPointer Property);
			
			virtual boost::shared_ptr<ConstitutiveLaw> Clone() const
			{
				boost::shared_ptr<ConstitutiveLaw> p_clone(new Isotropic_Damage_2D(
				mpFluencyCriteria->Clone(),
				mpSofteningBehavior->Clone(), 
				mpProperties));
				return p_clone;
			}

			/**
			 * Destructor.
			 */
			virtual ~Isotropic_Damage_2D();
			
			/**
			 * Operators 
			 */
			/**
			 * Operations
			 */
			bool Has( const Variable<double>& rThisVariable );
			bool Has( const Variable<Vector>& rThisVariable );
			bool Has( const Variable<Matrix>& rThisVariable );

			double& GetValue( const Variable<double>& rThisVariable, double& rValue );
			Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue );
			Matrix& GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue );

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
			
			void CalculateMaterialResponse (const Vector& StrainVector,
			const Matrix& DeformationGradient,
			Vector& StressVector,
			Matrix& AlgorithmicTangent,
			const ProcessInfo& CurrentProcessInfo,
			const Properties& props,
			const GeometryType& geom,
			const Vector& ShapeFunctionsValues,
			bool CalculateStresses,
			int CalculateTangent,
			bool SaveInternalVariables
                                                   );
						   
						   
	                 std::string Info() const;
			 int Check(const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo);

 
		
		private:

		  

		///@}
		///@name Serialization
		///@{	
		friend class Serializer;

		virtual void save(Serializer& rSerializer) const
		{
		rSerializer.save("Name","Isotropic_Damage_2D");
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
		}

		virtual void load(Serializer& rSerializer)
		{
		KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
		}
		  
		/// Internal Variables  
		double ml;
		double md_new;
		double md_old;
		double mr_old;
		double mr_new;
		bool mdamage_calculated;
		
		protected:
		
		
			



/*    enum myState
         {PlaneStress=1, PlaneStrain, Tri_D }; 
     
    myState State; */ 
    FluencyCriteriaPointer mpFluencyCriteria;
    SofteningHardeningCriteriaPointer mpSofteningBehavior;
    PropertiesPointer mpProperties;


   // Miembros Privados
   void CalculateNoDamageElasticMatrix(Matrix& C);
   void CalculateNoDamageStress(const Vector& StrainVector, Vector& rResult);
   void CalculateStressAndDamage(Vector& StressVector, const Vector& StrainVector);
   void CalculateStressAndTangentMatrix(const Vector& StrainVector, Matrix& algorithmicTangent);
   void CalculateStress( const Vector& StrainVector, Vector& StressVector);
   void CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult);
   void CalculateCauchyStresses( Vector& Cauchy_StressVector,const Matrix& F, const Vector& PK2_StressVector,const Vector& GreenLagrangeStrainVector);
   void Calculate( const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo);  
   void UpdateMaterial();
   //void CalculateDamage(const Matrix& ConstitutiveMatrix, const Vector& StressVector, double& d);
   
	}; // Class Isotropic_Damage_2D_2D
}  // namespace Kratos.
#endif // KRATOS_ISOTROPIC_2D_H_INCLUDED  defined 
