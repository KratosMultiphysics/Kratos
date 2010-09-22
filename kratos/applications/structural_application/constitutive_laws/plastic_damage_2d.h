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

#if !defined(KRATOS_PLASTIC_DAMAGE_2D_H_INCLUDED)
#define  KRATOS_PLASTIC_DAMAGE_2D_H_INCLUDED

// System includes 
#include <algorithm>


// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "structural_application.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "fluency_criteria/fluency_criteria.h"
#include "soft_hard_behavior/softening_hardening_criteria.h"
#include "includes/constitutive_law.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "custom_utilities/tensor_utils.h"
#include "includes/process_info.h"
#include "includes/properties.h"


namespace Kratos
{
	/**
	 * Defines a linear elastic isotropic constitutive law in 3D space.
	 * This material law is defined by the parameters E (Young's modulus) 
	 * and NU (Poisson ratio)
	 * As there are no further parameters the functionality is limited 
	 * to linear elasticity.
	 */
	class PlasticDamage2D : public ConstitutiveLaw
	{
		public:

		      
		        ///@name Type Definitions
		        typedef ConstitutiveLaw BaseType;
			/**
			 * Counted pointer of PlasticDamage3D
			 */
			typedef boost::shared_ptr<PlasticDamage2D> Pointer;

			typedef FluencyCriteria::Pointer FluencyCriteriaPointer;  

                        typedef SofteningHardeningCriteria::Pointer SofteningHardeningCriteriaPointer;   
                        
                        typedef Properties::Pointer PropertiesPointer;

                        enum Cases{Tracc_Pura, Comp_Pura, Mixto, None}; 




			
			/**
			 * Life Cycle 
			 */
			/**
			 * Default constructor.
			 */
			PlasticDamage2D();
			
                        virtual boost::shared_ptr<ConstitutiveLaw> Clone() const
			{
				boost::shared_ptr<ConstitutiveLaw> 
                                p_clone(new PlasticDamage2D(
                                mpFluencyCriteria->Clone(),
                                mpFluencyCriteria_Traction->Clone(), 
                                mpSofteningBehavior, 
                                mpProperties));
				
                            return p_clone;
			}
						

                        PlasticDamage2D(
                        FluencyCriteriaPointer FluencyCriteriaCompresion,
                        FluencyCriteriaPointer FluencyCriteriaTraction, 
                        SofteningHardeningCriteriaPointer SofteningBehavior, 
                        PropertiesPointer Property);
                      

                       	/**
			 * Destructor.
			 */
			virtual ~PlasticDamage2D();
			
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

		  void Calculate( const Variable<double>& rVariable, 
		  double& Output, 
		  const ProcessInfo& rCurrentProcessInfo);

                  void UpdateMaterial( const Vector& StrainVector,
                                     const Properties& props,
                                     const GeometryType& geom,
                                     const Vector& ShapeFunctionsValues,
                                     const ProcessInfo& CurrentProcessInfo);
                
                void CalculateMaterialResponse( const Vector& StrainVector,
                                            const Matrix& DeformationGradient,
                                            Vector& StressVector,
                                            Matrix& AlgorithmicTangent,
                                            const ProcessInfo& CurrentProcessInfo,
                                            const Properties& props, 
                                            const GeometryType& geom,
                                            const Vector& ShapeFunctionsValues,
                                            bool CalculateStresses = true,
                                            int CalculateTangent = true,
                                            bool SaveInternalVariables = true
                                          );


 
    private:

    protected:
 
    Cases mCase; 
   
 
    bool   mComputeTangentMatrix;
    double mNU;
    double mE;
    double mDE;
    double mlength;
    double mlocal_fail_factor;
    
    double mplastic_damage;
    double mcurrent_plastic_damage;
    double mcohesion;
    double mcurrent_cohesion;
    double mfriction_angle;
    double mmaxfriction_angle;
    double mcurrent_friction_angle;
    double mdilatancy_angle;  
    double mcurrent_dilatancy_angle;               
    double mmaxdilatancy_angle;
    double mdisipation;  
    double mcurrent_disipation; 
    double mefective_plastic_strain;
    double mcurrent_efective_plastic_strain;
   

    array_1d<double,3> mFt;  
    array_1d<double,3> mcurrent_Ft;  
    array_1d<double,4> mplastic_strain;
    array_1d<double,4> mcurrent_plastic_strain;
    Vector             mPrincipalStress;
    array_1d < array_1d < double,3 > ,3> mEigenVectors;
    //array_1d<unsigned int, 3 > mOrder; 

  

    vector<array_1d<double,4> > mDerivate_Potencial;
    vector<array_1d<double,4> > mDerivate_Fluency;
    std::vector<unsigned int>   mactive_surface;
    std::vector<unsigned int>   mactive_surface_T;
         

    FluencyCriteriaPointer mpFluencyCriteria;
    FluencyCriteriaPointer mpFluencyCriteria_Traction;

    SofteningHardeningCriteriaPointer mpSofteningBehavior;
    PropertiesPointer mpProperties;


   void CalculateElasticMatrix(boost::numeric::ublas::bounded_matrix<double,4,4>& C);
   void CalculateElasticStress(array_1d<double,4>& Strain, array_1d<double,4>& Stress);
   void Compute_Principal_Stress(const Vector& StressVector);
   void IdentifyMaximunAndMinumumPrincipalStres_CalculateOrder(array_1d<double, 3 >& PrincipalStress , array_1d<unsigned int,3>& order);
   void AssembleUpdateStressAndStrainTensor(
        const array_1d<double,3>& Sigma,  
        const array_1d<double,4>& StrainVector_aux,   
        const array_1d<unsigned int,3>& order,
	array_1d<double,4>& ElasticStrain,                                                
	array_1d<double,4>& ElasticStress);     



	}; // Class PlasticDamage2D 
}  // namespace Kratos.
#endif // KRATOS_ISOTROPIC_2D_H_INCLUDED  defined 
