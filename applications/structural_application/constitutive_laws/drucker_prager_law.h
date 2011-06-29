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
*   Date:                $Date: 2008-01-25 08:37:21 $
*   Revision:            $Revision: 1.3 $
*
* ***********************************************************/
#if !defined(KRATOS_DRUCKER_PRAGER_LAW_H_INCLUDED )
#define  KRATOS_DRUCKER_PRAGER_LAW_H_INCLUDED

// System includes 

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/constitutive_law.h"
#include "spaces/ublas_space.h"


namespace Kratos
{
    /**
     * Defines a linear elastic isotropic constitutive law in 3D space.
     * This material law is defined by the parameters E (Young's modulus) 
     * and NU (Poisson ratio)
     * As there are no further parameters the functionality is limited 
     * to linear elasticity.
     */
    class DruckerPragerLaw : public ConstitutiveLaw
    {
        public:
            /**
             * Type Definitions
             */
            typedef UblasSpace<double, CompressedMatrix, Vector> SpaceType;
            typedef ConstitutiveLaw BaseType;
            typedef array_1d<double, 81> MaterialTensorType;
            /**
             * Counted pointer of DruckerPragerLaw
             */
            typedef boost::shared_ptr<DruckerPragerLaw> Pointer;
            
            
            /**
             * Life Cycle 
             */
            
            /**
             * Default constructor.
             */
            DruckerPragerLaw();
            
            /**
             * Constructor
             * @param E the young's modulus of the specified material
             * @param NU the poisson ratio of the specified material
             */
            DruckerPragerLaw(const double E, const double NU);
            
            /**
             * Destructor.
             */
            virtual ~DruckerPragerLaw();
            
            boost::shared_ptr<ConstitutiveLaw> Clone() const;
            
            void InitializeMaterial( const Properties& props,
					const GeometryType& geom,
					const Vector& ShapeFunctionsValues);

			void InitializeSolutionStep( const Properties& props,
				const GeometryType& geom, //this is just to give the array of nodes
				const Vector& ShapeFunctionsValues ,
				const ProcessInfo& CurrentProcessInfo);
            
            void FinalizeSolutionStep( const Properties& props,
					const GeometryType& geom, const Vector& ShapeFunctionsValues ,const ProcessInfo& CurrentProcessInfo);
            
 			void SetValue( const Variable<Matrix >& rVariable, 
				const Matrix& Value, const ProcessInfo& rCurrentProcessInfo);

 			void SetValue( const Variable<Vector >& rVariable, 
				const Vector& rValue, const ProcessInfo& rCurrentProcessInfo);

 			Matrix& GetValue(const Variable<Matrix>& rVariable, Matrix& rValue);

    		Vector& GetValue(const Variable<Vector>& rVariable, Vector& rValue);

    		double& GetValue(const Variable<double>& rVariable, double& rValue);
            
//             template<class TVariableType> bool Has( const TVariableType& rThisVariable);
            /**
             * Operators 
             */
            
            
            /**
             * Operations
             */

            
            /**
             * calculates the current stress and the material tangent
             * NOTE: there are two versions of this function: one for a matrix representation
             * of the material tensor and one for a tensorial formulation. Each ConstitutiveLaw
             * HAS TO IMPLEMENT both of them (for convenience, there are conversation functions 
             * available in MathUtils for either of them)
             * Calculates the StressTensor and the algorithmic tangent Matrix for a given 
             * Elastic Left Cauchy Green Tensor in trial state after Simo 
             *[Comp. Meth. in Appl. Mech. and Eng. 99 (1992) 61-112] 
             * @param StressTensor 3times3 Kirchhoff Stress Tensor
             * @param LeftCauchyGreen_Trial elastic Left Cauchy Green Tensor in trial state =
             * delta[f_(n+1)]*b^e_n*delta[f_(n+1)]^T
             * @param algorithmicTangent \frac{\delta \tau_{n+1}}{\delta b_{n+1}}

             */
            virtual void CalculateStressAndTangentMatrix( Matrix& StressTensor,
                    const Matrix& LeftCauchyGreenTensor_trial,
                    MaterialTensorType& algorithmicTangent);
            
            /**
             * calculates the current stress and the material tangent
             * NOTE: there are two versions of this function: one for a matrix representation
             * of the material tensor and one for a tensorial formulation. Each ConstitutiveLaw
             * HAS TO IMPLEMENT both of them (for convenience, there are conversation functions 
             * available in MathUtils for either of them)
             * @param StressVector the calculated stress vector 
             * @param StrainVector the given strain vector
             * @param algorithmicTangent the calculated algorithmic tangent matrix
             */
            virtual void CalculateStressAndTangentMatrix( Vector& StressVector,
                    const Vector& StrainVector,
                    Matrix& algorithmicTangent);
            
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
                                      bool SaveInternalVariables = true );
            
				      
				      
		  private:

		  ///@}
		  ///@name Serialization
		  ///@{	
		  friend class Serializer;

		  virtual void save(Serializer& rSerializer) const
		  {
		     rSerializer.save("Name","DruckerPragerLaw");
		     KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw  );
		  }

		  virtual void load(Serializer& rSerializer)
		  {
		     KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw );
		  }
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
            
			Matrix mOldPlasticStrains;
			Matrix mCurrentPlasticStrains;
			Matrix mInsituStress;
			Matrix mCurrentStress;
			double mE;
			double mNU;
			double mFriction;
			double mCohesion;
			double mTOL1;
			double mTOL2;
			unsigned int mMaxIter;
			MaterialTensorType mI_dev;
			MaterialTensorType mElasticMaterialTensor;
			void SetOldPlasticStrains(Matrix& rValue);
            /**
             * Operations
             */

			Matrix CalculateStressTensor(Matrix& StressTensor,const Matrix& StrainTensor, Matrix& PlasticStrainTensor);

			Matrix CalculateElasticTangent(Matrix& tanC);

			Matrix InverseC(Matrix& InvC);
    }; // Class DruckerPragerLaw 
    
}  // namespace Kratos.

#endif // KRATOS_ISOTROPIC_LINEAR_ELASTIC_H_INCLUDED  defined 



