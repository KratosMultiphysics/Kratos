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
*   Date:                $Date: 2008-10-23 12:22:22 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/
#if !defined(KRATOS_HOOKS_LAW_H_INCLUDED )
#define  KRATOS_HOOKS_LAW_H_INCLUDED

// System includes 

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
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
    class GroutingMortar : public ConstitutiveLaw<Node<3> >
    {
        public:
            /**
             * Type Definitions
             */
            
            /**
             * Counted pointer of GroutingMortar
             */
            typedef boost::shared_ptr<GroutingMortar> Pointer;
            
            
            /**
             * Life Cycle 
             */
            
            /**
             * Default constructor.
             */
            GroutingMortar();
            
            /**
             * Constructor
             * @param E the young's modulus of the specified material
             * @param NU the poisson ratio of the specified material
             */
            GroutingMortar(const double E, const double NU);
            
            /**
             * Destructor.
             */
            virtual ~GroutingMortar();
            
            boost::shared_ptr<ConstitutiveLaw<Node<3> > > Clone() const;
            
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

 	          void SetValue( const Variable<double >& rVariable, 
					const double& rValue, const ProcessInfo& rCurrentProcessInfo);

 			Matrix GetValue(const Variable<Matrix>& rVariable);

    		Vector GetValue(const Variable<Vector>& rVariable);
            
			double GetValue(const Variable<double>& rVariable);
//             template<class TVariableType> bool Has( const TVariableType& rThisVariable);
            /**
             * Operators 
             */
            
            
            /**
             * Operations
             */

            /**
             * Calculates the StressTensor and the algorithmic tangent Matrix for a given 
             * Elastic Left Cauchy Green Tensor in trial state after Simo 
             *[Comp. Meth. in Appl. Mech. and Eng. 99 (1992) 61-112] 
             * @param StressTensor 3times3 Kirchhoff Stress Tensor
             * @param LeftCauchyGreen_Trial elastic Left Cauchy Green Tensor in trial state =
             * delta[f_(n+1)]*b^e_n*delta[f_(n+1)]^T
             * @param algorithmicTangent \frac{\delta \tau_{n+1}}{\delta b_{n+1}}
             */
    void CalculateStressAndTangentMatrix(Matrix& StressTensor, 
            const Matrix& StrainTensor, 
            array_1d<double,81>& algorithmicTangent);
  

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
            
			Matrix mElasticMaterialMatrix28days;
//             array_1d<double,81> mElasticMaterialTensor28days;
            double mE;
            double mTe;
            double mDTe;
            double mRatioE1E28;
            double mAe;
            double mBe;
            double mCe;
            double mDe;
            double mCurrentTime;
            double mDeltaTime;
            double mCurrentXi;

            Vector mlogEpsilon_n;
            Vector mlogEpsilon_current;
            Vector mlogEpsilon_t_n;
//             Matrix mEpsilon_n;
//             Matrix mEpsilon_current;
//             Matrix mEpsilon_t_n;

		    Matrix mStressTensor;
            Matrix mInsituStressTensor;

            Matrix mUpdatedLeftCauchyGreenTensor;
            Matrix mLeftCauchyGreenTensor_Old;
            Vector mMaterialParameters;

			void SetVariable( const Variable<Matrix>& rVariable, Matrix& rValue);
            /**
             * Operations
             */

             /**
             * Maps the algorithmic matrix, the henky strains and the kirchhoff stress tensor in 
             * principal state back to the actual configuration
             * @param  aep algorithmic tangent principal state
             * @param  principalStresses principal kirchhoff stresses
             * @param  stretches stretches in principal state
             * @param  henky strains in principal state
             * @param  LeftCauchyGreenTensor trial elastic Left Cauchy Green Tensor
             * @param  tanC algorithmic tangent in actual configuration
             * @param  UpdatedLeftCauchyGreenTensor elastic Left Cauchy Green Tensor
             */

            void InitializeMaterialDummy
                    (std::vector<std::vector<Matrix> >& C);

            double CalculateXi();
             /**
             * Spectral decomposition of the LeftCachyGreenTensor
             * @param LeftCauchyGreenTensor elastic Left Cauchy Green Tensor in trial state =
             * delta[f_(n+1)]*b^e_n*delta[f_(n+1)]^T
             * @param stretches principal stretches
             * @param henky henky strains in principal state
             */
            
            void SpectralDecomposition
                    (const Matrix& LeftCauchyGreenTensor, 
                     Vector& stretches, Vector& henky);
             /**
             * Perturbates the principal stretches if they are equal
             * @param  stretches principal stretches
             * @return perturbated principal stretches
             */
            
            Vector PerturbateLambda(const Vector& stretches);
            
             /**
             * Calculates the principal stresses
             * @param  principalStresses principal stresses
             * @param  aep algorithmic tangent in principal state
             * @param  logStrains henky strains in principal state
             */
            
            void CalculateStressAndTangentialStiffness_PrincipalState
                    (Vector& principalStresses, Matrix& aep,const Vector& logStrains);

             /**
             * Maps the algorithmic matrix, the henky strains and the kirchhoff stress tensor in 
             * principal state back to the actual configuration
             * @param  aep algorithmic tangent principal state
             * @param  principalStresses principal kirchhoff stresses
             * @param  stretches stretches in principal state
             * @param  henky strains in principal state
             * @param  LeftCauchyGreenTensor trial elastic Left Cauchy Green Tensor
             * @param  tanC algorithmic tangent in actual configuration
             * @param  UpdatedLeftCauchyGreenTensor elastic Left Cauchy Green Tensor
             */
            void InverseSpectralDecomposition(const Matrix aep, 
                                              const Vector& principalStresses, 
                                              const Vector& stretches, const Vector& henky, 
                                              const Matrix& LeftCauchyGreenTensor, 
                                              array_1d<double,81>& tanC,  
                                              Matrix& StressTensor, Matrix& UpdatedLeftCauchyGreenTensor);
            //GroutingMortar(const IsotropicPlaneStressWrinklingNew& rOther);
    }; // Class GroutingMortar 
    
}  // namespace Kratos.

#endif // KRATOS_ISOTROPIC_LINEAR_ELASTIC_H_INCLUDED  defined 



