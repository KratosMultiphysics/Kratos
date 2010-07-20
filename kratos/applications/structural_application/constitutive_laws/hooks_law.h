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
*   Last Modified by:    $Author: hurga $
*   Date:                $Date: 2009-03-05 12:01:22 $
*   Revision:            $Revision: 1.5 $
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
    class HooksLaw : public ConstitutiveLaw
    {
        public:
            /**
             * Type Definitions
             */
            typedef ConstitutiveLaw BaseType;
//             typedef BaseType::MaterialTensorType MaterialTensorType;
            
            /**
             * Counted pointer of HooksLaw
             */
            typedef boost::shared_ptr<HooksLaw> Pointer;
            
            
            /**
             * Life Cycle 
             */
            
            /**
             * Default constructor.
             */
            HooksLaw();
            
            /**
             * Constructor
             * @param E the young's modulus of the specified material
             * @param NU the poisson ratio of the specified material
             */
            HooksLaw(const double E, const double NU);
            
            /**
             * Destructor.
             */
            virtual ~HooksLaw();
            
            bool Has( const Variable<double>& rThisVariable );
            bool Has( const Variable<Vector>& rThisVariable );
            bool Has( const Variable<Matrix>& rThisVariable );
            
            double& GetValue( const Variable<double>& rThisVariable, double& rValue );
            Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue );
            
            void SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                              const ProcessInfo& rCurrentProcessInfo );
            void SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                              const ProcessInfo& rCurrentProcessInfo );
            
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
             * @param StressVector the calculated stress vector 
             * @param StrainVector the given strain vector
             * @param algorithmicTangent the calculated algorithmic tangent matrix
             */
             void UpdateMaterial(  const Vector& StrainVector,
                                      const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues,
                                      const ProcessInfo& CurrentProcessInfo );
                                      
             void CalculateStress(const Vector& StrainVector, Vector& StressVector);
    
             void CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult);
    
             void CalculateStressAndTangentMatrix( Vector& StressVector,
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
            
            double mE;
            double mNU;
            Matrix mC;
            Vector mMaterialParameters;
            Vector mCurrentStress;
            Vector mInsituStress;

            /**
             * Operations
             */

            
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


            //HooksLaw(const IsotropicPlaneStressWrinklingNew& rOther);
    }; // Class HooksLaw 
    
}  // namespace Kratos.

#endif // KRATOS_ISOTROPIC_LINEAR_ELASTIC_H_INCLUDED  defined 


