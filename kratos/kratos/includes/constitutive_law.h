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
*   Last Modified by:    $Author: nagel $
*   Date:                $Date: 2009-01-05 13:22:39 $
*   Revision:            $Revision: 1.4 $
*
* ***********************************************************/

#if !defined(KRATOS_CONSTITUTIVE_LAW )
#define  KRATOS_CONSTITUTIVE_LAW

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/node.h"
#include "includes/properties.h"
#include "geometries/geometry.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "containers/data_value_container.h"


namespace Kratos
{

    template< class TPointType = Node<3> > class ConstitutiveLaw
    {
        public:
            typedef ProcessInfo ProcessInfoType;
            typedef Geometry<TPointType> GeometryType;
            /**
             * Type of 4th order material tensor
             * needed for tensorial formulation of constitutive laws
             */
            typedef array_1d<double,81> MaterialTensorType;
            
            /** 
             * Counted pointer of ConstitutiveLaw 
             */
            KRATOS_CLASS_POINTER_DEFINITION( ConstitutiveLaw<TPointType> );
            
            /** 
             * Constructor.
             */
            ConstitutiveLaw()
            {
            }
            
            /** 
             * Destructor.
             */
            virtual ~ConstitutiveLaw()
            {
            }
            
            /**
             * Clone function
             */
            virtual boost::shared_ptr<ConstitutiveLaw<TPointType> > Clone() const
            {
                boost::shared_ptr<ConstitutiveLaw<TPointType> > p_clone(new ConstitutiveLaw());
                return p_clone;
            }
            
            /**
             * returns whether this constitutive Law has specified variable
             * @param rThisVariable the variable to be checked for
             * @return true if the variable is defined in the constitutive law
             */
            virtual bool Has( const Variable<double>& rThisVariable )
            {
                KRATOS_ERROR(std::logic_error,  "Called the virtual function for Has" , "");
            }
            /**
             * returns whether this constitutive Law has specified variable
             * @param rThisVariable the variable to be checked for
             * @return true if the variable is defined in the constitutive law
             */
            virtual bool Has( const Variable<Vector>& rThisVariable )
            {
                KRATOS_ERROR(std::logic_error,  "Called the virtual function for Has" , "");
            }
            /**
             * returns whether this constitutive Law has specified variable
             * @param rThisVariable the variable to be checked for
             * @return true if the variable is defined in the constitutive law
             */
            virtual bool Has( const Variable<Matrix>& rThisVariable )
            {
                KRATOS_ERROR(std::logic_error,  "Called the virtual function for Has" , "");
            }
            /**
             * returns the value of a specified variable
             * @param rThisVariable the variable to be returned
             * @return the value of the specified variable
             */
            virtual double GetValue( const Variable<double>& rThisVariable )
            {
                KRATOS_ERROR(std::logic_error,  "Called the virtual function for GetValue" , "");
            }
            /**
             * returns the value of a specified variable
             * @param rThisVariable the variable to be returned
             * @return the value of the specified variable
             */
            virtual Vector GetValue( const Variable<Vector>& rThisVariable )
            {
                KRATOS_ERROR(std::logic_error,  "Called the virtual function for GetValue" , "");
            }
            /**
             * returns the value of a specified variable
             * @param rThisVariable the variable to be returned
             * @return the value of the specified variable
             */
            virtual Matrix GetValue( const Variable<Matrix>& rThisVariable )
            {
                KRATOS_ERROR(std::logic_error,  "Called the virtual function for GetValue" , "");
            }
            /**
             * sets the value of a specified variable
             * @param rVariable the variable to be returned
             * @param Value new value of the specified variable
             * @param rCurrentProcessInfo the process info
             */
            virtual void SetValue( const Variable<double>& rVariable, 
                                   const double Value, 
                                   const ProcessInfo& rCurrentProcessInfo)
            {
                KRATOS_ERROR(std::logic_error,  "Called the virtual function for SetValue" , "");
            }
            /**
             * sets the value of a specified variable
             * @param rVariable the variable to be returned
             * @param Value new value of the specified variable
             * @param rCurrentProcessInfo the process info
             */
            virtual void SetValue( const Variable<array_1d<double,3> >& rVariable, 
                                   const array_1d<double,3>& Value, 
                                   const ProcessInfo& rCurrentProcessInfo)
            {
                KRATOS_ERROR(std::logic_error,  "Called the virtual function for SetValue" , "");
            }
            /**
             * sets the value of a specified variable
             * @param rVariable the variable to be returned
             * @param Value new value of the specified variable
             * @param rCurrentProcessInfo the process info
             */
            virtual void SetValue( const Variable<Vector >& rVariable, 
                                   const Vector& Value, 
                                   const ProcessInfo& rCurrentProcessInfo)
            {
                KRATOS_ERROR(std::logic_error,  "Called the virtual function for SetValue" , "");
            }
            /**
             * sets the value of a specified variable
             * @param rVariable the variable to be returned
             * @param Value new value of the specified variable
             * @param rCurrentProcessInfo the process info
             */
            virtual void SetValue( const Variable<Matrix >& rVariable, 
                                   const Matrix& Value, const ProcessInfo& rCurrentProcessInfo)
            {
                KRATOS_ERROR(std::logic_error,  "Called the virtual function for SetValue" , "");
            }
            /**
             * TODO: to be removed
             */
            template<class TVariableType> bool Has(const TVariableType& rThisVariable) const
            {
                KRATOS_ERROR(std::logic_error,  "Called the virtual function for GetValue" , "");
            }
            
            /**
             * This is to be called at the very beginning of the calculation
             * (e.g. from InitializeElement) in order to initialize all relevant
             * attributes of the constitutive law
             * @param props the Properties instance of the current element
             * @param geom the geometry of the current element
             * @param ShapeFunctionsValues the shape functions values in the current integration point
             */
            virtual void InitializeMaterial( const Properties& props,
                                             const GeometryType& geom,
                                             const Vector& ShapeFunctionsValues )
            {
            }
            
            /**
             * to be called at the beginning of each solution step
             * (e.g. from Element::InitializeSolutionStep)
             * @param props the Properties instance of the current element
             * @param geom the geometry of the current element
             * @param ShapeFunctionsValues the shape functions values in the current integration point
             * @param the current ProcessInfo instance
             */
            virtual void InitializeSolutionStep( const Properties& props,
                    const GeometryType& geom, //this is just to give the array of nodes
                    const Vector& ShapeFunctionsValues ,
                    const ProcessInfo& CurrentProcessInfo)
            {
            }
            
            /**
             * Updates the material (e.g. performs return-map algorithm)
             * This is to be called at each iteration right before the stresses
             * are requested
             * May be left out for linear calculations
             * @param StrainVector the current total strain vector $$\epsilon = \left[\epsilon_{11},
             * \epsilon_{22}, \epsilon_{33}, 2\,\epsilon_{33}, 2\,\epsilon_{33}, 2\,\epsilon_{33}\right]$$
             * @param props the Properties instance of the current element
             * @param geom the geometry of the current element
             * @param ShapeFunctionsValues the values of the shape functions in the current 
             * integration point
             * @param the current ProcessInfo instance
             */
            virtual void UpdateMaterial( const Vector& StrainVector,
                                         const Properties& props,
                                         const GeometryType& geom,
                                         const Vector& ShapeFunctionsValues ,
                                         const ProcessInfo& CurrentProcessInfo)
            {
            }

            /**
             * performs and calculates the stressVector and tangnet Matrix
             * This is to be called at each iteration right before the stresses
             * are requested
             * May be left out for linear calculations
             * @param StressTensor the calculated stress tensor
             * @param StrainTensor the given strain tensor
             * @param algorithmicTangent the 4th order algorithmic tangent tensor
             * @param the current ProcessInfo instance
              */

           virtual void CalculateMaterialResponse(
                            const Vector& StrainVector,
                            //const ProcessInfo& CurrentProcessInfo,
                            Vector& StressVector,
                            Matrix& algorithmicTangent,
                            bool calculate_stress_flag,
                            bool calculate_tangent_flag,
                            bool  save_internal_variables
                            )
          {
               if(calculate_stress_flag == true) this->CalculateStress(StrainVector,StressVector);
               if(calculate_tangent_flag == true) this->CalculateConstitutiveMatrix(StrainVector,algorithmicTangent);
          }

            
            /**
             * to be called at the end of each solution step
             * (e.g. from Element::FinalizeSolutionStep)
             * @param props the Properties instance of the current element
             * @param geom the geometry of the current element
             * @param ShapeFunctionsValues the shape functions values in the current integration point
             * @param the current ProcessInfo instance
             */
            virtual void FinalizeSolutionStep( const Properties& props,
                                               const GeometryType& geom, 
                                               const Vector& ShapeFunctionsValues ,
                                               const ProcessInfo& CurrentProcessInfo)
            {
            }
            
            /**
             * calculates the current stress and the material tangent
             * NOTE: there are two versions of this function: one for a matrix representation
             * of the material tensor and one for a tensorial formulation. Each ConstitutiveLaw
             * HAS TO IMPLEMENT both of them (for convenience, there are conversation functions 
             * available in MathUtils for either of them)
             * @param StressTensor the calculated stress vector 
             */
//             virtual void CalculateStressAndTangentMatrix(Matrix& StressTensor, 
//                     const Matrix& StrainTensor, 
//                     std::vector<std::vector<Matrix> >& algorithmicTangent)
//             {
//             }
            

            /**
             * calculates the current stress and the material tangent
             * NOTE: there are two versions of this function: one for a matrix representation
             * of the material tensor and one for a tensorial formulation. Each ConstitutiveLaw
             * HAS TO IMPLEMENT both of them (for convenience, there are conversation functions 
             * available in MathUtils for either of them)
             * @param StressTensor the calculated stress tensor
             * @param StrainTensor the given strain tensor
             * @param algorithmicTangent the 4th order algorithmic tangent tensor
             */
            virtual void CalculateStressAndTangentMatrix( Matrix& StressTensor,
                    const Matrix& StrainTensor,
                    MaterialTensorType& algorithmicTangent)
            {
                KRATOS_ERROR(std::logic_error,  
                             "Called the virtual function for CalculateStressAndTangentMatrix" , "");
            }
            
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
                    Matrix& algorithmicTangent)
            {
                KRATOS_ERROR(std::logic_error,  
                             "Called the virtual function for CalculateStressAndTangentMatrix" , "");
            }

            /**
             * this is to calculate the Cauchy stresses from a given Piola-Kirchhoff-2 stresses
             * vector
             * @param Cauchy_StressVector output: Cauchy stresses
             * @param F input: deformation gradient
             * @param PK2_StressVector input: Piola-Kirchhoff-2 stresses
             * @param GreenLagrangeStrainVector input Green-Lagrange strains
             */
            virtual void CalculateCauchyStresses( Vector& Cauchy_StressVector,
                    const Matrix& F,
                    const Vector& PK2_StressVector,
                    const Vector& GreenLagrangeStrainVector)
            {
            }
            
	/****/
		 virtual void CalculateStressAndTangentMatrix( Matrix& StressTensor,
			const Matrix& F,
                        const Matrix& StrainTensor,
                        Matrix& algorithmicTangent)
                {
                }


            /**
             * TO BE REMOVED
             */
            virtual void CalculateConstitutiveMatrix( const Vector& StrainVector, 
                    Matrix& ElasticityTensor )
            {
            }
            
            /**
             * calculates the current stresses (in vector representation).
             * This is intended for linear calculations where the algorithmic tangent
             * is not needed
             * @param StrainVector input: current total strains
             * @param StressVector output: the current stresses
             */
            virtual void CalculateStress( const Vector& StrainVector, 
                                          Vector& StressVector)
            {
            }
            
            /**
             * Calculates the value of a given variable on the current integration point
             * @param rVariable input: given variable
             * @param Output output: the calculated value of the variable
             * @param rCurrentProcessInfo input: the current ProcessInfo instance
             */
            virtual void Calculate( const Variable<double>& rVariable, 
                                    double& Output, 
                                    const ProcessInfo& rCurrentProcessInfo)
            {
            }
            
            /**
             * Calculates the value of a given variable on the current integration point
             * @param rVariable input: given variable
             * @param Output output: the calculated value of the variable
             * @param rCurrentProcessInfo input: the current ProcessInfo instance
             */
            virtual void Calculate( const Variable<Vector >& rVariable, 
                                    Vector& Output, 
                                    const ProcessInfo& rCurrentProcessInfo)
            {
            }
            
            /**
             * Calculates the value of a given variable on the current integration point
             * @param rVariable input: given variable
             * @param Output output: the calculated value of the variable
             * @param rCurrentProcessInfo input: the current ProcessInfo instance
             */
            virtual void Calculate( const Variable<Matrix >& rVariable, 
                                    Matrix& Output, 
                                    const ProcessInfo& rCurrentProcessInfo)
            {
            }
            
            /**
             * This can be used in order to reset all internal variables of the
             * constitutive law (e.g. if a model should be reset to its reference state)
             */
            virtual void ResetMaterial()
            {
            }
        
        protected:
 
 
        
        private:
    }; /* Class ConstitutiveLaw */
    
    /**
     * definition of CONSTITUTIVE_LAW variable
     */
    KRATOS_DEFINE_VARIABLE(ConstitutiveLaw<Node<3> >::Pointer, CONSTITUTIVE_LAW)
}  /* namespace Kratos.*/
#endif /* KRATOS_CONSTITUTIVE_LAW  defined */
