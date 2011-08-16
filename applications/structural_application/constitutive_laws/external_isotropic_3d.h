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
*   Date:                $Date: 2008-01-24 16:48:13 $
*   Revision:            $Revision: 1.2 $
*
* ***********************************************************/

#if !defined(KRATOS_EXTERNAL_ISOTROPIC_3D_H_INCLUDED )
#define  KRATOS_EXTERNAL_ISOTROPIC_3D_H_INCLUDED

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
    class ExternalIsotropic3D : public ConstitutiveLaw
    {
        public:
            /**
             * Type Definitions
             */
            typedef ConstitutiveLaw BaseType;
            
            /**
             * Counted pointer of ExternalIsotropic3D
             */
            typedef boost::shared_ptr<ExternalIsotropic3D> Pointer;
            
            /**
             * Life Cycle 
             * 
             * Default constructor.
             */
            ExternalIsotropic3D();
			
            virtual boost::shared_ptr<ConstitutiveLaw> Clone() const
            {
                boost::shared_ptr<ConstitutiveLaw> p_clone(new ExternalIsotropic3D());
                return p_clone;
            }
            
            /**
             * Destructor.
             */
            virtual ~ExternalIsotropic3D();
            
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
			
            void SetValue( const Variable<double>& rThisVariable, const double& rValue, 
                           const ProcessInfo& rCurrentProcessInfo );
            void SetValue( const Variable<array_1d<double, 3> >& rThisVariable, 
                           const array_1d<double, 3>& rValue, const ProcessInfo& rCurrentProcessInfo );
            void SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                           const ProcessInfo& rCurrentProcessInfo );
            void SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
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
            
            /**
             * As this constitutive law describes only linear elastic material properties
             * this function is rather useless and in fact does nothing
             */
            void InitializeSolutionStep( const Properties& props,
                                         const GeometryType& geom, //this is just to give the array of nodes
                                         const Vector& ShapeFunctionsValues ,
                                         const ProcessInfo& CurrentProcessInfo);
			
            void UpdateMaterial( const Vector& StrainVector,
                                 const Properties& props,
                                 const GeometryType& geom, //this is just to give the array of nodes
                                 const Vector& ShapeFunctionsValues ,
                                 const ProcessInfo& CurrentProcessInfo);
            
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
					  
              int Check(const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo);
            
//            /**
//             * converts a strain vector styled variable into its form, which the
//             * deviatoric parts are no longer multiplied by 2
//             */
            
            /**
             * Input and output
             * 
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
            
        /**
         * there are no protected class members
         */
        //protected:
        
        private:
            
            /**
             * calculates the linear elastic constitutive matrix in terms of Young's modulus and
             * Poisson ratio
             * @param E the Young's modulus
             * @param NU the Poisson ratio
             * @return the linear elastic constitutive matrix
             */
            void CalculateElasticMatrix(Matrix& C, const double E, const double NU);
            
            double mE,mNU;
            Vector mInSituStress;
            Matrix mCtangent;
            Vector mCurrentStress;
                      
	    
	    private:

	    ///@}
	    ///@name Serialization
	    ///@{	
	    friend class Serializer;

	    virtual void save(Serializer& rSerializer) const
	    {
	      rSerializer.save("Name","ExternalIsotropic3D");
	      KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw  );
	    }

	    virtual void load(Serializer& rSerializer)
	    {
	      KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw );
	    }
	    
            /**
             * Material properties to be assigned by Geostatistics
             * Here, these are used for testing purpose mainly. That is why
             * a full list of parameters is provided for a linear elastic material law
             * meaning of parameters:
             * 0 - Young's modulus E
             * 1 - Poisson's ratio NU
             * 2 - Friction angle PHI
             * 3 - Cohesion C
             * 4 - Compressive Strength
             * 5 - Tensile Strength
             * 6 - RMR Variance
             */
            Vector mMaterialParameters;
            
            
            /**
             * Assignment operator.
             */
            //ExternalIsotropic3D& operator=(const IsotropicPlaneStressWrinklingNew& rOther);
            /**
             * Copy constructor.
             */
            //ExternalIsotropic3D(const IsotropicPlaneStressWrinklingNew& rOther);
    }; // Class ExternalIsotropic3D 
}  // namespace Kratos.
#endif // KRATOS_EXTERNAL_ISOTROPIC_3D_H_INCLUDED  defined 
