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
 *   Date:                $Date: 2008-01-25 08:37:31 $
 *   Revision:            $Revision: 1.10 $
 *
 * ***********************************************************/

#if !defined(KRATOS_DUMMY_CONSTITUTIVE_LAW_H_INCLUDED )
#define  KRATOS_DUMMY_CONSTITUTIVE_LAW_H_INCLUDED

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

    class DummyConstitutiveLaw : public ConstitutiveLaw
    {
        public:
            /**
             * Type Definitions
             */
            typedef ConstitutiveLaw BaseType;
            /**
             * Counted pointer of DummyConstitutiveLaw
             */
            typedef boost::shared_ptr<DummyConstitutiveLaw> Pointer;

            /**
             * Life Cycle
             */
            /**
             * Default constructor.
             */
            DummyConstitutiveLaw();

            virtual boost::shared_ptr<ConstitutiveLaw> Clone() const
            {
                boost::shared_ptr<ConstitutiveLaw> p_clone( new DummyConstitutiveLaw() );
                return p_clone;
            }

            /**
             * Destructor.
             */
            virtual ~DummyConstitutiveLaw();

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
            void SetValue( const Variable<array_1d<double, 3 > >& rThisVariable,
                           const array_1d<double, 3 > & rValue, const ProcessInfo& rCurrentProcessInfo );
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
             * As this constitutive law describes only linear elastic material properties
             * this function is rather useless and in fact does nothing
             */
            void InitializeSolutionStep( const Properties& props,
                                         const GeometryType& geom, //this is just to give the array of nodes
                                         const Vector& ShapeFunctionsValues,
                                         const ProcessInfo& CurrentProcessInfo );

            void ResetMaterial( const Properties& props,
                                const GeometryType& geom,
                                const Vector& ShapeFunctionsValues );

            void FinalizeSolutionStep( const Properties& props,
                                       const GeometryType& geom, //this is just to give the array of nodes
                                       const Vector& ShapeFunctionsValues,
                                       const ProcessInfo& CurrentProcessInfo );


            /**
             * This function is designed to be called once to perform all the checks needed
             * on the input provided. Checks can be "expensive" as the function is designed
             * to catch user's errors.
             * @param props
             * @param geom
             * @param CurrentProcessInfo
             * @return
             */
            virtual int Check( const Properties& props,
                               const GeometryType& geom,
                               const ProcessInfo& CurrentProcessInfo );

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
             * returns the size of the strain vector of the current constitutive law
             * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
             */
            virtual SizeType GetStrainSize()
            {
                return 6;
            }

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

            virtual void save( Serializer& rSerializer ) const
            {
                KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw );
            }

            virtual void load( Serializer& rSerializer )
            {
                KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw );
            }
            
	    /**
             * Copy constructor.
             */
            //DummyConstitutiveLaw(const IsotropicPlaneStressWrinklingNew& rOther);
    }; // Class DummyConstitutiveLaw


} // namespace Kratos.
#endif // KRATOS_DUMMY_CONSTITUTIVE_LAW_H_INCLUDED  defined 
