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

#if !defined(KRATOS_DRUCKER_PRAGER_H_INCLUDED )
#define  KRATOS_DRUCKER_PRAGER_H_INCLUDED

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
     * Drucker-Prager plasticity model with linear hardening
     */

    class DruckerPrager : public ConstitutiveLaw
    {

        public:
            
            typedef ConstitutiveLaw BaseType;

            /**
             * Counted pointer of DruckerPrager
             */
            KRATOS_CLASS_POINTER_DEFINITION( DruckerPrager );

            /**
             * Life Cycle
             */
            DruckerPrager();

            /**
             * Destructor.
             */
            virtual ~DruckerPrager();

            /**
             * Clone function
             * @return a pointer to a new instance of this constitutive law
             */
            virtual boost::shared_ptr<BaseType> Clone() const
            {
                boost::shared_ptr<BaseType> p_clone( new DruckerPrager() );
                return p_clone;
            }

            /**
             * @return the working space dimension of the current constitutive law
             */
            virtual SizeType WorkingSpaceDimension()
            {
                return( 3 );
            }

            /**
             * returns the size of the strain vector of the current constitutive law
             */
            virtual SizeType GetStrainSize()
            {
                return( 6 );
            }

            /**
             * returns whether this constitutive Law has specified variable
             * @param rThisVariable the variable to be checked for
             * @return true if the variable is defined in the constitutive law
             */
            bool Has( const Variable<double>& rThisVariable );

            /**
             * returns whether this constitutive Law has specified variable
             * @param rThisVariable the variable to be checked for
             * @return true if the variable is defined in the constitutive law
             */
            bool Has( const Variable<Vector>& rThisVariable );

            /**
             * returns whether this constitutive Law has specified variable
             * @param rThisVariable the variable to be checked for
             * @return true if the variable is defined in the constitutive law
             */
            bool Has( const Variable<Matrix>& rThisVariable );

            /**
             * returns whether this constitutive Law has specified variable
             * @param rThisVariable the variable to be checked for
             * @return true if the variable is defined in the constitutive law
             * NOTE: fixed size array of 3 doubles (e.g. for 2D stresses, plastic strains, ...)
             */
            bool Has( const Variable<array_1d<double, 3 > >& rThisVariable );

            /**
             * returns whether this constitutive Law has specified variable
             * @param rThisVariable the variable to be checked for
             * @return true if the variable is defined in the constitutive law
             * NOTE: fixed size array of 6 doubles (e.g. for stresses, plastic strains, ...)
             */
            bool Has( const Variable<array_1d<double, 6 > >& rThisVariable );

            /**
             * returns the value of a specified variable
             * @param rThisVariable the variable to be returned
             * @param rValue a reference to the returned value
             * @param rValue output: the value of the specified variable
             */
            virtual double& GetValue( const Variable<double>& rThisVariable, double& rValue );

            /**
             * returns the value of a specified variable
             * @param rThisVariable the variable to be returned
             * @param rValue a reference to the returned value
             * @return the value of the specified variable
             */
            virtual Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue );

            /**
             * returns the value of a specified variable
             * @param rThisVariable the variable to be returned
             * @return the value of the specified variable
             */
            virtual Matrix& GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue );

            /**
             * returns the value of a specified variable
             * @param rThisVariable the variable to be returned
             * @param rValue a reference to the returned value
             * @return the value of the specified variable
             */
            virtual array_1d<double, 3 > & GetValue( const Variable<array_1d<double, 3 > >& rVariable,
                    array_1d<double, 3 > & rValue );

            /**
             * returns the value of a specified variable
             * @param rThisVariable the variable to be returned
             * @param rValue a reference to the returned value
             * @return the value of the specified variable
             */
            virtual array_1d<double, 6 > & GetValue( const Variable<array_1d<double, 6 > >& rVariable,
                    array_1d<double, 6 > & rValue );

            /**
             * sets the value of a specified variable
             * @param rVariable the variable to be returned
             * @param Value new value of the specified variable
             * @param rCurrentProcessInfo the process info
             */
            virtual void SetValue( const Variable<double>& rVariable,
                                   const double& Value,
                                   const ProcessInfo& rCurrentProcessInfo );

            /**
             * sets the value of a specified variable
             * @param rVariable the variable to be returned
             * @param Value new value of the specified variable
             * @param rCurrentProcessInfo the process info
             */
            virtual void SetValue( const Variable<Vector >& rVariable,
                                   const Vector& Value, const ProcessInfo& rCurrentProcessInfo );

            /**
             * sets the value of a specified variable
             * @param rVariable the variable to be returned
             * @param Value new value of the specified variable
             * @param rCurrentProcessInfo the process info
             */
            virtual void SetValue( const Variable<Matrix >& rVariable,
                                   const Matrix& Value, const ProcessInfo& rCurrentProcessInfo );

            /**
             * sets the value of a specified variable
             * @param rVariable the variable to be returned
             * @param Value new value of the specified variable
             * @param rCurrentProcessInfo the process info
             */
            virtual void SetValue( const Variable<array_1d<double, 3 > >& rVariable,
                                   const array_1d<double, 3 > & Value,
                                   const ProcessInfo& rCurrentProcessInfo );

            /**
             * sets the value of a specified variable
             * @param rVariable the variable to be returned
             * @param Value new value of the specified variable
             * @param rCurrentProcessInfo the process info
             */
            virtual void SetValue( const Variable<array_1d<double, 6 > >& rVariable,
                                   const array_1d<double, 6 > & Value,
                                   const ProcessInfo& rCurrentProcessInfo );

            /**
             * Is called to check whether the provided material parameters in the Properties
             * match the requirements of current constitutive model.
             * @param props the current Properties to be validated against.
             * @return true, if parameters are correct; false, if parameters are insufficient / faulty
             * NOTE: this has to implemented by each constitutive model. Returns false in base class since
             * no valid implementation is contained here.
             */
            virtual bool ValidateInput( const Properties& props );

            /**
             * returns the expected strain measure of this constitutive law (by default linear strains)
             * @return the expected strain measure
             */
            virtual StrainMeasure GetStrainMeasure()
            {
                return StrainMeasure_Linear;
            }

            /**
             * returns the stress measure of this constitutive law (by default 1st Piola-Kirchhoff stress in voigt notation)
             * @return the expected stress measure
             */
            virtual StressMeasure GetStressMeasure()
            {
                return StressMeasure_PK1;
            }

            /**
             * returns whether this constitutive model is formulated in incremental strains/stresses
             * NOTE: by default, all constitutive models should be formulated in total strains
             * @return true, if formulated in incremental strains/stresses, false otherwise
             */
            virtual bool IsIncremental()
            {
                return false;
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
                                             const Vector& ShapeFunctionsValues );

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
                                                 const Vector& ShapeFunctionsValues,
                                                 const ProcessInfo& CurrentProcessInfo );

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
                                               const Vector& ShapeFunctionsValues,
                                               const ProcessInfo& CurrentProcessInfo );

            virtual void InitializeNonLinearIteration( const Properties& props,
                    const GeometryType& geom,
                    const Vector& ShapeFunctionsValues,
                    const ProcessInfo& CurrentProcessInfo );

            /**
             * Computes the material response in terms of stresses and algorithmic tangent
             * @param StrainVector the current strains (total strains, input)
             * @param DeformationGradient the current deformation gradient (can be an empty matrix if a linear strain measure is used)
             * @param StressVector the computed stresses (output)
             * @param algorithmicTangent the material tangent matrix (output)
             * @param CurrentProcessInfo current ProcessInfo instance
             * @param props the material's Properties object
             * @param geom the element's geometry
             * @param ShapeFunctionsValues the shape functions values in the current integration pointer
             * @param CalculateStresses flag whether or not to compute the stress response
             * @param CalculateTangent flag to determine if to compute the material tangent
             * NOTE: the CalculateTangent flag is defined as int to allow for distinctive variants of the tangent
             * @param SaveInternalVariables flag whether or not to store internal (history) variables
             */
            virtual void CalculateMaterialResponse( const Vector& StrainVector,
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
             * Computes the volumetric part of the material response in terms of stresses and algorithmic tangent
             * @param StrainVector the current strains (total strains, input)
             * @param DeformationGradient the current deformation gradient (can be an empty matrix if a linear
             * @param StressVector the computed stresses (output)
             * @param algorithmicTangent the material tangent matrix (output)
             * @param CurrentProcessInfo current ProcessInfo instance
             * @param props the material's Properties object
             * @param geom the element's geometry
             * @param ShapeFunctionsValues the shape functions values in the current integration pointer
             * @param CalculateStresses flag whether or not to compute the stress response
             * @param CalculateTangent flag to determine if to compute the material tangent
             * NOTE: the CalculateTangent flag is defined as int to allow for distinctive variants of the tangent
             * @param SaveInternalVariables flag whether or not to store internal (history) variables
             */
            virtual void CalculateVolumetricResponse( const double VolumetricStrain,
                    const Matrix& DeformationGradient,
                    double& VolumetricStress,
                    double& AlgorithmicBulk,
                    const ProcessInfo& CurrentProcessInfo,
                    const Properties& props,
                    const GeometryType& geom,
                    const Vector& ShapeFunctionsValues,
                    bool CalculateStresses = true,
                    int CalculateTangent = true,
                    bool SaveInternalVariables = true )
            {
                KRATOS_ERROR(std::logic_error, "Volumetric response is not implemented in DruckerPrager", "");
            }

            /**
             * Computes the deviatoric part of the material response in terms of stresses and algorithmic tangent
             * @param StrainVector the current strains (total strains, input)
             * @param DeformationGradient the current deformation gradient (can be an empty matrix if a linear
             * @param StressVector the computed stresses (output)
             * @param algorithmicTangent the material tangent matrix (output)
             * @param CurrentProcessInfo current ProcessInfo instance
             * @param props the material's Properties object
             * TODO: add proper definition for algorithmic tangent
             */
            virtual void CalculateDeviatoricResponse( const Vector& StrainVector,
                    const Matrix& DeformationGradient,
                    Vector& StressVector,
                    Matrix& AlgorithmicTangent,
                    const ProcessInfo& CurrentProcessInfo,
                    const Properties& props,
                    const GeometryType& geom,
                    const Vector& ShapeFunctionsValues,
                    bool CalculateStresses = true,
                    int CalculateTangent = true,
                    bool SaveInternalVariables = true )
            {
                KRATOS_ERROR(std::logic_error, "Deviatoric response is not implemented in DruckerPrager", "");
            }

            /**
             * This can be used in order to reset all internal variables of the
             * constitutive law (e.g. if a model should be reset to its reference state)
             * @param props the Properties instance of the current element
             * @param geom the geometry of the current element
             * @param ShapeFunctionsValues the shape functions values in the current integration point
             * @param the current ProcessInfo instance
             */
            virtual void ResetMaterial( const Properties& props,
                                        const GeometryType& geom,
                                        const Vector& ShapeFunctionsValues );

            virtual void CalculateCauchyStresses( Vector& Cauchy_StressVector,
                                                  const Matrix& F,
                                                  const Vector& PK2_StressVector,
                                                  const Vector& GreenLagrangeStrainVector )
            {
                KRATOS_ERROR(std::logic_error, "Cauchy stresses are not implemented in DruckerPrager", "");
            }

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

        protected:

            /**
             * Member Variables
             */
            Vector mCurrentPlasticStrains;
            Vector mPrestress;
            double mAlpha;
            double mOldAlpha;
            Matrix mCtangent;
            Vector mCurrentStress;

            double mDensity, mE, mNU, mG, mK, mEta, mXi, mCohesion, mHardening, mPrestressFactor;
            double mdGamma;
            
            

            double unit4thSym3D[6][6];
            double unit2nd3D[6];
            Matrix mConsistentTangent;
            
            Vector mOldStrain;
            Vector mCurrentStrain;
            Vector mOldSElasticStrain;
            Vector mCurrentElasticStrain;
            Vector mCurrentStrainInc;
            Vector mOldStress;
            
            double mModelData[10];
            Vector mOldPlasticStrains;
            
            
            bool isYielded;


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
            
            void CalculateStress( const Vector& StrainVector, Vector& StressVector );
            
            void CalculateConstitutiveMatrix( const Vector& StrainVector, Matrix& rResult );

            void CalculateElasticMatrix( Matrix& C, const double E, const double NU );

            void CalculateUnit4thSym3D();

            void CalculateUnit2nd3D();


    }; // Class DruckerPrager
}  // namespace Kratos.
#endif // KRATOS_DRUCKER_PRAGER_H_INCLUDED  defined
