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
*   Last Modified by:    $Author: kazem $
*   Date:                $Date: 2009-01-16 10:49:53 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/

#if !defined(KRATOS_PLANE_STRAIN_H_INCLUDED )
#define  KRATOS_PLANE_STRAIN_H_INCLUDED

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
class PlaneStrain : public ConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ConstitutiveLaw BaseType;
    /**
     * Counted pointer of PlaneStrain
     */
    KRATOS_CLASS_POINTER_DEFINITION(PlaneStrain);

    /**
     * Life Cycle
     */
    /**
     * Default constructor.
     */
    PlaneStrain();

    virtual  ConstitutiveLaw::Pointer Clone() const
    {
         ConstitutiveLaw::Pointer p_clone(new PlaneStrain());
        return p_clone;
    }

    /**
     * Destructor.
     */
    virtual ~PlaneStrain();

    /**
     * Operators
     */
    /**
     * Operations
     */
    bool Has( const Variable<int>& rThisVariable );
    bool Has( const Variable<double>& rThisVariable );
    bool Has( const Variable<Vector>& rThisVariable );
    bool Has( const Variable<Matrix>& rThisVariable );

    int& GetValue( const Variable<int>& rThisVariable, int& rValue );
    double& GetValue( const Variable<double>& rThisVariable, double& rValue );
    Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue );
    Matrix& GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue );
    std::size_t GetStrainSize();

    void SetValue( const Variable<int>& rVariable,
                   const int& Value,
                   const ProcessInfo& rCurrentProcessInfo );
    void SetValue( const Variable<double>& rVariable,
                   const double& Value,
                   const ProcessInfo& rCurrentProcessInfo );
    void SetValue( const Variable<Vector>& rThisVariable,
                   const Vector& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
    void SetValue( const Variable<Matrix>& rThisVariable,
                   const Matrix& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
    /**
     * Material parameters are inizialized
     */
    virtual void InitializeMaterial( const Properties& props,
                                     const GeometryType& geom,
                                     const Vector& ShapeFunctionsValues );

    virtual void InitializeNonLinearIteration( const Properties& rMaterialProperties,
                                               const GeometryType& rElementGeometry,
                                               const Vector& rShapeFunctionsValues,
                                               const ProcessInfo& rCurrentProcessInfo );

    virtual void FinalizeNonLinearIteration( const Properties& rMaterialProperties,
                                             const GeometryType& rElementGeometry,
                                             const Vector& rShapeFunctionsValues,
                                             const ProcessInfo& rCurrentProcessInfo );

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
    /*            void InitializeSolutionStep( const Properties& props,
                        const GeometryType& geom, //this is just to give the array of nodes
                        const Vector& ShapeFunctionsValues ,
                        const ProcessInfo& CurrentProcessInfo);

                void FinalizeSolutionStep( const Properties& props,
                        const GeometryType& geom, //this is just to give the array of nodes
                        const Vector& ShapeFunctionsValues ,
                        const ProcessInfo& CurrentProcessInfo);
    */
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


    /**
     * converts a strain vector styled variable into its form, which the
     * deviatoric parts are no longer multiplied by 2
     */
    void Calculate(const Variable<Matrix >& rVariable, Matrix& rResult,
                   const ProcessInfo& rCurrentProcessInfo);

    void  CalculateMaterialResponse( const Vector& StrainVector,
                                     const Matrix& DeformationGradient,
                                     Vector& StressVector,
                                     Matrix& AlgorithmicTangent,
                                     const ProcessInfo& CurrentProcessInfo,
                                     const Properties& props,
                                     const GeometryType& geom,
                                     const Vector& ShapeFunctionsValues,
                                     bool CalculateStresses = true,
                                     int CalculateTangent = 1,
                                     bool SaveInternalVariables = true );

    int Check(const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo);

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

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
        rSerializer.save( "mE", mE );
        rSerializer.save( "mNU", mNU );
        rSerializer.save( "mDE", mDE );
        rSerializer.save( "mCurrentStress", mCurrentStress );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
        rSerializer.load( "mE", mE );
        rSerializer.load( "mNU", mNU );
        rSerializer.load( "mDE", mDE );
        rSerializer.load( "mCurrentStress", mCurrentStress );
    }

    /**
     * Static Member Variables
     */

    /**
     * calculates the linear elastic constitutive matrix in terms of Young's modulus and
     * Poisson ratio
     * @param E the Young's modulus
     * @param NU the Poisson ratio
     * @return the linear elastic constitutive matrix
     */
    void CalculateElasticMatrix(Matrix& C, const double E, const double NU);

    double mE, mNU, mDE;
    Vector mCurrentStress;



    /**
     * Un accessible methods
     */
    /**
     * Assignment operator.
     */
    //PlaneStrain& operator=(const IsotropicPlaneStressWrinklingNew& rOther);
    /**
     * Copy constructor.
     */
    //PlaneStrain(const IsotropicPlaneStressWrinklingNew& rOther);
}; // Class PlaneStrain
}  // namespace Kratos.
#endif // KRATOS_ISOTROPIC_2D_H_INCLUDED  defined 
