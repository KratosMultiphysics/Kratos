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

#if !defined(KRATOS_KONTINUOUS_DAMAGE_3D_H_INCLUDED )
#define  KRATOS_KONTINUOUS_DAMAGE_3D_H_INCLUDED

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
 * Defines an elastic-plastic constitutive law in 3D space with linear isotropic hardening.
 * This material law is defined by the parameters K (Kompression modulus)
 * and MU (Shear modulus).
 */
class ContinuousDamage3D : public ConstitutiveLaw<Node<3> >
{
public:
    /**
     * Type Definitions
     */
    typedef ConstitutiveLaw<Node<3> > BaseType;
    /**
     * Counted pointer of ContinuousDamage3D
     */
    KRATOS_CLASS_POINTER_DEFINITION(ContinuousDamage3D);

    /**
     * Life Cycle
     */
    /**
     * Default constructor.
     */
    ContinuousDamage3D();

    virtual ConstitutiveLaw<Node<3> >::Pointer Clone() const
    {
        ConstitutiveLaw<Node<3> >::Pointer p_clone(new ContinuousDamage3D());
        return p_clone;
    }

    /**
     * Destructor.
     */
    virtual ~ContinuousDamage3D();

    /**
     * Operators
     */
    /**
     * Operations
     */
    bool Has( const Variable<double>& rThisVariable );
    bool Has( const Variable<Vector>& rThisVariable );
    bool Has( const Variable<Matrix>& rThisVariable );

    double GetValue( const Variable<double>& rThisVariable );
    Vector GetValue( const Variable<Vector>& rThisVariable );
    Matrix GetValue( const Variable<Matrix>& rThisVariable );

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

    void CalculateStressAndTangentMatrix(Vector& StressVector,
                                         const Vector& StrainVector,
                                         Matrix& algorithmicTangent);


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

    int Check(const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo);

    /**
     * Calculates the cauchy stresses. For a given deformation and stress state
     * the cauchy stress vector is calculated
     * @param Cauchy_StressVector the result vector of cauchy stresses (will be overwritten)
     * @param F the current deformation gradient
     * @param PK2_StressVector the current second Piola-Kirchhoff-Stress vector
     * @param GreenLagrangeStrainVector the current Green-Lagrangian strains
     */
// 			void CalculateCauchyStresses( Vector& Cauchy_StressVector,
// 					const Matrix& F,
// 					const Vector& PK2_StressVector,
// 					const Vector& GreenLagrangeStrainVector);


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


    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw<Node<3> > );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw<Node<3> >);
    }


    /**
     * Static Member Variables
     */

    /**
     * calculates the linear elastic constitutive matrix in terms of Young's modulus and
     * Poisson ratio
     * @param K the Kompression Modulus
     * @param MU the Shear Modulus
     * @return the linear elastic constitutive matrix
     */
    void CalculateElasticMatrix(Matrix& C,  double K,  double MU);

    double mK,mMU,mDINF,mDSO,mALPHA,mALPHA0,mD,mDPartial;
    Matrix mC;
    Matrix mCtangent;
    Vector mInSituStress;
    Vector mCurrentStress;
    Vector mdevStrainVector;
    Vector mStressVector;


    /**
     * Un accessible methods
     */
    /**
     * Assignment operator.
     */
    //ContinuousDamage3D& operator=(const VonMisesPlaneStressWrinklingNew& rOther);
    /**
     * Copy constructor.
     */
    //ContinuousDamage3D(const VonMisesPlaneStressWrinklingNew& rOther);
}; // Class ContinuousDamage3D
}  // namespace Kratos.
#endif // KRATOS_VON_MISES_BILINEAR_3D_H_INCLUDED  defined 
