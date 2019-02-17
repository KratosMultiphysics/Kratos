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
 *   Last Modified by:    $Author: seyedali $
 *   Date:                $Date: 2014-05-15 00:00:00 $
 *   Revision:            $Revision: 0.54 $
 *
 * ***********************************************************/

// Remark: with reference to J. Oliver et al. (2008) and the 
// Java version of IMPL-EX algorithm by hbui

#if !defined(KRATOS_ISOTROPIC_DAMAGE_IMPLEX_H_INCLUDED )
#define  KRATOS_ISOTROPIC_DAMAGE_IMPLEX_H_INCLUDED

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
 * Defines an IMPL-EX isotropic damage constitutive law in 3D space.
 * This material law is defined by the parameters E (Young's modulus), 
 * NU (Poisson's ratio), Ft (Tensile Strength) and GF (Fracture Energy).
 */

class IsotropicDamageIMPLEX : public ConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ConstitutiveLaw BaseType;
    /**
     * Counted pointer of Isotropic_Damage_Implex
     */
    typedef boost::shared_ptr<IsotropicDamageIMPLEX> Pointer;

    /**
     * Life Cycle
     */
    /**
     * Default constructor.
     */
    IsotropicDamageIMPLEX();

    virtual boost::shared_ptr<ConstitutiveLaw> Clone() const
    {
        boost::shared_ptr<ConstitutiveLaw> p_clone( new IsotropicDamageIMPLEX() );
        return p_clone;
    }

    /**
     * Destructor.
     */
    virtual ~IsotropicDamageIMPLEX();

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

    void SetValue( const Variable<int>& rThisVariable, const int& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
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
                                 const GeometryType& geom,
                                 const Vector& ShapeFunctionsValues,
                                 const ProcessInfo& CurrentProcessInfo );

    void ResetMaterial( const Properties& props,
                        const GeometryType& geom,
                        const Vector& ShapeFunctionsValues );

    void FinalizeSolutionStep( const Properties& props,
                               const GeometryType& geom,
                               const Vector& ShapeFunctionsValues,
                               const ProcessInfo& CurrentProcessInfo );

    void FinalizeNonLinearIteration( const Properties& props,
                               const GeometryType& geom,
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
     * Calculates the elastic constitutive tensor
     * @param rResult elastic tangent operator
     */
    void CalculateElasticMatrix( Matrix& C );
    
    /**
     * Calculates the softening law based on an internal variable alpha
     * @param rResult softening function H
     */
    double SofteningLaw( const double& alpha );

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

    double mFt, mGF, mE, mNU, mE_0, ml, mE_f, mD;
    Vector mCurrentStress;
    double mAlpha, mAlpha_old, mAlpha_old_old, mdAlpha, mAlpha_alg;
    double mq, mq_old, mq_alg;
    double mDamage_alg;
    Matrix mC_alg;
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
}; // Class IsotropicDamageIMPLEX


} // namespace Kratos.
#endif // KRATOS_ISOTROPIC_DAMAGE_IMPLEX_H_INCLUDED  defined 
