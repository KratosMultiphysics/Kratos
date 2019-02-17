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

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes

#include "includes/define.h"
#include "constitutive_laws/isotropic_damage_implex.h"

#include "includes/constitutive_law.h"

#include "utilities/math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "structural_application.h"
#include "includes/properties.h"

namespace Kratos
{


/**
 * TO BE TESTED!!!
 */
IsotropicDamageIMPLEX::IsotropicDamageIMPLEX()
    : ConstitutiveLaw()
{
}

/**
 * TO BE TESTED!!!
 */
IsotropicDamageIMPLEX::~IsotropicDamageIMPLEX()
{
}

bool IsotropicDamageIMPLEX::Has( const Variable<double>& rThisVariable )
{
    if( rThisVariable == YOUNG_MODULUS || rThisVariable == POISSON_RATIO || rThisVariable == DAMAGE || rThisVariable == ALPHA)
        return true;
    return false;
}

bool IsotropicDamageIMPLEX::Has( const Variable<Vector>& rThisVariable )
{
    if( rThisVariable == STRESSES )
        return true;
    return false;
}

bool IsotropicDamageIMPLEX::Has( const Variable<Matrix>& rThisVariable )
{
    return false;
}

double& IsotropicDamageIMPLEX::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if( rThisVariable == YOUNG_MODULUS )
        return mE;
    if( rThisVariable == POISSON_RATIO )
        return mNU;
    if( rThisVariable == DAMAGE )
        return mD;
    if( rThisVariable == ALPHA )
        return mAlpha;
  
    KRATOS_THROW_ERROR( std::logic_error, "this variable is not supported", "" );
}

Vector& IsotropicDamageIMPLEX::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    if( rThisVariable == STRESSES )
        return mCurrentStress;
    
    KRATOS_THROW_ERROR( std::logic_error, "Vector Variable case not considered", "" );
}

Matrix& IsotropicDamageIMPLEX::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    KRATOS_THROW_ERROR( std::logic_error, "Vector Variable case not considered", "" );
}

void IsotropicDamageIMPLEX::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                                      const ProcessInfo& rCurrentProcessInfo )
{
}

void IsotropicDamageIMPLEX::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                                      const ProcessInfo& rCurrentProcessInfo )
{
    if( rThisVariable == YOUNG_MODULUS )
        mE = rValue;
    if( rThisVariable == POISSON_RATIO )
        mNU = rValue;
}

void IsotropicDamageIMPLEX::SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                                      const array_1d<double, 3>& rValue,
                                      const ProcessInfo& rCurrentProcessInfo )
{
}

void IsotropicDamageIMPLEX::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                                      const ProcessInfo& rCurrentProcessInfo )
{
    if( rThisVariable == STRESSES )
        mCurrentStress = rValue;
}

void IsotropicDamageIMPLEX::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                                      const ProcessInfo& rCurrentProcessInfo )
{
}

void IsotropicDamageIMPLEX::InitializeMaterial( const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues )
{
    mCurrentStress = ZeroVector(6);
    mFt = props[TENSILE_STRENGTH];
    mGF = props[FRACTURE_ENERGY];
    mE  = props[YOUNG_MODULUS];
    mNU = props[POISSON_RATIO];
    mE_0 = mFt / mE;
    ml   = pow( fabs( geom.Volume() ), 0.33333333333 )/2;
    mE_f = mGF / (ml*mFt);
    mD = 0.0;
    mAlpha = mE_0;
    mq = mAlpha;
    mAlpha_old = mAlpha;
    mAlpha_old_old = mAlpha_old;
    
    KRATOS_WATCH(mE_0);
    KRATOS_WATCH(mE_f);
    KRATOS_WATCH(ml);
}

void IsotropicDamageIMPLEX::ResetMaterial( const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues )
{
}

void IsotropicDamageIMPLEX::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues,
        const ProcessInfo& CurrentProcessInfo )
{ 
    // Initialize algorithmic IMPL-EX variables
    mAlpha_old_old = mAlpha_old;
    mAlpha_old = mAlpha;
    mq_old = mq;
    
//     KRATOS_WATCH(mAlpha_old_old);
//     KRATOS_WATCH(mAlpha_old);
//     KRATOS_WATCH(mAlpha);
    
    //==========================================================================================
    //-----------------------------------{ EXPLICIT ALGORITHM }---------------------------------
    //==========================================================================================
    
    // Compute elastic constitutive tensor
    Matrix C = ZeroMatrix(6,6);
    CalculateElasticMatrix(C);
    
    // Explicit extrapolation of the internal variable alpha
    mdAlpha = mAlpha_old - mAlpha_old_old;
    mAlpha_alg = mAlpha_old + mdAlpha;
    
    // Calculate softening function
    double H = SofteningLaw(mAlpha_old);
    
    // Compute stress-like internal variable
    mq_alg = mq_old + H*mdAlpha;
  
    // Compute extrapolated damage variable 
    mDamage_alg = 1.0 - ( mq_alg / mAlpha_alg );
    
    // Compute algorithmic tangent operator
    mC_alg = (1.0 - mDamage_alg)*C; 
}

void IsotropicDamageIMPLEX::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues,
        const ProcessInfo& CurrentProcessInfo )
{
  
}

void IsotropicDamageIMPLEX::FinalizeNonLinearIteration( const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues,
        const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_WATCH("NONLINEAR ITERATION ACTIVATED------------------------------------------------------------------------");
    // This function does not exist, yet...
}


void  IsotropicDamageIMPLEX::CalculateMaterialResponse( const Vector& StrainVector,
        const Matrix& DeformationGradient,
        Vector& StressVector,
        Matrix& AlgorithmicTangent,
        const ProcessInfo& CurrentProcessInfo,
        const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues,
        bool CalculateStresses,
        int CalculateTangent,
        bool SaveInternalVariables )
{
  
    //==========================================================================================
    //----------------------------------{ IMPLICIT ALGORITHM }----------------------------------
    //==========================================================================================
    
    // Compute elastic constitutive tensor
    Matrix C = ZeroMatrix(6,6);
    CalculateElasticMatrix(C);
    
    // Compute effective stresses
    StressVector = prod(C,StrainVector);
    
    // Equivalent strain based on the energy norm
    double eps_eq = sqrt((1.0/mE)*inner_prod(trans(StrainVector),StressVector));
   
    // Compute trial state
    double alpha_trial = mAlpha_old;
    
    // Determine loading function
    double f = eps_eq - alpha_trial;

    // Check loading function and compute damage multiplier
    double dlambda;
    if( f <= 0.0 )
    {
        // Elastic/unloading state
        dlambda = 0.0;
    }
    else
    {
        // Damage/loading state
        dlambda = f;
    }
    
//     KRATOS_WATCH(dlambda);
    
    // Update internal variable
    mAlpha = mAlpha_old + dlambda;
    
    // Calculate softening function
    double H = SofteningLaw(mAlpha_old);
    
    // Update stress-like internal variable
    mq = mq_old + H*dlambda;
    
    // Update damage variable
    mD = 1.0 - (mq/mAlpha);
    
//     KRATOS_WATCH(mD);
//     KRATOS_WATCH(mq);
//     KRATOS_WATCH(mAlpha);
    
//     if (mD > 1.0)
//     {
//        KRATOS_WATCH("DAMAGE VARIABLE IS NOT VALID.");
//        KRATOS_WATCH(mq/mAlpha);
//     }
//     
//     if (mD < 0)
//        KRATOS_WATCH("DAMAGE VARIABLE IS NEGATIVE ...");
    
    // Update corrected stresses
    Matrix C_t = (1.0 - mD)*C;
    noalias(AlgorithmicTangent) = C_t;
    mCurrentStress = prod( C_t, StrainVector );
    
    //==========================================================================================
    //----------------------------------{ EXPLICIT ALGORITHM }----------------------------------
    //==========================================================================================
    
    // Explicit extrapolation of the internal variable alpha
    mdAlpha = mAlpha_old - mAlpha_old_old;
    mAlpha_alg = mAlpha_old + mdAlpha;
    
    // Compute stress-like internal variable
    mq_alg = mq_old + H*mdAlpha;
  
    // Compute extrapolated damage variable 
    mDamage_alg = 1.0 - ( mq_alg / mAlpha_alg );
    
    // Compute algorithmic tangent operator
    mC_alg = (1.0 - mDamage_alg)*C; 
    
    // Set algorithmic tangent operator
    noalias(AlgorithmicTangent) = mC_alg;
    
    // Compute extrapolated stresses
    StressVector = prod( mC_alg, StrainVector );
    
}

void IsotropicDamageIMPLEX::CalculateElasticMatrix( Matrix& C )
{
    // Compute elastic constitutive tensor
    C.resize(6,6,false);
    double c1 = mE / (( 1.00 + mNU ) * ( 1 - 2 * mNU ) );
    double c2 = c1 * ( 1 - mNU );
    double c3 = c1 * mNU;
    double c4 = c1 * 0.5 * ( 1 - 2 * mNU );
    C(0,0) = c2;  C(0,1) = c3;  C(0,2) = c3;  C(0,3) = 0.0; C(0,4) = 0.0; C(0,5) = 0.0;
    C(1,0) = c3;  C(1,1) = c2;  C(1,2) = c3;  C(1,3) = 0.0; C(1,4) = 0.0; C(1,5) = 0.0;
    C(2,0) = c3;  C(2,1) = c3;  C(2,2) = c2;  C(2,3) = 0.0; C(2,4) = 0.0; C(2,5) = 0.0;
    C(3,0) = 0.0; C(3,1) = 0.0; C(3,2) = 0.0; C(3,3) = c4;  C(3,4) = 0.0; C(3,5) = 0.0;
    C(4,0) = 0.0; C(4,1) = 0.0; C(4,2) = 0.0; C(4,3) = 0.0; C(4,4) = c4;  C(4,5) = 0.0;
    C(5,0) = 0.0; C(5,1) = 0.0; C(5,2) = 0.0; C(5,3) = 0.0; C(5,4) = 0.0; C(5,5) = c4;
}

double IsotropicDamageIMPLEX::SofteningLaw( const double& alpha)
{
     double H = -(mE_0/mE_f)*exp(-(alpha-mE_0)/mE_f);
     
     return H;
}

//**********************************************************************
int IsotropicDamageIMPLEX::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY

    return 0;

    KRATOS_CATCH( "" );
}
} // Namespace Kratos
