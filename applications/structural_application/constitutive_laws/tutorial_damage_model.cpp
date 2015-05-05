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
*   Date:                $Date: 2009-01-14 17:14:12 $
*   Revision:            $Revision: 1.13 $
*
* ***********************************************************/

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes

#include "includes/define.h"
#include "constitutive_laws/tutorial_damage_model.h"

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
TutorialDamageModel::TutorialDamageModel()
    : ConstitutiveLaw()
{
}

/**
 * TO BE TESTED!!!
 */
TutorialDamageModel::~TutorialDamageModel()
{
}

bool TutorialDamageModel::Has( const Variable<double>& rThisVariable )
{
    if( rThisVariable == YOUNG_MODULUS || rThisVariable == POISSON_RATIO || rThisVariable == DAMAGE_E0 || rThisVariable == DAMAGE_EF )
        return true;
    return false;
}

bool TutorialDamageModel::Has( const Variable<Vector>& rThisVariable )
{
    if( rThisVariable == STRESSES )
        return true;
    return false;
}

bool TutorialDamageModel::Has( const Variable<Matrix>& rThisVariable )
{
    return false;
}

double& TutorialDamageModel::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if( rThisVariable == YOUNG_MODULUS )
        return mE;
    if( rThisVariable == POISSON_RATIO )
        return mNU;
    if( rThisVariable == DAMAGE_E0 )
        return mE_0;
    if( rThisVariable == DAMAGE_EF )
        return mE_f;

    KRATOS_THROW_ERROR( std::logic_error, "this variable is not supported", "" );
}

Vector& TutorialDamageModel::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    if( rThisVariable == STRESSES )
        return mCurrentStress;

    KRATOS_THROW_ERROR( std::logic_error, "Vector Variable case not considered", "" );
}

Matrix& TutorialDamageModel::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    KRATOS_THROW_ERROR( std::logic_error, "Vector Variable case not considered", "" );
}

void TutorialDamageModel::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                                    const ProcessInfo& rCurrentProcessInfo )
{
    if( rThisVariable == YOUNG_MODULUS )
        mE = rValue;
    if( rThisVariable == POISSON_RATIO )
        mNU = rValue;
    if( rThisVariable == DAMAGE_E0 )
        mE_0 = rValue;
    if( rThisVariable == DAMAGE_EF )
        mE_f = rValue;

}

void TutorialDamageModel::SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                                    const array_1d<double, 3>& rValue,
                                    const ProcessInfo& rCurrentProcessInfo )
{
}

void TutorialDamageModel::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                                    const ProcessInfo& rCurrentProcessInfo )
{
    if( rThisVariable == STRESSES )
        mCurrentStress = rValue;
}

void TutorialDamageModel::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                                    const ProcessInfo& rCurrentProcessInfo )
{
}


void TutorialDamageModel::InitializeMaterial( const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues )
{
    mCurrentStress = ZeroVector( 6 );
    mE = props[YOUNG_MODULUS];
    mNU = props[POISSON_RATIO];
    mE_0 = props[DAMAGE_E0];
    mE_f = props[DAMAGE_EF];
    mKappa = mE_0;
    mKappa_old = mKappa;
}

void TutorialDamageModel::ResetMaterial( const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues )
{
    mKappa = mE_0;
    mKappa_old = mKappa;
    mCurrentStress = ZeroVector( 6 );
}

void TutorialDamageModel::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

void TutorialDamageModel::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mKappa_old = mKappa;
}


void  TutorialDamageModel::CalculateMaterialResponse( const Vector& StrainVector,
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
    //equivalent strain
    double eps = sqrt(inner_prod( StrainVector, StrainVector ));
    //loading function
    double f = eps - mKappa;

    if( f > 0.0 )
    {
        mKappa = eps;
    }

    //damage parameter
    double d = 1.0 - mE_0/mKappa * exp(-(mKappa-mE_0)/mE_f);

    Matrix C = ZeroMatrix(6,6);
    double c1 = (1.0-d)*mE / (( 1.00 + mNU ) * ( 1 - 2 * mNU ) );
    double c2 = c1 * ( 1 - mNU );
    double c3 = c1 * mNU;
    double c4 = c1 * 0.5 * ( 1 - 2 * mNU );
    //filling material matrix
    C( 0, 0 ) = c2;
    C( 0, 1 ) = c3;
    C( 0, 2 ) = c3;
    C( 0, 3 ) = 0.0;
    C( 0, 4 ) = 0.0;
    C( 0, 5 ) = 0.0;
    C( 1, 0 ) = c3;
    C( 1, 1 ) = c2;
    C( 1, 2 ) = c3;
    C( 1, 3 ) = 0.0;
    C( 1, 4 ) = 0.0;
    C( 1, 5 ) = 0.0;
    C( 2, 0 ) = c3;
    C( 2, 1 ) = c3;
    C( 2, 2 ) = c2;
    C( 2, 3 ) = 0.0;
    C( 2, 4 ) = 0.0;
    C( 2, 5 ) = 0.0;
    C( 3, 0 ) = 0.0;
    C( 3, 1 ) = 0.0;
    C( 3, 2 ) = 0.0;
    C( 3, 3 ) = c4;
    C( 3, 4 ) = 0.0;
    C( 3, 5 ) = 0.0;
    C( 4, 0 ) = 0.0;
    C( 4, 1 ) = 0.0;
    C( 4, 2 ) = 0.0;
    C( 4, 3 ) = 0.0;
    C( 4, 4 ) = c4;
    C( 4, 5 ) = 0.0;
    C( 5, 0 ) = 0.0;
    C( 5, 1 ) = 0.0;
    C( 5, 2 ) = 0.0;
    C( 5, 3 ) = 0.0;
    C( 5, 4 ) = 0.0;
    C( 5, 5 ) = c4;

    noalias(AlgorithmicTangent) = C;
    StressVector = prod( C, StrainVector );
    mCurrentStress = StressVector;

}

//**********************************************************************
int TutorialDamageModel::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY

    return 0;

    KRATOS_CATCH( "" );
}
} // Namespace Kratos
