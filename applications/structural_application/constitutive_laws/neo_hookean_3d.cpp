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
#include "constitutive_laws/neo_hookean_3d.h"

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
NeoHookean3D::NeoHookean3D()
    : ConstitutiveLaw()
{
}

/**
 * TO BE TESTED!!!
 */
NeoHookean3D::~NeoHookean3D()
{
}

bool NeoHookean3D::Has( const Variable<int>& rThisVariable )
{
    return false;
}

bool NeoHookean3D::Has( const Variable<double>& rThisVariable )
{
    if ( rThisVariable == PRESTRESS_FACTOR || rThisVariable == YOUNG_MODULUS || rThisVariable == POISSON_RATIO )
        return true;

    return false;
}

bool NeoHookean3D::Has( const Variable<Vector>& rThisVariable )
{
    if ( rThisVariable == INSITU_STRESS )
        return true;

    if ( rThisVariable == PRESTRESS )
        return true;

    if ( rThisVariable == STRESSES )
        return true;

    return false;
}

bool NeoHookean3D::Has( const Variable<Matrix>& rThisVariable )
{
    return false;
}


int& NeoHookean3D::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    return rValue;
}

double& NeoHookean3D::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if ( rThisVariable == PRESTRESS_FACTOR ){
        rValue = mPrestressFactor;
    return rValue;
    }
    if(rThisVariable == YOUNG_MODULUS ){
       rValue = mE;
       return rValue;
     }


    if ( rThisVariable == POISSON_RATIO ){
        rValue = mNU;
        return rValue;
    }

    if(rThisVariable==DAMAGE){
        rValue = 0.00;
        return rValue;
    }

    if (rThisVariable==DELTA_TIME){
        rValue = sqrt(mE/mDE);
        return rValue;
    }

    rValue = 0.00;
    return rValue;
}

Vector& NeoHookean3D::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    if ( rThisVariable == INSITU_STRESS || rThisVariable == PRESTRESS )
    {
        const unsigned int size = mPrestress.size();
        if(rValue.size() != size)
            rValue.resize(size, false);
        noalias(rValue) = mPrestressFactor * mPrestress;
        return rValue;
    }

    if ( rThisVariable == STRESSES )
    {
        const unsigned int size = mCurrentStress.size();
        if(rValue.size() != size)
            rValue.resize(size, false);
        noalias(rValue)  = mCurrentStress;
        return rValue;
    }

    if ( rThisVariable == PLASTIC_STRAIN_VECTOR )
    {
        noalias(rValue) = ZeroVector( rValue.size() );
        return( rValue );
    }

    if ( rThisVariable == INTERNAL_VARIABLES )
    {
        rValue = ZeroVector( 1 );
        return( rValue );
    }

    KRATOS_THROW_ERROR( std::logic_error, "Vector Variable case not considered", "" );
}

Matrix& NeoHookean3D::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
        KRATOS_THROW_ERROR( std::logic_error, "Matrix Variable case not considered", "" );
}

void NeoHookean3D::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}

void NeoHookean3D::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == PRESTRESS_FACTOR )
        mPrestressFactor = rValue;
    if ( rThisVariable == YOUNG_MODULUS )
        mE = rValue;
    if ( rThisVariable == POISSON_RATIO )
        mNU = rValue;
}

void NeoHookean3D::SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                            const array_1d<double, 3>& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}

void NeoHookean3D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == INSITU_STRESS || rThisVariable == PRESTRESS )
    {
        noalias(mPrestress) = rValue;
    }
}

void NeoHookean3D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}


void NeoHookean3D::InitializeMaterial( const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues )
{
    mCurrentStress = ZeroVector( 6 );
    mPrestress = ZeroVector( 6 );
    mPrestressFactor = 1.0;
    mE = props[YOUNG_MODULUS];
    mNU = props[POISSON_RATIO];
    mDE = props[DENSITY];
}

void NeoHookean3D::ResetMaterial( const Properties& props,
                                 const GeometryType& geom,
                                 const Vector& ShapeFunctionsValues )
{
    mPrestress = ZeroVector( 6 );
    mPrestressFactor = 1.0;
}

void NeoHookean3D::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

void NeoHookean3D::InitializeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

void NeoHookean3D::FinalizeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

void NeoHookean3D::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

void  NeoHookean3D::CalculateMaterialResponse( const Vector& StrainVector,
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
    CalculateStress( StressVector, StrainVector );
    CalculateTangentMatrix( AlgorithmicTangent, StrainVector );
}

/**
 * TO BE TESTED!!!
 */
void NeoHookean3D::CalculateTangentMatrix( Matrix& C, const Vector& StrainVector )
{
    if ( C.size1() != 6 || C.size2() != 6 )
    {
        C.resize( 6, 6 );
    }

    double e_11 = StrainVector(0);
    double e_22 = StrainVector(1);
    double e_33 = StrainVector(2);
    double e_12 = StrainVector(3);
    double e_23 = StrainVector(4);
    double e_13 = StrainVector(5);

    double mu = mE / (2 * (1.0 + mNU));
    double lambda = mE * mNU / ((1.0 + mNU) * (1.0 - 2*mNU));

    double aux = 4*e_11*e_22 + 4*e_11*e_33 + 2*e_11 - pow(e_12, 2) - pow(e_13, 2) + 4*e_22*e_33 + 2*e_22 - pow(e_23, 2) + 2*e_33 + 1;

    C( 0, 0 ) = pow(2*e_22 + 2*e_33 + 1, 2)*(-lambda*log(aux) + lambda + 2*mu)/pow(aux, 2);
    C( 0, 1 ) = ((lambda*log(aux) - 2*mu)*(aux) + (2*e_11 + 2*e_33 + 1)*(2*e_22 + 2*e_33 + 1)*(-lambda*log(aux) + lambda + 2*mu))/pow(aux, 2);
    C( 0, 2 ) = ((lambda*log(aux) - 2*mu)*(aux) + (2*e_11 + 2*e_22 + 1)*(2*e_22 + 2*e_33 + 1)*(-lambda*log(aux) + lambda + 2*mu))/pow(aux, 2);
    C( 0, 3 ) = e_12*(2*e_22 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 0, 4 ) = e_23*(2*e_22 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 0, 5 ) = e_13*(2*e_22 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);

    C( 1, 0 ) = ((lambda*log(aux) - 2*mu)*(aux) + (2*e_11 + 2*e_33 + 1)*(2*e_22 + 2*e_33 + 1)*(-lambda*log(aux) + lambda + 2*mu))/pow(aux, 2);
    C( 1, 1 ) = pow(2*e_11 + 2*e_33 + 1, 2)*(-lambda*log(aux) + lambda + 2*mu)/pow(aux, 2);
    C( 1, 2 ) = ((lambda*log(aux) - 2*mu)*(aux) + (2*e_11 + 2*e_22 + 1)*(2*e_11 + 2*e_33 + 1)*(-lambda*log(aux) + lambda + 2*mu))/pow(aux, 2);
    C( 1, 3 ) = e_12*(2*e_11 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 1, 4 ) = e_23*(2*e_11 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 1, 5 ) = e_13*(2*e_11 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);

    C( 2, 0 ) = ((lambda*log(aux) - 2*mu)*(aux) + (2*e_11 + 2*e_22 + 1)*(2*e_22 + 2*e_33 + 1)*(-lambda*log(aux) + lambda + 2*mu))/pow(aux, 2);
    C( 2, 1 ) = ((lambda*log(aux) - 2*mu)*(aux) + (2*e_11 + 2*e_22 + 1)*(2*e_11 + 2*e_33 + 1)*(-lambda*log(aux) + lambda + 2*mu))/pow(aux, 2);
    C( 2, 2 ) = pow(2*e_11 + 2*e_22 + 1, 2)*(-lambda*log(aux) + lambda + 2*mu)/pow(aux, 2);
    C( 2, 3 ) = e_12*(2*e_11 + 2*e_22 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 2, 4 ) = e_23*(2*e_11 + 2*e_22 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 2, 5 ) = e_13*(2*e_11 + 2*e_22 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);

    C( 3, 0 ) = e_12*(2*e_22 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 3, 1 ) = e_12*(2*e_11 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 3, 2 ) = e_12*(2*e_11 + 2*e_22 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 3, 3 ) = (2*pow(e_12, 2)*(-lambda*log(aux) + lambda + 2*mu) - lambda*(aux)*log(aux) + 2*mu*(aux))/(2*pow(aux, 2));
    C( 3, 4 ) = e_12*e_23*(-lambda*log(aux) + lambda + 2*mu)/pow(aux, 2);
    C( 3, 5 ) = e_12*e_13*(-lambda*log(aux) + lambda + 2*mu)/pow(aux, 2);

    C( 4, 0 ) = e_23*(2*e_22 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 4, 1 ) = e_23*(2*e_11 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 4, 2 ) = e_23*(2*e_11 + 2*e_22 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 4, 3 ) = e_12*e_23*(-lambda*log(aux) + lambda + 2*mu)/pow(aux, 2);
    C( 4, 4 ) = (2*pow(e_23, 2)*(-lambda*log(aux) + lambda + 2*mu) - lambda*(aux)*log(aux) + 2*mu*(aux))/(2*pow(aux, 2));
    C( 4, 5 ) = e_13*e_23*(-lambda*log(aux) + lambda + 2*mu)/pow(aux, 2);

    C( 5, 0 ) = e_13*(2*e_22 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 5, 1 ) = e_13*(2*e_11 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 5, 2 ) = e_13*(2*e_11 + 2*e_22 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 5, 3 ) = e_12*e_13*(-lambda*log(aux) + lambda + 2*mu)/pow(aux, 2);
    C( 5, 4 ) = e_13*e_23*(-lambda*log(aux) + lambda + 2*mu)/pow(aux, 2);
    C( 5, 5 ) = (2*pow(e_13, 2)*(-lambda*log(aux) + lambda + 2*mu) - lambda*(aux)*log(aux) + 2*mu*(aux))/(2*pow(aux, 2));
}

/**
 * TO BE TESTED!!!
 */
void NeoHookean3D::CalculateStress( Vector& StressVector, const Vector& StrainVector )
{
    if ( StressVector.size() != 6 )
    {
        StressVector.resize( 6 );
    }

    double e_11 = StrainVector(0);
    double e_22 = StrainVector(1);
    double e_33 = StrainVector(2);
    double e_12 = StrainVector(3);
    double e_23 = StrainVector(4);
    double e_13 = StrainVector(5);

    double mu = mE / (2 * (1.0 + mNU));
    double lambda = mE * mNU / ((1.0 + mNU) * (1.0 - 2*mNU));

    double aux = 4*e_11*e_22 + 4*e_11*e_33 + 2*e_11 - pow(e_12, 2) - pow(e_13, 2) + 4*e_22*e_33 + 2*e_22 - pow(e_23, 2) + 2*e_33 + 1;

    StressVector(0) = (lambda*(2*e_22 + 2*e_33 + 1)*log(aux)/2 - mu*(2*e_22 + 2*e_33 + 1) + mu*aux) / aux;

    StressVector(1) = (lambda*(2*e_11 + 2*e_33 + 1)*log(aux)/2 - mu*(2*e_11 + 2*e_33 + 1) + mu*aux) / aux;

    StressVector(2) = (lambda*(2*e_11 + 2*e_22 + 1)*log(aux)/2 - mu*(2*e_11 + 2*e_22 + 1) + mu*aux) / aux;

    StressVector(3) = e_12*(-lambda*log(aux) + 2*mu) / (2*aux);

    StressVector(4) = e_23*(-lambda*log(aux) + 2*mu) / (2*aux);

    StressVector(5) = e_13*(-lambda*log(aux) + 2*mu) / (2*aux);

    noalias( StressVector ) -= mPrestressFactor * mPrestress;

    noalias(mCurrentStress) = StressVector;
}


//**********************************************************************
void NeoHookean3D::CalculateCauchyStresses(
    Vector& rCauchy_StressVector,
    const Matrix& rF,
    const Vector& rPK2_StressVector,
    const Vector& rGreenLagrangeStrainVector )
{
    Matrix S = MathUtils<double>::StressVectorToTensor( rPK2_StressVector );

    double J = MathUtils<double>::Det3( rF );
    boost::numeric::ublas::bounded_matrix<double, 3, 3> mstemp;
    boost::numeric::ublas::bounded_matrix<double, 3, 3> msaux;

    noalias( mstemp ) = prod( rF, S );
    noalias( msaux ) = prod( mstemp, trans( rF ) );
    msaux *= J;

    if ( rCauchy_StressVector.size() != 6 )
        rCauchy_StressVector.resize( 6 );

    rCauchy_StressVector[0] = msaux( 0, 0 );

    rCauchy_StressVector[1] = msaux( 1, 1 );

    rCauchy_StressVector[2] = msaux( 2, 2 );

    rCauchy_StressVector[3] = msaux( 0, 1 );

    rCauchy_StressVector[4] = msaux( 0, 2 );

    rCauchy_StressVector[5] = msaux( 1, 2 );
}

//**********************************************************************
int NeoHookean3D::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY

    if ( !props.Has( YOUNG_MODULUS ) || !props.Has(POISSON_RATIO) )
    {
        KRATOS_THROW_ERROR( std::logic_error, "this constitutive law requires YOUNG_MODULUS and POISSON_RATIO given as KRATOS variables", "" );
    }

    double nu = props[POISSON_RATIO];

    if ( nu > 0.499 && nu < 0.501 )
    {
        KRATOS_THROW_ERROR( std::logic_error, "invalid poisson ratio in input, close to incompressibility", "" );
        return -1;
    }

    return 0;

    KRATOS_CATCH( "" );
}
} // Namespace Kratos
