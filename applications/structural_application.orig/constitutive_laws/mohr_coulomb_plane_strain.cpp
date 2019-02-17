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
*   Date:                $Date: 2009-01-16 10:50:24 $
*   Revision:            $Revision: 1.14 $
*
* ***********************************************************/


// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes

#include "includes/define.h"
#include "constitutive_laws/mohr_coulomb_plane_strain.h"

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
MohrCoulombPlaneStrain::MohrCoulombPlaneStrain()
    : ConstitutiveLaw()
{
}

/**
 * TO BE TESTED!!!
 */
MohrCoulombPlaneStrain::~MohrCoulombPlaneStrain()
{
}


bool MohrCoulombPlaneStrain::Has( const Variable<double>& rThisVariable )
{
    return false;
}

bool MohrCoulombPlaneStrain::Has( const Variable<Vector>& rThisVariable )
{
    return false;
}

bool MohrCoulombPlaneStrain::Has( const Variable<Matrix>& rThisVariable )
{
    return false;
}

double& MohrCoulombPlaneStrain::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    return( rValue );
}

Vector& MohrCoulombPlaneStrain::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    return( rValue );
}

Matrix& MohrCoulombPlaneStrain::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    return( rValue );
}

void MohrCoulombPlaneStrain::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                                       const ProcessInfo& rCurrentProcessInfo )
{
}

void MohrCoulombPlaneStrain::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                                       const ProcessInfo& rCurrentProcessInfo )
{
}

void MohrCoulombPlaneStrain::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                                       const ProcessInfo& rCurrentProcessInfo )
{
}

void MohrCoulombPlaneStrain::Calculate( const Variable<Matrix >& rVariable, Matrix& rResult,
                                        const ProcessInfo& rCurrentProcessInfo )
{
}

/**
 * TO BE TESTED!!!
 */
void MohrCoulombPlaneStrain::InitializeMaterial( const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues )
{
    mE  = props[YOUNG_MODULUS];
    mNU = props[POISSON_RATIO];
    mDE = props[DENSITY];

// KRATOS_WATCH("INSIDE THE CONS LAW")
//  KRATOS_WATCH(mE)
//   KRATOS_WATCH(mNU)
}

/**
 * TO BE TESTED!!!
 */
void MohrCoulombPlaneStrain::CalculateElasticMatrix( Matrix& C, const double E, const double NU )
{
    double c1 = E / ( 1.00 - NU * NU );
    double c2 = c1 * NU;
    double c3 = 0.5 * E / ( 1 + NU );

    C( 0, 0 ) = c1;
    C( 0, 1 ) = c2;
    C( 0, 2 ) = 0.0;
    C( 1, 0 ) = c2;
    C( 1, 1 ) = c1;
    C( 1, 2 ) = 0.0;
    C( 2, 0 ) = 0.0;
    C( 2, 1 ) = 0.0;
    C( 2, 2 ) = c3;
}

/**
 * TO BE TESTED!!!
 */
void MohrCoulombPlaneStrain::CalculateStress( const Vector& StrainVector, Vector& StressVector )
{
    double c1 = mE / ( 1.00 - mNU * mNU );
    double c2 = c1 * mNU;
    double c3 = 0.5 * mE / ( 1 + mNU );


    StressVector[0] = c1 * StrainVector[0] + c2 * ( StrainVector[1] ) ;
    StressVector[1] = c1 * StrainVector[1] + c2 * ( StrainVector[0] ) ;
    StressVector[2] = c3 * StrainVector[2];

}

//  void MohrCoulombPlaneStrain::CalculatePlaneElasticMatrix(Matrix& C, const double E, const double NU)
// {
//  double c1 = E*(1.00-NU) / ((1.00 +NU)*(1.00 - 2*NU));
//  double c2 = E*NU / ((1.00 +NU)*(1.00 - 2*NU));
//  double c3 = 0.5*E / (1 + NU);

//  C(0,0) = c1; C(0,1) = c2; C(0,2) = 0.0;
//  C(1,0) = c2; C(1,1) = c1; C(1,2) = 0.0;
//  C(2,0) = 0.0; C(2,1) = 0.0; C(2,2) = c3;
// KRATOS_WATCH("inside D")
// }

void MohrCoulombPlaneStrain::CalculateStressAndTangentMatrix( Vector& StressVector,
        const Vector& StrainVector,
        Matrix& algorithmicTangent )

{
    MohrCoulombPlaneStrain::CalculateConstitutiveMatrix( StrainVector, algorithmicTangent );
}

/**
 * TO BE REVIEWED!!!
 */
void MohrCoulombPlaneStrain::CalculateConstitutiveMatrix( const Vector& StrainVector, Matrix& rResult )
{
    CalculateElasticMatrix( rResult, mE, mNU );
}

//  void MohrCoulombPlaneStrain::PlaneStrainConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
//  {
//   CalculatePlaneElasticMatrix( rResult, mE, mNU );
//  }



//**********************************************************************
void MohrCoulombPlaneStrain::CalculateCauchyStresses(
    Vector& rCauchy_StressVector,
    const Matrix& rF,
    const Vector& rPK2_StressVector,
    const Vector& rGreenLagrangeStrainVector )
{
    Matrix S = MathUtils<double>::StressVectorToTensor( rPK2_StressVector );

    double J = MathUtils<double>::Det2( rF );

    boost::numeric::ublas::bounded_matrix<double, 2, 2> temp;
    boost::numeric::ublas::bounded_matrix<double, 2, 2> aux;

    noalias( temp ) = prod( rF, S );
    noalias( aux ) = prod( temp, trans( rF ) );
    aux *= J;

    if ( rCauchy_StressVector.size() != 3 )
        rCauchy_StressVector.resize( 3 );

    rCauchy_StressVector[0] = aux( 0, 0 );

    rCauchy_StressVector[1] = aux( 1, 1 );

    rCauchy_StressVector[2] = aux( 0,1 );
}

//**********************************************************************
//**********************************************************************


int MohrCoulombPlaneStrain::Check(const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo)
{
    if(YOUNG_MODULUS.Key() == 0 || props[YOUNG_MODULUS]<= 0.00)
        KRATOS_THROW_ERROR(std::invalid_argument,"YOUNG_MODULUS has Key zero or invalid value ","");

    const double& nu = props[POISSON_RATIO];
    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );
    if(POISSON_RATIO.Key() == 0 || check==true) // props[POISSON_RATIO] == 1.00 || props[POISSON_RATIO] == -1.00)
        KRATOS_THROW_ERROR(std::invalid_argument,"POISSON_RATIO has Key zero invalid value ","");

    if(DENSITY.Key() == 0 || props[DENSITY]<0.00)
        KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY has Key zero or invalid value ","");

    return 0;
}

void MohrCoulombPlaneStrain::Calculate( const Variable<double>& rVariable,
                                        double& Output,
                                        const ProcessInfo& rCurrentProcessInfo )
{
    Output = sqrt( mE / mDE );
}


} // Namespace Kratos
