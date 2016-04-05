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
 *   Last Modified by:    $Author: hbui $
 *   Date:                $Date: 2013-05-22 17:20:00 $
 *   Revision:            $Revision: 1.12 $
*
* ***********************************************************/


// System includes
#include <iostream>

// External includes
#include <cmath>

// Project includes

#include "includes/define.h"
#include "plane_stress.h"

#include "includes/constitutive_law.h"

#include "utilities/math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "structural_application.h"
#include "includes/properties.h"

namespace Kratos
{



/**
 *	TO BE TESTED!!!
 */
PlaneStress::PlaneStress()
    : ConstitutiveLaw()
{
}
/**
 *	TO BE TESTED!!!
 */
PlaneStress::~PlaneStress()
{
}

bool PlaneStress::Has( const Variable<int>& rThisVariable )
{
    return false;
}

bool PlaneStress::Has( const Variable<double>& rThisVariable )
{
    return false;
}

bool PlaneStress::Has( const Variable<Vector>& rThisVariable )
{
    if ( rThisVariable == STRESSES )
        return true;

    return false;
}

bool PlaneStress::Has( const Variable<Matrix>& rThisVariable )
{
    return false;
}

void PlaneStress::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}

void PlaneStress::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    //    if ( rThisVariable == PRESTRESS_FACTOR )
//        mPrestressFactor = rValue;
    if ( rThisVariable == YOUNG_MODULUS )
        mE = rValue;
    if ( rThisVariable == POISSON_RATIO )
        mNU = rValue;
}

void PlaneStress::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}

void PlaneStress::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}

void PlaneStress::Calculate(const Variable<Matrix >& rVariable, Matrix& rResult,
                            const ProcessInfo& rCurrentProcessInfo)
{
}

double& PlaneStress::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if(rThisVariable == YOUNG_MODULUS )
    {
       rValue = mE;
       return rValue;
     }
 
    if ( rThisVariable == POISSON_RATIO )
    {
        rValue = mNU;
        return rValue;
    }
    
    return rValue;

}

Vector& PlaneStress::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    if ( rThisVariable == STRESSES )
    {
//        const unsigned int size = mCurrentStress.size();
//        rValue.resize(size, false );
//        noalias( rValue ) = mCurrentStress;

        rValue.resize(6, false);
        rValue(0) = mCurrentStress(0);
        rValue(1) = mCurrentStress(1);
        rValue(2) = 0.0;
        rValue(3) = mCurrentStress(2);
        rValue(4) = 0.0;
        rValue(5) = 0.0;
    }
    return rValue;
}

/**
 *	TO BE TESTED!!!
 */
void PlaneStress::InitializeMaterial( const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues )
{
    mCurrentStress = ZeroVector( 3 );
    mE  = props[YOUNG_MODULUS];
    mNU = props[POISSON_RATIO];
}


void PlaneStress::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

void PlaneStress::InitializeNonLinearIteration( const Properties& rMaterialProperties,
                                                const GeometryType& rElementGeometry,
                                                const Vector& rShapeFunctionsValues,
                                                const ProcessInfo& rCurrentProcessInfo )
{
}

void PlaneStress::FinalizeNonLinearIteration( const Properties& rMaterialProperties,
                                              const GeometryType& rElementGeometry,
                                              const Vector& rShapeFunctionsValues,
                                              const ProcessInfo& rCurrentProcessInfo )
{
}

void PlaneStress::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

void PlaneStress::CalculateMaterialResponse( const Vector& StrainVector,
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
    if(CalculateStresses == true)
    {
        if(StressVector.size() != 3)
            StressVector.resize(3, false);
        CalculateStress(StrainVector, StressVector);
    }
    if(CalculateTangent == 1)
    {
        if(AlgorithmicTangent.size1() != 3 || AlgorithmicTangent.size2() != 3)
            AlgorithmicTangent.resize(3, 3, false);
        CalculateConstitutiveMatrix(StrainVector, AlgorithmicTangent);
    }
}

/**
 *	TO BE TESTED!!!
 */
void PlaneStress::CalculateElasticMatrix(Matrix& C, const double E, const double NU)
{
    double c1 = E / (1.00 - NU*NU);
    double c2 = c1 * NU;
    double c3 = 0.5*E / (1 + NU);

    C(0,0) = c1;
    C(0,1) = c2;
    C(0,2) = 0.0;
    C(1,0) = c2;
    C(1,1) = c1;
    C(1,2) = 0.0;
    C(2,0) = 0.0;
    C(2,1) = 0.0;
    C(2,2) = c3;
}

/**
 *	TO BE TESTED!!!
 */
void PlaneStress::CalculateStress(const Vector& StrainVector, Vector& StressVector)
{
    double c1 = mE / (1.00 - mNU*mNU);
    double c2 = c1 * mNU;
    double c3 = 0.5*mE / (1 + mNU);

    StressVector[0] = c1 * StrainVector[0] + c2 * (StrainVector[1])	;
    StressVector[1] = c1 * StrainVector[1] + c2 * (StrainVector[0])	;
    StressVector[2] = c3 * StrainVector[2];

    noalias( mCurrentStress ) = StressVector;
}

/**
 *	TO BE REVIEWED!!!
 */
void PlaneStress::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
{
    CalculateElasticMatrix( rResult, mE, mNU );
}

//**********************************************************************
void PlaneStress::CalculateCauchyStresses(
    Vector& rCauchy_StressVector,
    const Matrix& rF,
    const Vector& rPK2_StressVector,
    const Vector& rGreenLagrangeStrainVector)
{
    Matrix S = MathUtils<double>::StressVectorToTensor( rPK2_StressVector );

    double J = MathUtils<double>::Det2( rF );
    boost::numeric::ublas::bounded_matrix<double,2,2> temp;
    boost::numeric::ublas::bounded_matrix<double,2,2> aux;
    noalias(temp) = prod(rF,S);
    noalias(aux) = prod(temp,trans(rF));
    aux *= J;

    if(rCauchy_StressVector.size() != 3)
        rCauchy_StressVector.resize(3,false);

    rCauchy_StressVector[0] = aux(0,0);
    rCauchy_StressVector[1] = aux(1,1);
    rCauchy_StressVector[2] = aux(0,1);
}

//**********************************************************************
//**********************************************************************
void PlaneStress::Calculate(const Variable<double>& rVariable,
                            double& Output,
                            const ProcessInfo& rCurrentProcessInfo)
{
}


/// Turn back information as a string.
std::string PlaneStress::Info() const
{
    std::stringstream buffer;
    buffer << "Plane stress for strong discontinuities" << std::endl;
    return buffer.str();
}

int PlaneStress::Check(const Properties& props,
                       const GeometryType& geom,
                       const ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    if(YOUNG_MODULUS.Key() == 0 || props[YOUNG_MODULUS]<= 0.00)
        KRATOS_THROW_ERROR(std::invalid_argument,"YOUNG_MODULUS has Key zero or invalid value ","");

    const double& nu = props[POISSON_RATIO];
    const bool check = bool( (nu >0.999 && nu<1.01 ) || (nu < -0.999 && nu > -1.01 ) );
    if(POISSON_RATIO.Key() == 0 || check==true) // props[POISSON_RATIO] == 1.00 || props[POISSON_RATIO] == -1.00)
        KRATOS_THROW_ERROR(std::invalid_argument,"POISSON_RATIO has Key zero invalid value ","");

    return 0;

    KRATOS_CATCH("");
}



} // Namespace Kratos
