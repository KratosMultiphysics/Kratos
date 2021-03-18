/*
==============================================================================
see soil_mechanics_appplication/LICENSE.txt
==============================================================================
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2 Mar 2021 $
//   Revision:            $Revision: 1.0 $
//
//



#include "univariate_phase_law.h"


namespace Kratos
{

UnivariatePhaseLaw::UnivariatePhaseLaw()
{}

UnivariatePhaseLaw::~UnivariatePhaseLaw()
{}

UnivariatePhaseLaw::Pointer UnivariatePhaseLaw::Clone() const
{
    return UnivariatePhaseLaw::Pointer(new UnivariatePhaseLaw());
}

bool UnivariatePhaseLaw::Has( const Variable<int>& rThisVariable )
{
}

bool UnivariatePhaseLaw::Has( const Variable<double>& rThisVariable )
{
}

bool UnivariatePhaseLaw::Has( const Variable<Vector>& rThisVariable )
{
}

bool UnivariatePhaseLaw::Has( const Variable<Matrix>& rThisVariable )
{
}

int& UnivariatePhaseLaw::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    return rValue;
}

double& UnivariatePhaseLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    return rValue;
}

Vector& UnivariatePhaseLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    return rValue;
}

Matrix& UnivariatePhaseLaw::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    return rValue;
}

void UnivariatePhaseLaw::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                           const ProcessInfo& rCurrentProcessInfo )
{
}

void UnivariatePhaseLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                           const ProcessInfo& rCurrentProcessInfo )
{
}

void UnivariatePhaseLaw::SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                           const array_1d<double, 3>& rValue,
                           const ProcessInfo& rCurrentProcessInfo )
{
}

void UnivariatePhaseLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                           const ProcessInfo& rCurrentProcessInfo )
{
}

void UnivariatePhaseLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                           const ProcessInfo& rCurrentProcessInfo )
{
}

double UnivariatePhaseLaw::GetValue(const double& v) const
{
    return 0.0;
}

double UnivariatePhaseLaw::GetDerivative(const double& v) const
{
    return 0.0;
}

double UnivariatePhaseLaw::GetSecondDerivative(const double& v) const
{
    return 0.0;
}

} /* namespace Kratos.*/
