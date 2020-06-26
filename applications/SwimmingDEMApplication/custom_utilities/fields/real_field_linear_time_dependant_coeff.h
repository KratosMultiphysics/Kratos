#ifndef KRATOS_REAL_FIELD_LINEAR_TIME_DEPENDANT_COEFF_H
#define KRATOS_REAL_FIELD_LINEAR_TIME_DEPENDANT_COEFF_H
// /* External includes */

// System includes

// Project includes
#include "includes/variables.h"

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

/* Project includes */
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "real_field.h"

namespace Kratos
{
class LinearRealField: public RealField
{
public:

KRATOS_CLASS_POINTER_DEFINITION(LinearRealField);

/// Default constructor.

LinearRealField(const double& a0, const double& b0, const double& c0,
                RealFunction& fa, RealFunction& fb, RealFunction& fc)
              : mX0(a0), mY0(b0), mZ0(c0), mFx(fa), mFy(fb), mFz(fc)
{}

/// Destructor.

virtual ~LinearRealField(){}


//***************************************************************************************************************
//***************************************************************************************************************

double Evaluate(const double time, const array_1d<double, 3>& coor) override
{
    return mX0 + mFx.Evaluate(time) * coor[0] + mY0 + mFy.Evaluate(time) * coor[1] + mZ0 + mFz.Evaluate(time) * coor[2];
}

//***************************************************************************************************************
//***************************************************************************************************************

double CalculateTimeDerivative(const double time, const array_1d<double, 3>& coor) override
{
    return mFx.CalculateDerivative(time) * coor[0] + mY0 + mFy.CalculateDerivative(time) * coor[1] + mZ0 + mFz.CalculateDerivative(time) * coor[2];;
}

//***************************************************************************************************************
//***************************************************************************************************************

void CalculateGradient(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& gradient) override
{
    gradient = ZeroVector(3);
}

//***************************************************************************************************************
//***************************************************************************************************************

void CalculateLaplacian(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& laplacian) override
{
    laplacian = ZeroVector(3);
}

//***************************************************************************************************************
//***************************************************************************************************************

///@}
///@name Inquiry
///@{


///@}
///@name Input and output
///@{

/// Turn back information as a stemplate<class T, std::size_t dim> tring.

virtual std::string Info() const override
{
return "";
}

/// Print information about this object.

virtual void PrintInfo(std::ostream& rOStream) const override
{
}

/// Print object's data.

virtual void PrintData(std::ostream& rOStream) const override
{
}


///@}
///@name Friends
///@{

///@}

protected:
///@name Protected static Member r_variables
///@{


///@}
///@name Protected member r_variables
///@{ template<class T, std::size_t dim>


///@}
///@name Protected Operators
///@{


///@}
///@name Protected Operations
///@{


///@}
///@name Protected  Access
///@{

///@}
///@name Protected Inquiry
///@{


///@}
///@name Protected LifeCycle
///@{


///@}

private:

///@name Static Member r_variables
///@{


///@}
///@name Member r_variables
///@{
double mX0;
double mY0;
double mZ0;
RealFunction& mFx;
RealFunction& mFy;
RealFunction& mFz;

///@}
///@name Private Operators
///@{

///@}
///@name Private Operations
///@{


///@}
///@name Private  Access
///@{


///@}
///@name Private Inquiry
///@{


///@}
///@name Un accessible methods
///@{

/// Assignment operator.
LinearRealField & operator=(LinearRealField const& rOther);


///@}

}; // Class LinearRealField

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_REAL_FIELD_LINEAR_TIME_DEPENDANT_COEFF_H
