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
class KRATOS_API(SWIMMING_DEM_APPLICATION) LinearRealField: public RealField
{
public:

KRATOS_CLASS_POINTER_DEFINITION(LinearRealField);

/// Default constructor.

LinearRealField(RealFunction::Pointer& rpA0,
                RealFunction::Pointer& rpA1,
                RealFunction::Pointer& rpA2,
                RealFunction::Pointer& rpB)
              : mA({rpA0, rpA1, rpA2}), mpB(rpB)
{
}

/// Destructor.

virtual ~LinearRealField(){}

double Evaluate(const double time,
                const array_1d<double, 3>& coor,
                const int i_thread) override;

double Evaluate(const double time,
                const DenseVector<double>& coor,
                const int i_thread) override;

double CalculateTimeDerivative(const double time,
                               const array_1d<double, 3>& coor,
                               const int i_thread) override;

double CalculateTimeDerivative(const double time,
                               const DenseVector<double>& coor,
                               const int i_thread) override;

void CalculateGradient(const double time,
                       const array_1d<double, 3>& coor,
                       array_1d<double, 3>& gradient,
                       const int i_thread) override;

void CalculateGradient(const double time,
                       const DenseVector<double>& coor,
                       DenseVector<double>& gradient,
                       const int i_thread) override;

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

std::vector<RealFunction::Pointer> mA;
RealFunction::Pointer mpB;

///@}
///@name Private Operators
///@{

///@}
///@name Private Operations
///@{
// template<TVector>
// double Evaluate(const double time,
//                 const TVector& coor,
//                 const int i_thread);

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
