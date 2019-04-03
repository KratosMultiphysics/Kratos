#if !defined(KRATOS_REAL_FIELD_H)
#define KRATOS_REAL_FIELD_H

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
#include "real_functions.h"

namespace Kratos
{
class RealField
{
public:

KRATOS_CLASS_POINTER_DEFINITION(RealField);

/// Default constructor.

RealField(){}

/// Destructor.

virtual ~RealField(){}

//***************************************************************************************************************
//***************************************************************************************************************

virtual double Evaluate(const double time, const array_1d<double, 3>& coor)
{
    return 0.0;
}

//***************************************************************************************************************
//***************************************************************************************************************

virtual double CalculateTimeDerivative(const double time, const array_1d<double, 3>& coor)
{
    return 0.0;
}

//***************************************************************************************************************
//***************************************************************************************************************

virtual void CalculateGradient(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& gradient)
{
    gradient[0] = 0.0;
    gradient[1] = 0.0;
    gradient[2] = 0.0;
}

//***************************************************************************************************************
//***************************************************************************************************************

virtual void CalculateLaplacian(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& laplacian)
{
    laplacian[0] = 0.0;
    laplacian[1] = 0.0;
    laplacian[2] = 0.0;
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

virtual std::string Info() const
{
return "";
}

/// Print information about this object.

virtual void PrintInfo(std::ostream& rOStream) const
{
}

/// Print object's data.

virtual void PrintData(std::ostream& rOStream) const
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
RealField & operator=(RealField const& rOther);


///@}

}; // Class RealField

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.
#endif // KRATOS_REAL_FIELD_H
