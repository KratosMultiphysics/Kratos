#ifndef KRATOS_CODINA_MANUFACTURED_FIELD_H
#define KRATOS_CODINA_MANUFACTURED_FIELD_H

#include "real_functions.h"
#include "real_field.h"

namespace Kratos
{
class CodinaManufacturedField:RealField
{
public:

KRATOS_CLASS_POINTER_DEFINITION(CodinaManufacturedField);

/// Default constructor.

CodinaManufacturedField(const RealField& alpha):RealField(), mAlpha(Alpha){}

/// Destructor.

virtual ~CodinaManufacturedField(){}

//***************************************************************************************************************
//***************************************************************************************************************

double Evaluate(const double time, const array_1d<double, 3>& coor)
{
    return 0.0;
}

//***************************************************************************************************************
//***************************************************************************************************************

double CalculateTimeDerivative(const double time, const array_1d<double, 3>& coor)
{
    return 0.0;
}

//***************************************************************************************************************
//***************************************************************************************************************

void CalculateGradient(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& gradient){}
void CalculateLaplacian(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& laplacian){}

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
RealFunction mF;
RealFunction mG;
RealField& mAlpha;
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
CodinaManufacturedField & operator=(CodinaManufacturedField const& rOther);


///@}

}; // Class CodinaManufacturedField

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.


#endif // KRATOS_CODINA_MANUFACTURED_FIELD_H
