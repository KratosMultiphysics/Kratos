#ifndef KRATOS_TIME_DEPENDANT_1D_POROSITY_FIELD_H
#define KRATOS_TIME_DEPENDANT_1D_POROSITY_FIELD_H

//#include "real_functions.h"
#include "real_field.h"

namespace Kratos
{
class TimeDependant1DPorosityField: public RealField
{
//friend class TimeDependant1DPorosityField;
public:

KRATOS_CLASS_POINTER_DEFINITION(TimeDependant1DPorosityField);

/// Default constructor.
using RealField::Evaluate;

TimeDependant1DPorosityField(const double& max_time): mC(max_time)
{
    mC = 2 * max_time;
}

/// Destructor.

virtual ~TimeDependant1DPorosityField(){}

//***************************************************************************************************************
//***************************************************************************************************************

double Evaluate(const double time, const array_1d<double, 3>& coor) override
{
    return ((coor[1] - 2) / (2 * time - mC));
}

//***************************************************************************************************************
//***************************************************************************************************************

double CalculateTimeDerivative(const double time, const array_1d<double, 3>& coor) override
{
    return (2 * (2 - coor[1]) / ((2 * time - 4) * (2 * time - mC)));
}

//***************************************************************************************************************
//***************************************************************************************************************

void CalculateGradient(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& gradient) override
{
    gradient[0] = 0.0;
    gradient[1] = 1.0 / (2 * time - mC);
    gradient[2] = 0.0;
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
double mC;
RealFunction mF;
RealFunction mG;

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
TimeDependant1DPorosityField & operator=(TimeDependant1DPorosityField const& rOther);


///@}

}; // Class TimeDependant1DPorosityField

class TimeDependantForceField:RealField
{

using RealField::Evaluate;


public:

TimeDependantForceField(const double max_time): mAlpha(TimeDependant1DPorosityField(max_time)){}


virtual ~TimeDependantForceField(){}

//***************************************************************************************************************
//***************************************************************************************************************

double Evaluate(const double time, const array_1d<double, 3>& coor)
{
    array_1d<double, 3> porosity_grad;

    double porosity = mAlpha.Evaluate(time, coor);
    mAlpha.CalculateGradient(time, coor, porosity_grad);

    return (- porosity * porosity_grad[1]);
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


virtual std::string Info() const
{
    return "";
}


virtual void PrintInfo(std::ostream& rOStream) const{}

virtual void PrintData(std::ostream& rOStream) const{}

protected:

private:

RealFunction mF;
RealFunction mG;
RealField& mAlpha;

TimeDependantForceField & operator=(TimeDependantForceField const& rOther);

}; // Class TimeDependantForceField
///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_TIME_DEPENDANT_1D_POROSITY_FIELD_H
