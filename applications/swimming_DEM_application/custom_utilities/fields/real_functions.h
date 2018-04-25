#ifndef KRATOS_REAL_FUNCTIONS_H
#define KRATOS_REAL_FUNCTIONS_H
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

namespace Kratos
{
class RealFunction
{
public:

KRATOS_CLASS_POINTER_DEFINITION(RealFunction);

/// Default constructor.

RealFunction(const double param1, const double param2):mA(param1), mB(param2){}

/// Destructor.

virtual ~RealFunction(){}


//***************************************************************************************************************
//***************************************************************************************************************

virtual double Evaluate(const double x)
{
    return 0.0;
}

//***************************************************************************************************************
//***************************************************************************************************************

virtual double CalculateDerivative(const double x)
{
    return 0.0;
}

//***************************************************************************************************************
//***************************************************************************************************************

virtual double CalculateSecondDerivative(const double x)
{
    return 0.0;
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
double mA;
double mB;

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
RealFunction & operator=(RealFunction const& rOther);

///@}

}; // Class RealFunction


//***************************************************************************************************************
//***************************************************************************************************************

class LinearFunction: public RealFunction
{
public:

LinearFunction(const double a, const double b): RealFunction(a, b){}

~LinearFunction(){}

double Evaluate(const double x) override
{
    return mA * x + mB;
}

double CalculateDerivative(const double x) override
{
    return mA;
}

double CalculateSecondDerivative(const double x) override
{
    return 0.0;
}

};

//***************************************************************************************************************
//***************************************************************************************************************

class PowerFunction: public RealFunction
{
public:

PowerFunction(const double a, const double b, const double c):RealFunction(a, b), mC(c){}

~PowerFunction(){}

double Evaluate(const double x) override
{
    return mA * pow(x, mB) + mC;
}

double CalculateDerivative(const double x) override
{
    return mA * mB * pow(x, mB - 1.0);
}

double CalculateSecondDerivative(const double x) override
{
    return mA * mB * (mB - 1.0) * pow(x, mB - 2.0);
}

private:

double mC;

};
//***************************************************************************************************************
//***************************************************************************************************************

class AdditionFunction: public RealFunction
{
public:

AdditionFunction(const double a, RealFunction& f, RealFunction& g):RealFunction(a, 1.0), mF(f), mG(g){}

~AdditionFunction(){}

double Evaluate(const double x) override
{
    return mA * (mF.Evaluate(x) + mG.Evaluate(x));
}

double CalculateDerivative(const double x) override
{
    return mA * (mF.CalculateDerivative(x) + mG.CalculateDerivative(x));
}

double CalculateSecondDerivative(const double x) override
{
    return mA * (mF.CalculateSecondDerivative(x) + mG.CalculateSecondDerivative(x));
}

private:

RealFunction& mF;
RealFunction& mG;
};

//***************************************************************************************************************
//***************************************************************************************************************

class ProductFunction: public RealFunction
{
public:

ProductFunction(const double a, RealFunction& f, RealFunction& g):RealFunction(a, 1.0), mF(f), mG(g){}

~ProductFunction(){}

double Evaluate(const double x) override
{
    return mA * mF.Evaluate(x) * mG.Evaluate(x);
}

double CalculateDerivative(const double x) override
{
    return mA * (mG.Evaluate(x) * mF.CalculateDerivative(x) + mF.Evaluate(x) * mG.CalculateDerivative(x));
}

double CalculateSecondDerivative(const double x) override
{
    return mA * (mF.CalculateSecondDerivative(x) * mG.Evaluate(x) + mF.Evaluate(x) * mG.CalculateSecondDerivative(x) + 2 * mF.CalculateDerivative(x) * mG.CalculateDerivative(x));
}

private:

RealFunction& mF;
RealFunction& mG;
};

//***************************************************************************************************************
//***************************************************************************************************************

class CompositionFunction: public RealFunction
{
public:

CompositionFunction(const double a, RealFunction& f, RealFunction& g):RealFunction(a, 1.0), mF(f), mG(g){}

~CompositionFunction(){}

double Evaluate(const double x) override
{
    return mA * mF.Evaluate(mG.Evaluate(x));
}

double CalculateDerivative(const double x) override
{
    return mA * mF.CalculateDerivative(mG.Evaluate(x)) * mG.CalculateDerivative(x);
}

double CalculateSecondDerivative(const double x) override
{
    return mA * (mF.CalculateSecondDerivative(x) * mG.CalculateDerivative(x) * mG.CalculateDerivative(x) + mF.CalculateDerivative(mG.Evaluate(x)) * mG.CalculateSecondDerivative(x));
}

private:

RealFunction& mF;
RealFunction& mG;
};

//***************************************************************************************************************
//***************************************************************************************************************

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.
#endif // KRATOS_REAL_FUNCTIONS_H
