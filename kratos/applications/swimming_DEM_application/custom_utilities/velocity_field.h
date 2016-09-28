#if !defined(KRATOS_VELOCITY_FIELD_H)
#define KRATOS_VELOCITY_FIELD_H

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

#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "real_functions.h"
#include "vector_field.h"
#include "swimming_DEM_application.h"


namespace Kratos
{
class VelocityField : public VectorField<3>
{
public:

KRATOS_CLASS_POINTER_DEFINITION(VelocityField);

/// Default constructor.

VelocityField():VectorField<3>(){}

/// Destructor.

virtual ~VelocityField(){}


//***************************************************************************************************************
//***************************************************************************************************************
void Evaluate(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& vector, const unsigned int n_threads = 1, const unsigned int i_thread = 1);

void CalculateTimeDerivative(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& deriv, const unsigned int n_threads = 1, const unsigned int i_thread = 1);

void CalculateGradient(const double time, const array_1d<double, 3>& coor, array_1d< array_1d<double, 3>, 3>& gradient, const unsigned int n_threads = 1, const unsigned int i_thread = 1);

void CalculateDivergence(const double time, const array_1d<double, 3>& coor, double& div, const unsigned int n_threads = 1, const unsigned int i_thread = 1);

void CalculateRotational(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& rot, const unsigned int n_threads = 1, const unsigned int i_thread = 1);

void CalculateLaplacian(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& lapl, const unsigned int n_threads = 1, const unsigned int i_thread = 1);

virtual void CalculateMaterialAcceleration(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& accel, const unsigned int n_threads = 1, const unsigned int i_thread = 1);

virtual void UpdateCoordinates(const double time, const array_1d<double, 3>& coor, const unsigned int n_threads = 1, const unsigned int i_thread = 1){}

void Evaluate(const double time, const vector<double>& coor, vector<double>& result, const unsigned int n_threads = 1, const unsigned int i_thread = 1);

void CalculateTimeDerivative(const double time, const vector<double>& coor, vector<double>& result, const unsigned int n_threads = 1, const unsigned int i_thread = 1);

double CalculateDivergence(const double time, const vector<double>& coor, const unsigned int n_threads = 1, const unsigned int i_thread = 1);

void CalculateRotational(const double time, const vector<double>& coor, vector<double>& result, const unsigned int n_threads = 1, const unsigned int i_thread = 1);

void CalculateLaplacian(const double time, const vector<double>& coor, vector<double>& result, const unsigned int n_threads = 1, const unsigned int i_thread = 1);

virtual void CalculateMaterialAcceleration(const double time, const vector<double>& coor, vector<double>& result, const unsigned int n_threads = 1, const unsigned int i_thread = 1);

void ImposeFieldOnNodes(ModelPart& r_model_part, const VariablesList& variables_to_be_imposed);
//virtual void UpdateCoordinates(const double time, const vector<double>& coor){}


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
// Values

virtual double U0(unsigned int i_thread = 1){return 0.0;}
virtual double U1(unsigned int i_thread = 1){return 0.0;}
virtual double U2(unsigned int i_thread = 1){return 0.0;}

// First-order derivatives

virtual double U0DT(unsigned int i_thread = 1){return 0.0;}
virtual double U0D0(unsigned int i_thread = 1){return 0.0;}
virtual double U0D1(unsigned int i_thread = 1){return 0.0;}
virtual double U0D2(unsigned int i_thread = 1){return 0.0;}

virtual double U1DT(unsigned int i_thread = 1){return 0.0;}
virtual double U1D0(unsigned int i_thread = 1){return 0.0;}
virtual double U1D1(unsigned int i_thread = 1){return 0.0;}
virtual double U1D2(unsigned int i_thread = 1){return 0.0;}

virtual double U2DT(unsigned int i_thread = 1){return 0.0;}
virtual double U2D0(unsigned int i_thread = 1){return 0.0;}
virtual double U2D1(unsigned int i_thread = 1){return 0.0;}
virtual double U2D2(unsigned int i_thread = 1){return 0.0;}

// Second-order derivatives

virtual double U0DTDT(unsigned int i_thread = 1){return 0.0;}
virtual double U0DTD0(unsigned int i_thread = 1){return 0.0;}
virtual double U0DTD1(unsigned int i_thread = 1){return 0.0;}
virtual double U0DTD2(unsigned int i_thread = 1){return 0.0;}
virtual double U0D0D0(unsigned int i_thread = 1){return 0.0;}
virtual double U0D0D1(unsigned int i_thread = 1){return 0.0;}
virtual double U0D0D2(unsigned int i_thread = 1){return 0.0;}
virtual double U0D1D1(unsigned int i_thread = 1){return 0.0;}
virtual double U0D1D2(unsigned int i_thread = 1){return 0.0;}
virtual double U0D2D2(unsigned int i_thread = 1){return 0.0;}

virtual double U1DTDT(unsigned int i_thread = 1){return 0.0;}
virtual double U1DTD0(unsigned int i_thread = 1){return 0.0;}
virtual double U1DTD1(unsigned int i_thread = 1){return 0.0;}
virtual double U1DTD2(unsigned int i_thread = 1){return 0.0;}
virtual double U1D0D0(unsigned int i_thread = 1){return 0.0;}
virtual double U1D0D1(unsigned int i_thread = 1){return 0.0;}
virtual double U1D0D2(unsigned int i_thread = 1){return 0.0;}
virtual double U1D1D1(unsigned int i_thread = 1){return 0.0;}
virtual double U1D1D2(unsigned int i_thread = 1){return 0.0;}
virtual double U1D2D2(unsigned int i_thread = 1){return 0.0;}

virtual double U2DTDT(unsigned int i_thread = 1){return 0.0;}
virtual double U2DTD0(unsigned int i_thread = 1){return 0.0;}
virtual double U2DTD1(unsigned int i_thread = 1){return 0.0;}
virtual double U2DTD2(unsigned int i_thread = 1){return 0.0;}
virtual double U2D0D0(unsigned int i_thread = 1){return 0.0;}
virtual double U2D0D1(unsigned int i_thread = 1){return 0.0;}
virtual double U2D0D2(unsigned int i_thread = 1){return 0.0;}
virtual double U2D1D1(unsigned int i_thread = 1){return 0.0;}
virtual double U2D1D2(unsigned int i_thread = 1){return 0.0;}
virtual double U2D2D2(unsigned int i_thread = 1){return 0.0;}

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

bool VariableIsInList(const VariablesList var_list, const VariableData& var);

/// Assignment operator.
VelocityField & operator=(VelocityField const& rOther);


///@}

}; // Class VelocityField

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_VELOCITY_FIELD_H defined
