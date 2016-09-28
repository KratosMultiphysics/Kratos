#if !defined(KRATOS_CELLULAR_FLOW_FIELD_H)
#define KRATOS_CELLULAR_FLOW_FIELD_H

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
#include "velocity_field.h"


namespace Kratos
{
class CellularFlowField : public VelocityField
{
public:

KRATOS_CLASS_POINTER_DEFINITION(CellularFlowField);

/// Default constructor.

CellularFlowField():VelocityField(),mL(1.0), mU(0.0), mK(2.72), mOmega(KRATOS_M_PI)
{
    mPiOverL = 0.5 * KRATOS_M_PI / mL;
    mOmegaUOverL = 0.5 * mOmega * mU / mL;

    unsigned int number_of_threads = omp_get_max_threads();
    ResizeVectorsForParallelism(number_of_threads);
}

CellularFlowField(const double wavelength, const double mean_flow_speed, const double oscillation_relative_amplitude, const double oscillation_angular_frequency)
                 :VelocityField(), mL(0.5 * wavelength), mU(mean_flow_speed), mK(oscillation_relative_amplitude), mOmega(oscillation_angular_frequency)
{
    mPiOverL = 0.5 * KRATOS_M_PI / mL;
    mOmegaUOverL = 0.5 * mOmega * mU / mL;

    unsigned int number_of_threads = omp_get_max_threads();
    ResizeVectorsForParallelism(number_of_threads);
}


/// Destructor.

virtual ~CellularFlowField(){}


//***************************************************************************************************************
//***************************************************************************************************************

void ResizeVectorsForParallelism(const unsigned int n_threads);

void UpdateCoordinates(const double time, const array_1d<double, 3>& coor, const unsigned int i_thread = 1);

void UpdateCoordinates(const double time, const vector<double>& coor, const unsigned int i_thread = 1);

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


// Values

double U0(unsigned int i_thread = 1);
double U1(unsigned int i_thread = 1);
double U2(unsigned int i_thread = 1);

// First-order derivatives

double U0DT(unsigned int i_thread = 1);
double U0D0(unsigned int i_thread = 1);
double U0D1(unsigned int i_thread = 1);
double U0D2(unsigned int i_thread = 1);

double U1DT(unsigned int i_thread = 1);
double U1D0(unsigned int i_thread = 1);
double U1D1(unsigned int i_thread = 1);
double U1D2(unsigned int i_thread = 1);

double U2DT(unsigned int i_thread = 1);
double U2D0(unsigned int i_thread = 1);
double U2D1(unsigned int i_thread = 1);
double U2D2(unsigned int i_thread = 1);

// Second-order derivatives

double U0DTDT(unsigned int i_thread = 1);
double U0DTD0(unsigned int i_thread = 1);
double U0DTD1(unsigned int i_thread = 1);
double U0DTD2(unsigned int i_thread = 1);
double U0D0D0(unsigned int i_thread = 1);
double U0D0D1(unsigned int i_thread = 1);
double U0D0D2(unsigned int i_thread = 1);
double U0D1D1(unsigned int i_thread = 1);
double U0D1D2(unsigned int i_thread = 1);
double U0D2D2(unsigned int i_thread = 1);

double U1DTDT(unsigned int i_thread = 1);
double U1DTD0(unsigned int i_thread = 1);
double U1DTD1(unsigned int i_thread = 1);
double U1DTD2(unsigned int i_thread = 1);
double U1D0D0(unsigned int i_thread = 1);
double U1D0D1(unsigned int i_thread = 1);
double U1D0D2(unsigned int i_thread = 1);
double U1D1D1(unsigned int i_thread = 1);
double U1D1D2(unsigned int i_thread = 1);
double U1D2D2(unsigned int i_thread = 1);

double U2DTDT(unsigned int i_thread = 1);
double U2DTD0(unsigned int i_thread = 1);
double U2DTD1(unsigned int i_thread = 1);
double U2DTD2(unsigned int i_thread = 1);
double U2D0D0(unsigned int i_thread = 1);
double U2D0D1(unsigned int i_thread = 1);
double U2D0D2(unsigned int i_thread = 1);
double U2D1D1(unsigned int i_thread = 1);
double U2D1D2(unsigned int i_thread = 1);
double U2D2D2(unsigned int i_thread = 1);

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
///

///@}
///@name Member r_variables
///@{
double mL;
double mU;
double mK;
double mPiOverL;
double mOmegaUOverL;
double mOmega;
std::vector<double> mSinOmegaT;
std::vector<double> mCosOmegaT;
std::vector<double> mSinPiX0;
std::vector<double> mCosPiX0;
std::vector<double> mSinPiX1;
std::vector<double> mCosPiX1;

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
CellularFlowField & operator=(CellularFlowField const& rOther);


///@}

}; // Class CellularFlowField

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_CELLULAR_FLOW_FIELD_H defined
