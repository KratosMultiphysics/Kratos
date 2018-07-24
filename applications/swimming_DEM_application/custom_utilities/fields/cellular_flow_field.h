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

/* Project includes */
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "real_functions.h"
#include "velocity_field.h"


namespace Kratos
{
class KRATOS_API(SWIMMING_DEM_APPLICATION) CellularFlowField : public VelocityField
{
public:

KRATOS_CLASS_POINTER_DEFINITION(CellularFlowField);

/// Default constructor.

CellularFlowField():VelocityField(),mL(1.0), mU(0.0), mK(2.72), mOmega(Globals::Pi)
{
    mOneOverL = 1.0 / mL;
    mOmegaUOverL = mOmega * mU / mL;
    unsigned int number_of_threads = OpenMPUtils::GetNumThreads();
    ResizeVectorsForParallelism(number_of_threads);
}

CellularFlowField(const double half_wavelength, const double max_flow_speed, const double oscillation_relative_amplitude, const double oscillation_angular_frequency)
                 :VelocityField(), mL(half_wavelength), mU(max_flow_speed), mK(oscillation_relative_amplitude), mOmega(oscillation_angular_frequency)
{

    mOneOverL = 1.0 / mL;
    mOmegaUOverL = mOmega * mU / mL;
    unsigned int maximum_number_of_threads = OpenMPUtils::GetNumThreads();
    ResizeVectorsForParallelism(maximum_number_of_threads);
}


/// Destructor.

virtual ~CellularFlowField(){}


//***************************************************************************************************************
//***************************************************************************************************************
void ResizeVectorsForParallelism(const int n_threads) override;

void UpdateCoordinates(const double time, const array_1d<double, 3>& coor, const int i_thread = 0) override;

void UpdateCoordinates(const double time, const DenseVector<double>& coor, const int i_thread = 0) override;

void LockCoordinates(const int i_thread = 0) override;

void UnlockCoordinates(const int i_thread = 0) override;
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


// Values

double U0(const int i_thread = 0) override;
double U1(const int i_thread = 0) override;
double U2(const int i_thread = 0) override;

// First-order derivatives

double U0DT(const int i_thread = 0) override;
double U0D0(const int i_thread = 0) override;
double U0D1(const int i_thread = 0) override;
double U0D2(const int i_thread = 0) override;

double U1DT(const int i_thread = 0) override;
double U1D0(const int i_thread = 0) override;
double U1D1(const int i_thread = 0) override;
double U1D2(const int i_thread = 0) override;

double U2DT(const int i_thread = 0) override;
double U2D0(const int i_thread = 0) override;
double U2D1(const int i_thread = 0) override;
double U2D2(const int i_thread = 0) override;

// Second-order derivatives

double U0DTDT(const int i_thread = 0) override;
double U0DTD0(const int i_thread = 0) override;
double U0DTD1(const int i_thread = 0) override;
double U0DTD2(const int i_thread = 0) override;
double U0D0D0(const int i_thread = 0) override;
double U0D0D1(const int i_thread = 0) override;
double U0D0D2(const int i_thread = 0) override;
double U0D1D1(const int i_thread = 0) override;
double U0D1D2(const int i_thread = 0) override;
double U0D2D2(const int i_thread = 0) override;

double U1DTDT(const int i_thread = 0) override;
double U1DTD0(const int i_thread = 0) override;
double U1DTD1(const int i_thread = 0) override;
double U1DTD2(const int i_thread = 0) override;
double U1D0D0(const int i_thread = 0) override;
double U1D0D1(const int i_thread = 0) override;
double U1D0D2(const int i_thread = 0) override;
double U1D1D1(const int i_thread = 0) override;
double U1D1D2(const int i_thread = 0) override;
double U1D2D2(const int i_thread = 0) override;

double U2DTDT(const int i_thread = 0) override;
double U2DTD0(const int i_thread = 0) override;
double U2DTD1(const int i_thread = 0) override;
double U2DTD2(const int i_thread = 0) override;
double U2D0D0(const int i_thread = 0) override;
double U2D0D1(const int i_thread = 0) override;
double U2D0D2(const int i_thread = 0) override;
double U2D1D1(const int i_thread = 0) override;
double U2D1D2(const int i_thread = 0) override;
double U2D2D2(const int i_thread = 0) override;

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
double mOneOverL;
double mOmegaUOverL;
double mOmega;
std::vector<int> mCoordinatesAreUpToDate;
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
