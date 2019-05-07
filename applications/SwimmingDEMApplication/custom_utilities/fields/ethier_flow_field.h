#if !defined(KRATOS_ETHIER_FLOW_FIELD_H)
#define KRATOS_ETHIER_FLOW_FIELD_H

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
#include "real_field.h"
#include "velocity_field.h"


namespace Kratos
{
class KRATOS_API(SWIMMING_DEM_APPLICATION) EthierVelocityField : public VelocityField
{
public:

KRATOS_CLASS_POINTER_DEFINITION(EthierVelocityField);

/// Default constructor.

EthierVelocityField():VelocityField(), mA(0.25 * Globals::Pi), mD(0.5 * Globals::Pi)
{
    unsigned int number_of_threads = OpenMPUtils::GetNumThreads();
    ResizeVectorsForParallelism(number_of_threads);
}

EthierVelocityField(const double a, const double b)
                 :VelocityField(), mA(a), mD(b)
{
    unsigned int number_of_threads = OpenMPUtils::GetNumThreads();
    ResizeVectorsForParallelism(number_of_threads);
}


/// Destructor.

virtual ~EthierVelocityField(){}


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

std::string Info() const override
{
    return "";
}

/// Print information about this object.

void PrintInfo(std::ostream& rOStream) const override
{
}

/// Print object's data.

void PrintData(std::ostream& rOStream) const override
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

double mA;
double mD;
std::vector<int> mCoordinatesAreUpToDate;
std::vector<double> mExpD2T;
std::vector<double> mExpAX;
std::vector<double> mExpAZ;
std::vector<double> mExpAY;
std::vector<double> mSinAXDY;
std::vector<double> mCosAXDY;
std::vector<double> mSinAYDZ;
std::vector<double> mCosAYDZ;
std::vector<double> mSinAZDX;
std::vector<double> mCosAZDX;


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
EthierVelocityField & operator=(EthierVelocityField const& rOther);


///@}

}; // Class EthierVelocityField


class KRATOS_API(SWIMMING_DEM_APPLICATION) EthierPressureField : public RealField
{
public:

KRATOS_CLASS_POINTER_DEFINITION(EthierPressureField);

/// Default constructor.

EthierPressureField():RealField(), mA(0.25 * Globals::Pi), mD(0.5 * Globals::Pi)
{
    unsigned int number_of_threads = OpenMPUtils::GetNumThreads();
    ResizeVectorsForParallelism(number_of_threads);
}

EthierPressureField(const double a, const double b)
                 :RealField(), mA(a), mD(b)
{
    unsigned int number_of_threads = OpenMPUtils::GetNumThreads();
    ResizeVectorsForParallelism(number_of_threads);
}


/// Destructor.

virtual ~EthierPressureField(){}

void ResizeVectorsForParallelism(const int n_threads);

//***************************************************************************************************************
//***************************************************************************************************************
double Evaluate(const double time, const array_1d<double, 3>& coor, const int i_thread) override;
double CalculateTimeDerivative(const double time, const array_1d<double, 3>& coor, const int i_thread) override;
void CalculateGradient(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& gradient, const int i_thread) override;
void CalculateLaplacian(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& laplacian, const int i_thread) override;

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

double mA;
double mD;
std::vector<int> mCoordinatesAreUpToDate;
std::vector<double> mExpD2T;
std::vector<double> mExpAX;
std::vector<double> mExpAZ;
std::vector<double> mExpAY;
std::vector<double> mSinAXDY;
std::vector<double> mCosAXDY;
std::vector<double> mSinAYDZ;
std::vector<double> mCosAYDZ;
std::vector<double> mSinAZDX;
std::vector<double> mCosAZDX;


///@}
///@name Private Operators
///@{

///@}
///@name Private Operations
///@{
double U(const int i_thread = 0) override;

// First-order derivatives
double UDT(const int i_thread = 0) override;
double UD1(const int i_thread = 0) override;
double UD2(const int i_thread = 0) override;
double UD3(const int i_thread = 0) override;

// Second-order derivatives
double UDTDT(const int i_thread = 0) override;
double UDTD0(const int i_thread = 0) override;
double UDTD1(const int i_thread = 0) override;
double UDTD2(const int i_thread = 0) override;
double UD0D0(const int i_thread = 0) override;
double UD0D1(const int i_thread = 0) override;
double UD0D2(const int i_thread = 0) override;
double UD1D1(const int i_thread = 0) override;
double UD1D2(const int i_thread = 0) override;
double UD2D2(const int i_thread = 0) override;


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
EthierPressureField & operator=(EthierPressureField const& rOther);


///@}

}; // Class EthierPressureField

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_ETHIER_FLOW_FIELD_H defined
