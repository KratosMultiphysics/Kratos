#if !defined(KRATOS_ETHIR_FLOW_FIELD_H)
#define KRATOS_ETHIR_FLOW_FIELD_H

// /* External includes */

// System includes

// Project includes
#include "includes/variables.h"
#include "swimming_DEM_application.h"

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
class ProductOfSines : public VelocityField
{
public:

KRATOS_CLASS_POINTER_DEFINITION(ProductOfSines);

/// Default constructor.

ProductOfSines():VelocityField()
{
    unsigned int number_of_threads = OpenMPUtils::GetNumThreads();
    ResizeVectorsForParallelism(number_of_threads);
    mOmega = Globals::Pi;
}

ProductOfSines(const double period)
{
    unsigned int number_of_threads = OpenMPUtils::GetNumThreads();
    ResizeVectorsForParallelism(number_of_threads);

    if (period == 0.0){
        KRATOS_THROW_ERROR(std::invalid_argument, "The period must be non-negative.", 0);
    }

    mOmega = Globals::Pi / period;
}


/// Destructor.

virtual ~ProductOfSines(){}


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

// First-order derivatives

double U0D0(const int i_thread = 0) override;
double U0D1(const int i_thread = 0) override;
double U0D2(const int i_thread = 0) override;

// Second-order derivatives

double U0D0D0(const int i_thread = 0) override;
double U0D0D1(const int i_thread = 0) override;
double U0D0D2(const int i_thread = 0) override;
double U0D1D1(const int i_thread = 0) override;
double U0D1D2(const int i_thread = 0) override;
double U0D2D2(const int i_thread = 0) override;

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
double mOmega;
std::vector<int> mCoordinatesAreUpToDate;
std::vector<double> mSin0;
std::vector<double> mCos0;
std::vector<double> mSin1;
std::vector<double> mCos1;
std::vector<double> mSin2;
std::vector<double> mCos2;

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
ProductOfSines & operator=(ProductOfSines const& rOther);


///@}

}; // Class ProductOfSines

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_ETHIR_FLOW_FIELD_H defined
