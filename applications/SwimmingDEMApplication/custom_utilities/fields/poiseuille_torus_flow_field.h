#if !defined(KRATOS_POISEUILLE_TORUS_FLOW_FIELD_H)
#define KRATOS_POISEUILLE_TORUS_FLOW_FIELD_H

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
#include "velocity_field.h"


namespace Kratos
{
class KRATOS_API(SWIMMING_DEM_APPLICATION) PoiseuilleTorus : public VelocityField
{

public:

KRATOS_CLASS_POINTER_DEFINITION(PoiseuilleTorus);

/// Default constructor.

PoiseuilleTorus():VelocityField(), mMajorRadius(1.), mMinorRadius(.5), mCenterVelocity(1.)
{
    unsigned int number_of_threads = ParallelUtilities::GetNumThreads();
    ResizeVectorsForParallelism(number_of_threads);
};

PoiseuilleTorus(double major_radius, double minor_radius, double center_velocity):VelocityField(), mMajorRadius(major_radius), mMinorRadius(minor_radius), mCenterVelocity(center_velocity)
{
    unsigned int number_of_threads = ParallelUtilities::GetNumThreads();
    ResizeVectorsForParallelism(number_of_threads);
}

/// Destructor.

virtual ~PoiseuilleTorus(){}

double getDistanceToCenter(const array_1d<double, 3>& coor);
double getVelocityModule();
double getThetaAngle(const array_1d<double, 3>& coor);

void UpdateCoordinates(const double time, const array_1d<double, 3>& coor, const int i_thread = 0) override;
void Evaluate(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& vector, const int i_thread = 0) override;

double U0(const int i_thread = 0) override;
double U1(const int i_thread = 0) override;
double U2(const int i_thread = 0) override;

virtual void CalculateMaterialAcceleration(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& accel, const int i_thread = 0) override;

private:
double mMajorRadius;
double mMinorRadius;
double mCenterVelocity;

double mTheta;
double mRho;

double mX;
double mY;
double mZ;


}; // Class PoiseuilleTorus

}  // Namespace Kratos

#endif