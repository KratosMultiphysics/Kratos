#if !defined(KRATOS_CONSTANT_VELOCITY_FIELD_H)
#define KRATOS_CONSTANT_VELOCITY_FIELD_H

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
class KRATOS_API(SWIMMING_DEM_APPLICATION) ConstantVelocityField : public VelocityField
{
public:

KRATOS_CLASS_POINTER_DEFINITION(ConstantVelocityField);

/// Default constructor.

ConstantVelocityField():VelocityField(), mVx(0.0), mVy(0.0), mVz(0.0) {}

ConstantVelocityField(const double vx, const double vy, const double vz):VelocityField(), mVx(vx), mVy(vy), mVz(vz) {}


/// Destructor.
virtual ~ConstantVelocityField(){}


void Evaluate(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& vector, const int i_thread) override;


virtual std::string Info() const override
{
    return "";
}

/// Print information about this object.
virtual void PrintInfo(std::ostream& rOStream) const override {}

/// Print object's data.
virtual void PrintData(std::ostream& rOStream) const override {}



protected:

    
    
private:

double mVx;
double mVy;
double mVz;

/// Assignment operator.
ConstantVelocityField & operator=(ConstantVelocityField const& rOther);

}; // Class ConstantVelocityField


} // namespace Kratos.

#endif // KRATOS_CONSTANT_VELOCITY_FIELD_H defined
