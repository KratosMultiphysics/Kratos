#if !defined(KRATOS_SHEAR_FLOW_1D_WITH_EXPONENTIAL_VISCOSITY_FIELD_H)
#define KRATOS_SHEAR_FLOW_1D_WITH_EXPONENTIAL_VISCOSITY_FIELD_H

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
class KRATOS_API(SWIMMING_DEM_APPLICATION) ShearFlow1DWithExponentialViscosityField : public VelocityField
{
public:

KRATOS_CLASS_POINTER_DEFINITION(ShearFlow1DWithExponentialViscosityField);

/// Default constructor.

ShearFlow1DWithExponentialViscosityField():VelocityField(), mUFarField(0.0), mZmax(1.0), mMaxSuspensionRelativeViscosity(10)
{}

ShearFlow1DWithExponentialViscosityField(const double u_far_field, const double z_max, const double max_relative_viscosity):
                                    VelocityField(), mUFarField(u_far_field), mZmax(z_max), mMaxSuspensionRelativeViscosity(max_relative_viscosity)
{}


/// Destructor.
virtual ~ShearFlow1DWithExponentialViscosityField(){}


void Evaluate(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& vector, const int i_thread) override;


virtual std::string Info() const override
{
    return "";
}

/// Print information about this object.
virtual void PrintInfo(std::ostream& rOStream) const override {}

/// Print object's data.
virtual void PrintData(std::ostream& rOStream) const override {}

void SetRimZoneThickness(const double z_max);
void SetViscosity(const double viscosity);

protected:

    
private:

double mUFarField;
double mZmax;
double mMaxSuspensionRelativeViscosity;

/// Assignment operator.
ShearFlow1DWithExponentialViscosityField & operator=(ShearFlow1DWithExponentialViscosityField const& rOther);

}; // Class ShearFlow1DWithExponentialViscosityField


} // namespace Kratos.

#endif // KRATOS_SHEAR_FLOW_1D_WITH_EXPONENTIAL_VISCOSITY_FIELD_H defined
