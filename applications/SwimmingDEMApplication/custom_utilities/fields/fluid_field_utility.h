#ifndef KRATOS_FLUID_FIELD_UTILITY_H
#define KRATOS_FLUID_FIELD_UTILITY_H
// /* External includes */

// System includes

/* Project includes */
#include "scalar_field_utility.h"
#include "vector_field_utility.h"
#include "velocity_field.h"
namespace Kratos
{
class KRATOS_API(SWIMMING_DEM_APPLICATION) FluidFieldUtility
{
public:

KRATOS_CLASS_POINTER_DEFINITION(FluidFieldUtility);

/// Default constructor.

FluidFieldUtility(SpaceTimeSet& rDomain,
                  RealField& PressureField,
                  VelocityField& rVelocityField,
                  const double FluidDensity = 1000.0,
                  const double FluidKinematicViscosity = 1e-6):
    mrDomain(rDomain),
    mrPressureField(PressureField),
    mrVelocityField(rVelocityField),
    mFluidDensity(FluidDensity),
    mFluidViscosity(FluidKinematicViscosity){}

/// Destructor.

virtual ~FluidFieldUtility(){}

virtual void ImposeFieldOnNodes(ModelPart& r_model_part, const VariablesList& variables_to_be_imposed);

virtual void ImposeFieldOnNodes(ModelPart& r_model_part, const Variable<array_1d<double, 3> >& fluid_variable_to_be_imposed, const Variable<double >& pressure_variable_to_be_imposed);

virtual void ImposeVelocityOnNodes(ModelPart& r_model_part, const Variable<array_1d<double, 3> >& container_variable)
{
    mrVelocityField.ImposeVelocityOnNodes(r_model_part, container_variable);
}

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
///@{

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
SpaceTimeSet& mrDomain;
RealField& mrPressureField;
VelocityField& mrVelocityField;
double mFluidDensity;
double mFluidViscosity;
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
FluidFieldUtility & operator=(FluidFieldUtility const& rOther);

///@}

}; // Class FluidFieldUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.
#endif // KRATOS_FLUID_FIELD_UTILITY_H
