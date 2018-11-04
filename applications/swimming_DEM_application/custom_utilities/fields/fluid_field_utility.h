#ifndef KRATOS_FLUID_FIELD_UTILITY_H
#define KRATOS_FLUID_FIELD_UTILITY_H
// /* External includes */

// System includes

/* Project includes */
#include "field_utility.h"
#include "velocity_field.h"
namespace Kratos
{
class KRATOS_API(SWIMMING_DEM_APPLICATION) FluidFieldUtility : public FieldUtility
{
public:

KRATOS_CLASS_POINTER_DEFINITION(FluidFieldUtility);

/// Default constructor.

FluidFieldUtility(): FieldUtility(){}

FluidFieldUtility(SpaceTimeSet::Pointer p_sts, VelocityField::Pointer p_vector_field, const double fluid_density = 1000.0, const double fluid_kinematic_viscosity = 1e-6):
    FieldUtility(p_sts, p_vector_field), mFluidDensity(fluid_density), mFluidViscosity(fluid_kinematic_viscosity){}

/// Destructor.

virtual ~FluidFieldUtility(){}

void ImposeFieldOnNodes(ModelPart& r_model_part, const VariablesList& variables_to_be_imposed) override;

void ImposeFieldOnNodes(ModelPart& r_model_part, const Variable<array_1d<double, 3> >& variable_to_be_imposed) override;

virtual void ImposeVelocityOnNodes(ModelPart& r_model_part, const Variable<array_1d<double, 3> >& container_variable)
{
    Kratos::shared_ptr<VelocityField> p_vel_field = Kratos::static_pointer_cast<VelocityField>(mpVectorField);
    p_vel_field->ImposeVelocityOnNodes(r_model_part, container_variable);
}


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
