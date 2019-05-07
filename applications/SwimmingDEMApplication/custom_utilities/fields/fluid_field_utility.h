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

void MarkNodesInside(ModelPart& r_model_part, const ProcessInfo& r_current_process_info)
{
    const int nnodes = r_model_part.Nodes().size();
    const double time = r_current_process_info[TIME];

    #pragma omp parallel for
    for (int i = 0; i < nnodes; ++i){
        ModelPart::NodeIterator node_it = r_model_part.NodesBegin() + i;
        double coor_x = node_it->X();
        double coor_y = node_it->Y();
        double coor_z = node_it->Z();
        bool is_in = mrDomain.IsIn(time, coor_x, coor_y, coor_z);
        node_it->Set(INSIDE, is_in);
    }
}

virtual void ImposeFieldOnNodes(ModelPart& r_model_part, const VariablesList& variables_to_be_imposed);

virtual void ImposeFieldOnNodes(ModelPart& r_model_part,
                                const Variable<array_1d<double, 3> >& fluid_variable_to_be_imposed,
                                const Variable<double >& pressure_variable_to_be_imposed);

virtual void ImposePressureFieldOnNodes(ModelPart& r_model_part, const Variable<double>& container_variable)
{

}

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
