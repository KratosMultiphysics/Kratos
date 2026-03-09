#include "fluid_field_utility.h"
#include "swimming_DEM_application.h"

namespace Kratos
{

void FluidFieldUtility::ImposeFieldOnNodes(ModelPart& r_model_part, const VariablesList& variables_to_be_imposed)
{
    const unsigned int nnodes = r_model_part.Nodes().size();

    MarkNodesInside(r_model_part, r_model_part.GetProcessInfo());

    #pragma omp parallel for
    for (int i = 0; i < (int)nnodes; ++i){
        ModelPart::NodeIterator node_it = r_model_part.NodesBegin() + i;
        double& fluid_viscosity = node_it->FastGetSolutionStepValue(FLUID_VISCOSITY_PROJECTED);
        double& fluid_density = node_it->FastGetSolutionStepValue(FLUID_DENSITY_PROJECTED);
        fluid_viscosity = mFluidViscosity;
        fluid_density = mFluidDensity;
    }
    mpVectorField->ImposeFieldOnNodes(r_model_part, variables_to_be_imposed);
}

void FluidFieldUtility::ImposeFieldOnNodes(ModelPart& r_model_part, const Variable<array_1d<double, 3> >& variable_to_be_imposed)
{
    const double time = r_model_part.GetProcessInfo()[TIME];

    MarkNodesInside(r_model_part, r_model_part.GetProcessInfo());

    #pragma omp parallel for firstprivate(time)
    for (int i = 0; i < (int)r_model_part.Nodes().size(); ++i){
        ModelPart::NodeIterator node_it = r_model_part.NodesBegin() + i;
        double& fluid_viscosity = node_it->FastGetSolutionStepValue(FLUID_VISCOSITY_PROJECTED);
        double& fluid_density = node_it->FastGetSolutionStepValue(FLUID_DENSITY_PROJECTED);
        fluid_viscosity = mFluidViscosity;
        fluid_density = mFluidDensity;
        const array_1d<double, 3>& coordinates = node_it->Coordinates();
        array_1d<double, 3>& vector = node_it->FastGetSolutionStepValue(variable_to_be_imposed);
        mpVectorField->Evaluate(time, coordinates, vector);
    }
}



} // namespace Kratos.




