#include "vector_field_utility.h"
#include "swimming_DEM_application.h"

namespace Kratos
{
void VectorFieldUtility::ImposeFieldOnNodes(ModelPart& r_model_part, const VariablesList& variables_to_be_imposed)
{
    MarkNodesInside(r_model_part, r_model_part.GetProcessInfo());
    mrVectorField.ImposeFieldOnNodes(r_model_part, variables_to_be_imposed);
}


void VectorFieldUtility::ImposeFieldOnNodes(ModelPart& r_model_part, const Variable<array_1d<double, 3> >& variable_to_be_imposed)
{
    MarkNodesInside(r_model_part, r_model_part.GetProcessInfo());
    const double time = r_model_part.GetProcessInfo()[TIME];

    #pragma omp parallel for firstprivate(time)
    for (int i = 0; i < (int)r_model_part.Nodes().size(); ++i){
        ModelPart::NodeIterator node_it = r_model_part.NodesBegin() + i;
        const array_1d<double, 3>& coordinates = node_it->Coordinates();
        array_1d<double, 3>& vector = node_it->FastGetSolutionStepValue(variable_to_be_imposed);
        mrVectorField.Evaluate(time, coordinates, vector);
    }
}


} // namespace Kratos.




