#include "scalar_field_utility.h"
#include "swimming_DEM_application.h"

namespace Kratos
{
void RealFieldUtility::ImposeFieldOnNodes(ModelPart& r_model_part, const VariablesList& variables_to_be_imposed)
{
    MarkNodesInside(r_model_part, r_model_part.GetProcessInfo());
    mrRealField.ImposeFieldOnNodes(r_model_part, variables_to_be_imposed);
}

} // namespace Kratos.




