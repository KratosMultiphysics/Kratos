#include "custom_utilities/cross_section_interpretation_utility.h"
#include "utilities/section_properties_utility.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

void CrossSectionInterpretationUtility::AssignSectionProperties(ModelPart& rModelPart)
{
    for (auto& r_prop : rModelPart.rProperties())
    {
        if (!r_prop.Has(CROSS_SECTION_TYPE)) {
            continue;
        }

        KRATOS_ERROR_IF_NOT(r_prop.Has(CROSS_SECTION_PARAMETERS))
            << "CROSS_SECTION_PARAMETERS must be defined for cross-section type: "
            << r_prop.GetValue(CROSS_SECTION_TYPE) << std::endl;

        const std::string& type = r_prop.GetValue(CROSS_SECTION_TYPE);
        const Vector& params = r_prop.GetValue(CROSS_SECTION_PARAMETERS);
        std::vector<double> geom_params(params.begin(), params.end());

        std::unordered_map<std::string, double> values;

        if (type == "RECT")
        {
            values = SectionPropertiesUtility::ComputeRectangularSection(geom_params);

        }
        else if (type == "I") {
            values = SectionPropertiesUtility::ComputeISection(geom_params);
        }
        else {
            KRATOS_ERROR << "Unsupported CROSS_SECTION_TYPE: " << type;
        }

        r_prop.SetValue(CROSS_AREA, values.at("AREA"));
        r_prop.SetValue(I22, values.at("I22"));
        r_prop.SetValue(I33, values.at("I33"));
        r_prop.SetValue(TORSIONAL_INERTIA, values.at("J"));
    }
}

} // namespace Kratos