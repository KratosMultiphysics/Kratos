#include "custom_utilities/section_properties_utility.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

void SectionPropertiesUtility::InterpretSections(ModelPart& rModelPart)
{
    for (auto& r_prop : rModelPart.rProperties())
    {
        // Skip if section already fully defined
        if (r_prop.Has(CROSS_AREA) &&
            r_prop.Has(I22) &&
            r_prop.Has(I33) &&
            r_prop.Has(TORSIONAL_INERTIA))
        {
            continue;
        }

        if (!r_prop.Has(CROSS_SECTION_TYPE))
            continue;

        const std::string& type = r_prop.GetValue(CROSS_SECTION_TYPE);

        if (type == "RECT")
        {
            KRATOS_ERROR_IF_NOT(r_prop.Has(CROSS_SECTION_PARAMETERS))
                << "RECT section requires CROSS_SECTION_PARAMETERS"
                << std::endl;

            const Vector& params =
                r_prop.GetValue(CROSS_SECTION_PARAMETERS);

            ComputeRectangularSection(r_prop, params);
        }
        else
        {
            KRATOS_ERROR << "Unsupported CROSS_SECTION_TYPE: "
                         << type << std::endl;
        }
    }
}

// ----------------------------------------------------------

void SectionPropertiesUtility::ComputeRectangularSection(
    Properties& rProperties,
    const Vector& rParameters)
{
    KRATOS_ERROR_IF(rParameters.size() != 2)
        << "RECT requires [b, h]" << std::endl;

    const double b = rParameters[0];
    const double h = rParameters[1];

    const double A = b * h;
    const double I22_value = b * std::pow(h,3) / 12.0;
    const double I33_value = h * std::pow(b,3) / 12.0;
    const double It = b * h * (b*b + h*h) / 12.0;

    rProperties.SetValue(CROSS_SECTION_PARAMETERS, rParameters);
    rProperties.SetValue(CROSS_AREA, A);
    rProperties.SetValue(I22, I22_value);
    rProperties.SetValue(I33, I33_value);
    rProperties.SetValue(TORSIONAL_INERTIA, It);
}

} // namespace Kratos
