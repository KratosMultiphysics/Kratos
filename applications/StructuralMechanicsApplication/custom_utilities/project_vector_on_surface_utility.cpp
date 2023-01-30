//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Klaus Sautter
//                   Philipp Bucher
//
//

// System includes

// External includes

// Project includes
#include "project_vector_on_surface_utility.h"
#include "structural_mechanics_application_variables.h"


namespace Kratos {

namespace {

typedef std::size_t SizeType;

void CheckVectorNorm(const Vector3& rVector)
{
    KRATOS_ERROR_IF(norm_2(rVector) < 1e-12) << "Vector has length of zero!" << std::endl;
}

void NormalizeVector(Vector3& rVector)
{
    CheckVectorNorm(rVector);
    rVector / norm_2(rVector);
}

Vector3 CheckAndReadNormalizedVector(Parameters VectorParam)
{
    const Vector vec = VectorParam.GetVector();
    KRATOS_ERROR_IF_NOT(vec.size() == 3) << "Vector is not of size 3!" << std::endl;

    Vector3 vec_return;
    vec_return[0] = vec[0];
    vec_return[1] = vec[1];
    vec_return[2] = vec[2];

    NormalizeVector(vec_return);

    return vec_return;
}

void ValidateParameters(Parameters ThisParameters)
{
    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"  : "Structure",
            "echo_level"       : 0,
            "projection_type"  : "planar",
            "global_direction" : [1,0,0],
            "variable_name"    : "PLEASE_SPECIFY",
            "visualize_in_vtk" : false,
            "method_specific_settings" : { },
            "check_local_space_dimension" : true
        })" );

    ThisParameters.ValidateAndAssignDefaults(default_parameters);
}

bool CheckElementLocalSpaceDimension(const bool& rCheckLocalSpaceDimension, const SizeType ElementID, const SizeType ElementLocalSpaceDimension)
{
    if (rCheckLocalSpaceDimension) {
        KRATOS_ERROR_IF_NOT(ElementLocalSpaceDimension==2) << "A projection plane must be provided for ProjectVectorOnSurfaceUtility for element " << ElementID << std::endl;
    } else {
        if (ElementLocalSpaceDimension!=2) return true;
    }
    return false;
}
} // helpers namespace


void ProjectVectorOnSurfaceUtility::Execute(ModelPart& rModelPart, Parameters ThisParameters)
{
    ValidateParameters(ThisParameters);
    const int echo_level = ThisParameters["echo_level"].GetInt();

    const std::string& r_variable_name = ThisParameters["variable_name"].GetString();
    KRATOS_ERROR_IF_NOT(KratosComponents<ArrayVariableType>::Has(r_variable_name)) << "Variable " << r_variable_name << " not known" << std::endl;
    const ArrayVariableType& r_variable = KratosComponents<ArrayVariableType>::Get(r_variable_name);

    const Vector3 global_direction = CheckAndReadNormalizedVector(ThisParameters["global_direction"]);

    // std::cout << std::endl << "Assigning " << r_variable_name << " orientation to elements using method: " << r_projection_type << std::endl;

    const std::string& r_projection_type = ThisParameters["projection_type"].GetString();
    const auto method_specific_settings = ThisParameters["method_specific_settings"];
    const bool check_local_space_dimension = ThisParameters["check_local_space_dimension"].GetBool();

    if (r_projection_type == "planar") {
        PlanarProjection(rModelPart, method_specific_settings, global_direction, r_variable, echo_level, check_local_space_dimension);
    } else if (r_projection_type == "radial") {
        RadialProjection(rModelPart, method_specific_settings, global_direction, r_variable, echo_level, check_local_space_dimension);
    } else if (r_projection_type == "spherical") {
        SphericalProjection(rModelPart, method_specific_settings, global_direction, r_variable, echo_level, check_local_space_dimension);
    } else {
        KRATOS_ERROR << "projection type: " << r_projection_type << " not available, please use planar,radial,spherical" << std::endl;
    }

}

void ProjectVectorOnSurfaceUtility::PlanarProjection(
        ModelPart& rModelPart,
        const Parameters ThisParameters,
        const Vector3& rGlobalDirection,
        const ArrayVariableType& rVariable,
        const int EchoLevel,
        const bool rCheckLocalSpaceDimension)
{
     auto& r_process_info = rModelPart.GetProcessInfo();

    // Declare working variables
    Matrix local_coordinate_orientation;

    // Loop over all elements in part
    for (auto &element : rModelPart.Elements()) {

        if (CheckElementLocalSpaceDimension(rCheckLocalSpaceDimension,element.Id(),element.GetGeometry().LocalSpaceDimension()))
        {
            continue;
        }
        // get local axis in cartesian coordinates
        element.Calculate(LOCAL_ELEMENT_ORIENTATION, local_coordinate_orientation, r_process_info);

        Vector local_axis_1 = ZeroVector(3);
        Vector local_axis_2 = ZeroVector(3);
        Vector local_axis_3 = ZeroVector(3);

        local_axis_1 = column(local_coordinate_orientation,0);
        local_axis_2 = column(local_coordinate_orientation,1);
        local_axis_3 = column(local_coordinate_orientation,2);

        // normalise local axis vectors (global cartesian)
        local_axis_1 /= norm_2(local_axis_1);
        local_axis_2 /= norm_2(local_axis_2);
        local_axis_3 /= norm_2(local_axis_3);

        // (Abaqus default projection)
        // http://130.149.89.49:2080/v6.8/books/gsa/default.htm?startat=ch05s03.html
        // Shell local axis 1 is the projection of Global X vector onto the shell surface.
        // If the Global X vector is normal to the shell surface,
        // the shell local 1-direction is the projection of the
        // Global Z vector onto the shell surface

        // First, check if specified global_vector is normal to the shell surface
        if (std::abs(inner_prod(rGlobalDirection, local_axis_1)) < std::numeric_limits<double>::epsilon() && std::abs(inner_prod(rGlobalDirection, local_axis_2)) < std::numeric_limits<double>::epsilon()) {
            KRATOS_ERROR << "Global direction is perpendicular to element " << element.GetId() << " please define a different projection plane or use another type of projection "
                << ", available: radial,spherical" << std::endl;
        } else {
            // Second, project the global vector onto the shell surface
            // http://www.euclideanspace.com/maths/geometry/elements/plane/lineOnPlane/index.htm
            // vector to be projected = vec_a
            // Surface normal = vec_b
            const Vector& vec_a = rGlobalDirection;
            const Vector& vec_b = local_axis_3;

            Vector a_cross_b = ZeroVector(3);
            Vector projected_result = ZeroVector(3);

            MathUtils<double>::CrossProduct(a_cross_b, vec_a, vec_b);
            MathUtils<double>::CrossProduct(projected_result, vec_b, a_cross_b);
            //noramlize projected result
            projected_result /= MathUtils<double>::Norm(projected_result);

            element.SetValue(rVariable, projected_result);
        }
    }
}

void ProjectVectorOnSurfaceUtility::RadialProjection(
        ModelPart& rModelPart,
        const Parameters ThisParameters,
        const Vector3& rGlobalDirection,
        const ArrayVariableType& rVariable,
        const int EchoLevel,
        const bool rCheckLocalSpaceDimension)
{
    const auto& r_process_info = rModelPart.GetProcessInfo();

    // Declare working variables
    Matrix local_coordinate_orientation;

    // Loop over all elements in part
    for (auto &element : rModelPart.Elements()) {

        if (CheckElementLocalSpaceDimension(rCheckLocalSpaceDimension,element.Id(),element.GetGeometry().LocalSpaceDimension()))
        {
            continue;
        }
        // get local axis in cartesian coordinates
        element.Calculate(LOCAL_ELEMENT_ORIENTATION, local_coordinate_orientation, r_process_info);

        Vector local_axis_1 = ZeroVector(3);
        Vector local_axis_2 = ZeroVector(3);
        Vector local_axis_3 = ZeroVector(3);

        local_axis_1 = column(local_coordinate_orientation,0);
        local_axis_2 = column(local_coordinate_orientation,1);
        local_axis_3 = column(local_coordinate_orientation,2);

        // normalise local axis vectors (global cartesian)
        local_axis_1 /= norm_2(local_axis_1);
        local_axis_2 /= norm_2(local_axis_2);
        local_axis_3 /= norm_2(local_axis_3);

        // (Abaqus default projection)
        // http://130.149.89.49:2080/v6.8/books/gsa/default.htm?startat=ch05s03.html
        // Shell local axis 1 is the projection of Global X vector onto the shell surface.
        // If the Global X vector is normal to the shell surface,
        // the shell local 1-direction is the projection of the
        // Global Z vector onto the shell surface

        // First, check if specified rGlobalDirection is normal to the shell surface
        if (std::abs(inner_prod(rGlobalDirection, local_axis_1)) < std::numeric_limits<double>::epsilon() && std::abs(inner_prod(rGlobalDirection, local_axis_2)) < std::numeric_limits<double>::epsilon()) {
            KRATOS_ERROR << "Global direction is perpendicular to element " << element.GetId() << " please define a different projection plane or use another type of projection "
                << ", available: planar,spherical" << std::endl;
        } else {
            Vector projected_result = ZeroVector(3);

            MathUtils<double>::CrossProduct(projected_result, rGlobalDirection, local_axis_3);

            //noramlize projected result
            projected_result /= MathUtils<double>::Norm(projected_result);

            element.SetValue(rVariable, projected_result);
        }
    }
}

void ProjectVectorOnSurfaceUtility::SphericalProjection(
        ModelPart& rModelPart,
        const Parameters ThisParameters,
        const Vector3& rGlobalDirection,
        const ArrayVariableType& rVariable,
        const int EchoLevel,
        const bool rCheckLocalSpaceDimension)
{
    KRATOS_ERROR << "SphericalProjection not implemented" << std::endl;
}

} // namespace Kratos.
