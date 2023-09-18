//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/checks.h"
#include "includes/model_part.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/variable_utils.h"

// Application includes
#include "gauss_point_error_utility.h"

namespace Kratos
{

/* Public functions *******************************************************/
Parameters GaussPointErrorUtility::GetDefaultSettings()
{
    Parameters default_settings(R"({
        "model_part_name" : "",
        "variable_name" : ""
    })");

    return default_settings;
}

GaussPointErrorUtility::GaussPointErrorUtility(
    ModelPart& rModelPart,
    Parameters rParameters)
    : mrModelPart(rModelPart)
{
}

GaussPointErrorUtility::GaussPointErrorUtility(
    Model& rModel,
    Parameters rParameters)
    : GaussPointErrorUtility(
        [&] (Model& x, Parameters& y) -> ModelPart& {
            y.ValidateAndAssignDefaults(GetDefaultSettings());
            KRATOS_ERROR_IF(y["model_part_name"].GetString() == "") << "\'model_part_name\' is empty. Please provide the origin model part name." << std::endl;
            return x.GetModelPart(y["model_part_name"].GetString());
        } (rModel, rParameters),
        rParameters)
{
}

double GaussPointErrorUtility::Execute()
{
    // Initialize values
    double disp_error_norm = 0.0;

    // Loop the elements to calculate the error
    const int n_elems = mrModelPart.NumberOfElements();
    for (int i_elem = 0; i_elem < n_elems; ++i_elem) {
        auto it_elem = mrModelPart.ElementsBegin() + i_elem;
        const auto p_geom = it_elem->pGetGeometry();
        const int n_nodes = p_geom->PointsNumber();

        // Initialize values
        Matrix r_N_container;
        const GeometryData::IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
        const auto& r_integrations_points = p_geom->IntegrationPoints(integration_method);

        array_1d<double,3> i_disp_exact;
        array_1d<double,3> i_gauss_disp_err;
        array_1d<double,3> i_gauss_disp_exact;

        r_N_container = p_geom->ShapeFunctionsValues(integration_method);
        Vector det_J_vector(r_integrations_points.size());
        p_geom->DeterminantOfJacobian(det_J_vector, integration_method);

        const unsigned int n_gauss = r_integrations_points.size();
        for (int i_gauss = 0; i_gauss < n_gauss; ++i_gauss ) {
            // Get Gauss pt. data
            const Vector& rN = row(r_N_container, i_gauss);
            const double w = r_integrations_points[i_gauss].Weight() * det_J_vector[i_gauss];

            // Calculate the Gauss pt. coordinates
            array_1d<double,3> coords = ZeroVector(3);
            coords[0] = ((*p_geom)[0]).X()*rN[0] + ((*p_geom)[1]).X()*rN[1] + ((*p_geom)[2]).X()*rN[2];
            coords[1] = ((*p_geom)[0]).Y()*rN[0] + ((*p_geom)[1]).Y()*rN[1] + ((*p_geom)[2]).Y()*rN[2];

            // Interpolate the Gauss pt. error
            i_gauss_disp_err = ZeroVector(3);
            i_gauss_disp_exact = ZeroVector(3);
            for (int i_node = 0; i_node < n_nodes; ++i_node) {
                const auto& r_node = (*p_geom)[i_node];
                CalculateExactSolution(r_node.Coordinates(), i_disp_exact);
                i_gauss_disp_err += rN[i_node] * r_node.FastGetSolutionStepValue(DISPLACEMENT);
                i_gauss_disp_exact += rN[i_node] * i_disp_exact;
            }
            i_gauss_disp_err -= i_gauss_disp_exact;

            // Add the current Gauss pt. error to the total norm
            disp_error_norm += w * std::pow(norm_2(i_gauss_disp_err), 2);
        }
    }

    disp_error_norm = std::sqrt(disp_error_norm);
    return disp_error_norm;
}

int GaussPointErrorUtility::Check()
{
    // Check that model part is not empty
    KRATOS_ERROR_IF(mrModelPart.NumberOfNodes() == 0) << "There are no nodes in the origin model part." << std::endl;
    KRATOS_ERROR_IF(mrModelPart.NumberOfElements() == 0) << "There are no elements in the origin model part." << std::endl;

    return 0;
}

/* Protected functions ****************************************************/

/* Private functions ******************************************************/

// Manufactured solution for the "annulus" plate test
void GaussPointErrorUtility::CalculateExactSolution(
    const array_1d<double,3>& rCoords,
    array_1d<double,3>& rExactSolution)
{
    const double x = rCoords[0];
    const double y = rCoords[1];

    // rExactSolution[0] = strain_factor*x;
    // rExactSolution[1] = strain_factor*y;
    // rExactSolution[2] = 0.0;

    // rExactSolution[0] = strain_factor*x*y;
    // rExactSolution[1] = strain_factor*x*y;
    // rExactSolution[2] = 0.0;

    rExactSolution[0] = 1.0e-3*x*y*std::sin(2.0*Kratos::Globals::Pi*x)*std::cos(2.0*Kratos::Globals::Pi*y);
    rExactSolution[1] = 1.0e-3*x*y*std::sin(2.0*Kratos::Globals::Pi*y)*std::cos(2.0*Kratos::Globals::Pi*x);
    rExactSolution[2] = 0.0;
}

};  // namespace Kratos.
