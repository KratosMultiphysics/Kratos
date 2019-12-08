// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Klaus B Sautter
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "geometries/point.h"
#include "custom_utilities/explicit_integration_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
namespace ExplicitIntegrationUtilities
{
double CalculateDeltaTime(
    ModelPart& rModelPart,
    Parameters ThisParameters
    )
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
    {
        "time_step_prediction_level" : 2.0,
        "max_delta_time"             : 1.0e0,
        "safety_factor"              : 0.80,
        "mass_factor"                : 1.0,
        "desired_delta_time"         : -1.0,
        "max_number_of_iterations"   : 10
    })" );


    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    const double time_step_prediction_level = ThisParameters["time_step_prediction_level"].GetDouble(); // The prediction level
    const double max_delta_time = ThisParameters["max_delta_time"].GetDouble(); // The prediction level
    const double safety_factor = ThisParameters["safety_factor"].GetDouble(); // The factor to not consider exactly the theoretical value
    double mass_factor = ThisParameters["mass_factor"].GetDouble(); // How the density of the solid is going to be multiplied (1.0 by default)
    const double desired_delta_time = ThisParameters["desired_delta_time"].GetDouble(); // The minimum delta time we want, if the value is negative not mass factor will be computed
    const bool compute_mass_factor = desired_delta_time < 0.0 ? false : true;
    const bool max_number_of_iterations = ThisParameters["max_number_of_iterations"].GetInt();

    // Getting process info
    ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

    // Initialize the value
    double stable_delta_time = 1000.0;

    // Actaully compute the value
    if (compute_mass_factor) {
        int iteration = 1;
        stable_delta_time = InnerCalculateDeltaTime(rModelPart, time_step_prediction_level, max_delta_time, safety_factor, mass_factor);
        if (stable_delta_time < desired_delta_time) {
            while (iteration < max_number_of_iterations) {
                mass_factor *= std::pow(desired_delta_time/stable_delta_time, 2);
                stable_delta_time = InnerCalculateDeltaTime(rModelPart, time_step_prediction_level, max_delta_time, safety_factor, mass_factor);
                KRATOS_INFO("ExplicitIntegrationUtilities") << "ITERATION NUMBER: " << iteration << "\tMass factor: " << mass_factor << "\nCurrent delta time: " << stable_delta_time << "\nDesired delta time:" << desired_delta_time << "\n RATIO: " << stable_delta_time/desired_delta_time << std::endl;
                if (stable_delta_time > desired_delta_time) {
                    break;
                }
                ++iteration;
            }
        }
    } else {
        stable_delta_time = InnerCalculateDeltaTime(rModelPart, time_step_prediction_level, max_delta_time, safety_factor, mass_factor);
    }

    if (stable_delta_time < max_delta_time) {
        r_current_process_info[DELTA_TIME] = stable_delta_time;
    }

    KRATOS_INFO_IF("ExplicitIntegrationUtilities", time_step_prediction_level > 1)
    << "  [EXPLICIT PREDICTION LEVEL " << time_step_prediction_level << " ] : (computed stable time step = " << stable_delta_time << " s)\n"
    << "  Using  = " << r_current_process_info[DELTA_TIME] << " s as time step DELTA_TIME)" << std::endl;

    return stable_delta_time;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

double InnerCalculateDeltaTime(
    ModelPart& rModelPart,
    const double TimeStepPredictionLevel,
    const double MaxDeltaTime,
    const double SafetyFactor,
    const double MassFactor
    )
{
    KRATOS_TRY

    // Initial delta time
    double delta_time = MaxDeltaTime / SafetyFactor;

    // Initialize the value
    double stable_delta_time = 1000.0;

    // Auxiliar parameters
    bool check_has_all_variables = true;
    double E(0.0), nu(0.0), rho(0.0), alpha(0.0), beta(0.0);

    // Iterate over elements
    ElementsArrayType& r_elements = rModelPart.Elements();
    const auto it_elem_begin = r_elements.begin();
    #pragma omp parallel for firstprivate(check_has_all_variables, stable_delta_time, E, nu, rho, alpha, beta)
    for (int i = 0; i < static_cast<int>(r_elements.size()); ++i) {
        auto it_elem = it_elem_begin + i;

        /* Get geometric and material properties */
        const Properties& r_properties = it_elem->GetProperties();
        auto& r_geometry = it_elem->GetGeometry();

        // Initialize
        check_has_all_variables = true;

        // Getting alpha Rayleigh parameter
        if (r_properties.Has(RAYLEIGH_ALPHA)) {
            alpha = r_properties[RAYLEIGH_ALPHA];
        } else {
            alpha = 0.0;
        }

        // Getting beta Rayleigh parameter
        if (r_properties.Has(RAYLEIGH_BETA)) {
            beta = r_properties[RAYLEIGH_BETA];
        } else {
            beta = 0.0;
        }

        // Getting Young modulus
        if (r_properties.Has(YOUNG_MODULUS)) {
            E = r_properties[YOUNG_MODULUS];
        } else {
            check_has_all_variables = false;
        }

        // Getting Poisson ratio
        if (r_properties.Has(POISSON_RATIO)) {
            nu = r_properties[POISSON_RATIO];
        } else {
            nu = -1.0;
        }

        // Getting density
        if (r_properties.Has(DENSITY)) {
            rho = MassFactor * r_properties[DENSITY];
        } else {
            check_has_all_variables = false;
        }

        if (check_has_all_variables) {
            // Computing length as the smallest side of the geometry
//             const double length = it_elem->GetGeometry().Length();
            double min_length = std::numeric_limits<double>::max();
            const auto edges = r_geometry.GenerateEdges();
            for (IndexType i_edge = 0; i_edge < r_geometry.EdgesNumber(); ++i_edge) {
                min_length = std::min(edges[i_edge].Length(), min_length);
            }

            // We compute the minimum height of the face too
            const auto faces = r_geometry.GenerateFaces();
            double max_length = 0.0;
            for (IndexType i_face = 0; i_face < r_geometry.FacesNumber(); ++i_face) {

                const auto sub_edges = faces[i_face].GenerateEdges();
                for (IndexType i_edge = 0; i_edge < faces[i_face].EdgesNumber(); ++i_edge) {
                    max_length = std::max(sub_edges[i_edge].Length(), max_length);
                }

                min_length = std::min(faces[i_face].Area()/max_length, min_length);
            }

            // Compute courant criterion
            const double bulk_modulus = (nu < 0.0) ? E : E / (3.0 * (1.0 - 2.0 * nu));
            const double wavespeed = std::sqrt(bulk_modulus / rho);
            const double w = 2.0 * wavespeed / min_length; // Frequency

            const double psi = 0.5 * (alpha / w + beta * w); // Critical ratio;
            stable_delta_time = (2.0 / w) * (std::sqrt(1.0 + psi * psi) - psi);

            if (stable_delta_time > 0.0) {
                #pragma omp critical
                if (stable_delta_time < delta_time) delta_time = stable_delta_time;
            }
        } else {
            KRATOS_ERROR << "Not enough parameters for prediction level " << TimeStepPredictionLevel << std::endl;
        }
    }

    stable_delta_time = delta_time * SafetyFactor;

    return stable_delta_time;

    KRATOS_CATCH("")
}

} // namespace ExplicitIntegrationUtilities
} // namespace Kratos
