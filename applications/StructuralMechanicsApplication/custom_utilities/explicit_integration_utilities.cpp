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
#include "custom_utilities/explicit_integration_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
namespace ExplicitIntegrationUtilities
{
double CalculateDeltaTime(
    ModelPart& rModelPart,
    const double PredictionLevel,
    const double MaximumDeltaTime,
    const double SafetyFactor
    )
{
    KRATOS_TRY

    // Getting process info
    ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
    ElementsArrayType& r_elements = rModelPart.Elements();

    // Initial delta time
    double delta_time = MaximumDeltaTime / SafetyFactor;

    // Initialize the value
    double stable_delta_time = 1000.0;

    // Auxiliar parameters
    bool check_has_all_variables = true;
    double E(0.0), nu(0.0), roh(0.0), alpha(0.0), beta(0.0);

    // Iterate over elements
    const auto it_elem_begin = rModelPart.ElementsBegin();
    #pragma omp parallel for firstprivate(check_has_all_variables, E, nu, roh, alpha, beta)
    for (int i = 0; i < static_cast<int>(r_elements.size()); ++i) {
        auto it_elem = it_elem_begin + i;

        /* Get geometric and material properties */
        const Properties& r_properties = it_elem->GetProperties();

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
            nu = 0.0;
        }

        // Getting density
        if (r_properties.Has(DENSITY)) {
            roh = r_properties[DENSITY];
        } else {
            check_has_all_variables = false;
        }

        if (check_has_all_variables) {
            const double length = it_elem->GetGeometry().Length();

            // Compute courant criterion
            const double bulk_modulus = E / (3.0 * (1.0 - 2.0 * nu));
            const double wavespeed = std::sqrt(bulk_modulus / roh);
            const double w = 2.0 * wavespeed / length; // Frequency

            const double psi = 0.5 * (alpha / w + beta * w); // Critical ratio;
            stable_delta_time = (2.0 / w) * (std::sqrt(1.0 + psi * psi) - psi);

            if (stable_delta_time > 0.0) {
                #pragma omp critical
                if (stable_delta_time < delta_time) delta_time = stable_delta_time;
            }
        } else {
            KRATOS_ERROR << "Not enough parameters for prediction level " << PredictionLevel << std::endl;
        }
    }

    stable_delta_time = delta_time * SafetyFactor;

    if (stable_delta_time < MaximumDeltaTime) {
        r_current_process_info[DELTA_TIME] = stable_delta_time;
    }

    KRATOS_INFO_IF("ExplicitIntegrationUtilities", PredictionLevel > 1)
    << "  [EXPLICIT PREDICTION LEVEL " << PredictionLevel << " ] : (computed stable time step = " << stable_delta_time << " s)\n"
    << "  Using  = " << r_current_process_info[DELTA_TIME] << " s as time step DELTA_TIME)" << std::endl;

    return stable_delta_time;

    KRATOS_CATCH("")
}

} // namespace ExplicitIntegrationUtilities
} // namespace Kratos
