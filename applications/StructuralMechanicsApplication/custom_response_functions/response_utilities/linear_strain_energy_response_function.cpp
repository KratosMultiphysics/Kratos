// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//

// System includes

// External includes

// Project includes
#include "linear_strain_energy_response_function.h"

namespace Kratos {

LinearStrainEnergyResponseFunction::LinearStrainEnergyResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
: mrModelPart(rModelPart), mResponseSettings(ResponseSettings)
{
    // TODO validate settings

    // Check if there are at primal elements, because the primal state is required
    ProcessInfo &r_current_process_info = mrModelPart.GetProcessInfo();
    KRATOS_ERROR_IF( r_current_process_info.Has(IS_ADJOINT) && r_current_process_info[IS_ADJOINT] )
            << "LinearStrainEnergyResponseFunction: Can not use adjoint model part!" << std::endl;
}

double LinearStrainEnergyResponseFunction::CalculateValue()
{
    KRATOS_TRY;

    ProcessInfo &r_current_process_info = mrModelPart.GetProcessInfo();
    double response_value = 0.0;

    // Sum all elemental strain energy values calculated as: W_e = u_e^T K_e u_e
    Matrix LHS;
    Vector RHS;
    Vector disp;

    for (auto& elem_i : mrModelPart.Elements())
    {
        elem_i.GetValuesVector(disp,0);

        elem_i.CalculateLocalSystem(LHS, RHS, r_current_process_info);

        // Compute linear strain energy 0.5*u*K*u
        response_value += 0.5 * inner_prod(disp, prod(LHS, disp));
        }

    return response_value;

    KRATOS_CATCH("");
}


} // namespace Kratos.

