// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Mahmoud Sesa, https://github.com/mahmoudsesa
//

// System includes

// External includes

// Project includes
#include "adjoint_nonlinear_strain_energy_response_function.h"

namespace Kratos 
{
    AdjointNonlinearStrainEnergyResponseFunction::AdjointNonlinearStrainEnergyResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    {
    }

    AdjointNonlinearStrainEnergyResponseFunction::~AdjointNonlinearStrainEnergyResponseFunction()
    {
    }

// Calculates the response value increment for one time step during the primal analysis
// It is not tested for problems involving self weight
// TODO changing the variable names to follow the guidline
void AdjointNonlinearStrainEnergyResponseFunction::CalculateResponseIncrement(ModelPart& rModelPart)
{
    KRATOS_TRY;

    ProcessInfo &r_current_process_info = rModelPart.GetProcessInfo();
    double response_increment_value = 0.0;

    // Check if there is adjoint elements, because response function calculation is done for primal analysis
    KRATOS_ERROR_IF( r_current_process_info.Has(IS_ADJOINT) && r_current_process_info[IS_ADJOINT] )
    << "Calculate value for strain energy response is only available when using primal elements" << std::endl;

    // sum all elemental strain energy increment values calculated by trapezoidal rule: E = 0.5 * (f_ext_i - f_ext_i-1) * (u_i - u_i-1)
    Matrix LHS;
    Vector RHS;
    Vector disp;
    Vector external_force;
    Vector external_force_previous_step;
    Vector disp_previous_step;
    Vector disp_increment;
    Vector average_load;

    for (auto& elem_i : rModelPart.Elements())
    {
        elem_i.GetValuesVector(disp,0);
        elem_i.GetValuesVector(disp_previous_step, 1);
        elem_i.CalculateLocalSystem(LHS, RHS, r_current_process_info);

        disp_increment = disp - disp_previous_step;
        external_force = -1.0 * RHS;
        external_force_previous_step = external_force - prod(LHS , disp_increment);
        average_load = 0.5 * (external_force + external_force_previous_step);

        response_increment_value += inner_prod(average_load , disp_increment);
    }

    m_response_value += response_increment_value; 

    KRATOS_CATCH("");
}

double AdjointNonlinearStrainEnergyResponseFunction::CalculateValue(ModelPart& rModelPart)
{
    KRATOS_TRY;

    this->CheckForBodyForces(rModelPart);

    return m_response_value;

    KRATOS_CATCH("");
}

void AdjointNonlinearStrainEnergyResponseFunction::CheckForBodyForces(ModelPart& rModelPart)
{
    const double numerical_limit = std::numeric_limits<double>::epsilon();

    for (auto& elem_i : rModelPart.Elements())
    {
    Vector acc = ZeroVector(3);
    if (elem_i.GetProperties().Has( VOLUME_ACCELERATION ))
        acc+= elem_i.GetProperties()[VOLUME_ACCELERATION];

    if( elem_i.GetGeometry()[0].SolutionStepsDataHas(VOLUME_ACCELERATION))
        acc += elem_i.GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION);

    KRATOS_ERROR_IF( norm_2(acc)>numerical_limit )
        << "Calculation of nonlinear strain energy response is untested for self weight!" << std::endl;
    }
}

// TODO, Discussion of the implementation of the response gradient
// Test the effieciency of calculatevalue() implementation,
// Implementation of the sensitivities, this should be independent from the response function implementation

} // namespace Kratos