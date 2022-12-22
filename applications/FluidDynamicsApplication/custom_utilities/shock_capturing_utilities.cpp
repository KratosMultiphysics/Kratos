//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Uxue Chasco
//
//

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h"
#include "utilities/element_size_calculator.h"

// Application includes
#include "shock_capturing_utilities.h"
#include "fluid_dynamics_application_variables.h"


namespace Kratos
{

    /* Public functions *******************************************************/

    /**
     * @brief Calculation of the Artificial viscosity according to the documentation
     */
    template <std::size_t TDim, std::size_t TNumNodes>
    void ShockCapturingUtilities<TDim, TNumNodes>::ExecuteRMethodType(const ModelPart &mrModelPart)
    {
        double alg_constant_shock_capturing = 0.8;
        for (auto it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); ++it_elem)
        {
            auto &r_geom = it_elem->GetGeometry();
            const auto& r_process_info = mrModelPart.GetProcessInfo();

            //  Velocity residual norm
            double residual_velocity_norm;
            it_elem->Calculate(VELOCITY_RESIDUAL_NORM, residual_velocity_norm, r_process_info);

            // Previous time step velocity gradient
            double Volume;
            array_1d<double, TNumNodes> N;
            BoundedMatrix<double, TNumNodes, TDim> DN_DX;
            GeometryUtils::CalculateGeometryData(r_geom, DN_DX, N, Volume);

            // Velocity gradient
            double velocity_gradient_norm;
            BoundedMatrix<double, TDim, TDim>
            velocity_gradient = ZeroMatrix(TDim, TDim);
            for (unsigned int n = 0; n < TNumNodes; n++)
            {
                const auto &r_vel = r_geom[n].FastGetSolutionStepValue(VELOCITY);
                for (unsigned int i = 0; i < TDim; i++){
                    for (unsigned int j = 0; j < TDim; j++)
                    {
                        velocity_gradient(i, j) += DN_DX(n, j) * r_vel[i];
                        velocity_gradient_norm += std::pow(velocity_gradient(i, j),2);
                    }
                }
            }
            // velocity gradient normal
            velocity_gradient_norm = sqrt(velocity_gradient_norm);
            const double h = ElementSizeCalculator<TDim, TNumNodes>::GradientsElementSize(DN_DX);

            // artificial viscosity according Ramon Method
            const double artificial_viscosity = 0.5 * (alg_constant_shock_capturing * h) * (residual_velocity_norm / velocity_gradient_norm);
            it_elem->SetValue(ARTIFICIAL_MASS_DIFFUSIVITY, velocity_gradient_norm);
            it_elem->SetValue(ARTIFICIAL_CONDUCTIVITY, residual_velocity_norm);
            it_elem->SetValue(ARTIFICIAL_DYNAMIC_VISCOSITY, artificial_viscosity);

            // TODO: We should use the GetValue in here. Do the SetValue first somewhere.
        }
    }

    template class ShockCapturingUtilities<2, 3>;
    template class ShockCapturingUtilities<3, 4>;
}
