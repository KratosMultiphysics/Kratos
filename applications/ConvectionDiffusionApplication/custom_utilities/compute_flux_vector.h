// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//           Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_COMPUTE_FLUX_VECTOR_UTILITY_INCLUDED)
#define KRATOS_COMPUTE_FLUX_VECTOR_UTILITY_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <unordered_map>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"

namespace Kratos
{

class ComputeFluxUtility
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ComputeFluxUtility);

    using NodeType = Node;

    /**
     * @brief Computes and stores the flux in the given model part.
     *
     * The nodal flux is computed as:
     * flux = -k * grad(phi)
     * where grad(phi) is computed from the unknown variable using shape function derivatives
     * and k is the conductivity from the diffusion_variable.
     *
     * The computed vector value is stored in the historical value of
     * @p rVectorFluxVar for each node.
     */
    static void ComputeVectorialFlux(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVectorFluxVar
    )
    {
        KRATOS_TRY;

        const auto& r_process_info = rModelPart.GetProcessInfo();
        const auto& r_settings = *r_process_info[CONVECTION_DIFFUSION_SETTINGS];
        const auto& r_unknown_var = r_settings.GetUnknownVariable();
        const auto& r_diffusion_variable = r_settings.GetDiffusionVariable();

        // Initialize nodal gradient accumulators
        std::unordered_map<std::size_t, array_1d<double, 3>> nodal_gradient;
        std::unordered_map<std::size_t, double> nodal_measure;
        
        for (auto& r_node : rModelPart.Nodes()) {
            nodal_gradient[r_node.Id()] = ZeroVector(3);
            nodal_measure[r_node.Id()] = 0.0;
        }

        // Iterate over elements and compute gradient contributions
        for (auto& r_element : rModelPart.Elements()) {
            auto& r_geometry = r_element.GetGeometry();
            const std::size_t num_nodes = r_geometry.size();
            const std::size_t dim = r_geometry.WorkingSpaceDimension();

            // Get integration method and points
            const auto integration_method = r_geometry.GetDefaultIntegrationMethod();
            const auto& r_integration_points = r_geometry.IntegrationPoints(integration_method);
            const auto& r_shape_functions = r_geometry.ShapeFunctionsValues(integration_method);
            
            // Compute det(J)
            Vector det_j = ZeroVector(r_integration_points.size());
            r_geometry.DeterminantOfJacobian(det_j, integration_method);
            
            // Sum contributions over integration points
            for (std::size_t i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
                const double gauss_weight = r_integration_points[i_gauss].Weight() * det_j[i_gauss];
                
                // Get shape function derivatives at Gauss point
                const auto& r_DN_De = r_geometry.ShapeFunctionDerivatives(1, i_gauss, integration_method);
                
                // Compute dN/dX from dN/dE and Jacobian
                Matrix DN_DX = ZeroMatrix(num_nodes, dim);
                Matrix J = ZeroMatrix(dim, dim);
                Matrix J_inv = ZeroMatrix(dim, dim);
                double det_J;
                
                r_geometry.Jacobian(J, i_gauss, integration_method);
                MathUtils<double>::GeneralizedInvertMatrix(J, J_inv, det_J);
                
                noalias(DN_DX) = prod(r_DN_De, J_inv);
                
                // Get nodal values of unknown variable
                Vector phi_values = ZeroVector(num_nodes);
                for (std::size_t i_node = 0; i_node < num_nodes; ++i_node) {
                    phi_values[i_node] = r_geometry[i_node].FastGetSolutionStepValue(r_unknown_var);
                }
                
                // Compute gradient at Gauss point: grad = dN/dX^T * phi
                array_1d<double, 3> gauss_grad = ZeroVector(3);
                for (std::size_t i = 0; i < dim; ++i) {
                    for (std::size_t j = 0; j < num_nodes; ++j) {
                        gauss_grad[i] += DN_DX(j, i) * phi_values[j];
                    }
                }
                
                // Accumulate gradient contribution to nodes (weighted by integration measure)
                for (std::size_t i_node = 0; i_node < num_nodes; ++i_node) {
                    const auto node_id = r_geometry[i_node].Id();
                    const double weight_per_node = gauss_weight / static_cast<double>(num_nodes);
                    nodal_gradient[node_id] += gauss_grad * weight_per_node;
                    nodal_measure[node_id] += weight_per_node;
                }
            }
        }

        // Normalize gradient and compute flux at each node
        for (auto& r_node : rModelPart.Nodes()) {
            const auto node_id = r_node.Id();
            
            // Normalize gradient by measure
            array_1d<double, 3> gradient = ZeroVector(3);
            if (nodal_measure[node_id] > 0.0) {
                gradient = nodal_gradient[node_id] / nodal_measure[node_id];
            }
            
            // Get conductivity and compute flux
            const double conductivity = r_node.FastGetSolutionStepValue(r_diffusion_variable);
            auto& r_flux = r_node.FastGetSolutionStepValue(rVectorFluxVar);
            r_flux = -conductivity * gradient;
        }

        KRATOS_CATCH("");
    }
};

} // namespace Kratos

#endif // KRATOS_COMPUTE_FLUX_VECTOR_UTILITY_INCLUDED
