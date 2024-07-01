// // KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
// //       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
// //      | (_| (_) | .` |\ V /___| |) | || _|| _|
// //       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
// //
// //  License:       BSD License
// //                 Kratos default license: kratos/license.txt
// //
// //  Main authors:  Riccardo Rossi
// //

#if !defined(KRATOS_SURFACE_INTEGRAL_INCLUDED)
#define KRATOS_SURFACE_INTEGRAL_INCLUDED

#define PRESSURE_ON_EULERIAN_MESH
#define USE_FEW_PARTICLES

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes
#ifdef _OPENMP
#include "omp.h"
#endif

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/variables.h"

namespace Kratos
{

    class ComputeFlux
    {
    public:
        KRATOS_CLASS_POINTER_DEFINITION(ComputeFlux);

        // typedef std::vector<const Variable<array_1d<double, 3>>*>& Vector3;

        typedef Node NodeType;
        typedef Geometry<NodeType> GeometryType;
        typedef ModelPart::NodesContainerType NodesArrayType;
        typedef ModelPart::ConditionsContainerType ConditionsArrayType;

        ComputeFlux(ModelPart &rModelPart, const Variable<array_1d<double, 3>> &rVar)
            : mrModelPart(rModelPart), mrVar(rVar)
        {
        }

        ~ComputeFlux()
        {
        }

        //**********************************************************************************************
        //**********************************************************************************************
        double ComputeSurfaceIntegral()
        {
            double flux = 0.0;
            // const ProcessInfo& r_current_process_info = mrModelPart.GetProcessInfo();

            block_for_each(mrModelPart.Conditions(), [&](Condition &rCond){
                GeometryType& r_geometry = rCond.GetGeometry();
                const size_t dim = r_geometry.WorkingSpaceDimension();

                const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod();
                const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method);
                const std::size_t number_of_integration_points = r_integration_points.size();
                const unsigned int number_of_nodes_element = r_geometry.size();

                // Values of the variable at the nodes
                Matrix nodal_var_values(number_of_nodes_element, dim);
                for (size_t i_node = 0; i_node < number_of_nodes_element; i_node++)
                {
                    Vector var_value = r_geometry[i_node].FastGetSolutionStepValue(mrVar);
                    for (size_t j = 0; j < dim; j++)
                    {
                        nodal_var_values(i_node, j) = var_value[j];
                    }
                }

                // Values of the shape functions
                const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(r_integration_method);
                // Value of the Jacobian
                Vector detJ0 = ZeroVector(number_of_integration_points);
                r_geometry.DeterminantOfJacobian(detJ0, r_integration_method);

                // Compute the values at the gauss points and compute the flux
                for ( IndexType i_gauss = 0; i_gauss < number_of_integration_points; ++i_gauss ) {

                    // Shape functions at the gauss point `i_point`
                    auto N = row(N_gausspoint, i_gauss);

                    // Gauss point values
                    const double gauss_point_volume = r_integration_points[i_gauss].Weight() * detJ0[i_gauss];
                    
                    // Unit normal
                    Vector unit_normal = r_geometry.UnitNormal(i_gauss);
                    unit_normal = unit_normal * (-1.0);
                    
                    // Sum the contribution of the flux of each gauss point for each space component
                    // detJ0mean += detJ0[i_gauss];
                    double element_flux = 0.;
                    for (size_t j = 0; j < dim; j++)
                    {
                        auto nodal_values_component_j = column(nodal_var_values, j);  // Vector with the j-th component of the nodes
                        double gauss_value_component_j = inner_prod(N, nodal_values_component_j);

                        flux += gauss_value_component_j * gauss_point_volume * unit_normal[j];
                        element_flux += gauss_value_component_j * gauss_point_volume * unit_normal[j];
                    }
                } }
            );

            return flux;
        }

    private:
        ModelPart &mrModelPart;
        const Variable<array_1d<double, 3>> &mrVar;
    };

} // namespace Kratos.

#endif // KRATOS_BFECC_LIMITER_CONVECTION_INCLUDED  defined
