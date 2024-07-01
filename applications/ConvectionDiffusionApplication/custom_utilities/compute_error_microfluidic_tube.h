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

#if !defined(ERROR_MICROFLUIDIC_TUBE_INCLUDED)
#define ERROR_MICROFLUIDIC_TUBE_INCLUDED

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

    class ComputeError
    {
    public:
        KRATOS_CLASS_POINTER_DEFINITION(ComputeError);

        // typedef std::vector<const Variable<array_1d<double, 3>>*>& Vector3;

        typedef Node<3> NodeType;
        typedef Geometry<NodeType> GeometryType;
        typedef ModelPart::NodesContainerType NodesArrayType;
        typedef ModelPart::ConditionsContainerType ConditionsArrayType;

        ComputeError(ModelPart &rModelPart)
            : mrModelPart(rModelPart)
        {
        }

        ~ComputeError()
        {
        }

        //**********************************************************************************************
        //**********************************************************************************************
        double ComputeErrorSurfaceIntegral()
        {
            double flux = 0.0;
            // const ProcessInfo& r_current_process_info = mrModelPart.GetProcessInfo();

            block_for_each(mrModelPart.Conditions(), [&](Condition &rCond)
                           {
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
                Vector vel_value = r_geometry[i_node].FastGetSolutionStepValue(VELOCITY);

                double temp_th_value = r_geometry[i_node].Z() > 0. ? 0. : 1.;
                double temp_node_value = r_geometry[i_node].FastGetSolutionStepValue(TEMPERATURE);
                temp_node_value = (temp_node_value - 1070.) / (1085. - 1070.);
                double diff_value = abs(temp_th_value - temp_node_value);
                // std::cout << "z = " << r_geometry[i_node].Z() << ", c_th = " << temp_th_value << ", c = " << temp_node_value << std::endl;
                for (size_t j = 0; j < dim; j++)
                {
                    // nodal_var_values(i_node, j) = diff_value * abs(vel_value[j]);
                    nodal_var_values(i_node, j) = 1.;
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
            } });

            return flux;
        }

    private:
        ModelPart &mrModelPart;
    };

} // namespace Kratos.

#endif //   defined
