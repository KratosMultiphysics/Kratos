//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Miguel Angel Celigueta
//                  Aditya Ghantasala
//

#include "fluid_post_process_utilities.h"


namespace Kratos {

    double FluidPostProcessUtilities::CalculateFlow(const ModelPart& rModelPart) {

        double flow = 0.0;
        auto conditions_begin = rModelPart.ConditionsBegin();
        const int num_conditions = rModelPart.NumberOfConditions();
        auto& r_data_comm = rModelPart.GetCommunicator().GetDataCommunicator();

        #pragma omp parallel for reduction(+:flow)
        for(int i_c=0; i_c<num_conditions;++i_c)
        {
            auto it_cond = conditions_begin + i_c;
            const auto& r_geom = it_cond->GetGeometry();
            const IndexType num_nodes = r_geom.size();

            const IntegrationMethodType integration_method = r_geom.GetDefaultIntegrationMethod();
            const auto& integration_points = r_geom.IntegrationPoints(integration_method);
            const Matrix& rNcontainer = r_geom.ShapeFunctionsValues(integration_method);
            const int num_integration_points = (int)integration_points.size();

            double average_velocity = 0.0;
            const double area = r_geom.Area();

            for (int point_number = 0; point_number < num_integration_points; ++point_number)
            {
                const auto point_local = integration_points[point_number];
                array_1d<double,3> unit_normal = r_geom.Normal(point_local);
                if( area < std::numeric_limits<double>::epsilon())
                    unit_normal *= 0.0;
                else
                    unit_normal *= 1/area;

                // Calculating the velocity on the gauss point
                Vector gauss_velocity = ZeroVector(3);
                for ( IndexType ii = 0; ii < num_nodes; ii++ ) {
                    gauss_velocity += rNcontainer( point_number, ii ) * r_geom[ii].FastGetSolutionStepValue( VELOCITY );
                }
                average_velocity += MathUtils<double>::Dot(gauss_velocity, unit_normal);
            }
            average_velocity /= num_integration_points;

            const double condition_flow = area * average_velocity;

            flow += condition_flow;
        } //for conditions

        // For MPI
        double total_flow = flow;
        if(r_data_comm.IsDistributed())
        {
            total_flow = r_data_comm.SumAll(flow);
        }

        return total_flow;
    }

} // namespace Kratos.
