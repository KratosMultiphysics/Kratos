//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Angel Celigueta

#include "post_process_utilities.h"


namespace Kratos {

    double PostProcessUtilities::ComputeFlow(const ModelPart& rModelPart) {

        double flow = 0.0;
        auto conditions_begin = rModelPart.ConditionsBegin();
        const int num_conditions = rModelPart.NumberOfConditions();

        #pragma omp parallel for reduction(+:flow)
        for(int i_c=0; i_c<num_conditions;++i_c) {

            auto it_cond = conditions_begin + i_c;
            const auto& r_geom = it_cond->GetGeometry();

            const double area = r_geom.Area();
            GeometryType::CoordinatesArrayType point_local;
            r_geom.PointLocalCoordinates(point_local, r_geom.Center()) ;
            array_1d<double,3> unit_normal = r_geom.Normal(point_local);

            if( area < std::numeric_limits<double>::epsilon())
                unit_normal *= 0.0;
            else
                unit_normal *= 1/area;

            double average_velocity = 0.0;

            for (int j=0; j<(int)r_geom.size(); j++) {
                average_velocity += MathUtils<double>::Dot(r_geom[j].FastGetSolutionStepValue(VELOCITY), unit_normal);
            }
            average_velocity *= (1.0 / (double)r_geom.size());

            const double condition_flow = area * average_velocity;

            flow += condition_flow;
        } //for conditions

        return flow;
    }

} // namespace Kratos.
