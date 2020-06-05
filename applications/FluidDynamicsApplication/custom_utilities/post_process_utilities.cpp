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
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/model_part.h"

namespace Kratos {
    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;

    void PostProcessUtilities::ComputeFlow(const ModelPart& rModelPart, double& flow) {

        flow = 0.0;
        for(auto i_cond = rModelPart.ConditionsBegin(); i_cond != rModelPart.ConditionsEnd(); ++i_cond) {
            const auto& geom = i_cond->GetGeometry();

            const double area = geom.Area();
            GeometryType::CoordinatesArrayType point_local;
            geom.PointLocalCoordinates(point_local, geom.Center()) ;
            array_1d<double,3> unitary_normal = geom.Normal(point_local);
            unitary_normal *= (1.0 / area);

            double average_velocity = 0.0;

            for (int j=0; j<(int)geom.size(); j++) {
                average_velocity += MathUtils<double>::Dot(geom[j].FastGetSolutionStepValue(VELOCITY), unitary_normal);
            }
            average_velocity *= (1.0 / (double)geom.size());

            const double condition_flow = area * average_velocity;

            flow += condition_flow;
        } //for conditions

    }

} // namespace Kratos.
