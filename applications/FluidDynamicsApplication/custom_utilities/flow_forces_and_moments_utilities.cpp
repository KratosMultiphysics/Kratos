//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "utilities/openmp_utils.h"
#include "utilities/variable_utils.h"

// Application includes
#include "flow_forces_and_moments_utilities.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
    /* Public functions *******************************************************/

    std::tuple<array_1d<double,3>, array_1d<double,3>> FlowForcesAndMomentsUtilities::CalculateBodyFittedFlowForcesAndMoments(ModelPart& rModelPart, const array_1d<double,3>& rReferencePoint){
        // Sum the reactions in the model part of interest.
        // Note that the reactions are assumed to be already computed.
        VariableUtils variable_utils;

        array_1d<double,3> flow_force =variable_utils.SumHistoricalVariable<array_1d<double,3>>(REACTION, rModelPart, 0);

        flow_force *= -1.0;

        array_1d<double,3> flow_moment = ZeroVector(3);

        for (auto& r_node : rModelPart.Nodes()) {

            const auto& r_reaction = r_node.FastGetSolutionStepValue(REACTION);

            // position vector from reference point
            array_1d<double,3> r;
            noalias(r) = r_node.Coordinates() - rReferencePoint;

            // moment = r x F
            array_1d<double,3> nodal_moment;
            MathUtils<double>::CrossProduct(nodal_moment, r, -r_reaction);

            flow_moment += nodal_moment;
        }

        return std::make_tuple(flow_force, flow_moment);
    }

    std::tuple<array_1d<double,3>, array_1d<double,3>> FlowForcesAndMomentsUtilities::CalculateEmbeddedFlowForcesAndMoments(
        ModelPart& rModelPart,
        const array_1d<double, 3>& rReferencePoint)
    {
        // Initialize total flow force
        array_1d<double, 3> flow_force = ZeroVector(3);
        double& flow_force_x = flow_force[0];
        double& flow_force_y = flow_force[1];
        double& flow_force_z = flow_force[2];

        // Initialize total flow moment
        array_1d<double, 3> flow_moment = ZeroVector(3);
        double& flow_moment_x = flow_moment[0];
        double& flow_moment_y = flow_moment[1];
        double& flow_moment_z = flow_moment[2];

        // Auxiliary vars
        array_1d<double, 3> elem_force;

        double force_x_red  = 0.0;
        double force_y_red  = 0.0;
        double force_z_red  = 0.0;

        double moment_x_red = 0.0;
        double moment_y_red = 0.0;
        double moment_z_red = 0.0;

        #pragma omp parallel for reduction(+:force_x_red,force_y_red,force_z_red,moment_x_red,moment_y_red,moment_z_red) private(elem_force) schedule(dynamic)
        for(int i = 0; i < static_cast<int>(rModelPart.Elements().size()); ++i){

            auto it_elem = rModelPart.ElementsBegin() + i;

            it_elem->Calculate(DRAG_FORCE, elem_force, rModelPart.GetProcessInfo());

            // ---- Geometry center ----
            const array_1d<double,3>& center =
                it_elem->GetGeometry().Center();

            // ---- Lever arm ----
            array_1d<double,3> lever_arm;
            noalias(lever_arm) = center - rReferencePoint;

            // ---- Moment contribution ----
            array_1d<double,3> elem_moment;
            MathUtils<double>::CrossProduct(elem_moment, lever_arm, elem_force);

            // ---- Reductions ----
            force_x_red  += elem_force[0];
            force_y_red  += elem_force[1];
            force_z_red  += elem_force[2];

            moment_x_red += elem_moment[0];
            moment_y_red += elem_moment[1];
            moment_z_red += elem_moment[2];
        }

        flow_force_x  += force_x_red;
        flow_force_y  += force_y_red;
        flow_force_z  += force_z_red;

        flow_moment_x += moment_x_red;
        flow_moment_y += moment_y_red;
        flow_moment_z += moment_z_red;

        // MPI synchronization
        flow_force  = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(flow_force);
        flow_moment = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(flow_moment);

        return std::make_tuple(flow_force, flow_moment);
    }

    array_1d<double, 3> FlowForcesAndMomentsUtilities::CalculateEmbeddedFlowForceCenter(const ModelPart& rModelPart)
    {
        // Initialize total flow force
        double tot_cut_area = 0.0;
        array_1d<double, 3> flow_force_center = ZeroVector(3);
        double& r_flow_force_center_x = flow_force_center[0];
        double& r_flow_force_center_y = flow_force_center[1];
        double& r_flow_force_center_z = flow_force_center[2];

        // Iterate the model part elements to compute the drag
        double elem_cut_area;
        array_1d<double, 3> elem_flow_force_center;

        // Auxiliary var to make the reduction
        double flow_force_x_center_red = 0.0;
        double flow_force_y_center_red = 0.0;
        double flow_force_z_center_red = 0.0;

        #pragma omp parallel for reduction(+:flow_force_x_center_red) reduction(+:flow_force_y_center_red) reduction(+:flow_force_z_center_red) reduction(+:tot_cut_area) private(elem_flow_force_center, elem_cut_area) schedule(dynamic)
        for(int i = 0; i < static_cast<int>(rModelPart.Elements().size()); ++i){
            auto it_elem = rModelPart.ElementsBegin() + i;
            it_elem->Calculate(CUTTED_AREA, elem_cut_area, rModelPart.GetProcessInfo());
            it_elem->Calculate(DRAG_FORCE_CENTER, elem_flow_force_center, rModelPart.GetProcessInfo());
            tot_cut_area += elem_cut_area;
            flow_force_x_center_red += elem_cut_area * elem_flow_force_center[0];
            flow_force_y_center_red += elem_cut_area * elem_flow_force_center[1];
            flow_force_z_center_red += elem_cut_area * elem_flow_force_center[2];
        }

        r_flow_force_center_x = flow_force_x_center_red;
        r_flow_force_center_y = flow_force_y_center_red;
        r_flow_force_center_z = flow_force_z_center_red;

        const double tol = 1.0e-12;
        if (tot_cut_area > tol) {
            flow_force_center /= tot_cut_area;
        }

        // Perform MPI synchronization
        flow_force_center = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(flow_force_center);

        return flow_force_center;
    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const FlowForcesAndMomentsUtilities& rThis) {

        rThis.PrintData(rOStream);
        return rOStream;
    }

}
