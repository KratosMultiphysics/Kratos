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
#include "drag_utilities.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
    /* Public functions *******************************************************/

    array_1d<double, 3> DragUtilities::CalculateBodyFittedDrag(ModelPart& rModelPart) {
        // Sum the reactions in the model part of interest.
        // Note that the reactions are assumed to be already computed.
        VariableUtils variable_utils;
        auto drag_force = variable_utils.SumHistoricalVariable<array_1d<double,3>>(REACTION, rModelPart, 0);
        drag_force *= -1.0;

        return drag_force;
    }

    array_1d<double, 3> DragUtilities::CalculateEmbeddedDrag(ModelPart& rModelPart) {

        // Initialize total drag force
        array_1d<double, 3> drag_force = ZeroVector(3);
        double& drag_x = drag_force[0];
        double& drag_y = drag_force[1];
        double& drag_z = drag_force[2];

        // Iterate the model part elements to compute the drag
        array_1d<double, 3> elem_drag;

        // Auxiliary var to make the reduction
        double drag_x_red = 0.0;
        double drag_y_red = 0.0;
        double drag_z_red = 0.0;

        #pragma omp parallel for reduction(+:drag_x_red) reduction(+:drag_y_red) reduction(+:drag_z_red) private(elem_drag) schedule(dynamic)
        for(int i = 0; i < static_cast<int>(rModelPart.Elements().size()); ++i){
            auto it_elem = rModelPart.ElementsBegin() + i;
            it_elem->Calculate(DRAG_FORCE, elem_drag, rModelPart.GetProcessInfo());
            drag_x_red += elem_drag[0];
            drag_y_red += elem_drag[1];
            drag_z_red += elem_drag[2];
        }

        drag_x += drag_x_red;
        drag_y += drag_y_red;
        drag_z += drag_z_red;

        // Perform MPI synchronization
        drag_force = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(drag_force);

        return drag_force;
    }

    array_1d<double, 3> DragUtilities::CalculateEmbeddedDragCenter(const ModelPart& rModelPart)
    {
        // Initialize total drag force
        double tot_cut_area = 0.0;
        array_1d<double, 3> drag_force_center = ZeroVector(3);
        double& r_drag_center_x = drag_force_center[0];
        double& r_drag_center_y = drag_force_center[1];
        double& r_drag_center_z = drag_force_center[2];

        // Iterate the model part elements to compute the drag
        double elem_cut_area;
        array_1d<double, 3> elem_drag_center;

        // Auxiliary var to make the reduction
        double drag_x_center_red = 0.0;
        double drag_y_center_red = 0.0;
        double drag_z_center_red = 0.0;

        #pragma omp parallel for reduction(+:drag_x_center_red) reduction(+:drag_y_center_red) reduction(+:drag_z_center_red) reduction(+:tot_cut_area) private(elem_drag_center, elem_cut_area) schedule(dynamic)
        for(int i = 0; i < static_cast<int>(rModelPart.Elements().size()); ++i){
            auto it_elem = rModelPart.ElementsBegin() + i;
            it_elem->Calculate(CUTTED_AREA, elem_cut_area, rModelPart.GetProcessInfo());
            it_elem->Calculate(DRAG_FORCE_CENTER, elem_drag_center, rModelPart.GetProcessInfo());
            tot_cut_area += elem_cut_area;
            drag_x_center_red += elem_cut_area * elem_drag_center[0];
            drag_y_center_red += elem_cut_area * elem_drag_center[1];
            drag_z_center_red += elem_cut_area * elem_drag_center[2];
        }

        r_drag_center_x = drag_x_center_red;
        r_drag_center_y = drag_y_center_red;
        r_drag_center_z = drag_z_center_red;

        const double tol = 1.0e-12;
        if (tot_cut_area > tol) {
            drag_force_center /= tot_cut_area;
        }

        // Perform MPI synchronization
        drag_force_center = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(drag_force_center);

        return drag_force_center;
    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const DragUtilities& rThis) {

        rThis.PrintData(rOStream);
        return rOStream;
    }

}
