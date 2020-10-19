//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Anoop Kodakkal, Mate Pentek
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
#include "base_moment_utilities.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
    /* Public functions *******************************************************/

    array_1d<double, 3> BaseMomentUtilities::CalculateBodyFittedBaseMoment(ModelPart& rModelPart, array_1d<double, 3> rReferencePoint) {

        array_1d<double, 3> base_moment = ZeroVector(3);
        double& moment_x = base_moment[0];
        double& moment_y = base_moment[1];
        double& moment_z = base_moment[2];

        // Auxiliary var to make the reduction
        double moment_x_red = 0.0;
        double moment_y_red = 0.0;
        double moment_z_red = 0.0;

        #pragma omp parallel for reduction(+:moment_x_red, moment_y_red, moment_z_red) schedule(dynamic)
        for(int i_node = 0; i_node < static_cast<int>(rModelPart.NumberOfNodes()); i_node++){
            auto it_node = rModelPart.NodesBegin() + i_node;
            auto drag = it_node->GetSolutionStepValue(REACTION,0);
            auto x = it_node->X() - rReferencePoint[0];
            auto y = it_node->Y() - rReferencePoint[1];
            auto z = it_node->Z() - rReferencePoint[2];
            moment_x_red +=  y * (-1) * drag[2] - z * (-1) * drag[1];
            moment_y_red +=  z * (-1) * drag[0] - x * (-1) * drag[2];
            moment_z_red +=  x * (-1) * drag[1] - y * (-1) * drag[0];
        }
        moment_x += moment_x_red;
        moment_y += moment_y_red;
        moment_z += moment_z_red;

        // // Perform MPI synchronization
        // noalias(base_moment) = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(drag_moment);
        // Perform MPI synchronization
        base_moment = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(base_moment);
        return base_moment;
    }

    array_1d<double, 3> BaseMomentUtilities::CalculateEmbeddedBaseMoment(ModelPart& rModelPart, array_1d<double, 3> rReferencePoint) {
        // Initialize total drag force
        array_1d<double, 3> base_moment = ZeroVector(3);
        double& moment_x = base_moment[0];
        double& moment_y = base_moment[1];
        double& moment_z = base_moment[2];

        // Iterate the model part elements to compute the moments
        array_1d<double, 3> elem_drag; 
        array_1d<double, 3> elem_drag_center;

        // Auxiliary var to make the reduction
        double moment_x_red = 0.0;
        double moment_y_red = 0.0;
        double moment_z_red = 0.0;

        #pragma omp parallel for reduction(+:moment_x_red, moment_y_red, moment_z_red) private(elem_drag) private(elem_drag_center) schedule(dynamic)
        for(int i = 0; i < static_cast<int>(rModelPart.Elements().size()); ++i){
            auto it_elem = rModelPart.ElementsBegin() + i;
            it_elem->Calculate(DRAG_FORCE_CENTER, elem_drag_center, rModelPart.GetProcessInfo());
            it_elem->Calculate(DRAG_FORCE, elem_drag, rModelPart.GetProcessInfo());           
            auto x = elem_drag_center[0] - rReferencePoint[0];
            if (elem_drag[0] == 0.0) {
               x = 0;
            }
            auto y = elem_drag_center[1] - rReferencePoint[1];
            if (elem_drag[1] == 0.0) {
               y = 0;
            }
            auto z = elem_drag_center[2] - rReferencePoint[2];
            if (elem_drag[2] == 0.0) {
               z = 0;
            }

            moment_x_red += y * elem_drag[2] - z * elem_drag[1];
            moment_y_red += z * elem_drag[0] - x * elem_drag[2];
            moment_z_red += x * elem_drag[1] - y * elem_drag[0];
        }

        moment_x += moment_x_red;
        moment_y += moment_y_red;
        moment_z += moment_z_red;

        // // Perform MPI synchronization
        // noalias(base_moment) = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(drag_moment);
        // Perform MPI synchronization
        base_moment = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(base_moment);
        return base_moment;    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const BaseMomentUtilities& rThis) {

        rThis.PrintData(rOStream);
        return rOStream;
    }

}
