//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Anoop Kodakkal
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
#include "drag_and_moment_utilities.h"

namespace Kratos
{
    /* Public functions *******************************************************/

    array_1d<double, 6> DragAndMomentUtilities::CalculateBodyFittedDragAndMoment(ModelPart& rModelPart, array_1d<double, 3> rReferencePoint) {

        array_1d<double, 6> drag_force_moment = ZeroVector(6);
        array_1d<double, 6> private_drag_force_moment = ZeroVector(6);

        //#pragma omp parallel for
        for(int i_node = 0; i_node < static_cast<int>(rModelPart.NumberOfNodes()); i_node++){
            auto it_node = rModelPart.NodesBegin() + i_node;
            auto drag = it_node->GetSolutionStepValue(REACTION,0);
            auto x = it_node->X() - rReferencePoint[0];
            auto y = it_node->Y() - rReferencePoint[1];
            auto z = it_node->Z() - rReferencePoint[2];
            private_drag_force_moment[0] += -1 * drag[0];
            private_drag_force_moment[1] += -1 * drag[1];
            private_drag_force_moment[2] += -1 * drag[2];
            private_drag_force_moment[3] +=  y * (-1) * drag[2] - z * (-1) * drag[1];
            private_drag_force_moment[4] +=  z * (-1) * drag[0] - x * (-1) * drag[2];
            private_drag_force_moment[5] +=  x * (-1) * drag[1] - y * (-1) * drag[0];
        }
        drag_force_moment[0] += private_drag_force_moment[0];
        drag_force_moment[1] += private_drag_force_moment[1];
        drag_force_moment[2] += private_drag_force_moment[2];
        drag_force_moment[3] += private_drag_force_moment[3];
        drag_force_moment[4] += private_drag_force_moment[4];
        drag_force_moment[5] += private_drag_force_moment[5];

        // Perform MPI synchronization
        //drag_force_moment = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(drag_force_moment);
        return drag_force_moment;
    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const DragAndMomentUtilities& rThis) {

        rThis.PrintData(rOStream);
        return rOStream;
    }

}
