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

namespace Kratos
{
    /* Public functions *******************************************************/

    array_1d<double, 3> DragUtilities::CalculateBodyFittedDrag(ModelPart& rModelPart) {
        // Sum the reactions in the model part of interest.
        // Note that the reactions are assumed to be already computed.
        VariableUtils variable_utils;
        array_1d<double, 3> drag_force = variable_utils.SumHistoricalNodeVectorVariable(REACTION, rModelPart, 0);
        drag_force *= -1.0;

        return drag_force;
    }
//rishith
    double DragUtilities::CalculateBodyShearForce(ModelPart& rModelPart,double centreX, double centreY, double centreZ) {
        // Sum the reactions in the model part of interest.
        // Note that the reactions are assumed to be already computed.
        KRATOS_INFO("we reached here");
        KRATOS_INFO("we reached here");
        VariableUtils variable_utils;
        
        double drag_force = variable_utils.SumHistoricalNodeVectorVariableDotWithNormal(REACTION, rModelPart,centreX ,centreY ,centreZ , 0); // to find shear force on the structure, change of name required
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
        rModelPart.GetCommunicator().SumAll(drag_force);

        return drag_force;
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
