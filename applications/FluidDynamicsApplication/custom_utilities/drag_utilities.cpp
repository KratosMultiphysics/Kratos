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
#include "utilities/math_utils.h"

// Application includes
#include "drag_utilities.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
    DragUtilities::DragUtilities()
    {
        mComputeBodyFittedDragOnNode = [](Node<3>& rNode) {
            return -rNode.GetSolutionStepValue(REACTION, 0);
        };

        mComputeEmbeddedDragOnElement = [this](Element& rElement) {
            array_1d<double,3> elem_drag;
            rElement.Calculate(DRAG_FORCE, elem_drag, mpModelPart->GetProcessInfo());
            return elem_drag;
        };


    }

    /* Public functions *******************************************************/

    array_1d<double, 3> DragUtilities::CalculateBodyFittedDrag(ModelPart& rModelPart) {
        mpModelPart = &rModelPart;
        
        auto result_tuple = OperateAndReduceOnComponents<array_1d<double,3>>(
            rModelPart,
            rModelPart.GetCommunicator().LocalMesh().Nodes(),
            std::make_tuple(mComputeBodyFittedDragOnNode)
        );

        return std::get<0>(result_tuple);
    }

    array_1d<double, 3> DragUtilities::CalculateEmbeddedDrag(ModelPart& rModelPart) {
        mpModelPart = &rModelPart;
        
        auto result_tuple = OperateAndReduceOnComponents<array_1d<double,3>>(
            rModelPart,
            rModelPart.Elements(), 
            std::make_tuple(mComputeEmbeddedDragOnElement)
        );
        
        return std::get<0>(result_tuple);
    }

    array_1d<double, 3> DragUtilities::CalculateEmbeddedDragCenter(const ModelPart& rModelPart)
    {
        auto compute_embedded_drag_center_utils = [&rModelPart](Element& rElement)
        {
            double cut_area; 
            array_1d<double,3> drag_center;
            rElement.Calculate(CUTTED_AREA, cut_area, rModelPart.GetProcessInfo());
            rElement.Calculate(DRAG_FORCE_CENTER, drag_center, rModelPart.GetProcessInfo());

            drag_center *= cut_area;
            return array_1d<double,4> {drag_center[0], drag_center[1], drag_center[2], cut_area};
        };

        auto reduced_objects = OperateAndReduceOnComponents<array_1d<double,4>>(
            rModelPart,
            rModelPart.Elements(),
            std::make_tuple(compute_embedded_drag_center_utils)
        );

        auto& r_reduced_utils = std::get<0>(reduced_objects);

        array_1d<double,3> drag_force_center {r_reduced_utils[0], r_reduced_utils[1], r_reduced_utils[2]};
        auto& r_tot_cut_area = r_reduced_utils[3];

        // Note: shouldn't this division happen AFTER summing up with MPI?
        // (with the total cut area as denominator also synchronized with MPI)
        const double tol = 1.0e-12;
        if (r_tot_cut_area > tol) {
            drag_force_center /= r_tot_cut_area;
        }

        // Perform MPI synchronization
        drag_force_center = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(drag_force_center);

        return drag_force_center;
    }

    DragUtilities::NodalVectorFunctor DragUtilities::MakeComputeMomentOnNodeFunctor(const array_1d<double,3>& rReferencePoint) const
    {
        return [rReferencePoint](Node<3>& rNode)
        {
            auto reaction = rNode.GetSolutionStepValue(REACTION, 0);
            return MathUtils<double>::CrossProduct<array_1d<double,3>>(reaction, rNode-rReferencePoint);
        };
    }

    DragUtilities::NodalVectorFunctor DragUtilities::GetBodyFittedDragFunctor() const
    {
        return mComputeBodyFittedDragOnNode;
    }

    DragUtilities::ElementVectorFunctor DragUtilities::GetEmbeddedDragFunctor() const
    {
        return mComputeEmbeddedDragOnElement;
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
