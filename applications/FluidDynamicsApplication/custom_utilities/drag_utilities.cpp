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

// Application includes
#include "drag_utilities.h"

namespace Kratos
{
    /* Public functions *******************************************************/

    array_1d<double, 3> DragUtilities::CalculateSlipDrag(ModelPart& rModelPart) {

        // Initialize total drag force
        array_1d<double, 3> drag_force = ZeroVector(3);

        // Initialize auxiliar arrays and partitioning
        const unsigned int num_threads = OpenMPUtils::GetNumThreads();
        std::vector<array_1d<double, 3>> thread_drag_force(num_threads);

        OpenMPUtils::PartitionVector conditions_partition;
        OpenMPUtils::DivideInPartitions(rModelPart.NumberOfConditions(), num_threads, conditions_partition);

        // Iterate the model part elements to compute the drag
        #pragma omp parallel shared(thread_drag_force)
        {
            // Compute each thread drag force values
            int thread_id = OpenMPUtils::ThisThread();
            ModelPart::ConditionIterator cond_begin = rModelPart.ConditionsBegin() + conditions_partition[thread_id];
            ModelPart::ConditionIterator cond_end = rModelPart.ConditionsBegin() + conditions_partition[thread_id + 1];

            array_1d<double, 3> aux_drag_force = ZeroVector(3);

            for (ModelPart::ConditionIterator it_cond = cond_begin; it_cond != cond_end; ++it_cond) {
                // Get condition geometry
                DragUtilities::GeometryType& r_geometry = it_cond->GetGeometry();
                const unsigned int n_nodes = r_geometry.PointsNumber();

                // Get Gauss pt. data
                const Matrix N_container = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
                const unsigned int n_gauss = r_geometry.IntegrationPointsNumber(GeometryData::GI_GAUSS_2);
                const IntegrationPointsArrayType gauss_points = r_geometry.IntegrationPoints(GeometryData::GI_GAUSS_2);

                // Get condition nodal pressure values
                Vector p_values(n_nodes);
                for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
                    p_values(i_node) = r_geometry[i_node].FastGetSolutionStepValue(PRESSURE);
                }

                // Pressure drag component integration
                array_1d<double, 3> cond_drag = ZeroVector(3);

                for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
                    // Gauss pt. area normal
                    array_1d<double, 3> area_normal = r_geometry.AreaNormal(gauss_points[i_gauss].Coordinates());

                    // Gauss pt. pressure interpolation
                    double p_gauss = 0.0;
                    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
                        p_gauss += N_container(i_gauss, i_node) * p_values(i_node);
                    }

                    // Add Gauss pt. pressure normal projection contribution
                    cond_drag -= gauss_points[i_gauss].Weight() * p_gauss * area_normal;
                }

                aux_drag_force += cond_drag;
            }

            thread_drag_force[thread_id] = aux_drag_force;
        }

        // Perform reduction
        for (unsigned int i_thread = 0; i_thread < num_threads; ++i_thread)
        {
            drag_force += thread_drag_force[i_thread];
        }

        // Perform MPI synchronization
        rModelPart.GetCommunicator().SumAll(drag_force);

        return drag_force;
    }

    array_1d<double, 3> DragUtilities::CalculateEmbeddedDrag(ModelPart& rModelPart) {
        
        // Initialize total drag force
        array_1d<double, 3> drag_force = ZeroVector(3);

        // Initialize auxiliar arrays and partitioning
        const unsigned int num_threads = OpenMPUtils::GetNumThreads();
        std::vector<array_1d<double, 3>> thread_drag_force(num_threads);

        OpenMPUtils::PartitionVector elements_partition;
        OpenMPUtils::DivideInPartitions(rModelPart.NumberOfElements(), num_threads, elements_partition);

        // Iterate the model part elements to compute the drag
        #pragma omp parallel shared(thread_drag_force)
        {
            // Compute each thread drag force values
            int thread_id = OpenMPUtils::ThisThread();
            ModelPart::ElementIterator elem_begin = rModelPart.ElementsBegin() + elements_partition[thread_id];
            ModelPart::ElementIterator elem_end = rModelPart.ElementsBegin() + elements_partition[thread_id + 1];

            array_1d<double, 3> aux_drag_force = ZeroVector(3);
            for (ModelPart::ElementIterator it_elem = elem_begin; it_elem != elem_end; ++it_elem) {
                array_1d<double, 3> elem_drag;
                it_elem->Calculate(DRAG_FORCE, elem_drag, rModelPart.GetProcessInfo());
                aux_drag_force += elem_drag;
            }

            thread_drag_force[thread_id] = aux_drag_force;
        }

        // Perform reduction
        for (unsigned int i_thread = 0; i_thread < num_threads; ++i_thread) {
            drag_force += thread_drag_force[i_thread];
        }

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
