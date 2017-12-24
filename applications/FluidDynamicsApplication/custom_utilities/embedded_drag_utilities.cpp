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
#include "utilities/openmp_utils.h"

// Application includes
#include "embedded_drag_utilities.h"

namespace Kratos
{
    /* Public functions *******************************************************/

    void EmbeddedDragUtilities::CalculateDrag(
        ModelPart& rModelPart,
        array_1d<double, 3>& rDragForce) {
        
        // Initialize total drag force
        rDragForce = ZeroVector(3);

        // Initialize auxiliar arrays and partitioning
        const unsigned int num_threads = OpenMPUtils::GetNumThreads();
        std::vector<array_1d<double,3>> thread_drag_force(num_threads);

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
            rDragForce += thread_drag_force[i_thread];
        }

        // Perform MPI synchronization
        rModelPart.GetCommunicator().SumAll(rDragForce);
    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const EmbeddedDragUtilities& rThis) {

        rThis.PrintData(rOStream);
        return rOStream;
    }

}
