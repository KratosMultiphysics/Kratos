// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Philipp Bucher
//                   Salman Yousaf
//

// System includes

// External includes

// Project includes
#include <math.h>
#include "custom_processes/compute_moment_of_inertia_process.h"
#include "custom_processes/total_structural_mass_process.h"
#include "structural_mechanics_application_variables.h"
namespace Kratos
{

void ComputeMomentOfInertiaProcess::Execute()
{
    KRATOS_TRY

     const std::size_t domain_size = mrThisModelPart.GetProcessInfo()[DOMAIN_SIZE];
     double moment_of_inertia = 0.0;
     
    // computing axis of rotation
     Vector3 axis_of_rotation = mrPoint2 - mrPoint1;
     const double axis_length =  std::sqrt(inner_prod(axis_of_rotation, axis_of_rotation));
     
    // Now we iterate over the elements
    auto& elements_array = mrThisModelPart.GetCommunicator().LocalMesh().Elements();
    // Making this loop omp-parallel requires locking all the geometries & nodes, which
    // is most probably not worth the effort
    for (auto& elem_i : elements_array) {
        const double elem_mass = TotalStructuralMassProcess::CalculateElementMass(elem_i, domain_size);
        Vector3 ACrossB;
        Vector3 B_vec = elem_i.GetGeometry().Center() - mrPoint1;
        MathUtils<double>::CrossProduct(ACrossB, axis_of_rotation,  B_vec);
        const double distance_from_axis = std::sqrt(inner_prod(ACrossB, ACrossB)) / axis_length ;
        moment_of_inertia += elem_mass * (distance_from_axis*distance_from_axis) ;
    }

    // sum up across partitions
    mrThisModelPart.GetCommunicator().SumAll(moment_of_inertia);

    std::stringstream info_stream;
    info_stream << "Moment of Inertia of ModelPart \"" << mrThisModelPart.Name() << "\"";

    KRATOS_INFO(info_stream.str()) << moment_of_inertia << std::endl;
    KRATOS_INFO("Hint")  << "Check variable MOMENT_OF_INERTIA in the process info in "
                         << "order to access to it at any moment" << std::endl;

    mrThisModelPart.GetProcessInfo()[MOMENT_OF_INERTIA] = moment_of_inertia;

    KRATOS_CATCH("")
} // class ComputeMomentOfInertiaProcess


} // namespace Kratos
