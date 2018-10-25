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
#include "custom_processes/compute_center_of_gravity_process.h"
#include "custom_processes/total_structural_mass_process.h"
#include "structural_mechanics_application_variables.h"
namespace Kratos
{

void ComputeCenterOfGravityProcess::Execute()
{
    KRATOS_TRY

    const std::size_t domain_size = mrThisModelPart.GetProcessInfo()[DOMAIN_SIZE];
    double total_mass = 0.0;
    array_1d<double, 3> center_of_gravity;

    // Now we iterate over the elements to calculate the total mass
    auto& elements_array = mrThisModelPart.GetCommunicator().LocalMesh().Elements();

    // Making this loop omp-parallel requires locking all the geometries & nodes, which
    // is most probably not worth the effort
    for (auto& elem_i : elements_array) {
        const double elem_mass = TotalStructuralMassProcess::CalculateElementMass(elem_i, domain_size);
        total_mass += elem_mass;

        center_of_gravity += elem_i.GetGeometry().Center() * elem_mass;
    }

    // sum up across partitions
    mrThisModelPart.GetCommunicator().SumAll(total_mass);
    mrThisModelPart.GetCommunicator().SumAll(center_of_gravity);

    center_of_gravity /= total_mass;

    std::stringstream info_stream;
    info_stream << "Center of Gravity of ModelPart \"" << mrThisModelPart.Name() << "\"";

    KRATOS_INFO(info_stream.str()) << center_of_gravity << std::endl;
    KRATOS_INFO("Hint")  << "Check variable CENTER_OF_GRAVITY in the process info in "
                         << "order to access to it in any moment" << std::endl;

    mrThisModelPart.GetProcessInfo()[CENTER_OF_GRAVITY] = center_of_gravity;

    KRATOS_CATCH("")
} // class ComputeCenterOfGravityProcess
} // namespace Kratos
