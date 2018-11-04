//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    L. Gracia, D.J. Vicente
//
//

// External includes

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "includes/kratos_parameters.h"

// Processes
#include "custom_processes/apply_component_table_process.hpp"
#include "custom_processes/dam_fix_temperature_condition_process.hpp"
#include "custom_processes/dam_bofang_condition_temperature_process.hpp"
#include "custom_processes/dam_reservoir_constant_temperature_process.hpp"
#include "custom_processes/dam_hydro_condition_load_process.hpp"
#include "custom_processes/dam_uplift_condition_load_process.hpp"
#include "custom_processes/dam_uplift_circular_condition_load_process.hpp"
#include "custom_processes/dam_westergaard_condition_load_process.hpp"
#include "custom_processes/dam_nodal_young_modulus_process.hpp"
#include "custom_processes/dam_chemo_mechanical_aging_young_process.hpp"
#include "custom_processes/dam_temperature_by_device_process.hpp"
#include "custom_processes/dam_added_mass_condition_process.hpp"
#include "custom_processes/dam_t_sol_air_heat_flux_process.hpp"
#include "custom_processes/dam_noorzai_heat_source_process.hpp"
#include "custom_processes/dam_azenha_heat_source_process.hpp"
#include "custom_processes/dam_nodal_reference_temperature_process.hpp"
#include "custom_processes/dam_grouting_reference_temperature_process.hpp"


namespace Kratos
{

namespace Python
{

void  AddCustomProcessesToPython(pybind11::module& m)
{

    using namespace pybind11;

    typedef Table<double,double> TableType;

    // Apply table values
    class_<ApplyComponentTableProcessDam, ApplyComponentTableProcessDam::Pointer, Process>
    (m, "ApplyComponentTableProcessDam")
    .def(init < ModelPart&, Parameters>());

    // Fix Temperature
    class_<DamFixTemperatureConditionProcess, DamFixTemperatureConditionProcess::Pointer, Process>
    (m, "DamFixTemperatureConditionProcess")
    .def(init < ModelPart&, Parameters&>());

    // Bofang Process
    class_<DamBofangConditionTemperatureProcess, DamBofangConditionTemperatureProcess::Pointer, Process>
    (m, "DamBofangConditionTemperatureProcess")
    .def(init < ModelPart&, Parameters&>());

    // Uniform Reservoir Temperature Process
    class_<DamReservoirConstantTemperatureProcess, DamReservoirConstantTemperatureProcess::Pointer, Process>
    (m, "DamReservoirConstantTemperatureProcess")
    .def(init < ModelPart&, Parameters&>());

    // Hydrostatic condition
    class_<DamHydroConditionLoadProcess, DamHydroConditionLoadProcess::Pointer, Process>
    (m, "DamHydroConditionLoadProcess")
    .def(init < ModelPart&, Parameters&>());

    // Uplift Condition
    class_<DamUpliftConditionLoadProcess, DamUpliftConditionLoadProcess::Pointer, Process>
    (m, "DamUpliftConditionLoadProcess")
    .def(init < ModelPart&, Parameters&>());

    // Uplift Condition for arch dams
    class_<DamUpliftCircularConditionLoadProcess, DamUpliftCircularConditionLoadProcess::Pointer, Process>
    (m, "DamUpliftCircularConditionLoadProcess")
    .def(init < ModelPart&, Parameters&>());

   // Westergaard Condition (for hydrostatic + hydrodynamic pressure)
    class_<DamWestergaardConditionLoadProcess, DamWestergaardConditionLoadProcess::Pointer, Process>
    (m, "DamWestergaardConditionLoadProcess")
    .def(init < ModelPart&, Parameters&>());

    // Nodal Young Modulus Process
    class_<DamNodalYoungModulusProcess, DamNodalYoungModulusProcess::Pointer, Process>
    (m, "DamNodalYoungModulusProcess")
    .def(init < ModelPart&, Parameters&>());

    // Chemo Mechanical Aging Young Modulus Process
    class_<DamChemoMechanicalAgingYoungProcess, DamChemoMechanicalAgingYoungProcess::Pointer, Process>
    (m, "DamChemoMechanicalAgingYoungProcess")
    .def(init < ModelPart&, Parameters&>());

    // Added Mass Distribution
    class_<DamAddedMassConditionProcess, DamAddedMassConditionProcess::Pointer, Process>
    (m, "DamAddedMassConditionProcess")
    .def(init < ModelPart&, Parameters&>());

    //Temperature by device
    class_<DamTemperaturebyDeviceProcess, DamTemperaturebyDeviceProcess::Pointer, Process>
    (m, "DamTemperaturebyDeviceProcess")
    .def(init < ModelPart&, Parameters&>());

    // Heat Flux by t_sol_air
    class_<DamTSolAirHeatFluxProcess, DamTSolAirHeatFluxProcess::Pointer, Process>
    (m, "DamTSolAirHeatFluxProcess")
    .def(init < ModelPart&, Parameters&>());

    // Heat Source According Noorzai (Adiabatic Hidratation)
    class_<DamNoorzaiHeatFluxProcess, DamNoorzaiHeatFluxProcess::Pointer, Process>
    (m, "DamNoorzaiHeatFluxProcess")
    .def(init < ModelPart&, Parameters&>());

    // Heat Source according Azenha (Arrhenius formulation NonAdiabatic Hidratation)
    class_<DamAzenhaHeatFluxProcess, DamAzenhaHeatFluxProcess::Pointer, Process>
    (m, "DamAzenhaHeatFluxProcess")
    .def(init < ModelPart&, Parameters&>());

    // Nodal Reference Temperature Process
    class_< DamNodalReferenceTemperatureProcess, DamNodalReferenceTemperatureProcess::Pointer, Process >
    (m, "DamNodalReferenceTemperatureProcess")
    .def(init < ModelPart&, TableType&, Parameters&>());

    // Grouting Reference Temperature Process
    class_< DamGroutingReferenceTemperatureProcess, DamGroutingReferenceTemperatureProcess::Pointer, Process >
    (m, "DamGroutingReferenceTemperatureProcess")
    .def(init < ModelPart&, Parameters&>());
    }

}  // namespace Python.
} // Namespace Kratos

