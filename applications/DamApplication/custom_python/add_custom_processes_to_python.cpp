//   
//   Project Name:           
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: $
//

// External includes

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "includes/kratos_parameters.h"

// Processes
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


namespace Kratos
{
	
namespace Python
{

void  AddCustomProcessesToPython(pybind11::module& m)
{
    using namespace pybind11;

    // Fix Temperature
    class_< DamFixTemperatureConditionProcess, Process >
    (m, "DamFixTemperatureConditionProcess")
    .def(init < ModelPart&, Parameters&>());
    
    // Bofang Process
    class_< DamBofangConditionTemperatureProcess, Process >
    (m, "DamBofangConditionTemperatureProcess")
    .def(init < ModelPart&, Parameters&>());

    // Uniform Reservoir Temperature Process
    class_< DamReservoirConstantTemperatureProcess, Process >
    (m, "DamReservoirConstantTemperatureProcess")
    .def(init < ModelPart&, Parameters&>());
        
    // Hydrostatic condition
    class_< DamHydroConditionLoadProcess, Process >
    (m, "DamHydroConditionLoadProcess")
    .def(init < ModelPart&, Parameters&>());
        
    // Uplift Condition
    class_< DamUpliftConditionLoadProcess, Process >
    (m, "DamUpliftConditionLoadProcess")
    .def(init < ModelPart&, Parameters&>());
    
    // Uplift Condition for arch dams   
    class_< DamUpliftCircularConditionLoadProcess, Process >
    (m, "DamUpliftCircularConditionLoadProcess")
    .def(init < ModelPart&, Parameters&>());
   
   // Westergaard Condition (for hydrostatic + hydrodynamic pressure)     
    class_< DamWestergaardConditionLoadProcess, Process >
    (m, "DamWestergaardConditionLoadProcess")
    .def(init < ModelPart&, Parameters&>());

    // Nodal Young Modulus Process     
    class_< DamNodalYoungModulusProcess, Process >
    (m, "DamNodalYoungModulusProcess")
    .def(init < ModelPart&, Parameters&>());

    // Chemo Mechanical Aging Young Modulus Process     
    class_< DamChemoMechanicalAgingYoungProcess, Process >
    (m, "DamChemoMechanicalAgingYoungProcess")
    .def(init < ModelPart&, Parameters&>());

    // Added Mass Distribution     
    class_< DamAddedMassConditionProcess, Process >
    (m, "DamAddedMassConditionProcess")
    .def(init < ModelPart&, Parameters&>());

    //Temperature by device     
    class_< DamTemperaturebyDeviceProcess, Process >
    (m, "DamTemperaturebyDeviceProcess")
    .def(init < ModelPart&, Parameters&>());

    // Heat Flux by t_sol_air      
    class_< DamTSolAirHeatFluxProcess, Process >
    (m, "DamTSolAirHeatFluxProcess")
    .def(init < ModelPart&, Parameters&>());

    // Heat Source According Noorzai (Adiabatic Hidratation)      
    class_< DamNoorzaiHeatFluxProcess, Process >
    (m, "DamNoorzaiHeatFluxProcess")
    .def(init < ModelPart&, Parameters&>());
    
    // Heat Source according Azenha (Arrhenius formulation NonAdiabatic Hidratation)
    class_< DamAzenhaHeatFluxProcess, Process >
    (m, "DamAzenhaHeatFluxProcess")
    .def(init < ModelPart&, Parameters&>());
    }

}  // namespace Python.
} // Namespace Kratos

