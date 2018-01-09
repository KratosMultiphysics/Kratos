//   
//   Project Name:           
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: $
//

// System includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "spaces/ublas_space.h"
#include "includes/kratos_parameters.h"

// Processes
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

using namespace boost::python;

void  AddCustomProcessesToPython() 
{    
    // Bofang Process
    class_< DamBofangConditionTemperatureProcess, bases< Process >, boost::noncopyable > ( "DamBofangConditionTemperatureProcess",
        init < ModelPart&, Parameters&>());

    // Uniform Reservoir Temperature Process
    class_< DamReservoirConstantTemperatureProcess, bases< Process >, boost::noncopyable > ( "DamReservoirConstantTemperatureProcess",
        init < ModelPart&, Parameters&>());
        
    // Hydrostatic condition
    class_< DamHydroConditionLoadProcess, bases< Process >, boost::noncopyable > ( "DamHydroConditionLoadProcess",
        init < ModelPart&, Parameters&>());
        
    // Uplift Condition
    class_< DamUpliftConditionLoadProcess, bases< Process >, boost::noncopyable > ( "DamUpliftConditionLoadProcess",
        init < ModelPart&, Parameters&>());
    
    // Uplift Condition for arch dams   
    class_< DamUpliftCircularConditionLoadProcess, bases< Process >, boost::noncopyable > ( "DamUpliftCircularConditionLoadProcess",
        init < ModelPart&, Parameters&>());
   
   // Westergaard Condition (for hydrostatic + hydrodynamic pressure)     
    class_< DamWestergaardConditionLoadProcess, bases< Process >, boost::noncopyable > ( "DamWestergaardConditionLoadProcess",
        init < ModelPart&, Parameters&>());

    // Nodal Young Modulus Process     
    class_< DamNodalYoungModulusProcess, bases< Process >, boost::noncopyable > ( "DamNodalYoungModulusProcess",
        init < ModelPart&, Parameters&>());

    // Chemo Mechanical Aging Young Modulus Process     
    class_< DamChemoMechanicalAgingYoungProcess, bases< Process >, boost::noncopyable > ( "DamChemoMechanicalAgingYoungProcess",
        init < ModelPart&, Parameters&>());

    // Added Mass Distribution     
    class_< DamAddedMassConditionProcess, bases< Process >, boost::noncopyable > ( "DamAddedMassConditionProcess",
        init < ModelPart&, Parameters&>());

    //Temperature by device     
    class_< DamTemperaturebyDeviceProcess, bases< Process >, boost::noncopyable > ( "DamTemperaturebyDeviceProcess",
        init < ModelPart&, Parameters&>());

    // Heat Flux by t_sol_air      
    class_< DamTSolAirHeatFluxProcess, bases< Process >, boost::noncopyable > ( "DamTSolAirHeatFluxProcess",
        init < ModelPart&, Parameters&>());

    // Heat Source According Noorzai (Adiabatic Hidratation)      
    class_< DamNoorzaiHeatFluxProcess, bases< Process >, boost::noncopyable > ( "DamNoorzaiHeatFluxProcess",
        init < ModelPart&, Parameters&>());
    
    // Heat Source according Azenha (Arrhenius formulation NonAdiabatic Hidratation)
    class_< DamAzenhaHeatFluxProcess, bases< Process >, boost::noncopyable > ( "DamAzenhaHeatFluxProcess",
    init < ModelPart&, Parameters&>());

    }

}  // namespace Python.
} // Namespace Kratos

