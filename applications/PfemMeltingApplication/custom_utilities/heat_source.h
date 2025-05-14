// KRATOS 
// _____   __               __  __      _ _   _             
//|  __ \ / _|             |  \/  |    | | | (_)            
//| |__) | |_ ___ _ __ ___ | \  / | ___| | |_ _ _ __   __ _ 
//|  ___/|  _/ _ \ '_ ` _ \| |\/| |/ _ \ | __| | '_ \ / _` |
//| |    | ||  __/ | | | | | |  | |  __/ | |_| | | | | (_| |
//|_|    |_| \___|_| |_| |_|_|  |_|\___|_|\__|_|_| |_|\__, |
//                                                     __/ |
//                                                    |___/ APPLICATION
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Julio Marti
//


#if !defined(KRATOS_PARTICLE_UTILITY_INCLUDED )
#define  KRATOS_PARTICLE_UTILITY_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "geometries/tetrahedra_3d_4.h"
#include "pfem_melting_application.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/timer.h"
//#include "processes/node_erase_process.h"
#include "utilities/binbased_fast_point_locator.h"
#include "includes/deprecated_variables.h"

#include <boost/timer.hpp>
#include "utilities/timer.h"

#ifdef _OPENMP
#include "omp.h"
#endif

namespace Kratos
{
  
 
  template<std::size_t TDim> class HeatSource
    {
    public:
      KRATOS_CLASS_POINTER_DEFINITION(HeatSource<TDim>);
      
      
      void Heat_Source(ModelPart & rLagrangianModelPart)
      {
        KRATOS_TRY
	  
        double density=0.0;
        double activation_energy=0.0;
        double arrhenius_coefficient=0.0;  
        double heat_of_vaporization=0.0;  
        double temperature=0.0;  
        double R=8.31; //universal gas constant
        double aux_var_polymer=0.0;

        for (ModelPart::NodesContainerType::iterator node_it = rLagrangianModelPart.NodesBegin();node_it != rLagrangianModelPart.NodesEnd(); ++node_it)
	  {
            node_it->FastGetSolutionStepValue(HEAT_FLUX) = 0.0; 

            density= node_it->FastGetSolutionStepValue(DENSITY);

            activation_energy= node_it->FastGetSolutionStepValue(ACTIVATION_ENERGY);

            arrhenius_coefficient= node_it->FastGetSolutionStepValue(ARRHENIUS_COEFFICIENT);

            heat_of_vaporization= node_it->FastGetSolutionStepValue(HEAT_OF_VAPORIZATION);

            temperature= node_it->FastGetSolutionStepValue(TEMPERATURE);

	    double E_over_R_polymer = activation_energy / R;

 	    if(temperature >800.0) temperature=800.0;

            aux_var_polymer= arrhenius_coefficient * exp(-E_over_R_polymer/temperature);


            node_it->FastGetSolutionStepValue(HEAT_FLUX) = (-1.0) * density * heat_of_vaporization * aux_var_polymer;
            

	  }
	
        KRATOS_CATCH("")
	  }
      
     


    };

} // namespace Kratos.

#endif // KRATOS_LAGRANGIAN_PARTICLES_UTILITIES_INCLUDED  defined


