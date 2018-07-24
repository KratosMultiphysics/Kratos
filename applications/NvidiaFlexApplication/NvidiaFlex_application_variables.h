/*
 * File:   NvidiaFlex_application_variables.h
 * Author: Salva Latorre
 *
 * Created on October 9, 2014, 10:54 AM
 */

#ifndef KRATOS_NVIDIAFLEX_APPLICATION_VARIABLES_H
#define KRATOS_NVIDIAFLEX_APPLICATION_VARIABLES_H

#include "includes/define.h"
#include "includes/variables.h"
//#include "includes/nvidia_flex_variables.h"

namespace Kratos {
    
    //KRATOS_DEFINE_APPLICATION_VARIABLE(NVIDIAFLEX_APPLICATION, double, MASS_FLOW)
    //KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(NVIDIAFLEX_APPLICATION, LINEAR_VELOCITY)

    class NvidiaFlexFlags {
        
        public:
            //KRATOS_DEFINE_LOCAL_APPLICATION_FLAG(NVIDIAFLEX_APPLICATION, HAS_ROTATION);
    };
}

#endif	/* KRATOS_NVIDIAFLEX_APPLICATION_VARIABLES_H */
