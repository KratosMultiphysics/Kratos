//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Eduardo Soudah
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"

// Application includes
#include "fluid_dynamics_biomedical_application.h"
#include "fluid_dynamics_biomedical_application_variables.h"

namespace Kratos
{

KratosFluidDynamicsBiomedicalApplication::KratosFluidDynamicsBiomedicalApplication():
    KratosApplication("FluidDynamicsBiomedicalApplication")
{}

void KratosFluidDynamicsBiomedicalApplication::Register() {
    KRATOS_INFO("") << "Initializing KratosFluidDynamicsBiomedicalApplication..." << std::endl;

    // Wall Shear Stress statistics (WSS)
    KRATOS_REGISTER_VARIABLE(WSS);   // Wall Shear Stress
	KRATOS_REGISTER_VARIABLE(TAWSS); // Time Averaged WSS
    KRATOS_REGISTER_VARIABLE(TWSS);  // Time WSS
	KRATOS_REGISTER_VARIABLE(ECAP);  // Endothelial cell activation potential (OSI/TAWSS)
	KRATOS_REGISTER_VARIABLE(RRT);   //  Relative Residence Time
	KRATOS_REGISTER_VARIABLE(OSI);   //  Oscillatory Shear Index
	KRATOS_REGISTER_VARIABLE(WALL_DISTANCE); // Signed distance to the wall
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(WSS_TANGENTIAL_STRESS) //  Tangential Shear Stress
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(WSS_NORMAL_STRESS)     //  Normal Shear Stress
	KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(TEMPORAL_OSI);         // Temporal OSI
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(OUTLET_NORMAL) // Outlet area normal
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(WALL_NORMAL)   // Wall area normal

}

}  // namespace Kratos.
