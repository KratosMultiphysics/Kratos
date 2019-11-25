//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "fluid_dynamics_application.h"
#include "includes/variables.h"

namespace Kratos
{

KratosFluidDynamicsBiomedicalApplication::KratosFluidDynamicsBiomedicalApplication():
    KratosApplication("FluidDynamicsBiomedicalApplication")
{}

void KratosFluidDynamicsBiomedicalApplication::Register() {
    // calling base class register to register Kratos components
    KratosApplication::Register();
    KRATOS_INFO("") << "Initializing KratosFluidDynamicsBiomedicalApplication..." << std::endl;

    // Wall Shear Stress statistics (WSS)
    KRATOS_REGISTER_VARIABLE(WSS);   // Wall Shear Stress
	KRATOS_REGISTER_VARIABLE(TAWSS); // Time Averaged WSS
    KRATOS_REGISTER_VARIABLE(TWSS);  // Time WSS
	KRATOS_REGISTER_VARIABLE(ECAP);  // Endothelial cell activation potential (OSI/TAWSS)
	KRATOS_REGISTER_VARIABLE(RRT);   //  Relative Residence Time
	KRATOS_REGISTER_VARIABLE(OSI);   //  Oscillatory Shear Index
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(WSS_TANGENTIAL_STRESS) //  Tangential Shear Stress
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(WSS_NORMAL_STRESS)     //  Normal Shear Stress
	KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(TEMPORAL_OSI);         // Temporal OSI

}

}  // namespace Kratos.
