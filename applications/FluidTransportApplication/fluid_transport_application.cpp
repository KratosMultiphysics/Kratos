//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand@{KRATOS_APP_AUTHOR}
//


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "includes/variables.h"

// Application includes
#include "fluid_transport_application.h"

namespace Kratos 
{

KratosFluidTransportApplication::KratosFluidTransportApplication():
    KratosApplication("FluidTransportApplication"),
    mSteadyConvectionDiffusionFICElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>>( Element::GeometryType::PointsArrayType (3) ) ) )
    {}

void KratosFluidTransportApplication::Register() 
{
 	// calling base class register to register Kratos components
 	KratosApplication::Register();
 	std::cout << "Initializing KratosFluidTransportApplication... " << std::endl;

    KRATOS_REGISTER_ELEMENT( "SteadyConvectionDiffusionFICElement2D3N", mSteadyConvectionDiffusionFICElement2D3N )
}

}  // namespace Kratos.
