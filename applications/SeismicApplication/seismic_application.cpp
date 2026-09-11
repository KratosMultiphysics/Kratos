//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Mahmoud Zidan
//


// System includes


// External includes


// Project includes
#include "includes/define.h"

#include "seismic_application.h"
#include "seismic_application_variables.h"
#include "includes/variables.h"

#include "geometries/line_gauss_lobatto_3d_2.h"


namespace Kratos {

KratosSeismicApplication::KratosSeismicApplication()
    : KratosApplication("SeismicApplication"),
    // Adding the fiber beam-column element
    mFiberBeamColumnElement3D2N(0, Element::GeometryType::Pointer(new LineGaussLobatto3D2<Node<3>>(Element::GeometryType::PointsArrayType(2))))
    {}

void KratosSeismicApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();

    KRATOS_INFO("") << "    KRATOS  ____       _               _      \n";
    KRATOS_INFO("") << "           / ___|  ___(_)___ _ __ ___ (_) ___ \n";
    KRATOS_INFO("") << "           \\___ \\ / _ \\ / __| '_ ` _ \\| |/ __|\n";
    KRATOS_INFO("") << "            ___) |  __/ \\__ \\ | | | | | | (__ \n";
    KRATOS_INFO("") << "           |____/ \\___|_|___/_| |_| |_|_|\\___|\n";
    KRATOS_INFO("") << "Initializing KratosSeismicApplication..." << std::endl;

    // Fiber Beam-Column Element Variables
    KRATOS_REGISTER_VARIABLE(ELEMENT_LOOP_TOLERANCE)
    KRATOS_REGISTER_VARIABLE(MAX_EQUILIBRIUM_ITERATIONS)
    KRATOS_REGISTER_VARIABLE(NUMBER_OF_SECTIONS)
    KRATOS_REGISTER_VARIABLE(CONCRETE_FIBERS_DATA)
    KRATOS_REGISTER_VARIABLE(STEEL_FIBERS_DATA)
    KRATOS_REGISTER_VARIABLE(BEAM_WIDTH)
    KRATOS_REGISTER_VARIABLE(BEAM_HEIGHT)

    // Fiber Beam-Column Constitutive parameters
    KRATOS_REGISTER_VARIABLE(CONCRETE_YIELD_STRENGTH)
    KRATOS_REGISTER_VARIABLE(CONCRETE_YIELD_STRAIN)
    KRATOS_REGISTER_VARIABLE(CONCRETE_CRUSHING_STRAIN)
    KRATOS_REGISTER_VARIABLE(STEEL_YOUNGS_MODULUS)
    KRATOS_REGISTER_VARIABLE(STEEL_HARDENING_RATIO)
    KRATOS_REGISTER_VARIABLE(STEEL_TRANSITION_VARIABLE)
    KRATOS_REGISTER_VARIABLE(STEEL_YIELD_STRENGTH)
    KRATOS_REGISTER_VARIABLE(STEEL_A1_COEFFICIENT)
    KRATOS_REGISTER_VARIABLE(STEEL_A2_COEFFICIENT)

    // Register the fiber beam-column element
    KRATOS_REGISTER_ELEMENT("FiberBeamColumnElement3D2N", mFiberBeamColumnElement3D2N)

}
}  // namespace Kratos.
