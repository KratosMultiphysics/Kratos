//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:      BSD License
//                Kratos default license: kratos/license.txt
//
//  Main authors: Riccardo Rossi
//                Raul Bravo
//


// System includes


// External includes


// Project includes
#include "rom_application.h"
#include "rom_application_variables.h"


namespace Kratos {

KratosRomApplication::KratosRomApplication():
    KratosApplication("RomApplication")
    {}

void KratosRomApplication::Register()
{
    KRATOS_INFO("") << "Initializing KratosRomApplication..." << std::endl;

    KRATOS_REGISTER_MODELER("HRomVisualizationMeshModeler", mHRomVisualizationMeshModeler);

    KRATOS_REGISTER_VARIABLE( ROM_BASIS )
    KRATOS_REGISTER_VARIABLE( ROM_LEFT_BASIS )
    KRATOS_REGISTER_VARIABLE ( HROM_WEIGHT )
    KRATOS_REGISTER_VARIABLE ( ROM_SOLUTION_INCREMENT )
    KRATOS_REGISTER_VARIABLE ( ROM_SOLUTION_BASE )
    KRATOS_REGISTER_VARIABLE ( ROM_SOLUTION_TOTAL )
    KRATOS_REGISTER_VARIABLE ( SOLUTION_BASE )
    KRATOS_REGISTER_VARIABLE ( SVD_PHI_MATRICES )
    KRATOS_REGISTER_VARIABLE ( NN_LAYERS )
    KRATOS_REGISTER_VARIABLE ( SOLUTION_REFERENCE )
    KRATOS_REGISTER_VARIABLE ( UPDATE_PHI_EFFECTIVE_BOOL )
}
}  // namespace Kratos.
