//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Sebastian Ares de Parga Regalado
//

// Application includes
#include "custom_processes/master_stiffness_process.h"

namespace Kratos
{
    /// This method is executed in order to initialize the current step
    void MasterStiffnessProcess::ExecuteInitializeSolutionStep()
    {
        KRATOS_TRY;

        this->ChangeVectorValues();
        
        KRATOS_CATCH("");
    }

    /// This method is executed in order to finalize the current step
    void MasterStiffnessProcess::ExecuteFinalizeSolutionStep()
    {
        KRATOS_TRY;
        
        this->MasterStiffnessVector();

        KRATOS_CATCH("");
    }


    ///------------------------------------------------------------------------------------

} // namespace Kratos.

