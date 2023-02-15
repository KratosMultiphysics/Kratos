// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Jonathan Nuttall
//

#pragma once

#include "dgeoapplication.h"

namespace Kratos
{
    KratosGeoApplication::KratosGeoApplication() {
        KRATOS_INFO("KratosGeoSuite") << "Setting Up Kratos" << std::endl;

        if (!kernel.IsImported("GeoMechanicsApplication"))
        {
            KRATOS_INFO("KratosGeoSuite") << "Importing GeoMechanicsApplication" << std::endl;
            geoApp = make_shared<KratosGeoMechanicsApplication>();
            kernel.ImportApplication(geoApp);
        }

        OpenMPUtils::SetNumThreads(1);
        if (this->GetEchoLevel() > 0)
        {
            OpenMPUtils::PrintOMPInfo();
        }

        this->SetEchoLevel(0);
    }

    Model& KratosGeoApplication::GetModelPointer()
    {
        return current_model;
    }

    int KratosGeoApplication::GetEchoLevel()
    {
        return echoLevel;
    }

    void KratosGeoApplication::SetEchoLevel(int level)
    {
        echoLevel = level;
    }

    void KratosGeoApplication::ResetModelParts()
    {
        KRATOS_INFO("Resetting Model") << "Setting Up Execution" << std::endl;
        current_model.Reset();
    }

    int KratosGeoApplication::execute_kratos_calculation(ModelPart& model_part, const std::vector<std::shared_ptr<Process>>& rProcesses,
        ImplicitSolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>::Pointer p_solving_strategy,
        double time, double delta_time, double number_iterations)
    {
        // Initialize
        for (auto& rProcess : rProcesses)
        {
            rProcess->ExecuteInitialize();
        }

        for (auto& rProcess : rProcesses)
        {
            rProcess->ExecuteBeforeSolutionLoop();
        }

        for (unsigned int iter = 0; iter < number_iterations; ++iter)
        {
            time += delta_time;
            model_part.CloneTimeStep(time);
            p_solving_strategy->Initialize();
            p_solving_strategy->InitializeSolutionStep();

            for (auto& rProcess : rProcesses)
            {
                rProcess->ExecuteInitializeSolutionStep();
            }

            p_solving_strategy->Predict();
            p_solving_strategy->SolveSolutionStep();

            for (auto& rProcess : rProcesses)
            {
                rProcess->ExecuteFinalizeSolutionStep();
            }

            p_solving_strategy->FinalizeSolutionStep();
        }

        for (auto process : rProcesses)
        {
            process->ExecuteFinalize();
        }

        return 0;
    }
}
