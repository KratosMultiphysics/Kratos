// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//					 license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher
//

/*
This is a dummy solver implementation that shows how the integration of a pure C++ client
into the CoSimulation framework works
*/

// System includes
#include <iostream>
#include <vector>
#include <array>
#include <stdexcept>

typedef std::vector<std::vector<std::array<double, 2>>> MeshType;
typedef std::vector<double> DataFieldType;

namespace { // helpers namespace

void Initialize(MeshType& rMesh, DataFieldType& rDataField, const int NumNodesPerDir)
{
    // defining (hard-coding for now) the size of the (2D) domain
    double domain_x[] = {0.0, 10,0};
    double domain_y[] = {0.0, 5,0};

    // creating the mesh
    rMesh.resize(NumNodesPerDir);

    for (int i_x=0; i_x<NumNodesPerDir; ++i_x) {
        rMesh[i_x].resize(NumNodesPerDir);
        for (int j_x=0; j_x<NumNodesPerDir; ++j_x) {

        }
    }


    // intializing the data-field
    const int size_data_field = NumNodesPerDir*NumNodesPerDir;
    rDataField.resize(size_data_field);
    for (int i=0; i<size_data_field; ++i) {
        rDataField[i] = 0.0;
    }

}


double AdvanceInTime(const double CurrentTime)
{
    std::cout << "  AdvanceInTime" << std::endl;
    return 0.0;
}
void InitializeSolutionStep()
{
    std::cout << "  InitializeSolutionStep" << std::endl;
}
void Predict()
{
    std::cout << "  Predict" << std::endl;
}
void SolveSolutionStep()
{
    std::cout << "  SolveSolutionStep" << std::endl;
}
void FinalizeSolutionStep()
{
    std::cout << "  FinalizeSolutionStep" << std::endl;
}
void OutputSolutionStep()
{
    std::cout << "  OutputSolutionStep" << std::endl;
}

void RunSolutionLoop()
{
    for (int i=0; i<3; ++i) {
        AdvanceInTime(0.0);
        InitializeSolutionStep();
        Predict();
        SolveSolutionStep();
        FinalizeSolutionStep();
        OutputSolutionStep();
        std::cout << std::endl;
    }
}

void ConductCoupling()
{
    // CoSimIO co_sim_io("dummy_solver_io_settings");

    // co_sim_io.Connect();

    // co_sim_io.Export(/*mesh*/);

    // // in this case the CoSimulation controls the solution loop / order of execution
    // while (true) {
    //     int command_id;
    //     std::string identifier;
    //     CoSimSolverWrapper::ReceiveCommand(command_id, identifier);

    //     if (command_id == 1) {
    //         std::vector<double> time_vec(1);
    //         co_sim_io.Import(/*data*/); // import current time
    //         const double current_time = time_vec[0];
    //         const double new_time = AdvanceInTime(current_time);
    //         time_vec[0] = new_time;
    //         co_sim_io.Export(/*data*/); // export new time
    //     } else if (command_id == 2) {
    //         InitializeSolutionStep();
    //     } else if (command_id == 3) {
    //         Predict();
    //     } else if (command_id == 4) {
    //         SolveSolutionStep();
    //     } else if (command_id == 5) {
    //         FinalizeSolutionStep();
    //     } else if (command_id == 11) {
    //         co_sim_io.Import(/*data*/);
    //     } else if (command_id == 12) {
    //         co_sim_io.Export(/*data*/);
    //     } else if (command_id == 6) {
    //        break;
    //     } else {
    //         error, unknown command!
    //     }
    // }

    // co_sim_io.Disconnect();
}

void Finalize()
{

}

void ParseInput(int argc, char **argv, int* Settings)
{
    if (argc > 3) {
        throw std::runtime_error("Max 2 input arguments accepted!");
    }

    for (int i=1; i<argc; ++i) {
        Settings[i-1] = std::atoi(argv[i]);
    }
    std::cout << "Using configuration:";
    std::cout << "\n    Number of nodes/dir: " << Settings[0];
    std::cout << "\n    Coupling is enabled: " << Settings[1] << std::endl << std::endl;
}

} // helpers namespace



int main(int argc, char **argv)
{
    // defining the default settings
    int settings[] = {
        10, // number of nodes/dir
        0   // enable coupling (true/false)
    };

    ParseInput(argc, argv, settings);

    MeshType mesh;
    DataFieldType data_field;

    Initialize(mesh, data_field, settings[0]);

    const bool do_coupling = settings[1];

    if (do_coupling) {
        std::cout << ">> Doing COUPLED simulation <<\n" << std::endl;
        ConductCoupling();
    } else {
        std::cout << ">> Doing STANDALONE simulation <<\n" << std::endl;
        RunSolutionLoop();
    }

    Finalize();

    return (0);
}

