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

// CoSimulation includes
#include "co_simulation_api/co_sim_io.h"

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
    return CurrentTime + 0.1;
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

void ImportGeometry(CoSim::CoSimIO& rCoSimIO, const std::string& rIdentifier="")
{
    throw std::runtime_error("not yet implemented");
}

void ExportGeometry(CoSim::CoSimIO& rCoSimIO, const std::string& rIdentifier="")
{
    throw std::runtime_error("not yet implemented");
}

void ImportMesh(CoSim::CoSimIO& rCoSimIO, MeshType& rMesh, const std::string& rIdentifier="")
{
    throw std::runtime_error("not yet implemented");
}

void ExportMesh(CoSim::CoSimIO& rCoSimIO, const MeshType& rMesh, const std::string& rIdentifier="")
{
    std::vector<double> node_coords(rMesh.size()*rMesh.size()*3, 0.0);
    // mesh has no cells, hence arguments are only dummy
    std::vector<int> connectivities;
    std::vector<int> cell_types;
    CoSim::DataContainers::Mesh mesh = {node_coords, connectivities, cell_types};
    rCoSimIO.Export(mesh, rIdentifier);
}

void ImportData(CoSim::CoSimIO& rCoSimIO, DataFieldType& rDataField, const std::string& rIdentifier="")
{
    CoSim::DataContainers::Data data = {rDataField};
    rCoSimIO.Import(data, rIdentifier);
}

void ExportData(CoSim::CoSimIO& rCoSimIO, const DataFieldType& rDataField, const std::string& rIdentifier="")
{
    CoSim::DataContainers::Data data = {rDataField};
    rCoSimIO.Export(data, rIdentifier);
}

void RunSolutionLoop(MeshType& rMesh, DataFieldType& rDataField)
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

void RunSolutionLoopWithWeakCoupling(MeshType& rMesh, DataFieldType& rDataField)
{
    // Note the following only works with one coupling inteface, requires more effort to make it work with multiple coupling interfaces.

    CoSim::CoSimIO co_sim_io("dummy_solver_cpp", "dummy_solver_io_settings");

    co_sim_io.Connect();
    ExportMesh(co_sim_io, rMesh); // send the coupling-interface mesh to be used for e.g. mapping

    for (int i=0; i<3; ++i) {
        AdvanceInTime(0.0);
        InitializeSolutionStep();
        Predict();

        ImportData(co_sim_io, rDataField);
        SolveSolutionStep();
        ExportData(co_sim_io, rDataField);

        FinalizeSolutionStep();
        OutputSolutionStep();
        std::cout << std::endl;
    }

    co_sim_io.Disconnect();
}

void RunSolutionLoopWithStrongCoupling(MeshType& rMesh, DataFieldType& rDataField)
{
    // Note the following only works with one coupling inteface, requires more effort to make it work with multiple coupling interfaces.

    CoSim::CoSimIO co_sim_io("dummy_solver_cpp", "dummy_solver_io_settings");

    co_sim_io.Connect();
    ExportMesh(co_sim_io, rMesh); // send the coupling-interface mesh to be used for e.g. mapping

    int control_signal;
    std::string identifier;

    for (int i=0; i<3; ++i) {
        AdvanceInTime(0.0);
        InitializeSolutionStep();
        Predict();
        while(true) {
            ImportData(co_sim_io, rDataField);
            SolveSolutionStep();
            ExportData(co_sim_io, rDataField);
            control_signal = co_sim_io.RecvControlSignal(identifier);
            if (control_signal == 51) { // convergence acheived
                break;
            }
        }

        FinalizeSolutionStep();
        OutputSolutionStep();
        std::cout << std::endl;
    }

    co_sim_io.Disconnect();
}

void RunSolutionCoSimulationOrchestrated(MeshType& rMesh, DataFieldType& rDataField)
{
    CoSim::CoSimIO co_sim_io("dummy_solver_cpp", "dummy_solver_io_settings");

    co_sim_io.Connect();

    int control_signal;
    std::string identifier;
    while(true) {
        // receive control signal to decide what to do
        // the signals are defined in KratosMultiphysics/applications/CoSimulationApplication/python_scripts/co_simulation_tools.py
        control_signal = co_sim_io.RecvControlSignal(identifier);
        if (control_signal == 1) {
            break; // coupled simulation is done
        }
        else if (control_signal == 11) {
            std::vector<double> time_vec(1);
            // co_sim_io.Import(/*data*/); // import current time
            const double current_time = time_vec[0];
            const double new_time = AdvanceInTime(current_time);
            time_vec[0] = new_time;
            // co_sim_io.Export(/*data*/); // export new time
        }
        else if (control_signal == 12) {
            InitializeSolutionStep();
        }
        else if (control_signal == 13) {
            SolveSolutionStep();
        }
        else if (control_signal == 14) {
            FinalizeSolutionStep();
        }
        else if (control_signal == 21) {
            ImportGeometry(co_sim_io, identifier);
        }
        else if (control_signal == 22) {
            ExportGeometry(co_sim_io, identifier);
        }
        else if (control_signal == 31) {
            ImportMesh(co_sim_io, rMesh, identifier);
        }
        else if (control_signal == 32) {
            ExportMesh(co_sim_io, rMesh, identifier);
        }
        else if (control_signal == 41) {
            ImportData(co_sim_io, rDataField, identifier);
        }
        else if (control_signal == 42) {
            ExportData(co_sim_io, rDataField, identifier);
        } else {
            throw std::runtime_error("Unknown control signal: " + std::to_string(control_signal));
        }
    }

    co_sim_io.Disconnect();
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
    std::cout << "\n    Number of nodes/dir: " << Settings[0] << std::endl;
}

} // helpers namespace



int main(int argc, char **argv)
{
    // defining the default settings
    int settings[] = {
        10, // number of nodes/dir
        0
    };

    ParseInput(argc, argv, settings);

    MeshType mesh;
    DataFieldType data_field;

    Initialize(mesh, data_field, settings[0]);

    if (settings[1] == 0) {
        std::cout << ">> Doing STANDALONE simulation <<\n" << std::endl;
        RunSolutionLoop(mesh, data_field);
    } else if (settings[1] == 1) {
        std::cout << ">> Doing COUPLED simulation (weakly coupled) <<\n" << std::endl;
        RunSolutionLoopWithWeakCoupling(mesh, data_field);
    } else if (settings[1] == 2) {
        std::cout << ">> Doing COUPLED simulation (strongly coupled) <<\n" << std::endl;
        RunSolutionLoopWithStrongCoupling(mesh, data_field);
    } else if (settings[1] == 3) {
        std::cout << ">> Doing COUPLED simulation orchestrated by CoSimulation <<\n" << std::endl;
        RunSolutionCoSimulationOrchestrated(mesh, data_field);
    }

    Finalize();

    std::cout << "finished simulation" << std::endl;

    return (0);
}

