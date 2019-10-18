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

CoSim::CoSimIO* p_co_sim_io; // "hack", this will be hidden in the future!

namespace { // helpers namespace

void Initialize(MeshType& rMesh, DataFieldType& rDataField, const int NumNodesPerDir)
{
    // defining (hard-coding for now) the size of the (2D) domain
    double domain_x[] = {0.0, 10,0};
    double domain_y[] = {0.0, 5,0};

    const double dx = (domain_x[1] - domain_x[0]) / (NumNodesPerDir-1);
    const double dy = (domain_y[1] - domain_y[0]) / (NumNodesPerDir-1);

    // creating the mesh
    rMesh.resize(NumNodesPerDir);

    for (int i_x=0; i_x<NumNodesPerDir; ++i_x) {
        rMesh[i_x].resize(NumNodesPerDir);
        for (int i_y=0; i_y<NumNodesPerDir; ++i_y) {
            std::array<double, 2> xy_coord = {i_x*dx, i_y*dy};
            rMesh[i_x][i_y] = xy_coord;
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
    std::cout << "\n\n  >>> AdvanceInTime: from: " << CurrentTime << " to: " << CurrentTime + 0.1 << std::endl;
    return CurrentTime + 0.1;
}
void InitializeSolutionStep()
{
    std::cout << "  >>> InitializeSolutionStep" << std::endl;
}
// void Predict()
// {
//     std::cout << "  >>> Predict" << std::endl;
// }
void SolveSolutionStep()
{
    std::cout << "  >>> SolveSolutionStep" << std::endl;
}
void FinalizeSolutionStep()
{
    std::cout << "  >>> FinalizeSolutionStep" << std::endl;
}
// void OutputSolutionStep()
// {
//     std::cout << "  >>> OutputSolutionStep" << std::endl;
// }

void ImportGeometry(const std::string& rIdentifier)
{
    throw std::runtime_error("not yet implemented");
}

void ExportGeometry(const std::string& rIdentifier)
{
    throw std::runtime_error("not yet implemented");
}

void ImportMesh(const std::string& rIdentifier)
{
    throw std::runtime_error("not yet implemented");
}

void ExportMesh(const std::string& rIdentifier)
{
    MeshType mesh;
    DataFieldType data_field;

    Initialize(mesh, data_field, 15);

    std::vector<double> node_coords(mesh.size()*mesh.size()*3, 0.0);
    int counter=0;
    for (int i_x=0; i_x<static_cast<int>(mesh.size()); ++i_x) {
        for (int i_y=0; i_y<static_cast<int>(mesh.size()); ++i_y) {
            node_coords[counter++] = mesh[i_x][i_y][0];
            node_coords[counter++] = mesh[i_x][i_y][1];
            node_coords[counter++] = 0.0; // for 3D
        }
    }

    // mesh has no cells, hence arguments are only dummy
    std::vector<int> connectivities;
    std::vector<int> cell_types;
    CoSim::DataContainers::Mesh data_mesh = {node_coords, connectivities, cell_types};
    p_co_sim_io->Export(data_mesh, rIdentifier);
}

void ImportData(const std::string& rIdentifier)
{
    DataFieldType data_field;

    CoSim::DataContainers::Data data = {data_field};

    p_co_sim_io->Import(data, rIdentifier);
}

void ExportData(const std::string& rIdentifier)
{
    MeshType mesh;
    DataFieldType data_field;

    Initialize(mesh, data_field, 15);

    CoSim::DataContainers::Data data = {data_field};

    p_co_sim_io->Export(data, rIdentifier);
}

void ExportMesh2(CoSim::CoSimIO& rCoSimIO, const MeshType& rMesh, const std::string& rIdentifier)
{
    std::vector<double> node_coords(rMesh.size()*rMesh.size()*3, 0.0);
    int counter=0;
    for (int i_x=0; i_x<static_cast<int>(rMesh.size()); ++i_x) {
        for (int i_y=0; i_y<static_cast<int>(rMesh.size()); ++i_y) {
            node_coords[counter++] = rMesh[i_x][i_y][0];
            node_coords[counter++] = rMesh[i_x][i_y][1];
            node_coords[counter++] = 0.0; // for 3D
        }
    }

    // mesh has no cells, hence arguments are only dummy
    std::vector<int> connectivities;
    std::vector<int> cell_types;
    CoSim::DataContainers::Mesh mesh = {node_coords, connectivities, cell_types};
    rCoSimIO.Export(mesh, rIdentifier);
}

void ImportData2(CoSim::CoSimIO& rCoSimIO, DataFieldType& rDataField, const std::string& rIdentifier)
{
    CoSim::DataContainers::Data data = {rDataField};
    rCoSimIO.Import(data, rIdentifier);
}

void ExportData2(CoSim::CoSimIO& rCoSimIO, DataFieldType& rDataField, const std::string& rIdentifier)
{
    CoSim::DataContainers::Data data = {rDataField};
    rCoSimIO.Export(data, rIdentifier);
}

void RunSolutionLoop(MeshType& rMesh, DataFieldType& rDataField)
{
    double current_time = 0.0;
    const double end_time = 0.49;
    while (current_time<end_time) {
        current_time = AdvanceInTime(current_time);
        InitializeSolutionStep();
        SolveSolutionStep();
        FinalizeSolutionStep();
        std::cout << std::endl;
    }
}

void RunSolutionLoopWithWeakCoupling(MeshType& rMesh, DataFieldType& rDataField)
{
    // Note the following only works with one coupling inteface, requires more effort to make it work with multiple coupling interfaces.

    CoSim::CoSimIO co_sim_io("dummy_solver_cpp", "dummy_solver_io_settings");

    co_sim_io.Connect();
    ExportMesh2(co_sim_io, rMesh, "interface");

    double current_time = 0.0;
    const double end_time = 0.49;
    while (current_time<end_time) {
        current_time = AdvanceInTime(current_time);
        InitializeSolutionStep();

        ImportData2(co_sim_io, rDataField, "interface_temp");
        SolveSolutionStep();
        ExportData2(co_sim_io, rDataField, "interface_pressure");

        FinalizeSolutionStep();
        std::cout << std::endl;
    }

    co_sim_io.Disconnect();
}

void RunSolutionLoopWithStrongCoupling(MeshType& rMesh, DataFieldType& rDataField)
{
    // Note the following only works with one coupling inteface, requires more effort to make it work with multiple coupling interfaces.

    CoSim::CoSimIO co_sim_io("dummy_solver_cpp", "dummy_solver_io_settings");

    co_sim_io.Connect();
    ExportMesh2(co_sim_io, rMesh, "interface");

    double current_time = 0.0;
    const double end_time = 0.49;
    while (current_time<end_time) {
        current_time = AdvanceInTime(current_time);
        InitializeSolutionStep();
        while(true) {
            ImportData2(co_sim_io, rDataField, "interface_temp");
            SolveSolutionStep();
            ExportData2(co_sim_io, rDataField, "interface_pressure");
            if (co_sim_io.IsConverged()) {
                break;
            }
        }

        FinalizeSolutionStep();
        std::cout << std::endl;
    }

    co_sim_io.Disconnect();
}

void RunSolutionCoSimulationOrchestrated(MeshType& rMesh, DataFieldType& rDataField)
{
    p_co_sim_io = new CoSim::CoSimIO("dummy_solver_cpp", "dummy_solver_io_settings");

    p_co_sim_io->Connect();

    p_co_sim_io->RegisterAdvanceInTime(&AdvanceInTime);
    p_co_sim_io->RegisterInitializeSolutionStep(&InitializeSolutionStep);
    p_co_sim_io->RegisterSolveSolutionStep(&SolveSolutionStep);
    p_co_sim_io->RegisterFinalizeSolutionStep(&FinalizeSolutionStep);

    p_co_sim_io->RegisterDataExchange(&ImportGeometry, "ImportGeometry");
    p_co_sim_io->RegisterDataExchange(&ExportGeometry, "ExportGeometry");
    p_co_sim_io->RegisterDataExchange(&ImportMesh, "ImportMesh");
    p_co_sim_io->RegisterDataExchange(&ExportMesh, "ExportMesh");
    p_co_sim_io->RegisterDataExchange(&ImportData, "ImportData");
    p_co_sim_io->RegisterDataExchange(&ExportData, "ExportData");

    p_co_sim_io->Run();

    p_co_sim_io->Disconnect();

    delete p_co_sim_io;
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

