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
This is a dummy solver implementation that shows how the integration of a pure C++ solver
into the CoSimulation framework works
*/

// System includes
#include <iostream>
#include <vector>
#include <array>
#include <stdexcept>

// CoSimulation includes
#include "co_simulation_io/co_sim_io.h"

typedef std::vector<std::vector<std::array<double, 2>>> MeshType;
typedef std::vector<double> DataFieldType;

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


void AdvanceInTime(double* pCurrentTime)
{
    std::cout << "\n\n  >>> AdvanceInTime: from: " << *pCurrentTime << " to: " << *pCurrentTime + 0.1 << std::endl;
    *pCurrentTime += 0.1;
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

void ImportGeometry(const std::string& rCommName, const std::string& rIdentifier)
{
    throw std::runtime_error("not yet implemented");
}

void ExportGeometry(const std::string& rCommName, const std::string& rIdentifier)
{
    throw std::runtime_error("not yet implemented");
}

void ImportMesh(const std::string& rCommName, const std::string& rIdentifier)
{
    throw std::runtime_error("not yet implemented");
}

void ExportMesh(const std::string& rCommName, const std::string& rIdentifier)
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
    // CoSim::DataContainers::Mesh data_mesh = {node_coords, connectivities, cell_types};
    // p_co_sim_io->Export(data_mesh, rIdentifier);
}

void ImportData(const std::string& rCommName, const std::string& rIdentifier)
{
    DataFieldType data_field;

    CoSimIO::ImportData(rCommName, rIdentifier, data_field);
}

void ExportData(const std::string& rCommName, const std::string& rIdentifier)
{
    MeshType mesh;
    DataFieldType data_field;

    Initialize(mesh, data_field, 15);

    CoSimIO::ExportData(rCommName, rIdentifier, data_field);
}

void ExportMeshToCoSim(const std::string& rCommName, const MeshType& rMesh, const std::string& rIdentifier)
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

    CoSimIO::ExportMesh(rCommName, rIdentifier, node_coords, connectivities, cell_types);
}

void ImportDataFromCoSim(const std::string& rCommName, DataFieldType& rDataField, const std::string& rIdentifier)
{
    CoSimIO::ImportData(rCommName, rIdentifier, rDataField);
}

void ExportDataToCoSim(const std::string& rCommName, DataFieldType& rDataField, const std::string& rIdentifier)
{
    CoSimIO::ExportData(rCommName, rIdentifier, rDataField);
}

void RunSolutionLoop(MeshType& rMesh, DataFieldType& rDataField)
{
    double current_time = 0.0;
    const double end_time = 0.49;
    while (current_time<end_time) {
        AdvanceInTime(&current_time);
        InitializeSolutionStep();
        SolveSolutionStep();
        FinalizeSolutionStep();
        std::cout << std::endl;
    }
}

void RunSolutionLoopWithWeakCoupling(MeshType& rMesh, DataFieldType& rDataField)
{
    // Note the following only works with one coupling interface, requires more effort to make it work with multiple coupling interfaces.

    const std::string comm_name("dummy_solver_weakly_coupled");

    CoSimIO::Connect(comm_name, "dummy_solver_io_settings");

    ExportMeshToCoSim(comm_name, rMesh, "interface");

    double current_time = 0.0;
    const double end_time = 0.49;
    while (current_time<end_time) {
        AdvanceInTime(&current_time);
        InitializeSolutionStep();

        ImportDataFromCoSim(comm_name, rDataField, "interface_temp");
        SolveSolutionStep();
        ExportDataToCoSim(comm_name, rDataField, "interface_pressure");

        FinalizeSolutionStep();
        std::cout << std::endl;
    }

    CoSimIO::Disconnect(comm_name);
}

void RunSolutionLoopWithStrongCoupling(MeshType& rMesh, DataFieldType& rDataField)
{
    // Note the following only works with one coupling interface, requires more effort to make it work with multiple coupling interfaces.

    const std::string comm_name("dummy_solver_strongly_coupled");

    CoSimIO::Connect(comm_name, "dummy_solver_io_settings");

    ExportMeshToCoSim(comm_name, rMesh, "interface");

    int convergence_signal;

    double current_time = 0.0;
    const double end_time = 0.49;
    while (current_time<end_time) {
        AdvanceInTime(&current_time);
        InitializeSolutionStep();

        ImportDataFromCoSim(comm_name, rDataField, "interface_temp");
        SolveSolutionStep();
        ExportDataToCoSim(comm_name, rDataField, "interface_pressure");

        CoSimIO::IsConverged(comm_name, &convergence_signal);
        if (convergence_signal) {break;}

        FinalizeSolutionStep();
        std::cout << std::endl;
    }

    CoSimIO::Disconnect(comm_name);
}

void RunSolutionCoSimulationOrchestrated(MeshType& rMesh, DataFieldType& rDataField)
{
    const std::string comm_name("dummy_solver_co_sim_controlled");

    CoSimIO::Connect(comm_name, "dummy_solver_io_settings");

    CoSimIO::Register(comm_name, "AdvanceInTime",          &AdvanceInTime);
    CoSimIO::Register(comm_name, "InitializeSolutionStep", &InitializeSolutionStep);
    CoSimIO::Register(comm_name, "SolveSolutionStep",      &SolveSolutionStep);
    CoSimIO::Register(comm_name, "FinalizeSolutionStep",   &FinalizeSolutionStep);

    CoSimIO::Register(comm_name, "ImportData",     &ImportData);
    CoSimIO::Register(comm_name, "ExportData",     &ExportData);
    CoSimIO::Register(comm_name, "ImportMesh",     &ImportMesh);
    CoSimIO::Register(comm_name, "ExportMesh",     &ExportMesh);
    CoSimIO::Register(comm_name, "ImportGeometry", &ImportGeometry);
    CoSimIO::Register(comm_name, "ExportGeometry", &ExportGeometry);

    CoSimIO::Run(comm_name);

    CoSimIO::Disconnect(comm_name);
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

