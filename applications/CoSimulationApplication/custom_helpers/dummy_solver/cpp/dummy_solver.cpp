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
#include "co_sim_io/co_sim_io.h"

typedef std::vector<std::vector<std::array<double, 2>>> MeshType;
typedef std::vector<double> DataFieldType;

namespace { // helpers namespace

static int sEchoLevel = 1;
static const double sDeltaTime = 0.1;
static const double sEndTime = 0.5;

static std::vector<double> sNodalCoords;
static std::vector<int> sElementConnectivities;
static std::vector<int> sElementTypes;

static std::vector<double> sFieldPressure;
static std::vector<double> sFieldVelocity;

#define solver_print(level) if(sEchoLevel>=level) std::cout << "Solver [C++]: "

void Initialize()
{
    // Defining the mesh

    /*    -- Mesh --
        0      2      3
        x------x------x
         \     |     /|\
          \  1 |  2 / | \
           \   |   /  |  \
            \  |  /   |   \
             \ | /  3 |  4 \
              \|/     |     \
               x------x-----x
               1      4     5
    */

    sNodalCoords = {
        0.0, 2.5, 1.0, // 0
        2.0, 0.0, 1.5, // 1
        2.0, 2.5, 1.5, // 2
        4.0, 2.5, 1.7, // 3
        4.0, 0.0, 1.7, // 4
        6.0, 0.0, 1.8 //  5
    };

    sElementConnectivities = {
        0, 1, 2, // 1
        1, 3, 2, // 2
        1, 4, 3, // 3
        3, 4, 5, // 4
    };

    sElementTypes = {
        5,5,5,5 // VTK_TRIANGLE
    };

    // Initializing the fields
    // sFieldPressure is coming from outside, we do nothing with it

    sFieldVelocity.resize(sNodalCoords.size()); // a value for each component on each node
    for (std::size_t i=0; i<sNodalCoords.size(); ++i) {
        sFieldVelocity[i] = i*0.125;
    }
}


void AdvanceInTime(double* pCurrentTime)
{
    *pCurrentTime += sDeltaTime;
    solver_print(2) << "AdvanceInTime; new time: " << *pCurrentTime << std::endl;
}

void InitializeSolutionStep()
{
    solver_print(2) << "InitializeSolutionStep" << std::endl;
}

void SolveSolutionStep()
{
    solver_print(2) << "SolveSolutionStep" << std::endl;
}

void FinalizeSolutionStep()
{
    solver_print(2) << "FinalizeSolutionStep" << std::endl;
}

void ImportMesh(const std::string& rCommName, const std::string& rIdentifier)
{
    // importing the mesh, but doing nothing with it to not overcomplicate this example
    std::vector<double> incoming_nodal_coords;
    std::vector<int> incoming_elem_connectivities;
    std::vector<int> incoming_elem_types;

    solver_print(3) << "\tBefore Importing Mesh from CoSim" << std::endl;
    CoSimIO::ImportMesh(rCommName, rIdentifier, incoming_nodal_coords, incoming_elem_connectivities, incoming_elem_types);
    solver_print(3) << "\tAfter Importing Mesh from CoSim" << std::endl;

    // Now copy the imported mesh to the solver-internal data structure to use it ...
    // "rIdentifier" is used to identify which mesh was imported
}

void ExportMesh(const std::string& rCommName, const std::string& rIdentifier)
{
    // usually use "rIdentifier" to select which mesh to export
    // (and then copy from solver-specific data structure to buffer which is passed to the IO)
    // in this example this is not necessary since we only export one mesh
    solver_print(3) << "\tBefore Exporting Mesh from CoSim" << std::endl;
    CoSimIO::ExportMesh(rCommName, rIdentifier, sNodalCoords, sElementConnectivities, sElementTypes);
    solver_print(3) << "\tAfter Exporting Mesh from CoSim" << std::endl;
}

void ImportDataFromCoSim(const std::string& rCommName, const std::string& rIdentifier)
{
    // usually use "rIdentifier" to select which data to import
    // (and then copy from solver-specific data structure to buffer which is passed to the IO)
    // in this example this is not necessary since we always import "field_pressure"

    solver_print(3) << "\tBefore Importing Data from CoSim" << std::endl;
    CoSimIO::ImportData(rCommName, rIdentifier, sFieldPressure);
    solver_print(3) << "\tAfter Importing Data from CoSim" << std::endl;
}

void ExportDataToCoSim(const std::string& rCommName, const std::string& rIdentifier)
{
    // usually use "rIdentifier" to select which data to export
    // (and then copy from solver-specific data structure to buffer which is passed to the IO)
    // in this example this is not necessary since we always export "field_velocity"

    solver_print(3) << "\tBefore Importing Data from CoSim" << std::endl;
    CoSimIO::ExportData(rCommName, rIdentifier, sFieldVelocity);
    solver_print(3) << "\tAfter Importing Data from CoSim" << std::endl;
}

void RunSolutionLoop()
{
    double current_time = 0.0;
    while (current_time < sEndTime) {
        AdvanceInTime(&current_time);
        InitializeSolutionStep();
        SolveSolutionStep();
        FinalizeSolutionStep();
        std::cout << std::endl;
    }
}

void RunSolutionLoopWithWeakCoupling()
{
    // Note the following only works with one coupling interface, requires more effort to make it work with multiple coupling interfaces.

    const std::string comm_name("external_dummy_solver");

    CoSimIO::Connect(comm_name, "unspecified");

    // ImportMesh(comm_name, "interface_mesh_quads");
    ExportMesh(comm_name, "interface_mesh_tri");

    double current_time = 0.0;
    while (current_time < sEndTime) {
        AdvanceInTime(&current_time);
        InitializeSolutionStep();

        ImportDataFromCoSim(comm_name, "field_pressure");
        SolveSolutionStep();
        ExportDataToCoSim(comm_name, "field_velocity");

        FinalizeSolutionStep();
        std::cout << std::endl;
    }

    CoSimIO::Disconnect(comm_name);
}

void RunSolutionLoopWithStrongCoupling()
{
    const std::string comm_name("external_dummy_solver");

    CoSimIO::Connect(comm_name, "unspecified");

    // ImportMesh(comm_name, "interface_mesh_quads");
    ExportMesh(comm_name, "interface_mesh_tri");

    int convergence_signal;

    double current_time = 0.0;
    while (current_time < sEndTime) {
        AdvanceInTime(&current_time);
        InitializeSolutionStep();

        while(true) {
            ImportDataFromCoSim(comm_name, "field_pressure");
            SolveSolutionStep();
            ExportDataToCoSim(comm_name, "field_velocity");

            CoSimIO::IsConverged(comm_name, convergence_signal);
            if (convergence_signal) {break;}
        }

        FinalizeSolutionStep();
        std::cout << std::endl;
    }

    CoSimIO::Disconnect(comm_name);
}

void RunSolutionCoSimulationOrchestrated()
{
    const std::string comm_name("external_dummy_solver");

    CoSimIO::Connect(comm_name, "unspecified");

    CoSimIO::Register(comm_name, "AdvanceInTime",          &AdvanceInTime);
    CoSimIO::Register(comm_name, "InitializeSolutionStep", &InitializeSolutionStep);
    CoSimIO::Register(comm_name, "SolveSolutionStep",      &SolveSolutionStep);
    CoSimIO::Register(comm_name, "FinalizeSolutionStep",   &FinalizeSolutionStep);

    CoSimIO::Register(comm_name, "ImportData",     &ImportDataFromCoSim);
    CoSimIO::Register(comm_name, "ExportData",     &ExportDataToCoSim);
    CoSimIO::Register(comm_name, "ImportMesh",     &ImportMesh);
    CoSimIO::Register(comm_name, "ExportMesh",     &ExportMesh);

    CoSimIO::Run(comm_name);

    CoSimIO::Disconnect(comm_name);
}

} // helpers namespace

int main(int argc, char **argv)
{
    // parsing input
    if (argc != 3) {
        throw std::runtime_error("Two input arguments required (level of coupling & echo-level)!");
    }

    const int level_of_coupling = std::atoi(argv[1]);
    sEchoLevel = std::atoi(argv[2]);

    Initialize();

    if (level_of_coupling == 0) {
        solver_print(1) << "Doing STANDALONE simulation" << std::endl;
        RunSolutionLoop();
    } else if (level_of_coupling == 1) {
        solver_print(1) << "Doing COUPLED simulation (weakly coupled)" << std::endl;
        RunSolutionLoopWithWeakCoupling();
    } else if (level_of_coupling == 2) {
        solver_print(1) << "Doing COUPLED simulation (strongly coupled)" << std::endl;
        RunSolutionLoopWithStrongCoupling();
    } else if (level_of_coupling == 3) {
        solver_print(1) << "Doing COUPLED simulation (orchestrated by CoSimulation)" << std::endl;
        RunSolutionCoSimulationOrchestrated();
    } else {
        throw std::runtime_error("ERROR, WRONG LEVEL OF COUPLING; CAN ONLY BE 0, 1, 2, 3, STOPPING");
    }

    solver_print(1) << "Exiting" << std::endl;

    return (0);
}

