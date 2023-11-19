//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes

// Project includes
#include "co_sim_io.hpp"
#include "includes/connect_impl.hpp"
#include "includes/connection.hpp"
#include "includes/data_communicator.hpp"
#include "includes/utilities.hpp"

// This file contains the implementation of the functions defined in "co_sim_io.hpp"

namespace CoSimIO {

Info Hello()
{
    std::cout << "Hello, this is the CoSimIO\n";
    std::cout << "The detached interface for coupled simulations together with the\n";
    std::cout << "CoSimulationApplication of KratosMultiphysics\n\"https://github.com/KratosMultiphysics/Kratos/tree/master/applications/CoSimulationApplication\"\n";
    std::cout << "Version:\n";
    std::cout << "    Major: " << GetMajorVersion() << "\n";
    std::cout << "    Minor: " << GetMinorVersion() << "\n";
    std::cout << "    Patch: " << GetPatchVersion() << "\n";
    std::cout << "For more information please visit \"https://github.com/KratosMultiphysics/CoSimIO\"" << std::endl;

    Info info;
    info.Set("major_version", GetMajorVersion());
    info.Set("minor_version", GetMinorVersion());
    info.Set("patch_version", GetPatchVersion());
    // TODO maybe add more things too here

    return info;
}

Info Connect(const Info& I_Settings)
{
    return Internals::ConnectImpl(
        I_Settings,
        std::make_shared<Internals::DataCommunicator>(),
        Internals::CommunicationFactory());
}


Info Disconnect(const Info& I_Info)
{
    using namespace Internals;
    const std::string connection_name = I_Info.Get<std::string>("connection_name");
    CO_SIM_IO_ERROR_IF_NOT(HasConnection(connection_name)) << "Trying to disconnect connection \"" << connection_name << "\" which does not exist!" << std::endl;

    auto info = GetConnection(connection_name).Disconnect(I_Info);
    RemoveConnection(connection_name);

    return info;
}


// Version for C++, there this input is a std::vector, which we have to wrap before passing it on
template<>
Info CO_SIM_IO_API ImportData(
    const Info& I_Info,
    std::vector<double>& rData)
{
    const std::string connection_name = I_Info.Get<std::string>("connection_name");
    using namespace CoSimIO::Internals;
    std::unique_ptr<DataContainer<double>> p_container(new DataContainerStdVector<double>(rData));
    return GetConnection(connection_name).ImportData(I_Info, *p_container);
}

// Version for C and fortran, there we already get a container
template<>
Info CO_SIM_IO_API ImportData(
    const Info& I_Info,
    CoSimIO::Internals::DataContainer<double>& rData)
{
    const std::string connection_name = I_Info.Get<std::string>("connection_name");
    return CoSimIO::Internals::GetConnection(connection_name).ImportData(I_Info, rData);
}

// Version for C++, there this input is a std::vector, which we have to wrap before passing it on
template<>
Info CO_SIM_IO_API ExportData(
    const Info& I_Info,
    const std::vector<double>& rData)
{
    const std::string connection_name = I_Info.Get<std::string>("connection_name");
    using namespace CoSimIO::Internals;
    std::unique_ptr<DataContainer<double>> p_container(new DataContainerStdVectorReadOnly<double>(rData));
    return GetConnection(connection_name).ExportData(I_Info, *p_container);
}

// Version for C and fortran, there we already get a container
template<>
Info CO_SIM_IO_API ExportData(
    const Info& I_Info,
    const CoSimIO::Internals::DataContainer<double>& rData)
{
    const std::string connection_name = I_Info.Get<std::string>("connection_name");
    return CoSimIO::Internals::GetConnection(connection_name).ExportData(I_Info, rData);
}

Info ImportMesh(
    const Info& I_Info,
    ModelPart& O_ModelPart)
{
    const std::string connection_name = I_Info.Get<std::string>("connection_name");
    return CoSimIO::Internals::GetConnection(connection_name).ImportMesh(I_Info, O_ModelPart);
}

Info ExportMesh(
    const Info& I_Info,
    const ModelPart& I_ModelPart)
{
    const std::string connection_name = I_Info.Get<std::string>("connection_name");
    return CoSimIO::Internals::GetConnection(connection_name).ExportMesh(I_Info, I_ModelPart);
}

Info ImportInfo(
    const Info& I_Info)
{
    const std::string connection_name = I_Info.Get<std::string>("connection_name");
    return CoSimIO::Internals::GetConnection(connection_name).ImportInfo(I_Info);
}

Info ExportInfo(
    const Info& I_Info)
{
    const std::string connection_name = I_Info.Get<std::string>("connection_name");
    return CoSimIO::Internals::GetConnection(connection_name).ExportInfo(I_Info);
}

Info Run(const Info& I_Info)
{
    const std::string connection_name = I_Info.Get<std::string>("connection_name");
    return CoSimIO::Internals::GetConnection(connection_name).Run(I_Info);
}


template<>
Info CO_SIM_IO_API Register(
    const Info& I_Info,
    std::function<Info(const Info&)> I_FunctionPointer)
{
    const std::string connection_name = I_Info.Get<std::string>("connection_name");
    const std::string function_name = I_Info.Get<std::string>("function_name");
    return CoSimIO::Internals::GetConnection(connection_name).Register(function_name, I_FunctionPointer);
}

template<>
Info CO_SIM_IO_API Register(
    const Info& I_Info,
    Info (*I_FunctionPointer)(const Info&))
{
    auto fct_callback = [I_FunctionPointer](const Info& I_Info)
    {
        Info info = I_FunctionPointer(I_Info);
        return info;
    };

    const std::string connection_name = I_Info.Get<std::string>("connection_name");
    const std::string function_name = I_Info.Get<std::string>("function_name");
    return CoSimIO::Internals::GetConnection(connection_name).Register(function_name, fct_callback);
}

} // namespace CoSimIO
