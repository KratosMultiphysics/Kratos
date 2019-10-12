// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher
//

#ifndef KRATOS_CO_SIM_FILE_COMM_H_INCLUDED
#define KRATOS_CO_SIM_FILE_COMM_H_INCLUDED

// System includes
#include <chrono>
#include <thread>
#include <iomanip>

// Project includes
#include "co_sim_comm.h"

namespace CoSim {

namespace { // helpers namespace

static double ElapsedSeconds(const std::chrono::steady_clock::time_point& rStartTime)
{
    using namespace std::chrono;
    return duration_cast<duration<double>>(steady_clock::now() - rStartTime).count();
}

static bool FileExists(const std::string& rFileName)
{
    std::ifstream infile(rFileName);
    return infile.good(); // no need to close manually
}

static std::string GetTempFileName(const std::string& rFileName)
{
    // return std::string(rFileName).insert(CommDir.length()+1, ".");
}

static std::string GetFullPath(const std::string& rFileName)
{
    // return CommDir + "/" + rFileName; // TODO check if this work in Win
}

static void WaitForFile(const std::string& rFileName)
{
    // EMPIRE_API_LOG(1) << "Waiting for file: \"" << rFileName << "\"" << std::endl;
    while(!FileExists(rFileName)) {
        std::this_thread::sleep_for(std::chrono::milliseconds(500)); // wait 0.5s before next check
        // EMPIRE_API_LOG(3) << "    Waiting" << std::endl;
    }
    // EMPIRE_API_LOG(1) << "Found file: \"" << rFileName << "\"" << std::endl;
}

static void RemoveFile(const std::string& rFileName)
{
    if (std::remove(rFileName.c_str()) != 0) {
        // EMPIRE_API_LOG(0) << "Warning: \"" << rFileName << "\" could not be deleted!" << std::endl;
    }
}

static void MakeFileVisible(const std::string& rFinalFileName)
{
    if (std::rename(GetTempFileName(rFinalFileName).c_str(), rFinalFileName.c_str()) != 0) {
        // EMPIRE_API_LOG(0) << "Warning: \"" << rFinalFileName << "\" could not be made visible!" << std::endl;
    }
}

template <typename T>
static void CheckStream(const T& rStream, const std::string& rFileName)
{
    if (!rStream.is_open()) {
        std::stringstream err_msg;
        err_msg << rFileName << " could not be opened!";
        throw std::runtime_error(err_msg.str());
    }
}

static void SendArray(const std::string& rFileName, const int sizeOfArray, const double *data)
{
    // EMPIRE_API_LOG(2) << "Attempting to send array \"" << rFileName << "\" with size: " << sizeOfArray << " ..." << std::endl;

    const auto start_time(std::chrono::steady_clock::now());

    std::ofstream output_file;
    output_file.open(GetTempFileName(rFileName));
    CheckStream(output_file, rFileName);

    output_file << std::scientific << std::setprecision(14); // TODO maybe this should be configurable

    output_file << sizeOfArray << "\n";

    for (int i=0; i<sizeOfArray-1; ++i) {
        output_file << data[i] << " ";
    }
    output_file << data[sizeOfArray-1]; // outside to not have trailing whitespace

    output_file.close();
    MakeFileVisible(rFileName);

    // EMPIRE_API_LOG(2) << "Finished sending array" << std::endl;

    // if (PrintTiming) {
    //     // EMPIRE_API_LOG(0) << "Sending Array \"" << rFileName << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;
    // }
}

static void ReceiveArray(const std::string& rFileName, const int sizeOfArray, double *data)
{
    // EMPIRE_API_LOG(2) << "Attempting to receive array \"" << rFileName << "\" with size: " << sizeOfArray << " ..." << std::endl;

    WaitForFile(rFileName);

    const auto start_time(std::chrono::steady_clock::now());

    std::ifstream input_file(rFileName);
    CheckStream(input_file, rFileName);

    input_file >> std::setprecision(14); // TODO maybe this should be configurable

    int size_read;
    input_file >> size_read; // the first number in the file is the size of the array

    if (size_read != sizeOfArray) {
        std::stringstream err_msg;
        err_msg << "The received size for array \"" << rFileName << "\" is different from what is expected:";
        err_msg << "\n    Expected size: " << sizeOfArray;
        err_msg << "\n    Received size: " << size_read;
        throw std::runtime_error(err_msg.str());
    }

    for (int i=0; i<sizeOfArray; ++i) {
        input_file >> data[i];
    }

    RemoveFile(rFileName);

    // EMPIRE_API_LOG(2) << "Finished receiving array" << std::endl;

    // if (PrintTiming) {
    //     // EMPIRE_API_LOG(0) << "Receiving Array \"" << rFileName << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;
    // }
}

static int GetVtkCellType(const int NumberOfNodes)
{
    if (NumberOfNodes == 1) {
        return 1;
    } else if (NumberOfNodes == 2) {
        return 3;
    } else if (NumberOfNodes == 3) {
        return 5;
    } else if (NumberOfNodes == 4) {
        return 9;
    } else {
        std::stringstream err_msg;
        err_msg << "Unsupported number of nodes/element: " << NumberOfNodes;
        throw std::runtime_error(err_msg.str());
    }
}

} // helpers namespace


class FileComm : public CoSimComm
{
public:
    explicit FileComm(const std::string& rName, SettingsType& rSettings)
        : CoSimComm(rName, rSettings)
    {
        const SettingsType default_settings = {
            {"communication_folder_name_suffix", ""},
            {"use_folder_for_communication" , "0"}
        };
        Tools::AddMissingSettings(default_settings, CoSimComm::mrSettings);

        mCommFolderSuffix = CoSimComm::mrSettings.at("communication_folder_name_suffix");
        mCommInFolder = (CoSimComm::mrSettings.at("use_folder_for_communication") == "1");
    }

private:

    std::string mCommFolderSuffix = "";
    bool mCommInFolder = false;

    bool ConnectDetail() override
    {
        return true; // nothing needed here for file-based communication
    }

    bool DisconnectDetail() override
    {
        return true; // nothing needed here for file-based communication
    }

    bool ImportDetail(DataContainers::Mesh& rDataContainer, const std::string& rIdentifier) override
    {
        return true;
    }

    bool ExportDetail(const DataContainers::Mesh& rDataContainer, const std::string& rIdentifier) override
    {
        return true;
    }

    bool ImportDetail(DataContainers::Data& rDataContainer, const std::string& rIdentifier) override
    {
        return true;
    }

    bool ExportDetail(const DataContainers::Data& rDataContainer, const std::string& rIdentifier) override
    {
        return true;
    }

    bool ImportDetail(int& rDataContainer, const std::string& rIdentifier) override
    {
        return true;
    }

    bool ExportDetail(const int& rDataContainer, const std::string& rIdentifier) override
    {
        return true;
    }

};

} // namespace CoSim

#endif /* KRATOS_CO_SIM_FILE_COMM_H_INCLUDED */