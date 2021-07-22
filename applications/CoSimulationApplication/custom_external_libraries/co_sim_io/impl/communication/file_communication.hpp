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

#ifndef CO_SIM_IO_FILE_COMMUNICATION_INCLUDED
#define CO_SIM_IO_FILE_COMMUNICATION_INCLUDED

// System includes
#include <chrono>
#include <thread>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <system_error>

// Project includes
#include "communication.hpp"
#include "../vtk_utilities.hpp"
#include "../filesystem_inc.hpp"

namespace CoSimIO {
namespace Internals {

namespace { // helpers namespace

template <typename T>
static void CheckStream(const T& rStream, const fs::path& rPath)
{
    CO_SIM_IO_ERROR_IF_NOT(rStream.is_open()) << rPath << " could not be opened!" << std::endl;
}

} // helpers namespace


class FileCommunication : public Communication
{
public:
    explicit FileCommunication(const Info& I_Settings) : Communication(I_Settings)
    {
        mCommFolder = GetWorkingDirectory();
        mCommFolder /= ".CoSimIOFileComm_" + GetConnectionName();
        mCommInFolder = I_Settings.Get<bool>("use_folder_for_communication", true);
    }

    ~FileCommunication() override
    {
        if (GetIsConnected()) {
            CO_SIM_IO_INFO("CoSimIO") << "Warning: Disconnect was not performed, attempting automatic disconnection!" << std::endl;
            Info tmp;
            Disconnect(tmp);
        }
    }

private:

    fs::path mCommFolder;
    bool mCommInFolder = true;
    mutable int mFileIndex = 0;

    Info ConnectDetail(const Info& I_Info) override
    {
        if (mCommInFolder) {
            if (GetIsPrimaryConnection()) {
                // delete and recreate directory to remove potential leftovers
                std::error_code ec;
                fs::remove_all(mCommFolder, ec);
                if (ec) {
                    CO_SIM_IO_INFO("CoSimIO") << "Warning, communication directory (" << mCommFolder << ")could not be deleted!\nError code: " << ec.message() << std::endl;
                }
                if (!fs::exists(mCommFolder)) {
                    fs::create_directory(mCommFolder);
                }
            }
        }

        ExchangeSyncFileWithPartner("connect");

        Info info;
        info.Set("is_connected", true);
        return info;
    }

    Info DisconnectDetail(const Info& I_Info) override
    {
        ExchangeSyncFileWithPartner("disconnect");

        if (mCommInFolder && GetIsPrimaryConnection()) {
            // delete directory to remove potential leftovers
            std::error_code ec;
            fs::remove_all(mCommFolder, ec);
            if (ec) {
                CO_SIM_IO_INFO("CoSimIO") << "Warning, communication directory (" << mCommFolder << ")could not be deleted!\nError code: " << ec.message() << std::endl;
            }
        }

        Info info;
        info.Set("is_connected", false);
        return info;
    }

    void ExchangeSyncFileWithPartner(const std::string& rIdentifier) const
    {
        const fs::path file_name_primary(GetFileName("CoSimIO_primary_" + rIdentifier + "_" + GetConnectionName(), "sync"));
        const fs::path file_name_secondary(GetFileName("CoSimIO_secondary_" + rIdentifier + "_" + GetConnectionName(), "sync"));

        if (GetIsPrimaryConnection()) {
            std::ofstream sync_file;
            sync_file.open(GetTempFileName(file_name_primary));
            sync_file.close();
            CO_SIM_IO_ERROR_IF_NOT(fs::exists(GetTempFileName(file_name_primary))) << "Primary sync file " << file_name_primary << " could not be created!" << std::endl;
            MakeFileVisible(file_name_primary);

            WaitForPath(file_name_secondary);
            RemovePath(file_name_secondary);

            WaitUntilFileIsRemoved(file_name_primary);
        } else {
            WaitForPath(file_name_primary);
            RemovePath(file_name_primary);

            std::ofstream sync_file;
            sync_file.open(GetTempFileName(file_name_secondary));
            sync_file.close();
            CO_SIM_IO_ERROR_IF_NOT(fs::exists(GetTempFileName(file_name_secondary))) << "Secondary sync file " << file_name_secondary << " could not be created!" << std::endl;
            MakeFileVisible(file_name_secondary);

            WaitUntilFileIsRemoved(file_name_secondary);
        }
    }

    Info ImportInfoImpl(const Info& I_Info) override
    {
        const std::string identifier = I_Info.Get<std::string>("identifier");
        CheckEntry(identifier, "identifier");

        const fs::path file_name(GetFileName("CoSimIO_info_" + GetConnectionName() + "_" + identifier, "dat"));

        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to import Info in file " << file_name << " ..." << std::endl;

        WaitForPath(file_name);

        std::ifstream input_file(file_name);
        CheckStream(input_file, file_name);

        Info imported_info;
        imported_info.Load(input_file);

        input_file.close(); // TODO check return value?
        RemovePath(file_name);

        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished importing Info" << std::endl;

        return imported_info;
    }

    Info ExportInfoImpl(const Info& I_Info) override
    {
        const std::string identifier = I_Info.Get<std::string>("identifier");
        CheckEntry(identifier, "identifier");

        const fs::path file_name(GetFileName("CoSimIO_info_" + GetConnectionName() + "_" + identifier, "dat"));

        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to export Info in file " << file_name << " ..." << std::endl;

        WaitUntilFileIsRemoved(file_name); // TODO maybe this can be queued somehow ... => then it would not block the sender

        std::ofstream output_file;
        output_file.open(GetTempFileName(file_name));
        CheckStream(output_file, file_name);

        I_Info.Save(output_file);

        output_file.close();
        MakeFileVisible(file_name);

        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished exporting Info" << std::endl;

        return Info(); // TODO use
    }

    Info ImportDataImpl(
        const Info& I_Info,
        Internals::DataContainer<double>& rData) override
    {
        const std::string identifier = I_Info.Get<std::string>("identifier");
        CheckEntry(identifier, "identifier");

        const fs::path file_name(GetFileName("CoSimIO_data_" + GetConnectionName() + "_" + identifier, "dat"));

        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to import array \"" << identifier << "\" in file " << file_name << " ..." << std::endl;

        WaitForPath(file_name);

        const auto start_time(std::chrono::steady_clock::now());

        std::ifstream input_file(file_name);
        CheckStream(input_file, file_name);

        input_file >> std::setprecision(14); // TODO maybe this should be configurable

        int size_read;
        input_file >> size_read; // the first number in the file is the size of the array

        rData.resize(size_read);

        for (int i=0; i<size_read; ++i) {
            input_file >> rData[i];
        }

        input_file.close();
        RemovePath(file_name);

        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished importing array with size: " << size_read << std::endl;

        CO_SIM_IO_INFO_IF("CoSimIO", GetPrintTiming()) << "Importing Array \"" << identifier << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;

        return Info(); // TODO use
    }

    Info ExportDataImpl(
        const Info& I_Info,
        const Internals::DataContainer<double>& rData) override
    {
        const std::string identifier = I_Info.Get<std::string>("identifier");
        CheckEntry(identifier, "identifier");

        const fs::path file_name(GetFileName("CoSimIO_data_" + GetConnectionName() + "_" + identifier, "dat"));

        WaitUntilFileIsRemoved(file_name); // TODO maybe this can be queued somehow ... => then it would not block the sender

        const std::size_t size = rData.size();
        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to export array \"" << identifier << "\" with size: " << size << " in file " << file_name << " ..." << std::endl;

        const auto start_time(std::chrono::steady_clock::now());

        std::ofstream output_file;
        output_file.open(GetTempFileName(file_name));
        CheckStream(output_file, file_name);

        output_file << std::scientific << std::setprecision(14); // TODO maybe this should be configurable

        output_file << size << "\n";

        for (std::size_t i=0; i<size-1; ++i) {
            output_file << rData[i] << " ";
        }
        // TODO check if size == 0!
        output_file << rData[size-1]; // outside to not have trailing whitespace

        output_file.close();
        MakeFileVisible(file_name);

        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished exporting array" << std::endl;

        CO_SIM_IO_INFO_IF("CoSimIO", GetPrintTiming()) << "Exporting Array \"" << identifier << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;

        return Info(); // TODO use
    }

    Info ImportMeshImpl(
        const Info& I_Info,
        ModelPart& O_ModelPart) override
    {
        const std::string identifier = I_Info.Get<std::string>("identifier");
        CheckEntry(identifier, "identifier");

        const fs::path file_name(GetFileName("CoSimIO_mesh_" + GetConnectionName() + "_" + identifier, "vtk"));

        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to import mesh \"" << identifier << "\" in file " << file_name << " ..." << std::endl;

        WaitForPath(file_name);

        const auto start_time(std::chrono::steady_clock::now());

        std::ifstream input_file(file_name);
        CheckStream(input_file, file_name);

        // reading file
        std::string current_line;
        std::vector<double> nodal_coords;
        std::vector<IdType> nodal_ids;
        std::vector<IdType> element_ids;
        std::vector<ElementType> element_types;
        std::vector<ConnectivitiesType> element_connectivities;

        while (std::getline(input_file, current_line)) {
            // reading nodes
            if (current_line.find("POINTS") != std::string::npos) {
                std::size_t num_nodes;
                current_line = current_line.substr(current_line.find("POINTS") + 7); // removing "POINTS"
                std::istringstream line_stream(current_line);
                line_stream >> num_nodes;

                CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Mesh contains " << num_nodes << " Nodes" << std::endl;

                nodal_coords.resize(3*num_nodes);
                nodal_ids.resize(num_nodes);

                for (std::size_t i=0; i<num_nodes*3; ++i) {
                    input_file >> nodal_coords[i];
                }
            }

            // reading connectivities
            if (current_line.find("CELLS") != std::string::npos) {
                std::size_t num_elems, num_nodes_per_elem;
                current_line = current_line.substr(current_line.find("CELLS") + 6); // removing "CELLS"
                std::istringstream line_stream(current_line);
                line_stream >> num_elems;

                element_ids.resize(num_elems);
                element_types.resize(num_elems);
                element_connectivities.resize(num_elems);

                CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Mesh contains " << num_elems << " Elements" << std::endl;

                for (std::size_t i=0; i<num_elems; ++i) {
                    input_file >> num_nodes_per_elem;
                    element_connectivities[i].resize(num_nodes_per_elem);
                    for (std::size_t j=0; j<num_nodes_per_elem; ++j) {
                        input_file >> element_connectivities[i][j];
                    }
                }
            }

            // reading node Ids
            if (current_line.find("NODE_ID") != std::string::npos) {
                for (std::size_t i=0; i<nodal_ids.size(); ++i) { // nodal_ids was resized to correct size above
                    input_file >> nodal_ids[i];
                }
            }

            // reading element Ids
            if (current_line.find("ELEMENT_ID") != std::string::npos) {
                for (std::size_t i=0; i<element_ids.size(); ++i) { // element_ids was resized to correct size above
                    input_file >> element_ids[i];
                }
            }

            // reading element types
            if (current_line.find("ELEMENT_TYPE") != std::string::npos) {
                int enum_temp;
                for (std::size_t i=0; i<element_types.size(); ++i) { // element_types was resized to correct size above
                    input_file >> enum_temp; // using a temp variable as enums cannot be read directly
                    element_types[i] = static_cast<CoSimIO::ElementType>(enum_temp);
                }
            }
        }

        // filling ModelPart with read information
        for (std::size_t i=0; i<nodal_ids.size(); ++i) {
            O_ModelPart.CreateNewNode(
                nodal_ids[i],
                nodal_coords[i*3],
                nodal_coords[i*3+1],
                nodal_coords[i*3+2]);
        }
        for (std::size_t i=0; i<element_ids.size(); ++i) {
            for (auto& conn : element_connectivities[i]) {
                conn = nodal_ids[conn]; // transforming vtk Ids back to original Ids
            }
            O_ModelPart.CreateNewElement(
                element_ids[i],
                element_types[i],
                element_connectivities[i]);
        }

        input_file.close();
        RemovePath(file_name);

        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished importing mesh" << std::endl;

        CO_SIM_IO_INFO_IF("CoSimIO", GetPrintTiming()) << "Importing Mesh \"" << identifier << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;

        return Info(); // TODO use
    }

    Info ExportMeshImpl(
        const Info& I_Info,
        const ModelPart& I_ModelPart) override
    {
        const std::string identifier = I_Info.Get<std::string>("identifier");
        CheckEntry(identifier, "identifier");

        const fs::path file_name(GetFileName("CoSimIO_mesh_" + GetConnectionName() + "_" + identifier, "vtk"));

        WaitUntilFileIsRemoved(file_name); // TODO maybe this can be queued somehow ... => then it would not block the sender

        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to export mesh \"" << identifier << "\" with " << I_ModelPart.NumberOfNodes() << " Nodes | " << I_ModelPart.NumberOfElements() << " Elements in file " << file_name << " ..." << std::endl;

        const auto start_time(std::chrono::steady_clock::now());

        std::ofstream output_file;
        output_file.open(GetTempFileName(file_name));
        CheckStream(output_file, file_name);

        output_file << std::scientific << std::setprecision(7); // TODO maybe this should be configurable

        // write file header
        output_file << "# vtk DataFile Version 4.0\n";
        output_file << "CoSimIO FileCommunication\n";
        output_file << "ASCII\n";
        output_file << "DATASET UNSTRUCTURED_GRID\n\n";

        // write nodes and create Id map
        std::unordered_map<IdType, IdType> id_map;
        IdType vtk_id = 0;
        output_file << "POINTS " << I_ModelPart.NumberOfNodes() << " float\n";
        for (auto node_it=I_ModelPart.NodesBegin(); node_it!=I_ModelPart.NodesEnd(); ++node_it) {
            output_file << (*node_it)->X() << " " << (*node_it)->Y() << " " << (*node_it)->Z() << "\n";
            id_map[(*node_it)->Id()] = vtk_id++;
        }
        output_file << "\n";

        // get cells size information
        std::size_t cell_list_size = 0;
        for (auto elem_it=I_ModelPart.ElementsBegin(); elem_it!=I_ModelPart.ElementsEnd(); ++elem_it) {
            cell_list_size += (*elem_it)->NumberOfNodes() + 1; // +1 for size of connectivity
        }

        // write cells connectivity
        const auto const_id_map = id_map; // const reference to not accidentially modify the map
        output_file << "CELLS " << I_ModelPart.NumberOfElements() << " " << cell_list_size << "\n";
        for (auto elem_it=I_ModelPart.ElementsBegin(); elem_it!=I_ModelPart.ElementsEnd(); ++elem_it) {
            const std::size_t num_nodes_cell = (*elem_it)->NumberOfNodes();
            output_file << num_nodes_cell << " ";
            std::size_t node_counter = 0;
            for (auto node_it=(*elem_it)->NodesBegin(); node_it!=(*elem_it)->NodesEnd(); ++node_it) {
                const IdType node_id = (*node_it)->Id();
                auto id_iter = const_id_map.find(node_id);
                CO_SIM_IO_ERROR_IF(id_iter == const_id_map.end()) << "The node with Id " << node_id << " is not part of the ModelPart but used for Element with Id " << (*elem_it)->Id() << std::endl;
                output_file << id_iter->second;
                if (node_counter++<num_nodes_cell-1) output_file << " "; // not adding a whitespace after last number
            }
            output_file << "\n";
        }
        output_file << "\n";

        // write cell types
        output_file << "CELL_TYPES " << I_ModelPart.NumberOfElements() << "\n";
        for (auto elem_it=I_ModelPart.ElementsBegin(); elem_it!=I_ModelPart.ElementsEnd(); ++elem_it) {
            output_file << static_cast<int>(GetVtkCellTypeForElementType((*elem_it)->Type())) << "\n";
        }
        output_file << "\n";

        // writing node Ids
        output_file << "POINT_DATA " << I_ModelPart.NumberOfNodes() << "\n";
        output_file << "FIELD FieldData 1" << "\n";
        output_file << "NODE_ID 1 " << I_ModelPart.NumberOfNodes() << " int\n";
        for (auto node_it=I_ModelPart.NodesBegin(); node_it!=I_ModelPart.NodesEnd(); ++node_it) {
            output_file << (*node_it)->Id() << "\n";
        }
        output_file << "\n";

        // writing element Ids
        output_file << "CELL_DATA " << I_ModelPart.NumberOfElements() << "\n";
        output_file << "FIELD FieldData 1" << "\n";
        output_file << "ELEMENT_ID 1 " << I_ModelPart.NumberOfElements() << " int\n";
        for (auto elem_it=I_ModelPart.ElementsBegin(); elem_it!=I_ModelPart.ElementsEnd(); ++elem_it) {
            output_file << (*elem_it)->Id() << "\n";
        }
        output_file << "\n";

        // writing element types
        output_file << "CELL_DATA " << I_ModelPart.NumberOfElements() << "\n";
        output_file << "FIELD FieldData 1" << "\n";
        output_file << "ELEMENT_TYPE 1 " << I_ModelPart.NumberOfElements() << " int\n";
        for (auto elem_it=I_ModelPart.ElementsBegin(); elem_it!=I_ModelPart.ElementsEnd(); ++elem_it) {
            output_file << static_cast<int>((*elem_it)->Type()) << "\n";
        }

        output_file.close();
        MakeFileVisible(file_name);

        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished exporting mesh" << std::endl;

        CO_SIM_IO_INFO_IF("CoSimIO", GetPrintTiming()) << "Exporting Mesh \"" << identifier << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;

        return Info(); // TODO use
    }

    fs::path GetTempFileName(const fs::path& rPath) const
    {
        if (mCommInFolder) {
            return rPath.string().insert(mCommFolder.string().length()+1, ".");
        } else {
            return "." + rPath.string();
        }
    }

    fs::path GetFileName(const fs::path& rPath, const std::string& rExtension) const
    {
        fs::path local_copy(rPath);
        local_copy += "_" + std::to_string((mFileIndex++)%100) + "." + rExtension;

        if (mCommInFolder) {
            return mCommFolder / local_copy;
        } else {
            return local_copy;
        }
    }

    void WaitForPath(const fs::path& rPath) const
    {
        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>0) << "Waiting for: " << rPath << std::endl;
        while(!fs::exists(rPath)) {
            std::this_thread::sleep_for(std::chrono::milliseconds(5)); // wait 0.001s before next check
        }
        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>0) << "Found: " << rPath << std::endl;
    }

    void WaitUntilFileIsRemoved(const fs::path& rPath) const
    {
        if (fs::exists(rPath)) { // only issue the wating message if the file exists initially
            CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>0) << "Waiting for: " << rPath << " to be removed" << std::endl;
            while(fs::exists(rPath)) {
                std::this_thread::sleep_for(std::chrono::milliseconds(5)); // wait 0.001s before next check
            }
            CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>0) << rPath << " was removed" << std::endl;
        }
    }

    void MakeFileVisible(const fs::path& rPath) const
    {
        std::error_code ec;
        fs::rename(GetTempFileName(rPath), rPath, ec);
        CO_SIM_IO_ERROR_IF(ec) << rPath << " could not be made visible!\nError code: " << ec.message() << std::endl;
    }

    void RemovePath(const fs::path& rPath) const
    {
        // In windows the file cannot be removed if another file handle is using it
        // this can be the case here if the partner checks if the file (still) exists
        // hence we try multiple times to delete it
        std::error_code ec;
        for (std::size_t i=0; i<5; ++i) {
            if (fs::remove(rPath, ec)) {
                return; // if file could be removed succesfully then return
            }
        }
        CO_SIM_IO_ERROR << rPath << " could not be deleted!\nError code: " << ec.message() << std::endl;
    }

};

} // namespace Internals
} // namespace CoSimIO

#endif // CO_SIM_IO_FILE_COMMUNICATION_INCLUDED
