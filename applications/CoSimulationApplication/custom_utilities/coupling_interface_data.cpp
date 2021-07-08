//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes


// External includes


// Project includes
#include "coupling_interface_data.h"


namespace Kratos {

CouplingInterfaceData::CouplingInterfaceData(
    Parameters Settings,
    Model& rModel,
    const std::string& rName,
    const std::string& rSolverName)
    : mName(rName), mSolverName(rSolverName)
{
    Settings.ValidateAndAssignDefaults(CouplingInterfaceData::GetDefaultParameters());

    // checking name of CouplingInterfaceData
    KRATOS_ERROR_IF(rName == "") << "No \"name\" was specified!" << std::endl;
    const char disallowed_chars[] = {'.', ',', ':', ';', '>', '<', '/', '\'', '|', '*', '!', '"', ' '};
    for (const auto ch : disallowed_chars) {
        KRATOS_ERROR_IF_NOT(rName.find(ch) == std::string::npos) << "Name contains a character that is not allowed: \"" << std::string(1,ch) << "\"!" << std::endl;
    }

    // checking name of ModelPart
    mModelPartName = Settings["model_part_name"].GetString();
    KRATOS_ERROR_IF(mModelPartName == "") << "No \"model_part_name\" was specified!" << std::endl;

    mpModelpart = &rModel.GetModelPart(mModelPartName);


    // getting data location
    const std::string location = Settings["location"].GetString();
    const std::unordered_map<std::string, DataLocation> name_data_location_map {
        {"node_historical", DataLocation::NodeHistorical},
        {"node_non_historical", DataLocation::NodeNonHistorical},
        {"element", DataLocation::Element},
        {"condition", DataLocation::Condition},
        {"model_part", DataLocation::ModelPart},
    };

    auto location_iter = name_data_location_map.find(location);
    if (location_iter == name_data_location_map.end()) {
        std::vector<std::string> admissible_locations;
        admissible_locations.reserve(name_data_location_map.size());
        for (const auto& r_loc : name_data_location_map) {
            admissible_locations.push_back(r_loc.first);
        }
        KRATOS_ERROR << "\"" << location << "\" is not allowed as \"location\", only the following options are possible:\n" << admissible_locations << std::endl;
    }
    mDataLocation = location_iter->second;
}


Parameters CouplingInterfaceData::GetDefaultParameters()
{
    return Parameters(R"({
        "model_part_name" : "",
        "variable_name"   : "",
        "location"        : "node_historical",
        "dimension"       : -1
    })");
}


std::size_t CouplingInterfaceData::GetBufferSize() const
{
    if (mDataLocation == DataLocation::NodeHistorical) {
        return mpModelpart->GetBufferSize();
    } else {
        return 1;
    }
}

std::size_t CouplingInterfaceData::Size() const
{
    return 1;
}

void CouplingInterfaceData::GetData() const
{

}

void CouplingInterfaceData::SetData()
{

}



std::string CouplingInterfaceData::Info() const
{
    return "";
}

void CouplingInterfaceData::PrintInfo(std::ostream& rOStream) const
{

}

void CouplingInterfaceData::PrintData(std::ostream& rOStream) const
{

}


}  // namespace Kratos.
