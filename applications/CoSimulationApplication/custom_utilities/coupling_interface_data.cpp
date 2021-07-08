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
    mModelPartName = Settings["model_part_name"].GetString();
    mpModelpart = &rModel.GetModelPart(mModelPartName);
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
