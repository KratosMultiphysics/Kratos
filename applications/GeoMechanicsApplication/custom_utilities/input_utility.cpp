// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "input_utility.h"
#include "includes/model_part_io.h"
#include "utilities/read_materials_utility.h"

namespace Kratos {

Parameters InputUtility::ProjectParametersFrom(const std::string &rProjectFilePath) const {
    std::ifstream t(rProjectFilePath);
    std::stringstream buffer;
    buffer << t.rdbuf();
    Parameters projFile{buffer.str()};
    return projFile;
}

void InputUtility::ReadModelFromFile(const std::filesystem::path &rModelPartFilePath,
                                            Kratos::ModelPart &rModelPart) const
{
    // Note that the file extension of the model part file must be excluded, since that is automatically appended by the
    // constructor of class ModelPartIO
    ModelPartIO reader{rModelPartFilePath.generic_string()};
    reader.ReadModelPart(rModelPart);
    KRATOS_INFO("KratosGeoSettlement") << "Read the mesh data from " << rModelPartFilePath << std::endl;
}

void InputUtility::AddMaterialsFrom(const std::string &rMaterialFilePath, Model &rModel) const
{
    const std::string parameters = R"({ "Parameters" : { "materials_filename" :")" + rMaterialFilePath + R"("}})";
    Parameters material_file{parameters};
    ReadMaterialsUtility reader{material_file, rModel};
}

}