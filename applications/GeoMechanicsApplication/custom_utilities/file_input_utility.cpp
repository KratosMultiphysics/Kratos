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

#include "file_input_utility.h"
#include "includes/model_part_io.h"
#include "utilities/read_materials_utility.h"

namespace Kratos {

Parameters FileInputUtility::ProjectParametersFromFile(const std::filesystem::path& rProjectFilePath) const {
    std::ifstream t(rProjectFilePath);
    std::stringstream buffer;
    buffer << t.rdbuf();
    Parameters projFile{buffer.str()};
    return projFile;
}

void FileInputUtility::ReadModelFromFile(const std::filesystem::path& rModelPartFilePath,
                                         Kratos::ModelPart& rModelPart) const
{
    // Note that the file extension of the model part file must be excluded, since that is automatically appended by the
    // constructor of class ModelPartIO
    ModelPartIO reader{rModelPartFilePath.generic_string()};
    reader.ReadModelPart(rModelPart);
    KRATOS_INFO("FileInputUtility::ReadModelFromFile") << "Read the mesh data from " << rModelPartFilePath << std::endl;
}

void FileInputUtility::AddMaterialsFromFile(const std::filesystem::path& rMaterialFilePath, Model& rModel) const
{
    const std::string parameters = R"({ "Parameters" : { "materials_filename" :")" + rMaterialFilePath.generic_string() + R"("}})";
    Parameters material_file{parameters};
    ReadMaterialsUtility reader{material_file, rModel};
}

}