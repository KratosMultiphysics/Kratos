// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//

#include "input_utilities.h"
#include "containers/model.h"
#include "utilities/read_materials_utility.h"


namespace Kratos
{

Parameters InputUtilities::ProjectParametersFrom(const std::string& rProjectFilePath)
{
    std::ifstream t(rProjectFilePath);
    std::stringstream buffer;
    buffer << t.rdbuf();
    Parameters projFile{buffer.str()};
    return projFile;
}

void InputUtilities::AddMaterialsFrom(const std::string& rMaterialFilePath,
                                      Kratos::Model&     rModel)
{
    const std::string parameters = R"({ "Parameters" : { "materials_filename" :")" + rMaterialFilePath + R"("}})";
    Parameters material_file{parameters};
    ReadMaterialsUtility reader{material_file, rModel};
}

}
