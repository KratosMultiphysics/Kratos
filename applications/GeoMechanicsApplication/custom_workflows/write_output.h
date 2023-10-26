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

#pragma once

#include <string>

#include "includes/gid_io.h"
#include "includes/kratos_parameters.h"

namespace Kratos {

class Model;
class ModelPart;

// A collection of functions to write output to files. According to Kratos' coding style we should use static member
// functions rather than nonmember (i.e. free) functions. See
// https://github.com/KratosMultiphysics/Kratos/wiki/Namespaces-vs-Static-Classes for details.
class GeoOutputWriter {
public:
    GeoOutputWriter(Parameters Settings,
                    const std::string& rWorkingDirectory,
                    ModelPart& rModelPart);

    void WriteGiDOutput(ModelPart& rModelPart,
                        Parameters Settings,
                        const std::string& rWorkingDirectory,
                        bool WriteHydraulicHeadToNodes = true);

    void FinalizeResults();

private:
    void WriteNodalOutput(const std::vector<std::string>& rOutputItemNames,
                          GidIO<>& rGidIO,
                          const ModelPart& rModelPart);
    void WriteIntegrationPointOutput(const std::vector<std::string>& rOutputItemNames,
                                     GidIO<>& rGidIO,
                                     const ModelPart& rModelPart);
    void CalculateNodalHydraulicHead(GidIO<>& rGidIO, ModelPart& rModelPart);

    GiD_PostMode GetGiDPostModeFrom(const Parameters& rGiDPostFlags);
    MultiFileFlag GetMultiFileFlagFrom(const Parameters& rGiDPostFlags);
    WriteDeformedMeshFlag GetWriteDeformedMeshFlagFrom(const Parameters& rGiDPostFlags);
    WriteConditionsFlag GetWriteConditionsFlagFrom(const Parameters& rGiDPostFlags);

    GidIO<> mGidIO;
};

} // namespace Kratos
