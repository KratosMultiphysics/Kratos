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

#include "includes/kratos_parameters.h"
#include "includes/gid_io.h"


namespace Kratos
{

class Model;
class ModelPart;

// A collection of functions to write output to files. According to Kratos' coding style we should use static member
// functions rather than nonmember (i.e. free) functions. See
// https://github.com/KratosMultiphysics/Kratos/wiki/Namespaces-vs-Static-Classes for details.
class GeoOutputWriter
{
public:
    GeoOutputWriter() = delete;
    ~GeoOutputWriter() = delete;
    GeoOutputWriter(const GeoOutputWriter&) = delete;
    GeoOutputWriter& operator=(const GeoOutputWriter&) = delete;

    static void WriteGiDOutput(ModelPart&         rModelPart,
                               Parameters         Settings,
                               const std::string& rWorkingDirectory,
                               bool WriteHydraulicHeadToNodes = true);

private:
    static void WriteNodalOutput(const std::vector<std::string>& rOutputItemNames,
                                 GidIO<>&                        rGidIO,
                                 ModelPart&                      rModelPart);
    static void WriteIntegrationPointOutput(const std::vector<std::string>& rOutputItemNames,
                                            GidIO<>&                        rGidIO,
                                            ModelPart&                      rModelPart);
    static void CalculateNodalHydraulicHead(GidIO<>& rGidIO, ModelPart& rModelPart);
    static void PrintGaussVariable(const std::any& input, GidIO<>& rGidIO, const ModelPart& rModelPart);

    static GiD_PostMode GetGiDPostModeFrom(const Parameters& rGiDPostFlags);
    static MultiFileFlag GetMultiFileFlagFrom(const Parameters& rGiDPostFlags);
    static WriteDeformedMeshFlag GetWriteDeformedMeshFlagFrom(const Parameters& rGiDPostFlags);
    static WriteConditionsFlag GetWriteConditionsFlagFrom(const Parameters& rGiDPostFlags);

    template<class T>
    static bool PrintType(const std::any& input, GidIO<> &rGidIO, const ModelPart& rModelPart, double time)
    {
        try
        {
            const auto variable = std::any_cast<T>(input);
            rGidIO.PrintOnGaussPoints(variable, rModelPart, time, 0);
        }
        catch (const std::bad_any_cast&)
        {
            return false;
        }

        return true;
    }

};

}
