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
#include "stub_process_info_parser.h"

namespace Kratos
{

const std::string vector_parameter_string =
R"({
    "model_part_name": "PorousDomain.Body_Acceleration-auto-1",
    "variable_name":   "VOLUME_ACCELERATION",
    "active":          [false,true,false],
    "value":           [0.0,-10,0.0],
    "table":           [0,0,0]
})";

const std::string parameter_field_string =
R"({
    "model_part_name": "PorousDomain.Initial_OCR_field",
    "variable_name":   "OCR",
    "func_type":       "input",
    "function":        "1.5",
    "dataset":         "empty"
})";

const std::string excavation_string =
R"({
    "model_part_name": "PorousDomain.Excavation-auto-1",
    "variable_name": "EXCAVATION",
    "deactivate_soil_part": false
})";

const std::string k0_string =
R"({
    "model_part_name": "PorousDomain.porous_computational_model_part",
    "variable_name": "CAUCHY_STRESS_TENSOR"
})";

const std::vector<ProcessParameters> process_list =
{ProcessParameters{"ApplyVectorConstraintTableProcess", Parameters{vector_parameter_string}},
 ProcessParameters{"SetParameterFieldProcess", Parameters{parameter_field_string}},
 ProcessParameters{"ApplyExcavationProcess", Parameters{excavation_string}},
 ProcessParameters{"ApplyK0ProcedureProcess", Parameters{k0_string}}};

std::vector<ProcessParameters> StubProcessInfoParser::GetProcessList(const Kratos::Parameters& rProcessParameters) const
{
    return process_list;
}

std::size_t StubProcessInfoParser::NumberOfProcesses()
{
    return process_list.size();
}

}