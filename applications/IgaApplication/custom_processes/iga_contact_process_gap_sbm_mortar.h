//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Andrea Gorgi

#pragma once

// System includes
#include <vector>

// Project includes
#include "containers/model.h"
#include "geometries/geometry.h"
#include "processes/process.h"

namespace Kratos
{

class KRATOS_API(IGA_APPLICATION) IgaContactProcessGapSbmMortar
    : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(IgaContactProcessGapSbmMortar);

    using NodeType = Node;
    using GeometryType = Geometry<NodeType>;
    using GeometryPointerType = GeometryType::Pointer;
    using GeometriesArrayType = GeometryType::GeometriesArrayType;
    using CoordinatesArrayType = GeometryType::CoordinatesArrayType;

    using IndexType = std::size_t;
    using SizeType = std::size_t;

    using PropertiesPointerType = Properties::Pointer;

    IgaContactProcessGapSbmMortar(
        Model& rModel,
        Parameters ThisParameters);

    ~IgaContactProcessGapSbmMortar() override = default;

    void Execute() override;

    void ExecuteInitializeSolutionStep() override
    {
        Execute();
    }

    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters = Parameters(R"(
        {
            "analysis_model_part_name" : "ModelPart",
            "contact_sub_model_part_name" : "contact",
            "contact_condition_name" : "GapSbmContactCondition",
            "echo_level" : 0,
            "numbered_considered_neighbours" : 4,
            "shape_function_derivatives_order" : 5,
            "number_of_integration_points_per_span" : 5,
            "projection_distance_scale" : 2.0,
            "projection_distance_fallback_scale" : 0.5,
            "contact_parameters" : {
                "slave_model_part" : {
                    "sub_model_part_name" : "",
                    "layer_name" : "",
                    "property_id" : 1
                },
                "master_model_part" : {
                    "sub_model_part_name" : "",
                    "layer_name" : "",
                    "property_id" : 1
                }
            }
        })");

        return default_parameters;
    }

    std::string Info() const override
    {
        return "IgaContactProcessGapSbmMortar";
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IgaContactProcessGapSbmMortar";
    }

    void PrintData(std::ostream& rOStream) const override
    {
    }

private:
    void ClearCreatedPairingConditions();

    Model* mpModel = nullptr;
    Parameters mParameters;
    SizeType mEchoLevel = 0;

    ModelPart* mrSlaveModelPart = nullptr;
    ModelPart* mrMasterModelPart = nullptr;
    ModelPart* mpContactModelPart = nullptr;

    PropertiesPointerType mpPropMaster;
    PropertiesPointerType mpPropSlave;
    PropertiesPointerType mpContactConditionProperties;

    std::vector<IndexType> mCreatedPairingConditionIds;
};

inline std::istream& operator >> (
    std::istream& rIStream,
    IgaContactProcessGapSbmMortar& rThis);

inline std::ostream& operator << (
    std::ostream& rOStream,
    const IgaContactProcessGapSbmMortar& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos
