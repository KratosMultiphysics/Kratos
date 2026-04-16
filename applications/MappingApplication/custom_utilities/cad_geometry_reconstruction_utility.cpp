//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan Ignacio Camarotti
//

// System includes
#include <fstream>
#include <sstream>

// Project includes
#include "custom_utilities/cad_geometry_reconstruction_utility.h"

namespace Kratos
{
namespace CadGeometryReconstructionUtility
{

namespace
{

Parameters ReadParametersFile(
    const std::string& rDataFileName,
    const int EchoLevel)
{
    const std::string data_file_name =
        (rDataFileName.size() >= 9 && rDataFileName.compare(rDataFileName.size() - 9, 9, ".cad.json") == 0)
            ? rDataFileName
            : rDataFileName + ".cad.json";

    std::ifstream infile(data_file_name);
    KRATOS_ERROR_IF_NOT(infile.good())
        << "CAD geometry file: " << data_file_name << " cannot be found." << std::endl;

    KRATOS_INFO_IF("CadGeometryReconstructionUtility", EchoLevel > 0)
        << "Reading CAD geometry file: \"" << data_file_name << "\"" << std::endl;

    std::stringstream buffer;
    buffer << infile.rdbuf();

    return Parameters(buffer.str());
}

template<class TGeometry>
void SetIdOrName(
    const Parameters rParameters,
    typename TGeometry::Pointer pGeometry)
{
    if (rParameters.Has("brep_id")) {
        pGeometry->SetId(rParameters["brep_id"].GetInt());
    } else if (rParameters.Has("brep_name")) {
        pGeometry->SetId(rParameters["brep_name"].GetString());
    }
}

std::string GetIdOrName(const Parameters rParameters)
{
    if (rParameters.Has("brep_id")) {
        return std::to_string(rParameters["brep_id"].GetInt());
    } else if (rParameters.Has("brep_name")) {
        return rParameters["brep_name"].GetString();
    } else {
        return "no_id_assigned";
    }
}

Point::Pointer ReadPoint(const Parameters rParameters)
{
    const SizeType number_of_entries = rParameters.size();

    KRATOS_ERROR_IF((number_of_entries != 1) && (number_of_entries != 2))
        << "Control points as Point need to be provided in the following structure: "
        << "[[x, y, z, weight]] or [id, [x, y, z, weight]]." << std::endl;

    const Vector cp = rParameters[number_of_entries - 1].GetVector();

    return Kratos::make_shared<Point>(cp[0], cp[1], cp[2]);
}

Vector ReadControlPointWeightVector(const Parameters rParameters)
{
    Vector control_point_weights = ZeroVector(rParameters.size());

    KRATOS_ERROR_IF(rParameters.size() == 0)
        << "Length of control point list is zero!" << std::endl;

    const SizeType number_of_entries = rParameters[0].size();

    KRATOS_ERROR_IF(rParameters[0][number_of_entries - 1].size() != 4)
        << "Control points need to be provided in the following structure: "
        << "[[x, y, z, weight]] or [id, [x, y, z, weight]]." << std::endl;

    for (IndexType cp_idx = 0; cp_idx < rParameters.size(); ++cp_idx) {
        control_point_weights[cp_idx] = rParameters[cp_idx][number_of_entries - 1][3].GetDouble();
    }

    return control_point_weights;
}

void ReadControlPointVectorFromExistingNodes(
    ContainerNodeType& rControlPoints,
    const Parameters rParameters,
    ModelPart& rModelPart,
    const int EchoLevel)
{
    KRATOS_ERROR_IF_NOT(rParameters.IsArray())
        << "\"control_points\" section needs to be an array." << std::endl;

    KRATOS_INFO_IF("CadGeometryReconstructionUtility", EchoLevel > 1)
        << "Using already existing nodes as control points. CAD cp count = "
        << rParameters.size() << ", model part node count = "
        << rModelPart.NumberOfNodes() << std::endl;

    KRATOS_ERROR_IF(rParameters.size() != rModelPart.NumberOfNodes())
        << "Mismatch between CAD control point count (" << rParameters.size()
        << ") and existing model part node count (" << rModelPart.NumberOfNodes() << ")." << std::endl;

    rControlPoints.reserve(rModelPart.NumberOfNodes());

    for (auto it_node = rModelPart.NodesBegin(); it_node != rModelPart.NodesEnd(); ++it_node) {
        rControlPoints.push_back(*(it_node.base()));
    }
}

NurbsSurfacePointerType ReadNurbsSurfaceWithExistingNodes(
    const Parameters rParameters,
    ModelPart& rModelPart,
    const int EchoLevel)
{
    bool is_rational = true;
    if (rParameters.Has("is_rational")) {
        is_rational = rParameters["is_rational"].GetBool();
    }

    KRATOS_ERROR_IF_NOT(rParameters.Has("knot_vectors"))
        << "Missing 'knot_vectors' in nurbs surface." << std::endl;
    KRATOS_ERROR_IF(rParameters["knot_vectors"].size() != 2)
        << "'knot_vectors' need to be of size two." << std::endl;

    const Vector knot_vector_u = rParameters["knot_vectors"][0].GetVector();
    const Vector knot_vector_v = rParameters["knot_vectors"][1].GetVector();

    KRATOS_ERROR_IF_NOT(rParameters.Has("degrees"))
        << "Missing 'degrees' in nurbs surface." << std::endl;
    KRATOS_ERROR_IF(rParameters["degrees"].size() != 2)
        << "'degrees' need to be of size two." << std::endl;

    const int p = rParameters["degrees"][0].GetInt();
    const int q = rParameters["degrees"][1].GetInt();

    KRATOS_ERROR_IF_NOT(rParameters.Has("control_points"))
        << "Missing 'control_points' in nurbs surface." << std::endl;

    ContainerNodeType control_points;
    ReadControlPointVectorFromExistingNodes(
        control_points,
        rParameters["control_points"],
        rModelPart,
        EchoLevel);

    if (is_rational) {
        const Vector control_point_weights = ReadControlPointWeightVector(rParameters["control_points"]);
        return Kratos::make_shared<NurbsSurfaceType>(
            control_points,
            p,
            q,
            knot_vector_u,
            knot_vector_v,
            control_point_weights);
    }

    return Kratos::make_shared<NurbsSurfaceType>(
        control_points,
        p,
        q,
        knot_vector_u,
        knot_vector_v);
}

typename BrepCurveOnSurfaceType::Pointer ReadTrimmingCurve(
    const Parameters rParameters,
    NurbsSurfacePointerType pNurbsSurface,
    ModelPart& rModelPart,
    const int EchoLevel)
{
    KRATOS_ERROR_IF_NOT(rParameters.Has("curve_direction"))
        << "Missing 'curve_direction' in trimming curve." << std::endl;
    const bool curve_direction = rParameters["curve_direction"].GetBool();

    KRATOS_ERROR_IF_NOT(rParameters.Has("parameter_curve"))
        << "Missing 'parameter_curve' in trimming curve." << std::endl;

    const auto parameter_curve = rParameters["parameter_curve"];

    const bool is_rational = parameter_curve.Has("is_rational")
        ? parameter_curve["is_rational"].GetBool()
        : true;

    KRATOS_ERROR_IF_NOT(parameter_curve.Has("knot_vector"))
        << "Missing 'knot_vector' in parameter_curve." << std::endl;
    const Vector knot_vector = parameter_curve["knot_vector"].GetVector();

    KRATOS_ERROR_IF_NOT(parameter_curve.Has("degree"))
        << "Missing 'degree' in parameter_curve." << std::endl;
    const int curve_degree = parameter_curve["degree"].GetInt();

    KRATOS_ERROR_IF_NOT(parameter_curve.Has("control_points"))
        << "Missing 'control_points' in parameter_curve." << std::endl;

    ContainerEmbeddedNodeType curve_points;
    const auto control_points = parameter_curve["control_points"];
    for (IndexType cp_idx = 0; cp_idx < control_points.size(); ++cp_idx) {
        curve_points.push_back(ReadPoint(control_points[cp_idx]));
    }

    NurbsTrimmingCurvePointerType p_trimming_curve;
    if (is_rational) {
        const Vector curve_weights = ReadControlPointWeightVector(parameter_curve["control_points"]);
        p_trimming_curve = Kratos::make_shared<NurbsTrimmingCurveType>(
            curve_points,
            curve_degree,
            knot_vector,
            curve_weights);
    } else {
        p_trimming_curve = Kratos::make_shared<NurbsTrimmingCurveType>(
            curve_points,
            curve_degree,
            knot_vector);
    }

    KRATOS_ERROR_IF_NOT(parameter_curve.Has("active_range"))
        << "Missing 'active_range' in parameter_curve." << std::endl;

    const Vector active_range_vector = parameter_curve["active_range"].GetVector();
    KRATOS_ERROR_IF(active_range_vector.size() != 2)
        << "'active_range' must have size 2." << std::endl;

    const NurbsInterval brep_active_range(active_range_vector[0], active_range_vector[1]);

    auto p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(
        pNurbsSurface,
        p_trimming_curve,
        brep_active_range,
        curve_direction);

    if (rParameters.Has("trim_index")) {
        p_brep_curve_on_surface->SetId(rParameters["trim_index"].GetInt());
    }

    return p_brep_curve_on_surface;
}

BrepCurveOnSurfaceLoopType ReadTrimmingCurveVector(
    const Parameters rParameters,
    NurbsSurfacePointerType pNurbsSurface,
    ModelPart& rModelPart,
    const int EchoLevel)
{
    KRATOS_ERROR_IF(rParameters.size() < 1)
        << "Trimming curve list has no element." << std::endl;

    BrepCurveOnSurfaceLoopType trimming_brep_curve_vector(rParameters.size());

    for (IndexType tc_idx = 0; tc_idx < rParameters.size(); ++tc_idx) {
        trimming_brep_curve_vector[tc_idx] = ReadTrimmingCurve(
            rParameters[tc_idx],
            pNurbsSurface,
            rModelPart,
            EchoLevel);
    }

    return trimming_brep_curve_vector;
}

std::tuple<BrepCurveOnSurfaceLoopArrayType, BrepCurveOnSurfaceLoopArrayType> ReadBoundaryLoops(
    const Parameters rParameters,
    NurbsSurfacePointerType pNurbsSurface,
    ModelPart& rModelPart,
    const int EchoLevel)
{
    BrepCurveOnSurfaceLoopArrayType outer_loops;
    BrepCurveOnSurfaceLoopArrayType inner_loops;

    for (IndexType bl_idx = 0; bl_idx < rParameters.size(); ++bl_idx)
    {
        KRATOS_ERROR_IF_NOT(rParameters[bl_idx].Has("loop_type"))
            << "Missing 'loop_type' in boundary loop " << bl_idx << "." << std::endl;

        const std::string loop_type = rParameters[bl_idx]["loop_type"].GetString();

        KRATOS_ERROR_IF_NOT(rParameters[bl_idx].Has("trimming_curves"))
            << "Missing 'trimming_curves' in boundary loop " << bl_idx << "." << std::endl;

        auto trimming_curves = ReadTrimmingCurveVector(
            rParameters[bl_idx]["trimming_curves"],
            pNurbsSurface,
            rModelPart,
            EchoLevel);

        if (loop_type == "outer") {
            outer_loops.resize(outer_loops.size() + 1, false);
            outer_loops[outer_loops.size() - 1] = trimming_curves;
        } else if (loop_type == "inner") {
            inner_loops.resize(inner_loops.size() + 1, false);
            inner_loops[inner_loops.size() - 1] = trimming_curves;
        } else {
            KRATOS_ERROR << "Loop type \"" << loop_type << "\" is not supported." << std::endl;
        }
    }

    return std::make_tuple(outer_loops, inner_loops);
}

void ReadAndAddEmbeddedEdges(
    typename BrepSurfaceType::Pointer pBrepSurface,
    const Parameters rParameters,
    NurbsSurfacePointerType pNurbsSurface,
    ModelPart& rModelPart,
    const int EchoLevel)
{
    if (rParameters.Has("embedded_edges")) {
        if (rParameters["embedded_edges"].size() > 0) {
            BrepCurveOnSurfaceArrayType embedded_edges = ReadTrimmingCurveVector(
                rParameters["embedded_edges"],
                pNurbsSurface,
                rModelPart,
                EchoLevel);

            pBrepSurface->AddEmbeddedEdges(embedded_edges);
        }
    }
}

void ReadBrepSurfaceWithExistingNodes(
    const Parameters rParameters,
    ModelPart& rModelPart,
    const int EchoLevel)
{
    KRATOS_INFO_IF("CadGeometryReconstructionUtility", EchoLevel > 0)
        << "Reading BrepSurface \"" << GetIdOrName(rParameters)
        << "\" using existing nodes." << std::endl;

    KRATOS_ERROR_IF_NOT(rParameters.Has("surface"))
        << "Missing 'surface' in brep face." << std::endl;

    auto p_surface = ReadNurbsSurfaceWithExistingNodes(
        rParameters["surface"],
        rModelPart,
        EchoLevel);

    const bool is_trimmed = rParameters["surface"].Has("is_trimmed")
        ? rParameters["surface"]["is_trimmed"].GetBool()
        : true;

    if (rParameters.Has("boundary_loops"))
    {
        BrepCurveOnSurfaceLoopArrayType outer_loops;
        BrepCurveOnSurfaceLoopArrayType inner_loops;
        std::tie(outer_loops, inner_loops) = ReadBoundaryLoops(
            rParameters["boundary_loops"],
            p_surface,
            rModelPart,
            EchoLevel);

        auto p_brep_surface = Kratos::make_shared<BrepSurfaceType>(
            p_surface,
            outer_loops,
            inner_loops,
            is_trimmed);

        p_surface->SetGeometryParent(p_brep_surface.get());

        SetIdOrName<BrepSurfaceType>(rParameters, p_brep_surface);

        ReadAndAddEmbeddedEdges(
            p_brep_surface,
            rParameters,
            p_surface,
            rModelPart,
            EchoLevel);

        rModelPart.AddGeometry(p_brep_surface);
    }
    else
    {
        auto p_brep_surface = Kratos::make_shared<BrepSurfaceType>(p_surface);

        SetIdOrName<BrepSurfaceType>(rParameters, p_brep_surface);

        ReadAndAddEmbeddedEdges(
            p_brep_surface,
            rParameters,
            p_surface,
            rModelPart,
            EchoLevel);

        rModelPart.AddGeometry(p_brep_surface);
    }
}

void ReadBrepSurfacesWithExistingNodes(
    const Parameters rParameters,
    ModelPart& rModelPart,
    const int EchoLevel)
{
    KRATOS_ERROR_IF_NOT(rParameters.IsArray())
        << "\"faces\" section needs to be an array of BrepSurfaces." << std::endl;

    for (IndexType brep_surface_i = 0; brep_surface_i < rParameters.size(); ++brep_surface_i) {
        ReadBrepSurfaceWithExistingNodes(
            rParameters[brep_surface_i],
            rModelPart,
            EchoLevel);
    }
}

void ReadBrepsWithExistingNodes(
    const Parameters rParameters,
    ModelPart& rModelPart,
    const int EchoLevel)
{
    for (IndexType brep_index = 0; brep_index < rParameters.size(); ++brep_index)
    {
        KRATOS_INFO_IF("CadGeometryReconstructionUtility", EchoLevel > 0)
            << "Reading Brep \"" << GetIdOrName(rParameters[brep_index])
            << "\" - faces using existing nodes." << std::endl;

        if (rParameters[brep_index].Has("faces")) {
            ReadBrepSurfacesWithExistingNodes(
                rParameters[brep_index]["faces"],
                rModelPart,
                EchoLevel);
        }
    }
}

void ReadGeometryModelPartWithExistingNodes(
    const Parameters rCadJsonParameters,
    ModelPart& rModelPart,
    const int EchoLevel)
{
    KRATOS_ERROR_IF_NOT(rCadJsonParameters.Has("breps"))
        << "Missing \"breps\" section." << std::endl;

    ReadBrepsWithExistingNodes(
        rCadJsonParameters["breps"],
        rModelPart,
        EchoLevel);
}

} // unnamed namespace

void ReconstructModelPartBrepGeometryFromCadJson(
    const std::string& rDataFileName,
    ModelPart& rModelPart,
    const int EchoLevel)
{
    const Parameters cad_json_parameters = ReadParametersFile(rDataFileName, EchoLevel);
    ReconstructModelPartBrepGeometryFromCadJson(cad_json_parameters, rModelPart, EchoLevel);
}

void ReconstructModelPartBrepGeometryFromCadJson(
    const Parameters CadJsonParameters,
    ModelPart& rModelPart,
    const int EchoLevel)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rModelPart.NumberOfNodes() == 0)
        << "Trying to reconstruct BRep geometry in model part \""
        << rModelPart.FullName() << "\" but it has no nodes." << std::endl;

    KRATOS_INFO_IF("CadGeometryReconstructionUtility", EchoLevel > 0)
        << "Reconstructing BRep geometry in model part \""
        << rModelPart.FullName() << "\" from CAD json using existing nodes." << std::endl;

    ReadGeometryModelPartWithExistingNodes(
        CadJsonParameters,
        rModelPart,
        EchoLevel);

    KRATOS_CATCH("")
}

} // namespace CadGeometryReconstructionUtility
} // namespace Kratos