//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolo' Antonelli
//                   Andrea Gorgi
//


// System includes
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "iga_modeler_sbm.h"
#include "geometries/nurbs_curve_geometry.h"
#include "integration/integration_point_utilities.h"
#include "iga_application_variables.h"


namespace Kratos
{
namespace
{
using ProjectionCoordinates = array_1d<double, 3>;
using NurbsSkinCurve = NurbsCurveGeometry<2, PointerVector<Node>>;

struct NurbsProjectionResult
{
    Geometry<Node>::Pointer pGeometry;
    ProjectionCoordinates LocalCoordinates = ZeroVector(3);
    ProjectionCoordinates SkinCoordinates = ZeroVector(3);
    double SquaredDistance = std::numeric_limits<double>::max();
    std::string LayerName;
    IndexType InputOrder = 0;
    bool IsConvergedProjection = false;
};

// Finds the closest valid projection of a point onto the NURBS skin curves.
NurbsProjectionResult ProjectToNurbsSkin(
    const ProjectionCoordinates& rPoint,
    Model& rModel,
    const double SearchRadius,
    const double Tolerance,
    const std::string& rRequiredLayer = "")
{
    NurbsProjectionResult result;
    std::vector<Geometry<Node>::Pointer> curves;

    // Collects registered NURBS skin curves from the model's main model parts and sorts them by ID.
    for (const auto& r_model_part_name : rModel.GetModelPartNames()) {
        ModelPart& r_model_part = rModel.GetModelPart(r_model_part_name);
        if (r_model_part.IsSubModelPart()) {
            continue;
        }
        for (auto it = r_model_part.Geometries().ptr_begin();
             it != r_model_part.Geometries().ptr_end(); ++it) {
            auto* p_curve = dynamic_cast<NurbsSkinCurve*>((*it).get());
            if (!p_curve || !p_curve->Has(IDENTIFIER)) {
                continue;
            }
            if (!p_curve->Has(CONDITION_NAME)) {
                continue;
            }
            const std::string& r_condition_name =
                p_curve->GetValue(CONDITION_NAME);
            if (!KratosComponents<Condition>::Has(r_condition_name)) {
                continue;
            }
            curves.push_back(*it);
        }
    }
    std::sort(curves.begin(), curves.end(), [](const auto& a, const auto& b) {
        return a->Id() < b->Id();
    });

    const double admissible_radius_squared = std::pow(SearchRadius + Tolerance, 2);
    const auto evaluate = [&](const Geometry<Node>::Pointer& pGeometry,
                              const ProjectionCoordinates& rLocalCoordinates,
                              const IndexType InputOrder,
                              const bool IsConvergedProjection) {
        if (!std::isfinite(rLocalCoordinates[0])) {
            return;
        }
        ProjectionCoordinates position;
        pGeometry->GlobalCoordinates(position, rLocalCoordinates);
        const ProjectionCoordinates distance = position - rPoint;
        const double squared_distance = inner_prod(distance, distance);
        if (!std::isfinite(squared_distance) ||
            squared_distance > admissible_radius_squared) {
            return;
        }
        const double distance_norm = std::sqrt(squared_distance);
        const double best_distance_norm = std::sqrt(result.SquaredDistance);

        // Keeps the closest projection, breaking distance ties by input order.
        if (!result.pGeometry ||
            distance_norm < best_distance_norm - Tolerance ||
            (std::abs(distance_norm - best_distance_norm) <= Tolerance &&
             InputOrder < result.InputOrder)) {
            result.pGeometry = pGeometry;
            result.LocalCoordinates = rLocalCoordinates;
            result.SkinCoordinates = position;
            result.SquaredDistance = squared_distance;
            result.LayerName = pGeometry->GetValue(IDENTIFIER);
            result.InputOrder = InputOrder;
            result.IsConvergedProjection = IsConvergedProjection;
        }
    };

    IndexType input_order = 0;
    for (auto p_geometry : curves) {
        auto* p_curve = static_cast<NurbsSkinCurve*>(p_geometry.get());
        const std::string& r_layer = p_curve->GetValue(IDENTIFIER);
        if (!rRequiredLayer.empty() && r_layer != rRequiredLayer) {
            ++input_order;
            continue;
        }

        for (const auto& r_interval : p_curve->KnotSpanIntervals()) {
            const double lower = r_interval.GetT0();
            const double upper = r_interval.GetT1();
            if (upper <= lower) {
                continue;
            }

            // Evaluates the final projection iterate from the span midpoint and both endpoints.
            ProjectionCoordinates local_coordinates = ZeroVector(3);
            local_coordinates[0] = 0.5 * (lower + upper);
            const bool converged = p_curve->ProjectionPointGlobalToLocalSpace(
                rPoint, local_coordinates, Tolerance) != 0;
            evaluate(p_geometry, local_coordinates, input_order, converged);

            ProjectionCoordinates lower_coordinates = ZeroVector(3);
            lower_coordinates[0] = lower;
            evaluate(p_geometry, lower_coordinates, input_order, false);

            ProjectionCoordinates upper_coordinates = ZeroVector(3);
            upper_coordinates[0] = upper;
            evaluate(p_geometry, upper_coordinates, input_order, false);
        }
        ++input_order;
    }

    KRATOS_ERROR_IF_NOT(result.pGeometry)
        << "::[IgaModelerSbm]:: No valid NURBS projection for quadrature point "
        << rPoint << " within search radius " << SearchRadius
        << "." << std::endl;
    return result;
}

} // namespace
///@name Stages
///@{

void IgaModelerSbm::SetupModelPart()
{
    KRATOS_ERROR_IF_NOT(mParameters.Has("analysis_model_part_name"))
        << "Missing \"analysis_model_part_name\" section" << std::endl;
    ModelPart& analysis_model_part = mpModel->GetModelPart(mParameters["analysis_model_part_name"].GetString());

    KRATOS_ERROR_IF_NOT(mParameters.Has("element_condition_list"))
        << "Missing \"element_condition_list\" section" << std::endl;

    const Parameters iga_physics_parameters = mParameters["element_condition_list"];

    CreateIntegrationDomain(
        analysis_model_part,
        iga_physics_parameters);
    
    ActivateNodesInElementsAndCleanRoot(analysis_model_part);
}

///@}

void IgaModelerSbm::CreateIntegrationDomain(
    ModelPart& rModelPart,
    const Parameters rPhysicsParameters) const
{

    KRATOS_ERROR_IF_NOT(rPhysicsParameters.IsArray())
        << "\"element_condition_list\" needs to be an array." << std::endl;

    for (SizeType i = 0; i < rPhysicsParameters.size(); ++i)
    {
        CreateIntegrationDomainPerUnit(
            rModelPart,
            rPhysicsParameters[i]);
    }
}

void IgaModelerSbm::CreateIntegrationDomainPerUnit(
    ModelPart& rModelPart,
    const Parameters rPhysicsParameters) const
{
    KRATOS_ERROR_IF_NOT(rPhysicsParameters.Has("iga_model_part"))
        << "::[IgaModelerSbm]:: \"iga_model_part\" needs to be specified." << std::endl;

    KRATOS_ERROR_IF_NOT(rPhysicsParameters.Has("shape_function_derivatives_order"))
        << "::[IgaModelerSbm]:: \"shape_function_derivatives_order\" needs to be specified." << std::endl;

    std::string sub_model_part_name = rPhysicsParameters["iga_model_part"].GetString();

    ModelPart& sub_model_part = rModelPart.HasSubModelPart(sub_model_part_name)
        ? rModelPart.GetSubModelPart(sub_model_part_name)
        : rModelPart.CreateSubModelPart(sub_model_part_name);

    // Generate the list of geometries, which are needed, here.
    GeometriesArrayType geometry_list;
    GetGeometryList(geometry_list, rModelPart, rPhysicsParameters);
    
    KRATOS_ERROR_IF_NOT(rPhysicsParameters.Has("geometry_type")) << 
         "::[IgaModelerSbm]:: Missing \"geometry_type\" parameter." << rPhysicsParameters << std::endl;
                
    std::string geometry_type = rPhysicsParameters["geometry_type"].GetString();

    KRATOS_ERROR_IF(
        geometry_type != "GeometrySurface" &&
        geometry_type != "SurfaceEdge")
        << "::[IgaModelerSbm]:: Unsupported \"geometry_type\": \""
        << geometry_type
        << "\". Available options: GeometrySurface, SurfaceEdge."
        << std::endl;

    if (!rPhysicsParameters.Has("sbm_parameters"))
        CreateQuadraturePointGeometries(
            geometry_list, sub_model_part, rPhysicsParameters, geometry_type); 
    else 
        CreateQuadraturePointGeometriesSbm(
            geometry_list, sub_model_part, rPhysicsParameters, geometry_type);

    KRATOS_INFO_IF("CreateIntegrationDomainElementCondition", mEchoLevel > 3)
        << "Creation of elements/ conditions finished in: " << sub_model_part << std::endl;
}

///@}
///@name CAD functionalities
///@{

void IgaModelerSbm::GetGeometryList(
    GeometriesArrayType& rGeometryList,
    ModelPart& rModelPart,
    const Parameters rPhysicsParameters) const
{
    /* we have three cases:
        1) background geometry: i.e. the surface or volume. they are type = "element"
        2) outer_loop -> is_inner = false : type = "condition" (2D brep curves & 3D brep surfaces)
        2) inner_loop -> is_inner = true  : type = "condition" (2D brep curves & 3D brep surfaces)
    */ 

    const std::string type = rPhysicsParameters["type"].GetString();

    // get the dimension of the geometry
    SizeType dim;
    // Get the space dimension 
    if (rModelPart.GetValue(KNOT_VECTOR_W).size() == 0) {
        dim = 2;
    } else {
        dim = 3;
    }

    if (type == "element")
    {
        int surface_brep_id = 1; 
        rGeometryList.push_back(rModelPart.pGetGeometry(surface_brep_id));
    } 
    else if (type == "condition")
    {
        if (rPhysicsParameters.Has("sbm_parameters"))
        {
            KRATOS_ERROR_IF_NOT(rPhysicsParameters["sbm_parameters"].Has("is_inner")) << 
                "::[IgaModelerSbm]:: Missing \"is_inner\" parameter in the \"sbm_parameters\"." << rPhysicsParameters["sbm_parameters"] 
                << std::endl;
                        
            if (rPhysicsParameters["sbm_parameters"]["is_inner"].GetBool()) // inner loop
            {
                // Surface id is 1
                int inner_brep_id = 2;
                ModelPart& surrogate_model_part_outer = rModelPart.GetSubModelPart("surrogate_outer");
                if (surrogate_model_part_outer.NumberOfConditions() == 0)
                    if (dim == 2) // 2D case
                        inner_brep_id += 4; // if there is no outer we use the 4 sides of the rectangle
                    else if (dim == 3) // 3D case
                        inner_brep_id += 6; // if there is no outer we use the 6 sides of the parallelepiped   
                    else
                        KRATOS_ERROR << "::[IgaModelerSbm]:: The dimension of the geometry is not 2 or 3." << std::endl;                 
                else 
                    // if outer loop is present take the number of conditions
                    inner_brep_id += surrogate_model_part_outer.NumberOfConditions();

                // INNER   
                ModelPart& surrogate_model_part_inner = rModelPart.GetSubModelPart("surrogate_inner");

                KRATOS_ERROR_IF(surrogate_model_part_inner.NumberOfElements() == 0) 
                    << "::[IgaModelerSbm]:: The surrogate_model_part_inner has zero elements (no inner loop/boundary defined)."
                    << "Something might be missing in the NurbsModelerSbm." << std::endl;

                for (const auto& rElem : surrogate_model_part_inner.Elements()) {
                    /*
                    Each element in the surrogate_model_part_inner represents a surrogate boundary loop. First "node.Id()" is the id of the first condition and
                        the second "node.Id()" is the last condition of that loop. (Essential for multiple inner loops)
                    */
                    const auto& r_geometry = rElem.GetGeometry();
                    KRATOS_ERROR_IF(r_geometry.PointsNumber() < 2)
                        << "Surrogate loop element " << rElem.Id() << " has <2 geometry points." << std::endl;

                    // First/last condition IDs encoded as the first two geometry nodes
                    const IndexType first_condition_id = r_geometry[0].Id();
                    const IndexType last_condition_id  = r_geometry[1].Id();

                    SizeType size_surrogate_loop = last_condition_id - first_condition_id + 1;

                    for (SizeType j = 0; j < size_surrogate_loop; ++j) {
                        rGeometryList.push_back(rModelPart.pGetGeometry(inner_brep_id));
                        inner_brep_id++;
                    }
                }
            } else // outer loop
            {
                int outer_brep_id = 2;
                // OUTER
                ModelPart& surrogate_model_part_outer = rModelPart.GetSubModelPart("surrogate_outer");
                
                if (surrogate_model_part_outer.NumberOfConditions() > 0) {
                    // both for 2D and 3D
                    const int size_surrogate_loop_outer = surrogate_model_part_outer.NumberOfConditions();
                    for (int j = 0; j < size_surrogate_loop_outer; ++j) {
                        rGeometryList.push_back(rModelPart.pGetGeometry(outer_brep_id));
                        outer_brep_id++;
                    }

                }
            }
        }
        // else -> body-fitted case
        else {
            KRATOS_ERROR_IF_NOT(rPhysicsParameters.Has("brep_ids")) << 
                "::[IgaModelerSbm]:: Missing \"brep_ids\" parameter in body-fitted boundary condition\"." << rPhysicsParameters 
                << std::endl;
            if (rPhysicsParameters.Has("brep_ids")) {
                for (SizeType i = 0; i < rPhysicsParameters["brep_ids"].size(); ++i) {
                    rGeometryList.push_back(rModelPart.pGetGeometry(rPhysicsParameters["brep_ids"][i].GetInt()));
                }
            }
        }
    } else 
    {
        KRATOS_ERROR << "::[IgaModelerSbm]:: type " << type << " not defined in the IgaModelerSbm."
                     << " Available types are: elements, conditions." << std::endl;
    }

    KRATOS_ERROR_IF(rGeometryList.size() == 0)
        << "::[IgaModelerSbm]:: Empty geometry list in GetGeometryList. Physics parameter is: " << rPhysicsParameters << std::endl;
}

void IgaModelerSbm::CreateQuadraturePointGeometries(
    GeometriesArrayType& rGeometryList,
    ModelPart& rModelPart,
    const Parameters rParameters,
    std::string GeometryType) const
{
    KRATOS_ERROR_IF_NOT(rParameters.Has("type"))
        << "::[IgaModelerSbm]:: \"type\" need to be specified." << std::endl;
    std::string type = rParameters["type"].GetString();
    KRATOS_ERROR_IF_NOT(rParameters.Has("name"))
        << "::[IgaModelerSbm]:: \"name\" need to be specified." << std::endl;
    std::string name = rParameters["name"].GetString();

    const int shape_function_derivatives_order =
    rParameters["shape_function_derivatives_order"].GetInt();

   KRATOS_ERROR_IF(shape_function_derivatives_order < 1)
    << "::[IgaModelerSbm]:: \"shape_function_derivatives_order\" "
    << "must be at least 1, but received "
    << shape_function_derivatives_order << std::endl;

    std::string quadrature_method = rParameters.Has("quadrature_method")
        ? rParameters["integration_rule"].GetString()
        : "GAUSS";

    KRATOS_INFO_IF("CreateQuadraturePointGeometries", mEchoLevel > 0)
        << "Creating " << name << "s of type: " << type
        << " for " << rGeometryList.size() << " geometries"
        << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

    
    for (SizeType i = 0; i < rGeometryList.size(); ++i)
    {
        GeometriesArrayType geometries;
        IntegrationInfo integration_info = rGeometryList[i].GetDefaultIntegrationInfo();
        for (IndexType i = 0; i < integration_info.LocalSpaceDimension(); ++i) {
            if (quadrature_method == "GAUSS") {
                integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GAUSS);
            }
            else if (quadrature_method == "GRID") {
                integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GRID);
            }
            else {
                KRATOS_INFO("CreateQuadraturePointGeometries") << "Quadrature method: " << quadrature_method
                    << " is not available. Available options are \"GAUSS\" and \"GRID\". Default quadrature method is being considered." << std::endl;
            }
        }

        if (rParameters.Has("number_of_integration_points_per_span")) {
            for (IndexType i = 0; i < integration_info.LocalSpaceDimension(); ++i) {
                integration_info.SetNumberOfIntegrationPointsPerSpan(i, rParameters["number_of_integration_points_per_span"].GetInt());
            }
        }
        if (GeometryType == "SurfaceEdge"
            && rGeometryList[i].GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Coupling_Geometry)
        {
            rGeometryList[i].GetGeometryPart(0).CreateQuadraturePointGeometries(
                geometries, shape_function_derivatives_order, integration_info);
        }
        else
        {
            rGeometryList[i].CreateQuadraturePointGeometries(
                geometries, shape_function_derivatives_order, integration_info);
        }

        KRATOS_INFO_IF("CreateQuadraturePointGeometries", mEchoLevel > 1)
            << geometries.size() << " quadrature point geometries have been created." << std::endl;

        if (type == "element") {
            // Get the mesh sizes from the iga model part (fallback if not on parent).
            Vector knot_span_sizes;
            if (rModelPart.GetParentModelPart().Has(KNOT_SPAN_SIZES)) {
                knot_span_sizes = rModelPart.GetParentModelPart().GetValue(KNOT_SPAN_SIZES);
            } else if (rModelPart.Has(KNOT_SPAN_SIZES)) {
                knot_span_sizes = rModelPart.GetValue(KNOT_SPAN_SIZES);
            } else if (rModelPart.GetRootModelPart().Has(KNOT_SPAN_SIZES)) {
                knot_span_sizes = rModelPart.GetRootModelPart().GetValue(KNOT_SPAN_SIZES);
            }

            SizeType id = 1;
            if (rModelPart.GetRootModelPart().Elements().size() > 0)
                id = rModelPart.GetRootModelPart().Elements().back().Id() + 1;

            this->CreateElements(
                geometries.ptr_begin(), geometries.ptr_end(),
                rModelPart, name, id, PropertiesPointerType(), knot_span_sizes);
        }
        else if (type == "condition") {
            // Get the mesh sizes from the iga model part (fallback if not on parent).
            Vector knot_span_sizes;
            if (rModelPart.GetParentModelPart().Has(KNOT_SPAN_SIZES)) {
                knot_span_sizes = rModelPart.GetParentModelPart().GetValue(KNOT_SPAN_SIZES);
            } else {
                KRATOS_ERROR << "KNOT_SPAN_SIZES not found in parent model part." << std::endl;
            }

            SizeType id = 1;
            if (rModelPart.GetRootModelPart().Conditions().size() > 0)
                id = rModelPart.GetRootModelPart().Conditions().back().Id() + 1;
            this->CreateConditions(
                geometries.ptr_begin(), geometries.ptr_end(),
                rModelPart, name, id, PropertiesPointerType(), knot_span_sizes);
        }
        else {
            KRATOS_ERROR << "\"type\" does not exist: " << type
                << ". Possible types are \"element\" and \"condition\"." << std::endl;
        }
    }
}

void IgaModelerSbm::CreateQuadraturePointGeometriesSbm(
    GeometriesArrayType& rGeometryList,
    ModelPart& rModelPart,
    const Parameters rParameters,
    std::string GeometryType) const
{
    KRATOS_ERROR_IF_NOT(rParameters.Has("name"))
        << "\"name\" needs to be specified." << std::endl;
                            
    const std::string name = rParameters["name"].GetString();

    if (name == "SbmCondition") {
        CreateQuadraturePointGeometriesSbmByProjectionLayer(
            rGeometryList, rModelPart, rParameters);
    } else {
        CreateQuadraturePointGeometriesSbmByFixedConditionName(
            rGeometryList, rModelPart, rParameters, GeometryType, name);
    }
}

void IgaModelerSbm::CreateQuadraturePointGeometriesSbmByProjectionLayer(
    GeometriesArrayType& rGeometryList,
    ModelPart& rModelPart,
    const Parameters rParameters) const
{
    KRATOS_ERROR_IF(rParameters["type"].GetString() != "condition")
        << "::[IgaModelerSbm]:: SBM operators require type \"condition\"."
        << std::endl;

    const bool is_inner = rParameters["sbm_parameters"]["is_inner"].GetBool();
    const std::string loop_name = is_inner ? "inner" : "outer";

    KRATOS_ERROR_IF_NOT(mParameters.Has("skin_model_part_name"))
        << "::[IgaModelerSbm]:: Missing \"skin_model_part_name\" in "
        << "modeler parameters." << std::endl;

    const std::string skin_model_part_name =
        mParameters["skin_model_part_name"].GetString();

    ModelPart& r_skin_loop = mpModel->GetModelPart(
        skin_model_part_name).GetSubModelPart(loop_name);

    const Vector& r_knot_span_sizes =
        rModelPart.GetParentModelPart().GetValue(KNOT_SPAN_SIZES);

    KRATOS_ERROR_IF(r_knot_span_sizes.size() < 2)
        << "::[IgaModelerSbm]:: Two KNOT_SPAN_SIZES are required in 2D."
        << std::endl;

    const double h = std::max(r_knot_span_sizes[0], r_knot_span_sizes[1]);
    const double search_radius = std::sqrt(2.0) * h;
    const double tolerance = std::max(1e-6, 1e-3 * h);

    const int derivatives_order = rParameters["shape_function_derivatives_order"].GetInt();

    std::string quadrature_method = rParameters.Has("quadrature_method")
        ? rParameters["integration_rule"].GetString()
        : "GAUSS";

    KRATOS_ERROR_IF(quadrature_method != "GAUSS" &&
                    quadrature_method != "GRID")
        << "::[IgaModelerSbm]:: Unsupported quadrature method: "
        << quadrature_method << "." << std::endl;

    //Creation or retrieval of Id identifiers for conditions and nodes
    SizeType condition_id = rModelPart.GetRootModelPart().NumberOfConditions() == 0
        ? 1 : rModelPart.GetRootModelPart().Conditions().back().Id() + 1;
    SizeType node_id = r_skin_loop.GetRootModelPart().NumberOfNodes() == 0
        ? 1 : r_skin_loop.GetRootModelPart().Nodes().back().Id() + 1;

    for (auto& r_surrogate_geometry : rGeometryList) {
        const SizeType required_derivatives_order =
            r_surrogate_geometry.pGetGeometryPart(
                Geometry<Node>::BACKGROUND_GEOMETRY_INDEX)->
                    PolynomialDegree(0) + 1;
        KRATOS_ERROR_IF(derivatives_order < 0 || static_cast<SizeType>(derivatives_order) < required_derivatives_order)
            << "::[IgaModelerSbm]:: shape_function_derivatives_order must be "
            << "at least " << required_derivatives_order << "." << std::endl;


        GeometriesArrayType quadrature_geometries;
        IntegrationInfo integration_info = r_surrogate_geometry.GetDefaultIntegrationInfo();

        for (IndexType d = 0; d < integration_info.LocalSpaceDimension(); ++d)
        {
            integration_info.SetQuadratureMethod(
                d, quadrature_method == "GRID"
                ? IntegrationInfo::QuadratureMethod::GRID
                : IntegrationInfo::QuadratureMethod::GAUSS);

            if (rParameters.Has("number_of_integration_points_per_span")) {
                integration_info.SetNumberOfIntegrationPointsPerSpan(
                    d, rParameters[
                        "number_of_integration_points_per_span"].GetInt());
            }
        }

        r_surrogate_geometry.CreateQuadraturePointGeometries
        (quadrature_geometries, derivatives_order, integration_info);
        KRATOS_ERROR_IF(quadrature_geometries.empty())
            << "::[IgaModelerSbm]:: No quadrature point geometries were created."<< std::endl;


        std::vector<NurbsProjectionResult> projections;
        std::vector<std::string> layers;
        std::vector<SizeType> votes;

        projections.reserve(quadrature_geometries.size());
        for (const auto& r_quadrature_geometry : quadrature_geometries)
        {
            auto projection = ProjectToNurbsSkin(
                r_quadrature_geometry.Center().Coordinates(),
                *mpModel,
                search_radius,
                tolerance);

            //Layer voting process
            const auto layer_it = std::find(layers.begin(), layers.end(), projection.LayerName);
            if (layer_it == layers.end()) {
                layers.push_back(projection.LayerName);
                votes.push_back(1);
            } else {
                ++votes[std::distance(layers.begin(), layer_it)];
            }
            projections.push_back(std::move(projection));
        }

        const IndexType winner = std::distance(
            votes.begin(), std::max_element(votes.begin(), votes.end()));

        const std::string winning_layer = layers[winner];
        for (IndexType j = 0; j < projections.size(); ++j)
        // Re-project quadrature points to the winning layer closest point
        {
            if (projections[j].LayerName != winning_layer) {
                projections[j] = ProjectToNurbsSkin(
                    quadrature_geometries[j].Center().Coordinates(),
                    *mpModel,
                    search_radius,
                    tolerance,
                    winning_layer);
            }
        }
        // Create or retrieve the submodelparts for the winning layer and the projection data.
        ModelPart& r_condition_layer =
            rModelPart.HasSubModelPart(winning_layer)
            ? rModelPart.GetSubModelPart(winning_layer)
            : rModelPart.CreateSubModelPart(winning_layer);
        ModelPart& r_skin_layer =
            r_skin_loop.HasSubModelPart(winning_layer)
            ? r_skin_loop.GetSubModelPart(winning_layer)
            : r_skin_loop.CreateSubModelPart(winning_layer);
        // Keep projection nodes in the skin hierarchy for nodal boundary-value processes.
        ModelPart& r_projection_data = r_skin_layer.HasSubModelPart("projection_data")
            ? r_skin_layer.GetSubModelPart("projection_data")
            : r_skin_layer.CreateSubModelPart("projection_data");

        for (IndexType j = 0; j < projections.size(); ++j) {
            const auto& r_projection = projections[j];
            KRATOS_WARNING_IF("IgaModelerSbm", !r_projection.IsConvergedProjection)
                << "Closest-point convergence was not verified for quadrature point "
                << quadrature_geometries[j].Center().Coordinates()
                << ". Using the closest evaluated point on the NURBS skin "
                << r_projection.SkinCoordinates << " (distance "
                << std::sqrt(r_projection.SquaredDistance) << ")." << std::endl;
            KRATOS_ERROR_IF_NOT(r_projection.pGeometry->Has(CONDITION_NAME))
                << "::[IgaModelerSbm]:: Missing CONDITION_NAME on NURBS geometry "
                << r_projection.pGeometry->Id() << "." << std::endl;
            const std::string& condition_name = r_projection.pGeometry->GetValue(CONDITION_NAME);
            KRATOS_ERROR_IF(condition_name.empty())
                << "::[IgaModelerSbm]:: Empty condition name for NURBS geometry "
                << r_projection.pGeometry->Id() << "." << std::endl;
            KRATOS_ERROR_IF_NOT(KratosComponents<Condition>::Has(condition_name))
                << condition_name << " not registered." << std::endl;

            const Condition& r_reference_condition =
                KratosComponents<Condition>::Get(condition_name);
            auto p_condition = r_reference_condition.Create(
                condition_id++, quadrature_geometries(j),
                PropertiesPointerType());

            p_condition->SetValue(IDENTIFIER, loop_name);
            p_condition->SetValue(KNOT_SPAN_SIZES, r_knot_span_sizes);

            const auto& r_position = r_projection.SkinCoordinates;
            // Nodal boundary-value processes on the skin hierarchy include this projection.
            auto p_projection_node = r_projection_data.CreateNewNode(
                node_id++, r_position[0], r_position[1], r_position[2]);

            std::vector<ProjectionCoordinates> derivatives;

            r_projection.pGeometry->GlobalSpaceDerivatives(
                derivatives, r_projection.LocalCoordinates, 1);

            const double tangent_norm = norm_2(derivatives[1]);
            KRATOS_ERROR_IF_NOT(std::isfinite(tangent_norm) && tangent_norm > 0.0)
                << "::[IgaModelerSbm]:: Degenerate tangent on NURBS geometry "
                << r_projection.pGeometry->Id() << "." << std::endl;

            const ProjectionCoordinates tangent = derivatives[1] / tangent_norm;
            ProjectionCoordinates normal = ZeroVector(3);

            const double normal_sign = is_inner ? -1.0 : 1.0;
            normal[0] = normal_sign * tangent[1];
            normal[1] = -normal_sign * tangent[0];

            // Assign the normal to the projection node and set NEIGHBOUR_NODES.
            p_projection_node->SetValue(NORMAL, normal);
            p_condition->SetValue(NEIGHBOUR_NODES,
                GlobalPointersVector<Node>({p_projection_node}));
            r_condition_layer.AddCondition(p_condition);
        }
    }
}

void IgaModelerSbm::CreateQuadraturePointGeometriesSbmByFixedConditionName(
    GeometriesArrayType& rGeometryList,
    ModelPart& rModelPart,
    const Parameters rParameters,
    std::string GeometryType,
    std::string ConditionName) const
{
    KRATOS_ERROR_IF_NOT(rParameters.Has("type"))
        << "\"type\" needs to be specified." << std::endl;

    // Only conditions should call CreateQuadraturePointGeometriesSbm
    std::string type = rParameters["type"].GetString();
    bool check_input_type = (type == "condition");
    KRATOS_ERROR_IF_NOT(check_input_type) << ":::[IgaModelerSbm]::: type != \"condition\" in CreateQuadraturePointGeometriesSbm. "
                                          << "It must be a condition to apply the sbm operators. type: " << type << std::endl;

    const int shape_function_derivatives_order =
    rParameters["shape_function_derivatives_order"].GetInt();

    const SizeType required_shape_function_derivatives_order =
    rGeometryList[0].pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)-> PolynomialDegree(0) + 1;

    KRATOS_ERROR_IF(
        shape_function_derivatives_order < 0 ||
        static_cast<SizeType>(shape_function_derivatives_order) < required_shape_function_derivatives_order)
        << "::[IgaModelerSbm]:: \"shape_function_derivatives_order\" must be at least " << required_shape_function_derivatives_order << ", but received " << shape_function_derivatives_order << std::endl;

    std::string quadrature_method = rParameters.Has("quadrature_method")
        ? rParameters["integration_rule"].GetString()
        : "GAUSS";

    KRATOS_INFO_IF("CreateQuadraturePointGeometries", mEchoLevel > 0)
        << "Creating " << ConditionName << "s of type: " << type
        << " for " << rGeometryList.size() << " geometries"
        << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

    // Check if the sbm projection operation is needed (there is no need for background domain and for body-fitted boundary conditions)
    PointVector points;   
    const std::string skin_model_part_name = mParameters.Has("skin_model_part_name")
            ? mParameters["skin_model_part_name"].GetString()
            : "skin_model_part";

    ModelPart& skin_model_part = mpModel->HasModelPart(skin_model_part_name)
            ? mpModel->GetModelPart(skin_model_part_name)
            : KRATOS_ERROR << "::[CreateQuadraturePointGeometriesSbm]::: Sbm case -> skin_model_part has not been defined before. "
                            << "Maybe you are not calling the nurbs_modeler_sbm" << std::endl;

    // inner & outer are defaulf sub model part names
    auto& skin_sub_model_part_in = skin_model_part.GetSubModelPart("inner");
    auto& skin_sub_model_part_out = skin_model_part.GetSubModelPart("outer");

    const bool is_inner = rParameters["sbm_parameters"]["is_inner"].GetBool();
    
    const std::string surrogate_sub_model_part_name = is_inner ? "surrogate_inner" : "surrogate_outer";
    
    ModelPart& surrogate_sub_model_part = rModelPart.GetParentModelPart().HasSubModelPart(surrogate_sub_model_part_name)
            ? rModelPart.GetParentModelPart().GetSubModelPart(surrogate_sub_model_part_name)
            : KRATOS_ERROR << "::[CreateQuadraturePointGeometriesSbm]::: Sbm case -> surrogate_sub_model_part has not been defined before."
                            << "Maybe you are not calling the nurbs_modeler_sbm" << std::endl;

    if (is_inner) { // INNER
        for (auto &i_cond : skin_sub_model_part_in.Conditions()) {
            points.push_back(PointTypePointer(new PointType(i_cond.Id(), i_cond.GetGeometry().Center().X(), i_cond.GetGeometry().Center().Y(), i_cond.GetGeometry().Center().Z())));
        }
    } 
    else { // OUTER
        for (auto &i_cond : skin_sub_model_part_out.Conditions()) {
            points.push_back(PointTypePointer(new PointType(i_cond.Id(), i_cond.GetGeometry().Center().X(), i_cond.GetGeometry().Center().Y(), i_cond.GetGeometry().Center().Z())));
        }
    }
    
    // Get the mesh sizes from the surrogate model part
    const Vector& knot_span_sizes = surrogate_sub_model_part.GetParentModelPart().GetValue(KNOT_SPAN_SIZES);

    double knot_span_reference_size = knot_span_sizes[0];
    if (knot_span_sizes[1] > knot_span_reference_size) {knot_span_reference_size = knot_span_sizes[1];}
    if (knot_span_sizes.size() > 2) {if (knot_span_sizes[2] > knot_span_reference_size) {knot_span_reference_size = knot_span_sizes[2];}}

    const int domain_size = rModelPart.GetProcessInfo()[DOMAIN_SIZE];
    double search_radius;
    if (domain_size == 2) {
        search_radius = 2*std::sqrt(2.0) * knot_span_reference_size;
    } else {
        search_radius = 3*std::sqrt(3.0) * knot_span_reference_size;
    }

    DynamicBins testBins(points.begin(), points.end());
    
    // Maximum number of results to be found in the search in radius
    const int number_of_results = 1e6; 

    ModelPart::NodesContainerType::ContainerType results(number_of_results);
    std::vector<double> list_of_distances(number_of_results);
    for (SizeType i = 0; i < rGeometryList.size(); ++i)
    {
        GeometriesArrayType geometries;
        IntegrationInfo integration_info = rGeometryList[i].GetDefaultIntegrationInfo();
        for (IndexType i = 0; i < integration_info.LocalSpaceDimension(); ++i) {
            if (quadrature_method == "GAUSS") {
                integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GAUSS);
            }
            else if (quadrature_method == "GRID") {
                integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GRID);
            }
            else {
                KRATOS_INFO("CreateQuadraturePointGeometries") << "Quadrature method: " << quadrature_method
                    << " is not available. Available options are \"GAUSS\" and \"GRID\". Default quadrature method is being considered." << std::endl;
            }
        }

        if (rParameters.Has("number_of_integration_points_per_span")) {
            for (IndexType i = 0; i < integration_info.LocalSpaceDimension(); ++i) {
                integration_info.SetNumberOfIntegrationPointsPerSpan(i, rParameters["number_of_integration_points_per_span"].GetInt());
            }
        }
    
        rGeometryList[i].CreateQuadraturePointGeometries(geometries, shape_function_derivatives_order, integration_info);

        KRATOS_INFO_IF("CreateQuadraturePointGeometries", mEchoLevel > 1)
            << geometries.size() << " quadrature point geometries have been created." << std::endl;

        SizeType id = 1;
        if (rModelPart.GetRootModelPart().NumberOfConditions() > 0)
            id = rModelPart.GetRootModelPart().Conditions().back().Id() + 1;
        
        std::vector<int> list_id_closest_condition(geometries.size());

        for (SizeType j= 0; j < geometries.size() ; j++) {  

            const Point integration_point = geometries[j].Center(); 
            PointerType p_integration_point = PointerType(new PointType(1, integration_point.X(), integration_point.Y(), integration_point.Z()));

            // Use the search in radius to find the closest point
            SizeType obtained_results = testBins.SearchInRadius(*p_integration_point, search_radius, results.begin(), list_of_distances.begin(), number_of_results);

            double minimum_distance = 1e14;

            // Find the nearest node
            IndexType nearest_node_id;
            for (IndexType k = 0; k < obtained_results; k++) {
                double current_distance = list_of_distances[k];   
                if (current_distance < minimum_distance) { 
                    minimum_distance = current_distance;
                    nearest_node_id = k;
                }
            }
            KRATOS_ERROR_IF(obtained_results == 0) << "::[IgaModelerSbm]:: Zero points found in serch for projection of point: " <<
                integration_point << std::endl;
            
            // store id closest condition
            list_id_closest_condition[j] = results[nearest_node_id]->Id();
        }

        if (is_inner) {
            // pass the skin_sub_model_part_in
            this->CreateConditions( geometries.ptr_begin(), geometries.ptr_end(),
                rModelPart, skin_sub_model_part_in, list_id_closest_condition, ConditionName, id, PropertiesPointerType(), is_inner, knot_span_sizes);
        }
        else{
            // pass the skin_sub_model_part_out
            this->CreateConditions(geometries.ptr_begin(), geometries.ptr_end(),
                rModelPart, skin_sub_model_part_out, list_id_closest_condition, ConditionName, id, PropertiesPointerType(), is_inner, knot_span_sizes);
        }
    }
}

///@}
///@name Generate Elements and Conditions
///@{

void IgaModelerSbm::CreateElements(
    typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
    typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
    ModelPart& rModelPart,
    std::string& rElementName,
    SizeType& rIdCounter,
    PropertiesPointerType pProperties,
    const Vector KnotSpanSizes) const
{
    KRATOS_ERROR_IF(!KratosComponents<Element>::Has(rElementName))
        << rElementName << " not registered." << std::endl;

    const Element& rReferenceElement = KratosComponents<Element>::Get(rElementName);

    ElementsContainerType new_element_list;

    KRATOS_INFO_IF("CreateElements", mEchoLevel > 2)
        << "Creating elements of type " << rElementName
        << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

    SizeType num_elements = std::distance(rGeometriesBegin, rGeometriesEnd);
    new_element_list.reserve(num_elements);

    int count = 0;
    for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it)
    {
        new_element_list.push_back(
            rReferenceElement.Create(rIdCounter, (*it), pProperties));


        // Set knot span sizes to the condition
        new_element_list.GetContainer()[count]->SetValue(KNOT_SPAN_SIZES, KnotSpanSizes);

        rIdCounter++;
        count++;
    }

    rModelPart.AddElements(new_element_list.begin(), new_element_list.end());
}

void IgaModelerSbm::CreateConditions(
    typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
    typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
    ModelPart& rModelPart,
    std::string& rConditionName,
    SizeType& rIdCounter,
    PropertiesPointerType pProperties,
    const Vector KnotSpanSizes) const
{
    const Condition& reference_condition = KratosComponents<Condition>::Get(rConditionName);

    ModelPart::ConditionsContainerType new_condition_list;

    KRATOS_INFO_IF("CreateConditions", mEchoLevel > 2)
        << "Creating conditions of type " << rConditionName
        << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

    int count_list_closest_condition = 0;
    for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it)
    {
        new_condition_list.push_back(
            reference_condition.Create(rIdCounter, (*it), pProperties));
        
        // Set knot span sizes to the condition
        new_condition_list.GetContainer()[count_list_closest_condition]->SetValue(KNOT_SPAN_SIZES, KnotSpanSizes);

        rIdCounter++;
        count_list_closest_condition++;
    }

    rModelPart.AddConditions(new_condition_list.begin(), new_condition_list.end());
}


void IgaModelerSbm::CreateConditions(
    typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
    typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
    ModelPart& rModelPart,
    ModelPart& rSkinModelPart,
    std::vector<int>& listIdClosestCondition,
    std::string& rConditionName,
    SizeType& rIdCounter,
    PropertiesPointerType pProperties,
    bool IsInner,
    const Vector KnotSpanSizes) const
{
    const Condition& reference_condition = KratosComponents<Condition>::Get(rConditionName);

    ModelPart::ConditionsContainerType new_condition_list;

    KRATOS_INFO_IF("CreateConditions", mEchoLevel > 2)
        << "Creating conditions of type " << rConditionName
        << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

    int count_list_closest_condition = 0;

    // 2D case
    if (rSkinModelPart.ConditionsBegin()->GetGeometry().size() == 2) {

        for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it) {
            new_condition_list.push_back(reference_condition.Create(rIdCounter, (*it), pProperties));

            IndexType condId = listIdClosestCondition[count_list_closest_condition];

            Condition::Pointer cond1 = &rSkinModelPart.GetCondition(condId);

            IndexType condId2;
            if (condId == rSkinModelPart.ConditionsBegin()->Id()) condId2 = (rSkinModelPart.ConditionsEnd()-1)->Id();
            else condId2 = condId-1;

            Condition::Pointer cond2 = &rSkinModelPart.GetCondition(condId2);

            new_condition_list.GetContainer()[count_list_closest_condition]->SetValue(NEIGHBOUR_CONDITIONS, GlobalPointersVector<Condition>({cond1,cond2}));
            if (IsInner) {
                new_condition_list.GetContainer()[count_list_closest_condition]->SetValue(IDENTIFIER, "inner");
            } else {
                new_condition_list.GetContainer()[count_list_closest_condition]->SetValue(IDENTIFIER, "outer");
            }
            new_condition_list.GetContainer()[count_list_closest_condition]->SetValue(KNOT_SPAN_SIZES, KnotSpanSizes);
                        
            rIdCounter++;
            count_list_closest_condition++;
        }
    } else {
        // 3D case
        for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it) {
            new_condition_list.push_back(reference_condition.Create(rIdCounter, (*it), pProperties));

            IndexType condId = listIdClosestCondition[count_list_closest_condition];
            Condition::Pointer cond1 = &rSkinModelPart.GetCondition(condId);

            new_condition_list.GetContainer()[count_list_closest_condition]->SetValue(NEIGHBOUR_CONDITIONS, GlobalPointersVector<Condition>({cond1}));
            if (IsInner) {
                new_condition_list.GetContainer()[count_list_closest_condition]->SetValue(IDENTIFIER, "inner");
            } else {
                new_condition_list.GetContainer()[count_list_closest_condition]->SetValue(IDENTIFIER, "outer");
            }
            new_condition_list.GetContainer()[count_list_closest_condition]->SetValue(KNOT_SPAN_SIZES, KnotSpanSizes);
            rIdCounter++;
            count_list_closest_condition++;
        }
    }
    
    rModelPart.AddConditions(new_condition_list.begin(), new_condition_list.end());
}

void IgaModelerSbm::CreateConditions(
    typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
    typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
    ModelPart& rModelPart,
    ModelPart& rSkinModelPart,
    std::vector<int>& rListIdClosestCondition,
    std::vector<int>& rListIdSecondClosestCondition,
    SizeType& rIdCounter,
    PropertiesPointerType pProperties,
    const bool IsInner,
    const Vector KnotSpanSizes) const
{

    ModelPart::ConditionsContainerType new_condition_list;

    int count_list_closest_condition = 0;

    // count how many conditions for each layer
    std::vector<std::string> bc_on_skin_projection_type;
    std::vector<SizeType> count_bc_on_skin_projection_type;
    for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it) {
        std::string condition_layer_name = rSkinModelPart.GetCondition(rListIdClosestCondition[count_list_closest_condition]).GetValue(LAYER_NAME);

        // Find the condition name in bc_on_skin_projection_type
        auto it_name = std::find(bc_on_skin_projection_type.begin(), bc_on_skin_projection_type.end(), condition_layer_name);
        if (it_name != bc_on_skin_projection_type.end()) {
            // Increment the count for the existing condition name
            size_t index = std::distance(bc_on_skin_projection_type.begin(), it_name);
            count_bc_on_skin_projection_type[index]++;
        } else {
            // Add new condition name and initialize its count
            bc_on_skin_projection_type.push_back(condition_layer_name);
            count_bc_on_skin_projection_type.push_back(1); // Initialize count to 1
        }
        count_list_closest_condition++;
    }
    //--------------------
    std::string max_layer_condition_name;
    std::string max_condition_name;
    std::string layer_name;
    SizeType max_count = 0;
    for (size_t i = 0; i < count_bc_on_skin_projection_type.size(); ++i) {
        if (count_bc_on_skin_projection_type[i] > max_count) {
            max_count = count_bc_on_skin_projection_type[i];
            max_layer_condition_name = bc_on_skin_projection_type[i];
        }
    }
    // create a pool of conditions of the right type
    std::vector<int> list_id_closest_condition_of_correct_bc;
    for (IndexType i = 0; i < rListIdClosestCondition.size(); i++) {
        IndexType condId = rListIdClosestCondition[i];
        std::string condition_layer_name = rSkinModelPart.GetCondition(rListIdClosestCondition[i]).GetValue(LAYER_NAME);
        if (condition_layer_name == max_layer_condition_name) 
        {
            list_id_closest_condition_of_correct_bc.push_back(condId);
            layer_name = max_layer_condition_name;
            max_condition_name = rSkinModelPart.GetCondition(rListIdClosestCondition[i]).GetValue(CONDITION_NAME);
        }
    }

    KRATOS_ERROR_IF(list_id_closest_condition_of_correct_bc.size() != max_count) << "ERROR in list_id_closest_condition_of_correct_bc" << std::endl;
    // correct the projections
    count_list_closest_condition = 0;
    for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it) {

        std::string condition_layer_name = rSkinModelPart.GetCondition(rListIdClosestCondition[count_list_closest_condition]).GetValue(LAYER_NAME);

        auto gp_coord = (*it)->Center();
        int best_cond_id = -1;
        if (condition_layer_name != max_layer_condition_name) 
        {
            bool use_projections_of_the_others_quadrature_points = false;
            // search for the second closest condition
            if (rListIdSecondClosestCondition[count_list_closest_condition] != -1)
            {
                // if there is a second closest condition, we need to check if it is of the correct layer
                rListIdClosestCondition[count_list_closest_condition] = rListIdSecondClosestCondition[count_list_closest_condition];
                std::string second_closest_condition_layer_name = rSkinModelPart.GetCondition(rListIdSecondClosestCondition[count_list_closest_condition]).GetValue(LAYER_NAME);
                
                if (second_closest_condition_layer_name != max_layer_condition_name)
                    use_projections_of_the_others_quadrature_points = true;
            } else {
                use_projections_of_the_others_quadrature_points = true;
            }

            // if there is no second closest condition, we need to find the closest condition of the correct layer
            if (use_projections_of_the_others_quadrature_points)
            {
                double best_distance = 1e16;
                for (IndexType i = 0; i < max_count; i++)
                {
                    int cond_id = list_id_closest_condition_of_correct_bc[i];
                    auto cond_center = (&rSkinModelPart.GetCondition(cond_id))->GetGeometry().Center();
                    double curr_distance = norm_2(cond_center-gp_coord);
                    if (curr_distance < best_distance) 
                    {
                        best_distance = curr_distance;
                        best_cond_id = cond_id;
                    }
                }
                rListIdClosestCondition[count_list_closest_condition] = best_cond_id;
            }
        }
        count_list_closest_condition++;
    }

    KRATOS_ERROR_IF(!KratosComponents<Condition>::Has(max_condition_name))
            << max_condition_name << " not registered." << std::endl;

    const Condition& rReferenceCondition = KratosComponents<Condition>::Get(max_condition_name);
    count_list_closest_condition = 0;

    ModelPart& r_layer_model_part = rModelPart.HasSubModelPart(layer_name) ? 
                                    rModelPart.GetSubModelPart(layer_name) : 
                                    rModelPart.CreateSubModelPart(layer_name);
    // 2D case
    if (rSkinModelPart.ConditionsBegin()->GetGeometry().size() == 2) {

        for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it) {
            new_condition_list.push_back(
                rReferenceCondition.Create(rIdCounter, (*it), pProperties));

            IndexType condId = rListIdClosestCondition[count_list_closest_condition];

            Condition::Pointer cond1 = &rSkinModelPart.GetCondition(condId);
            IndexType condId2;  
            if (condId == rSkinModelPart.ConditionsBegin()->Id()) {
                condId2 = (rSkinModelPart.ConditionsEnd()-1)->Id();
            }
            else condId2 = condId-1;
            Condition::Pointer cond2 = &rSkinModelPart.GetCondition(condId2);

            // Add closest projection node
            NodePointerVector empty_vector;
            PointTypePointer projection_node = cond1->GetGeometry()(0);
            if (norm_2(projection_node->GetValue(NORMAL)) > 1e-13)
            {
                empty_vector.push_back(cond1->GetGeometry()(0)); // Just it_node-plane neighbours
                (*it)->SetValue(NEIGHBOUR_NODES, empty_vector);
            }
            
            new_condition_list.GetContainer()[count_list_closest_condition]->SetValue(NEIGHBOUR_CONDITIONS, GlobalPointersVector<Condition>({cond1,cond2}));
            if (IsInner) {
                new_condition_list.GetContainer()[count_list_closest_condition]->SetValue(IDENTIFIER, "inner");
            } else {
                new_condition_list.GetContainer()[count_list_closest_condition]->SetValue(IDENTIFIER, "outer");
            }
            new_condition_list.GetContainer()[count_list_closest_condition]->SetValue(KNOT_SPAN_SIZES, KnotSpanSizes);                            

            rIdCounter++;
            count_list_closest_condition++;
        }
    } else {
        // TODO: 3D case
        KRATOS_ERROR << "CreateConditions: 3D case not implemented yet." << std::endl;
    }
    
    r_layer_model_part.AddConditions(new_condition_list.begin(), new_condition_list.end());
}

void IgaModelerSbm::ActivateNodesInElementsAndCleanRoot(ModelPart& rAnalysisModelPart) const
{
    for (auto& r_node : rAnalysisModelPart.Nodes()) {
        r_node.Set(ACTIVE, false);
    }

    for (auto& r_elem : rAnalysisModelPart.Elements()) {
        auto& r_geom = r_elem.GetGeometry();
        for (auto& r_node : r_geom) {
            r_node.Set(ACTIVE, true);
        }
    }

    ModelPart& r_root = rAnalysisModelPart.GetRootModelPart();
    std::vector<ModelPart::IndexType> node_ids_to_remove;
    node_ids_to_remove.reserve(r_root.NumberOfNodes());

    for (auto& r_node : r_root.Nodes()) {
        if (r_node.IsDefined(ACTIVE) && r_node.IsNot(ACTIVE)) {
            node_ids_to_remove.push_back(r_node.Id());
        }
    }

    for (const auto node_id : node_ids_to_remove) {
        r_root.RemoveNode(node_id);
    }
}

///@}
}
