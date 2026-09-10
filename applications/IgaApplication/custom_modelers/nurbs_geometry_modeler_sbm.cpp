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

// Project includes
#include "nurbs_geometry_modeler_sbm.h"
#include "custom_utilities/create_breps_sbm_utilities.h"
#include "custom_processes/snake_sbm_process.h"
#include "iga_application_variables.h"

// System includes
#include <algorithm>
#include <limits>

namespace Kratos
{

namespace
{

using BrepCurveOnSurfaceType =
    BrepCurveOnSurface<PointerVector<Node>, true, PointerVector<Point>>;

std::string MapLocalRefinementLayerName(const std::string& rSkinLayerName)
{
    if (rSkinLayerName == "COUPLING_SIDE_OUTER") {
        return "COUPLING_CONDITION_OUTER";
    }
    if (rSkinLayerName == "COUPLING_SIDE_INNER") {
        return "COUPLING_CONDITION_INNER";
    }
    return rSkinLayerName;
}

double SquaredDistanceToSegment(
    const array_1d<double, 3>& rPoint,
    const array_1d<double, 3>& rSegmentStart,
    const array_1d<double, 3>& rSegmentEnd)
{
    const array_1d<double, 3> segment = rSegmentEnd - rSegmentStart;
    const double length_squared = inner_prod(segment, segment);
    if (length_squared <= std::numeric_limits<double>::epsilon()) {
        const array_1d<double, 3> delta = rPoint - rSegmentStart;
        return inner_prod(delta, delta);
    }

    const array_1d<double, 3> point_delta = rPoint - rSegmentStart;
    const double projection = std::clamp(
        inner_prod(point_delta, segment) / length_squared,
        0.0,
        1.0);
    const array_1d<double, 3> closest_point = rSegmentStart + projection * segment;
    const array_1d<double, 3> delta = rPoint - closest_point;
    return inner_prod(delta, delta);
}

void AddLayerConditionMetadata(
    Node& rNode,
    const std::string& rLayerName,
    const std::string& rConditionName)
{
    auto connected_layers = rNode.GetValue(CONNECTED_LAYERS);
    auto connected_conditions = rNode.GetValue(CONNECTED_CONDITIONS);
    if (connected_conditions.size() < connected_layers.size()) {
        connected_conditions.resize(connected_layers.size());
    }

    const auto layer_it = std::find(
        connected_layers.begin(), connected_layers.end(), rLayerName);
    if (layer_it == connected_layers.end()) {
        connected_layers.push_back(rLayerName);
        connected_conditions.push_back(rConditionName);
    } else {
        const auto layer_index = static_cast<std::size_t>(
            std::distance(connected_layers.begin(), layer_it));
        if (connected_conditions[layer_index].empty()) {
            connected_conditions[layer_index] = rConditionName;
        }
    }

    rNode.SetValue(CONNECTED_LAYERS, connected_layers);
    rNode.SetValue(CONNECTED_CONDITIONS, connected_conditions);
}

template<class TDataContainerType>
void AddNeighbourGeometry(
    TDataContainerType& rContainer,
    const Geometry<Node>::Pointer& pGeometry)
{
    auto neighbour_geometries = rContainer.GetValue(NEIGHBOUR_GEOMETRIES);
    const bool already_present = std::any_of(
        neighbour_geometries.begin(),
        neighbour_geometries.end(),
        [&](const auto& pExistingGeometry) {
            return pExistingGeometry.get() == pGeometry.get();
        });
    if (!already_present) {
        neighbour_geometries.push_back(pGeometry);
        rContainer.SetValue(NEIGHBOUR_GEOMETRIES, neighbour_geometries);
    }
}

void PrepareLocalRefinementSurrogateData(
    ModelPart& rIgaModelPart,
    ModelPart& rSkinModelPart,
    const NurbsSurfaceGeometry<3, PointerVector<Node>>::Pointer& pSurface)
{
    KRATOS_ERROR_IF_NOT(pSurface)
        << "NurbsGeometryModelerSbm: local-refinement preparation requires a valid NURBS surface."
        << std::endl;

    auto prepare_loop = [&](const std::string& rLoopName, const std::string& rSurrogateName) {
        if (!rSkinModelPart.HasSubModelPart(rLoopName) ||
            !rIgaModelPart.HasSubModelPart(rSurrogateName)) {
            return;
        }

        ModelPart& r_skin_loop = rSkinModelPart.GetSubModelPart(rLoopName);
        ModelPart& r_surrogate_loop = rIgaModelPart.GetSubModelPart(rSurrogateName);
        if (r_skin_loop.NumberOfConditions() == 0) {
            return;
        }

        const std::size_t max_surface_degree = std::max(
            pSurface->PolynomialDegree(0), pSurface->PolynomialDegree(1));
        const std::size_t shape_function_derivatives_order =
            2 * max_surface_degree + 1;

        for (auto& r_surrogate_condition : r_surrogate_loop.Conditions()) {
            const auto surrogate_center = r_surrogate_condition.GetGeometry().Center();
            const Condition* p_closest_skin_condition = nullptr;
            double closest_distance_squared = std::numeric_limits<double>::max();

            for (const auto& r_skin_condition : r_skin_loop.Conditions()) {
                const auto& r_skin_geometry = r_skin_condition.GetGeometry();
                if (r_skin_geometry.PointsNumber() < 2 || !r_skin_condition.Has(LAYER_NAME)) {
                    continue;
                }

                const double distance_squared = SquaredDistanceToSegment(
                    surrogate_center.Coordinates(),
                    r_skin_geometry[0].Coordinates(),
                    r_skin_geometry[1].Coordinates());
                if (distance_squared < closest_distance_squared) {
                    closest_distance_squared = distance_squared;
                    p_closest_skin_condition = &r_skin_condition;
                }
            }

            KRATOS_ERROR_IF_NOT(p_closest_skin_condition)
                << "NurbsGeometryModelerSbm: could not associate surrogate condition #"
                << r_surrogate_condition.Id() << " in '" << r_surrogate_loop.FullName()
                << "' with a skin layer." << std::endl;

            const std::string skin_layer_name =
                p_closest_skin_condition->GetValue(LAYER_NAME);
            const std::string surrogate_layer_name =
                MapLocalRefinementLayerName(skin_layer_name);
            const std::string condition_name =
                p_closest_skin_condition->Has(CONDITION_NAME)
                    ? p_closest_skin_condition->GetValue(CONDITION_NAME)
                    : "";
            r_surrogate_condition.SetValue(LAYER_NAME, surrogate_layer_name);

            ModelPart& r_layer_model_part =
                rIgaModelPart.HasSubModelPart(skin_layer_name)
                    ? rIgaModelPart.GetSubModelPart(skin_layer_name)
                    : rIgaModelPart.CreateSubModelPart(skin_layer_name);
            r_layer_model_part.AddCondition(
                r_surrogate_loop.pGetCondition(r_surrogate_condition.Id()));

            auto& r_surrogate_geometry = r_surrogate_condition.GetGeometry();
            for (IndexType i = 0; i < r_surrogate_geometry.PointsNumber(); ++i) {
                r_layer_model_part.AddNode(r_surrogate_geometry.pGetPoint(i));
                AddLayerConditionMetadata(
                    r_surrogate_geometry[i], skin_layer_name, condition_name);
            }

            KRATOS_ERROR_IF_NOT(r_surrogate_condition.Has(BREP_ID))
                << "NurbsGeometryModelerSbm: surrogate condition #"
                << r_surrogate_condition.Id() << " in '" << r_surrogate_loop.FullName()
                << "' has no BREP_ID after CreateSurrogateBoundary." << std::endl;

            const IndexType brep_id = static_cast<IndexType>(
                r_surrogate_condition.GetValue(BREP_ID));
            if (r_surrogate_condition.Has(BREP_MODEL_PART_FULL_NAME)) {
                const std::string& r_brep_model_part_name =
                    r_surrogate_condition.GetValue(BREP_MODEL_PART_FULL_NAME);
                KRATOS_ERROR_IF(
                    !r_brep_model_part_name.empty() &&
                    r_brep_model_part_name != rIgaModelPart.FullName())
                    << "NurbsGeometryModelerSbm: surrogate condition #"
                    << r_surrogate_condition.Id() << " references BREP model part '"
                    << r_brep_model_part_name << "', expected '"
                    << rIgaModelPart.FullName() << "'." << std::endl;
            }
            KRATOS_ERROR_IF_NOT(rIgaModelPart.HasGeometry(brep_id))
                << "NurbsGeometryModelerSbm: BREP geometry #" << brep_id
                << " referenced by surrogate condition #" << r_surrogate_condition.Id()
                << " was not found in '" << rIgaModelPart.FullName() << "'." << std::endl;

            auto p_brep_geometry = rIgaModelPart.pGetGeometry(brep_id);
            r_layer_model_part.AddGeometry(p_brep_geometry);
            auto p_brep_curve =
                std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_brep_geometry);
            KRATOS_ERROR_IF_NOT(p_brep_curve)
                << "NurbsGeometryModelerSbm: local-refinement BREP #" << brep_id
                << " is not a BrepCurveOnSurface." << std::endl;

            const NurbsInterval domain_interval = p_brep_curve->DomainInterval();
            IntegrationPoint<1> center_integration_point(
                0.5 * (domain_interval.GetT0() + domain_interval.GetT1()));
            Geometry<Node>::IntegrationPointsArrayType integration_points;
            integration_points.push_back(center_integration_point);
            Geometry<Node>::GeometriesArrayType quadrature_points;
            IntegrationInfo integration_info =
                p_brep_curve->GetDefaultIntegrationInfo();
            p_brep_curve->CreateQuadraturePointGeometries(
                quadrature_points,
                shape_function_derivatives_order,
                integration_points,
                integration_info);

            KRATOS_ERROR_IF(quadrature_points.size() != 1)
                << "NurbsGeometryModelerSbm: expected one central quadrature geometry for BREP #"
                << brep_id << ", got " << quadrature_points.size() << "." << std::endl;

            const auto p_central_geometry = quadrature_points(0);
            AddNeighbourGeometry(*p_brep_geometry, p_central_geometry);
            AddNeighbourGeometry(r_surrogate_condition, p_central_geometry);
            for (IndexType i = 0; i < r_surrogate_geometry.PointsNumber(); ++i) {
                AddNeighbourGeometry(r_surrogate_geometry[i], p_central_geometry);
            }
        }
    };

    prepare_loop("inner", "surrogate_inner");
    prepare_loop("outer", "surrogate_outer");
}

} // namespace

///@name Stages
///@{

///@}
///@name Private Operations
///@{
void NurbsGeometryModelerSbm::CreateAndAddRegularGrid2D(
    ModelPart& rModelPart, 
    const Point& A_xyz, 
    const Point& B_xyz,
    const Point& A_uvw, 
    const Point& B_uvw, 
    const SizeType OrderU, 
    const SizeType OrderV, 
    const SizeType NumKnotSpansU, 
    const SizeType NumKnotSpansV, 
    const bool AddSurfaceToModelPart)
{   

    // Call the CreateAndAddRegularGrid2D method of the base class NurbsGeometryModeler
    NurbsGeometryModeler::CreateAndAddRegularGrid2D(rModelPart, A_xyz, B_xyz,
        A_uvw, B_uvw, OrderU, OrderV, NumKnotSpansU, NumKnotSpansV, false);
        
    // Create the Domain/Iga Model Part
    const std::string iga_model_part_name = mParameters["model_part_name"].GetString();
    ModelPart& r_iga_model_part = mpModel->HasModelPart(iga_model_part_name)
                                ? mpModel->GetModelPart(iga_model_part_name)
                                : mpModel->CreateModelPart(iga_model_part_name);

    // compute unique_knot_vector_u
    Vector unique_knot_vector_u(2+(NumKnotSpansU-1));
    unique_knot_vector_u[0] = mKnotVectorU[0]; 
    unique_knot_vector_u[NumKnotSpansU] = mKnotVectorU[mKnotVectorU.size()-1];
    for (SizeType i_knot_insertion = 0; i_knot_insertion < NumKnotSpansU-1; i_knot_insertion++) {
        unique_knot_vector_u[i_knot_insertion+1] = mInsertKnotsU[i_knot_insertion];
    }
    // compute unique_knot_vector_v
    Vector unique_knot_vector_v(2+(NumKnotSpansV-1));
    unique_knot_vector_v[0] = mKnotVectorV[0]; 
    unique_knot_vector_v[NumKnotSpansV] = mKnotVectorV[mKnotVectorV.size()-1];
    for (SizeType i_knot_insertion = 0; i_knot_insertion < NumKnotSpansV-1; i_knot_insertion++) {
        unique_knot_vector_v[i_knot_insertion+1] = mInsertKnotsV[i_knot_insertion];
    }
    // Set the value of the knot vectors
    r_iga_model_part.SetValue(KNOT_VECTOR_U, unique_knot_vector_u);
    r_iga_model_part.SetValue(KNOT_VECTOR_V, unique_knot_vector_v);

    // If neither skin_inner nor skin_outer exists
    if (!(mParameters.Has("skin_model_part_inner_initial_name") || mParameters.Has("skin_model_part_outer_initial_name"))){
        
        Vector knot_step_uv= ZeroVector(2);
        knot_step_uv[0] = std::abs(unique_knot_vector_u[std::ceil(unique_knot_vector_u.size()/2) +1] - unique_knot_vector_u[std::ceil(unique_knot_vector_u.size()/2)] ) ;
        knot_step_uv[1] = std::abs(unique_knot_vector_v[std::ceil(unique_knot_vector_v.size()/2) +1] - unique_knot_vector_v[std::ceil(unique_knot_vector_v.size()/2)] ) ;

        // saving the knot span sizes
        r_iga_model_part.SetValue(KNOT_SPAN_SIZES, knot_step_uv);

        // Create the breps for the outer sbm boundary
        CreateBrepsSbmUtilities<Node, Point, false> CreateBrepsSbmUtilities(mEchoLevel);
        CreateBrepsSbmUtilities.CreateSurrogateBoundary(mpSurface, A_uvw, B_uvw, rModelPart);

        //TODO: This must be turned to an error once we finish the ongoing SBM BCs development
        KRATOS_WARNING("None of the 'skin_model_part_name' have not been defined ") << 
                        "in the nurbs_geometry_modeler_sbm in the project paramer json" << std::endl;
        return;
    }

    // Create the True Model Part -> contains all the true boundary features
    std::string skin_model_part_name;

    // Retrieve skin_model_part_inner_initial_name if it exists
    std::string skin_model_part_inner_initial_name = "skin_model_part_outer_initial_name";
    if (mParameters.Has("skin_model_part_inner_initial_name")) {
        skin_model_part_inner_initial_name = mParameters["skin_model_part_inner_initial_name"].GetString();
    }

    // Retrieve skin_model_part_outer_initial_name if it exists;
    std::string skin_model_part_outer_initial_name = "skin_model_part_outer_initial_name";
    if (mParameters.Has("skin_model_part_outer_initial_name")) {
        skin_model_part_outer_initial_name = mParameters["skin_model_part_outer_initial_name"].GetString();
    }

    // Create the surrogate sub model parts inner and outer
    ModelPart& surrogate_sub_model_part_inner = r_iga_model_part.HasSubModelPart("surrogate_inner")
        ? r_iga_model_part.GetSubModelPart("surrogate_inner")
        : r_iga_model_part.CreateSubModelPart("surrogate_inner");

    ModelPart& surrogate_sub_model_part_outer = r_iga_model_part.HasSubModelPart("surrogate_outer")
        ? r_iga_model_part.GetSubModelPart("surrogate_outer")
        : r_iga_model_part.CreateSubModelPart("surrogate_outer");

    if (mParameters.Has("skin_model_part_name"))
        skin_model_part_name = mParameters["skin_model_part_name"].GetString();
    else
        KRATOS_ERROR << "The skin_model_part name '" << skin_model_part_name << "' was not defined in the project parameters.\n" << std::endl;

    // inner
    mpModel->HasModelPart(skin_model_part_inner_initial_name)
        ? mpModel->GetModelPart(skin_model_part_inner_initial_name)
        : mpModel->CreateModelPart(skin_model_part_inner_initial_name);
    // outer
    mpModel->HasModelPart(skin_model_part_outer_initial_name)
        ? mpModel->GetModelPart(skin_model_part_outer_initial_name)
        : mpModel->CreateModelPart(skin_model_part_outer_initial_name);
    
    // Skin model part refined after Snake Process — get or create
    ModelPart& skin_model_part = mpModel->HasModelPart(skin_model_part_name)
        ? mpModel->GetModelPart(skin_model_part_name)
        : mpModel->CreateModelPart(skin_model_part_name);

    // Ensure "inner" submodel part exists
    ModelPart& skin_inner = skin_model_part.HasSubModelPart("inner")
        ? skin_model_part.GetSubModelPart("inner")
        : skin_model_part.CreateSubModelPart("inner");

    // Ensure "outer" submodel part exists
    ModelPart& skin_outer = skin_model_part.HasSubModelPart("outer")
        ? skin_model_part.GetSubModelPart("outer")
        : skin_model_part.CreateSubModelPart("outer");


    // Create the parameters for the SnakeSbmProcess
    Kratos::Parameters snake_parameters;
    snake_parameters.AddString("model_part_name", iga_model_part_name);
    snake_parameters.AddString("skin_model_part_name", skin_model_part_name);
    snake_parameters.AddDouble("echo_level", mEchoLevel);
    snake_parameters.AddString("skin_model_part_inner_initial_name", skin_model_part_inner_initial_name);
    snake_parameters.AddString("skin_model_part_outer_initial_name", skin_model_part_outer_initial_name);
    if (mParameters.Has("lambda_inner"))
        snake_parameters.AddDouble("lambda_inner", mParameters["lambda_inner"].GetDouble());
    if (mParameters.Has("lambda_outer"))
        snake_parameters.AddDouble("lambda_outer", mParameters["lambda_outer"].GetDouble());
    if (mParameters.Has("number_of_inner_loops"))
        snake_parameters.AddDouble("number_of_inner_loops", mParameters["number_of_inner_loops"].GetInt());
    if (mParameters.Has("number_initial_points_if_importing_nurbs"))
        snake_parameters.AddInt("number_initial_points_if_importing_nurbs", mParameters["number_initial_points_if_importing_nurbs"].GetInt());
    if (mParameters.Has("create_surr_outer_from_surr_inner"))
        snake_parameters.AddBool("create_surr_outer_from_surr_inner", mParameters["create_surr_outer_from_surr_inner"].GetBool());
    if (mParameters.Has("create_surr_inner_from_surr_outer"))
        snake_parameters.AddBool("create_surr_inner_from_surr_outer", mParameters["create_surr_inner_from_surr_outer"].GetBool());

    // Create the surrogate_sub_model_part for inner and outer
    SnakeSbmProcess snake_sbm_process(*mpModel, snake_parameters);
    snake_sbm_process.Execute();

    // Create the breps for the outer sbm boundary
    CreateBrepsSbmUtilities<Node, Point, true> CreateBrepsSbmUtilities(mEchoLevel);
    CreateBrepsSbmUtilities.CreateSurrogateBoundary(mpSurface, surrogate_sub_model_part_inner, surrogate_sub_model_part_outer, A_uvw, B_uvw, r_iga_model_part);

    if (mParameters["use_for_local_refinement"].GetBool()) {
        PrepareLocalRefinementSurrogateData(
            r_iga_model_part,
            skin_model_part,
            mpSurface);
    }


}

// 3D 
    void NurbsGeometryModelerSbm::CreateAndAddRegularGrid3D( 
        ModelPart& rModelPart,
        const Point& A_xyz,
        const Point& B_xyz,
        const Point& A_uvw,
        const Point& B_uvw,
        const SizeType OrderU,
        const SizeType OrderV,
        const SizeType OrderW,
        const SizeType NumKnotSpansU,
        const SizeType NumKnotSpansV,
        const SizeType NumKnotSpansW,
        const bool AddVolumeToModelPart)
    {   

        // Call the CreateAndAddRegularGrid3D method of the base class NurbsGeometryModeler
        NurbsGeometryModeler::CreateAndAddRegularGrid3D(rModelPart, A_xyz, B_xyz,
            A_uvw, B_uvw, OrderU, OrderV, OrderW, NumKnotSpansU, NumKnotSpansV, NumKnotSpansW, false);
                 
        // Create the Domain/Iga Model Part
        const std::string iga_model_part_name = mParameters["model_part_name"].GetString();
        ModelPart& iga_model_part = mpModel->HasModelPart(iga_model_part_name)
                                    ? mpModel->GetModelPart(iga_model_part_name)
                                    : mpModel->CreateModelPart(iga_model_part_name);

        // Create the True Model Part -> contains all the true boundary features
        std::string skin_model_part_inner_initial_name = "SkinModelPartInnerInitial";
        std::string skin_model_part_outer_initial_name = "SkinModelPartOuterInitial";
        std::string skin_model_part_name;
        if (mParameters.Has("skin_model_part_inner_initial_name")) {
            skin_model_part_inner_initial_name = mParameters["skin_model_part_inner_initial_name"].GetString();

            KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(skin_model_part_inner_initial_name)) 
                         << "The skin_model_part '" << skin_model_part_inner_initial_name << "' was not created in the model.\n" 
                         << "Check the reading of the mdpa file in the import mdpa modeler."<< std::endl;
        }
        if (mParameters.Has("skin_model_part_outer_initial_name")) {
            skin_model_part_outer_initial_name = mParameters["skin_model_part_outer_initial_name"].GetString();

            KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(skin_model_part_outer_initial_name)) 
                         << "The skin_model_part '" << skin_model_part_outer_initial_name << "' was not created in the model.\n" 
                         << "Check the reading of the mdpa file in the import mdpa modeler."<< std::endl;
        }

        // Create the surrogate sub model parts inner and outer
        ModelPart& surrogate_sub_model_part_inner = iga_model_part.CreateSubModelPart("surrogate_inner");
        ModelPart& surrogate_sub_model_part_outer = iga_model_part.CreateSubModelPart("surrogate_outer");

        // compute unique_knot_vector_u
        Vector unique_knot_vector_u(2+(NumKnotSpansU-1));
        unique_knot_vector_u[0] = mKnotVectorU[0]; unique_knot_vector_u[NumKnotSpansU] = mKnotVectorU[mKnotVectorU.size()-1];
        for (SizeType i_knot_insertion = 0; i_knot_insertion < NumKnotSpansU-1; i_knot_insertion++) {
            unique_knot_vector_u[i_knot_insertion+1] = mInsertKnotsU[i_knot_insertion];
        }

        // compute unique_knot_vector_v
        Vector unique_knot_vector_v(2+(NumKnotSpansV-1));
        unique_knot_vector_v[0] = mKnotVectorV[0]; unique_knot_vector_v[NumKnotSpansV] = mKnotVectorV[mKnotVectorV.size()-1];
        for (SizeType i_knot_insertion = 0; i_knot_insertion < NumKnotSpansV-1; i_knot_insertion++) {
            unique_knot_vector_v[i_knot_insertion+1] = mInsertKnotsV[i_knot_insertion];
        }

        // compute unique_knot_vector_w
        Vector unique_knot_vector_w(2+(NumKnotSpansW-1));
        unique_knot_vector_w[0] = mKnotVectorW[0]; unique_knot_vector_w[NumKnotSpansW] = mKnotVectorW[mKnotVectorW.size()-1];
        for (SizeType i_knot_insertion = 0; i_knot_insertion < NumKnotSpansW-1; i_knot_insertion++) {
            unique_knot_vector_w[i_knot_insertion+1] = mInsertKnotsW[i_knot_insertion];
        }

        // Set the value of the knot vectors
        iga_model_part.SetValue(KNOT_VECTOR_U, unique_knot_vector_u);
        iga_model_part.SetValue(KNOT_VECTOR_V, unique_knot_vector_v);
        iga_model_part.SetValue(KNOT_VECTOR_W, unique_knot_vector_w);

        // Save knot span sizes for 3D.
        Vector knot_step_uvw = ZeroVector(3);
        const SizeType mid_u = static_cast<SizeType>(std::ceil(unique_knot_vector_u.size() / 2.0));
        const SizeType mid_v = static_cast<SizeType>(std::ceil(unique_knot_vector_v.size() / 2.0));
        const SizeType mid_w = static_cast<SizeType>(std::ceil(unique_knot_vector_w.size() / 2.0));
        knot_step_uvw[0] = std::abs(unique_knot_vector_u[mid_u + 1] - unique_knot_vector_u[mid_u]);
        knot_step_uvw[1] = std::abs(unique_knot_vector_v[mid_v + 1] - unique_knot_vector_v[mid_v]);
        knot_step_uvw[2] = std::abs(unique_knot_vector_w[mid_w + 1] - unique_knot_vector_w[mid_w]);
        iga_model_part.SetValue(KNOT_SPAN_SIZES, knot_step_uvw);

        // If there is not neither skin_inner nor skin_outer throw an error since you are using the sbm modeler
        if (!(mParameters.Has("skin_model_part_inner_initial_name") || mParameters.Has("skin_model_part_outer_initial_name"))){
            
            // Create the breps for the outer sbm boundary
            CreateBrepsSbmUtilities<Node, Point> CreateBrepsSbmUtilities(mEchoLevel);
            CreateBrepsSbmUtilities.CreateSurrogateBoundary(mpVolume, A_uvw, B_uvw, rModelPart);

            KRATOS_WARNING("None of the 'skin_model_part_name' have not been defined ") << 
                            "in the nurbs_geometry_modeler_sbm in the project paramer json" << std::endl;
            
            return;
        }
        
        if (mParameters.Has("skin_model_part_name"))
            skin_model_part_name = mParameters["skin_model_part_name"].GetString();
        else
            KRATOS_ERROR << "The skin_model_part name '" << skin_model_part_name << "' was not defined in the project parameters.\n" << std::endl;
 
        // inner
        mpModel->HasModelPart(skin_model_part_inner_initial_name)
            ? mpModel->GetModelPart(skin_model_part_inner_initial_name)
            : mpModel->CreateModelPart(skin_model_part_inner_initial_name);
        // outer
        mpModel->HasModelPart(skin_model_part_outer_initial_name)
            ? mpModel->GetModelPart(skin_model_part_outer_initial_name)
            : mpModel->CreateModelPart(skin_model_part_outer_initial_name);
        
        // Skin model part refined after Snake Process
        ModelPart& skin_model_part = mpModel->CreateModelPart(skin_model_part_name);
        skin_model_part.CreateSubModelPart("inner");
        skin_model_part.CreateSubModelPart("outer");
        
        
        // Create the parameters for the SnakeSbmProcess
        Kratos::Parameters snake_parameters;
        snake_parameters.AddString("model_part_name", iga_model_part_name);
        snake_parameters.AddString("skin_model_part_name", skin_model_part_name);
        snake_parameters.AddDouble("echo_level", mEchoLevel);
        snake_parameters.AddString("skin_model_part_inner_initial_name", skin_model_part_inner_initial_name);
        snake_parameters.AddString("skin_model_part_outer_initial_name", skin_model_part_outer_initial_name);
        if (mParameters.Has("lambda_inner"))
            snake_parameters.AddDouble("lambda_inner", mParameters["lambda_inner"].GetDouble());
        if (mParameters.Has("lambda_outer"))
            snake_parameters.AddDouble("lambda_outer", mParameters["lambda_outer"].GetDouble());
        if (mParameters.Has("number_of_inner_loops"))
            snake_parameters.AddDouble("number_of_inner_loops", mParameters["number_of_inner_loops"].GetInt());
        
        // Create the surrogate_sub_model_part for inner and outer // TODO: extend this in 3D
        SnakeSbmProcess snake_sbm_process(*mpModel, snake_parameters);
        snake_sbm_process.Execute();

        // Create the breps for the outer sbm boundary // TODO: extend this in 3D
        CreateBrepsSbmUtilities<Node, Point> CreateBrepsSbmUtilities(mEchoLevel);
        CreateBrepsSbmUtilities.CreateSurrogateBoundary(mpVolume, surrogate_sub_model_part_inner, surrogate_sub_model_part_outer, A_uvw, B_uvw, iga_model_part);
    }


const Parameters NurbsGeometryModelerSbm::GetDefaultParameters() const
{
    return Parameters(R"(
    {
        "echo_level": 1,
        "model_part_name" : "IgaModelPart",
        "lower_point_xyz": [0.0, 0.0, 0.0],
        "upper_point_xyz": [1.0, 1.0, 0.0],
        "lower_point_uvw": [0.0, 0.0, 0.0],
        "upper_point_uvw": [1.0, 1.0, 0.0],
        "polynomial_order" : [2, 2],
        "number_of_knot_spans" : [10, 10],
        "create_surr_outer_from_surr_inner": false,
        "create_surr_inner_from_surr_outer": false,
        "lambda_inner": 0.5,
        "lambda_outer": 0.5,
        "number_of_inner_loops": 0,
        "use_for_local_refinement": false
    })");
}

const Parameters NurbsGeometryModelerSbm::GetValidParameters() const
{
    return Parameters(R"(
    {
        "echo_level": 1,
        "model_part_name" : "IgaModelPart",
        "lower_point_xyz": [0.0, 0.0, 0.0],
        "upper_point_xyz": [1.0, 1.0, 0.0],
        "lower_point_uvw": [0.0, 0.0, 0.0],
        "upper_point_uvw": [1.0, 1.0, 0.0],
        "polynomial_order" : [2, 2],
        "number_of_knot_spans" : [10, 10],
        "create_surr_outer_from_surr_inner": false,
        "create_surr_inner_from_surr_outer": false,
        "lambda_inner": 0.5,
        "lambda_outer": 0.5,
        "number_of_inner_loops": 0,
        "number_initial_points_if_importing_nurbs": 1,
        "skin_model_part_inner_initial_name": "skin_model_part_inner_initial",
        "skin_model_part_outer_initial_name": "skin_model_part_outer_initial",
        "skin_model_part_name": "skin_model_part",
        "use_for_local_refinement": false
    })");
}

} // end namespace kratos
