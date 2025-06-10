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

namespace Kratos
{

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
    ModelPart& surrogate_sub_model_part_inner = r_iga_model_part.CreateSubModelPart("surrogate_inner");
    ModelPart& surrogate_sub_model_part_outer = r_iga_model_part.CreateSubModelPart("surrogate_outer");
    
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

    // Create the surrogate_sub_model_part for inner and outer
    SnakeSbmProcess snake_sbm_process(*mpModel, snake_parameters);
    snake_sbm_process.Execute();

    // Create the breps for the outer sbm boundary
    CreateBrepsSbmUtilities<Node, Point, true> CreateBrepsSbmUtilities(mEchoLevel);
    CreateBrepsSbmUtilities.CreateSurrogateBoundary(mpSurface, surrogate_sub_model_part_inner, surrogate_sub_model_part_outer, A_uvw, B_uvw, r_iga_model_part);
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
        "lambda_inner": 0.5,
        "lambda_outer": 0.5,
        "number_of_inner_loops": 0
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
        "lambda_inner": 0.5,
        "lambda_outer": 0.5,
        "number_of_inner_loops": 0,
        "skin_model_part_inner_initial_name": "skin_model_part_inner_initial",
        "skin_model_part_outer_initial_name": "skin_model_part_outer_initial",
        "skin_model_part_name": "skin_model_part"
    })");
}

} // end namespace kratos
