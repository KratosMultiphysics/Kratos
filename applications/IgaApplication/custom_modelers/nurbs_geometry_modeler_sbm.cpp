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
#include "includes/define.h"
#include "nurbs_geometry_modeler_sbm.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities.h"
#include "custom_utilities/create_breps_sbm_utilities.h"

namespace Kratos
{

    ///@name Stages
    ///@{

    void NurbsGeometryModelerSbm::SetupGeometryModel(){

        //----------------------------------------------------------------------------------------------------------------
        KRATOS_INFO_IF("NurbsGeometryModelerSbm", mEchoLevel > 1) << "[NURBS MODELER SBM]:: calling NurbsGeometryModeler" << std::endl;
        
        // Call the SetupGeometryModel method of the base class NurbsGeometryModeler
        NurbsGeometryModeler::SetupGeometryModel();
    }

    ///@}
    ///@name Private Operations
    ///@{
    void NurbsGeometryModelerSbm::CreateAndAddRegularGrid2D( ModelPart& r_model_part, const Point& A_xyz, const Point& B_xyz,
        const Point& A_uvw, const Point& B_uvw, SizeType OrderU, SizeType OrderV,SizeType NumKnotSpansU, SizeType NumKnotSpansV, bool add_surface_to_model_part)
    {   

        // Call the CreateAndAddRegularGrid2D method of the base class NurbsGeometryModeler
        NurbsGeometryModeler::CreateAndAddRegularGrid2D(r_model_part, A_xyz, B_xyz,
            A_uvw, B_uvw, OrderU, OrderV, NumKnotSpansU, NumKnotSpansV, false);
         
        // Create the Domain/Iga Model Part
        const std::string iga_model_part_name = mParameters["model_part_name"].GetString();
        ModelPart& iga_model_part = mpModel->HasModelPart(iga_model_part_name)
                                    ? mpModel->GetModelPart(iga_model_part_name)
                                    : mpModel->CreateModelPart(iga_model_part_name);

        // Create the True Model Part -> contains all the true boundary features
        std::string skin_model_part_inner_initial_name = "skin_model_part_inner_initial_name";
        std::string skin_model_part_outer_initial_name = "skin_model_part_outer_initial_name";
        std::string skin_model_part_name;
        if (mParameters.Has("skin_model_part_inner_initial_name")) {
            skin_model_part_inner_initial_name = mParameters["skin_model_part_inner_initial_name"].GetString();

            if (!mpModel->HasModelPart(skin_model_part_inner_initial_name)) 
            KRATOS_ERROR << "The skin_model_part '" << skin_model_part_inner_initial_name << "' was not created in the model.\n" 
                         << "Check the reading of the mdpa file in the import mdpa modeler."<< std::endl;
        }
        if (mParameters.Has("skin_model_part_outer_initial_name")) {
            skin_model_part_outer_initial_name = mParameters["skin_model_part_outer_initial_name"].GetString();

            if (!mpModel->HasModelPart(skin_model_part_outer_initial_name)) 
            KRATOS_ERROR << "The skin_model_part '" << skin_model_part_outer_initial_name << "' was not created in the model.\n" 
                         << "Check the reading of the mdpa file in the import mdpa modeler."<< std::endl;
        }
        // If there is not neither skin_inner nor skin_outer throw an error since you are using the sbm modeler
        if (!(mParameters.Has("skin_model_part_inner_initial_name") || mParameters.Has("skin_model_part_outer_initial_name"))){
        
            // Create the breps for the outer sbm boundary
            CreateBrepsSBMUtilities<Node, Point> CreateBrepsSBMUtilities(mEchoLevel);
            CreateBrepsSBMUtilities.CreateSurrogateBoundary(mpSurface, r_model_part, A_uvw, B_uvw);
            
            KRATOS_WARNING("None of the 'skin_model_part_name' have not been defined ") << 
                            "in the nurbs_geometry_modeler_sbm in the project paramer json" << std::endl;
            return;
        }
        
        if (mParameters.Has("skin_model_part_name"))
            skin_model_part_name = mParameters["skin_model_part_name"].GetString();
        else
            KRATOS_ERROR << "The skin_model_part name '" << skin_model_part_name << "' was not defined in the project parameters.\n" << std::endl;
 
        // inner
        ModelPart& skin_model_part_inner_initial = mpModel->HasModelPart(skin_model_part_inner_initial_name)
            ? mpModel->GetModelPart(skin_model_part_inner_initial_name)
            : mpModel->CreateModelPart(skin_model_part_inner_initial_name);
        // outer
        ModelPart& skin_model_part_outer_initial = mpModel->HasModelPart(skin_model_part_outer_initial_name)
            ? mpModel->GetModelPart(skin_model_part_outer_initial_name)
            : mpModel->CreateModelPart(skin_model_part_outer_initial_name);

        // Create the surrogate sub model parts inner and outer
        ModelPart& surrogate_sub_model_part_inner = iga_model_part.CreateSubModelPart("surrogate_inner");
        ModelPart& surrogate_sub_model_part_outer = iga_model_part.CreateSubModelPart("surrogate_outer");

        // Skin model part refined after Snake Process
        ModelPart& skin_model_part = mpModel->CreateModelPart(skin_model_part_name);
        skin_model_part.CreateSubModelPart("inner");
        skin_model_part.CreateSubModelPart("outer");
        
        
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
        SnakeSBMUtilities::CreateTheSnakeCoordinates(iga_model_part, skin_model_part_inner_initial, skin_model_part_outer_initial, skin_model_part, mEchoLevel,
                                                     unique_knot_vector_u, unique_knot_vector_v, mParameters) ;
        // Create the breps for the outer sbm boundary
        CreateBrepsSBMUtilities<Node, Point> CreateBrepsSBMUtilities(mEchoLevel);
        CreateBrepsSBMUtilities.CreateSurrogateBoundary(mpSurface, r_model_part, surrogate_sub_model_part_inner, surrogate_sub_model_part_outer, A_uvw, B_uvw);
    }

} // end namespace kratos
