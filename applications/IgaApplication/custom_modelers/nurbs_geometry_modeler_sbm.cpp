//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
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
        KRATOS_INFO_IF("NurbsGeometryModelerSbm", mEchoLevel > 1) << "[NURBS MODELER]:: STARTING" << std::endl;
        // Get bounding box physical space.
        KRATOS_ERROR_IF_NOT(mParameters.Has("lower_point_xyz"))
            << "NurbsGeometryModelerSbm: Missing \"lower_point_xyz\" section" << std::endl;

        KRATOS_ERROR_IF_NOT( mParameters["lower_point_xyz"].GetVector().size() == 3 )
            << "NurbsGeometryModelerSbm: \"lower_point_xyz\" must be defined in 3 dimensions." << std::endl;

        Point point_a_xyz(mParameters["lower_point_xyz"].GetVector());

        KRATOS_ERROR_IF_NOT(mParameters.Has("upper_point_xyz"))
            << "NurbsGeometryModelerSbm: Missing \"upper_point_xyz\" section" << std::endl;

        KRATOS_ERROR_IF_NOT( mParameters["upper_point_xyz"].GetVector().size() == 3 )
            << "NurbsGeometryModelerSbm: \"upper_point_xyz\" must be defined in 3 dimensions." << std::endl;

        Point point_b_xyz(mParameters["upper_point_xyz"].GetVector());

        // Get bounding box parametric space.
        KRATOS_ERROR_IF_NOT(mParameters.Has("lower_point_uvw"))
            << "NurbsGeometryModelerSbm: Missing \"lower_point_uvw\" section" << std::endl;

        KRATOS_ERROR_IF_NOT( mParameters["lower_point_uvw"].GetVector().size() == 3 )
            << "NurbsGeometryModelerSbm: \"lower_point_uvw\" must be defined in 3 dimensions." << std::endl;

        Point point_a_uvw(mParameters["lower_point_uvw"].GetVector());

        KRATOS_ERROR_IF_NOT(mParameters.Has("upper_point_uvw"))
            << "NurbsGeometryModelerSbm: Missing \"upper_point_uvw\" section" << std::endl;

        KRATOS_ERROR_IF_NOT( mParameters["upper_point_uvw"].GetVector().size() == 3 )
            << "NurbsGeometryModelerSbm: \"upper_point_uvw\" must be defined in 3 dimensions." << std::endl;

        Point point_b_uvw(mParameters["upper_point_uvw"].GetVector());

        // Get polynomial order and number of knotspans
        KRATOS_ERROR_IF_NOT(mParameters.Has("polynomial_order"))
            << "NurbsGeometryModelerSbm: Missing \"polynomial_order\" section" << std::endl;

        KRATOS_ERROR_IF_NOT(mParameters.Has("number_of_knot_spans"))
            << "NurbsGeometryModelerSbm: Missing \"number_of_knot_spans\" section" << std::endl;

        /// local space dimension is defined by the number of given polynomial orders.
        SizeType local_space_dimension =  mParameters["polynomial_order"].size();
        SizeType size_number_of_knot_spans =  mParameters["number_of_knot_spans"].size();

        KRATOS_ERROR_IF( local_space_dimension != size_number_of_knot_spans )
            << "Size of given Vectors: \"polynomial_order\" and \"number_of_knot_spans\" do not match." << std::endl;

        // Create model part in case it does not exist.
        KRATOS_ERROR_IF_NOT(mParameters.Has("model_part_name"))
            << "NurbsGeometryModelerSbm: Missing \"model_part_name\" section" << std::endl;

        ModelPart& model_part = mpModel->HasModelPart(mParameters["model_part_name"].GetString()) ?
            mpModel->GetModelPart(mParameters["model_part_name"].GetString()) :
            mpModel->CreateModelPart(mParameters["model_part_name"].GetString());

        if( local_space_dimension == 2) {
            SizeType p_u =  mParameters["polynomial_order"].GetArrayItem(0).GetInt();
            SizeType p_v =  mParameters["polynomial_order"].GetArrayItem(1).GetInt();
            SizeType num_knot_span_u =  mParameters["number_of_knot_spans"].GetArrayItem(0).GetInt();
            SizeType num_knot_span_v =  mParameters["number_of_knot_spans"].GetArrayItem(1).GetInt();

            CreateAndAddRegularGrid2D(model_part, point_a_xyz, point_b_xyz, point_a_uvw, point_b_uvw, p_u, p_v, num_knot_span_u, num_knot_span_v);
            KRATOS_INFO_IF("NurbsGeometryModelerSbm", mEchoLevel > 1) << "[NURBS MODELER]:: grid 2D created" << std::endl;
        }
        else if( local_space_dimension == 3) {
            SizeType p_u =  mParameters["polynomial_order"].GetArrayItem(0).GetInt();
            SizeType p_v =  mParameters["polynomial_order"].GetArrayItem(1).GetInt();
            SizeType p_w =  mParameters["polynomial_order"].GetArrayItem(2).GetInt();
            SizeType num_knot_span_u =  mParameters["number_of_knot_spans"].GetArrayItem(0).GetInt();
            SizeType num_knot_span_w =  mParameters["number_of_knot_spans"].GetArrayItem(2).GetInt();
            SizeType num_knot_span_v =  mParameters["number_of_knot_spans"].GetArrayItem(1).GetInt();

            // TODO -> 3D
        }
        else {
            KRATOS_ERROR << "Nurbs Geometry Modeler is only available for surfaces and volumes." << std::endl;
        }

    }

    ///@}
    ///@name Private Operations
    ///@{
    void NurbsGeometryModelerSbm::CreateAndAddRegularGrid2D( ModelPart& r_model_part, const Point& A_xyz, const Point& B_xyz,
        const Point& A_uvw, const Point& B_uvw, SizeType OrderU, SizeType OrderV,SizeType NumKnotSpansU, SizeType NumKnotSpansV)
    {
        KRATOS_ERROR_IF( B_xyz.X() <= A_xyz.X() || B_xyz.Y() <= A_xyz.Y() ) << "NurbsGeometryModelerSbm: "
            << "The two Points A_xyz and B_xyz must meet the following requirement: (B_xyz-A_xyz) > (0,0,0). However, (B_xyz-A_xyz)=" << B_xyz-A_xyz << std::endl;

        KRATOS_ERROR_IF( B_uvw.X() <= A_uvw.X() || B_uvw.Y() <= A_uvw.Y() ) << "NurbsGeometryModelerSbm: "
            << "The two Points A_uvw and B_uvw must meet the following requirement: (B_uvw-A_uvw) > (0,0,0). However, (B_uvw-A_uvw)=" << B_uvw-A_uvw << std::endl;


        PointerVector<NodeType> points;
        double delta_x = std::abs(B_xyz.X()-A_xyz.X())/OrderU;
        double delta_y = std::abs(B_xyz.Y()-A_xyz.Y())/OrderV;

        double x = A_xyz.X();
        double y = A_xyz.Y();

        for( IndexType j=0; j < OrderV+1; ++j){
            for( IndexType i = 0; i < OrderU+1; ++i){
                // Proper node ID is defined later.
                points.push_back(Kratos::make_intrusive<NodeType>(0, x, y, 0.0 ));
                x += delta_x;
            }
            x = A_xyz.X();
            y += delta_y;
        }

        // Set up knot vectors according to the given polynomial degree p.
        SizeType number_of_knots_u = 2*OrderU;
        Vector knot_vector_u(number_of_knots_u);
        for( IndexType i=0; i < number_of_knots_u; ++i){
            if( i < OrderU ){
                knot_vector_u[i] = A_uvw[0];
            } else {
                knot_vector_u[i] = B_uvw[0];
            }
        }
        SizeType number_of_knots_v = 2*OrderV;
        Vector knot_vector_v(number_of_knots_v);
        for( IndexType i=0; i < number_of_knots_v; ++i){
            if( i < OrderV ){
                knot_vector_v[i] = A_uvw[1];
            } else {
                knot_vector_v[i] = B_uvw[1];
            }
        }

        // Create bivariant nurbs surface.
        auto p_surface_geometry = Kratos::make_shared<NurbsSurfaceGeometryType>(
            points, OrderU, OrderV, knot_vector_u, knot_vector_v);


        // Set up knots for knot refinement according to the given number of elements in each direction.
        double delta_knot_u = (B_uvw[0]-A_uvw[0]) / NumKnotSpansU;
        double knot_u = A_uvw[0];
        std::vector<double> insert_knots_u(NumKnotSpansU-1);
        for( IndexType i = 0; i < NumKnotSpansU-1; i++){
            knot_u += delta_knot_u;
            insert_knots_u[i] = knot_u;
        }

        double delta_knot_v = (B_uvw[1]-A_uvw[1]) / NumKnotSpansV;
        double knot_v = A_uvw[1];
        std::vector<double> insert_knots_v(NumKnotSpansV-1);
        for( IndexType i = 0; i < NumKnotSpansV-1; i++){
            knot_v += delta_knot_v;
            insert_knots_v[i] = knot_v;
        }

        // Set up weights. Required in case no refinement is performed.
        Vector WeightsRefined(points.size(),1.0);


        KRATOS_INFO_IF("::[NurbsGeometryModelerSbm]::", mEchoLevel > 0) << "Ending the CreateTheSnakeCoordinates" << std::endl;

        KRATOS_ERROR_IF_NOT(mParameters.Has("model_part_name"))
            << "Missing \"domain_model_part_name\" in NurbsGeometryModelerSbm Parameters." << std::endl;
            
        
        // Create the Domain/Iga Model Part
        const std::string iga_model_part_name = mParameters["model_part_name"].GetString();
        ModelPart& iga_model_part = mpModel->HasModelPart(iga_model_part_name)
                                    ? mpModel->GetModelPart(iga_model_part_name)
                                    : mpModel->CreateModelPart(iga_model_part_name);

        // Create the True Model Part -> contains all the true boundary features
        std::string skin_model_part_name;
        if (mParameters.Has("skin_model_part_name")) {
            skin_model_part_name = mParameters["skin_model_part_name"].GetString();
        }
        else {
            KRATOS_ERROR << "The 'skin_model_part_name' has not been defined" << 
                            "in the nurbs_geometry_modeler_sbm, please define it" << 
                            " in the project paramer json" << std::endl;
        }
        // const std::string skin_model_part_name = mParameters["skin_model_part_name"].GetString();
        if (!mpModel->HasModelPart(skin_model_part_name)) 
            KRATOS_ERROR << "The skin_model_part '" << skin_model_part_name << "' was not created in the model.\n" 
                         << "Check the reading of the mdpa file in the import mdpa modeler."<< std::endl;

        ModelPart& skin_model_part = mpModel->GetModelPart(skin_model_part_name);

        // TODO: to be deleted
        //// Read refinement parameters
        double knot_step_u; double knot_step_v;
        const Parameters refinements_parameters = ReadParamatersFile("refinements.iga.json");
        Parameters refinements_parameters_model_part;
        bool isDefined = false;
        for (int i = 0; i < refinements_parameters["refinements"].size(); i++) {
            if (refinements_parameters["refinements"][i]["model_part_name"].GetString() == iga_model_part_name) {
                refinements_parameters_model_part = refinements_parameters["refinements"][i];
                isDefined = true;
                break;
            }
        }

        ModelPart& surrogate_sub_model_part_inner = iga_model_part.CreateSubModelPart("surrogate_inner");
        ModelPart& surrogate_sub_model_part_outer = iga_model_part.CreateSubModelPart("surrogate_outer");

        // Skin model part refined after Snake Process
        ModelPart& skin_sub_model_part_in = skin_model_part.CreateSubModelPart("inner");
        ModelPart& skin_sub_model_part_out = skin_model_part.CreateSubModelPart("outer");
        
        // Suppose to be in sbm case
        SnakeSBMUtilitiesNew::CreateTheSnakeCoordinates(iga_model_part, skin_model_part, mEchoLevel,
                                                     knot_vector_u, knot_vector_v, mParameters) ;

        KRATOS_WATCH("TO THE END")
        exit(0);
        
        // Create the breps for the outer sbm boundary
        // CreateBrepsSBMUtilities<Node, Point> CreateBrepsSBMUtilities(mEchoLevel);

        // CreateBrepsSBMUtilities.CreateSurrogateBoundary(p_surface_geometry, r_model_part, surrogate_model_part_inner, surrogate_model_part_outer, A_uvw, B_uvw);
        // // Perform knot refinement.
        // PointerVector<NodeType> PointsRefined = p_surface_geometry->Points();
        // if( NumKnotSpansU > 1) {
        //     Vector KnotsURefined;
        //     PointsRefined = PointerVector<NodeType>(0);
        //     WeightsRefined.clear();

        //     NurbsSurfaceRefinementUtilities::KnotRefinementU( *p_surface_geometry, insert_knots_u,
        //         PointsRefined, KnotsURefined, WeightsRefined);

        //     p_surface_geometry->SetInternals(PointsRefined,
        //         p_surface_geometry->PolynomialDegreeU(), p_surface_geometry->PolynomialDegreeV(),
        //         KnotsURefined, p_surface_geometry->KnotsV(),
        //         WeightsRefined);
        // }

        // if( NumKnotSpansV > 1) {

        //     Vector KnotsVRefined;
        //     PointsRefined = PointerVector<NodeType>(0);
        //     WeightsRefined.clear();

        //     NurbsSurfaceRefinementUtilities::KnotRefinementV( *p_surface_geometry, insert_knots_v,
        //         PointsRefined, KnotsVRefined, WeightsRefined);

        //     p_surface_geometry->SetInternals(PointsRefined,
        //         p_surface_geometry->PolynomialDegreeU(), p_surface_geometry->PolynomialDegreeV(),
        //         p_surface_geometry->KnotsU(), KnotsVRefined,
        //         WeightsRefined);
        // }

        // IndexType node_id = 1;
        // // if( r_model_part.NumberOfNodes() > 0 ){
        // //     node_id = (r_model_part.NodesEnd() - 1)->Id() + 1;
        // // }

        // if(r_model_part.GetParentModelPart().NumberOfNodes() > 0 ){
        //     node_id = (r_model_part.GetParentModelPart().NodesEnd() - 1)->Id() + 1;
        // }

        // for (IndexType i = 0; i < PointsRefined.size(); ++i) {
        //     if (PointsRefined(i)->Id() == 0) {
        //         PointsRefined(i) = r_model_part.CreateNewNode(node_id, PointsRefined[i][0], PointsRefined[i][1], PointsRefined[i][2]);
        //         node_id++;
        //     }
        // }

        // // SetInternals again to make geometry point to the same nodes, as constructed in the model part.
        // p_surface_geometry->SetInternals(PointsRefined,
        //     p_surface_geometry->PolynomialDegreeU(), p_surface_geometry->PolynomialDegreeV(),
        //     p_surface_geometry->KnotsU(), p_surface_geometry->KnotsV(), WeightsRefined);
    }


    Parameters NurbsGeometryModelerSbm::ReadParamatersFile(
        const std::string& rDataFileName) const
    {
        // Check if rDataFileName ends with ".cad.json" and add it if needed.
        const std::string data_file_name = (rDataFileName.compare(rDataFileName.size() - 9, 9, ".iga.json") != 0)
            ? rDataFileName + ".iga.json"
            : rDataFileName;

        std::ifstream infile(data_file_name);
        KRATOS_ERROR_IF_NOT(infile.good()) << "Physics fil: "
            << data_file_name << " cannot be found." << std::endl;
        KRATOS_INFO_IF("ReadParamatersFile", mEchoLevel > 3)
            << "Reading file: \"" << data_file_name << "\"" << std::endl;

        std::stringstream buffer;
        buffer << infile.rdbuf();

        return Parameters(buffer.str());
    }
} // end namespace kratos
