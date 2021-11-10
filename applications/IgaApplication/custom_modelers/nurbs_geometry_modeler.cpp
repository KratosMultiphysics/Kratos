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
#include "nurbs_geometry_modeler.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities.h"

namespace Kratos
{

    ///@name Stages
    ///@{

    void NurbsGeometryModeler::SetupGeometryModel(){
        // Get bounding box
        KRATOS_ERROR_IF_NOT(mParameters.Has("lower_point"))
            << "NurbsGeometryModeler: Missing \"lower_point\" section" << std::endl;

        KRATOS_ERROR_IF_NOT( mParameters["lower_point"].GetVector().size() == 3 )
            << "NurbsGeometryModeler: \"lower_point\" must be defined in 3 dimensions." << std::endl;

        Point point_a(mParameters["lower_point"].GetVector());

        KRATOS_ERROR_IF_NOT(mParameters.Has("upper_point"))
            << "NurbsGeometryModeler: Missing \"upper_point\" section" << std::endl;

        KRATOS_ERROR_IF_NOT( mParameters["upper_point"].GetVector().size() == 3 )
            << "NurbsGeometryModeler: \"upper_point\" must be defined in 3 dimensions." << std::endl;

        Point point_b(mParameters["upper_point"].GetVector());

        // Get polynomial order and number of knotspans
        KRATOS_ERROR_IF_NOT(mParameters.Has("polynomial_order"))
            << "NurbsGeometryModeler: Missing \"polynomial_order\" section" << std::endl;

        KRATOS_ERROR_IF_NOT(mParameters.Has("number_of_knot_spans"))
            << "NurbsGeometryModeler: Missing \"number_of_knot_spans\" section" << std::endl;

        /// local space dimension is defined by the number of given polynomial orders.
        SizeType local_space_dimension =  mParameters["polynomial_order"].size();
        SizeType size_number_of_knot_spans =  mParameters["number_of_knot_spans"].size();

        KRATOS_ERROR_IF( local_space_dimension != size_number_of_knot_spans )
            << "Size of given Vectors: \"polynomial_order\" and \"number_of_knot_spans\" do not match." << std::endl;

        // Create model part in case it does not exist.
        KRATOS_ERROR_IF_NOT(mParameters.Has("model_part_name"))
            << "NurbsGeometryModeler: Missing \"model_part_name\" section" << std::endl;

        ModelPart& model_part = mpModel->HasModelPart(mParameters["model_part_name"].GetString()) ?
            mpModel->GetModelPart(mParameters["model_part_name"].GetString()) :
            mpModel->CreateModelPart(mParameters["model_part_name"].GetString());

        if( local_space_dimension == 2) {
            SizeType p_u =  mParameters["polynomial_order"].GetArrayItem(0).GetInt();
            SizeType p_v =  mParameters["polynomial_order"].GetArrayItem(1).GetInt();
            SizeType num_knot_span_u =  mParameters["number_of_knot_spans"].GetArrayItem(0).GetInt();
            SizeType num_knot_span_v =  mParameters["number_of_knot_spans"].GetArrayItem(1).GetInt();

            CreateAndAddRegularGrid2D(model_part, point_a, point_b, p_u, p_v, num_knot_span_u, num_knot_span_v);
        }
        else if( local_space_dimension == 3) {
            SizeType p_u =  mParameters["polynomial_order"].GetArrayItem(0).GetInt();
            SizeType p_v =  mParameters["polynomial_order"].GetArrayItem(1).GetInt();
            SizeType p_w =  mParameters["polynomial_order"].GetArrayItem(2).GetInt();
            SizeType num_knot_span_u =  mParameters["number_of_knot_spans"].GetArrayItem(0).GetInt();
            SizeType num_knot_span_w =  mParameters["number_of_knot_spans"].GetArrayItem(2).GetInt();
            SizeType num_knot_span_v =  mParameters["number_of_knot_spans"].GetArrayItem(1).GetInt();

            CreateAndAddRegularGrid3D(model_part, point_a, point_b, p_u, p_v, p_w, num_knot_span_u, num_knot_span_v, num_knot_span_w);
        }
        else {
            KRATOS_ERROR << "Nurbs Geometry Modeler is only available for surfaces and volumes." << std::endl;
        }

    }

    ///@}
    ///@name Private Operations
    ///@{
    void NurbsGeometryModeler::CreateAndAddRegularGrid2D( ModelPart& r_model_part, const Point& A, const Point& B,
        SizeType OrderU, SizeType OrderV,SizeType NumKnotSpansU, SizeType NumKnotSpansV)
    {
        KRATOS_ERROR_IF( B.X() <= A.X() || B.Y() <= A.Y() ) << "NurbsGeometryModeler: "
            << "The two Points A and B must meet the following requirement: (B-A) > (0,0,0). However, (B-A)=" << B-A << std::endl;

        PointerVector<NodeType> points;
        double delta_x = std::abs(B.X()-A.X())/OrderU;
        double delta_y = std::abs(B.Y()-A.Y())/OrderV;

        double x = A.X();
        double y = A.Y();

        for( IndexType j=0; j < OrderV+1; ++j){
            for( IndexType i = 0; i < OrderU+1; ++i){
                // Proper node ID is defined later.
                points.push_back(Kratos::make_intrusive<NodeType>(0, x, y, 0.0 ));
                x += delta_x;
            }
            x = A.X();
            y += delta_y;
        }

        // Set up knot vectors according to the given polynomial degree p.
        SizeType number_of_knots_u = 2*OrderU;
        Vector knot_vector_u(number_of_knots_u);
        for( IndexType i=0; i < number_of_knots_u; ++i){
            if( i < OrderU ){
                knot_vector_u[i] = 0.0;
            } else {
                knot_vector_u[i] = 1.0;
            }
        }
        SizeType number_of_knots_v = 2*OrderV;
        Vector knot_vector_v(number_of_knots_v);
        for( IndexType i=0; i < number_of_knots_v; ++i){
            if( i < OrderV ){
                knot_vector_v[i] = 0.0;
            } else {
                knot_vector_v[i] = 1.0;
            }
        }

        // Create bivariant nurbs surface.
        auto p_surface_geometry = Kratos::make_shared<NurbsSurfaceGeometryType>(
            points, OrderU, OrderV, knot_vector_u, knot_vector_v);

        // Set up knots for knot refinement according to the given number of elements in each direction.
        double delta_knot_u = 1.0 / NumKnotSpansU;
        double knot_u = 0.0;
        std::vector<double> insert_knots_u(NumKnotSpansU-1);
        for( IndexType i = 0; i < NumKnotSpansU-1; i++){
            knot_u += delta_knot_u;
            insert_knots_u[i] = knot_u;
        }

        double delta_knot_v = 1.0 / NumKnotSpansV;
        double knot_v = 0.0;
        std::vector<double> insert_knots_v(NumKnotSpansV-1);
        for( IndexType i = 0; i < NumKnotSpansV-1; i++){
            knot_v += delta_knot_v;
            insert_knots_v[i] = knot_v;
        }

        // Set up weights. Required in case no refinement is performed.
        Vector WeightsRefined(points.size(),1.0);

        // Add geometry to model part
        if( mParameters.Has("geometry_name") ){
            p_surface_geometry->SetId(mParameters["geometry_name"].GetString());
        } else {
            const SizeType number_of_geometries = r_model_part.NumberOfGeometries();
            SizeType last_geometry_id = 0;
            if( number_of_geometries > 0 ){
                for( auto it = r_model_part.GeometriesBegin(); it!= r_model_part.GeometriesEnd(); ++it){
                    last_geometry_id = it->Id();
                }
            }
            p_surface_geometry->SetId(last_geometry_id+1);
        }
        r_model_part.AddGeometry(p_surface_geometry);

        // Perform knot refinement.
        PointerVector<NodeType> PointsRefined = p_surface_geometry->Points();
        if( NumKnotSpansU > 1) {
            Vector KnotsURefined;
            PointsRefined = PointerVector<NodeType>(0);
            WeightsRefined.clear();

            NurbsSurfaceRefinementUtilities::KnotRefinementU( *p_surface_geometry, insert_knots_u,
                PointsRefined, KnotsURefined, WeightsRefined);

            p_surface_geometry->SetInternals(PointsRefined,
                p_surface_geometry->PolynomialDegreeU(), p_surface_geometry->PolynomialDegreeV(),
                KnotsURefined, p_surface_geometry->KnotsV(),
                WeightsRefined);
        }
        if( NumKnotSpansV > 1) {

            Vector KnotsVRefined;
            PointsRefined = PointerVector<NodeType>(0);
            WeightsRefined.clear();

            NurbsSurfaceRefinementUtilities::KnotRefinementV( *p_surface_geometry, insert_knots_v,
                PointsRefined, KnotsVRefined, WeightsRefined);

            p_surface_geometry->SetInternals(PointsRefined,
                p_surface_geometry->PolynomialDegreeU(), p_surface_geometry->PolynomialDegreeV(),
                p_surface_geometry->KnotsU(), KnotsVRefined,
                WeightsRefined);
        }

        IndexType node_id = 1;
        if( r_model_part.NumberOfNodes() > 0 ){
            node_id = (r_model_part.NodesEnd() - 1)->Id() + 1;
        }

        for (IndexType i = 0; i < PointsRefined.size(); ++i) {
            if (PointsRefined(i)->Id() == 0) {
                PointsRefined(i) = r_model_part.CreateNewNode(node_id, PointsRefined[i][0], PointsRefined[i][1], PointsRefined[i][2]);
                node_id++;
            }
        }

        // SetInternals again to make geometry point to the same nodes, as constructed in the model part.
        p_surface_geometry->SetInternals(PointsRefined,
            p_surface_geometry->PolynomialDegreeU(), p_surface_geometry->PolynomialDegreeV(),
            p_surface_geometry->KnotsU(), p_surface_geometry->KnotsV(), WeightsRefined);
    }


    void NurbsGeometryModeler::CreateAndAddRegularGrid3D( ModelPart& r_model_part, const Point& A, const Point& B, SizeType OrderU, SizeType OrderV, SizeType OrderW,
        SizeType NumKnotSpansU, SizeType NumKnotSpansV, SizeType NumKnotSpansW )
    {
        KRATOS_ERROR_IF( B.X() <= A.X() || B.Y() <= A.Y() || B.Z() <= A.Z() ) << "NurbsGeometryModeler: "
            << "The two Points A and B must meet the following requirement: (B-A) > (0,0,0). However, (B-A)=" << B-A << std::endl;

        // Set up control points.
        // Note: The CP's are uniformly subdivided with increasing p.
        //       This is allowed, as all lines of the cartesian grid are straight.
        PointerVector<NodeType> points;
        double delta_x = std::abs(B.X()-A.X())/OrderU;
        double delta_y = std::abs(B.Y()-A.Y())/OrderV;
        double delta_z = std::abs(B.Z()-A.Z())/OrderW;
        double x = A.X();
        double y = A.Y();
        double z = A.Z();
        for( IndexType k=0; k < OrderW+1; ++k){
            for( IndexType j=0; j < OrderV+1; ++j){
                for( IndexType i = 0; i < OrderU+1; ++i){
                    points.push_back(Kratos::make_intrusive<NodeType>(0, x, y, z ));
                    x += delta_x;
                }
                x = A.X();
                y += delta_y;
            }
            y = A.Y();
            z += delta_z;
        }
        // Set up knot vectors according to the given polynomial degree p.
        SizeType number_of_knots_u = 2*OrderU;
        Vector knot_vector_u(number_of_knots_u);
        for( IndexType i=0; i < number_of_knots_u; ++i){
            if( i < OrderU ){
                knot_vector_u[i] = 0.0;
            } else {
                knot_vector_u[i] = 1.0;
            }
        }
        SizeType number_of_knots_v = 2*OrderV;
        Vector knot_vector_v(number_of_knots_v);
        for( IndexType i=0; i < number_of_knots_v; ++i){
            if( i < OrderV ){
                knot_vector_v[i] = 0.0;
            } else {
                knot_vector_v[i] = 1.0;
            }
        }
        SizeType number_of_knots_w = 2*OrderW;
        Vector knot_vector_w(number_of_knots_w);
        for( IndexType i=0; i < number_of_knots_w; ++i){
            if( i < OrderW ){
                knot_vector_w[i] = 0.0;
            } else {
                knot_vector_w[i] = 1.0;
            }
        }

        // Create trivariant nurbs cube.
        auto p_volume_geometry = Kratos::make_shared<NurbsVolumeGeometryType>(
            points, OrderU, OrderV, OrderW, knot_vector_u, knot_vector_v, knot_vector_w);

        // Set up knots for knot refinement according to the given number of elements in each direction.
        double delta_knot_u = 1.0 / NumKnotSpansU;
        double knot_u = 0.0;
        std::vector<double> insert_knots_u(NumKnotSpansU-1);
        for( IndexType i = 0; i < NumKnotSpansU-1; i++){
            knot_u += delta_knot_u;
            insert_knots_u[i] = knot_u;
        }

        double delta_knot_v = 1.0 / NumKnotSpansV;
        double knot_v = 0.0;
        std::vector<double> insert_knots_v(NumKnotSpansV-1);
        for( IndexType i = 0; i < NumKnotSpansV-1; i++){
            knot_v += delta_knot_v;
            insert_knots_v[i] = knot_v;
        }

        double delta_knot_w = 1.0 / NumKnotSpansW;
        double knot_w = 0.0;
        std::vector<double> insert_knots_w(NumKnotSpansW-1);
        for( IndexType i = 0; i < NumKnotSpansW-1; i++){
            knot_w += delta_knot_w;
            insert_knots_w[i] = knot_w;
        }

        // Add geometry to model part
        if( mParameters.Has("geometry_name") ){
            p_volume_geometry->SetId(mParameters["geometry_name"].GetString());
        } else {
            const SizeType number_of_geometries = r_model_part.NumberOfGeometries();
            SizeType last_geometry_id = 0;
            if( number_of_geometries > 0 ){
                for( auto it = r_model_part.GeometriesBegin(); it!= r_model_part.GeometriesEnd(); ++it){
                    last_geometry_id = it->Id();
                }
            }
            p_volume_geometry->SetId(last_geometry_id+1);
        }
        r_model_part.AddGeometry(p_volume_geometry);

        // Perform knot refinement.
        PointerVector<NodeType> PointsRefined = p_volume_geometry->Points();

        if( NumKnotSpansU > 1) {
            Vector KnotsURefined;
            PointsRefined = PointerVector<NodeType>(0);

            NurbsVolumeRefinementUtilities::KnotRefinementU( *p_volume_geometry, insert_knots_u,
                PointsRefined, KnotsURefined);

            p_volume_geometry->SetInternals(PointsRefined,
                p_volume_geometry->PolynomialDegreeU(), p_volume_geometry->PolynomialDegreeV(), p_volume_geometry->PolynomialDegreeW(),
                KnotsURefined, p_volume_geometry->KnotsV(), p_volume_geometry->KnotsW());
        }
        if( NumKnotSpansV > 1) {
            Vector KnotsVRefined;
            PointsRefined = PointerVector<NodeType>(0);

            NurbsVolumeRefinementUtilities::KnotRefinementV( *p_volume_geometry, insert_knots_v,
                PointsRefined, KnotsVRefined);

            p_volume_geometry->SetInternals(PointsRefined,
                p_volume_geometry->PolynomialDegreeU(), p_volume_geometry->PolynomialDegreeV(), p_volume_geometry->PolynomialDegreeW(),
                p_volume_geometry->KnotsU(), KnotsVRefined, p_volume_geometry->KnotsW());
        }
        if( NumKnotSpansW > 1) {
            Vector KnotsWRefined;
            PointsRefined = PointerVector<NodeType>(0);

            NurbsVolumeRefinementUtilities::KnotRefinementW( *p_volume_geometry, insert_knots_w,
                PointsRefined, KnotsWRefined);

            p_volume_geometry->SetInternals(PointsRefined,
                p_volume_geometry->PolynomialDegreeU(), p_volume_geometry->PolynomialDegreeV(), p_volume_geometry->PolynomialDegreeW(),
                p_volume_geometry->KnotsU(), p_volume_geometry->KnotsV(), KnotsWRefined);
        }

        // Add nodes to model part
        IndexType node_id = 1;
        if( r_model_part.NumberOfNodes() > 0 ){
            node_id = (r_model_part.NodesEnd() - 1)->Id() + 1;
        }
        for (IndexType i = 0; i < PointsRefined.size(); ++i) {
            if (PointsRefined(i)->Id() == 0) {
                PointsRefined(i) = r_model_part.CreateNewNode(node_id, PointsRefined[i][0], PointsRefined[i][1], PointsRefined[i][2]);
                node_id++;
            }
        }
        // SetInternals again to make geometry point to the same nodes, as constructed in the model part.
        p_volume_geometry->SetInternals(PointsRefined,
            p_volume_geometry->PolynomialDegreeU(), p_volume_geometry->PolynomialDegreeV(), p_volume_geometry->PolynomialDegreeW(),
            p_volume_geometry->KnotsU(), p_volume_geometry->KnotsV(), p_volume_geometry->KnotsW());
    }
    ///@}
} // end namespace kratos
