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
#include "nurbs_geometry_modeler.h"
#include "includes/define.h"

namespace Kratos
{
    ///@name Stages
    ///@{

    void NurbsGeometryModeler::SetupGeometryModel(){
        // Get bounding box
        KRATOS_ERROR_IF_NOT(mParameters.Has("lower_point"))
            << "NurbsGeometryModeler: Missing \"lower_point\" section" << std::endl;
        Point point_a(mParameters["lower_point"].GetVector());

        KRATOS_ERROR_IF_NOT(mParameters.Has("upper_point"))
            << "NurbsGeometryModeler: Missing \"upper_point\" section" << std::endl;
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

        // if( local_space_dimension == 2) {
        //     SizeType p_u =  mParameters["polynomial_order"].GetArrayItem(0).GetInt();
        //     SizeType p_v =  mParameters["polynomial_order"].GetArrayItem(1).GetInt();
        //     SizeType num_knot_span_u =  mParameters["number_of_knot_spans"].GetArrayItem(0).GetInt();
        //     SizeType num_knot_span_v =  mParameters["number_of_knot_spans"].GetArrayItem(1).GetInt();
        //     //CreateGeometry2D(..)
        // }
        if( local_space_dimension == 3) {
            SizeType p_u =  mParameters["polynomial_order"].GetArrayItem(0).GetInt();
            SizeType p_v =  mParameters["polynomial_order"].GetArrayItem(1).GetInt();
            SizeType p_w =  mParameters["polynomial_order"].GetArrayItem(2).GetInt();
            SizeType num_knot_span_u =  mParameters["number_of_knot_spans"].GetArrayItem(0).GetInt();
            SizeType num_knot_span_w =  mParameters["number_of_knot_spans"].GetArrayItem(2).GetInt();
            SizeType num_knot_span_v =  mParameters["number_of_knot_spans"].GetArrayItem(1).GetInt();

            CreateGeometry3D(point_a, point_b, p_u, p_v, p_w, num_knot_span_u, num_knot_span_v, num_knot_span_w);
        }
        else {
            KRATOS_ERROR << "Murbs Geometry Modeler is not yet implemented for 1D and 2D geometries." << std::endl;
        }

        // Here add geometry and nodes to model part.
        KRATOS_ERROR_IF_NOT(mParameters.Has("model_part_name"))
            << "NurbsGeometryModeler: Missing \"model_part_name\" section" << std::endl;

        ModelPart& model_part = mpModel->HasModelPart(mParameters["model_part_name"].GetString()) ?
            mpModel->GetModelPart(mParameters["model_part_name"].GetString()) :
            mpModel->CreateModelPart(mParameters["model_part_name"].GetString());

        const SizeType number_of_geometries = model_part.NumberOfGeometries();
        mpGeometry->SetId(number_of_geometries+1);
        model_part.AddGeometry(mpGeometry);

        for( IndexType i = 0; i < mpGeometry->size(); ++i){
            mpGeometry->pGetPoint(i)->SetSolutionStepVariablesList(model_part.pGetNodalSolutionStepVariablesList());
            model_part.AddNode(mpGeometry->pGetPoint(i),0);
        }
    }

    ///@}
    ///@name Private Operations
    ///@{

    void NurbsGeometryModeler::CreateGeometry3D( const Point& A, const Point& B, SizeType OrderU, SizeType OrderV, SizeType OrderW,
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
        IndexType node_id = 1;
        for( IndexType k=0; k < OrderW+1; ++k){
            for( IndexType j=0; j < OrderV+1; ++j){
                for( IndexType i = 0; i < OrderU+1; ++i){
                    points.push_back(Kratos::make_intrusive<NodeType>(node_id, x, y, z ));
                    x += delta_x;
                    node_id++;
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
        mpGeometry = Kratos::make_shared<NurbsVolumeGeometryType>(
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

        // Perform knot refinement.
        if( NumKnotSpansU > 1) {
            mpGeometry = NurbsVolumeUtilities::KnotRefinementU(*mpGeometry, insert_knots_u);
        }
        if( NumKnotSpansV > 1) {
            mpGeometry = NurbsVolumeUtilities::KnotRefinementV(*mpGeometry, insert_knots_v);
        }
        if( NumKnotSpansW > 1) {
            mpGeometry = NurbsVolumeUtilities::KnotRefinementW(*mpGeometry, insert_knots_w);
        }
    }
    ///@}

} // end namespace kratos
