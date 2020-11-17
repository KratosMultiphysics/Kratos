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

#if !defined(KRATOS_NURBS_VOLUME_GRID_MODELER_H_INCLUDED )
#define  KRATOS_NURBS_VOLUME_GRID_MODELER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/nurbs_volume_geometry.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_volume_utilities.h"
#include "includes/node.h"
#include "includes/model_part.h"

namespace Kratos {

class NurbsVolumeGridModeler
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( NurbsVolumeGridModeler );

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef Node<3> NodeType;

    typedef NurbsVolumeGeometry<PointerVector<NodeType>> NurbsVolumeGeometryType;
    typedef typename NurbsVolumeGeometryType::Pointer NurbsVolumeGeometryPointerType;

    ///@}
    ///@name Operations
    ///@{
    
    /**
     * @brief Creates a cartesian grid composed out of trivariant B-spline cubes.
     * @param PointA Lower point of bounding box.
     * @param PointB Upper point of bounding box.
     * @param Order  Polynomial degree in each direction u,v,w.
     * @param NumElements Number of equidistant elements/knot spans in each direction u,v,w.
     * @return Pointer to NurbsVolumeGeometry
     * @note The CP'S are defined as nodes and added to the rModelPart.
     **/ 
    static NurbsVolumeGeometryPointerType CreateGrid( ModelPart& rModelPart, const Point A, const Point B, SizeType OrderU, SizeType OrderV, SizeType OrderW,
        SizeType NumElementsU, SizeType NumElementsV, SizeType NumElementsW )
    {
        KRATOS_ERROR_IF( B.X() <= A.X() || B.Y() <= A.Y() || B.Z() <= A.Z() ) << "NurbsVolumeGridModeler: " 
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
                    points.push_back(NodeType::Pointer(new NodeType(node_id, x, y, z )));
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
        NurbsVolumeGeometryPointerType geometry
            = NurbsVolumeGeometryPointerType( new NurbsVolumeGeometryType(points, OrderU, OrderV, OrderW, knot_vector_u, knot_vector_v, knot_vector_w));

        // Set up knots for knot refinement according to the given number of elements in each direction.
        double delta_knot_u = 1.0 / NumElementsU;
        double knot_u = 0.0;
        std::vector<double> insert_knots_u(NumElementsU-1);
        for( IndexType i = 0; i < NumElementsU-1; i++){
            knot_u += delta_knot_u;
            insert_knots_u[i] = knot_u;
        }

        double delta_knot_v = 1.0 / NumElementsV;
        double knot_v = 0.0;
        std::vector<double> insert_knots_v(NumElementsV-1);
        for( IndexType i = 0; i < NumElementsV-1; i++){
            knot_v += delta_knot_v;
            insert_knots_v[i] = knot_v;
        }

        double delta_knot_w = 1.0 / NumElementsW;
        double knot_w = 0.0;
        std::vector<double> insert_knots_w(NumElementsW-1);
        for( IndexType i = 0; i < NumElementsW-1; i++){
            knot_w += delta_knot_w;
            insert_knots_w[i] = knot_w;
        }

        // Perform knot refinement.
        if( NumElementsU > 1)
            geometry = NurbsVolumeUtilities::KnotRefinementU(*geometry, insert_knots_u);
        if( NumElementsV > 1)
            geometry = NurbsVolumeUtilities::KnotRefinementV(*geometry, insert_knots_v);
        if( NumElementsW > 1)
            geometry = NurbsVolumeUtilities::KnotRefinementW(*geometry, insert_knots_w);

        // Add CP's to model part.
        for( int i = 0; i < geometry->size(); ++i){
            geometry->pGetPoint(i)->SetSolutionStepVariablesList(rModelPart.pGetNodalSolutionStepVariablesList());
            rModelPart.AddNode(geometry->pGetPoint(i),0);
        }
        return geometry;
    }

};

} // End namesapce Kratos
#endif // KRATOS_NURBS_VOLUME_GRID_MODELER_H_INCLUDED