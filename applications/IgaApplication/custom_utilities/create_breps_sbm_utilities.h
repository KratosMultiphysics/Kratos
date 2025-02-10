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

# pragma once

// System includes

// External includes

// Project includes
#include "includes/io.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

// Geometries
#include "geometries/nurbs_curve_geometry.h"
#include "geometries/brep_surface.h"
#include "geometries/brep_curve.h"
#include "geometries/brep_curve_on_surface.h"


namespace Kratos
{

///@name Kratos Classes
///@{
template<class TNodeType = Node, class TEmbeddedNodeType = Point>
class CreateBrepsSbmUtilities : public IO
{
    public:

    ///@}
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CreateBrepsSbmUtilities
    KRATOS_CLASS_POINTER_DEFINITION(CreateBrepsSbmUtilities);

    using SizeType = std::size_t;
    using IndexType = std::size_t;

    using GeometryType = Geometry<TNodeType>;
    using GeometryPointerType = typename GeometryType::Pointer;
    using GeometrySurrogateArrayType = DenseVector<GeometryPointerType>;

    using ContainerNodeType = PointerVector<TNodeType>;
    using ContainerEmbeddedNodeType = PointerVector<TEmbeddedNodeType>;
    
    using NurbsSurfaceGeometryType = NurbsSurfaceGeometry<3, PointerVector<NodeType>>;
    using NurbsSurfaceGeometryPointerType = typename NurbsSurfaceGeometryType::Pointer;

    using BrepSurfaceType = BrepSurface<ContainerNodeType, true, ContainerEmbeddedNodeType>;
    using BrepCurveOnSurfaceType = BrepCurveOnSurface<ContainerNodeType, true, ContainerEmbeddedNodeType>;

    using BrepCurveOnSurfaceLoopType = DenseVector<typename BrepCurveOnSurfaceType::Pointer>;
    using BrepCurveOnSurfaceLoopArrayType = DenseVector<DenseVector<typename BrepCurveOnSurfaceType::Pointer>>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with path to input file.
    CreateBrepsSbmUtilities(
        SizeType EchoLevel = 0)
        : mEchoLevel(EchoLevel)
    {
    }

    /// Destructor.
    ~CreateBrepsSbmUtilities() = default;

    ///@}
    ///@name Python exposed Functions
    ///@{

    
    /**
     * @brief Create a Surrogate Boundary object
     * 
     * @param pSurface 
     * @param rModelPart 
     * @param rSurrogateModelPartInner 
     * @param rSurrogateModelPartOuter 
     * @param A_uvw 
     * @param B_uvw 
     * Adds the surface geometry to the herin provided model_part and create the boundary breps for the SBM case.
     */
    void CreateSurrogateBoundary(NurbsSurfaceGeometryPointerType& pSurface, 
                                 ModelPart& rModelPart, 
                                 ModelPart& rSurrogateModelPartInner, 
                                 ModelPart& rSurrogateModelPartOuter, 
                                 const Point& A_uvw, 
                                 const Point& B_uvw)
    {
        CreateBrepSurface(pSurface, rModelPart, rSurrogateModelPartInner, rSurrogateModelPartOuter, mEchoLevel);
        CreateBrepCurveOnSurfaces(pSurface, rModelPart, rSurrogateModelPartInner, rSurrogateModelPartOuter, A_uvw, B_uvw, mEchoLevel);
    }
    
    /**
     * @brief Create a Surrogate Boundary object
     * 
     * @param pSurface 
     * @param rModelPart 
     * @param A_uvw 
     * @param B_uvw 
     * Adds the surface geometry to the herin provided model_part and create the boundary breps when SBM is not needed.
     */
    void CreateSurrogateBoundary(NurbsSurfaceGeometryPointerType& pSurface, 
                                 ModelPart& rModelPart, 
                                 const Point& A_uvw, 
                                 const Point& B_uvw)
    {
        CreateBrepSurface(pSurface, rModelPart, mEchoLevel);
        IndexType id_brep_curve_on_surface = 2; // because id 1 is the brep surface
        CreateBrepCurvesOnRectangle(rModelPart, pSurface, A_uvw, B_uvw, id_brep_curve_on_surface);
    }


private:

    /**
     * @brief Create a Brep Surface object
     * 
     * @param pSurface 
     * @param rModelPart 
     * @param EchoLevel 
     */
    static void CreateBrepSurface(
        NurbsSurfaceGeometryPointerType pSurface,
        ModelPart& rModelPart,
        SizeType EchoLevel = 0)
    {
        KRATOS_INFO_IF("ReadBrepSurface", (EchoLevel > 3))
            << "Creating BrepSurface \""<< std::endl;

        BrepCurveOnSurfaceLoopArrayType outer_loops, inner_loops;

        auto p_brep_surface =
            Kratos::make_shared<BrepSurfaceType>(
                pSurface, 
                outer_loops,
                inner_loops,
                false);

        // Sets the brep as geometry parent of the nurbs surface.
        pSurface->SetGeometryParent(p_brep_surface.get());
        p_brep_surface->SetId(1);
        rModelPart.AddGeometry(p_brep_surface);
    }

    /**
     * @brief Create a Brep Surface object for the SBM case
     * 
     * @param pSurface 
     * @param rModelPart 
     * @param rSurrogateModelPartInner 
     * @param rSurrogateModelPartOuter 
     * @param EchoLevel 
     */
    static void CreateBrepSurface(
        NurbsSurfaceGeometryPointerType pSurface,
        ModelPart& rModelPart,
        ModelPart& rSurrogateModelPartInner, 
        ModelPart& rSurrogateModelPartOuter,
        SizeType EchoLevel = 0)
    {
        KRATOS_INFO_IF("ReadBrepSurface", (EchoLevel > 3))
            << "Creating BrepSurface \""<< std::endl;

        BrepCurveOnSurfaceLoopArrayType outer_loops, inner_loops;
        
        GeometrySurrogateArrayType surrogate_outer_loop_geometries(rSurrogateModelPartOuter.NumberOfConditions());
        GeometrySurrogateArrayType surrogate_inner_loop_geometries(rSurrogateModelPartInner.NumberOfConditions());
        int count = 0;
        for (auto i_cond : rSurrogateModelPartOuter.Conditions())
        {
            surrogate_outer_loop_geometries[count] = i_cond.pGetGeometry();
            count++;
        }

        count = 0;
        for (auto i_cond : rSurrogateModelPartInner.Conditions())
        {
            surrogate_inner_loop_geometries[count] = i_cond.pGetGeometry();
            count++;
        }

        auto p_brep_surface =
            Kratos::make_shared<BrepSurfaceType>(
                pSurface, 
                outer_loops,
                inner_loops);

        p_brep_surface->SetSurrogateInnerLoopGeometries(surrogate_inner_loop_geometries);
        p_brep_surface->SetSurrogateOuterLoopGeometries(surrogate_outer_loop_geometries);

        // Sets the brep as geometry parent of the nurbs surface.
        pSurface->SetGeometryParent(p_brep_surface.get());
        p_brep_surface->SetId(1);
        rModelPart.AddGeometry(p_brep_surface);
    }

    /**
     * @brief Create a Brep Curve On Surfaces object
     * 
     * @param pSurface 
     * @param rModelPart 
     * @param rSurrogateModelPartInner 
     * @param rSurrogateModelPartOuter 
     * @param A_uvw 
     * @param B_uvw 
     * @param EchoLevel 
     */
    static void CreateBrepCurveOnSurfaces(
        const NurbsSurfaceGeometryPointerType pSurface,
        ModelPart& rModelPart,
        const ModelPart& rSurrogateModelPartInner, 
        const ModelPart& rSurrogateModelPartOuter,
        const Point& A_uvw, 
        const Point& B_uvw,
        const SizeType EchoLevel = 0) {
    
        // OUTER 
        IndexType id_brep_curve_on_surface = 2; // because id 1 is the brep surface

        if (rSurrogateModelPartOuter.Nodes().size() > 0) {
            SizeType size_surrogate_loop_outer = rSurrogateModelPartOuter.Nodes().size();
            std::vector<double> surrogate_coords_x_outer(size_surrogate_loop_outer);
            std::vector<double> surrogate_coords_y_outer(size_surrogate_loop_outer);
            int count_surrogate_loop_outer = 0;
            for (auto i_node = rSurrogateModelPartOuter.NodesEnd()-1; i_node != rSurrogateModelPartOuter.NodesBegin()-1; i_node--) {
                surrogate_coords_x_outer[count_surrogate_loop_outer] = i_node->X();
                surrogate_coords_y_outer[count_surrogate_loop_outer] = i_node->Y();
                count_surrogate_loop_outer++;
            }
            std::vector<NurbsCurveGeometry<2, PointerVector<Point>>::Pointer> trimming_curves_outer;
            for (std::size_t i = 0; i < surrogate_coords_x_outer.size(); ++i) {
                Vector active_range_knot_vector = ZeroVector(2);
                
                Point::Pointer firstBrepPoint = Kratos::make_shared<Point>(surrogate_coords_x_outer[i], surrogate_coords_y_outer[i], 0.0);
                Point::Pointer second_brep_point = Kratos::make_shared<Point>(
                    surrogate_coords_x_outer[(i + 1) % surrogate_coords_x_outer.size()],  // Wrap around for the last point
                    surrogate_coords_y_outer[(i + 1) % surrogate_coords_y_outer.size()],  // Wrap around for the last point
                    0.0);
                // Compute the knot vector needed
                if (surrogate_coords_x_outer[(i + 1) % surrogate_coords_x_outer.size()]==surrogate_coords_x_outer[i]) {
                    active_range_knot_vector[0] = surrogate_coords_y_outer[i];
                    active_range_knot_vector[1] = surrogate_coords_y_outer[(i + 1) % surrogate_coords_y_outer.size()];
                }
                else {
                    active_range_knot_vector[0] = surrogate_coords_x_outer[i];
                    active_range_knot_vector[1] = surrogate_coords_x_outer[(i + 1) % surrogate_coords_x_outer.size()];
                }
                //// Order the active_range_knot_vector
                if (active_range_knot_vector[0] > active_range_knot_vector[1]) {
                    double temp = active_range_knot_vector[1];
                    active_range_knot_vector[1] = active_range_knot_vector[0] ;
                    active_range_knot_vector[0] = temp ;
                }
                NurbsCurveGeometry<2, PointerVector<Point>>::Pointer p_trimming_curve = CreateSingleBrep(firstBrepPoint, second_brep_point, active_range_knot_vector);
                trimming_curves_outer.push_back(p_trimming_curve);
            }

            BrepCurveOnSurfaceLoopType trimming_brep_curve_vector_outer(surrogate_coords_x_outer.size());

            for (std::size_t i = 0; i < trimming_curves_outer.size(); ++i) {
                Vector active_range_vector = ZeroVector(2);
                if (surrogate_coords_x_outer[(i + 1) % surrogate_coords_x_outer.size()] == surrogate_coords_x_outer[i]) {
                    active_range_vector[0] = surrogate_coords_y_outer[i];
                    active_range_vector[1] = surrogate_coords_y_outer[(i + 1) % surrogate_coords_x_outer.size()];
                }
                else {
                    active_range_vector[0] = surrogate_coords_x_outer[i];
                    active_range_vector[1] = surrogate_coords_x_outer[(i + 1) % surrogate_coords_x_outer.size()];
                }
                bool curve_direction = true;

                // Re-order
                if (active_range_vector[0] > active_range_vector[1]) {
                    double temp = active_range_vector[1];
                    active_range_vector[1] = active_range_vector[0] ;
                    active_range_vector[0] = temp ;
                }

                NurbsInterval brep_active_range(active_range_vector[0], active_range_vector[1]);

                auto p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(
                    pSurface, trimming_curves_outer[i], brep_active_range, curve_direction);
                p_brep_curve_on_surface->SetId(id_brep_curve_on_surface);

                rModelPart.AddGeometry(p_brep_curve_on_surface);
                id_brep_curve_on_surface++;

            }

        } else {
            CreateBrepCurvesOnRectangle(rModelPart, pSurface, A_uvw, B_uvw, id_brep_curve_on_surface);
        }

        // INNER
        for (IndexType iel = 1; iel < rSurrogateModelPartInner.Elements().size()+1; iel++) {
            int first_surrogate_node_id = rSurrogateModelPartInner.pGetElement(iel)->GetGeometry()[0].Id(); // Element 1 because is the only surrogate loop
            int last_surrogate_node_id = rSurrogateModelPartInner.pGetElement(iel)->GetGeometry()[1].Id();  // Element 1 because is the only surrogate loop
            int size_surrogate_loop_inner = last_surrogate_node_id - first_surrogate_node_id + 1;

            int count_surrogate_loop_inner = 0;
            std::vector<double> surrogate_coords_x_inner(size_surrogate_loop_inner);
            std::vector<double> surrogate_coords_y_inner(size_surrogate_loop_inner);
            for (int id_node = first_surrogate_node_id; id_node < last_surrogate_node_id+1; id_node++) {
                surrogate_coords_x_inner[count_surrogate_loop_inner] = rSurrogateModelPartInner.GetNode(id_node).X();
                surrogate_coords_y_inner[count_surrogate_loop_inner] = rSurrogateModelPartInner.GetNode(id_node).Y();
                count_surrogate_loop_inner++;
            }
            std::vector<NurbsCurveGeometry<2, PointerVector<Point>>::Pointer> trimming_curves_inner;

            for (std::size_t i = 0; i < surrogate_coords_x_inner.size(); ++i) {
                Vector active_range_knot_vector = ZeroVector(2);
                
                Point::Pointer first_brep_point = Kratos::make_shared<Point>(surrogate_coords_x_inner[i], surrogate_coords_y_inner[i], 0.0);
                Point::Pointer second_brep_point = Kratos::make_shared<Point>(
                    surrogate_coords_x_inner[(i + 1) % surrogate_coords_x_inner.size()],  // Wrap around for the last point
                    surrogate_coords_y_inner[(i + 1) % surrogate_coords_y_inner.size()],  // Wrap around for the last point
                    0.0);
                // Compute the knot vector needed
                if (surrogate_coords_x_inner[(i + 1) % surrogate_coords_x_inner.size()]==surrogate_coords_x_inner[i]) {
                    active_range_knot_vector[0] = surrogate_coords_y_inner[i];
                    active_range_knot_vector[1] = surrogate_coords_y_inner[(i + 1) % surrogate_coords_y_inner.size()];
                }
                else {
                    active_range_knot_vector[0] = surrogate_coords_x_inner[i];
                    active_range_knot_vector[1] = surrogate_coords_x_inner[(i + 1) % surrogate_coords_x_inner.size()];
                }
                // Order the active_range_knot_vector
                if (active_range_knot_vector[0] > active_range_knot_vector[1]) {
                    double temp = active_range_knot_vector[1];
                    active_range_knot_vector[1] = active_range_knot_vector[0] ;
                    active_range_knot_vector[0] = temp ;
                }
                NurbsCurveGeometry<2, PointerVector<Point>>::Pointer p_trimming_curve = CreateSingleBrep(first_brep_point, second_brep_point, active_range_knot_vector);
                trimming_curves_inner.push_back(p_trimming_curve);
            }

            BrepCurveOnSurfaceLoopType trimming_brep_curve_vector(surrogate_coords_x_inner.size());


            for (std::size_t i = 0; i < trimming_curves_inner.size(); ++i) {
                Vector active_range_vector = ZeroVector(2);
                if (surrogate_coords_x_inner[(i + 1) % surrogate_coords_x_inner.size()]==surrogate_coords_x_inner[i]) {
                    active_range_vector[0] = surrogate_coords_y_inner[i];
                    active_range_vector[1] = surrogate_coords_y_inner[(i + 1) % surrogate_coords_x_inner.size()];
                }
                else {
                    active_range_vector[0] = surrogate_coords_x_inner[i];
                    active_range_vector[1] = surrogate_coords_x_inner[(i + 1) % surrogate_coords_x_inner.size()];
                }
                bool curve_direction = true;

                // re-order
                if (active_range_vector[0] > active_range_vector[1]) {
                    double temp = active_range_vector[1];
                    active_range_vector[1] = active_range_vector[0] ;
                    active_range_vector[0] = temp ;
                }

                NurbsInterval brep_active_range(active_range_vector[0], active_range_vector[1]);

                auto p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(
                    pSurface, trimming_curves_inner[i], brep_active_range, curve_direction);
                p_brep_curve_on_surface->SetId(id_brep_curve_on_surface);
                rModelPart.AddGeometry(p_brep_curve_on_surface);
                id_brep_curve_on_surface++;
                trimming_brep_curve_vector[i] = p_brep_curve_on_surface ;

            }

        }
    } 
    
    /**
     * @brief Create a Single Brep object
     * 
     * @param firstBrepPoint 
     * @param secondBrepPoint 
     * @param activeRangeKnotVector 
     * @return NurbsCurveGeometry<2, PointerVector<Point>>::Pointer 
     */
    static typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer CreateSingleBrep(
        const Point::Pointer firstBrepPoint,
        const Point::Pointer secondBrepPoint, 
        const Vector activeRangeKnotVector)
        {
            // Create the data for the trimming curves
            PointerVector<Point> control_points;
            control_points.push_back(firstBrepPoint);
            control_points.push_back(secondBrepPoint);
            const int polynomial_degree = 1;
            Vector knot_vector = ZeroVector(4) ;
            knot_vector[0] = activeRangeKnotVector[0] ;
            knot_vector[1] = activeRangeKnotVector[0] ;
            knot_vector[2] = activeRangeKnotVector[1] ;
            knot_vector[3] = activeRangeKnotVector[1] ;
            // Create the trimming curves
            typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer p_trimming_curve(
                new NurbsCurveGeometry<2, PointerVector<Point>>(
                    control_points,
                    polynomial_degree,
                    knot_vector));   
            return p_trimming_curve;
        }
    
    /**
     * @brief Create a Brep Curves On Rectangle object
     * 
     * @param rModelPart 
     * @param pSurfaceGeometry 
     * @param A_uvw 
     * @param B_uvw 
     * @param lastGeometryId 
     */
    static void CreateBrepCurvesOnRectangle(ModelPart& rModelPart, 
                                            const NurbsSurfaceGeometryPointerType pSurfaceGeometry, 
                                            const Point& A_uvw, 
                                            const Point& B_uvw, 
                                            IndexType &lastGeometryId) {
        Vector knot_vector = ZeroVector(2);
        knot_vector[0] = 0.0;
        knot_vector[1] = std::abs(B_uvw[0] - A_uvw[0]);
        const int p = 1;

        NurbsCurveGeometry<2, PointerVector<Point>>::PointsArrayType segment1;
        segment1.push_back(Point::Pointer(new Point(A_uvw[0], A_uvw[1])));
        segment1.push_back(Point::Pointer(new Point(B_uvw[0], A_uvw[1])));

        NurbsCurveGeometry<2, PointerVector<Point>>::PointsArrayType segment2;
        segment2.push_back(Point::Pointer(new Point(B_uvw[0], A_uvw[1])));
        segment2.push_back(Point::Pointer(new Point(B_uvw[0], B_uvw[1])));
        
        NurbsCurveGeometry<2, PointerVector<Point>>::PointsArrayType segment3;
        segment3.push_back(Point::Pointer(new Point(B_uvw[0], B_uvw[1])));
        segment3.push_back(Point::Pointer(new Point(A_uvw[0], B_uvw[1])));
        
        NurbsCurveGeometry<2, PointerVector<Point>>::PointsArrayType segment4;
        segment4.push_back(Point::Pointer(new Point(A_uvw[0], B_uvw[1])));
        segment4.push_back(Point::Pointer(new Point(A_uvw[0], A_uvw[1])));

        auto p_curve_1 = Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(segment1, p, knot_vector);
        auto p_curve_2 = Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(segment2, p, knot_vector);
        auto p_curve_3 = Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(segment3, p, knot_vector);
        auto p_curve_4 = Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(segment4, p, knot_vector);
        
        auto brep_curve_on_surface = BrepCurveOnSurface< PointerVector<NodeType>, true, PointerVector<Point>>(pSurfaceGeometry, p_curve_1);
        auto p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(pSurfaceGeometry, p_curve_1);
        p_brep_curve_on_surface->SetId(lastGeometryId);
        rModelPart.AddGeometry(p_brep_curve_on_surface);
        
        brep_curve_on_surface = BrepCurveOnSurface< PointerVector<NodeType>, true, PointerVector<Point>>(pSurfaceGeometry, p_curve_2);
        p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(pSurfaceGeometry, p_curve_2);
        p_brep_curve_on_surface->SetId(++lastGeometryId);
        rModelPart.AddGeometry(p_brep_curve_on_surface);

        brep_curve_on_surface = BrepCurveOnSurface< PointerVector<NodeType>, true, PointerVector<Point>>(pSurfaceGeometry, p_curve_3);
        p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(pSurfaceGeometry, p_curve_3);
        p_brep_curve_on_surface->SetId(++lastGeometryId);
        rModelPart.AddGeometry(p_brep_curve_on_surface);
        
        brep_curve_on_surface = BrepCurveOnSurface< PointerVector<NodeType>, true, PointerVector<Point>>(pSurfaceGeometry, p_curve_4);
        p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(pSurfaceGeometry, p_curve_4);
        p_brep_curve_on_surface->SetId(++lastGeometryId);
        rModelPart.AddGeometry(p_brep_curve_on_surface);

        lastGeometryId++;
    }
    
    ///@}
    ///@name Utility functions
    ///@{

    int mEchoLevel;

    ///@}
}; // Class CreateBrepsSbmUtilities
}  // namespace Kratos.
