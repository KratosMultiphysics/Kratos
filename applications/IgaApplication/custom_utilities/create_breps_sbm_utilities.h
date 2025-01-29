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
#include "geometries/coupling_geometry.h"
#include "geometries/point_on_geometry.h"

#include "geometries/nurbs_curve_geometry.h"
#include "geometries/nurbs_surface_geometry.h"

#include "geometries/brep_surface.h"
#include "geometries/brep_curve.h"

#include "geometries/nurbs_curve_on_surface_geometry.h"
#include "geometries/brep_curve_on_surface.h"


namespace Kratos
{

///@name Kratos Classes
///@{
template<class TNodeType = Node, class TEmbeddedNodeType = Point>
class CreateBrepsSBMUtilities : public IO
{
    public:

    ///@}
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CreateBrepsSBMUtilities
    KRATOS_CLASS_POINTER_DEFINITION(CreateBrepsSBMUtilities);

    using SizeType = std::size_t;
    using IndexType = std::size_t;

    using GeometryType = Geometry<TNodeType>;
    using GeometryPointerType = typename GeometryType::Pointer;

    using ContainerNodeType = PointerVector<TNodeType>;
    using ContainerEmbeddedNodeType = PointerVector<TEmbeddedNodeType>;

    using CouplingGeometryType = CouplingGeometry<TNodeType>;

    using NurbsSurfaceType = NurbsSurfaceGeometry<3, ContainerNodeType>;
    using NurbsTrimmingCurveType = NurbsCurveGeometry<2, ContainerEmbeddedNodeType>;

    using NurbsSurfacePointerType = typename NurbsSurfaceType::Pointer;
    using NurbsTrimmingCurvePointerType = typename NurbsTrimmingCurveType::Pointer;

    using NurbsSurfaceGeometryType = NurbsSurfaceGeometry<3, PointerVector<NodeType>>;
    using NurbsSurfaceGeometryPointerType = typename NurbsSurfaceGeometryType::Pointer;

    using BrepSurfaceType = BrepSurface<ContainerNodeType, true, ContainerEmbeddedNodeType>;
    using BrepCurveOnSurfaceType = BrepCurveOnSurface<ContainerNodeType, true, ContainerEmbeddedNodeType>;

    using BrepCurveType = BrepCurve<ContainerNodeType, ContainerEmbeddedNodeType>;
    using PointOnGeometryOnSurfaceType = PointOnGeometry<ContainerNodeType, 3, 2>;
    using PointOnGeometryOnCurveType = PointOnGeometry<ContainerNodeType, 3, 1>;

    using BrepCurveOnSurfaceArrayType = DenseVector<typename BrepCurveOnSurfaceType::Pointer>;
    using BrepCurveOnSurfaceLoopType = DenseVector<typename BrepCurveOnSurfaceType::Pointer>;
    using BrepCurveOnSurfaceLoopArrayType = DenseVector<DenseVector<typename BrepCurveOnSurfaceType::Pointer>>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with path to input file.
    CreateBrepsSBMUtilities(
        SizeType EchoLevel = 0)
        : mEchoLevel(EchoLevel)
    {
    }

    /// Destructor.
    ~CreateBrepsSBMUtilities() = default;

    ///@}
    ///@name Python exposed Functions
    ///@{

    /// Adds the surface geometry to the herin provided model_part and create the boundary breps for the SBM case.
    void CreateSurrogateBoundary(NurbsSurfaceGeometryPointerType& p_surface, ModelPart& rModelPart, ModelPart& rSurrogateModelPartInner, ModelPart& rSurrogateModelPartOuter, const Point& A_uvw, const Point& B_uvw)
    {
        CreateBrepSurface(p_surface, rModelPart, rSurrogateModelPartInner, rSurrogateModelPartOuter, mEchoLevel);
        CreateBrepCurveOnSurfaces(p_surface, rModelPart, rSurrogateModelPartInner, rSurrogateModelPartOuter, A_uvw, B_uvw, mEchoLevel);
    }

    /// Adds the surface geometry to the herin provided model_part and create the boundary breps when SBM is not needed.
    void CreateSurrogateBoundary(NurbsSurfaceGeometryPointerType& p_surface, ModelPart& rModelPart, const Point& A_uvw, const Point& B_uvw)
    {
        CreateBrepSurface(p_surface, rModelPart, mEchoLevel);
        int id_brep_curve_on_surface = 2;
        CreateBrepCurvesOnRectangle(rModelPart, p_surface, A_uvw, B_uvw, id_brep_curve_on_surface);
    }


private:

    /**
     * @brief Create a Brep Surface object
     * 
     * @param p_surface 
     * @param rModelPart 
     * @param EchoLevel 
     */
    static void CreateBrepSurface(
        NurbsSurfaceGeometryPointerType p_surface,
        ModelPart& rModelPart,
        SizeType EchoLevel = 0)
    {
        KRATOS_INFO_IF("ReadBrepSurface", (EchoLevel > 3))
            << "Creating BrepSurface \""<< std::endl;

        BrepCurveOnSurfaceLoopArrayType outer_loops, inner_loops;

        auto p_brep_surface =
            Kratos::make_shared<BrepSurfaceType>(
                p_surface, 
                outer_loops,
                inner_loops,
                false);

        // Sets the brep as geometry parent of the nurbs surface.
        p_surface->SetGeometryParent(p_brep_surface.get());
        p_brep_surface->SetId(1);
        rModelPart.AddGeometry(p_brep_surface);
    }

    /**
     * @brief Create a Brep Surface object for the SBM case
     * 
     * @param p_surface 
     * @param rModelPart 
     * @param rSurrogateModelPartInner 
     * @param rSurrogateModelPartOuter 
     * @param EchoLevel 
     */
    static void CreateBrepSurface(
        NurbsSurfaceGeometryPointerType p_surface,
        ModelPart& rModelPart,
        ModelPart& rSurrogateModelPartInner, 
        ModelPart& rSurrogateModelPartOuter,
        SizeType EchoLevel = 0)
    {
        KRATOS_INFO_IF("ReadBrepSurface", (EchoLevel > 3))
            << "Creating BrepSurface \""<< std::endl;

        BrepCurveOnSurfaceLoopArrayType outer_loops, inner_loops;

        auto p_brep_surface =
            Kratos::make_shared<BrepSurfaceType>(
                p_surface, 
                outer_loops,
                inner_loops,
                rSurrogateModelPartInner,
                rSurrogateModelPartOuter);

        // Sets the brep as geometry parent of the nurbs surface.
        p_surface->SetGeometryParent(p_brep_surface.get());
        p_brep_surface->SetId(1);
        rModelPart.AddGeometry(p_brep_surface);
    }


    static void CreateBrepCurveOnSurfaces(
        NurbsSurfaceGeometryPointerType p_surface,
        ModelPart& rModelPart,
        ModelPart& rSurrogateModelPartInner, 
        ModelPart& rSurrogateModelPartOuter,
        const Point& A_uvw, const Point& B_uvw,
        SizeType EchoLevel = 0) {
    
        // OUTER :
        // Each element in the surrogate_model_part represents a surrogate boundary loop. First "node" is the initial ID of the first surrogate node and
        // the second "node" is the last surrogate node of that loop. (We have done this in the case we have multiple surrogate boundaries and 1 model part)

        int id_brep_curve_on_surface = 2;

        if (rSurrogateModelPartOuter.Nodes().size() > 0) {
            int sizeSurrogateLoop_outer = rSurrogateModelPartOuter.Nodes().size();
            std::vector<double> surrogatecoord_x_outer(sizeSurrogateLoop_outer);
            std::vector<double> surrogatecoord_y_outer(sizeSurrogateLoop_outer);
            int countSurrogateLoop_outer = 0;
            for (auto i_node = rSurrogateModelPartOuter.NodesEnd()-1; i_node != rSurrogateModelPartOuter.NodesBegin()-1; i_node--) {
                surrogatecoord_x_outer[countSurrogateLoop_outer] = i_node->X();
                surrogatecoord_y_outer[countSurrogateLoop_outer] = i_node->Y();
                countSurrogateLoop_outer++;
            }
            // OUTER 
            std::vector<NurbsCurveGeometry<2, PointerVector<Point>>::Pointer> trimming_curves_GPT_outer;
            for (std::size_t i = 0; i < surrogatecoord_x_outer.size(); ++i) {
                Vector active_range_knot_vector = ZeroVector(2);
                
                Point::Pointer point1 = Kratos::make_shared<Point>(surrogatecoord_x_outer[i], surrogatecoord_y_outer[i], 0.0);
                Point::Pointer point2 = Kratos::make_shared<Point>(
                    surrogatecoord_x_outer[(i + 1) % surrogatecoord_x_outer.size()],  // Wrap around for the last point
                    surrogatecoord_y_outer[(i + 1) % surrogatecoord_y_outer.size()],  // Wrap around for the last point
                    0.0);
                // Compute the knot vector needed
                if (surrogatecoord_x_outer[(i + 1) % surrogatecoord_x_outer.size()]==surrogatecoord_x_outer[i]) {
                    active_range_knot_vector[0] = surrogatecoord_y_outer[i];
                    active_range_knot_vector[1] = surrogatecoord_y_outer[(i + 1) % surrogatecoord_y_outer.size()];
                }
                else {
                    active_range_knot_vector[0] = surrogatecoord_x_outer[i];
                    active_range_knot_vector[1] = surrogatecoord_x_outer[(i + 1) % surrogatecoord_x_outer.size()];
                }
                //// Order the active_range_knot_vector
                if (active_range_knot_vector[0] > active_range_knot_vector[1]) {
                    double temp = active_range_knot_vector[1];
                    active_range_knot_vector[1] = active_range_knot_vector[0] ;
                    active_range_knot_vector[0] = temp ;
                }
                NurbsCurveGeometry<2, PointerVector<Point>>::Pointer p_trimming_curve = CreateSingleBrep(point1, point2, active_range_knot_vector);
                trimming_curves_GPT_outer.push_back(p_trimming_curve);
            }

            BrepCurveOnSurfaceLoopType trimming_brep_curve_vector_outer(surrogatecoord_x_outer.size());

            for (std::size_t i = 0; i < trimming_curves_GPT_outer.size(); ++i) {
                Vector active_range_vector = ZeroVector(2);
                if (surrogatecoord_x_outer[(i + 1) % surrogatecoord_x_outer.size()]==surrogatecoord_x_outer[i]) {
                    active_range_vector[0] = surrogatecoord_y_outer[i];
                    active_range_vector[1] = surrogatecoord_y_outer[(i + 1) % surrogatecoord_x_outer.size()];
                }
                else {
                    active_range_vector[0] = surrogatecoord_x_outer[i];
                    active_range_vector[1] = surrogatecoord_x_outer[(i + 1) % surrogatecoord_x_outer.size()];
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
                    p_surface, trimming_curves_GPT_outer[i], brep_active_range, curve_direction);
                p_brep_curve_on_surface->SetId(id_brep_curve_on_surface);

                rModelPart.AddGeometry(p_brep_curve_on_surface);
                id_brep_curve_on_surface++;

            }

        } else {
            CreateBrepCurvesOnRectangle(rModelPart, p_surface, A_uvw, B_uvw, id_brep_curve_on_surface);
        }

        // INNER
        for (IndexType iel = 1; iel < rSurrogateModelPartInner.Elements().size()+1; iel++) {
            int firstSurrogateNodeId = rSurrogateModelPartInner.pGetElement(iel)->GetGeometry()[0].Id(); // Element 1 because is the only surrogate loop
            int lastSurrogateNodeId = rSurrogateModelPartInner.pGetElement(iel)->GetGeometry()[1].Id();  // Element 1 because is the only surrogate loop
            int sizeSurrogateLoop = lastSurrogateNodeId - firstSurrogateNodeId + 1;

            int countSurrogateLoop = 0;
            std::vector<double> surrogatecoord_x(sizeSurrogateLoop);
            std::vector<double> surrogatecoord_y(sizeSurrogateLoop);
            for (int id_node = firstSurrogateNodeId; id_node < lastSurrogateNodeId+1; id_node++) {
                surrogatecoord_x[countSurrogateLoop] = rSurrogateModelPartInner.GetNode(id_node).X();
                surrogatecoord_y[countSurrogateLoop] = rSurrogateModelPartInner.GetNode(id_node).Y();
                countSurrogateLoop++;
            }
            std::vector<NurbsCurveGeometry<2, PointerVector<Point>>::Pointer> trimming_curves_GPT;

            for (std::size_t i = 0; i < surrogatecoord_x.size(); ++i) {
                Vector active_range_knot_vector = ZeroVector(2);
                
                Point::Pointer point1 = Kratos::make_shared<Point>(surrogatecoord_x[i], surrogatecoord_y[i], 0.0);
                Point::Pointer point2 = Kratos::make_shared<Point>(
                    surrogatecoord_x[(i + 1) % surrogatecoord_x.size()],  // Wrap around for the last point
                    surrogatecoord_y[(i + 1) % surrogatecoord_y.size()],  // Wrap around for the last point
                    0.0);
                // Compute the knot vector needed
                if (surrogatecoord_x[(i + 1) % surrogatecoord_x.size()]==surrogatecoord_x[i]) {
                    active_range_knot_vector[0] = surrogatecoord_y[i];
                    active_range_knot_vector[1] = surrogatecoord_y[(i + 1) % surrogatecoord_y.size()];
                }
                else {
                    active_range_knot_vector[0] = surrogatecoord_x[i];
                    active_range_knot_vector[1] = surrogatecoord_x[(i + 1) % surrogatecoord_x.size()];
                }
                // Order the active_range_knot_vector
                if (active_range_knot_vector[0] > active_range_knot_vector[1]) {
                    double temp = active_range_knot_vector[1];
                    active_range_knot_vector[1] = active_range_knot_vector[0] ;
                    active_range_knot_vector[0] = temp ;
                }
                NurbsCurveGeometry<2, PointerVector<Point>>::Pointer p_trimming_curve = CreateSingleBrep(point1, point2, active_range_knot_vector);
                trimming_curves_GPT.push_back(p_trimming_curve);
            }

            BrepCurveOnSurfaceLoopType trimming_brep_curve_vector(surrogatecoord_x.size());


            for (std::size_t i = 0; i < trimming_curves_GPT.size(); ++i) {
                Vector active_range_vector = ZeroVector(2);
                if (surrogatecoord_x[(i + 1) % surrogatecoord_x.size()]==surrogatecoord_x[i]) {
                    active_range_vector[0] = surrogatecoord_y[i];
                    active_range_vector[1] = surrogatecoord_y[(i + 1) % surrogatecoord_x.size()];
                }
                else {
                    active_range_vector[0] = surrogatecoord_x[i];
                    active_range_vector[1] = surrogatecoord_x[(i + 1) % surrogatecoord_x.size()];
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
                    p_surface, trimming_curves_GPT[i], brep_active_range, curve_direction);
                p_brep_curve_on_surface->SetId(id_brep_curve_on_surface);
                rModelPart.AddGeometry(p_brep_curve_on_surface);
                id_brep_curve_on_surface++;
                trimming_brep_curve_vector[i] = p_brep_curve_on_surface ;

            }

        }
    } 
    


    ///@}
    ///@name Utility functions
    ///@{

    int mEchoLevel;

    static typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer CreateSingleBrep(
        Point::Pointer point1, Point::Pointer point2, Vector active_range_knot_vector)
        {
            // Create the data for the trimming curves
            PointerVector<Point> control_points;
            control_points.push_back(point1);
            control_points.push_back(point2);
            int polynomial_degree = 1;
            Vector knot_vector = ZeroVector(4) ;
            knot_vector[0] = active_range_knot_vector[0] ;
            knot_vector[1] = active_range_knot_vector[0] ;
            knot_vector[2] = active_range_knot_vector[1] ;
            knot_vector[3] = active_range_knot_vector[1] ;
            // Create the trimming curves
            typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer p_trimming_curve(
                new NurbsCurveGeometry<2, PointerVector<Point>>(
                    control_points,
                    polynomial_degree,
                    knot_vector));   
            return p_trimming_curve;
        }
    

    
    static void CreateBrepCurvesOnRectangle(ModelPart& r_model_part, NurbsSurfaceGeometryPointerType p_surface_geometry, const Point& A_uvw, const Point& B_uvw, int &last_geometry_id) {
        Vector knot_vector = ZeroVector(2);
        knot_vector[0] = 0.0;
        knot_vector[1] = std::abs(B_uvw[0] - A_uvw[0]);
        int p = 1;

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
        
        
        
        auto brep_curve_on_surface = BrepCurveOnSurface< PointerVector<NodeType>, true, PointerVector<Point>>(p_surface_geometry, p_curve_1);
        auto p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(p_surface_geometry, p_curve_1);
        p_brep_curve_on_surface->SetId(last_geometry_id);
        r_model_part.AddGeometry(p_brep_curve_on_surface);
        
        brep_curve_on_surface = BrepCurveOnSurface< PointerVector<NodeType>, true, PointerVector<Point>>(p_surface_geometry, p_curve_2);
        p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(p_surface_geometry, p_curve_2);
        p_brep_curve_on_surface->SetId(++last_geometry_id);
        r_model_part.AddGeometry(p_brep_curve_on_surface);

        brep_curve_on_surface = BrepCurveOnSurface< PointerVector<NodeType>, true, PointerVector<Point>>(p_surface_geometry, p_curve_3);
        p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(p_surface_geometry, p_curve_3);
        p_brep_curve_on_surface->SetId(++last_geometry_id);
        r_model_part.AddGeometry(p_brep_curve_on_surface);
        
        brep_curve_on_surface = BrepCurveOnSurface< PointerVector<NodeType>, true, PointerVector<Point>>(p_surface_geometry, p_curve_4);
        p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(p_surface_geometry, p_curve_4);
        p_brep_curve_on_surface->SetId(++last_geometry_id);
        r_model_part.AddGeometry(p_brep_curve_on_surface);

        last_geometry_id++;
    }

    ///@}
}; // Class CreateBrepsSBMUtilities
}  // namespace Kratos.
