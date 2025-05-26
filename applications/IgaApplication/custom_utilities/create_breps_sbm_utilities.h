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
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

// Geometries
#include "geometries/nurbs_curve_geometry.h"
#include "geometries/brep_surface.h"
#include "geometries/brep_curve.h"
#include "geometries/brep_curve_on_surface.h"
#include "includes/io.h"

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
    
    using NurbsSurfaceGeometryType = NurbsSurfaceGeometry<3, PointerVector<TNodeType>>;
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
     * @param rSurrogateModelPartInner 
     * @param rSurrogateModelPartOuter 
     * @param rCoordsA 
     * @param rCoordsB 
     * @param rModelPart 
     */
    void CreateSurrogateBoundary(NurbsSurfaceGeometryPointerType& pSurface, 
                                 const ModelPart& rSurrogateModelPartInner, 
                                 const ModelPart& rSurrogateModelPartOuter, 
                                 const Point& rCoordsA, 
                                 const Point& rCoordsB,
                                 ModelPart& rModelPart)
    {
        CreateBrepSurface(pSurface, rSurrogateModelPartInner, rSurrogateModelPartOuter, rModelPart, mEchoLevel);
        CreateBrepCurveOnSurfaces(pSurface, rSurrogateModelPartInner, rSurrogateModelPartOuter, rCoordsA, rCoordsB, rModelPart, mEchoLevel);
    }
    
    /**
     * @brief Create a Surrogate Boundary object
     * 
     * @param pSurface 
     * @param rCoordsA 
     * @param rCoordsB 
     * @param rModelPart 
     */
    void CreateSurrogateBoundary(NurbsSurfaceGeometryPointerType& pSurface, 
                                const Point& rCoordsA, 
                                const Point& rCoordsB,
                                ModelPart& rModelPart)
    {
        CreateBrepSurface(pSurface, rModelPart, mEchoLevel);
        IndexType id_brep_curve_on_surface = 2; // because id 1 is the brep surface
        CreateBrepCurvesOnRectangle(pSurface, rCoordsA, rCoordsB, id_brep_curve_on_surface, rModelPart);
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
        const SizeType EchoLevel = 0)
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
     * @brief Create a Brep Surface object
     * 
     * @param pSurface 
     * @param rSurrogateModelPartInner 
     * @param rSurrogateModelPartOuter 
     * @param rModelPart 
     * @param EchoLevel 
     */
    static void CreateBrepSurface(
        NurbsSurfaceGeometryPointerType pSurface,
        const ModelPart& rSurrogateModelPartInner, 
        const ModelPart& rSurrogateModelPartOuter,
        ModelPart& rModelPart,
        const SizeType EchoLevel = 0)
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
     * @param rSurrogateModelPartInner 
     * @param rSurrogateModelPartOuter 
     * @param rCoordsA 
     * @param rCoordsB 
     * @param rModelPart 
     * @param EchoLevel 
     */
    static void CreateBrepCurveOnSurfaces(
        const NurbsSurfaceGeometryPointerType pSurface,
        const ModelPart& rSurrogateModelPartInner, 
        const ModelPart& rSurrogateModelPartOuter,
        const Point& rCoordsA, 
        const Point& rCoordsB,
        ModelPart& rModelPart,
        const SizeType EchoLevel = 0) {
        // OUTER 
        IndexType id_brep_curve_on_surface = 2; // because id 1 is the brep surface

        if (rSurrogateModelPartOuter.NumberOfConditions() > 0)
        {
            for (auto &i_cond : rSurrogateModelPartOuter.Conditions()) {
                Geometry<Node>::PointsArrayType points;
                Point first_point  = i_cond.GetGeometry()[0];
                Point second_point = i_cond.GetGeometry()[1];
                Vector knot_vector = ZeroVector(2);

                // define active range: the first point and second point are in order
                Vector active_range_knot_vector = ZeroVector(2);
                // Compute the knot vector needed
                if (first_point[0] == second_point[0]) {
                    // the brep curve is vertical
                    active_range_knot_vector[0] = first_point[1];
                    active_range_knot_vector[1] = second_point[1];
                } else {
                    // the brep curve is horizontal
                    active_range_knot_vector[0] = first_point[0];
                    active_range_knot_vector[1] = second_point[0];
                }
                // Always sort the range, regardless of direction
                std::sort(active_range_knot_vector.begin(), active_range_knot_vector.end());

                // check if the brep is entering or exiting
                if (i_cond.Is(BOUNDARY))
                {
                    // if true is exiting -> reverse the order of the 
                    Point temp = first_point;
                    first_point = second_point;
                    second_point = temp;
                }
                
                Point::Pointer first_brep_point = Kratos::make_shared<Point>(first_point[0], first_point[1], 0.0);
                Point::Pointer p_second_brep_point = Kratos::make_shared<Point>(second_point[0], second_point[1], 0.0);
                
                NurbsCurveGeometry<2, PointerVector<Point>>::Pointer p_trimming_curve = CreateBrepCurve(first_brep_point, p_second_brep_point, active_range_knot_vector);

                NurbsInterval brep_active_range(active_range_knot_vector[0], active_range_knot_vector[1]);

                bool curve_direction = true;
                auto p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(
                    pSurface, p_trimming_curve, brep_active_range, curve_direction);

                p_brep_curve_on_surface->SetId(id_brep_curve_on_surface);

                rModelPart.AddGeometry(p_brep_curve_on_surface);
                id_brep_curve_on_surface++;
            }
            
        } else {
            CreateBrepCurvesOnRectangle(pSurface, rCoordsA, rCoordsB, id_brep_curve_on_surface, rModelPart);
        }


        // INNER
        for (IndexType i_element = 1; i_element < rSurrogateModelPartInner.Elements().size()+1; i_element++) {
            IndexType first_surrogate_cond_id = rSurrogateModelPartInner.pGetElement(i_element)->GetGeometry()[0].Id(); // Element 1 because is the only surrogate loop
            IndexType last_surrogate_cond_id = rSurrogateModelPartInner.pGetElement(i_element)->GetGeometry()[1].Id();  // Element 1 because is the only surrogate loop
    
            for (IndexType cond_id = first_surrogate_cond_id; cond_id <= last_surrogate_cond_id; ++cond_id) {
                
                auto p_cond = rSurrogateModelPartInner.pGetCondition(cond_id);
                Geometry<Node>::PointsArrayType points;
                Point first_point  = p_cond->GetGeometry()[0];
                Point second_point = p_cond->GetGeometry()[1];
                Vector knot_vector = ZeroVector(2);

                // define active range: the first point and second point are in order
                Vector active_range_knot_vector = ZeroVector(2);
                // Compute the knot vector needed
                if (first_point[0] == second_point[0]) {
                    // the brep curve is vertical
                    active_range_knot_vector[0] = first_point[1];
                    active_range_knot_vector[1] = second_point[1];
                } else {
                    // the brep curve is horizontal
                    active_range_knot_vector[0] = first_point[0];
                    active_range_knot_vector[1] = second_point[0];
                }
                // Always sort the range, regardless of direction
                std::sort(active_range_knot_vector.begin(), active_range_knot_vector.end());

                // check if the brep is entering or exiting
                if (p_cond->Is(BOUNDARY))
                {
                    // if true is exiting -> reverse the order of the 
                    Point temp = first_point;
                    first_point = second_point;
                    second_point = temp;
                }
                
                Point::Pointer p_first_brep_point = Kratos::make_shared<Point>(first_point[0], first_point[1], 0.0);
                Point::Pointer p_second_brep_point = Kratos::make_shared<Point>(second_point[0], second_point[1], 0.0);
                NurbsCurveGeometry<2, PointerVector<Point>>::Pointer p_trimming_curve = CreateBrepCurve(p_first_brep_point, p_second_brep_point, active_range_knot_vector);
                NurbsInterval brep_active_range(active_range_knot_vector[0], active_range_knot_vector[1]);

                bool curve_direction = true;
                auto p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(
                    pSurface, p_trimming_curve, brep_active_range, curve_direction);

                p_brep_curve_on_surface->SetId(id_brep_curve_on_surface);
                rModelPart.AddGeometry(p_brep_curve_on_surface);
                id_brep_curve_on_surface++;
            }
        }
    } 
    
    /**
     * @brief Create a Single Brep object
     * 
     * @param pFirstBrepPoint 
     * @param pSecondBrepPoint 
     * @param rActiveRangeKnotVector 
     * @return NurbsCurveGeometry<2, PointerVector<Point>>::Pointer 
     */
    static typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer CreateBrepCurve(
        const Point::Pointer pFirstBrepPoint,
        const Point::Pointer pSecondBrepPoint, 
        const Vector& rActiveRangeKnotVector)
        {
            // Create the data for the trimming curves
            PointerVector<Point> control_points;
            control_points.push_back(pFirstBrepPoint);
            control_points.push_back(pSecondBrepPoint);
            const int polynomial_degree = 1;
            Vector knot_vector = ZeroVector(4) ;
            knot_vector[0] = rActiveRangeKnotVector[0] ;
            knot_vector[1] = rActiveRangeKnotVector[0] ;
            knot_vector[2] = rActiveRangeKnotVector[1] ;
            knot_vector[3] = rActiveRangeKnotVector[1] ;
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
     * @param pSurfaceGeometry 
     * @param rCoordsA 
     * @param rCoordsB 
     * @param rLastGeometryId 
     * @param rModelPart 
     */
    static void CreateBrepCurvesOnRectangle(const NurbsSurfaceGeometryPointerType pSurfaceGeometry, 
                                            const Point& rCoordsA, 
                                            const Point& rCoordsB, 
                                            IndexType& rLastGeometryId,
                                            ModelPart& rModelPart) {
        Vector knot_vector = ZeroVector(2);
        knot_vector[0] = 0.0;
        knot_vector[1] = std::abs(rCoordsB[0] - rCoordsA[0]);
        const int p = 1;

        NurbsCurveGeometry<2, PointerVector<Point>>::PointsArrayType segment1;
        segment1.push_back(Point::Pointer(new Point(rCoordsA[0], rCoordsA[1])));
        segment1.push_back(Point::Pointer(new Point(rCoordsB[0], rCoordsA[1])));

        NurbsCurveGeometry<2, PointerVector<Point>>::PointsArrayType segment2;
        segment2.push_back(Point::Pointer(new Point(rCoordsB[0], rCoordsA[1])));
        segment2.push_back(Point::Pointer(new Point(rCoordsB[0], rCoordsB[1])));
        
        NurbsCurveGeometry<2, PointerVector<Point>>::PointsArrayType segment3;
        segment3.push_back(Point::Pointer(new Point(rCoordsB[0], rCoordsB[1])));
        segment3.push_back(Point::Pointer(new Point(rCoordsA[0], rCoordsB[1])));
        
        NurbsCurveGeometry<2, PointerVector<Point>>::PointsArrayType segment4;
        segment4.push_back(Point::Pointer(new Point(rCoordsA[0], rCoordsB[1])));
        segment4.push_back(Point::Pointer(new Point(rCoordsA[0], rCoordsA[1])));

        auto p_curve_1 = Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(segment1, p, knot_vector);
        auto p_curve_2 = Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(segment2, p, knot_vector);
        auto p_curve_3 = Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(segment3, p, knot_vector);
        auto p_curve_4 = Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(segment4, p, knot_vector);
        
        auto brep_curve_on_surface = BrepCurveOnSurface< PointerVector<TNodeType>, true, PointerVector<Point>>(pSurfaceGeometry, p_curve_1);
        auto p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(pSurfaceGeometry, p_curve_1);
        p_brep_curve_on_surface->SetId(rLastGeometryId);
        rModelPart.AddGeometry(p_brep_curve_on_surface);
        
        brep_curve_on_surface = BrepCurveOnSurface< PointerVector<TNodeType>, true, PointerVector<Point>>(pSurfaceGeometry, p_curve_2);
        p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(pSurfaceGeometry, p_curve_2);
        p_brep_curve_on_surface->SetId(++rLastGeometryId);
        rModelPart.AddGeometry(p_brep_curve_on_surface);

        brep_curve_on_surface = BrepCurveOnSurface< PointerVector<TNodeType>, true, PointerVector<Point>>(pSurfaceGeometry, p_curve_3);
        p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(pSurfaceGeometry, p_curve_3);
        p_brep_curve_on_surface->SetId(++rLastGeometryId);
        rModelPart.AddGeometry(p_brep_curve_on_surface);
        
        brep_curve_on_surface = BrepCurveOnSurface< PointerVector<TNodeType>, true, PointerVector<Point>>(pSurfaceGeometry, p_curve_4);
        p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(pSurfaceGeometry, p_curve_4);
        p_brep_curve_on_surface->SetId(++rLastGeometryId);
        rModelPart.AddGeometry(p_brep_curve_on_surface);

        rLastGeometryId++;
    }
    
    ///@}
    ///@name Utility functions
    ///@{

    int mEchoLevel;

    ///@}
}; // Class CreateBrepsSbmUtilities
}  // namespace Kratos.
