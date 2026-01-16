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
#include "iga_application_variables.h"

// Geometries
#include "geometries/nurbs_curve_geometry.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/brep_surface.h"
#include "geometries/brep_volume.h"
#include "geometries/brep_curve.h"
#include "geometries/brep_curve_on_surface.h"
#include "geometries/brep_surface_on_volume.h"

#include "includes/io.h"

#include "geometries/quadrilateral_3d_4.h"

namespace Kratos
{

///@name Kratos Classes
///@{
template<class TNodeType = Node, class TEmbeddedNodeType = Point, bool TShiftedBoundary = true>
class CreateBrepsSbmUtilities : public IO
{
    public:

    ///@}
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CreateBrepsSbmUtilities
    KRATOS_CLASS_POINTER_DEFINITION(CreateBrepsSbmUtilities);

    using GeometryType = Geometry<TNodeType>;
    using GeometryPointerType = typename GeometryType::Pointer;
    using GeometrySurrogateArrayType = DenseVector<GeometryPointerType>;

    using ContainerNodeType = PointerVector<TNodeType>;
    using ContainerEmbeddedNodeType = PointerVector<TEmbeddedNodeType>;
    
    using NurbsSurfaceGeometryType = NurbsSurfaceGeometry<3, PointerVector<TNodeType>>;
    using NurbsSurfaceGeometryPointerType = typename NurbsSurfaceGeometryType::Pointer;
    using NurbsVolumeGeometryType = NurbsVolumeGeometry<PointerVector<TNodeType>>;
    using NurbsVolumeGeometryPointerType = typename NurbsVolumeGeometryType::Pointer;

    using BrepSurfaceType = BrepSurface<ContainerNodeType, TShiftedBoundary, ContainerEmbeddedNodeType>;
    using BrepCurveOnSurfaceType = BrepCurveOnSurface<ContainerNodeType, TShiftedBoundary, ContainerEmbeddedNodeType>;
    using BrepVolumeType = BrepVolume<ContainerNodeType, false, ContainerEmbeddedNodeType>;
    using BrepSurfaceOnVolumeType = BrepSurfaceOnVolume<ContainerNodeType, false, ContainerEmbeddedNodeType>;
    using BrepVolumeTypeSbm = BrepVolume<ContainerNodeType, true, ContainerEmbeddedNodeType>;
    using BrepSurfaceOnVolumeTypeSbm = BrepSurfaceOnVolume<ContainerNodeType, true, ContainerEmbeddedNodeType>;

    using BrepCurveOnSurfaceLoopArrayType = DenseVector<DenseVector<typename BrepCurveOnSurfaceType::Pointer>>;
    using BrepSurfaceOnVolumeLoopArrayType = DenseVector<DenseVector<typename BrepSurfaceOnVolumeType::Pointer>>;
    using BrepSurfaceOnVolumeLoopArrayTypeSbm = DenseVector<DenseVector<typename BrepSurfaceOnVolumeTypeSbm::Pointer>>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with path to input file.
    CreateBrepsSbmUtilities(
        std::size_t EchoLevel = 0)
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
        std::size_t id_brep_curve_on_surface = 2; // because id 1 is the brep surface
        CreateBrepCurvesOnRectangle(pSurface, rCoordsA, rCoordsB, id_brep_curve_on_surface, rModelPart);
    }

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
    void CreateSurrogateBoundary(NurbsVolumeGeometryPointerType& pVolume, 
        const ModelPart& rSurrogateModelPartInner, 
        const ModelPart& rSurrogateModelPartOuter, 
        const Point& rCoordsA, 
        const Point& rCoordsB,
        ModelPart& rModelPart)
    {
        CreateBrepVolume(pVolume, rSurrogateModelPartInner, rSurrogateModelPartOuter, rModelPart, mEchoLevel);
        CreateBrepSurfaceOnVolume(pVolume, rSurrogateModelPartInner, rSurrogateModelPartOuter, rCoordsA, rCoordsB, rModelPart, mEchoLevel);
    }

    /**
     * @brief Create a Surrogate Boundary object 3D
     * 
     * @param pVolume 
     * @param rCoordsA 
     * @param rCoordsB 
     * @param rModelPart 
     */
    void CreateSurrogateBoundary(
        NurbsVolumeGeometryPointerType& pVolume, 
        const Point& rCoordsA, 
        const Point& rCoordsB,
        ModelPart& rModelPart)
    {
        CreateBrepVolume(pVolume, rModelPart, mEchoLevel);
        IndexType id_brep_surface_on_volume = 2; // because id 1 is the brep volume
        CreateBrepSurfaceOnParallelepiped(pVolume, rCoordsA, rCoordsB, id_brep_surface_on_volume, rModelPart);
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
        const std::size_t EchoLevel = 0)
    {
        KRATOS_INFO_IF("ReadBrepSurface", (EchoLevel > 3))
            << "Creating BrepSurface \""<< std::endl;

        BrepCurveOnSurfaceLoopArrayType outer_loops, inner_loops;

        // Set to FALSE the TShiftedBoundary flag
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
        const std::size_t EchoLevel = 0)
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

        auto p_surrogate_outer_loop_geometries = Kratos::make_shared<GeometrySurrogateArrayType>(surrogate_outer_loop_geometries);
        auto p_surrogate_inner_loop_geometries = Kratos::make_shared<GeometrySurrogateArrayType>(surrogate_inner_loop_geometries);

        p_brep_surface->SetSurrogateOuterLoopGeometries(p_surrogate_outer_loop_geometries);
        p_brep_surface->SetSurrogateInnerLoopGeometries(p_surrogate_inner_loop_geometries);

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
        const std::size_t EchoLevel = 0) {
        // OUTER 
        std::size_t id_brep_curve_on_surface = 2; // because id 1 is the brep surface

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
                i_cond.SetValue(BREP_ID, static_cast<int>(id_brep_curve_on_surface));
                auto& r_geom_cond = i_cond.GetGeometry();
                if (r_geom_cond.size() >= 2) {
                    r_geom_cond[0].SetValue(BREP_ID, static_cast<int>(id_brep_curve_on_surface));
                    r_geom_cond[1].SetValue(BREP_ID, static_cast<int>(id_brep_curve_on_surface));
                }

                rModelPart.AddGeometry(p_brep_curve_on_surface);
                id_brep_curve_on_surface++;
            }
            
        } else {
            CreateBrepCurvesOnRectangle(pSurface, rCoordsA, rCoordsB, id_brep_curve_on_surface, rModelPart);
        }


        // INNER
        const auto& r_elems = rSurrogateModelPartInner.Elements();

        for (const auto& rElem : r_elems) {
            const auto& r_geom = rElem.GetGeometry();
            std::size_t first_surrogate_cond_id = r_geom[0].Id();
            std::size_t last_surrogate_cond_id  = r_geom[1].Id();

            for (std::size_t cond_id = first_surrogate_cond_id; cond_id <= last_surrogate_cond_id; ++cond_id) {
                
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
                                            std::size_t& rLastGeometryId,
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



    // 3D
    /**
     * @brief Create a Brep Surface object
     * 
     * @param pVolume 
     * @param rModelPart 
     * @param EchoLevel 
     */
    static void CreateBrepVolume(
        NurbsVolumeGeometryPointerType pVolume,
        ModelPart& rModelPart,
        const SizeType EchoLevel)
    {

        auto p_nurbs_cube = Kratos::make_shared<NurbsVolumeGeometry<PointerVector<NodeType>>>(*pVolume);
        
        KRATOS_INFO_IF("ReadBrepVolume", (EchoLevel > 3))
            << "Creating BrepVolume \""<< std::endl;

        BrepSurfaceOnVolumeLoopArrayType outer_loops, inner_loops;

        KRATOS_WATCH("we are setting false--> BodyFitted")
        auto p_brep_volume =
            Kratos::make_shared<BrepVolumeType>(
                pVolume, 
                outer_loops,
                inner_loops,
                false);

        // Sets the brep as geometry parent of the nurbs surface.
        pVolume->SetGeometryParent(p_brep_volume.get());
        p_brep_volume->SetId(1);
        rModelPart.AddGeometry(p_brep_volume);
    }


    /**
     * @brief Create a Brep Volume object
     * 
     * @param pSurface 
     * @param rSurrogateModelPartInner 
     * @param rSurrogateModelPartOuter 
     * @param rModelPart 
     * @param EchoLevel 
     */
    static void CreateBrepVolume(
        NurbsVolumeGeometryPointerType pVolume,
        const ModelPart& rSurrogateModelPartInner, 
        const ModelPart& rSurrogateModelPartOuter,
        ModelPart& rModelPart,
        const SizeType EchoLevel = 0)
    {
        KRATOS_INFO_IF("ReadBrepVolume ", (EchoLevel > 3))
            << "Creating BrepVolume \""<< std::endl;

        BrepSurfaceOnVolumeLoopArrayTypeSbm outer_loops, inner_loops;
        
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

        auto p_brep_volume =
            Kratos::make_shared<BrepVolumeTypeSbm>(
                pVolume, 
                outer_loops,
                inner_loops);

        auto p_surrogate_outer_loop_geometries = Kratos::make_shared<GeometrySurrogateArrayType>(surrogate_outer_loop_geometries);
        auto p_surrogate_inner_loop_geometries = Kratos::make_shared<GeometrySurrogateArrayType>(surrogate_inner_loop_geometries);

        p_brep_volume->SetSurrogateOuterLoopGeometries(p_surrogate_outer_loop_geometries);
        p_brep_volume->SetSurrogateInnerLoopGeometries(p_surrogate_inner_loop_geometries);

        // Sets the brep as geometry parent of the nurbs surface.
        pVolume->SetGeometryParent(p_brep_volume.get());
        p_brep_volume->SetId(1);
        rModelPart.AddGeometry(p_brep_volume);
    }


    /**
     * @brief Create a Brep Surface On Cube object
     * 
     * @param pVolumeGeometry 
     * @param rCoordsA 
     * @param rCoordsB 
     * @param rLastGeometryId 
     * @param rModelPart 
     */
    static void CreateBrepSurfaceOnParallelepiped(const NurbsVolumeGeometryPointerType pVolumeGeometry, 
            const Point& rCoordsA, 
            const Point& rCoordsB, 
            IndexType& rLastGeometryId,
            ModelPart& rModelPart) 
    {
        Vector knot_vector_u = ZeroVector(2); 
        Vector knot_vector_v = ZeroVector(2); 
        Vector normal = ZeroVector(3);
        const SizeType p = 1;
        const SizeType q = 1;
        
        // TO DO: DOES NOT WORK IF THE KNOT SPANS ARE RECTANGLES
        // lower
        Geometry<NodeType>::PointsArrayType points_lower;
        knot_vector_u[0] = 0.0; knot_vector_u[1] = std::abs(rCoordsA[0]-rCoordsB[0]);  
        knot_vector_v[0] = 0.0; knot_vector_v[1] = std::abs(rCoordsA[1]-rCoordsB[1]); 
        points_lower.push_back(NodeType::Pointer(new NodeType(1, rCoordsA[0], rCoordsA[1], rCoordsA[2])));
        points_lower.push_back(NodeType::Pointer(new NodeType(2, rCoordsA[0], rCoordsB[1], rCoordsA[2])));
        points_lower.push_back(NodeType::Pointer(new NodeType(3, rCoordsB[0], rCoordsA[1], rCoordsA[2])));
        points_lower.push_back(NodeType::Pointer(new NodeType(4, rCoordsB[0], rCoordsB[1], rCoordsA[2])));
        normal = ZeroVector(3);
        normal[2] = -1;
        auto p_surface_1 = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(points_lower, p, q, knot_vector_u, knot_vector_v);
        auto brep_p_surface_1 = Kratos::make_shared<BrepSurface<PointerVector<NodeType>,  false, PointerVector<Point>>>(p_surface_1);
        auto p_brep_surface_on_volume_1 = Kratos::make_shared<BrepSurfaceOnVolume< PointerVector<NodeType>, false, PointerVector<NodeType>>>(pVolumeGeometry, p_surface_1);
        p_brep_surface_on_volume_1->SetId(rLastGeometryId);
        p_brep_surface_on_volume_1->SetNormalSbm(normal);
        rModelPart.AddGeometry(p_brep_surface_on_volume_1);

        // upper
        Geometry<NodeType>::PointsArrayType points_upper;
        knot_vector_u[0] = 0.0; knot_vector_u[1] = std::abs(rCoordsA[0]-rCoordsB[0]);  
        knot_vector_v[0] = 0.0; knot_vector_v[1] = std::abs(rCoordsA[1]-rCoordsB[1]); 
        points_upper.push_back(NodeType::Pointer(new NodeType(1, rCoordsA[0], rCoordsA[1], rCoordsB[2])));
        points_upper.push_back(NodeType::Pointer(new NodeType(2, rCoordsA[0], rCoordsB[1], rCoordsB[2])));
        points_upper.push_back(NodeType::Pointer(new NodeType(3, rCoordsB[0], rCoordsA[1], rCoordsB[2])));
        points_upper.push_back(NodeType::Pointer(new NodeType(4, rCoordsB[0], rCoordsB[1], rCoordsB[2])));
        normal = ZeroVector(3);
        normal[2] = 1;
        auto p_surface_2 = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(points_upper, p, q, knot_vector_u, knot_vector_v);
        auto brep_p_surface_2 = Kratos::make_shared<BrepSurface<PointerVector<NodeType>,  false, PointerVector<Point>>>(p_surface_2);
        auto p_brep_surface_on_volume_2 = Kratos::make_shared<BrepSurfaceOnVolume< PointerVector<NodeType>, false, PointerVector<NodeType>>>(pVolumeGeometry, p_surface_2);
        p_brep_surface_on_volume_2->SetId(++rLastGeometryId);
        p_brep_surface_on_volume_2->SetNormalSbm(normal);
        rModelPart.AddGeometry(p_brep_surface_on_volume_2);

        // front
        Geometry<NodeType>::PointsArrayType points_front;
        knot_vector_u[0] = 0.0; knot_vector_u[1] = std::abs(rCoordsA[0]-rCoordsB[0]);  
        knot_vector_v[0] = 0.0; knot_vector_v[1] = std::abs(rCoordsA[2]-rCoordsB[2]); 
        points_front.push_back(NodeType::Pointer(new NodeType(1, rCoordsA[0], rCoordsA[1], rCoordsA[2])));
        points_front.push_back(NodeType::Pointer(new NodeType(2, rCoordsB[0], rCoordsA[1], rCoordsA[2])));
        points_front.push_back(NodeType::Pointer(new NodeType(3, rCoordsA[0], rCoordsA[1], rCoordsB[2]))); 
        points_front.push_back(NodeType::Pointer(new NodeType(4, rCoordsB[0], rCoordsA[1], rCoordsB[2])));
        normal = ZeroVector(3);
        normal[1] = -1;
        auto p_surface_3 = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(points_front, p, q, knot_vector_u, knot_vector_v);
        auto brep_p_surface_3 = Kratos::make_shared<BrepSurface<PointerVector<NodeType>,  false, PointerVector<Point>>>(p_surface_3);
        auto p_brep_surface_on_volume_3 = Kratos::make_shared<BrepSurfaceOnVolume< PointerVector<NodeType>, false, PointerVector<NodeType>>>(pVolumeGeometry, p_surface_3);
        p_brep_surface_on_volume_3->SetId(++rLastGeometryId);
        p_brep_surface_on_volume_3->SetNormalSbm(normal);
        rModelPart.AddGeometry(p_brep_surface_on_volume_3);

        // back
        Geometry<NodeType>::PointsArrayType points_back;
        knot_vector_u[0] = 0.0; knot_vector_u[1] = std::abs(rCoordsA[0]-rCoordsB[0]);  
        knot_vector_v[0] = 0.0; knot_vector_v[1] = std::abs(rCoordsA[2]-rCoordsB[2]); 
        points_back.push_back(NodeType::Pointer(new NodeType(1, rCoordsA[0], rCoordsB[1], rCoordsA[2])));
        points_back.push_back(NodeType::Pointer(new NodeType(2, rCoordsA[0], rCoordsB[1], rCoordsB[2])));
        points_back.push_back(NodeType::Pointer(new NodeType(3, rCoordsB[0], rCoordsB[1], rCoordsA[2])));
        points_back.push_back(NodeType::Pointer(new NodeType(4, rCoordsB[0], rCoordsB[1], rCoordsB[2])));
        normal = ZeroVector(3);
        normal[1] = 1;
        auto p_surface_4 = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(points_back, p, q, knot_vector_u, knot_vector_v);
        auto brep_p_surface_4 = Kratos::make_shared<BrepSurface<PointerVector<NodeType>,  false, PointerVector<Point>>>(p_surface_4);
        auto p_brep_surface_on_volume_4 = Kratos::make_shared<BrepSurfaceOnVolume< PointerVector<NodeType>, false, PointerVector<NodeType>>>(pVolumeGeometry, p_surface_4);
        p_brep_surface_on_volume_4->SetId(++rLastGeometryId);
        p_brep_surface_on_volume_4->SetNormalSbm(normal);
        rModelPart.AddGeometry(p_brep_surface_on_volume_4);

        // left
        Geometry<NodeType>::PointsArrayType points_left;
        knot_vector_u[0] = 0.0; knot_vector_u[1] = std::abs(rCoordsA[1]-rCoordsB[1]);  
        knot_vector_v[0] = 0.0; knot_vector_v[1] = std::abs(rCoordsA[2]-rCoordsB[2]); 
        points_left.push_back(NodeType::Pointer(new NodeType(1, rCoordsA[0], rCoordsA[1], rCoordsA[2])));
        points_left.push_back(NodeType::Pointer(new NodeType(2, rCoordsA[0], rCoordsB[1], rCoordsA[2])));
        points_left.push_back(NodeType::Pointer(new NodeType(3, rCoordsA[0], rCoordsA[1], rCoordsB[2])));
        points_left.push_back(NodeType::Pointer(new NodeType(4, rCoordsA[0], rCoordsB[1], rCoordsB[2])));
        normal = ZeroVector(3);
        normal[0] = -1;
        auto p_surface_5 = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(points_left, p, q, knot_vector_u, knot_vector_v);
        auto brep_p_surface_5 = Kratos::make_shared<BrepSurface<PointerVector<NodeType>,  false, PointerVector<Point>>>(p_surface_5);
        auto p_brep_surface_on_volume_5 = Kratos::make_shared<BrepSurfaceOnVolume< PointerVector<NodeType>, false, PointerVector<NodeType>>>(pVolumeGeometry, p_surface_5);
        // auto p_brep_surface_on_volume = Kratos::make_shared<BrepSurfaceOnVolumeType>(pVolumeGeometry, brep_p_surface_1);
        p_brep_surface_on_volume_5->SetId(++rLastGeometryId);
        p_brep_surface_on_volume_5->SetNormalSbm(normal);
        rModelPart.AddGeometry(p_brep_surface_on_volume_5);

        // right
        Geometry<NodeType>::PointsArrayType points_right;
        knot_vector_u[0] = 0.0; knot_vector_u[1] = std::abs(rCoordsA[1]-rCoordsB[1]);  
        knot_vector_v[0] = 0.0; knot_vector_v[1] = std::abs(rCoordsA[2]-rCoordsB[2]); 
        points_right.push_back(NodeType::Pointer(new NodeType(1, rCoordsB[0], rCoordsA[1], rCoordsA[2])));
        points_right.push_back(NodeType::Pointer(new NodeType(2 ,rCoordsB[0], rCoordsA[1], rCoordsB[2])));
        points_right.push_back(NodeType::Pointer(new NodeType(3 ,rCoordsB[0], rCoordsB[1], rCoordsA[2])));
        points_right.push_back(NodeType::Pointer(new NodeType(4 ,rCoordsB[0], rCoordsB[1], rCoordsB[2])));
        normal = ZeroVector(3);
        normal[0] = 1;
        auto p_surface_6 = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(points_right, p, q, knot_vector_u, knot_vector_v);
        auto brep_p_surface_6 = Kratos::make_shared<BrepSurface<PointerVector<NodeType>,  false, PointerVector<Point>>>(p_surface_6);
        auto p_brep_surface_on_volume_6 = Kratos::make_shared<BrepSurfaceOnVolume< PointerVector<NodeType>, false, PointerVector<NodeType>>>(pVolumeGeometry, p_surface_6);
        p_brep_surface_on_volume_6->SetId(++rLastGeometryId);
        p_brep_surface_on_volume_6->SetNormalSbm(normal);
        rModelPart.AddGeometry(p_brep_surface_on_volume_6);

        rLastGeometryId++;
    }


    /**
     * @brief Create a Brep Surface On Volume object
     * 
     * @param pVolume 
     * @param rSurrogateModelPartInner 
     * @param rSurrogateModelPartOuter 
     * @param rCoordsA 
     * @param rCoordsB 
     * @param rModelPart 
     * @param EchoLevel 
     */
    static void CreateBrepSurfaceOnVolume(
        const NurbsVolumeGeometryPointerType pVolume,
        const ModelPart& rSurrogateModelPartInner, 
        const ModelPart& rSurrogateModelPartOuter,
        const Point& rCoordsA, 
        const Point& rCoordsB,
        ModelPart& rModelPart,
        const SizeType EchoLevel = 0) 
    {

        // OUTER 
        IndexType id_brep_surface_on_volume = 2; // because id 1 is the brep surface

        if (rSurrogateModelPartOuter.NumberOfConditions() > 0)
        {
            Vector normal = ZeroVector(3);
            for (auto &i_cond : rSurrogateModelPartOuter.Conditions()) {
                Geometry<Node>::PointsArrayType points;
                normal = ZeroVector(3);

                Point A_uvw_sbm = i_cond.GetGeometry()[0];
                Point B_uvw_sbm = i_cond.GetGeometry()[2];
                Point C_uvw_sbm = i_cond.GetGeometry()[1];
                Point D_uvw_sbm = i_cond.GetGeometry()[3];

                Vector knot_vector_u = ZeroVector(2); 
                knot_vector_u[0] = 0.0;
                Vector knot_vector_v = ZeroVector(2); 
                knot_vector_v[0] = 0.0;

                // Check the orientation of the sbm face
                array_1d<double, 3> diagonalAB(3);
                diagonalAB[0] = B_uvw_sbm[0] - A_uvw_sbm[0];
                diagonalAB[1] = B_uvw_sbm[1] - A_uvw_sbm[1];
                diagonalAB[2] = B_uvw_sbm[2] - A_uvw_sbm[2];
                array_1d<double, 3> x_unit({1.0, 0.0, 0.0});
                array_1d<double, 3> y_unit({0.0, 1.0, 0.0});
                array_1d<double, 3> z_unit({0.0, 0.0, 1.0});
                int perpendicular_direction = -1;
                if (std::abs(inner_prod(diagonalAB, x_unit)) < 1e-13 ) {
                    // the normal is parallel to x
                    knot_vector_u[1] = std::abs(A_uvw_sbm[1]-B_uvw_sbm[1]);
                    knot_vector_v[1] = std::abs(A_uvw_sbm[2]-B_uvw_sbm[2]);
                    perpendicular_direction = 0;
                    normal[0] = 1; 
                }
                else if (std::abs(inner_prod(diagonalAB, y_unit)) < 1e-13 ){
                    // the normal is parallel to y
                    knot_vector_u[1] = std::abs(A_uvw_sbm[0]-B_uvw_sbm[0]);
                    knot_vector_v[1] = std::abs(A_uvw_sbm[2]-B_uvw_sbm[2]);
                    perpendicular_direction = 1;
                    normal[1] = 1; 
                }
                else if (std::abs(inner_prod(diagonalAB, z_unit)) < 1e-13 ){
                    // the normal is parallel to z
                    knot_vector_u[1] = std::abs(A_uvw_sbm[0]-B_uvw_sbm[0]);
                    knot_vector_v[1] = std::abs(A_uvw_sbm[1]-B_uvw_sbm[1]);
                    perpendicular_direction = 2;
                    normal[2] = 1; 
                } else {
                    KRATOS_ERROR << "Surrogate face not parallel to any unit direction" << std::endl;
                }
                
                const SizeType p = 1;
                const SizeType q = 1;

                points.push_back(NodeType::Pointer(new NodeType(1, A_uvw_sbm[0], A_uvw_sbm[1], A_uvw_sbm[2])));
                points.push_back(NodeType::Pointer(new NodeType(2, D_uvw_sbm[0], D_uvw_sbm[1], D_uvw_sbm[2])));
                points.push_back(NodeType::Pointer(new NodeType(3, C_uvw_sbm[0], C_uvw_sbm[1], C_uvw_sbm[2])));
                points.push_back(NodeType::Pointer(new NodeType(4, B_uvw_sbm[0], B_uvw_sbm[1], B_uvw_sbm[2])));
                
                auto p_surface_sbm = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(points, p, q, knot_vector_u, knot_vector_v);
                auto brep_p_surface_sbm = Kratos::make_shared<BrepSurface<PointerVector<NodeType>, true, PointerVector<Point>>>(p_surface_sbm);
                auto p_brep_surface_on_volume_sbm = Kratos::make_shared<BrepSurfaceOnVolume< PointerVector<NodeType>, true, PointerVector<NodeType>>>(pVolume, p_surface_sbm);
                p_brep_surface_on_volume_sbm->SetId(id_brep_surface_on_volume);
                
                // Set if the surface on volume is entering of exiting. 
                // Note we save BOUNDARY true to the condition of the surrogate model part when in entering.
                bool isExiting = !i_cond.Is(BOUNDARY);
                int is_exiting_direction = -1;
                /* is_exiting_direction :
                                -1 -> entering surface/external body fitted surface
                                0 -> exiting surface in x direction
                                1 -> exiting surface in y direction
                                2 -> exiting surface in z direction
                                */
                if (!isExiting) {
                    is_exiting_direction = perpendicular_direction;
                    normal = normal * (-1);
                }
                p_brep_surface_on_volume_sbm->SetIsExitingDirectionSbm(is_exiting_direction);
                // KRATOS_WATCH("\n")
                // KRATOS_WATCH(is_exiting_direction)
                // KRATOS_WATCH(isExiting)
                // KRATOS_WATCH(i_cond)
                // KRATOS_WATCH(normal)
                p_brep_surface_on_volume_sbm->SetNormalSbm(normal);
                
                rModelPart.AddGeometry(p_brep_surface_on_volume_sbm);
                id_brep_surface_on_volume++;

            }
            
        } else {
            CreateBrepSurfaceOnParallelepiped(pVolume, rCoordsA, rCoordsB, id_brep_surface_on_volume, rModelPart);
        }


        // INNER
        for (IndexType i_element = 1; i_element < rSurrogateModelPartInner.Elements().size()+1; i_element++) {
            IndexType first_surrogate_cond_id = rSurrogateModelPartInner.pGetElement(i_element)->GetGeometry()[0].Id(); // Element 1 because is the only surrogate loop
            IndexType last_surrogate_cond_id = rSurrogateModelPartInner.pGetElement(i_element)->GetGeometry()[1].Id();  // Element 1 because is the only surrogate loop
            
            Vector normal = ZeroVector(3);
            for (IndexType cond_id = first_surrogate_cond_id; cond_id <= last_surrogate_cond_id; ++cond_id) {
                
                auto p_cond = rSurrogateModelPartInner.pGetCondition(cond_id);

                Geometry<Node>::PointsArrayType points;
                normal = ZeroVector(3);

                Point A_uvw_sbm = p_cond->GetGeometry()[0];
                Point B_uvw_sbm = p_cond->GetGeometry()[2];
                Point C_uvw_sbm = p_cond->GetGeometry()[1];
                Point D_uvw_sbm = p_cond->GetGeometry()[3];

                Vector knot_vector_u = ZeroVector(2); 
                knot_vector_u[0] = 0.0;
                Vector knot_vector_v = ZeroVector(2); 
                knot_vector_v[0] = 0.0;

                // Check the orientation of the sbm face
                array_1d<double, 3> diagonalAB(3);
                diagonalAB[0] = B_uvw_sbm[0] - A_uvw_sbm[0];
                diagonalAB[1] = B_uvw_sbm[1] - A_uvw_sbm[1];
                diagonalAB[2] = B_uvw_sbm[2] - A_uvw_sbm[2];
                array_1d<double, 3> x_unit({1.0, 0.0, 0.0});
                array_1d<double, 3> y_unit({0.0, 1.0, 0.0});
                array_1d<double, 3> z_unit({0.0, 0.0, 1.0});
                int perpendicular_direction = -1;
                if (std::abs(inner_prod(diagonalAB, x_unit)) < 1e-13 ) {
                    // the normal is parallel to x
                    knot_vector_u[1] = std::abs(A_uvw_sbm[1]-B_uvw_sbm[1]);
                    knot_vector_v[1] = std::abs(A_uvw_sbm[2]-B_uvw_sbm[2]);
                    perpendicular_direction = 0;
                    normal[0] = 1; 
                }
                else if (std::abs(inner_prod(diagonalAB, y_unit)) < 1e-13 ){
                    // the normal is parallel to y
                    knot_vector_u[1] = std::abs(A_uvw_sbm[0]-B_uvw_sbm[0]);
                    knot_vector_v[1] = std::abs(A_uvw_sbm[2]-B_uvw_sbm[2]);
                    perpendicular_direction = 1;
                    normal[1] = 1; 
                }
                else if (std::abs(inner_prod(diagonalAB, z_unit)) < 1e-13 ){
                    // the normal is parallel to z
                    knot_vector_u[1] = std::abs(A_uvw_sbm[0]-B_uvw_sbm[0]);
                    knot_vector_v[1] = std::abs(A_uvw_sbm[1]-B_uvw_sbm[1]);
                    perpendicular_direction = 2;
                    normal[2] = 1; 
                } else {
                    KRATOS_ERROR << "Surrogate face not parallel to any unit direction" << std::endl;
                }
                
                const SizeType p = 1;
                const SizeType q = 1;

                points.push_back(NodeType::Pointer(new NodeType(1, A_uvw_sbm[0], A_uvw_sbm[1], A_uvw_sbm[2])));
                points.push_back(NodeType::Pointer(new NodeType(2, D_uvw_sbm[0], D_uvw_sbm[1], D_uvw_sbm[2])));
                points.push_back(NodeType::Pointer(new NodeType(3, C_uvw_sbm[0], C_uvw_sbm[1], C_uvw_sbm[2])));
                points.push_back(NodeType::Pointer(new NodeType(4, B_uvw_sbm[0], B_uvw_sbm[1], B_uvw_sbm[2])));
                
                auto p_surface_sbm = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(points, p, q, knot_vector_u, knot_vector_v);
                auto brep_p_surface_sbm = Kratos::make_shared<BrepSurface<PointerVector<NodeType>, true, PointerVector<Point>>>(p_surface_sbm);
                auto p_brep_surface_on_volume_sbm = Kratos::make_shared<BrepSurfaceOnVolume< PointerVector<NodeType>, true, PointerVector<NodeType>>>(pVolume, p_surface_sbm);
                p_brep_surface_on_volume_sbm->SetId(id_brep_surface_on_volume);
                
                // Set if the surface on volume is entering of exiting. 
                // Note we save BOUNDARY true to the condition of the surrogate model part when in entering.
                bool isExiting = !p_cond->Is(BOUNDARY);
                int is_exiting_direction = -1;
                /* is_exiting_direction :
                                -1 -> entering surface/external body fitted surface
                                0 -> exiting surface in x direction
                                1 -> exiting surface in y direction
                                2 -> exiting surface in z direction
                                */
                if (isExiting) {
                    is_exiting_direction = perpendicular_direction;
                    normal = normal * (-1);
                }
                p_brep_surface_on_volume_sbm->SetIsExitingDirectionSbm(is_exiting_direction);
                p_brep_surface_on_volume_sbm->SetNormalSbm(normal);
                
                rModelPart.AddGeometry(p_brep_surface_on_volume_sbm);
                id_brep_surface_on_volume++;
            }
        }
        
    } 
    
    ///@}
    ///@name Utility functions
    ///@{

    int mEchoLevel;

    ///@}
}; // Class CreateBrepsSbmUtilities
}  // namespace Kratos.
