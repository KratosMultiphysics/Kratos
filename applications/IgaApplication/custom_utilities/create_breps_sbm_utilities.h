//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_CREATE_BREPS_SBM_UTILITIES_INCLUDED )
#define  KRATOS_CREATE_BREPS_SBM_UTILITIES_INCLUDED


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
#include "geometries/nurbs_volume_geometry.h"

#include "geometries/brep_volume.h"
#include "geometries/brep_surface.h"
#include "geometries/brep_curve.h"

#include "geometries/nurbs_surface_on_volume_geometry.h"
#include "geometries/nurbs_curve_on_surface_geometry.h"

#include "geometries/brep_surface_on_volume.h"
#include "geometries/brep_curve_on_surface.h"


namespace Kratos
{

///@name Kratos Classes
///@{
/// Input for CAD-files.
/** Gives IO capabilities for Nurbs based Brep models in the JSON format defined in
https://amses-journal.springeropen.com/articles/10.1186/s40323-018-0109-4. */
template<class TNodeType = Node, class TEmbeddedNodeType = Point>
class CreateBrepsSBMUtilities : public IO
{
    public:

    ///@}
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CreateBrepsSBMUtilities
    KRATOS_CLASS_POINTER_DEFINITION(CreateBrepsSBMUtilities);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    typedef Geometry<TNodeType> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointerType;

    typedef PointerVector<TNodeType> ContainerNodeType;
    typedef PointerVector<TEmbeddedNodeType> ContainerEmbeddedNodeType;

    typedef CouplingGeometry<TNodeType> CouplingGeometryType;

    typedef NurbsSurfaceGeometry<3, ContainerNodeType> NurbsSurfaceType;
    typedef NurbsCurveGeometry<2, ContainerEmbeddedNodeType> NurbsTrimmingCurveType;

    typedef typename NurbsSurfaceType::Pointer NurbsSurfacePointerType;
    typedef typename NurbsTrimmingCurveType::Pointer NurbsTrimmingCurvePointerType;

    // MODIFIED--------------------------------------------------------------------
    typedef NurbsSurfaceGeometry<3, PointerVector<NodeType>> NurbsSurfaceGeometryType;
    typedef typename NurbsSurfaceGeometryType::Pointer NurbsSurfaceGeometryPointerType;

    typedef NurbsVolumeGeometry<PointerVector<NodeType>> NurbsVolumeGeometryType;
    typedef typename NurbsVolumeGeometryType::Pointer NurbsVolumeGeometryPointerType;
    //-----------------------------------------------------------------------------

    typedef BrepVolume<ContainerNodeType, ContainerEmbeddedNodeType> BrepVolumeType;
    typedef BrepSurfaceOnVolume<ContainerNodeType, ContainerEmbeddedNodeType> BrepSurfaceOnVolumeType;

    typedef BrepSurface<ContainerNodeType, ContainerEmbeddedNodeType> BrepSurfaceType;
    typedef BrepCurveOnSurface<ContainerNodeType, ContainerEmbeddedNodeType> BrepCurveOnSurfaceType;
    
    typedef BrepCurve<ContainerNodeType, ContainerEmbeddedNodeType> BrepCurveType;
    typedef PointOnGeometry<ContainerNodeType, 3, 2> PointOnGeometryOnSurfaceType;
    typedef PointOnGeometry<ContainerNodeType, 3, 1> PointOnGeometryOnCurveType;

    typedef DenseVector<typename BrepCurveOnSurfaceType::Pointer> BrepCurveOnSurfaceArrayType;
    typedef DenseVector<typename BrepCurveOnSurfaceType::Pointer> BrepCurveOnSurfaceLoopType;
    typedef DenseVector<DenseVector<typename BrepCurveOnSurfaceType::Pointer>> BrepCurveOnSurfaceLoopArrayType;

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

    /// Adds all CAD geometries to the herin provided model_part.
    void CreateSurrogateBoundary(NurbsSurfaceGeometryPointerType& p_surface, ModelPart& rModelPart, ModelPart& rSurrogateModelPart_inner, ModelPart& rSurrogateModelPart_outer, const Point& A_uvw, const Point& B_uvw)
    {
        CreateBrepSurface(p_surface, rModelPart, rSurrogateModelPart_inner, rSurrogateModelPart_outer, mEchoLevel);

        CreateBrepCurveOnSurfaces(p_surface, rModelPart, rSurrogateModelPart_inner, rSurrogateModelPart_outer, A_uvw, B_uvw, mEchoLevel);

    }
    
    // 3D !!!!!
    /// Adds all CAD geometries to the herin provided model_part.
    void CreateSurrogateBoundary(NurbsVolumeGeometryPointerType& p_volume, ModelPart& rModelPart, ModelPart& rSurrogateModelPart_inner, ModelPart& rSurrogateModelPart_outer, const Point& A_uvw, const Point& B_uvw)
    {
        // Create brep volume
        CreateBrepVolume(p_volume, rModelPart, rSurrogateModelPart_inner, rSurrogateModelPart_outer, mEchoLevel);

        CreateBrepSurfaceOnVolume(p_volume, rModelPart, rSurrogateModelPart_inner, rSurrogateModelPart_outer, A_uvw, B_uvw, mEchoLevel);
    }


private:

    // 2D
    static void CreateBrepSurface(
        NurbsSurfaceGeometryPointerType p_surface,
        ModelPart& rModelPart,
        ModelPart& rSurrogateModelPart_inner, 
        ModelPart& rSurrogateModelPart_outer,
        SizeType EchoLevel = 0)
    {
        KRATOS_INFO_IF("ReadBrepSurface", (EchoLevel > 3))
            << "Creating BrepSurface \""<< std::endl;

        bool case1 = rSurrogateModelPart_outer.Nodes().size() > 0;
        bool case2 = rSurrogateModelPart_inner.Nodes().size() > 0;
        BrepCurveOnSurfaceLoopArrayType outer_loops, inner_loops;

        if (case1 || case2) p_surface->SetValue(IS_SBM, true);

        // WITHOUT CLIPPER
        auto p_brep_surface =
            Kratos::make_shared<BrepSurfaceType>(
                p_surface, 
                outer_loops,
                inner_loops,
                rSurrogateModelPart_inner,
                rSurrogateModelPart_outer);

        p_surface->SetValue(IS_SBM, true);
        
        /// Sets the brep as geometry parent of the nurbs surface.
        p_surface->SetGeometryParent(p_brep_surface.get());

        SizeType last_geometry_id = rModelPart.GetParentModelPart().Geometries().size();
        p_brep_surface->SetId(1);
        rModelPart.AddGeometry(p_brep_surface);

    }

    // 3D
    static void CreateBrepVolume(
        NurbsVolumeGeometryPointerType& p_volume,
        ModelPart& rModelPart,
        ModelPart& rSurrogateModelPart_inner, 
        ModelPart& rSurrogateModelPart_outer,
        SizeType EchoLevel = 0)
    {
        KRATOS_INFO_IF("ReadBrepVolume", (EchoLevel > 3))
            << "Creating BrepVolume "<< std::endl;

        // WITHOUT CLIPPER
        auto p_brep_volume =
            Kratos::make_shared<BrepVolumeType>(
                p_volume, 
                rSurrogateModelPart_inner,
                rSurrogateModelPart_outer);

        KRATOS_INFO_IF("p_brep_volume has been created", (EchoLevel > 3))<< std::endl;

        /// Sets the brep as geometry parent of the nurbs surface.
        p_volume->SetGeometryParent(p_brep_volume.get());

        SizeType last_geometry_id = 0;
        p_brep_volume->SetId(1);

        rModelPart.AddGeometry(p_brep_volume);
        KRATOS_INFO_IF("p_brep_volume -> AddGeometry", (EchoLevel > 3))<< std::endl;
    }


    static void CreateBrepCurveOnSurfaces(
        NurbsSurfaceGeometryPointerType p_surface,
        ModelPart& rModelPart,
        ModelPart& rSurrogateModelPart_inner, 
        ModelPart& rSurrogateModelPart_outer,
        const Point& A_uvw, const Point& B_uvw,
        SizeType EchoLevel = 0) {
    
        //_____________________________________OUTER
        // Each element in the surrogate_model_part represents a surrogate boundary loop. First "node" is the initial ID of the first surrogate node and
        // the second "node" is the last surrogate node of that loop. (We have done this in the case we have multiple surrogate boundaries and 1 model part)

        int Id_brep_curve_on_surface = 2;

        if (rSurrogateModelPart_outer.Nodes().size() > 0) {
            int sizeSurrogateLoop_outer = rSurrogateModelPart_outer.Nodes().size();
            std::vector<double> surrogatecoord_x_outer(sizeSurrogateLoop_outer);
            std::vector<double> surrogatecoord_y_outer(sizeSurrogateLoop_outer);
            int countSurrogateLoop_outer = 0;
            for (auto i_node = rSurrogateModelPart_outer.NodesEnd()-1; i_node != rSurrogateModelPart_outer.NodesBegin()-1; i_node--) {
                surrogatecoord_x_outer[countSurrogateLoop_outer] = i_node->X();
                surrogatecoord_y_outer[countSurrogateLoop_outer] = i_node->Y();
                countSurrogateLoop_outer++;
            }
            // //********************************************** *OUTER 
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

                // Metti sempre in ordine crescente
                if (active_range_vector[0] > active_range_vector[1]) {
                    double temp = active_range_vector[1];
                    active_range_vector[1] = active_range_vector[0] ;
                    active_range_vector[0] = temp ;
                }

                NurbsInterval brep_active_range(active_range_vector[0], active_range_vector[1]);

                auto p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(
                    p_surface, trimming_curves_GPT_outer[i], brep_active_range, curve_direction);
                p_brep_curve_on_surface->SetId(Id_brep_curve_on_surface);

                rModelPart.AddGeometry(p_brep_curve_on_surface);
                Id_brep_curve_on_surface++;

            }

        } else {
            CreateBrepCurvesOnRectangle(rModelPart, p_surface, A_uvw, B_uvw, Id_brep_curve_on_surface);
        }

        //********************************************** INNER
        // int sizeSurrogateLoop = rSurrogateModelPart_inner.Nodes().size();
        for (int iel = 1; iel < rSurrogateModelPart_inner.Elements().size()+1; iel++) {
            int firstSurrogateNodeId = rSurrogateModelPart_inner.pGetElement(iel)->GetGeometry()[0].Id(); // Element 1 because is the only surrogate loop
            int lastSurrogateNodeId = rSurrogateModelPart_inner.pGetElement(iel)->GetGeometry()[1].Id();  // Element 1 because is the only surrogate loop
            int sizeSurrogateLoop = lastSurrogateNodeId - firstSurrogateNodeId + 1;

            int countSurrogateLoop = 0;
            std::vector<double> surrogatecoord_x(sizeSurrogateLoop);
            std::vector<double> surrogatecoord_y(sizeSurrogateLoop);
            for (int id_node = firstSurrogateNodeId; id_node < lastSurrogateNodeId+1; id_node++) {
                surrogatecoord_x[countSurrogateLoop] = rSurrogateModelPart_inner.GetNode(id_node).X();
                surrogatecoord_y[countSurrogateLoop] = rSurrogateModelPart_inner.GetNode(id_node).Y();
                countSurrogateLoop++;
            }
            //**********************************************
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
                //// Order the active_range_knot_vector
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

                // Metti sempre in ordine crescente
                if (active_range_vector[0] > active_range_vector[1]) {
                    double temp = active_range_vector[1];
                    active_range_vector[1] = active_range_vector[0] ;
                    active_range_vector[0] = temp ;
                }

                NurbsInterval brep_active_range(active_range_vector[0], active_range_vector[1]);

                auto p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(
                    p_surface, trimming_curves_GPT[i], brep_active_range, curve_direction);
                p_brep_curve_on_surface->SetId(Id_brep_curve_on_surface);
                rModelPart.AddGeometry(p_brep_curve_on_surface);
                Id_brep_curve_on_surface++;
                trimming_brep_curve_vector[i] = p_brep_curve_on_surface ;

            }

        }
    } 
    

    static void CreateBrepSurfaceOnVolume(
        NurbsVolumeGeometryPointerType &p_volume,
        ModelPart& rModelPart,
        ModelPart& rSurrogateModelPart_inner, 
        ModelPart& rSurrogateModelPart_outer,
        const Point& A_uvw, const Point& B_uvw,
        SizeType EchoLevel = 0) {
    
        //_____________________________________OUTER
        // Each element in the surrogate_model_part represents a surrogate boundary loop. First "node" is the initial ID of the first surrogate node and
        // the second "node" is the last surrogate node of that loop. (We have done this in the case we have multiple surrogate boundaries and 1 model part)

        int Id_brep_surface_on_volume = 2;

        if (rSurrogateModelPart_outer.Nodes().size() > 0) {

            Vector normal = ZeroVector(3);
            for (auto &i_cond : rSurrogateModelPart_outer.Conditions()) {
                normal = ZeroVector(3);
                Geometry<NodeType>::PointsArrayType points;

                Point A_uvw_sbm = i_cond.GetGeometry()[0];
                Point B_uvw_sbm = i_cond.GetGeometry()[2];
                Point C_uvw_sbm = i_cond.GetGeometry()[1];
                Point D_uvw_sbm = i_cond.GetGeometry()[3];

                Vector knot_vector_u = ZeroVector(2); knot_vector_u[0] = 0.0;
                Vector knot_vector_v = ZeroVector(2); knot_vector_v[0] = 0.0;

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
                    KRATOS_ERROR << "Surrofate face not parallel to any unit direction" << std::endl;
                }
                
                const SizeType p = 1;
                const SizeType q = 1;

                points.push_back(NodeType::Pointer(new NodeType(1, A_uvw_sbm[0], A_uvw_sbm[1], A_uvw_sbm[2])));
                points.push_back(NodeType::Pointer(new NodeType(2, D_uvw_sbm[0], D_uvw_sbm[1], D_uvw_sbm[2])));
                points.push_back(NodeType::Pointer(new NodeType(3, C_uvw_sbm[0], C_uvw_sbm[1], C_uvw_sbm[2])));
                points.push_back(NodeType::Pointer(new NodeType(4, B_uvw_sbm[0], B_uvw_sbm[1], B_uvw_sbm[2])));
                
                auto p_surface_sbm = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(points, p, q, knot_vector_u, knot_vector_v);
                auto brep_p_surface_sbm = Kratos::make_shared<BrepSurface<PointerVector<NodeType>, PointerVector<Point>>>(p_surface_sbm);
                auto p_brep_surface_on_volume_sbm = Kratos::make_shared<BrepSurfaceOnVolume< PointerVector<NodeType>, PointerVector<NodeType>>>(p_volume, p_surface_sbm);
                p_brep_surface_on_volume_sbm->SetId(Id_brep_surface_on_volume);
                
                // Set if the surface on volume is entering of exiting. 
                // Note we save BOUNDARY true to the condition of the surrogate model part when in entering.
                bool isExiting = !i_cond.Is(BOUNDARY);
                int isExitingDirection = -1;
                /* isExitingDirection :
                                -1 -> entering surface/external body fitted surface
                                0 -> exiting surface in x direction
                                1 -> exiting surface in y direction
                                2 -> exiting surface in z direction
                                */
                if (!isExiting) {
                    isExitingDirection = perpendicular_direction;
                    normal = normal * (-1);
                }
                p_brep_surface_on_volume_sbm->SetIsExitingDirectionSBM(isExitingDirection);
                p_brep_surface_on_volume_sbm->SetNormalSBM(normal);
                
                rModelPart.AddGeometry(p_brep_surface_on_volume_sbm);

                Id_brep_surface_on_volume++;
                
            }

        } else {
            CreateBrepSurfaceOnParallelepiped(rModelPart, p_volume, A_uvw, B_uvw, Id_brep_surface_on_volume);
        }

        KRATOS_WATCH("numero outer surfaces on Volume")
        KRATOS_WATCH(Id_brep_surface_on_volume-1)

        // //********************************************** INNER
        Vector normal = ZeroVector(3);
        for (auto &i_cond : rSurrogateModelPart_inner.Conditions()) {
            normal = ZeroVector(3);
            Geometry<NodeType>::PointsArrayType points;

            Point A_uvw_sbm = i_cond.GetGeometry()[0];
            Point B_uvw_sbm = i_cond.GetGeometry()[2];
            Point C_uvw_sbm = i_cond.GetGeometry()[1];
            Point D_uvw_sbm = i_cond.GetGeometry()[3];

            Vector knot_vector_u = ZeroVector(2); knot_vector_u[0] = 0.0;
            Vector knot_vector_v = ZeroVector(2); knot_vector_v[0] = 0.0;

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
                KRATOS_ERROR << "Surrofate face not parallel to any unit direction" << std::endl;
            }
            
            const SizeType p = 1;
            const SizeType q = 1;

            points.push_back(NodeType::Pointer(new NodeType(1, A_uvw_sbm[0], A_uvw_sbm[1], A_uvw_sbm[2])));
            points.push_back(NodeType::Pointer(new NodeType(2, D_uvw_sbm[0], D_uvw_sbm[1], D_uvw_sbm[2])));
            points.push_back(NodeType::Pointer(new NodeType(3, C_uvw_sbm[0], C_uvw_sbm[1], C_uvw_sbm[2])));
            points.push_back(NodeType::Pointer(new NodeType(4, B_uvw_sbm[0], B_uvw_sbm[1], B_uvw_sbm[2])));
            
            auto p_surface_sbm = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(points, p, q, knot_vector_u, knot_vector_v);
            auto brep_p_surface_sbm = Kratos::make_shared<BrepSurface<PointerVector<NodeType>, PointerVector<Point>>>(p_surface_sbm);
            auto p_brep_surface_on_volume_sbm = Kratos::make_shared<BrepSurfaceOnVolume< PointerVector<NodeType>, PointerVector<NodeType>>>(p_volume, p_surface_sbm);
            p_brep_surface_on_volume_sbm->SetId(Id_brep_surface_on_volume);
            
            // Set if the surface on volume is entering of exiting. 
            // Note we save BOUNDARY true to the condition of the surrogate model part when in entering.
            bool isExiting = !i_cond.Is(BOUNDARY);
            int isExitingDirection = -1;
            /* isExitingDirection :
                            -1 -> entering surface/external body fitted surface
                             0 -> exiting surface in x direction
                             1 -> exiting surface in y direction
                             2 -> exiting surface in z direction
                            */
            if (isExiting) {
                isExitingDirection = perpendicular_direction;
                normal = normal * (-1);
            }
            p_brep_surface_on_volume_sbm->SetIsExitingDirectionSBM(isExitingDirection);
            p_brep_surface_on_volume_sbm->SetNormalSBM(normal);
            
            rModelPart.AddGeometry(p_brep_surface_on_volume_sbm);

            Id_brep_surface_on_volume++;
            
        }

        KRATOS_WATCH("numero outer+inner surfaces on Volume")
        KRATOS_WATCH(Id_brep_surface_on_volume-1)


        // // int sizeSurrogateLoop = rSurrogateModelPart_inner.Nodes().size();
        // for (int iel = 1; iel < rSurrogateModelPart_inner.Elements().size()+1; iel++) {
        //     int firstSurrogateNodeId = rSurrogateModelPart_inner.pGetElement(iel)->GetGeometry()[0].Id(); // Element 1 because is the only surrogate loop
        //     int lastSurrogateNodeId = rSurrogateModelPart_inner.pGetElement(iel)->GetGeometry()[1].Id();  // Element 1 because is the only surrogate loop
        //     int sizeSurrogateLoop = lastSurrogateNodeId - firstSurrogateNodeId + 1;

        //     int countSurrogateLoop = 0;
        //     std::vector<double> surrogatecoord_x(sizeSurrogateLoop);
        //     std::vector<double> surrogatecoord_y(sizeSurrogateLoop);
        //     for (int id_node = firstSurrogateNodeId; id_node < lastSurrogateNodeId+1; id_node++) {
        //         surrogatecoord_x[countSurrogateLoop] = rSurrogateModelPart_inner.GetNode(id_node).X();
        //         surrogatecoord_y[countSurrogateLoop] = rSurrogateModelPart_inner.GetNode(id_node).Y();
        //         countSurrogateLoop++;
        //     }
        //     //**********************************************
        //     std::vector<NurbsCurveGeometry<2, PointerVector<Point>>::Pointer> trimming_curves_GPT;

        //     for (std::size_t i = 0; i < surrogatecoord_x.size(); ++i) {
        //         Vector active_range_knot_vector = ZeroVector(2);
                
        //         Point::Pointer point1 = Kratos::make_shared<Point>(surrogatecoord_x[i], surrogatecoord_y[i], 0.0);
        //         Point::Pointer point2 = Kratos::make_shared<Point>(
        //             surrogatecoord_x[(i + 1) % surrogatecoord_x.size()],  // Wrap around for the last point
        //             surrogatecoord_y[(i + 1) % surrogatecoord_y.size()],  // Wrap around for the last point
        //             0.0);
        //         // Compute the knot vector needed
        //         if (surrogatecoord_x[(i + 1) % surrogatecoord_x.size()]==surrogatecoord_x[i]) {
        //             active_range_knot_vector[0] = surrogatecoord_y[i];
        //             active_range_knot_vector[1] = surrogatecoord_y[(i + 1) % surrogatecoord_y.size()];
        //         }
        //         else {
        //             active_range_knot_vector[0] = surrogatecoord_x[i];
        //             active_range_knot_vector[1] = surrogatecoord_x[(i + 1) % surrogatecoord_x.size()];
        //         }
        //         //// Order the active_range_knot_vector
        //         if (active_range_knot_vector[0] > active_range_knot_vector[1]) {
        //             double temp = active_range_knot_vector[1];
        //             active_range_knot_vector[1] = active_range_knot_vector[0] ;
        //             active_range_knot_vector[0] = temp ;
        //         }
        //         NurbsCurveGeometry<2, PointerVector<Point>>::Pointer p_trimming_curve = CreateSingleBrep(point1, point2, active_range_knot_vector);
        //         trimming_curves_GPT.push_back(p_trimming_curve);
        //     }

        //     BrepCurveOnSurfaceLoopType trimming_brep_curve_vector(surrogatecoord_x.size());


        //     for (std::size_t i = 0; i < trimming_curves_GPT.size(); ++i) {
        //         Vector active_range_vector = ZeroVector(2);
        //         if (surrogatecoord_x[(i + 1) % surrogatecoord_x.size()]==surrogatecoord_x[i]) {
        //             active_range_vector[0] = surrogatecoord_y[i];
        //             active_range_vector[1] = surrogatecoord_y[(i + 1) % surrogatecoord_x.size()];
        //         }
        //         else {
        //             active_range_vector[0] = surrogatecoord_x[i];
        //             active_range_vector[1] = surrogatecoord_x[(i + 1) % surrogatecoord_x.size()];
        //         }
        //         bool curve_direction = true;

        //         // Metti sempre in ordine crescente
        //         if (active_range_vector[0] > active_range_vector[1]) {
        //             double temp = active_range_vector[1];
        //             active_range_vector[1] = active_range_vector[0] ;
        //             active_range_vector[0] = temp ;
        //         }

        //         NurbsInterval brep_active_range(active_range_vector[0], active_range_vector[1]);

        //         auto p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(
        //             p_volume, trimming_curves_GPT[i], brep_active_range, curve_direction);
        //         p_brep_curve_on_surface->SetId(Id_brep_curve_on_surface);
        //         rModelPart.AddGeometry(p_brep_curve_on_surface);
        //         Id_brep_curve_on_surface++;
        //         trimming_brep_curve_vector[i] = p_brep_curve_on_surface ;

        //     }

        // }
    }


    ///@}
    ///@name Utility functions
    ///@{


    Parameters mCadJsonParameters;
    int mEchoLevel;

    //// MODIFIED
    static typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer CreateSingleBrep(
        Point::Pointer point1, Point::Pointer point2, Vector active_range_knot_vector)
        {
            // Create the data for the trimming curves
            PointerVector<Point> control_points;
            control_points.push_back(point1);
            control_points.push_back(point2);
            int polynomial_degree = 1; // probably correct
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
        knot_vector[1] = 2.0;
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
        
        
        
        auto brep_curve_on_surface = BrepCurveOnSurface< PointerVector<NodeType>, PointerVector<Point>>(p_surface_geometry, p_curve_1);
        auto p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(p_surface_geometry, p_curve_1);
        p_brep_curve_on_surface->SetId(last_geometry_id);
        r_model_part.AddGeometry(p_brep_curve_on_surface);
        
        brep_curve_on_surface = BrepCurveOnSurface< PointerVector<NodeType>, PointerVector<Point>>(p_surface_geometry, p_curve_2);
        p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(p_surface_geometry, p_curve_2);
        p_brep_curve_on_surface->SetId(++last_geometry_id);
        r_model_part.AddGeometry(p_brep_curve_on_surface);

        brep_curve_on_surface = BrepCurveOnSurface< PointerVector<NodeType>, PointerVector<Point>>(p_surface_geometry, p_curve_3);
        p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(p_surface_geometry, p_curve_3);
        p_brep_curve_on_surface->SetId(++last_geometry_id);
        r_model_part.AddGeometry(p_brep_curve_on_surface);
        
        brep_curve_on_surface = BrepCurveOnSurface< PointerVector<NodeType>, PointerVector<Point>>(p_surface_geometry, p_curve_4);
        p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(p_surface_geometry, p_curve_4);
        p_brep_curve_on_surface->SetId(++last_geometry_id);
        r_model_part.AddGeometry(p_brep_curve_on_surface);

        last_geometry_id++;
    }
    
    // 3D
    static void CreateBrepSurfaceOnParallelepiped(ModelPart& r_model_part, NurbsVolumeGeometryPointerType p_volume_geometry, const Point& A_uvw, const Point& B_uvw, int &last_geometry_id) {
        // Vector knot_vector_u = ZeroVector(2); knot_vector_u[0] = 0.0; knot_vector_u[1] = 1.0;  
        // Vector knot_vector_v = ZeroVector(2); knot_vector_v[0] = 0.0; knot_vector_v[1] = 1.0; 
        
        Vector knot_vector_u = ZeroVector(2); 
        Vector knot_vector_v = ZeroVector(2); 
        Vector Normal = ZeroVector(3);
        const SizeType p = 1;
        const SizeType q = 1;
        
        

        // TO DO: DOES NOT WORK IF THE KNOT SPANS ARE RECTANGLES
        // lower
        Geometry<NodeType>::PointsArrayType points_lower;
        knot_vector_u[0] = 0.0; knot_vector_u[1] = std::abs(A_uvw[0]-B_uvw[0]);  
        knot_vector_v[0] = 0.0; knot_vector_v[1] = std::abs(A_uvw[1]-B_uvw[1]); 
        points_lower.push_back(NodeType::Pointer(new NodeType(1, A_uvw[0], A_uvw[1], A_uvw[2])));
        points_lower.push_back(NodeType::Pointer(new NodeType(2, A_uvw[0], B_uvw[1], A_uvw[2])));// CAMBIATO
        points_lower.push_back(NodeType::Pointer(new NodeType(3, B_uvw[0], A_uvw[1], A_uvw[2])));// CAMBIATO
        points_lower.push_back(NodeType::Pointer(new NodeType(4, B_uvw[0], B_uvw[1], A_uvw[2])));
        
        Normal = ZeroVector(3);
        Normal[2] = -1;
        auto p_surface_1 = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(points_lower, p, q, knot_vector_u, knot_vector_v);
        auto brep_p_surface_1 = Kratos::make_shared<BrepSurface<PointerVector<NodeType>, PointerVector<Point>>>(p_surface_1);
        auto p_brep_surface_on_volume_1 = Kratos::make_shared<BrepSurfaceOnVolume< PointerVector<NodeType>, PointerVector<NodeType>>>(p_volume_geometry, p_surface_1);
        // auto p_brep_surface_on_volume = Kratos::make_shared<BrepSurfaceOnVolumeType>(p_volume_geometry, brep_p_surface_1);
        p_brep_surface_on_volume_1->SetId(last_geometry_id);
        p_brep_surface_on_volume_1->SetNormalSBM(Normal);
        r_model_part.AddGeometry(p_brep_surface_on_volume_1);

        // upper
        Geometry<NodeType>::PointsArrayType points_upper;
        knot_vector_u[0] = 0.0; knot_vector_u[1] = std::abs(A_uvw[0]-B_uvw[0]);  
        knot_vector_v[0] = 0.0; knot_vector_v[1] = std::abs(A_uvw[1]-B_uvw[1]); 

        points_upper.push_back(NodeType::Pointer(new NodeType(1, A_uvw[0], A_uvw[1], B_uvw[2])));
        points_upper.push_back(NodeType::Pointer(new NodeType(2, A_uvw[0], B_uvw[1], B_uvw[2]))); // CAMBIATO
        points_upper.push_back(NodeType::Pointer(new NodeType(3, B_uvw[0], A_uvw[1], B_uvw[2]))); // CAMBIATO
        points_upper.push_back(NodeType::Pointer(new NodeType(4, B_uvw[0], B_uvw[1], B_uvw[2])));
        Normal = ZeroVector(3);
        Normal[2] = 1;
        auto p_surface_2 = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(points_upper, p, q, knot_vector_u, knot_vector_v);
        auto brep_p_surface_2 = Kratos::make_shared<BrepSurface<PointerVector<NodeType>, PointerVector<Point>>>(p_surface_2);
        auto p_brep_surface_on_volume_2 = Kratos::make_shared<BrepSurfaceOnVolume< PointerVector<NodeType>, PointerVector<NodeType>>>(p_volume_geometry, p_surface_2);
        p_brep_surface_on_volume_2->SetId(++last_geometry_id);
        p_brep_surface_on_volume_2->SetNormalSBM(Normal);
        r_model_part.AddGeometry(p_brep_surface_on_volume_2);

        // front
        Geometry<NodeType>::PointsArrayType points_front;
        knot_vector_u[0] = 0.0; knot_vector_u[1] = std::abs(A_uvw[0]-B_uvw[0]);  
        knot_vector_v[0] = 0.0; knot_vector_v[1] = std::abs(A_uvw[2]-B_uvw[2]); 
        
        points_front.push_back(NodeType::Pointer(new NodeType(1, A_uvw[0], A_uvw[1], A_uvw[2])));
        points_front.push_back(NodeType::Pointer(new NodeType(2, B_uvw[0], A_uvw[1], A_uvw[2])));
        points_front.push_back(NodeType::Pointer(new NodeType(3, A_uvw[0], A_uvw[1], B_uvw[2]))); 
        points_front.push_back(NodeType::Pointer(new NodeType(4, B_uvw[0], A_uvw[1], B_uvw[2])));
        Normal = ZeroVector(3);
        Normal[1] = -1;
        auto p_surface_3 = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(points_front, p, q, knot_vector_u, knot_vector_v);
        auto brep_p_surface_3 = Kratos::make_shared<BrepSurface<PointerVector<NodeType>, PointerVector<Point>>>(p_surface_3);
        auto p_brep_surface_on_volume_3 = Kratos::make_shared<BrepSurfaceOnVolume< PointerVector<NodeType>, PointerVector<NodeType>>>(p_volume_geometry, p_surface_3);
        p_brep_surface_on_volume_3->SetId(++last_geometry_id);
        p_brep_surface_on_volume_3->SetNormalSBM(Normal);
        r_model_part.AddGeometry(p_brep_surface_on_volume_3);

        // back
        Geometry<NodeType>::PointsArrayType points_back;
        knot_vector_u[0] = 0.0; knot_vector_u[1] = std::abs(A_uvw[0]-B_uvw[0]);  
        knot_vector_v[0] = 0.0; knot_vector_v[1] = std::abs(A_uvw[2]-B_uvw[2]); 
        points_back.push_back(NodeType::Pointer(new NodeType(1, A_uvw[0], B_uvw[1], A_uvw[2])));
        points_back.push_back(NodeType::Pointer(new NodeType(2, A_uvw[0], B_uvw[1], B_uvw[2]))); // CAMBIATO
        points_back.push_back(NodeType::Pointer(new NodeType(3, B_uvw[0], B_uvw[1], A_uvw[2]))); // CAMBIATO
        points_back.push_back(NodeType::Pointer(new NodeType(4, B_uvw[0], B_uvw[1], B_uvw[2])));
        Normal = ZeroVector(3);
        Normal[1] = 1;
        auto p_surface_4 = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(points_back, p, q, knot_vector_u, knot_vector_v);
        auto brep_p_surface_4 = Kratos::make_shared<BrepSurface<PointerVector<NodeType>, PointerVector<Point>>>(p_surface_4);
        auto p_brep_surface_on_volume_4 = Kratos::make_shared<BrepSurfaceOnVolume< PointerVector<NodeType>, PointerVector<NodeType>>>(p_volume_geometry, p_surface_4);
        p_brep_surface_on_volume_4->SetId(++last_geometry_id);
        p_brep_surface_on_volume_4->SetNormalSBM(Normal);
        r_model_part.AddGeometry(p_brep_surface_on_volume_4);

        // left
        Geometry<NodeType>::PointsArrayType points_left;
        knot_vector_u[0] = 0.0; knot_vector_u[1] = std::abs(A_uvw[1]-B_uvw[1]);  
        knot_vector_v[0] = 0.0; knot_vector_v[1] = std::abs(A_uvw[2]-B_uvw[2]); 
        points_left.push_back(NodeType::Pointer(new NodeType(1, A_uvw[0], A_uvw[1], A_uvw[2])));
        points_left.push_back(NodeType::Pointer(new NodeType(2, A_uvw[0], B_uvw[1], A_uvw[2])));
        points_left.push_back(NodeType::Pointer(new NodeType(3, A_uvw[0], A_uvw[1], B_uvw[2])));
        points_left.push_back(NodeType::Pointer(new NodeType(4, A_uvw[0], B_uvw[1], B_uvw[2])));
        Normal = ZeroVector(3);
        Normal[0] = -1;
        auto p_surface_5 = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(points_left, p, q, knot_vector_u, knot_vector_v);
        auto brep_p_surface_5 = Kratos::make_shared<BrepSurface<PointerVector<NodeType>, PointerVector<Point>>>(p_surface_5);
        auto p_brep_surface_on_volume_5 = Kratos::make_shared<BrepSurfaceOnVolume< PointerVector<NodeType>, PointerVector<NodeType>>>(p_volume_geometry, p_surface_5);
        // auto p_brep_surface_on_volume = Kratos::make_shared<BrepSurfaceOnVolumeType>(p_volume_geometry, brep_p_surface_1);
        p_brep_surface_on_volume_5->SetId(++last_geometry_id);
        p_brep_surface_on_volume_5->SetNormalSBM(Normal);
        r_model_part.AddGeometry(p_brep_surface_on_volume_5);

        // right
        Geometry<NodeType>::PointsArrayType points_right;
        knot_vector_u[0] = 0.0; knot_vector_u[1] = std::abs(A_uvw[1]-B_uvw[1]);  
        knot_vector_v[0] = 0.0; knot_vector_v[1] = std::abs(A_uvw[2]-B_uvw[2]); 
        points_right.push_back(NodeType::Pointer(new NodeType(1, B_uvw[0], A_uvw[1], A_uvw[2])));
        points_right.push_back(NodeType::Pointer(new NodeType(2 ,B_uvw[0], A_uvw[1], B_uvw[2]))); // CAMBIATO
        points_right.push_back(NodeType::Pointer(new NodeType(3 ,B_uvw[0], B_uvw[1], A_uvw[2]))); // CAMBIATO
        points_right.push_back(NodeType::Pointer(new NodeType(4 ,B_uvw[0], B_uvw[1], B_uvw[2])));
        Normal = ZeroVector(3);
        Normal[0] = 1;
        auto p_surface_6 = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(points_right, p, q, knot_vector_u, knot_vector_v);
        auto brep_p_surface_6 = Kratos::make_shared<BrepSurface<PointerVector<NodeType>, PointerVector<Point>>>(p_surface_6);
        auto p_brep_surface_on_volume_6 = Kratos::make_shared<BrepSurfaceOnVolume< PointerVector<NodeType>, PointerVector<NodeType>>>(p_volume_geometry, p_surface_6);
        p_brep_surface_on_volume_6->SetId(++last_geometry_id);
        p_brep_surface_on_volume_6->SetNormalSBM(Normal);
        r_model_part.AddGeometry(p_brep_surface_on_volume_6);

        last_geometry_id++;
    }

    ///@}
}; // Class CreateBrepsSBMUtilities
}  // namespace Kratos.

#endif // KRATOS_CAD_JSON_INPUT_INCLUDED  defined
