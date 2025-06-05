//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_CAD_JSON_INPUT_INCLUDED )
#define  KRATOS_CAD_JSON_INPUT_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/io.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

// Geometries
#include "geometries/coupling_geometry.h"

#include "geometries/nurbs_curve_geometry.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/nurbs_curve_on_surface_geometry.h"

#include "geometries/brep_surface.h"
#include "geometries/brep_curve_on_surface.h"
#include "geometries/brep_curve.h"

#include "geometries/point_on_geometry.h"

#include "geometries/thb_surface_geometry.h"
#include "geometries/thb_brep_surface.h"

namespace Kratos
{

///@name Kratos Classes
///@{
/// Input for CAD-files.
/** Gives IO capabilities for Nurbs based Brep models in the JSON format defined in
https://amses-journal.springeropen.com/articles/10.1186/s40323-018-0109-4. */
template<class TNodeType = Node, class TEmbeddedNodeType = Point>
class CadJsonInput : public IO
{
    public:

    ///@}
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CadJsonInput
    KRATOS_CLASS_POINTER_DEFINITION(CadJsonInput);

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

    typedef BrepSurface<ContainerNodeType, false, ContainerEmbeddedNodeType> BrepSurfaceType;
    typedef BrepCurveOnSurface<ContainerNodeType, false, ContainerEmbeddedNodeType> BrepCurveOnSurfaceType;
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
    CadJsonInput(
        const std::string& rDataFileName,
        SizeType EchoLevel = 0)
        : mEchoLevel(EchoLevel)
    {
        mCadJsonParameters = ReadParamatersFile(rDataFileName, EchoLevel);
    }

    /// Constructor with KratosParameters.
    CadJsonInput(
        Parameters CadJsonParameters,
        SizeType EchoLevel = 0)
        : mCadJsonParameters(CadJsonParameters)
        , mEchoLevel(EchoLevel)
    {
    }

    /// Destructor.
    ~CadJsonInput() = default;

    ///@}
    ///@name Python exposed Functions
    ///@{

    /// Adds all CAD geometries to the herin provided model_part.
    void ReadModelPart(ModelPart& rModelPart) override
    {
        ReadGeometryModelPart(mCadJsonParameters, rModelPart, mEchoLevel);
    }

    ///@}

private:
    ///@name Static Functions
    ///@{

    /// Allows static access without own memory.
    static void ReadGeometryModelPart(
        const Parameters rCadJsonParameters,
        ModelPart& rModelPart,
        SizeType EchoLevel = 0)
    {
        KRATOS_ERROR_IF_NOT(rCadJsonParameters.Has("breps"))
            << "Missing \"breps\" section" << std::endl;

        ReadBreps(rCadJsonParameters["breps"], rModelPart, EchoLevel);
    }

    ///@}
    ///@name Read in Brep
    ///@{

    static void ReadBreps(
        const Parameters rParameters,
        ModelPart& rModelPart,
        SizeType EchoLevel = 0)
    {
        for (IndexType brep_index = 0; brep_index < rParameters.size(); brep_index++)
        {
            KRATOS_INFO_IF("ReadBreps", (EchoLevel > 0))
                << "Reading Brep \"" << GetIdOrName(rParameters[brep_index])
                << "\" - faces." << std::endl;

            if (rParameters[brep_index].Has("faces"))
            {
                ReadBrepSurfaces(rParameters[brep_index]["faces"], rModelPart, EchoLevel);
            }
        }

        for (IndexType brep_index = 0; brep_index < rParameters.size(); brep_index++)
        {
            KRATOS_INFO_IF("ReadBreps", (EchoLevel > 0))
                << "Reading Brep \"" << GetIdOrName(rParameters[brep_index])
                << "\" - edges." << std::endl;

            if (rParameters[brep_index].Has("edges"))
            {
                ReadBrepCurveOnSurfaces(rParameters[brep_index]["edges"], rModelPart, EchoLevel);
            }
        }

        for (IndexType brep_index = 0; brep_index < rParameters.size(); brep_index++)
        {
            KRATOS_INFO_IF("ReadBreps", (EchoLevel > 0))
                << "Reading Brep \"" << GetIdOrName(rParameters[brep_index])
                << "\" - points." << std::endl;

            if (rParameters[brep_index].Has("vertices"))
            {
                ReadPointsOnGeometries(rParameters[brep_index]["vertices"], rModelPart, EchoLevel);
            }
        }
    }

    ///@}
    ///@name Read in Brep Geometries
    ///@{

    static void ReadBrepSurfaces(
        const Parameters rParameters,
        ModelPart& rModelPart,
        SizeType EchoLevel = 0)
    {
        KRATOS_ERROR_IF_NOT(rParameters.IsArray())
            << "\"faces\" section needs to be an array of BrepSurfaces." << std::endl;

        KRATOS_INFO_IF("ReadBrepSurfaces", EchoLevel > 2)
            << "Reading " << rParameters.size() << " BrepSurfaces..." << std::endl;

        for (IndexType brep_surface_i = 0; brep_surface_i < rParameters.size(); ++brep_surface_i)
        {
            ReadBrepSurface(rParameters[brep_surface_i], rModelPart, EchoLevel);
        }
    }

    static void ReadBrepSurface(
        const Parameters rParameters,
        ModelPart& rModelPart,
        SizeType EchoLevel = 0)
    {
        KRATOS_INFO_IF("ReadBrepSurface", (EchoLevel > 3))
            << "Reading BrepSurface \"" << GetIdOrName(rParameters) << "\"" << std::endl;

        KRATOS_ERROR_IF_NOT(HasIdOrName(rParameters))
            << "Missing 'brep_id' or 'brep_name' in brep face." << std::endl;

        KRATOS_ERROR_IF_NOT(rParameters.Has("surface"))
            << "Missing 'surface' in brep face." << std::endl;

        auto p_surface(ReadNurbsSurface<3, TNodeType>(
            rParameters["surface"], rModelPart, EchoLevel));

        const bool is_trimmed = (rParameters["surface"].Has("is_trimmed"))
            ? rParameters["surface"]["is_trimmed"].GetBool()
            : true;
        KRATOS_INFO_IF("ReadBrepSurface", (EchoLevel > 4) && !rParameters["surface"].Has("is_trimmed"))
            << "For BrepSurface \"" << GetIdOrName(rParameters) << "\""
            << "\", is_trimmed is not provided in the input."
            << " is_trimmed = true is considered." << std::endl;

        if (rParameters.Has("boundary_loops"))
        {
            BrepCurveOnSurfaceLoopArrayType outer_loops, inner_loops;
            tie(outer_loops, inner_loops) =
                ReadBoundaryLoops(rParameters["boundary_loops"], p_surface, rModelPart, EchoLevel);

            auto p_brep_surface =
                Kratos::make_shared<BrepSurfaceType>(
                    p_surface,
                    outer_loops,
                    inner_loops,
                    is_trimmed);

            /// Sets the brep as geometry parent of the nurbs surface.
            p_surface->SetGeometryParent(p_brep_surface.get());

            SetIdOrName<BrepSurfaceType>(rParameters, p_brep_surface);

            ReadAndAddEmbeddedEdges(p_brep_surface, rParameters, p_surface, rModelPart, EchoLevel);

            rModelPart.AddGeometry(p_brep_surface);
        }
        else
        {
            KRATOS_INFO_IF("ReadBrepSurface", (EchoLevel > 4))
                << "For BrepSurface \"" << GetIdOrName(rParameters) << "\""
                << "\", boundary_loops are not provided in the input."
                << " It will be considered as untrimmed." << std::endl;

            auto p_brep_surface =
                Kratos::make_shared<BrepSurfaceType>(
                    p_surface);

            SetIdOrName<BrepSurfaceType>(rParameters, p_brep_surface);

            ReadAndAddEmbeddedEdges(p_brep_surface, rParameters, p_surface, rModelPart, EchoLevel);

            rModelPart.AddGeometry(p_brep_surface);
        }
    }

    ///@}
    ///@name Read in Surface Trimming
    ///@{

    static BrepCurveOnSurfaceLoopType
        ReadTrimmingCurveVector(
            const Parameters rParameters,
            typename NurbsSurfaceType::Pointer pNurbsSurface,
            ModelPart& rModelPart,
            SizeType EchoLevel = 0)
    {
        KRATOS_ERROR_IF(rParameters.size() < 1)
            << "Trimming curve list has no element." << std::endl;

        BrepCurveOnSurfaceLoopType
            trimming_brep_curve_vector(rParameters.size());

        for (IndexType tc_idx = 0; tc_idx < rParameters.size(); tc_idx++)
        {
            trimming_brep_curve_vector[tc_idx] = ReadTrimmingCurve(
                rParameters[tc_idx], pNurbsSurface, rModelPart, EchoLevel);
        }

        return trimming_brep_curve_vector;
    }

    static typename BrepCurveOnSurfaceType::Pointer
        ReadTrimmingCurve(
            const Parameters rParameters,
            typename NurbsSurfaceType::Pointer pNurbsSurface,
            ModelPart& rModelPart,
            SizeType EchoLevel = 0)
    {
        KRATOS_ERROR_IF_NOT(rParameters.Has("curve_direction"))
            << "Missing 'curve_direction' in nurbs curve" << std::endl;
        bool curve_direction = rParameters["curve_direction"].GetBool();

        KRATOS_ERROR_IF_NOT(rParameters.Has("parameter_curve"))
            << "Missing 'parameter_curve' in nurbs curve" << std::endl;

        auto p_trimming_curve(ReadNurbsCurve<2, TEmbeddedNodeType>(
            rParameters["parameter_curve"], rModelPart, EchoLevel));

        KRATOS_ERROR_IF_NOT(rParameters["parameter_curve"].Has("active_range"))
            << "Missing 'active_range' in parameter_curve, in trimming curve." << std::endl;
        Vector active_range_vector = rParameters["parameter_curve"]["active_range"].GetVector();
        NurbsInterval brep_active_range(active_range_vector[0], active_range_vector[1]);

        auto p_brep_curve_on_surface
            = Kratos::make_shared<BrepCurveOnSurfaceType>(
                pNurbsSurface, p_trimming_curve, brep_active_range, curve_direction);

        if (rParameters.Has("trim_index")) {
            p_brep_curve_on_surface->SetId(rParameters["trim_index"].GetInt());
        }

        return p_brep_curve_on_surface;
    }

    static std::tuple<BrepCurveOnSurfaceLoopArrayType, BrepCurveOnSurfaceLoopArrayType>
        ReadBoundaryLoops(
            const Parameters rParameters,
            typename NurbsSurfaceType::Pointer pNurbsSurface,
            ModelPart& rModelPart,
            SizeType EchoLevel = 0)
    {
        BrepCurveOnSurfaceLoopArrayType outer_loops;
        BrepCurveOnSurfaceLoopArrayType inner_loops;

        for (IndexType bl_idx = 0; bl_idx < rParameters.size(); bl_idx++)
        {
            KRATOS_ERROR_IF_NOT(rParameters[bl_idx].Has("loop_type"))
                << "Missing 'loop_type' in boundary loops, in "
                << bl_idx << " loop." << std::endl;
            std::string loop_type = rParameters[bl_idx]["loop_type"].GetString();

            KRATOS_ERROR_IF_NOT(rParameters[bl_idx].Has("trimming_curves"))
                << "Missing 'trimming_curves' in boundary loops"
                << bl_idx << " loop." << std::endl;
            auto trimming_curves(ReadTrimmingCurveVector(
                rParameters[bl_idx]["trimming_curves"], pNurbsSurface, rModelPart, EchoLevel));

            if (loop_type == "outer")
            {
                outer_loops.resize(outer_loops.size() + 1);
                outer_loops[outer_loops.size() - 1] = trimming_curves;
            }
            else if (loop_type == "inner")
            {
                inner_loops.resize(inner_loops.size() + 1);
                inner_loops[inner_loops.size() - 1] = trimming_curves;
            }
            else
            {
                KRATOS_ERROR << "Loop type: " << loop_type
                    << " is not supported." << std::endl;
            }
        }

        return std::make_tuple(outer_loops, inner_loops);
    }

    static void ReadAndAddEmbeddedEdges(
            typename BrepSurfaceType::Pointer pBrepSurface,
            const Parameters rParameters,
            typename NurbsSurfaceType::Pointer pNurbsSurface,
            ModelPart& rModelPart,
            SizeType EchoLevel = 0)
    {
        if (rParameters.Has("embedded_edges")) {
            if (rParameters["embedded_edges"].size() > 0) {
                BrepCurveOnSurfaceArrayType embedded_edges(ReadTrimmingCurveVector(
                    rParameters["embedded_edges"], pNurbsSurface, rModelPart, EchoLevel));

                pBrepSurface->AddEmbeddedEdges(embedded_edges);
            }
        }
    }

    ///@}
    ///@name Read in Nurbs Geometries
    ///@{

    static void ReadBrepCurveOnSurfaces(
        const Parameters rParameters,
        ModelPart& rModelPart,
        SizeType EchoLevel = 0)
    {
        KRATOS_ERROR_IF_NOT(rParameters.IsArray())
            << "\"faces\" section needs to be an array of BrepSurfaces." << std::endl;

        KRATOS_INFO_IF("ReadBrepCurveOnSurfaces", EchoLevel > 2)
            << "Reading " << rParameters.size() << " BrepEdge..." << std::endl;

        for (IndexType i = 0; i < rParameters.size(); i++)
        {
            ReadBrepEdge(rParameters[i], rModelPart, EchoLevel);
        }
    }

    static void ReadBrepEdge(
        const Parameters rParameters,
        ModelPart& rModelPart,
        SizeType EchoLevel = 0)
    {
        KRATOS_ERROR_IF_NOT(HasIdOrName(rParameters))
            << "Missing 'brep_id' or 'brep_name' in brep edge" << std::endl;

        if (rParameters.Has("topology"))
        {
            if (rParameters["topology"].size() == 0)
            {
                ReadBrepCurve(rParameters, rModelPart, EchoLevel);
            }
            else if (rParameters["topology"].size() == 1)
            {
                ReadBrepEdgeBrepCurveOnSurface(rParameters, rModelPart, EchoLevel);
            }
            else { // More than one topology means that a coupling geometry is required.
                ReadCouplingGeometry(rParameters, rModelPart, EchoLevel);
            }
        }
    }

    static void ReadBrepCurve(
        const Parameters rParameters,
        ModelPart& rModelPart,
        SizeType EchoLevel = 0)
    {
        KRATOS_ERROR_IF_NOT(HasIdOrName(rParameters))
            << "Missing 'brep_id' or 'brep_name' in brep curve." << std::endl;

        KRATOS_INFO_IF("ReadBrepCurve", (EchoLevel > 3))
            << "Reading BrepCurve \"" << GetIdOrName(rParameters) << "\"" << std::endl;

        KRATOS_ERROR_IF_NOT(rParameters.Has("3d_curve"))
            << "Missing '3d_curve' in brep curve." << std::endl;

        auto p_curve = ReadNurbsCurve<3, TNodeType>(
            rParameters["3d_curve"], rModelPart, EchoLevel);

        auto p_brep_curve = Kratos::make_shared<BrepCurveType>(p_curve);

        SetIdOrName<BrepCurveType>(rParameters, p_brep_curve);

        rModelPart.AddGeometry(p_brep_curve);
    }


    static void ReadBrepEdgeBrepCurveOnSurface(
        const Parameters & rParameters,
        ModelPart & rModelPart,
        SizeType EchoLevel = 0)
    {
        KRATOS_INFO_IF("ReadBrepEdge", (EchoLevel > 3))
            << "Reading BrepEdge \"" << GetIdOrName(rParameters) << "\"" << std::endl;

        KRATOS_ERROR_IF_NOT(HasIdOrName(rParameters["topology"][0]))
            << "Missing 'brep_id' or 'brep_name' in topology" << std::endl;

        KRATOS_INFO_IF("ReadBrepEdge", (EchoLevel > 4))
            << "Getting trim: \"" << rParameters["topology"][0]["trim_index"].GetInt()
            << "\" from geometry: \"" << GetIdOrName(rParameters["topology"][0])
            << "\"." << std::endl;

        GeometryPointerType p_geometry = GetGeometry(rParameters["topology"][0], rModelPart);
        GeometryPointerType p_brep_trim =
            p_geometry->pGetGeometryPart(rParameters["topology"][0]["trim_index"].GetInt());

        auto p_brep_curve_on_surface
            = dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_brep_trim);
        KRATOS_ERROR_IF(p_brep_curve_on_surface == nullptr)
            << "dynamic_cast from Geometry to BrepCurveOnSurface not successfull. Brep Id: "
            << GetIdOrName(rParameters["topology"][0]) << " and trim index: "
            << rParameters["topology"][0]["trim_index"].GetInt() << std::endl;

        bool relative_direction = true;
        if (rParameters["topology"][0].Has("relative_direction")) {
            relative_direction = rParameters["topology"][0]["relative_direction"].GetBool();
        }
        else {
            KRATOS_INFO_IF("ReadBrepEdge", (EchoLevel > 4))
                << "For trim: \"" << rParameters["topology"][0]["trim_index"].GetInt()
                << "\" from geometry: \"" << GetIdOrName(rParameters["topology"][0])
                << "\", no relative_direction is provided in the input." << std::endl;
        }

        auto p_nurbs_curve_on_surface = p_brep_curve_on_surface->pGetCurveOnSurface();

        auto brep_nurbs_interval = p_brep_curve_on_surface->DomainInterval();
        auto p_bre_edge_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(
            p_nurbs_curve_on_surface, brep_nurbs_interval, relative_direction);

        SetIdOrName<BrepCurveOnSurfaceType>(rParameters, p_bre_edge_brep_curve_on_surface);

        rModelPart.AddGeometry(p_bre_edge_brep_curve_on_surface);
    }

    static void ReadCouplingGeometry(
        const Parameters rParameters,
        ModelPart& rModelPart,
        SizeType EchoLevel = 0)
    {
        KRATOS_INFO_IF("ReadCouplingGeometry", (EchoLevel > 3))
            << "Reading CouplingGeometry \"" << GetIdOrName(rParameters) << "\"" << std::endl;

        typename CouplingGeometryType::GeometryPointerVector geometry_vector;

        for (IndexType i = 0; i < rParameters["topology"].size(); i++)
        {
            KRATOS_ERROR_IF_NOT(HasIdOrName(rParameters["topology"][0]))
                << "Missing 'brep_id' or 'brep_name' in topology of Coupling - BrepEdge" << std::endl;

            auto p_geometry = GetGeometry(rParameters["topology"][i], rModelPart);

            geometry_vector.push_back(
                p_geometry->pGetGeometryPart(rParameters["topology"][i]["trim_index"].GetInt()));
        }

        auto p_coupling_geometry = Kratos::make_shared<CouplingGeometryType>(
            geometry_vector);

        SetIdOrName<CouplingGeometryType>(rParameters, p_coupling_geometry);

        rModelPart.AddGeometry(p_coupling_geometry);
    }

    static void ReadPointsOnGeometries(
        const Parameters rParameters,
        ModelPart& rModelPart,
        SizeType EchoLevel = 0)
    {
        KRATOS_ERROR_IF_NOT(rParameters.IsArray())
            << "\"vertices\" section needs to be an array of PointsOnGeometries." << std::endl;

        KRATOS_INFO_IF("ReadPointsOnGeometries", EchoLevel > 2)
            << "Reading " << rParameters.size() << " PointsOnGeometries..." << std::endl;

        for (IndexType brep_point_i = 0; brep_point_i < rParameters.size(); ++brep_point_i)
        {
            KRATOS_INFO_IF("ReadPointsOnGeometries", (EchoLevel > 3))
                << "Reading PointOnGeometry \"" << GetIdOrName(rParameters[brep_point_i]) << "\"" << std::endl;

            if (rParameters[brep_point_i]["topology"].size() == 1) {
                if (rParameters[brep_point_i]["topology"][0].Has("local_coordinates"))
                {
                    auto p_geometry = GetGeometry(rParameters[brep_point_i]["topology"][0], rModelPart);

                    GeometryPointerType p_point_on_geometry = ReadPointOnGeometry(rParameters[brep_point_i]["topology"][0], rModelPart, p_geometry, EchoLevel);

                    SetIdOrName(rParameters[brep_point_i], p_point_on_geometry);
                    rModelPart.AddGeometry(p_point_on_geometry);
                }
                else {
                    KRATOS_ERROR << "Topology of brep point: " << GetIdOrName(rParameters[brep_point_i])
                        << " does not provide local coordinates. No other import option provided. Topology as following: "
                        << rParameters[brep_point_i] << std::endl;
                }
            }
            else if (rParameters[brep_point_i]["topology"].size() > 1) {
                std::vector<GeometryPointerType> coupling_point_geometries(rParameters[brep_point_i]["topology"].size());

                for (IndexType i = 0; i < rParameters[brep_point_i]["topology"].size(); ++i) {
                    auto p_geometry = GetGeometry(rParameters[brep_point_i]["topology"][i], rModelPart);

                    coupling_point_geometries[i] = ReadPointOnGeometry(rParameters[brep_point_i]["topology"][i], rModelPart, p_geometry, EchoLevel);
                }

                auto p_coupling_geometry = Kratos::make_shared<CouplingGeometryType>(
                    coupling_point_geometries);

                SetIdOrName<CouplingGeometryType>(rParameters[brep_point_i], p_coupling_geometry);
                rModelPart.AddGeometry(p_coupling_geometry);
            }
            else {
                KRATOS_ERROR << "PointOnGeometry did not provide any topology, which is not allowed.\n Parameters are read as following: \n"
                    << rParameters[brep_point_i] << std::endl;
            }
        }
    }

    static GeometryPointerType ReadPointOnGeometry(
        const Parameters rParameters,
        ModelPart& rModelPart,
        GeometryPointerType pBackgroundGeometry,
        SizeType EchoLevel = 0)
    {
        auto coordinates_vector = rParameters["local_coordinates"].GetVector();
        array_1d<double, 3> local_coordinates;
        local_coordinates[0] = coordinates_vector[0];
        local_coordinates[1] = coordinates_vector[1];
        local_coordinates[2] = coordinates_vector[2];

        if (pBackgroundGeometry->LocalSpaceDimension() == 2)
            return Kratos::make_shared<PointOnGeometryOnSurfaceType>(local_coordinates, pBackgroundGeometry);
        else if (pBackgroundGeometry->LocalSpaceDimension() == 1)
            return Kratos::make_shared<PointOnGeometryOnCurveType>(local_coordinates, pBackgroundGeometry);
        else {
            KRATOS_ERROR << "Local space dimension of: " << pBackgroundGeometry->LocalSpaceDimension()
                << " is not supported for point on geometry." << std::endl;
        }
    }

    ///@}
    ///@name Read in Nurbs Geometries
    ///@{

    /* @brief read NurbsCurves from the given parameter input.
    *
    * The input needs to be provided in the following shape:
    * {
    *     "is_rational": bool,
    *     "degree": p:int,
    *     "knot_vector": [ double, ... ],
    *     "active_range": [ double, double ],
    *     "control_points": [
    *         [ id: int, [ x, y, z, weight ] ],
    *         ...
    *     ]
    * }
    */
    template<int TWorkingSpaceDimension, class TThisNodeType>
    static typename NurbsCurveGeometry<TWorkingSpaceDimension, PointerVector<TThisNodeType>>::Pointer
        ReadNurbsCurve(
            const Parameters rParameters,
            ModelPart& rModelPart,
            SizeType EchoLevel = 0)
    {
        bool is_rational = true;
        if (rParameters.Has("is_rational")) {
            is_rational = rParameters["is_rational"].GetBool();
        }
        else {
            KRATOS_INFO_IF("ReadNurbsCurve", (EchoLevel > 4))
                << "\"is_rational\" is not provided within \"surface\". Thus, it is considered as rational. "
                << "If this curve is non-rational the computation of the shape functions is less optimized."
                << std::endl;
        }

        KRATOS_ERROR_IF_NOT(rParameters.Has("knot_vector"))
            << "Missing 'knot_vector' in nurbs curve" << std::endl;
        Vector knot_vector = rParameters["knot_vector"].GetVector();

        KRATOS_ERROR_IF_NOT(rParameters.Has("degree"))
            << "Missing 'degree' in nurbs curve" << std::endl;
        int polynomial_degree = rParameters["degree"].GetInt();

        PointerVector<TThisNodeType> control_points;

        ReadControlPointVector(control_points,
            rParameters["control_points"], rModelPart, EchoLevel);

        if (is_rational)
        {
            Vector control_point_weights = ReadControlPointWeightVector(
                rParameters["control_points"]);

            return Kratos::make_shared<NurbsCurveGeometry<TWorkingSpaceDimension, PointerVector<TThisNodeType>>>(
                NurbsCurveGeometry<TWorkingSpaceDimension, PointerVector<TThisNodeType>>(
                    control_points,
                    polynomial_degree,
                    knot_vector,
                    control_point_weights));
        }
        typename NurbsCurveGeometry<TWorkingSpaceDimension, PointerVector<TThisNodeType>>::Pointer p_nurbs_curve(
            new NurbsCurveGeometry<TWorkingSpaceDimension, PointerVector<TThisNodeType>>(
                control_points,
                polynomial_degree,
                knot_vector));

        return p_nurbs_curve;
    }

    /* @brief read NurbsSurfaces from the given parameter input.
    *
    * The input needs to be provided in the following shape:
    * {
    *     "is_trimmed": bool,
    *     "is_rational" : bool,
    *     "degrees" : [int, int] ,
    *     "knot_vectors" : [
    *         [ double, ... ],
    *             [double, ...]
    *     ],
    *     "control_points" : [
    *         [ id:p, [0, 10, 0, 1] ],
    *         ...
    *     ]
    * }
    */
    template<int TWorkingSpaceDimension, class TThisNodeType>
    static typename NurbsSurfaceGeometry<TWorkingSpaceDimension, PointerVector<TThisNodeType>>::Pointer
        ReadNurbsSurface(
            const Parameters rParameters,
            ModelPart& rModelPart,
            SizeType EchoLevel = 0)
    {
        bool is_rational = true;
        if (rParameters.Has("is_rational")) {
            is_rational = rParameters["is_rational"].GetBool();
        }
        else {
            KRATOS_INFO_IF("ReadNurbsSurface", (EchoLevel > 4))
                << "\"is_rational\" is not provided within \"surface\". Thus, it is considered as rational. "
                << "If this surface is non-rational the computation of the shape functions is less optimized."
                << std::endl;
        }

        KRATOS_ERROR_IF_NOT(rParameters.Has("knot_vectors"))
            << "Missing 'knot_vector' in nurbs surface" << std::endl;
        KRATOS_ERROR_IF(rParameters["knot_vectors"].size() != 2)
            << "'knot_vectors' need to be of size two, knot_vector_u and knot_vector_v" << std::endl;
        Vector knot_vector_u = rParameters["knot_vectors"][0].GetVector();
        Vector knot_vector_v = rParameters["knot_vectors"][1].GetVector();

        KRATOS_ERROR_IF_NOT(rParameters.Has("degrees"))
            << "Missing 'degrees' in nurbs surface" << std::endl;
        KRATOS_ERROR_IF(rParameters["degrees"].size() != 2)
            << "'degrees' need to be of size two, p and q" << std::endl;
        int p = rParameters["degrees"][0].GetInt();
        int q = rParameters["degrees"][1].GetInt();

        PointerVector<TThisNodeType> control_points;

        ReadControlPointVector(control_points,
            rParameters["control_points"], rModelPart, EchoLevel);

        if (is_rational)
        {
            Vector control_point_weights = ReadControlPointWeightVector(
                rParameters["control_points"]);

            return Kratos::make_shared<NurbsSurfaceGeometry<TWorkingSpaceDimension, PointerVector<TThisNodeType>>>(
                control_points,
                p,
                q,
                knot_vector_u,
                knot_vector_v,
                control_point_weights);
        }
        typename NurbsSurfaceGeometry<TWorkingSpaceDimension, PointerVector<TThisNodeType>>::Pointer p_nurbs_surface(
            new NurbsSurfaceGeometry<TWorkingSpaceDimension, PointerVector<TThisNodeType>>(
                control_points,
                p,
                q,
                knot_vector_u,
                knot_vector_v));

        return p_nurbs_surface;
    }

    ///@}
    ///@name Read in Control Points
    ///@{

    /// Reads the weights of all control points and provides them in a Vector.
    static Vector ReadControlPointWeightVector(
        const Parameters rParameters,
        SizeType EchoLevel = 0)
    {
        Vector control_point_weights = ZeroVector(rParameters.size());
        KRATOS_ERROR_IF(rParameters.size() == 0)
            << "Length of control point list is zero!" << std::endl;

        SizeType number_of_entries = rParameters[0].size();
        KRATOS_ERROR_IF(rParameters[0][number_of_entries - 1].size() != 4)
            << "Control points need to be provided in following structure: [[x, y, z, weight]] or [id, [x, y, z, weight]]"
            << "Size of inner vector incorrect!"
            << std::endl;

        for (IndexType cp_idx = 0; cp_idx < rParameters.size(); cp_idx++)
        {
            control_point_weights[cp_idx] = rParameters[cp_idx][number_of_entries - 1][3].GetDouble();
        }

        return control_point_weights;
    }

    /// Reads a Node::Pointer-vector of control points.
    static void ReadControlPointVector(
        PointerVector<Node>& rControlPoints,
        const Parameters rParameters,
        ModelPart& rModelPart,
        SizeType EchoLevel = 0)
    {
        KRATOS_ERROR_IF_NOT(rParameters.IsArray())
            << "\"control_points\" section needs to be an array." << std::endl;

        KRATOS_INFO_IF("ReadControlPointVector", EchoLevel > 4)
            << "Reading " << rParameters.size() << " control points of type Node." << std::endl;

        for (IndexType cp_idx = 0; cp_idx < rParameters.size(); cp_idx++)
        {
            rControlPoints.push_back(ReadAndCreateNode(rParameters[cp_idx], rModelPart));
        }
    }

    /// Reads a Point::Pointer-vector of control points.
    static void ReadControlPointVector(
        PointerVector<Point>& rControlPoints,
        const Parameters rParameters,
        ModelPart& rModelPart,
        SizeType EchoLevel = 0)
    {
        KRATOS_ERROR_IF_NOT(rParameters.IsArray())
            << "\"control_points\" section needs to be an array." << std::endl;

        KRATOS_INFO_IF("ReadControlPointVector", EchoLevel > 4)
            << "Reading " << rParameters.size() << " control points of type Point." << std::endl;

        for (IndexType cp_idx = 0; cp_idx < rParameters.size(); cp_idx++)
        {
            rControlPoints.push_back(ReadPoint(rParameters[cp_idx]));
        }
    }

    ///@}
    ///@name Read in Nodes/ Points
    ///@{

    /* Reads, and returns a Pointer to Node.
    * Input needs to be a Parameter object:
    * [id, [x, y, z, weight]]
    */
    static Node::Pointer ReadAndCreateNode(
        const Parameters rParameters,
        ModelPart& rModelPart,
        SizeType EchoLevel = 0)
    {
        SizeType number_of_entries = rParameters.size();
        KRATOS_ERROR_IF((number_of_entries != 2))
            << "Control points as Node need to be provided in following structure: [id, [x, y, z, weight]]"
            << std::endl;

        IndexType id = rParameters[0].GetInt();
        Vector cp = rParameters[1].GetVector();

        return rModelPart.CreateNewNode(id, cp[0], cp[1], cp[2]);
    }

    /* Reads, and returns a Pointer to Point.
    * Input needs to be a Parameter object:
    * [[x, y, z, weight]] or [id, [x, y, z, weight]]
    */
    static Point::Pointer ReadPoint(
        const Parameters rParameters,
        SizeType EchoLevel = 0)
    {
        SizeType number_of_entries = rParameters.size();
        KRATOS_ERROR_IF((number_of_entries != 1) && (number_of_entries != 2))
            << "Control points as Point need to be provided in following structure: "
            << "[[x, y, z, weight]] or [id, [x, y, z, weight]]" << std::endl;

        Vector cp = rParameters[number_of_entries - 1].GetVector();

        return Kratos::make_shared<Point>(cp[0], cp[1], cp[2]);
    }

    ///@}
    ///@name Utility functions
    ///@{

    /// Sets the geometry Id with either the 'brep_id' or the 'brep_name'.
    template<class TGeometry>
    static void SetIdOrName(
        const Parameters rParameters,
        typename TGeometry::Pointer pGeometry)
    {
        if (rParameters.Has("brep_id")) {
            pGeometry->SetId(rParameters["brep_id"].GetInt());
        }
        else if (rParameters.Has("brep_name")) {
            pGeometry->SetId(rParameters["brep_name"].GetString());
        }
    }

    /// Sets the geometry Id with either the 'brep_id' or the 'brep_name'.
    static void SetIdOrName(
        const Parameters rParameters,
        GeometryPointerType pGeometry)
    {
        if (rParameters.Has("brep_id")) {
            pGeometry->SetId(rParameters["brep_id"].GetInt());
        }
        else if (rParameters.Has("brep_name")) {
            pGeometry->SetId(rParameters["brep_name"].GetString());
        }
    }

    /// Returns the string of either the 'brep_id' or the 'brep_name'. Used for output massages.
    static std::string GetIdOrName(
        const Parameters rParameters)
    {
        if (rParameters.Has("brep_id")) {
            return std::to_string(rParameters["brep_id"].GetInt());
        }
        else if (rParameters.Has("brep_name")) {
            return rParameters["brep_name"].GetString();
        }
        else {
            return "no_id_assigned";
        }
    }

    /// Checks if one of the 'brep_id' or the 'brep_name' is provided.
    static bool HasIdOrName(
        const Parameters rParameters)
    {
        return (rParameters.Has("brep_id") || rParameters.Has("brep_name"));
    }

    /// Returns the geometry with either the 'brep_id' or the 'brep_name'.
    static typename GeometryType::Pointer GetGeometry(
        const Parameters rParameters,
        ModelPart& rModelPart)
    {
        if (rParameters.Has("brep_id")) {
            return rModelPart.pGetGeometry(rParameters["brep_id"].GetInt());
        }
        else { // if (rParameters["topology"][i].Has("brep_name"))
            return rModelPart.pGetGeometry(rParameters["brep_name"].GetString());
        }
    }

    /// Reads in a json formatted file and returns its KratosParameters instance.
    static Parameters ReadParamatersFile(
        const std::string& rDataFileName,
        SizeType EchoLevel = 0)
    {
        // Check if rDataFileName ends with ".cad.json" and add it if needed.
        const std::string data_file_name = (rDataFileName.compare(rDataFileName.size() - 9, 9, ".cad.json") != 0)
            ? rDataFileName + ".cad.json"
            : rDataFileName;

        std::ifstream infile(data_file_name);
        KRATOS_ERROR_IF_NOT(infile.good()) << "CAD geometry file: "
            << data_file_name << " cannot be found." << std::endl;
        KRATOS_INFO_IF("ReadParamatersFile", EchoLevel > 3)
            << "Reading file: \"" << data_file_name << "\"" << std::endl;

        std::stringstream buffer;
        buffer << infile.rdbuf();

        return Parameters(buffer.str());
    }

    ///@}
    ///@name Members
    ///@{

    Parameters mCadJsonParameters;
    int mEchoLevel;

    ///@}
}; // Class CadJsonInput
}  // namespace Kratos.

#endif // KRATOS_CAD_JSON_INPUT_INCLUDED  defined
