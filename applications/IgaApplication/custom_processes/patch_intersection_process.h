//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli

#pragma once

#include "processes/process.h"
#include "includes/model_part.h"
#include "geometries/geometry.h"
#include "integration/integration_info.h"
#include "geometries/nurbs_curve_geometry.h"
#include "geometries/brep_curve_on_surface.h"
#include "geometries/brep_surface.h"
#include "geometries/nurbs_surface_geometry.h"
#include "includes/properties.h"

namespace Kratos {

class PatchIntersectionProcess : public Process {
public:
    KRATOS_CLASS_POINTER_DEFINITION(PatchIntersectionProcess);

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef Node NodeType;

    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointerType;
    typedef typename Properties::Pointer PropertiesPointerType;

    // typedef NurbsSurfaceGeometry<3, PointerVector<NodeType>> NurbsSurfaceGeometryType;
    // typedef typename NurbsSurfaceGeometryType::Pointer NurbsSurfaceGeometryPointerType;

    // Param curve in UV space controlled by Nodes (Nodes store (u,v,0))
    using ParamCurveType        = NurbsCurveGeometry<2, PointerVector<Node>>;
    using BrepCurveOnSurfaceType= BrepCurveOnSurface<PointerVector<Node>, false, PointerVector<Node>>;

    // Your surface type alias (as you already had)
    using NurbsSurfaceGeometryType        = NurbsSurfaceGeometry<3, PointerVector<Node>>;
    using NurbsSurfaceGeometryPointerType = NurbsSurfaceGeometryType::Pointer;

    using NurbsCurveUVType        = Kratos::NurbsCurveGeometry<2, Kratos::PointerVector<NodeType>>;
    using IPArray = GeometryType::IntegrationPointsArrayType;



    // Constructor knobs:
    // - echo_level: 0 quiet, 1 key matches, 2 verbose debug
    // - tolerance: geometric comparisons in param-space
    // - patch_prefix: prefix of patch submodel parts (e.g. "Patch")
    // - internal_sub/external_sub: submodel part names used inside each patch
    PatchIntersectionProcess(
        ModelPart& rModelPart,
        int echo_level = 0,
        double tolerance = 1e-12,
        std::string patch_prefix = "Patch",
        std::string internal_sub = "internal_boundaries",
        std::string external_sub = "external_boundaries");

    void Execute() override;

private:
    using GeometriesArray   = GeometryType::GeometriesArrayType;

    ModelPart&   mrModelPart;
    int          mEchoLevel{0};
    double       mTol{1e-12};
    std::string  mPatchPrefix{"Patch"};
    std::string  mInternalSubName{"internal_boundaries"};
    std::string  mExternalSubName{"external_boundaries"};

    void ComputeIntersections();

    // Create body-fitted conditions (quadrature-point geometries) on a set of
    // parametric sub-intervals along an existing BREP edge.
    static void CreateBodyFittedConditionsOnSegments(
        GeometryType::Pointer pEdge,
        const std::vector<std::pair<double,double>>& segments, // local param intervals on the edge
        ModelPart& rTargetSubModelPart,                        // where to store the created conditions
        const std::string& rConditionName,                     // registered Kratos condition name
        SizeType& rIdCounter,                                  // running id (updated in-place)
        const Vector& rKnotSpanSizes,                          // to store mesh size info on the condition
        SizeType shape_func_deriv_order = 1,
        const std::string& quadrature_method = "GAUSS",
        int nips_per_span = -1);

    // Creates the two BrepCurveOnSurface pieces on a rectangle edge (helper).
    void CreateAndAddBrepCurve(
        const NurbsSurfaceGeometryPointerType pSurfaceGeometry,
        const Point& rCoordsA,
        const Point& rCoordsB,
        IndexType& rLastGeometryId,
        ModelPart& rModelPart
    );

    // Creates body-fitted conditions from the given BrepCurveOnSurface geometries.
    // - Uses GAUSS by default (1D), derivative order = 1 (same as IgaModelerSbm default)
    // - Sets KNOT_SPAN_SIZES on each created condition
    void CreateConditionsFromBrepCurves(
        const std::vector<GeometryPointerType>& rBrepCurves,
        ModelPart& rTargetSubModelPart,
        const std::string& rConditionName,
        SizeType ip_per_span = -1 /*optional: if >0, uses GRID with that many IP/span*/);

    void CreateConditionsFromBrepCurvesWithMirroredNeighbours(
        const std::vector<GeometryPointerType>& rBrepCurvesPrimary,
        const std::vector<GeometryPointerType>& rBrepCurvesMirror,
        ModelPart& rTargetSubModelPart,
        const std::string& rConditionName,
        SizeType ip_per_span,          // force same #IPs on both sides
        SizeType deriv_order = 1       // usually 1 for interface lines
    );


};

} // namespace Kratos
