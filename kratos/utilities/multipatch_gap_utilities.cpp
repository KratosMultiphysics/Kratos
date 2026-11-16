//     |  /           |
//     ' /   __| _` | __|  _ \   __|
//     . \  |   (   | |   (   |\__ `
//    _|\_\_|  \__,_|\__|\___/ ____/
//                    Multi-Physics
//
//   License:         BSD License
//                    Kratos default license: kratos/license.txt
//
//   Main author:     Nicolò Antonelli
//

#include <chrono>
#include <limits>
#include <map>
#include "utilities/multipatch_gap_utilities.h"
#include "utilities/gap_quadrature_utility.h"
#include "geometries/line_2d_2.h"
#include "geometries/triangle_2d_3.h"
#include "integration/triangle_gauss_legendre_integration_points.h"

namespace Kratos {
namespace {
// TODO: to remove temporary timing instrumentation after profiling.
thread_local bool profiling_enabled = false;
thread_local double projection_time = 0.0, triangulation_time = 0.0;
using ProfileClock = std::chrono::steady_clock;
using Utility = FemMultipatchGapUtilities;
using C = Utility::Coordinates;
using Polygon = std::vector<C>;
using Key = std::array<ModelPart::IndexType, 2>;
struct Edge { Node::Pointer pStart, pEnd; Element::Pointer pParent; };
struct EdgePair { std::vector<Edge> Sides; };
using Edges = std::map<Key, EdgePair>;

double cross(const C& rFirst, const C& rSecond, const C& rThird)
{
    return (rSecond[0]-rFirst[0])*(rThird[1]-rFirst[1]) - (rSecond[1]-rFirst[1])*(rThird[0]-rFirst[0]);
}
double distance(const C& rFirst, const C& rSecond) { return std::hypot(rFirst[0]-rSecond[0], rFirst[1]-rSecond[1]); }
double area(const Polygon& rPolygon)
{
    double value = 0.0;
    for (std::size_t i = 1; i+1<rPolygon.size(); ++i) value += cross(rPolygon[0], rPolygon[i], rPolygon[i+1]);
    return value * 0.5;
}
Edges edges(ModelPart& rPart, const bool Cartesian, const bool AllowTriangles = false)
{
    Edges result;
    for (auto& r_element : rPart.Elements()) {
        auto& r_geometry = r_element.GetGeometry();
        const bool triangle = r_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3;
        KRATOS_ERROR_IF((r_geometry.GetGeometryType() != GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4 &&
                         !(AllowTriangles && triangle)) || r_geometry.Area() <= 0.0)
            << "Gap-SBM requires positively oriented Quad4 parents or supported Triangle2D3 parents." << std::endl;
        for (const auto& r_node : r_geometry) {
            KRATOS_ERROR_IF(std::abs(r_node.Z()) > 1e-12) << "Gap-SBM requires z=0." << std::endl;
        }
        if (triangle) {
            KRATOS_ERROR_IF(cross(r_geometry[0].Coordinates(), r_geometry[1].Coordinates(),
                                  r_geometry[2].Coordinates()) <= 0.0)
                << "Gap-SBM requires counterclockwise triangular parents." << std::endl;
            const double scale = std::sqrt(r_geometry.Area());
            bool horizontal = false, vertical = false;
            for (const auto& r_edge : r_geometry.GenerateEdges()) {
                horizontal = horizontal || std::abs(r_edge[1].Y() - r_edge[0].Y()) <= 1e-10 * scale;
                vertical = vertical || std::abs(r_edge[1].X() - r_edge[0].X()) <= 1e-10 * scale;
            }
            KRATOS_ERROR_IF_NOT(horizontal && vertical)
                << "Gap-SBM triangular parents must be halves of Cartesian cells." << std::endl;
        }
        if (Cartesian && !triangle) {
            const double scale = std::sqrt(r_geometry.Area());
            KRATOS_ERROR_IF(std::abs(r_geometry[1].Y()-r_geometry[0].Y()) > 1e-10*scale ||
                            std::abs(r_geometry[2].X()-r_geometry[1].X()) > 1e-10*scale ||
                            std::abs(r_geometry[3].Y()-r_geometry[2].Y()) > 1e-10*scale ||
                            std::abs(r_geometry[0].X()-r_geometry[3].X()) > 1e-10*scale)
                << "Gap-SBM background must be Cartesian Quad4." << std::endl;
        }
        for (auto& r_edge : r_geometry.GenerateEdges()) {
            Key key{r_edge[0].Id(), r_edge[1].Id()};
            std::sort(key.begin(), key.end());
            auto& r_sides = result[key].Sides;
            r_sides.push_back({r_edge.pGetPoint(0), r_edge.pGetPoint(1), rPart.pGetElement(r_element.Id())});
            KRATOS_ERROR_IF(r_sides.size()>2) << "Non-manifold patch edge." << std::endl;
        }
    }
    return result;
}
std::vector<Edge> contour(ModelPart& rInner, const double Tol)
{
    const auto all = edges(rInner, false);
    std::map<ModelPart::IndexType, Edge> next;
    for (const auto& [r_key, r_pair] : all) {
        if (r_pair.Sides.size()==1) {
            const auto& r_edge = r_pair.Sides[0];
            KRATOS_ERROR_IF_NOT(next.emplace(r_edge.pStart->Id(), r_edge).second)
                << "Inner contour must be a single closed manifold." << std::endl;
        }
    }
    KRATOS_ERROR_IF(next.empty()) << "Empty fitted inner boundary." << std::endl;
    std::vector<Edge> result;
    auto id = next.begin()->first;
    const auto first = id;
    do {
        const auto it = next.find(id);
        KRATOS_ERROR_IF(it==next.end() || result.size()>=next.size())
            << "Invalid fitted boundary loop." << std::endl;
        result.push_back(it->second);
        id = it->second.pEnd->Id();
    } while (id!=first);
    KRATOS_ERROR_IF(result.size()!=next.size()) << "Multiple inner boundary loops are not supported." << std::endl;
    Polygon points;
    for (const auto& r_edge: result) points.push_back(r_edge.pStart->Coordinates());
    KRATOS_ERROR_IF(area(points)<=Tol*Tol) << "Inner contour must be counterclockwise." << std::endl;
    for (const auto& r_edge:result) {
        const double length = distance(r_edge.pStart->Coordinates(), r_edge.pEnd->Coordinates());
        KRATOS_ERROR_IF(length<=Tol) << "Zero-length inner segment." << std::endl;
        for (const auto& r_point:points) {
            KRATOS_ERROR_IF(cross(r_edge.pStart->Coordinates(), r_edge.pEnd->Coordinates(), r_point)<-Tol*length)
                << "First Gap-SBM implementation requires a convex fitted contour." << std::endl;
        }
    }
    return result;
}
Polygon clip(Polygon polygon, const std::vector<Edge>& rBoundary)
{
    for (const auto& r_edge:rBoundary) {
        Polygon out;
        if (polygon.empty()) break;
        const auto& r_edge_start = r_edge.pStart->Coordinates();
        const auto& r_edge_end = r_edge.pEnd->Coordinates();
        for (std::size_t i = 0; i<polygon.size(); ++i) {
            const C& r_start = polygon[i]; const C& r_end = polygon[(i+1)%polygon.size()];
            const double ds = cross(r_edge_start, r_edge_end, r_start), dt = cross(r_edge_start, r_edge_end, r_end);
            if (ds>=0.0) out.push_back(r_start);
            if ((ds<0.0 && dt>0.0)||(ds>0.0 && dt<0.0)) {
                const double w = ds/(ds-dt);
                C point = r_start+w*(r_end-r_start);
                out.push_back(point);
            }
        }
        polygon = std::move(out);
    }
    return polygon;
}
std::set<ModelPart::IndexType> cut_cells(ModelPart& rBackground,
                                        const std::vector<Edge>& rBoundary, const double Tol)
{
    std::set<ModelPart::IndexType> result;
    for (const auto& r_element:rBackground.Elements()) {
        Polygon polygon;
        for (const auto& r_node:r_element.GetGeometry()) polygon.push_back(r_node.Coordinates());
        const Polygon intersection = clip(polygon, rBoundary);
        if (intersection.size()<3 || area(intersection)<=0.0) continue;
        // A corner cut has area proportional to penetration squared: an
        // area cutoff Tol*h can miss cuts much deeper than Tol. Test an
        // interior point against both contours using distance tolerances.
        // This still excludes roundoff-sized slivers on fitted edges.
        C interior = ZeroVector(3);
        for (const auto& r_point:intersection) interior+=r_point;
        interior/=static_cast<double>(intersection.size());
        bool inside = true;
        for (std::size_t i = 0; i<polygon.size(); ++i) {
            const auto& r_start = polygon[i];
            const auto& r_end = polygon[(i+1)%polygon.size()];
            inside = inside && cross(r_start, r_end, interior)>Tol*distance(r_start, r_end);
        }
        for (const auto& r_edge:rBoundary) {
            const auto& r_start = r_edge.pStart->Coordinates();
            const auto& r_end = r_edge.pEnd->Coordinates();
            inside = inside && cross(r_start, r_end, interior)>Tol*distance(r_start, r_end);
        }
        if (inside) result.insert(r_element.Id());
    }
    return result;
}
double triangle_overlap(Polygon polygon, const std::array<C, 3>& rTriangle)
{
    for (std::size_t k = 0; k<3 && !polygon.empty(); ++k) {
        Polygon out;
        for (std::size_t i = 0; i<polygon.size(); ++i) {
            const C& r_start = polygon[i]; const C& r_end = polygon[(i+1)%polygon.size()];
            const double da = cross(rTriangle[k], rTriangle[(k+1)%3], r_start);
            const double db = cross(rTriangle[k], rTriangle[(k+1)%3], r_end);
            if (da>=0.0) out.push_back(r_start);
            if ((da<0.0 && db>0.0)||(da>0.0 && db<0.0)) {
                const C point = r_start+(da/(da-db))*(r_end-r_start); out.push_back(point);
            }
        }
        polygon = std::move(out);
    }
    return std::abs(area(polygon));
}
struct Projection { C Point; double S = 0.0; };
// Avoid the name "project": ADL also finds Boost uBLAS matrix projections
// through array_1d, causing invalid matrix-proxy template instantiations.
Projection project_to_contour(const C& rPoint, const std::vector<Edge>& rBoundary, const std::vector<double>& rLengths)
{
    double best = std::numeric_limits<double>::max();
    Projection result;
    for (std::size_t i = 0; i<rBoundary.size(); ++i) {
        const C a = rBoundary[i].pStart->Coordinates();
        const C d = rBoundary[i].pEnd->Coordinates()-a;
        const double length = rLengths[i+1]-rLengths[i];
        const double t = std::clamp(((rPoint[0]-a[0])*d[0]+(rPoint[1]-a[1])*d[1])/(length*length), 0.0, 1.0);
        const C q = a+t*d;
        const double dist = distance(rPoint, q);
        if (dist<best) { best = dist; result = {q, rLengths[i]+t*length}; }
    }
    if (result.S>=rLengths.back()) result.S = 0.0;
    return result;
}
std::vector<std::array<C, 3>> triangulate(Polygon polygon, const double Tol)
{
    // Remove duplicate/collinear points before ordinary ear clipping.
    bool changed = true;
    while (changed && polygon.size()>2) {
        changed = false;
        for (std::size_t i = 0; i<polygon.size(); ++i) {
            const auto& r_previous = polygon[(i+polygon.size()-1)%polygon.size()];
            const auto& r_next = polygon[(i+1)%polygon.size()];
            if (distance(r_previous, polygon[i])<=Tol || std::abs(cross(r_previous, polygon[i], r_next))<=Tol*distance(r_previous, r_next)) {
                polygon.erase(polygon.begin()+i); changed = true; break;
            }
        }
    }
    std::vector<std::array<C, 3>> result;
    if (polygon.size()<3 || std::abs(area(polygon))<=Tol*Tol) return result;
    if (area(polygon)<0.0) std::reverse(polygon.begin(), polygon.end());
    while (polygon.size()>3) {
        bool found = false;
        for (std::size_t i = 0; i<polygon.size(); ++i) {
            const std::size_t a = (i+polygon.size()-1)%polygon.size(), b = i, c = (i+1)%polygon.size();
            if (cross(polygon[a], polygon[b], polygon[c])<=Tol*Tol) continue;
            bool contains = false;
            for (std::size_t j = 0; j<polygon.size(); ++j) {
                if (j==a||j==b||j==c) continue;
                if (cross(polygon[a], polygon[b], polygon[j])>=-Tol*Tol && cross(polygon[b], polygon[c], polygon[j])>=-Tol*Tol &&
                    cross(polygon[c], polygon[a], polygon[j])>=-Tol*Tol) { contains = true; break; }
            }
            if (contains) continue;
            result.push_back({polygon[a], polygon[b], polygon[c]}); polygon.erase(polygon.begin()+b); found = true; break;
        }
        KRATOS_ERROR_IF_NOT(found) << "Gap polygon is not simple; invalid projection topology." << std::endl;
    }
    if (cross(polygon[0], polygon[1], polygon[2])>Tol*Tol) result.push_back({polygon[0], polygon[1], polygon[2]});
    return result;
}

// Try a fan from either endpoint of the surrogate edge. Each accepted
// triangle then has a parent node; interface edges remain polygon edges.
// A height-function gap need not be star-shaped from either endpoint.
std::vector<std::array<C, 3>> triangulate_level_set_gap(const Polygon& rPolygon, const double Tol)
{
    if (rPolygon.size() < 3 || std::abs(area(rPolygon)) <= Tol * Tol) {
        return triangulate(rPolygon, Tol);
    }
    const double orientation = area(rPolygon) > 0.0 ? 1.0 : -1.0;
    std::vector<std::array<C, 3>> best;
    double best_quality = -1.0;
    for (std::size_t anchor = 0; anchor < 2; ++anchor) {
        std::vector<std::array<C, 3>> candidate;
        double minimum_quality = std::numeric_limits<double>::max();
        double candidate_area = 0.0;
        bool valid = true;
        for (std::size_t offset = 1; offset + 1 < rPolygon.size(); ++offset) {
            const C& r_a = rPolygon[anchor];
            const C& r_b = rPolygon[(anchor + offset) % rPolygon.size()];
            const C& r_c = rPolygon[(anchor + offset + 1) % rPolygon.size()];
            const double determinant = orientation * cross(r_a, r_b, r_c);
            if (determinant < -Tol * Tol) {
                valid = false;
                break;
            }
            if (determinant <= Tol * Tol) continue;
            const double ab = distance(r_a, r_b), bc = distance(r_b, r_c), ca = distance(r_c, r_a);
            minimum_quality = std::min(minimum_quality, determinant / (ab * ab + bc * bc + ca * ca));
            candidate_area += 0.5 * determinant;
            candidate.push_back(orientation > 0.0 ? std::array<C, 3>{r_a, r_b, r_c}
                                                   : std::array<C, 3>{r_a, r_c, r_b});
        }
        // Near-degenerate pieces are treated as in the existing ear clipper;
        // reject a fan if they would change the polygon area appreciably.
        const double area_tolerance = rPolygon.size() * Tol * Tol + 1e-12 * std::abs(area(rPolygon));
        if (valid && !candidate.empty() && std::abs(candidate_area - std::abs(area(rPolygon))) <= area_tolerance &&
            minimum_quality > best_quality) {
            best_quality = minimum_quality;
            best = std::move(candidate);
        }
    }
    return best.empty() ? triangulate(rPolygon, Tol) : best;
}

Geometry<Node>::Pointer quadrature(Element& rParent, const std::array<C, 3>& rTriangle)
{
    Geometry<Node>::IntegrationPointsArrayType points;
    const double determinant = cross(rTriangle[0], rTriangle[1], rTriangle[2]);
    // The reference triangle rule integrates polynomials through degree two.
    for (const auto& r_point : TriangleGaussLegendreIntegrationPoints2::IntegrationPoints()) {
        const C global_point = rTriangle[0] + r_point.X() * (rTriangle[1] - rTriangle[0])
            + r_point.Y() * (rTriangle[2] - rTriangle[0]);
        points.emplace_back(global_point[0], global_point[1], 0.0, determinant * r_point.Weight());
    }
    return GapQuadratureUtility::Create(rParent.GetGeometry(), points);
}
} // unnamed namespace

FemMultipatchUtilities::OwnerMapType FemMultipatchGapUtilities::UpdateActivation(
    ModelPart& rBackground, const std::vector<ModelPart*>& rInnerPatches, const double Tolerance)
{
    KRATOS_TRY
    KRATOS_ERROR_IF(!std::isfinite(Tolerance)||Tolerance<=0.0) << "Invalid gap tolerance." << std::endl;
    const auto all_edges = edges(rBackground, true);
    KRATOS_ERROR_IF(rBackground.NumberOfElements()==0) << "Empty background." << std::endl;
    double xmin = std::numeric_limits<double>::max(), ymin = xmin;
    double xmax = -xmin, ymax = -ymin;
    for (const auto& r_node:rBackground.Nodes()) {
        xmin = std::min(xmin, r_node.X()); ymin = std::min(ymin, r_node.Y());
        xmax = std::max(xmax, r_node.X()); ymax = std::max(ymax, r_node.Y());
    }
    FemMultipatchUtilities::OwnerMapType owners;
    std::set<ModelPart::IndexType> independent;
    for (const auto& r_node:rBackground.Nodes()) independent.insert(r_node.Id());
    for (std::size_t i = 0; i<rInnerPatches.size(); ++i) {
        KRATOS_ERROR_IF(&rInnerPatches[i]->GetRootModelPart()!=&rBackground.GetRootModelPart())
            << "Patches must share a root ModelPart." << std::endl;
        for (const auto& r_node:rInnerPatches[i]->Nodes()) {
            KRATOS_ERROR_IF(r_node.X()<=xmin+Tolerance || r_node.X()>=xmax-Tolerance ||
                            r_node.Y()<=ymin+Tolerance || r_node.Y()>=ymax-Tolerance)
                << "Inner patches must lie strictly inside the background." << std::endl;
            KRATOS_ERROR_IF_NOT(independent.insert(r_node.Id()).second) << "Patches must have independent nodes." << std::endl;
        }
        const auto boundary = contour(*rInnerPatches[i], Tolerance);
        for (const auto id:cut_cells(rBackground, boundary, Tolerance)) {
            KRATOS_ERROR_IF_NOT(owners.emplace(id, i+1).second)
                << "Inner patches or their conservative surrogate holes overlap." << std::endl;
        }
    }
    for (const auto& [r_key, r_pair]:all_edges) {
        if (r_pair.Sides.size()==1) {
            KRATOS_ERROR_IF(owners.count(r_pair.Sides[0].pParent->Id()))
                << "Surrogate holes must lie strictly inside the background." << std::endl;
        } else {
            const auto a = owners.find(r_pair.Sides[0].pParent->Id()), b = owners.find(r_pair.Sides[1].pParent->Id());
            KRATOS_ERROR_IF(a!=owners.end() && b!=owners.end() && a->second!=b->second)
                << "Surrogate holes must not touch." << std::endl;
        }
    }
    // Reject vertex contacts too, including non-manifold diagonal holes.
    std::map<ModelPart::IndexType, std::set<std::size_t>> node_owners;
    for (const auto& r_element:rBackground.Elements()) {
        const auto it = owners.find(r_element.Id());
        if (it!=owners.end()) for (const auto& r_node:r_element.GetGeometry()) node_owners[r_node.Id()].insert(it->second);
    }
    for (const auto& [r_id, r_owners]:node_owners) KRATOS_ERROR_IF(r_owners.size()>1) << "Surrogate holes touch at a vertex." << std::endl;
    // Commit flags only after classification/ownership validation succeeds.
    for (auto& r_element:rBackground.Elements()) r_element.Set(ACTIVE, owners.count(r_element.Id())==0);
    return owners;
    KRATOS_CATCH("")
}

void FemMultipatchGapUtilities::ResetProfilingTimes() { projection_time = triangulation_time = 0.0; }
std::array<double, 2> FemMultipatchGapUtilities::GetProfilingTimes() { return {projection_time, triangulation_time}; }
void FemMultipatchGapUtilities::EnableProfiling(bool Enabled) { profiling_enabled = Enabled; }

FemMultipatchGapUtilities::GapData FemMultipatchGapUtilities::BuildLevelSetPhase(
    ModelPart& rPhase, const Matrix& rContour, const bool NegativeSide, const double Tolerance)
{
    KRATOS_TRY
    KRATOS_ERROR_IF(!std::isfinite(Tolerance) || Tolerance <= 0.0)
        << "Invalid level-set gap tolerance." << std::endl;
    KRATOS_ERROR_IF(rPhase.NumberOfElements() == 0 || rContour.size1() < 2 || rContour.size2() != 2)
        << "Level-set gap requires phase elements and an Nx2 contour." << std::endl;
    Polygon contour_points;
    contour_points.reserve(rContour.size1());
    for (std::size_t i = 0; i < rContour.size1(); ++i) {
        KRATOS_ERROR_IF(!std::isfinite(rContour(i, 0)) || !std::isfinite(rContour(i, 1)))
            << "Non-finite level-set contour point." << std::endl;
        KRATOS_ERROR_IF(i && rContour(i, 0) <= rContour(i - 1, 0) + Tolerance)
            << "Level-set contour must be a strictly increasing-x height function." << std::endl;
        contour_points.push_back(C{rContour(i, 0), rContour(i, 1), 0.0});
    }
    const auto height = [&](const double x) {
        KRATOS_ERROR_IF(x < contour_points.front()[0] - Tolerance || x > contour_points.back()[0] + Tolerance)
            << "Level-set contour does not span the phase patch." << std::endl;
        const auto it = std::lower_bound(contour_points.begin(), contour_points.end(), x,
            [](const C& rPoint, double value) { return rPoint[0] < value; });
        if (it == contour_points.begin()) return (*it)[1];
        if (it == contour_points.end()) return contour_points.back()[1];
        const auto& r_first = *(it - 1); const auto& r_second = *it;
        return r_first[1] + (x - r_first[0]) * (r_second[1] - r_first[1]) / (r_second[0] - r_first[0]);
    };
    const auto all_edges = edges(rPhase, true, true);
    std::vector<Edge> boundary;
    for (const auto& [r_key, r_pair] : all_edges) {
        if (r_pair.Sides.size() != 1) continue;
        const auto& r_edge = r_pair.Sides.front();
        const double dx = r_edge.pEnd->X() - r_edge.pStart->X();
        if (NegativeSide ? dx >= -Tolerance : dx <= Tolerance) continue;
        const double xmin = std::min(r_edge.pStart->X(), r_edge.pEnd->X());
        const double xmax = std::max(r_edge.pStart->X(), r_edge.pEnd->X());
        const auto edge_height = [&](double x) {
            return r_edge.pStart->Y() + (x - r_edge.pStart->X()) *
                (r_edge.pEnd->Y() - r_edge.pStart->Y()) / dx;
        };
        const auto valid_height = [&](double x, double level) {
            const double y = edge_height(x);
            KRATOS_ERROR_IF(NegativeSide ? y > level + Tolerance : y < level - Tolerance)
                << "Phase surrogate edge crosses the level-set interface at x=" << x << "." << std::endl;
        };
        valid_height(xmin, height(xmin));
        valid_height(xmax, height(xmax));
        for (const auto& r_point : contour_points)
            if (r_point[0] > xmin && r_point[0] < xmax) valid_height(r_point[0], r_point[1]);
        boundary.push_back(r_edge);
    }
    std::sort(boundary.begin(), boundary.end(), [](const Edge& rFirst, const Edge& rSecond) {
        return std::min(rFirst.pStart->X(), rFirst.pEnd->X()) < std::min(rSecond.pStart->X(), rSecond.pEnd->X());
    });
    KRATOS_ERROR_IF(boundary.empty()) << "No level-set surrogate boundary found." << std::endl;
    KRATOS_ERROR_IF(std::abs(std::min(boundary.front().pStart->X(), boundary.front().pEnd->X()) - contour_points.front()[0]) > Tolerance ||
                    std::abs(std::max(boundary.back().pStart->X(), boundary.back().pEnd->X()) - contour_points.back()[0]) > Tolerance)
        << "Level-set contour must connect both lateral walls." << std::endl;
    GapData result;
    double gap_area = 0.0;
    for (std::size_t i = 0; i < boundary.size(); ++i) {
        const auto& r_edge = boundary[i];
        const double xmin = std::min(r_edge.pStart->X(), r_edge.pEnd->X());
        const double xmax = std::max(r_edge.pStart->X(), r_edge.pEnd->X());
        if (i) {
            const auto& r_previous = boundary[i - 1];
            KRATOS_ERROR_IF(std::abs(xmin - std::max(r_previous.pStart->X(), r_previous.pEnd->X())) > Tolerance)
                << "Level-set phase requires one contiguous surrogate edge per column." << std::endl;
            const double y_left = r_previous.pStart->X() > r_previous.pEnd->X() ?
                r_previous.pStart->Y() : r_previous.pEnd->Y();
            const double y_right = r_edge.pStart->X() < r_edge.pEnd->X() ?
                r_edge.pStart->Y() : r_edge.pEnd->Y();
            const double y_interface = height(xmin);
            const double low = std::min(y_left, y_right), high = std::max(y_left, y_right);
            const double overlap_low = NegativeSide ? high : y_interface;
            const double overlap_high = NegativeSide ? y_interface : low;
            if (overlap_high - overlap_low > Tolerance) {
                result.Connectors.push_back({r_previous.pParent, r_edge.pParent,
                    C{xmin, overlap_low, 0.0}, C{xmin, overlap_high, 0.0}});
            }
            if (high - low > Tolerance) {
                // The taller surrogate column is bulk on one side and a gap
                // extension on the other. Split the seam at every bulk edge.
                const bool bulk_left = NegativeSide ? y_left > y_right : y_left < y_right;
                double covered = 0.0;
                for (auto& r_element : rPhase.Elements()) {
                    const auto& r_geometry = r_element.GetGeometry();
                    const auto center = r_geometry.Center();
                    if (bulk_left ? center.X() >= xmin - Tolerance : center.X() <= xmin + Tolerance) continue;
                    for (const auto& r_face : r_geometry.GenerateEdges()) {
                        if (std::abs(r_face[0].X() - xmin) > Tolerance ||
                            std::abs(r_face[1].X() - xmin) > Tolerance) continue;
                        const double bottom = std::max(low, std::min(r_face[0].Y(), r_face[1].Y()));
                        const double top = std::min(high, std::max(r_face[0].Y(), r_face[1].Y()));
                        if (top - bottom <= Tolerance) continue;
                        auto p_bulk = rPhase.pGetElement(r_element.Id());
                        result.Connectors.push_back({bulk_left ? p_bulk : r_previous.pParent,
                            bulk_left ? r_edge.pParent : p_bulk, C{xmin, bottom, 0.0}, C{xmin, top, 0.0}});
                        covered += top - bottom;
                    }
                }
                KRATOS_ERROR_IF(std::abs(covered - (high - low)) > Tolerance)
                    << "Level-set surrogate staircase has missing or overlapping bulk neighbours." << std::endl;
            }
        }
        std::vector<C> line{C{xmin, height(xmin), 0.0}};
        for (const auto& r_point : contour_points)
            if (r_point[0] > xmin + Tolerance && r_point[0] < xmax - Tolerance) line.push_back(r_point);
        line.push_back(C{xmax, height(xmax), 0.0});
        result.Surrogate.push_back({r_edge.pParent, nullptr, r_edge.pStart->Coordinates(), r_edge.pEnd->Coordinates()});
        const double y_at_min = r_edge.pStart->X() < r_edge.pEnd->X() ? r_edge.pStart->Y() : r_edge.pEnd->Y();
        const double y_at_max = r_edge.pStart->X() > r_edge.pEnd->X() ? r_edge.pStart->Y() : r_edge.pEnd->Y();
        Polygon polygon{C{xmin, y_at_min, 0.0}, C{xmax, y_at_max, 0.0}};
        for (auto it = line.rbegin(); it != line.rend(); ++it) polygon.push_back(*it);
        for (const auto& r_triangle : triangulate_level_set_gap(polygon, Tolerance)) {
            gap_area += 0.5 * cross(r_triangle[0], r_triangle[1], r_triangle[2]);
            result.Cells.push_back({r_edge.pParent, r_triangle, quadrature(*r_edge.pParent, r_triangle)});
        }
        for (std::size_t j = 1; j < line.size(); ++j) {
            result.Physical.push_back({r_edge.pParent, nullptr,
                NegativeSide ? line[j] : line[j - 1], NegativeSide ? line[j - 1] : line[j]});
        }
    }
    double bulk_area = 0.0;
    double outer_y = NegativeSide ? std::numeric_limits<double>::max() : -std::numeric_limits<double>::max();
    for (const auto& r_element : rPhase.Elements()) {
        bulk_area += r_element.GetGeometry().Area();
        for (const auto& r_node : r_element.GetGeometry())
            outer_y = NegativeSide ? std::min(outer_y, r_node.Y()) : std::max(outer_y, r_node.Y());
    }
    double physical_area = 0.0;
    for (std::size_t i = 1; i < contour_points.size(); ++i) {
        const double width = contour_points[i][0] - contour_points[i - 1][0];
        const double mean_height = 0.5 * (contour_points[i][1] + contour_points[i - 1][1]);
        physical_area += width * (NegativeSide ? mean_height - outer_y : outer_y - mean_height);
    }
    KRATOS_ERROR_IF(std::abs(bulk_area + gap_area - physical_area) >
                    100.0 * Tolerance * std::max(1.0, std::abs(physical_area)))
        << "Level-set phase bulk and gap do not cover their physical area." << std::endl;
    return result;
    KRATOS_CATCH("")
}

FemMultipatchGapUtilities::GapData FemMultipatchGapUtilities::BuildHorizontalPhase(
    ModelPart& rPhase, const double InterfaceHeight, const bool NegativeSide, const double Tolerance)
{
    KRATOS_TRY
    KRATOS_ERROR_IF(!std::isfinite(InterfaceHeight) || !std::isfinite(Tolerance) || Tolerance <= 0.0)
        << "Invalid horizontal interface or tolerance." << std::endl;
    KRATOS_ERROR_IF(rPhase.NumberOfElements() == 0) << "Empty phase patch." << std::endl;
    const auto all_edges = edges(rPhase, true);
    double boundary_y = NegativeSide ? -std::numeric_limits<double>::max() : std::numeric_limits<double>::max();
    for (const auto& r_element : rPhase.Elements()) {
        for (const auto& r_node : r_element.GetGeometry()) {
            boundary_y = NegativeSide ? std::max(boundary_y, r_node.Y()) : std::min(boundary_y, r_node.Y());
            KRATOS_ERROR_IF(NegativeSide ? r_node.Y() > InterfaceHeight + Tolerance : r_node.Y() < InterfaceHeight - Tolerance)
                << "Phase elements must lie entirely on their own side of the interface." << std::endl;
        }
    }
    std::vector<Edge> boundary;
    for (const auto& [r_key, r_pair] : all_edges) {
        if (r_pair.Sides.size() != 1) continue;
        const auto& r_edge = r_pair.Sides.front();
        if (std::abs(r_edge.pStart->Y() - boundary_y) <= Tolerance &&
            std::abs(r_edge.pEnd->Y() - boundary_y) <= Tolerance) boundary.push_back(r_edge);
    }
    std::sort(boundary.begin(), boundary.end(), [](const Edge& rFirst, const Edge& rSecond) {
        return std::min(rFirst.pStart->X(), rFirst.pEnd->X()) < std::min(rSecond.pStart->X(), rSecond.pEnd->X());
    });
    KRATOS_ERROR_IF(boundary.empty()) << "No horizontal surrogate boundary." << std::endl;
    GapData result;
    for (std::size_t i = 0; i < boundary.size(); ++i) {
        const auto& r_edge = boundary[i];
        const auto& r_start = r_edge.pStart->Coordinates();
        const auto& r_end = r_edge.pEnd->Coordinates();
        KRATOS_ERROR_IF((NegativeSide && r_start[0] <= r_end[0]) || (!NegativeSide && r_start[0] >= r_end[0]))
            << "Surrogate boundary orientation is inconsistent with phase side." << std::endl;
        if (i > 0) {
            const auto& r_previous = boundary[i-1];
            const double previous_x = std::max(r_previous.pStart->X(), r_previous.pEnd->X());
            const double current_x = std::min(r_start[0], r_end[0]);
            KRATOS_ERROR_IF(std::abs(current_x - previous_x) > Tolerance)
                << "Surrogate boundary must be one contiguous edge chain." << std::endl;
            C a = ZeroVector(3), b = ZeroVector(3);
            a[0] = b[0] = current_x; a[1] = boundary_y; b[1] = InterfaceHeight;
            if (distance(a, b) > Tolerance && r_previous.pParent != r_edge.pParent) {
                // The right-hand normal points from the left parent to the right.
                if (b[1] < a[1]) std::swap(a, b);
                result.Connectors.push_back({r_previous.pParent, r_edge.pParent, a, b});
            }
        }
        C projected_start = r_start, projected_end = r_end;
        projected_start[1] = projected_end[1] = InterfaceHeight;
        result.Surrogate.push_back({r_edge.pParent, nullptr, r_start, r_end});
        result.Physical.push_back({r_edge.pParent, nullptr, projected_start, projected_end});
        const Polygon polygon{r_start, r_end, projected_end, projected_start};
        for (const auto& r_triangle : triangulate(polygon, Tolerance)) {
            result.Cells.push_back({r_edge.pParent, r_triangle, quadrature(*r_edge.pParent, r_triangle)});
        }
    }
    return result;
    KRATOS_CATCH("")
}

FemMultipatchGapUtilities::GapData FemMultipatchGapUtilities::Build(
    ModelPart& rBackground, ModelPart& rInner, const double Tolerance)
{
    KRATOS_TRY
    KRATOS_ERROR_IF(!std::isfinite(Tolerance)||Tolerance<=0.0) << "Invalid gap tolerance." << std::endl;
    const auto boundary = contour(rInner, Tolerance);
    const auto inactive = cut_cells(rBackground, boundary, Tolerance);
    const auto all_edges = edges(rBackground, true);
    GapData result;
    std::vector<double> lengths{0.0};
    for (const auto& r_edge:boundary) lengths.push_back(lengths.back()+distance(r_edge.pStart->Coordinates(), r_edge.pEnd->Coordinates()));
    const double perimeter = lengths.back();
    std::map<ModelPart::IndexType, Projection> projections;
    std::map<ModelPart::IndexType, std::vector<std::size_t>> incident;
    std::vector<C> centers;
    double removed_area = 0.0;
    for (const auto id:inactive) {
        const auto& r_element = rBackground.GetElement(id);
        KRATOS_ERROR_IF(r_element.Is(ACTIVE)) << "Run conservative activation before building the gap." << std::endl;
        removed_area+=r_element.GetGeometry().Area();
    }
    for (const auto& [r_key, r_pair]:all_edges) {
        if (r_pair.Sides.size()!=2) continue;
        const bool a = inactive.count(r_pair.Sides[0].pParent->Id()), b = inactive.count(r_pair.Sides[1].pParent->Id());
        if (a==b) continue;
        const auto& r_edge = r_pair.Sides[a?1:0];
        KRATOS_ERROR_IF_NOT(r_edge.pParent->Is(ACTIVE)) << "Surrogate holes overlap or touch." << std::endl;
        const std::size_t index = result.Surrogate.size();
        result.Surrogate.push_back({r_edge.pParent, nullptr, r_edge.pStart->Coordinates(), r_edge.pEnd->Coordinates()});
        for (const auto& p_node:{r_edge.pStart, r_edge.pEnd}) {
            if (projections.find(p_node->Id()) == projections.end()) {
                const auto profile_start = profiling_enabled ? ProfileClock::now() : ProfileClock::time_point{};
                projections.emplace(p_node->Id(), project_to_contour(p_node->Coordinates(), boundary, lengths));
                if (profiling_enabled) projection_time += std::chrono::duration<double>(ProfileClock::now()-profile_start).count();
            }
            incident[p_node->Id()].push_back(index);
        }
        const auto p0 = projections.at(r_edge.pStart->Id()), p1 = projections.at(r_edge.pEnd->Id());
        double delta = p0.S-p1.S;
        if (delta < -Tolerance) delta+=perimeter;
        if (std::abs(delta)<Tolerance) delta = 0.0;
        Polygon polygon{r_edge.pStart->Coordinates(), r_edge.pEnd->Coordinates(), p1.Point};
        double current = p1.S;
        C previous = p1.Point;
        // Follow the actual CCW inner contour, not a chord across a corner.
        for (std::size_t j = 0; j<boundary.size()*2 && current < p1.S+delta-Tolerance; ++j) {
            std::size_t k = 0;
            double s = current;
            if (s>=perimeter-Tolerance) s-=perimeter;
            while (k+1<boundary.size() && lengths[k+1]<=s+Tolerance) ++k;
            const double offset = current-s;
            const double end = std::min(p1.S+delta, offset+lengths[k+1]);
            KRATOS_ERROR_IF(end<=current) << "Invalid contour traversal." << std::endl;
            const double fraction = (end-offset-lengths[k])/(lengths[k+1]-lengths[k]);
            const C point = boundary[k].pStart->Coordinates()+fraction*(boundary[k].pEnd->Coordinates()-boundary[k].pStart->Coordinates());
            polygon.push_back(point);
            if (distance(previous, point)>Tolerance) {
                // Reverse the inner contour: normal points from background into inner patch.
                result.Physical.push_back({r_edge.pParent, boundary[k].pParent, point, previous});
            }
            previous = point; current = end;
        }
        if (distance(polygon.back(), p0.Point)>Tolerance) polygon.push_back(p0.Point);
        C center = ZeroVector(3);
        for (const auto& r_point:polygon) center+=r_point;
        center/=static_cast<double>(polygon.size()); centers.push_back(center);
        // Tolerance is a distance. A coincident edge can acquire a tiny
        // signed area from projection roundoff, so scale the area tolerance
        // by the local edge length rather than squaring the distance.
        const double area_tolerance = Tolerance * std::max(
            distance(r_edge.pStart->Coordinates(), r_edge.pEnd->Coordinates()),
            distance(p0.Point, p1.Point));
        KRATOS_ERROR_IF(area(polygon)>area_tolerance)
            << "Gap projection reverses orientation; invalid additive gap." << std::endl;
        const auto profile_start = profiling_enabled ? ProfileClock::now() : ProfileClock::time_point{};
        const auto triangles = triangulate(polygon, Tolerance);
        if (profiling_enabled) triangulation_time += std::chrono::duration<double>(ProfileClock::now()-profile_start).count();
        for (const auto& r_triangle:triangles) {
            result.Cells.push_back({r_edge.pParent, r_triangle, quadrature(*r_edge.pParent, r_triangle)});
        }
    }
    KRATOS_ERROR_IF(result.Surrogate.empty()) << "No closed internal surrogate boundary found." << std::endl;
    for (const auto& [r_id, r_faces]:incident) {
        KRATOS_ERROR_IF(r_faces.size()!=2) << "Surrogate boundary must be a closed manifold." << std::endl;
        const auto& r_first = result.Surrogate[r_faces[0]]; const auto& r_second = result.Surrogate[r_faces[1]];
        C start = rBackground.GetNode(r_id).Coordinates(), end = projections.at(r_id).Point;
        if (r_first.pFirst==r_second.pFirst || distance(start, end)<=Tolerance) continue;
        const C difference = centers[r_faces[1]]-centers[r_faces[0]];
        if ((end[1]-start[1])*difference[0]-(end[0]-start[0])*difference[1]<0.0) std::swap(start, end);
        result.Connectors.push_back({r_first.pFirst, r_second.pFirst, start, end});
    }
    double gap_area = 0.0, interface_length = 0.0;
    for (const auto& r_cell:result.Cells) gap_area+=0.5*cross(r_cell.Vertices[0], r_cell.Vertices[1], r_cell.Vertices[2]);
    for (const auto& r_segment:result.Physical) interface_length+=distance(r_segment.Start, r_segment.End);
    Polygon inner_points;
    for (const auto& r_edge:boundary) inner_points.push_back(r_edge.pStart->Coordinates());
    const double scale = std::max(1.0, removed_area);
    std::vector<std::array<double, 4>> bounds;
    bounds.reserve(result.Cells.size());
    for (const auto& r_cell : result.Cells) {
        const auto& r_triangle_vertices = r_cell.Vertices;
        bounds.push_back({std::min({r_triangle_vertices[0][0], r_triangle_vertices[1][0], r_triangle_vertices[2][0]}),
                          std::min({r_triangle_vertices[0][1], r_triangle_vertices[1][1], r_triangle_vertices[2][1]}),
                          std::max({r_triangle_vertices[0][0], r_triangle_vertices[1][0], r_triangle_vertices[2][0]}),
                          std::max({r_triangle_vertices[0][1], r_triangle_vertices[1][1], r_triangle_vertices[2][1]})});
    }
    for (std::size_t i = 0; i<result.Cells.size(); ++i) {
        const auto& r_triangle_vertices = result.Cells[i].Vertices;
        Polygon polygon(r_triangle_vertices.begin(), r_triangle_vertices.end());
        KRATOS_ERROR_IF(std::abs(area(clip(polygon, boundary)))>Tolerance*scale)
            << "Gap region overlaps the fitted inner patch." << std::endl;
        for (std::size_t j = 0; j<i; ++j) {
            // Disjoint bounding boxes cannot have positive-area overlap.
            if (bounds[i][2] < bounds[j][0] || bounds[j][2] < bounds[i][0] ||
                bounds[i][3] < bounds[j][1] || bounds[j][3] < bounds[i][1]) continue;
            KRATOS_ERROR_IF(triangle_overlap(polygon, result.Cells[j].Vertices)>Tolerance*scale)
                << "Gap integration regions overlap." << std::endl;
        }
    }
    KRATOS_ERROR_IF(std::abs(gap_area+area(inner_points)-removed_area)>100*Tolerance*scale)
        << "Gap partition does not cover the removed background area: gap area = " << gap_area
        << ", inner area = " << area(inner_points) << ", removed area = " << removed_area
        << ", mismatch = " << gap_area+area(inner_points)-removed_area << "." << std::endl;
    KRATOS_ERROR_IF(std::abs(interface_length-perimeter)>100*Tolerance*std::max(1.0, perimeter))
        << "Gap partition does not cover the fitted interface." << std::endl;
    return result;
    KRATOS_CATCH("")
}

namespace {
void materialize_gap_geometry(const Utility::GapData& rData, ModelPart& rInterface,
                              const double Tolerance, ModelPart* pInner)
{
    KRATOS_TRY
    KRATOS_ERROR_IF(rInterface.NumberOfGeometries()!=0 || rInterface.NumberOfSubModelParts()!=0)
        << "Gap geometry destination must be empty; create it only once." << std::endl;
    auto& r_root = rInterface.GetRootModelPart();
    ModelPart::IndexType id = 1;
    for (const auto& r_geometry:r_root.Geometries()) id = std::max(id, r_geometry.Id()+1);
    const auto make_line = [&](ModelPart& rPart, const Utility::Segment& rSegment) {
        // Preserve actual mesh vertices on true/surrogate boundaries. Split
        // integration-only endpoints need no registered solution node.
        Node::Pointer p_start, p_end;
        const auto& r_parent = rSegment.pFirst->GetGeometry();
        for (std::size_t i = 0; i < r_parent.size(); ++i) {
            if (distance(r_parent[i].Coordinates(), rSegment.Start) <= Tolerance) p_start = r_parent.pGetPoint(i);
            if (distance(r_parent[i].Coordinates(), rSegment.End) <= Tolerance) p_end = r_parent.pGetPoint(i);
        }
        if (!p_start) p_start = Kratos::make_intrusive<Node>(0, rSegment.Start[0], rSegment.Start[1], 0.0);
        if (!p_end) p_end = Kratos::make_intrusive<Node>(0, rSegment.End[0], rSegment.End[1], 0.0);
        if (p_start->Id() != 0) rPart.AddNode(p_start);
        if (p_end->Id() != 0) rPart.AddNode(p_end);
        auto p_line = Kratos::make_shared<Line2D2<Node>>(p_start, p_end); p_line->SetId(id++);
        GlobalPointersVector<Element> parents;
        parents.push_back(GlobalPointer<Element>(rSegment.pFirst.get()));
        if (rSegment.pSecond) parents.push_back(GlobalPointer<Element>(rSegment.pSecond.get()));
        p_line->SetValue(NEIGHBOUR_ELEMENTS, parents); rPart.AddGeometry(p_line);
    };
    auto& r_truth = rInterface.CreateSubModelPart("TrueBoundary");
    if (pInner) {
        for (const auto& r_edge:contour(*pInner, Tolerance)) {
            make_line(r_truth, {r_edge.pParent, nullptr, r_edge.pStart->Coordinates(), r_edge.pEnd->Coordinates()});
        }
    } else {
        for (const auto& r_segment:rData.Physical) make_line(r_truth, r_segment);
    }
    auto& r_surrogate = rInterface.CreateSubModelPart("SurrogateBoundary");
    for (const auto& r_segment:rData.Surrogate) make_line(r_surrogate, r_segment);
    auto& r_volumes = rInterface.CreateSubModelPart("GapVolumes");
    auto& r_triangles = rInterface.CreateSubModelPart("GapTriangles");
    for (const auto& r_cell:rData.Cells) {
        r_cell.pQuadrature->SetId(id++);
        GlobalPointersVector<Element> parents;
        parents.push_back(GlobalPointer<Element>(r_cell.pParent.get()));
        r_cell.pQuadrature->SetValue(NEIGHBOUR_ELEMENTS, parents); r_volumes.AddGeometry(r_cell.pQuadrature);
        Geometry<Node>::PointsArrayType points;
        for (const auto& r_point:r_cell.Vertices) points.push_back(Kratos::make_intrusive<Node>(0, r_point[0], r_point[1], 0.0));
        auto p_triangle = Kratos::make_shared<Triangle2D3<Node>>(points); p_triangle->SetId(id++);
        r_triangles.AddGeometry(p_triangle);
    }
    auto& r_physical = rInterface.CreateSubModelPart("PhysicalInterfaces");
    for (const auto& r_segment:rData.Physical) make_line(r_physical, r_segment);
    auto& r_connectors = rInterface.CreateSubModelPart("GapInterfaces");
    for (const auto& r_segment:rData.Connectors) make_line(r_connectors, r_segment);
    KRATOS_CATCH("")
}
} // unnamed namespace

void FemMultipatchGapUtilities::CreateGeometry(ModelPart& rBackground, ModelPart& rInner,
                                               ModelPart& rInterface, const double Tolerance)
{
    KRATOS_ERROR_IF(rInterface.NumberOfGeometries()!=0 || rInterface.NumberOfSubModelParts()!=0)
        << "Gap geometry destination must be empty; create it only once." << std::endl;
    const auto data = Build(rBackground, rInner, Tolerance);
    materialize_gap_geometry(data, rInterface, Tolerance, &rInner);
}

void FemMultipatchGapUtilities::CreateHorizontalPhaseGeometry(
    ModelPart& rPhase, ModelPart& rInterface, const double InterfaceHeight,
    const bool NegativeSide, const double Tolerance)
{
    KRATOS_ERROR_IF(&rPhase.GetRootModelPart() != &rInterface.GetRootModelPart())
        << "Phase and interface must share a root ModelPart." << std::endl;
    const auto data = BuildHorizontalPhase(rPhase, InterfaceHeight, NegativeSide, Tolerance);
    materialize_gap_geometry(data, rInterface, Tolerance, nullptr);
}

void FemMultipatchGapUtilities::CreateLevelSetPhaseGeometry(
    ModelPart& rPhase, ModelPart& rInterface, const Matrix& rContour,
    const bool NegativeSide, const double Tolerance)
{
    KRATOS_ERROR_IF(&rPhase.GetRootModelPart() != &rInterface.GetRootModelPart())
        << "Phase and interface must share a root ModelPart." << std::endl;
    const auto data = BuildLevelSetPhase(rPhase, rContour, NegativeSide, Tolerance);
    materialize_gap_geometry(data, rInterface, Tolerance, nullptr);
}
} // namespace Kratos
