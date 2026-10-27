// License: BSD License (kratos/license.txt)
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include "utilities/multipatch_gap_utilities.h"
#include "utilities/gap_quadrature_utility.h"
#include "geometries/line_2d_2.h"
#include "geometries/triangle_2d_3.h"
#include "integration/line_gauss_legendre_integration_points.h"

namespace Kratos {
namespace {
using Utility = FemMultipatchGapUtilities;
using C = Utility::Coordinates;
using Polygon = std::vector<C>;
using Key = std::array<ModelPart::IndexType, 2>;
struct Edge { Node::Pointer pStart, pEnd; Element::Pointer pParent; };
struct EdgePair { std::vector<Edge> Sides; };
using Edges = std::map<Key, EdgePair>;

double cross(const C& a, const C& b, const C& c)
{
    return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);
}
double distance(const C& a, const C& b) { return std::hypot(a[0]-b[0], a[1]-b[1]); }
double area(const Polygon& p)
{
    double value = 0.0;
    for (std::size_t i = 1; i+1<p.size(); ++i) value += cross(p[0], p[i], p[i+1]);
    return value * 0.5;
}
Edges edges(ModelPart& rPart, const bool Cartesian)
{
    Edges result;
    for (auto& r_element : rPart.Elements()) {
        auto& g = r_element.GetGeometry();
        KRATOS_ERROR_IF(g.GetGeometryType() != GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4 || g.Area() <= 0.0)
            << "Gap-SBM requires positively oriented Quad4 parents." << std::endl;
        for (const auto& n : g) {
            KRATOS_ERROR_IF(std::abs(n.Z()) > 1e-12) << "Gap-SBM requires z=0." << std::endl;
        }
        if (Cartesian) {
            const double scale = std::sqrt(g.Area());
            KRATOS_ERROR_IF(std::abs(g[1].Y()-g[0].Y()) > 1e-10*scale ||
                            std::abs(g[2].X()-g[1].X()) > 1e-10*scale ||
                            std::abs(g[3].Y()-g[2].Y()) > 1e-10*scale ||
                            std::abs(g[0].X()-g[3].X()) > 1e-10*scale)
                << "Gap-SBM background must be Cartesian Quad4." << std::endl;
        }
        for (auto& edge : g.GenerateEdges()) {
            Key key{edge[0].Id(), edge[1].Id()};
            std::sort(key.begin(), key.end());
            auto& pair = result[key].Sides;
            pair.push_back({edge.pGetPoint(0), edge.pGetPoint(1), rPart.pGetElement(r_element.Id())});
            KRATOS_ERROR_IF(pair.size()>2) << "Non-manifold patch edge." << std::endl;
        }
    }
    return result;
}
std::vector<Edge> contour(ModelPart& rInner, const double Tol)
{
    const auto all = edges(rInner, false);
    std::map<ModelPart::IndexType, Edge> next;
    for (const auto& [key, pair] : all) {
        if (pair.Sides.size()==1) {
            const auto& e = pair.Sides[0];
            KRATOS_ERROR_IF_NOT(next.emplace(e.pStart->Id(), e).second)
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
    for (const auto& e: result) points.push_back(e.pStart->Coordinates());
    KRATOS_ERROR_IF(area(points)<=Tol*Tol) << "Inner contour must be counterclockwise." << std::endl;
    for (const auto& e:result) {
        const double length = distance(e.pStart->Coordinates(), e.pEnd->Coordinates());
        KRATOS_ERROR_IF(length<=Tol) << "Zero-length inner segment." << std::endl;
        for (const auto& p:points) {
            KRATOS_ERROR_IF(cross(e.pStart->Coordinates(), e.pEnd->Coordinates(), p)<-Tol*length)
                << "First Gap-SBM implementation requires a convex fitted contour." << std::endl;
        }
    }
    return result;
}
Polygon clip(Polygon p, const std::vector<Edge>& boundary)
{
    for (const auto& e:boundary) {
        Polygon out;
        if (p.empty()) break;
        const auto& a = e.pStart->Coordinates();
        const auto& b = e.pEnd->Coordinates();
        for (std::size_t i = 0; i<p.size(); ++i) {
            const C& s = p[i]; const C& t = p[(i+1)%p.size()];
            const double ds = cross(a, b, s), dt = cross(a, b, t);
            if (ds>=0.0) out.push_back(s);
            if ((ds<0.0 && dt>0.0)||(ds>0.0 && dt<0.0)) {
                const double w = ds/(ds-dt);
                C point = s+w*(t-s);
                out.push_back(point);
            }
        }
        p = std::move(out);
    }
    return p;
}
std::set<ModelPart::IndexType> cut_cells(ModelPart& rBackground,
                                        const std::vector<Edge>& boundary, const double Tol)
{
    std::set<ModelPart::IndexType> result;
    for (const auto& e:rBackground.Elements()) {
        Polygon polygon;
        for (const auto& n:e.GetGeometry()) polygon.push_back(n.Coordinates());
        const Polygon intersection = clip(polygon, boundary);
        if (intersection.size()<3 || area(intersection)<=0.0) continue;
        // A corner cut has area proportional to penetration squared: an
        // area cutoff Tol*h can miss cuts much deeper than Tol. Test an
        // interior point against both contours using distance tolerances.
        // This still excludes roundoff-sized slivers on fitted edges.
        C interior = ZeroVector(3);
        for (const auto& point:intersection) interior+=point;
        interior/=static_cast<double>(intersection.size());
        bool inside = true;
        for (std::size_t i = 0; i<polygon.size(); ++i) {
            const auto& a = polygon[i];
            const auto& b = polygon[(i+1)%polygon.size()];
            inside = inside && cross(a, b, interior)>Tol*distance(a, b);
        }
        for (const auto& edge:boundary) {
            const auto& a = edge.pStart->Coordinates();
            const auto& b = edge.pEnd->Coordinates();
            inside = inside && cross(a, b, interior)>Tol*distance(a, b);
        }
        if (inside) result.insert(e.Id());
    }
    return result;
}
double triangle_overlap(Polygon p, const std::array<C, 3>& triangle)
{
    for (std::size_t k = 0; k<3 && !p.empty(); ++k) {
        Polygon out;
        for (std::size_t i = 0; i<p.size(); ++i) {
            const C& a = p[i]; const C& b = p[(i+1)%p.size()];
            const double da = cross(triangle[k], triangle[(k+1)%3], a);
            const double db = cross(triangle[k], triangle[(k+1)%3], b);
            if (da>=0.0) out.push_back(a);
            if ((da<0.0 && db>0.0)||(da>0.0 && db<0.0)) {
                const C point = a+(da/(da-db))*(b-a); out.push_back(point);
            }
        }
        p = std::move(out);
    }
    return std::abs(area(p));
}
struct Projection { C Point; double S = 0.0; };
// Avoid the name "project": ADL also finds Boost uBLAS matrix projections
// through array_1d, causing invalid matrix-proxy template instantiations.
Projection project_to_contour(const C& p, const std::vector<Edge>& boundary, const std::vector<double>& lengths)
{
    double best = std::numeric_limits<double>::max();
    Projection result;
    for (std::size_t i = 0; i<boundary.size(); ++i) {
        const C a = boundary[i].pStart->Coordinates();
        const C d = boundary[i].pEnd->Coordinates()-a;
        const double length = lengths[i+1]-lengths[i];
        const double t = std::clamp(((p[0]-a[0])*d[0]+(p[1]-a[1])*d[1])/(length*length), 0.0, 1.0);
        const C q = a+t*d;
        const double dist = distance(p, q);
        if (dist<best) { best = dist; result = {q, lengths[i]+t*length}; }
    }
    if (result.S>=lengths.back()) result.S = 0.0;
    return result;
}
std::vector<std::array<C, 3>> triangulate(Polygon p, const double Tol)
{
    // Remove duplicate/collinear points before ordinary ear clipping.
    bool changed = true;
    while (changed && p.size()>2) {
        changed = false;
        for (std::size_t i = 0; i<p.size(); ++i) {
            const auto& prev = p[(i+p.size()-1)%p.size()];
            const auto& next = p[(i+1)%p.size()];
            if (distance(prev, p[i])<=Tol || std::abs(cross(prev, p[i], next))<=Tol*distance(prev, next)) {
                p.erase(p.begin()+i); changed = true; break;
            }
        }
    }
    std::vector<std::array<C, 3>> result;
    if (p.size()<3 || std::abs(area(p))<=Tol*Tol) return result;
    if (area(p)<0.0) std::reverse(p.begin(), p.end());
    while (p.size()>3) {
        bool found = false;
        for (std::size_t i = 0; i<p.size(); ++i) {
            const std::size_t a = (i+p.size()-1)%p.size(), b = i, c = (i+1)%p.size();
            if (cross(p[a], p[b], p[c])<=Tol*Tol) continue;
            bool contains = false;
            for (std::size_t j = 0; j<p.size(); ++j) {
                if (j==a||j==b||j==c) continue;
                if (cross(p[a], p[b], p[j])>=-Tol*Tol && cross(p[b], p[c], p[j])>=-Tol*Tol &&
                    cross(p[c], p[a], p[j])>=-Tol*Tol) { contains = true; break; }
            }
            if (contains) continue;
            result.push_back({p[a], p[b], p[c]}); p.erase(p.begin()+b); found = true; break;
        }
        KRATOS_ERROR_IF_NOT(found) << "Gap polygon is not simple; invalid projection topology." << std::endl;
    }
    if (cross(p[0], p[1], p[2])>Tol*Tol) result.push_back({p[0], p[1], p[2]});
    return result;
}
Geometry<Node>::Pointer quadrature(Element& rParent, const std::array<C, 3>& t)
{
    Geometry<Node>::IntegrationPointsArrayType points;
    const double determinant = cross(t[0], t[1], t[2]);
    // Duffy 3x3 rule integrates degree-four Cartesian Q1 products exactly.
    for (const auto& a:LineGaussLegendreIntegrationPoints3::IntegrationPoints()) {
        for (const auto& b:LineGaussLegendreIntegrationPoints3::IntegrationPoints()) {
            const double u = 0.5*(1+a.X()), v = 0.5*(1+b.X());
            const C p = t[0]+u*(t[1]-t[0])+(1-u)*v*(t[2]-t[0]);
            points.emplace_back(p[0], p[1], 0.0, 0.25*determinant*(1-u)*a.Weight()*b.Weight());
        }
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
    for (const auto& n:rBackground.Nodes()) {
        xmin = std::min(xmin, n.X()); ymin = std::min(ymin, n.Y());
        xmax = std::max(xmax, n.X()); ymax = std::max(ymax, n.Y());
    }
    FemMultipatchUtilities::OwnerMapType owners;
    std::set<ModelPart::IndexType> independent;
    for (const auto& n:rBackground.Nodes()) independent.insert(n.Id());
    for (std::size_t i = 0; i<rInnerPatches.size(); ++i) {
        KRATOS_ERROR_IF(&rInnerPatches[i]->GetRootModelPart()!=&rBackground.GetRootModelPart())
            << "Patches must share a root ModelPart." << std::endl;
        for (const auto& n:rInnerPatches[i]->Nodes()) {
            KRATOS_ERROR_IF(n.X()<=xmin+Tolerance || n.X()>=xmax-Tolerance ||
                            n.Y()<=ymin+Tolerance || n.Y()>=ymax-Tolerance)
                << "Inner patches must lie strictly inside the background." << std::endl;
            KRATOS_ERROR_IF_NOT(independent.insert(n.Id()).second) << "Patches must have independent nodes." << std::endl;
        }
        const auto boundary = contour(*rInnerPatches[i], Tolerance);
        for (const auto id:cut_cells(rBackground, boundary, Tolerance)) {
            KRATOS_ERROR_IF_NOT(owners.emplace(id, i+1).second)
                << "Inner patches or their conservative surrogate holes overlap." << std::endl;
        }
    }
    for (const auto& [key, pair]:all_edges) {
        if (pair.Sides.size()==1) {
            KRATOS_ERROR_IF(owners.count(pair.Sides[0].pParent->Id()))
                << "Surrogate holes must lie strictly inside the background." << std::endl;
        } else {
            const auto a = owners.find(pair.Sides[0].pParent->Id()), b = owners.find(pair.Sides[1].pParent->Id());
            KRATOS_ERROR_IF(a!=owners.end() && b!=owners.end() && a->second!=b->second)
                << "Surrogate holes must not touch." << std::endl;
        }
    }
    // Reject vertex contacts too, including non-manifold diagonal holes.
    std::map<ModelPart::IndexType, std::set<std::size_t>> node_owners;
    for (const auto& e:rBackground.Elements()) {
        const auto it = owners.find(e.Id());
        if (it!=owners.end()) for (const auto& n:e.GetGeometry()) node_owners[n.Id()].insert(it->second);
    }
    for (const auto& [id, set]:node_owners) KRATOS_ERROR_IF(set.size()>1) << "Surrogate holes touch at a vertex." << std::endl;
    // Commit flags only after classification/ownership validation succeeds.
    for (auto& e:rBackground.Elements()) e.Set(ACTIVE, owners.count(e.Id())==0);
    return owners;
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
    for (const auto& e:boundary) lengths.push_back(lengths.back()+distance(e.pStart->Coordinates(), e.pEnd->Coordinates()));
    const double perimeter = lengths.back();
    std::map<ModelPart::IndexType, Projection> projections;
    std::map<ModelPart::IndexType, std::vector<std::size_t>> incident;
    std::vector<C> centers;
    double removed_area = 0.0;
    for (const auto id:inactive) {
        const auto& e = rBackground.GetElement(id);
        KRATOS_ERROR_IF(e.Is(ACTIVE)) << "Run conservative activation before building the gap." << std::endl;
        removed_area+=e.GetGeometry().Area();
    }
    for (const auto& [key, pair]:all_edges) {
        if (pair.Sides.size()!=2) continue;
        const bool a = inactive.count(pair.Sides[0].pParent->Id()), b = inactive.count(pair.Sides[1].pParent->Id());
        if (a==b) continue;
        const auto& e = pair.Sides[a?1:0];
        KRATOS_ERROR_IF_NOT(e.pParent->Is(ACTIVE)) << "Surrogate holes overlap or touch." << std::endl;
        const std::size_t index = result.Surrogate.size();
        result.Surrogate.push_back({e.pParent, nullptr, e.pStart->Coordinates(), e.pEnd->Coordinates()});
        for (const auto& node:{e.pStart, e.pEnd}) {
            projections.emplace(node->Id(), project_to_contour(node->Coordinates(), boundary, lengths));
            incident[node->Id()].push_back(index);
        }
        const auto p0 = projections.at(e.pStart->Id()), p1 = projections.at(e.pEnd->Id());
        double delta = p0.S-p1.S;
        if (delta < -Tolerance) delta+=perimeter;
        if (std::abs(delta)<Tolerance) delta = 0.0;
        Polygon polygon{e.pStart->Coordinates(), e.pEnd->Coordinates(), p1.Point};
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
                result.Physical.push_back({e.pParent, boundary[k].pParent, point, previous});
            }
            previous = point; current = end;
        }
        if (distance(polygon.back(), p0.Point)>Tolerance) polygon.push_back(p0.Point);
        C center = ZeroVector(3);
        for (const auto& p:polygon) center+=p;
        center/=static_cast<double>(polygon.size()); centers.push_back(center);
        // Tolerance is a distance. A coincident edge can acquire a tiny
        // signed area from projection roundoff, so scale the area tolerance
        // by the local edge length rather than squaring the distance.
        const double area_tolerance = Tolerance * std::max(
            distance(e.pStart->Coordinates(), e.pEnd->Coordinates()),
            distance(p0.Point, p1.Point));
        KRATOS_ERROR_IF(area(polygon)>area_tolerance)
            << "Gap projection reverses orientation; invalid additive gap." << std::endl;
        for (const auto& triangle:triangulate(polygon, Tolerance)) {
            result.Cells.push_back({e.pParent, triangle, quadrature(*e.pParent, triangle)});
        }
    }
    KRATOS_ERROR_IF(result.Surrogate.empty()) << "No closed internal surrogate boundary found." << std::endl;
    for (const auto& [id, faces]:incident) {
        KRATOS_ERROR_IF(faces.size()!=2) << "Surrogate boundary must be a closed manifold." << std::endl;
        const auto& a = result.Surrogate[faces[0]]; const auto& b = result.Surrogate[faces[1]];
        C start = rBackground.GetNode(id).Coordinates(), end = projections.at(id).Point;
        if (a.pFirst==b.pFirst || distance(start, end)<=Tolerance) continue;
        const C difference = centers[faces[1]]-centers[faces[0]];
        if ((end[1]-start[1])*difference[0]-(end[0]-start[0])*difference[1]<0.0) std::swap(start, end);
        result.Connectors.push_back({a.pFirst, b.pFirst, start, end});
    }
    double gap_area = 0.0, interface_length = 0.0;
    for (const auto& cell:result.Cells) gap_area+=0.5*cross(cell.Vertices[0], cell.Vertices[1], cell.Vertices[2]);
    for (const auto& segment:result.Physical) interface_length+=distance(segment.Start, segment.End);
    Polygon inner_points;
    for (const auto& e:boundary) inner_points.push_back(e.pStart->Coordinates());
    const double scale = std::max(1.0, removed_area);
    for (std::size_t i = 0; i<result.Cells.size(); ++i) {
        const auto& t = result.Cells[i].Vertices;
        Polygon polygon(t.begin(), t.end());
        KRATOS_ERROR_IF(std::abs(area(clip(polygon, boundary)))>Tolerance*scale)
            << "Gap region overlaps the fitted inner patch." << std::endl;
        for (std::size_t j = 0; j<i; ++j) {
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

void FemMultipatchGapUtilities::CreateGeometry(ModelPart& rBackground, ModelPart& rInner,
                                               ModelPart& rInterface, const double Tolerance)
{
    KRATOS_TRY
    KRATOS_ERROR_IF(rInterface.NumberOfGeometries()!=0 || rInterface.NumberOfSubModelParts()!=0)
        << "Gap geometry destination must be empty; create it only once." << std::endl;
    const auto data = Build(rBackground, rInner, Tolerance);
    auto& root = rInterface.GetRootModelPart();
    ModelPart::IndexType id = 1;
    for (const auto& g:root.Geometries()) id = std::max(id, g.Id()+1);
    const auto make_line = [&](ModelPart& part, const Segment& segment) {
        // Preserve actual mesh vertices on true/surrogate boundaries. Split
        // integration-only endpoints need no registered solution node.
        Node::Pointer a, b;
        const auto& parent = segment.pFirst->GetGeometry();
        for (std::size_t i = 0; i < parent.size(); ++i) {
            if (distance(parent[i].Coordinates(), segment.Start) <= Tolerance) a = parent.pGetPoint(i);
            if (distance(parent[i].Coordinates(), segment.End) <= Tolerance) b = parent.pGetPoint(i);
        }
        if (!a) a = Kratos::make_intrusive<Node>(0, segment.Start[0], segment.Start[1], 0.0);
        if (!b) b = Kratos::make_intrusive<Node>(0, segment.End[0], segment.End[1], 0.0);
        if (a->Id() != 0) part.AddNode(a);
        if (b->Id() != 0) part.AddNode(b);
        auto line = Kratos::make_shared<Line2D2<Node>>(a, b); line->SetId(id++);
        GlobalPointersVector<Element> parents;
        parents.push_back(GlobalPointer<Element>(segment.pFirst.get()));
        if (segment.pSecond) parents.push_back(GlobalPointer<Element>(segment.pSecond.get()));
        line->SetValue(NEIGHBOUR_ELEMENTS, parents); part.AddGeometry(line);
    };
    auto& truth = rInterface.CreateSubModelPart("TrueBoundary");
    for (const auto& edge:contour(rInner, Tolerance)) {
        make_line(truth, {edge.pParent, nullptr, edge.pStart->Coordinates(), edge.pEnd->Coordinates()});
    }
    auto& surrogate = rInterface.CreateSubModelPart("SurrogateBoundary");
    for (const auto& segment:data.Surrogate) make_line(surrogate, segment);
    auto& volumes = rInterface.CreateSubModelPart("GapVolumes");
    auto& triangles = rInterface.CreateSubModelPart("GapTriangles");
    for (const auto& cell:data.Cells) {
        cell.pQuadrature->SetId(id++);
        GlobalPointersVector<Element> parents;
        parents.push_back(GlobalPointer<Element>(cell.pParent.get()));
        cell.pQuadrature->SetValue(NEIGHBOUR_ELEMENTS, parents); volumes.AddGeometry(cell.pQuadrature);
        Geometry<Node>::PointsArrayType points;
        for (const auto& p:cell.Vertices) points.push_back(Kratos::make_intrusive<Node>(0, p[0], p[1], 0.0));
        auto triangle = Kratos::make_shared<Triangle2D3<Node>>(points); triangle->SetId(id++);
        triangles.AddGeometry(triangle);
    }
    auto& physical = rInterface.CreateSubModelPart("PhysicalInterfaces");
    for (const auto& segment:data.Physical) make_line(physical, segment);
    auto& connectors = rInterface.CreateSubModelPart("GapInterfaces");
    for (const auto& segment:data.Connectors) make_line(connectors, segment);
    KRATOS_CATCH("")
}
} // namespace Kratos
