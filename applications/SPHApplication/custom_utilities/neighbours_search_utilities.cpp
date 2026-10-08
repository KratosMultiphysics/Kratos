#include "custom_utilities/neighbours_search_utilities.h"

// System includes
#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>
#include <vector>

// Project includes
#include "sph_application_variables.h"
#include "spatial_containers/bins_static.h"
#include "spatial_containers/point_object.h"

namespace Kratos
{

namespace
{

/**
 * @brief Wrap nodes/elements as PointObject instances so they can be used by the spatial bins
 * @param TContainerType The type of the container holding the objects.
 */
template<class TContainerType>
auto CreateSearchPoints(
    const TContainerType& rContainer,
    const std::size_t DomainSize)
{
    using ObjectType = typename TContainerType::value_type;
    using SearchPointType = PointObject<ObjectType>;
    using SearchPointPointerType = typename SearchPointType::Pointer;

    std::vector<SearchPointPointerType> points;
    points.reserve(rContainer.size());

    for (auto it_object = rContainer.begin(); it_object != rContainer.end(); ++it_object) {
        auto p_point = SearchPointPointerType(new SearchPointType(*(it_object.base())));

        if constexpr (std::is_same_v<ObjectType, Element>) {
            noalias(p_point->Coordinates()) = it_object->GetGeometry()[0].Coordinates();
        }

        if (DomainSize == 2) (*p_point)[2] = 0.0;

        points.push_back(p_point);
    }

    return points;
}

/**
 * @brief Compute the squared distance between two points in 3D space.
 */
template<class TPointType>
double SquaredDistance(
    const TPointType& rPointA,
    const TPointType& rPointB)
{
    const double dx = rPointA[0] - rPointB[0];
    const double dy = rPointA[1] - rPointB[1];
    const double dz = rPointA[2] - rPointB[2];
    return dx * dx + dy * dy + dz * dz;
}

/**
 * @brief Search for PointObject instances within a given radius of a search point, including the search point itself.
 * @tparam TSearchPointPointerType The type of the pointer to the search point.
 * @tparam TBinsType The type of the spatial bins used for searching.
 * @param pSearchPoint Pointer to the central search point.
 * @param Radius The radius within which to search for neighbours.
 */
template<class TSearchPointPointerType, class TBinsType>
std::vector<TSearchPointPointerType> SearchInRadiusInclusive(
    TSearchPointPointerType pSearchPoint,
    const double Radius,
    const std::size_t NumberOfPoints,
    TBinsType& rBins)
{
    // BinsStatic uses a strict comparison with the squared search radius.
    // std::nexafter takes the next representable value after Radius towards infinity
    double bins_radius = std::nextafter(Radius, std::numeric_limits<double>::infinity());

    // Allocates memory for a maximum of 128 neighbours, which should be enough for all cases. 
    std::size_t upper_bound = std::min<std::size_t>(128, NumberOfPoints);
    std::vector<TSearchPointPointerType> candidates(upper_bound);

    std::size_t number_of_candidates;
    number_of_candidates = rBins.SearchInRadius(*pSearchPoint, bins_radius, candidates.begin(), upper_bound);

    candidates.resize(number_of_candidates);

    return candidates;
}

} // namespace

double NeighboursSearchUtilities::ComputeSmoothingLength(const ModelPart& rThisModelPart, double Coeff)
{
    return Coeff * ComputeInterparticleMinDist(rThisModelPart);
}

double NeighboursSearchUtilities::ComputeInterparticleMinDist(const ModelPart& rThisModelPart)
{
    const auto& rnodes = rThisModelPart.Nodes();
    const std::size_t domain_size = rThisModelPart.GetProcessInfo()[DOMAIN_SIZE];

    double min_dist = std::numeric_limits<double>::max();

    auto points = CreateSearchPoints(rnodes, domain_size);

    using SearchPointType = PointObject<Node>;
    using SearchPointVectorType = std::vector<SearchPointType::Pointer>;
    using StaticBinsType = Bins<3, SearchPointType, SearchPointVectorType>;

    StaticBinsType bins(points.begin(), points.end());

    for (auto& p_point : points) {
        auto p_nearest_point = bins.SearchNearestPointInner(p_point); // returns the nearest point to p_point
        min_dist = std::min(min_dist, SquaredDistance(*p_point, *p_nearest_point));
    }

    return std::sqrt(domain_size * min_dist);
}

void NeighboursSearchUtilities::SearchNeighbours(
    ModelPart& rThisModelPart,
    double Radius)
{
    auto& r_elem = rThisModelPart.Elements();
    const std::size_t number_of_elements = r_elem.size();

    const std::size_t domain_size = rThisModelPart.GetProcessInfo()[DOMAIN_SIZE];
    auto points = CreateSearchPoints(r_elem, domain_size);

    using SearchPointType = PointObject<Element>;
    using SearchPointPointerType = SearchPointType::Pointer;
    using SearchPointVectorType = std::vector<SearchPointPointerType>;
    using StaticBinsType = Bins<3, SearchPointType, SearchPointVectorType>;

    StaticBinsType bins(points.begin(), points.end());

    for (auto& p_point : points) {
        auto candidate_points = SearchInRadiusInclusive(p_point, Radius, number_of_elements, bins);

        std::vector<Element::Pointer> neighbours;
        neighbours.reserve(candidate_points.size() + 1);

        const auto p_element = p_point->pGetObject();

        for (const auto& p_candidate_point : candidate_points) {
            neighbours.push_back(p_candidate_point->pGetObject());
        }

        auto& r_element = *p_element;
        r_element.SetValue(NEIGHBOURS, neighbours);
    }
}

}