//
//   Project Name:        Kratos
//   Last Modified by:    $Author: clabra $
//   Date:                $Date: 2007-03-29 19:37:47 $
//   Revision:            $Revision: 1.2 $
//
//

#if !defined(KRATOS_POINT_POINT_SEARCH_H_INCLUDED)
#define  KRATOS_POINT_POINT_SEARCH_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// include kratos definitions
#include "includes/define.h"

// Project includes
#include "utilities/openmp_utils.h"

// Configures
#include "spatial_containers/spatial_search.h"
#include "point_configure.h"
// Search
#include "spatial_containers/bins_dynamic_objects.h"
#include "spatial_containers/bins_dynamic.h"

// External includes

/* Timer defines */
#include "utilities/timer.h"
#ifdef CUSTOMTIMER
#define KRATOS_TIMER_START(t) Timer::Start(t);
#define KRATOS_TIMER_STOP(t) Timer::Stop(t);
#else
#define KRATOS_TIMER_START(t)
#define KRATOS_TIMER_STOP(t)
#endif

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/

class PointPointSearch: public SpatialSearch
{
public:
///@name Type Definitions
///@{

/// Pointer definition of PointPointSearch
KRATOS_CLASS_POINTER_DEFINITION(PointPointSearch);

typedef PointType*                                  PointPointerType;
typedef std::vector<PointPointerType>*              PointVector;
typedef std::vector<PointPointerType>::iterator     PointIterator;

typedef double*                                     DistanceVector;
typedef double*                                     DistanceIterator;

//Configure Types
typedef PointConfigure<3>                           PointConfigureType;
//Bin Types
typedef BinsObjectDynamic<PointConfigureType>       PointBinsType;
typedef PointerVectorSet<Point, IndexedObject>      PointSetType;


///@}
///@name Life Cycle
///@{

/// Default constructor.
PointPointSearch(){}

/// Destructor.
~PointPointSearch(){}

void SearchPointsImplementation(
        NodesContainerType const& r_nodes,
        NodesContainerType const& r_nodes_to_find,
        RadiusArrayType const& radius,
        VectorResultNodesContainerType& r_results,
        VectorDistanceType& r_results_distances)
{
    KRATOS_TRY

    int max_n_of_neigh_nodes = r_nodes_to_find.size();

    NodesContainerType::ContainerType& nodes         = const_cast <NodesContainerType::ContainerType&> (r_nodes.GetContainer());
    NodesContainerType::ContainerType& nodes_to_find = const_cast <NodesContainerType::ContainerType&> (r_nodes_to_find.GetContainer());

    PointSetType::ContainerType nodes_temp;
    PointSetType::ContainerType nodes_to_find_temp;

    std::map<Point::Pointer, Node<3>::Pointer> map_point_to_node;

    nodes_temp.reserve(nodes.size());

    for (NodesContainerType::ContainerType::iterator it = nodes.begin(); it != nodes.end(); ++it){
        nodes_temp.push_back(std::make_shared<Point>((*it)->Coordinates()));
    }

    nodes_to_find_temp.reserve(nodes_to_find.size());

    for (NodesContainerType::ContainerType::iterator it = nodes_to_find.begin(); it != nodes_to_find.end(); ++it){
        nodes_to_find_temp.push_back(std::make_shared<Point>((*it)->Coordinates()));
    }

    PointBinsType bins(nodes_to_find_temp.begin(), nodes_to_find_temp.end());

    #pragma omp parallel
    {
        PointSetType::ContainerType local_results(max_n_of_neigh_nodes);
        DistanceType                local_results_distances(max_n_of_neigh_nodes);
        std::size_t                 n_of_results = 0;

        #pragma omp for
        for (int i = 0; i < static_cast<int>(nodes.size()); ++i){
            PointSetType::ContainerType::iterator i_results_begin = local_results.begin();
            DistanceType::iterator i_distances_results_begin      = local_results_distances.begin();

            n_of_results = bins.SearchObjectsInRadiusExclusive(nodes_temp[i], radius[i], i_results_begin, i_distances_results_begin, max_n_of_neigh_nodes);

            r_results[i].reserve(n_of_results);

            for (PointSetType::ContainerType::iterator it = local_results.begin(); it != local_results.begin() + n_of_results; ++it){
                r_results[i].push_back(map_point_to_node[ *(it.base() )]);
            }

            r_results_distances[i].insert(r_results_distances[i].begin(), local_results_distances.begin(), local_results_distances.begin() + n_of_results);
        }

            KRATOS_WATCH(r_results_distances.size())
    }

    KRATOS_CATCH("")
}

/// Turn back information as a string.
virtual std::string Info() const override
{
std::stringstream buffer;
buffer << "PointPointSearch" ;

return buffer.str();
}

/// Print information about this object.
virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "PointPointSearch";}

/// Print object's data.
virtual void PrintData(std::ostream& rOStream) const override {}


///@}
///@name Friends
///@{


///@}

protected:
///@name Protected static Member Variables
///@{


///@}
///@name Protected member Variables
///@{



///@}
///@name Protected Operators
///@{


///@}
///@name Protected Operations
///@{


///@}
///@name Protected  Access
///@{


///@}
///@name Protected Inquiry
///@{


///@}
///@name Protected LifeCycle
///@{


///@}

private:
///@name Static Member Variables
///@{


///@}
///@name Member Variables
///@{


///@}
///@name Private Operators
///@{


///@}
///@name Private Operations
///@{


///@}
///@name Private  Access
///@{


///@}
///@name Private Inquiry
///@{


///@}
///@name Un accessible methods
///@{

/// Assignment operator.
PointPointSearch& operator=(PointPointSearch const& rOther)
{
    return *this;
}

/// Copy constructor.
PointPointSearch(PointPointSearch const& rOther)
{
    *this = rOther;
}


}; // Class PointPointSearch


}  // namespace Kratos.

#endif // KRATOS_POINT_POINT_SEARCH_H_INCLUDED  defined


