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
typedef PointerVectorSet<Point<3>, IndexedObject>   PointSetType;


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
        const RadiusArrayType & radius,
        VectorResultNodesContainerType& r_results,
        VectorDistanceType& r_results_distances)
{
    KRATOS_TRY

    int max_number_of_nodes = r_nodes.size();

    NodesContainerType::ContainerType& nodes      = const_cast <NodesContainerType::ContainerType&> (r_nodes.GetContainer());
    NodesContainerType::ContainerType& nodes_bins = const_cast <NodesContainerType::ContainerType&> (r_nodes_to_find.GetContainer());

    PointSetType::ContainerType nodes_to_geometrical_object_pointer_vector;
    PointSetType::ContainerType bins_of_nodes_to_geometrical_object_pointer_vector;

    nodes_to_geometrical_object_pointer_vector.reserve(nodes.size());
    bins_of_nodes_to_geometrical_object_pointer_vector.reserve(nodes_bins.size());

    for (NodesContainerType::ContainerType::iterator it = nodes.begin(); it != nodes.end(); ++it){
        nodes_to_geometrical_object_pointer_vector.push_back(*it);
    }

    for (NodesContainerType::ContainerType::iterator it = nodes_bins.begin(); it != nodes_bins.end(); ++it){
        bins_of_nodes_to_geometrical_object_pointer_vector.push_back(*it);
    }

    PointBinsType bins(bins_of_nodes_to_geometrical_object_pointer_vector.begin(), bins_of_nodes_to_geometrical_object_pointer_vector.end());

    #pragma omp parallel
    {
        PointSetType::ContainerType local_results(max_number_of_nodes);
        DistanceType                         local_results_distances(max_number_of_nodes);
        std::size_t                          number_of_results = 0;

        #pragma omp for
        for (int i = 0; i < static_cast<int>(nodes.size()); ++i){
            PointSetType::ContainerType::iterator i_results = local_results.begin();
            DistanceType::iterator i_distances_results      = local_results_distances.begin();

            number_of_results = bins.SearchObjectsInRadiusExclusive(nodes_to_geometrical_object_pointer_vector[i], radius[i], i_results, i_distances_results, max_number_of_nodes);
            r_results[i].reserve(number_of_results);

            for (PointSetType::ContainerType::iterator it = local_results.begin(); it != local_results.begin() + number_of_results; ++it){
                Node<3>::Pointer p_node = boost::dynamic_pointer_cast<Node<3> >(*it);
                r_results[i].push_back(p_node);
                r_results_distances[i].insert(r_results_distances[i].begin(), local_results_distances.begin(), local_results_distances.begin() + number_of_results);
            }
        }
    }

    KRATOS_CATCH("")
}

/// Turn back information as a string.
virtual std::string Info() const
{
std::stringstream buffer;
buffer << "PointPointSearch" ;

return buffer.str();
}

/// Print information about this object.
virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "PointPointSearch";}

/// Print object's data.
virtual void PrintData(std::ostream& rOStream) const {}


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


