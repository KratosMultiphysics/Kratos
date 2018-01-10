//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//



// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define_python.h"
#include "python/add_search_strategies_to_python.h"
#include "spatial_containers/spatial_search.h"

namespace Kratos
{

namespace Python
{
  
void  AddSearchStrategiesToPython()
{
    using namespace boost::python;
  
    class_<SpatialSearch, boost::noncopyable >
             ("SpatialSearch", init< >())
             ;
}

}  // namespace Python.

} // Namespace Kratos

