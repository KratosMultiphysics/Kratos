//
//   Project Name:        Kratos
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2008-12-09 20:20:55 $
//   Revision:            $Revision: 1.5 $
//
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>

//Trilinos includes
#include "mpi.h"

// Project includes
#include "custom_utilities/zoltan_partition_utility.h"
#include "custom_python/add_zoltan_processes_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

void AddZoltanProcessesToPython()
{
    class_<ZoltanPartitionUtility, boost::noncopyable >
    ("ZoltanPartitionUtility",
     init< >() )
    .def("CalculatePartition", &ZoltanPartitionUtility::CalculatePartition )
    ;
}


} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
