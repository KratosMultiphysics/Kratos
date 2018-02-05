//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "mapping_application.h"
#include "mapping_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_mappers_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;



BOOST_PYTHON_MODULE(KratosMappingApplication)
{

    class_<KratosMappingApplication,
           KratosMappingApplication::Pointer,
           bases<KratosApplication>, boost::noncopyable >("KratosMappingApplication")
           ;

    AddCustomStrategiesToPython();
    AddCustomUtilitiesToPython();
    AddCustomProcessesToPython();
    AddCustomMappersToPython();

    //registering variables in python

//	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_AREA);


}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
