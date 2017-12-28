//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//

#if !defined(KRATOS_ADD_COUPLING_STRATEGIES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_COUPLING_STRATEGIES_TO_PYTHON_H_INCLUDED

// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"


namespace Kratos
{

    namespace Python
    {

      void  AddCustomCouplingStrategiesToPython()
      {
		  using namespace boost::python;
	  }

    }  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_STRATEGIES_PYTHON_H_INCLUDED  defined
