// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____ 
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _ 
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Nelson Lafontaine
//                   Jordi Cotela Dalmau
//                   Riccardo Rossi
//                   Vicente Mataix Ferrándiz
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "meshing_application.h"
#include "custom_python/add_meshers_to_python.h"
#include "custom_python/add_processes_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_strategies_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;



BOOST_PYTHON_MODULE(KratosMeshingApplication)
{

    class_<KratosMeshingApplication,
           KratosMeshingApplication::Pointer,
           bases<KratosApplication>, boost::noncopyable >("KratosMeshingApplication")
           ;
    AddMeshersToPython();
    AddProcessesToPython();
    AddCustomUtilitiesToPython();
    AddCustomStrategiesToPython();

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(AVERAGE_NODAL_ERROR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ANISOTROPIC_RATIO)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(AUXILIAR_GRADIENT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(AUXILIAR_HESSIAN)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MMG_METRIC)
        
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(COUNTER)

    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(WEIGHT_FATHER_NODES) //used in the cutting planes app

}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined



