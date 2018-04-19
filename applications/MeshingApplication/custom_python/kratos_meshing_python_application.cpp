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

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"
#include "meshing_application.h"
#include "custom_python/add_meshers_to_python.h"
#include "custom_python/add_processes_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_strategies_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace pybind11;



PYBIND11_MODULE(KratosMeshingApplication,m)
{

    class_<KratosMeshingApplication,
           KratosMeshingApplication::Pointer,
           KratosApplication >(m, "KratosMeshingApplication")
           .def(init<>())
           ;
    AddMeshersToPython(m);
    AddProcessesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomStrategiesToPython(m);

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AVERAGE_NODAL_ERROR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ANISOTROPIC_RATIO)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AUXILIAR_GRADIENT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AUXILIAR_HESSIAN)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MMG_METRIC)
        
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(COUNTER)

    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(WEIGHT_FATHER_NODES) //used in the cutting planes app

}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined



