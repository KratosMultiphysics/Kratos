//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//   License:        BSD License
//   Kratos default license: kratos/license.txt
//
//   Project Name:        $AdjointFluidApplication        $
//   Last modified by:    $Author: michael.andre@tum.de   $
//   Date:                $Date:         November  2016   $
//   Revision:            $Revision:                0.0   $

// System includes

// External includes
#include "pybind11/pybind11.h"

// Project includes
#include "processes/process.h"
#include "custom_processes/output_primal_solution_process.h"
#include "custom_processes/input_primal_solution_process.h"
#include "includes/model_part.h"

namespace Kratos
{

namespace Python
{

void AddCustomProcessesToPython(pybind11::module& m)
{
    using namespace pybind11;

    class_< OutputPrimalSolutionProcess >(m,"OutputPrimalSolutionProcess")
    .def(init<ModelPart&, Parameters&>())
    ;

    class_< InputPrimalSolutionProcess >(m,"InputPrimalSolutionProcess")
    .def(init<ModelPart&, Parameters&>())
    ;
}

} /* namespace Python */

} /* namespace Kratos */
