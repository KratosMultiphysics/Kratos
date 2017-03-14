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
#include <boost/python.hpp>

// Project includes
#include "processes/process.h"
#include "custom_processes/output_primal_solution_process.h"
#include "custom_processes/input_primal_solution_process.h"
#include "includes/model_part.h"

namespace Kratos
{

namespace Python
{

void AddCustomProcessesToPython()
{
    using namespace boost::python;

    class_< OutputPrimalSolutionProcess, bases<Process> >
    ("OutputPrimalSolutionProcess", init<ModelPart&, Parameters&>())
    ;

    class_< InputPrimalSolutionProcess, bases<Process> >
    ("InputPrimalSolutionProcess", init<ModelPart&, Parameters&>())
    ;
}

} /* namespace Python */

} /* namespace Kratos */
