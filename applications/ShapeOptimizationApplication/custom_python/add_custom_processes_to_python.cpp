//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//   License:        BSD License
//   Kratos default license: kratos/license.txt
//
//   Project Name:        $ShapeOptimizationApplication     $
//   Last modified by:    $Author: martin.fusseder@tum.de   $
//   Date:                $Date:         August  2017       $
//   Revision:            $Revision:                0.0     $

// System includes

// fusseder TODO: this file was just copied, maybe rename it

// External includes
#include <boost/python.hpp>

// Project includes
#include "processes/process.h"
#include "custom_processes/output_primal_solution_process.h"
#include "custom_processes/input_primal_solution_process.h"
#include "custom_processes/replace_elements_and_conditions_for_adjoint_problem_process.h"
#include "includes/model_part.h"
#include "custom_python/add_custom_processes_to_python.h"

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

    class_<ReplaceElementsAndConditionsForAdjointProblemProcess , bases<Process>, boost::noncopyable >("ReplaceElementsAndConditionsForAdjointProblemProcess",
            init<ModelPart&, Parameters>())
    ;
}

} /* namespace Python */

} /* namespace Kratos */
