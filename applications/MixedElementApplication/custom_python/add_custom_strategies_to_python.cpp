//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes


// External includes
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_strategies/convergencecriterias/mixed_element_criteria.h"

// Project includes
#include "includes/define_python.h"

#include "spaces/ublas_space.h"



namespace Kratos
{

namespace Python
{
using namespace pybind11;

void AddCustomStrategiesToPython(pybind11::module& m)
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaBaseType;

    class_< MixedElementConvergenceCriteria< SparseSpaceType,  LocalSpaceType > ,  ConvergenceCriteriaBaseType  >(m,"MixedElementConvergenceCriteria")
            .def(init<double, double >())
            ;

}

} // namespace Python.

} // Namespace Kratos

