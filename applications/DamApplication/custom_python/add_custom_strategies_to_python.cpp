//
//   Project Name:
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: $
//

// External includes
#include "spaces/ublas_space.h"

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "includes/kratos_parameters.h"

//linear solvers
#include "linear_solvers/linear_solver.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"

//builders and solvers

//schemes
#include "custom_strategies/schemes/incrementalupdate_static_smoothing_scheme.hpp"
#include "custom_strategies/schemes/incrementalupdate_static_damped_smoothing_scheme.hpp"
#include "custom_strategies/schemes/bossak_displacement_smoothing_scheme.hpp"
#include "custom_strategies/schemes/dam_UP_scheme.hpp"
#include "custom_strategies/schemes/dam_P_scheme.hpp"


namespace Kratos
{

namespace Python
{

namespace py = pybind11;

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

    //custom scheme types
    typedef IncrementalUpdateStaticSmoothingScheme< SparseSpaceType, LocalSpaceType >  IncrementalUpdateStaticSmoothingSchemeType;
    typedef IncrementalUpdateStaticDampedSmoothingScheme< SparseSpaceType, LocalSpaceType >  IncrementalUpdateStaticDampedSmoothingSchemeType;
    typedef BossakDisplacementSmoothingScheme< SparseSpaceType, LocalSpaceType >  BossakDisplacementSmoothingSchemeType;
    typedef DamUPScheme< SparseSpaceType, LocalSpaceType >  DamUPSchemeType;
    typedef DamPScheme< SparseSpaceType, LocalSpaceType >  DamPSchemeType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// Schemes
    py::class_< IncrementalUpdateStaticSmoothingSchemeType, typename IncrementalUpdateStaticSmoothingSchemeType::Pointer, BaseSchemeType >
    (m, "IncrementalUpdateStaticSmoothingScheme")
    .def(py::init< >());

    py::class_< IncrementalUpdateStaticDampedSmoothingSchemeType, typename IncrementalUpdateStaticDampedSmoothingSchemeType::Pointer, BaseSchemeType >
    (m, "IncrementalUpdateStaticDampedSmoothingScheme")
    .def(py::init< double, double >());

    py::class_< BossakDisplacementSmoothingSchemeType, typename BossakDisplacementSmoothingSchemeType::Pointer, BaseSchemeType >
    (m, "BossakDisplacementSmoothingScheme")
    .def(py::init< double, double, double >());

	py::class_< DamUPSchemeType, typename DamUPSchemeType::Pointer, BaseSchemeType >
    (m, "DamUPScheme")
    .def(py::init< double, double, double, double >());

    py::class_< DamPSchemeType, typename DamPSchemeType::Pointer, BaseSchemeType >
    (m, "DamPScheme")
    .def(py::init< double, double >());
}

}  // namespace Python.
} // Namespace Kratos

