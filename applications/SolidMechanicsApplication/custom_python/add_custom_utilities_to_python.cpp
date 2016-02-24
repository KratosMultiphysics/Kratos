//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// External includes

// Project includes
#include "includes/node.h"
#include "includes/define.h"
#include "processes/process.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

//Application includes
#include "custom_python/add_custom_utilities_to_python.h"

//Utilities
#include "custom_utilities/line_search_calculation_utilities.hpp"
#include "custom_utilities/comparison_utilities.hpp"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"

#include "custom_utilities/energy_utilities.h"
#include "custom_utilities/isotropic_damage_utilities.hpp"

namespace Kratos
{

namespace Python
{

// //Boundary Skin Generator
// void GenerateSkin(ModelPart& model_part,char* ConditionName,unsigned int dimension,unsigned int preserve )
// {
//   GenerateBoundarySkin(model_part,KratosComponents<Condition>::Get(ConditionName),dimension,preserve);
// }



void  AddCustomUtilitiesToPython()
{

    using namespace boost::python;


    typedef UblasSpace<double, CompressedMatrix, Vector>    SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector>               LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    typedef Scheme< SparseSpaceType, LocalSpaceType >            SchemeType;
    typedef SchemeType::Pointer                           SchemePointerType;

    typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
    typedef BuilderAndSolverType::Pointer                   BuilderAndSolverPointerType;

    typedef Process                                         ProcessBaseType;

    class_<EnergyUtilities>("EnergyUtilities",init<>())
    .def("GetTotalKineticEnergy",&EnergyUtilities::GetTotalKineticEnergy)
    .def("CalculateNodalMass",&EnergyUtilities::CalculateNodalMass)
    .def("GetTotalStrainEnergy",&EnergyUtilities::GetTotalStrainEnergy)
    .def("GetGravitationalEnergy",&EnergyUtilities::GetGravitationalEnergy)
    .def("GetExternallyAppliedEnergy",&EnergyUtilities::GetExternallyAppliedEnergy)
    ;


}

}  // namespace Python.

} // Namespace Kratos

