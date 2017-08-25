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

#include "custom_utilities/eigenvector_to_solution_step_variable_transfer_utility.hpp"


namespace Kratos
{

  namespace Python
  {

  inline
  void TransferEigenvector1(
        EigenvectorToSolutionStepVariableTransferUtility& rThisUtil,
        ModelPart& rModelPart,
        int iEigenMode)
  {
    rThisUtil.Transfer(rModelPart,iEigenMode);
  }

  inline
  void TransferEigenvector2(
        EigenvectorToSolutionStepVariableTransferUtility& rThisUtil,
        ModelPart& rModelPart,
        int iEigenMode,
        int step)
  {
    rThisUtil.Transfer(rModelPart,iEigenMode,step);
  }

  void  AddCustomUtilitiesToPython()
  {

      using namespace boost::python;

      class_<EnergyUtilities>("EnergyUtilities",init<>())
	.def("GetTotalKinematicEnergy",&EnergyUtilities::GetTotalKinematicEnergy)
	.def("CalculateNodalMass",&EnergyUtilities::CalculateNodalMass)
	.def("GetTotalStrainEnergy",&EnergyUtilities::GetTotalStrainEnergy)
	.def("GetGravitationalEnergy",&EnergyUtilities::GetGravitationalEnergy)
	.def("GetExternallyAppliedEnergy",&EnergyUtilities::GetExternallyAppliedEnergy)
	;

      class_<EigenvectorToSolutionStepVariableTransferUtility>("EigenvectorToSolutionStepVariableTransferUtility")
	.def("Transfer",TransferEigenvector1)
	.def("Transfer",TransferEigenvector2)
	;

  }

  }  // namespace Python.

} // Namespace Kratos

