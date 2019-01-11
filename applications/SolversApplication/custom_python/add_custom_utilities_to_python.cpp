//
//   Project Name:        KratosSolversApplication    $
//   Created by:          $Author:        JMCarbonell $
//   Last modified by:    $Co-Author:                 $
//   Date:                $Date:         January 2019 $
//   Revision:            $Revision:              0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"

// Utilities
#include "custom_utilities/eigenvector_to_solution_step_variable_transfer_utility.hpp"


namespace Kratos
{

namespace Python
{

inline
void TransferEigenvector1(EigenvectorToSolutionStepVariableTransferUtility& rThisUtil,
                          ModelPart& rModelPart,
                          int iEigenMode)
{
  rThisUtil.Transfer(rModelPart,iEigenMode);
}

inline
void TransferEigenvector2(EigenvectorToSolutionStepVariableTransferUtility& rThisUtil,
                          ModelPart& rModelPart,
                          int iEigenMode,
                          int step)
{
  rThisUtil.Transfer(rModelPart,iEigenMode,step);
}

void  AddCustomUtilitiesToPython(pybind11::module& m)
{

  namespace py = pybind11;

  py::class_<EigenvectorToSolutionStepVariableTransferUtility>(m,"EigenvectorToSolutionStepVariableTransferUtility")
      .def("Transfer",TransferEigenvector1)
      .def("Transfer",TransferEigenvector2)
      ;

}

}  // namespace Python.

} // Namespace Kratos
