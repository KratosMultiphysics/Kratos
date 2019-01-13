//
//   Project Name:        KratosUmatApplication      $
//   Created by:          $Author:      JMCarbonell  $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:      September 2017 $
//   Revision:            $Revision:             0.0 $
//
//

// System includes
#include <pybind11/stl.h>

// External includes

// Project includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

//models
#include "custom_models/hypoplastic_umat_small_strain_model.hpp"
#include "custom_models/von_mises_umat_small_strain_model.hpp"
#include "custom_models/von_mises_umat_large_strain_model.hpp"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;


void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
{

  // models
  py::class_<VonMisesSmallStrainUmatModel, typename VonMisesSmallStrainUmatModel::Pointer, ConstitutiveModel>
      (m, "VonMisesSmallStrainUmatModel")
      .def( py::init<>() )
      ;
  py::class_<VonMisesLargeStrainUmatModel, typename VonMisesLargeStrainUmatModel::Pointer, ConstitutiveModel>
      (m, "VonMisesLargeStrainUmatModel")
      .def( py::init<>() )
      ;
  py::class_<HypoplasticSmallStrainUmatModel, typename HypoplasticSmallStrainUmatModel::Pointer, ConstitutiveModel>
      (m, "HypoplasticSmallStrainUmatModel")
      .def( py::init<>() )
      ;
}

}  // namespace Python.
} // Namespace Kratos
