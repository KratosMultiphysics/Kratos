//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

// System includes
#include <pybind11/pybind11.h>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/constitutive_law.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/properties.h"

#include "python/variable_indexing_python.h"
#include "python/add_mesh_to_python.h"


//Application includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

//constitutive laws
#include "custom_constitutive/zarate_law.hpp"
#include "custom_constitutive/fem_dem_elastic_law.hpp"

namespace Kratos
{
namespace Python
{
	void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
	{
		py::class_<ZarateLaw, typename ZarateLaw::Pointer, ConstitutiveLaw >
			(m, "ZarateLaw").def(py::init<>() )
			;

		py::class_<FemDemElasticLaw, typename FemDemElasticLaw::Pointer, ConstitutiveLaw >
			(m, "FemDemElasticLaw").def(py::init<>() )
			;
	}
}  // namespace Python.
}  // namespace Kratos.