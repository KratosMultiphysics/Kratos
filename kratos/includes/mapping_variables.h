//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_MAPPING_VARIABLES_H_INCLUDED)
#define KRATOS_MAPPING_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "containers/variable.h"
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"
#include "includes/kratos_components.h"
#include "includes/mortar_classes.h"

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_API

namespace Kratos
{
   KRATOS_DEFINE_VARIABLE( boost::shared_ptr<ConditionMap>, MAPPING_PAIRS ) // An unordened map of which contains the structure

} // namespace Kratos

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT

#endif // KRATOS_MAPPING_VARIABLES_H_INCLUDED defined
