//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_DEFINE_DEPRECATED_H_INCLUDED )
#define  KRATOS_DEFINE_DEPRECATED_H_INCLUDED

/* System includes */
#include <stdexcept>
#include <sstream>

/* External includes */

/* Project includes */
#include "includes/kratos_version.h"

#if __cplusplus >= 201402L
#define KRATOS_DEPRECATED [[deprecated]]
#define KRATOS_DEPRECATED_MESSAGE(deprecated_message) [[deprecated(deprecated_message)]]
#elif __GNUC__
#define KRATOS_DEPRECATED __attribute__((deprecated))
#define KRATOS_DEPRECATED_MESSAGE(deprecated_message) KRATOS_DEPRECATED
#elif defined(_MSC_VER)
#define KRATOS_DEPRECATED __declspec(deprecated)
#define KRATOS_DEPRECATED_MESSAGE(deprecated_message) KRATOS_DEPRECATED
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define KRATOS_DEPRECATED
#define KRATOS_DEPRECATED_MESSAGE(deprecated_message)
#endif

#define KRATOS_VERSION_AS_INT (KRATOS_MAJOR_VERSION*100+KRATOS_MINOR_VERSION*10+KRATOS_PATCH_VERSION)
#define KRATOS_TO_BE_REMOVED_IN_VERSION(major_version, minor_version, patch_version) static_assert((major_version*100+minor_version*10+patch_version) > KRATOS_VERSION_AS_INT, "This method is deprecated for the current version. Please remove");

#endif /* KRATOS_DEFINE_DEPRECATED_H_INCLUDED  defined */
