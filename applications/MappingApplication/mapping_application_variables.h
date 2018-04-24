//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"


#if !defined(KRATOS_MAPPING_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_MAPPING_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_application.h"
#include "custom_utilities/mapper_utilities.h"

namespace Kratos
{
    typedef std::vector<MapperInterfaceInfo::Pointer> InterfaceInfoPointerVector;
    KRATOS_DEFINE_APPLICATION_VARIABLE( MAPPING_APPLICATION, int, INTERFACE_EQUATION_ID )
    KRATOS_DEFINE_APPLICATION_VARIABLE( MAPPING_APPLICATION, double/*InterfaceInfoPointerVector*/, INTERFACE_INFO )

    // TODO
    KRATOS_DEFINE_APPLICATION_VARIABLE( MAPPING_APPLICATION, double, NEIGHBOR_RANK)
    KRATOS_DEFINE_APPLICATION_VARIABLE( MAPPING_APPLICATION, NEIGHBOR_COORDINATES)
}

#endif	/* KRATOS_MAPPING_APPLICATION_VARIABLES_H_INCLUDED */