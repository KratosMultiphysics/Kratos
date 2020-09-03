//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
//

#ifndef KRATOS_MPM_MPI_SEARCH_IMPORT_H
#define KRATOS_MPM_MPI_SEARCH_IMPORT_H

// Project inlcudes
#include "includes/model_part.h"

namespace Kratos {

typedef ModelPart::ElementsContainerType ElementsArrayType;

namespace MPM_MPI_SEARCH {

    void SearchElements(ModelPart& rMPMModelPart, ModelPart& rBackgroundGridModelPart,
                        std::vector<Element::Pointer>& rMissingElements,
                        std::vector<Condition::Pointer>& rMissingConditions,
                        const std::size_t MaxNumberOfResults, const double Tolerance);

} // namespace MPM_MPI_SEARCH
} // namespace Kratos

#endif // KRATOS_MPM_MPI_SEARCH_IMPORT_H