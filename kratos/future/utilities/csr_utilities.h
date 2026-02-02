//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project include
#include "containers/nd_data.h"

namespace Kratos
{

///@addtogroup KratosCore
///@{

///@name Kratos Globals
///@{


///@}
///@name Type Definitions
///@{


///@}
///@name  Enum's
///@{


///@}
///@name  Functions
///@{


///@}
///@name Kratos Classes
///@{

class CsrUtilities final
{
public:
    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CsrUtilities() = delete;

    /// Destructor.
    ~CsrUtilities() = delete;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    template<class TContainerType, class TCsrMatrixType , class TIndexVectorType>
    static void GetCsrEquationIdIndices(
        const TContainerType& rContainer,
        const TCsrMatrixType& rCsrMatrix,
        NDData<IndexType>& rNdData)
    {
        // Resize NDData
        // Note that we assume the entities in the container to be of the same type
        DenseVector<unsigned int> nd_data_shape(3);
        nd_data_shape[0] = rContainer.size();
        nd_data_shape[1] = rContainer.begin()->GetDofs().size();
        nd_data_shape[2] = rContainer.begin()->GetDofs().size();
        rNdData = NdData<IndexType>(nd_data_shape);

        //
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    ///@}
    ///@name Friends
    ///@{


    ///@}
}; // Class CsrUtilities

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.
