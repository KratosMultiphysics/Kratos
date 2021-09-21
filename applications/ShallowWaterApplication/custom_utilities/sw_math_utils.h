//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_SW_MATH_UTILS_H_INCLUDED
#define KRATOS_SW_MATH_UTILS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class SwMathUtils
 * @ingroup ShallowWaterApplication
 * @brief Various mathematical utilities functions
 * @author Miguel Maso Sotomayor
 */
template<class TDataType>
class SwMathUtils
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    template<class TMatrixTypeA, class TMatrixTypeB>
    static void AddMatrix(
        TMatrixTypeA& rDestination,
        const TMatrixTypeB& rInput,
        const IndexType InitialRow,
        const IndexType InitialCol)
    {
        KRATOS_TRY
        for(IndexType i = 0; i < rInput.size1(); ++i) {
            for(IndexType j = 0; j < rInput.size2(); ++j) {
                rDestination(InitialRow+i, InitialCol+j) += rInput(i,j);
            }
        }
        KRATOS_CATCH("")
    }

    template<class TVectorTypeA, class TVectorTypeB>
    static void AddVector(
        TVectorTypeA& rDestination,
        const TVectorTypeB& rInput,
        const IndexType InitialIndex)
    {
        KRATOS_TRY
        for(IndexType i = 0; i < rInput.size(); ++i) {
            rDestination[InitialIndex+i] += rInput[i];
        }
        KRATOS_CATCH("")
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

private:
    ///@name Un accessible methods
    ///@{

    SwMathUtils(void);

    SwMathUtils(SwMathUtils& rSource);

    ///@}

}; // Class SwMathUtils

///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // KRATOS_SW_MATH_UTILS_H_INCLUDED  defined


