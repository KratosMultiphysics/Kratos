//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "new_utility.h"

namespace Kratos
{
    NewUtility::NewUtility()
    {
    }

    void NewUtility::PrintSomething()
    {
        KRATOS_WATCH("caio")
    }

    std::vector<int> NewUtility::CreateBinaryArray(std::size_t length)
    {

        std::vector<int> binary_array(length, 0);

        return binary_array;
    }

    std::vector<int> NewUtility::ModifyArray(std::vector<int>& my_binary_array, std::vector<int>& EquationId)
    {
        for (int i_eq = 0; i_eq < EquationId.size(); i_eq++ ){
            if (my_binary_array[EquationId[i_eq]] == 0) {
                my_binary_array[EquationId[i_eq]] = 1 ;
            }
        }
        // KRATOS_WATCH(my_binary_array)
        return my_binary_array;
    }

}  // namespace Kratos.
