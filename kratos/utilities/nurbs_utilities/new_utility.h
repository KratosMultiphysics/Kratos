//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


#if !defined(KRATOS_NEW_UTILITY_H_INCLUDED )
#define  KRATOS_NEW_UTILITY_H_INCLUDED

// System includes
#include "includes/define.h"

namespace Kratos
{
    ///@name Kratos Classes
    ///@{

    class KRATOS_API(IGA_APPLICATION) NewUtility
    {
    public:
        // Costructor
        NewUtility();

        void PrintSomething();

        std::vector<int> CreateBinaryArray(std::size_t length);

        std::vector<int> ModifyArray(std::vector<int>& my_binary_array, std::vector<int>& EquationId);

    private:

    }; // Class NewUtility

}  // namespace Kratos

#endif // KRATOS_DIRECTOR_UTILITIES_H_INCLUDED  defined