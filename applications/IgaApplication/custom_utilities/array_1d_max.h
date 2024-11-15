//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


#if !defined(KRATOS_ARRAY_1D_MAX_H_INCLUDED )
#define  KRATOS_ARRAY_1D_MAX_H_INCLUDED

#include "includes/define.h" // For array_1d
#include <algorithm>         // For std::max_element
#include "includes/kratos_flags.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "utilities/geometry_utilities.h"


namespace Kratos
{

    class KRATOS_API(IGA_APPLICATION) Array1DMax
    {
        public:
            /// Constructor
            Array1DMax();
            // Function to find the maximum value in array_1d
            static double FindMaxInArray1D(const array_1d<double, 3>& array) {
            return *std::max_element(array.begin(), array.end());
            }
    };

} // namespace Kratos.
#endif // KRATOS_ARRAY_1D_MAX_H_INCLUDED  defined