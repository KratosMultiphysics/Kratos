// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "principal_stresses.hpp"

namespace Kratos::Geo
{

PrincipalStresses::PrincipalStresses(const std::initializer_list<double>& rValues)
{
    KRATOS_DEBUG_ERROR_IF_NOT(rValues.size() == msVectorSize)
        << "PrincipalStresses can only be initialized with a vector of size " << msVectorSize
        << ", got " << rValues.size() << std::endl;
    std::ranges::copy(rValues, mValues.begin());
}

const PrincipalStresses::InternalVectorType& PrincipalStresses::Values() const { return mValues; }

PrincipalStresses::InternalVectorType& PrincipalStresses::Values() { return mValues; }

} // namespace Kratos::Geo