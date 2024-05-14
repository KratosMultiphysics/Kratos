// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//                   John van Esch
//                   Gennady Markelov
//

#include "custom_constitutive/pw_filter_law.h"
#include "geo_mechanics_application_variables.h"
#include "includes/global_variables.h"

namespace Kratos
{

GeoPwFilterLaw::GeoPwFilterLaw() : mNumberOfDimensions{2} {}

GeoPwFilterLaw::GeoPwFilterLaw(std::size_t NumberOfDimensions)
    : mNumberOfDimensions{NumberOfDimensions}
{
    KRATOS_ERROR_IF(mNumberOfDimensions != 1)
        << "Got invalid number of dimensions. The dimension has to be 1, but got: " << mNumberOfDimensions
        << std::endl;
}

ConstitutiveLaw::Pointer GeoPwFilterLaw::Clone() const
{
    return Kratos::make_shared<GeoPwFilterLaw>(*this);
}

SizeType GeoPwFilterLaw::WorkingSpaceDimension() { return mNumberOfDimensions; }

Matrix GeoPwFilterLaw::CalculatePwFilterMatrix(const Properties& rProp, const ProcessInfo& rProcessInfo) const
{
    const double equivalent_radius_square = rProp[CROSS_AREA] / Globals::Pi;
    return ScalarMatrix(1, 1, equivalent_radius_square * 0.125);
}

} // Namespace Kratos
