// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#pragma once

// Project includes
#include "geo_aliases.h"
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/node.h"
#include "includes/properties.h"
#include "includes/variables.h"

#include <optional>
#include <string>

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) CheckUtilities
{
public:
    static void CheckDomainSize(double                            DomainSize,
                                std::size_t                       Id,
                                const std::optional<std::string>& PrintName = std::nullopt);

    static void CheckHasNodalSolutionStepData(const Geometry<Node>&             rGeometry,
                                              const Geo::ConstVariableDataRefs& rVariableRefs);
    static void CheckHasDofs(const Geometry<Node>& rGeometry, const Geo::ConstVariableDataRefs& rVariableRefs);

private:
    static std::string PrintVectorContent(const std::vector<size_t>& rVector);
}; /* Class CheckUtilities*/

class KRATOS_API(GEO_MECHANICS_APPLICATION) CheckProperties
{
public:
    CheckProperties(const std::string& rPrintName,
                    const Properties&  rProperties,
                    bool               IncludeLower = false,
                    bool               IncludeUpper = false)
        : mrPrintName(rPrintName),
          mrProperties(rProperties),
          mId(rProperties.Id()),
          mIncludeLower(IncludeLower),
          mIncludeUpper(IncludeUpper)
    {
    }

    CheckProperties(const std::string& rPrintName,
                    const Properties&  rProperties,
                    std::size_t        Id,
                    bool               IncludeLower = false,
                    bool               IncludeUpper = false)
        : mrPrintName(rPrintName), mrProperties(rProperties), mId(Id), mIncludeLower(IncludeLower), mIncludeUpper(IncludeUpper)
    {
    }

    std::unique_ptr<CheckProperties> AsExclusive() const
    {
        return std::make_unique<CheckProperties>(mrPrintName, mrProperties, mId, false, false);
    }

    std::unique_ptr<CheckProperties> AsInclusive() const
    {
        return std::make_unique<CheckProperties>(mrPrintName, mrProperties, mId, true, true);
    }

    template <typename T>
    void Check(const Variable<T>& rVariable) const
    {
        CheckAvailabilityOnly(rVariable);
        CheckRangeBounds(rVariable, mLowerBound, mUpperBound);
    }

    template <typename T>
    void Check(const Variable<T>& rVariable, double UpperBound) const
    {
        CheckAvailabilityOnly(rVariable);
        CheckRangeBounds(rVariable, mLowerBound, UpperBound);
    }

    template <typename T>
    void Check(const Variable<T>& rVariable, double LowerBound, double UpperBound) const
    {
        CheckAvailabilityOnly(rVariable);
        CheckRangeBounds(rVariable, LowerBound, UpperBound);
    }

    template <typename T>
    void CheckAvailabilityOnly(const Variable<T>& rVariable) const
    {
        if (!mrProperties.Has(rVariable))
            KRATOS_ERROR << rVariable.Name() << " does not exist in the " << mrPrintName << " "
                         << mId << "." << std::endl;
    }

    template <typename T>
    void CheckAvailabilityAndEmpty(const Variable<T>& rVariable) const
    {
        CheckAvailabilityOnly(rVariable);
        if (mrProperties[rVariable].empty())
            KRATOS_ERROR << rVariable.Name() << " is empty in the " << mrPrintName << " " << mId
                         << "." << std::endl;
    }

private:
    const std::string mrPrintName;
    const Properties& mrProperties;
    const std::size_t mId;
    const bool        mIncludeLower;
    const bool        mIncludeUpper;
    const double      mLowerBound = 0.0;
    const double      mUpperBound = std::numeric_limits<double>::max();

    template <typename T>
    void CheckRangeBounds(const Variable<T>& rVariable, double LowerBound, double UpperBound) const
    {
        const auto value = mrProperties[rVariable];
        bool       in_range;
        if (!mIncludeLower && !mIncludeUpper) { // Open interval
            in_range = (value > LowerBound && value < UpperBound);
        } else if (mIncludeLower && mIncludeUpper) { // Closed interval
            in_range = (value >= LowerBound && value <= UpperBound);
        } else if (mIncludeLower) { // Left-closed right-open
            in_range = (value >= LowerBound && value < UpperBound);
        } else { // Right-closed left-open
            in_range = (value > LowerBound && value <= UpperBound);
        }
        if (!in_range) {
            std::ostringstream print_range;
            print_range << (mIncludeLower ? "[" : "(") << LowerBound << "; "
                        << ((UpperBound == std::numeric_limits<double>::max()) ? "-" : std::to_string(UpperBound))
                        << (mIncludeUpper ? "]" : ")");
            KRATOS_ERROR << rVariable.Name() << " in the " << mrPrintName << " " << mId
                         << " has an invalid value: " << value << " out of the range "
                         << print_range.str() << "." << std::endl;
        }
    }
};

} /* namespace Kratos.*/
