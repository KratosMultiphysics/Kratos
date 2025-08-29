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
    enum class Bounds {
        AllInclusive,
        AllExclusive,
        InclusiveLowerAndExclusiveUpper,
        ExclusiveLowerAndInclusiveUpper
    };

    CheckProperties(const std::string& rPrintName, const Properties& rProperties, Bounds Option)
        : mrPrintName(rPrintName), mrProperties(rProperties), mId(rProperties.Id()), mOption(Option)
    {
    }

    CheckProperties(const std::string& rPrintName, const Properties& rProperties, std::size_t Id, Bounds Option)
        : mrPrintName(rPrintName), mrProperties(rProperties), mId(Id), mOption(Option)
    {
    }

    std::unique_ptr<CheckProperties> SingleUseBounds(Bounds Option) const
    {
        return std::make_unique<CheckProperties>(mrPrintName, mrProperties, mId, Option);
    }

    void SetNewBounds(Bounds Option) { mOption = Option; }

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
    Bounds            mOption;
    const double      mLowerBound = 0.0;
    const double      mUpperBound = std::numeric_limits<double>::max();

    template <typename T>
    void CheckRangeBounds(const Variable<T>& rVariable, double LowerBound, double UpperBound) const
    {
        const auto value = mrProperties[rVariable];
        bool       in_range;
        switch (mOption) {
        case Bounds::AllExclusive:
            in_range = (value > LowerBound && value < UpperBound);
            break;
        case Bounds::AllInclusive:
            in_range = (value >= LowerBound && value <= UpperBound);
            break;
        case Bounds::InclusiveLowerAndExclusiveUpper:
            in_range = (value >= LowerBound && value < UpperBound);
            break;
        case Bounds::ExclusiveLowerAndInclusiveUpper:
            in_range = (value > LowerBound && value <= UpperBound);
        }
        if (!in_range) {
            std::ostringstream print_range;
            const auto         include_lower_bound = (mOption == Bounds::AllInclusive) ||
                                             (mOption == Bounds::InclusiveLowerAndExclusiveUpper);
            const auto include_upper_bound = (mOption == Bounds::AllInclusive) ||
                                             (mOption == Bounds::ExclusiveLowerAndInclusiveUpper);
            print_range << (include_lower_bound ? "[" : "(") << LowerBound << "; "
                        << ((UpperBound == std::numeric_limits<double>::max()) ? "-" : std::to_string(UpperBound))
                        << (include_upper_bound ? "]" : ")");
            KRATOS_ERROR << rVariable.Name() << " in the " << mrPrintName << " " << mId
                         << " has an invalid value: " << value << " out of the range "
                         << print_range.str() << "." << std::endl;
        }
    }
};

} /* namespace Kratos.*/
