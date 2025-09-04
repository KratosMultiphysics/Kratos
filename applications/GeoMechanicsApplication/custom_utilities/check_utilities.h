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

    CheckProperties(const Properties& rProperties, const std::string& rPrintName, Bounds RangeBoundsType)
        : mrProperties(rProperties), mrPrintName(rPrintName), mId(rProperties.Id()), mRangeBoundsType(RangeBoundsType)
    {
    }

    CheckProperties(const Properties& rProperties, const std::string& rPrintName, std::size_t Id, Bounds RangeBoundsType)
        : mrProperties(rProperties), mrPrintName(rPrintName), mId(Id), mRangeBoundsType(RangeBoundsType)
    {
    }

    CheckProperties(const CheckProperties&) = delete;

    CheckProperties(CheckProperties&&) = default;

    CheckProperties SingleUseBounds(Bounds RangeBoundsType) const
    {
        return CheckProperties(mrProperties, mrPrintName, mId, RangeBoundsType);
    }

    void SetNewRangeBounds(Bounds RangeBoundsType) const
    {
        mHistoryOfTypes.push_back(mRangeBoundsType);
        mRangeBoundsType = RangeBoundsType;
    }

    void RestorePreviousBounds() const;

    template <typename T>
    void Check(const Variable<T>& rVariable) const
    {
        CheckAvailabilityOnly(rVariable);
        CheckRangeBounds(rVariable, mDefaultLowerBound, mDefaultUpperBound);
    }

    template <typename T>
    void Check(const Variable<T>& rVariable, double UpperBound) const
    {
        CheckAvailabilityOnly(rVariable);
        CheckRangeBounds(rVariable, mDefaultLowerBound, UpperBound);
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

    void CheckPermeabilityProperties(size_t Dimension) const;

private:
    const Properties&           mrProperties;
    const std::string           mrPrintName;
    const std::size_t           mId;
    mutable Bounds              mRangeBoundsType;
    const double                mDefaultLowerBound = 0.0;
    const double                mDefaultUpperBound = std::numeric_limits<double>::max();
    mutable std::vector<Bounds> mHistoryOfTypes;

    template <typename T>
    void CheckRangeBounds(const Variable<T>& rVariable, double LowerBound, double UpperBound) const
    {
        const auto value = mrProperties[rVariable];
        bool       in_range;
        using enum CheckProperties::Bounds;
        switch (mRangeBoundsType) {
        case AllExclusive:
            in_range = (value > LowerBound && value < UpperBound);
            break;
        case AllInclusive:
            in_range = (value >= LowerBound && value <= UpperBound);
            break;
        case InclusiveLowerAndExclusiveUpper:
            in_range = (value >= LowerBound && value < UpperBound);
            break;
        case ExclusiveLowerAndInclusiveUpper:
            in_range = (value > LowerBound && value <= UpperBound);
            break;
        default:
            KRATOS_ERROR << " Unknown type of range bounds";
            break;
        }
        if (!in_range) {
            std::ostringstream print_range;
            const auto         include_lower_bound = (mRangeBoundsType == AllInclusive) ||
                                             (mRangeBoundsType == InclusiveLowerAndExclusiveUpper);
            const auto include_upper_bound = (mRangeBoundsType == AllInclusive) ||
                                             (mRangeBoundsType == ExclusiveLowerAndInclusiveUpper);
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
