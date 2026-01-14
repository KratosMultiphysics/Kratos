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

    static void CheckForNonZeroZCoordinateIn2D(const Geometry<Node>& rGeometry);

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

    CheckProperties(const Properties& rProperties, const std::string& rPrintName, Bounds RangeBoundsType);

    CheckProperties(const Properties& rProperties, const std::string& rPrintName, std::size_t ElementId, Bounds RangeBoundsType);

    CheckProperties(const CheckProperties&) = delete;

    CheckProperties(CheckProperties&&) = default;

    CheckProperties SingleUseBounds(Bounds RangeBoundsType) const;

    void SetNewRangeBounds(Bounds RangeBoundsType) const;

    template <typename T>
    void Check(const Variable<T>& rVariable) const
    {
        CheckAvailability(rVariable);
        CheckRangeBounds(rVariable, mDefaultLowerBound, std::nullopt);
    }

    template <typename T>
    void Check(const Variable<T>& rVariable, double UpperBound) const
    {
        CheckAvailability(rVariable);
        CheckRangeBounds(rVariable, mDefaultLowerBound, UpperBound);
    }

    template <typename T>
    void Check(const Variable<T>& rVariable, double LowerBound, double UpperBound) const
    {
        CheckAvailability(rVariable);
        CheckRangeBounds(rVariable, LowerBound, UpperBound);
    }

    template <typename T>
    void CheckAvailability(const Variable<T>& rVariable) const
    {
        if (!mrProperties.Has(rVariable))
            KRATOS_ERROR << rVariable.Name() << " does not exist" << print_property_id()
                         << print_element_id() << "." << std::endl;
    }

    template <typename T>
    void CheckAvailabilityAndNotEmpty(const Variable<T>& rVariable) const
    {
        CheckAvailability(rVariable);
        if (mrProperties[rVariable].empty())
            KRATOS_ERROR << rVariable.Name() << " is empty" << print_property_id()
                         << print_element_id() << "." << std::endl;
    }

    template <typename T, typename Eq>
    void CheckAvailabilityAndEquality(const Kratos::Variable<T>& rVariable, const Eq& rName) const
    {
        CheckAvailability(rVariable);
        std::string quotes = "";
        if constexpr (std::is_same_v<T, std::string>) {
            quotes = "\"";
        }
        KRATOS_ERROR_IF_NOT(mrProperties[rVariable] == rName)
            << rVariable.Name() << " has a value of " << quotes << mrProperties[rVariable] << quotes
            << " instead of " << quotes << rName << quotes << print_property_id()
            << print_element_id() << "." << std::endl;
    }

    template <typename T>
    requires(!std::is_same_v<T, bool>) void CheckAvailabilityAndSpecified(const Variable<T>& rVariable) const
    {
        CheckAvailability(rVariable);
        if (!mrProperties[rVariable])
            KRATOS_ERROR << rVariable.Name() << " needs to be specified" << print_property_id()
                         << print_element_id() << "." << std::endl;
    }

    void CheckPermeabilityProperties(size_t Dimension) const;

private:
    const Properties&                mrProperties;
    const std::string                mrPrintName;
    const std::optional<std::size_t> mElementId = std::nullopt;
    mutable Bounds                   mRangeBoundsType;
    const double                     mDefaultLowerBound = 0.0;

    template <typename T>
    void CheckRangeBounds(const Variable<T>& rVariable, double LowerBound, std::optional<double> UpperBound) const
    {
        const auto value = mrProperties[rVariable];
        bool       in_range;
        using enum CheckProperties::Bounds;

        switch (mRangeBoundsType) {
        case AllExclusive:
            in_range = (value > LowerBound && (!UpperBound || value < *UpperBound));
            break;
        case AllInclusive:
            in_range = (value >= LowerBound && (!UpperBound || value <= *UpperBound));
            break;
        case InclusiveLowerAndExclusiveUpper:
            in_range = (value >= LowerBound && (!UpperBound || value < *UpperBound));
            break;
        case ExclusiveLowerAndInclusiveUpper:
            in_range = (value > LowerBound && (!UpperBound || value <= *UpperBound));
            break;
        default:
            KRATOS_ERROR << " Unknown type of range bounds";
            break;
        }
        if (!in_range) {
            const auto include_lower_bound = (mRangeBoundsType == AllInclusive) ||
                                             (mRangeBoundsType == InclusiveLowerAndExclusiveUpper);
            const auto include_upper_bound =
                UpperBound && ((mRangeBoundsType == AllInclusive) ||
                               (mRangeBoundsType == ExclusiveLowerAndInclusiveUpper));
            std::ostringstream print_range;
            print_range << (include_lower_bound ? "[" : "(") << LowerBound << ", "
                        << (UpperBound ? double_to_string(*UpperBound) : "-")
                        << (include_upper_bound ? "]" : ")");
            KRATOS_ERROR << rVariable.Name() << print_property_id() << print_element_id()
                         << " has an invalid value: " << value << " is out of the range "
                         << print_range.str() << "." << std::endl;
        }
    }

    std::string double_to_string(double value) const;

    std::string print_property_id() const;

    std::string print_element_id() const;
};

} /* namespace Kratos.*/
