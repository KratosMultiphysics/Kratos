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

    std::unique_ptr<CheckProperties> AsExclusive()
    {
        return std::make_unique<CheckProperties>(mrPrintName, mrProperties, mId, false, false);
    }

    std::unique_ptr<CheckProperties> AsInclusive()
    {
        return std::make_unique<CheckProperties>(mrPrintName, mrProperties, mId, true, true);
    }

    void Check(const Variable<double>& rVariable);
    void Check(const Variable<double>& rVariable, double UpperBound);
    void Check(const Variable<double>& rVariable, double LowerBound, double UpperBound);

    template <typename T>
    void CheckAvailabilityOnly(const Variable<T>& rVariable)
    {
        if (!mrProperties.Has(rVariable))
            KRATOS_ERROR << rVariable.Name() << " does not exist in the " << mrPrintName << " "
                         << mId << "." << std::endl;
    }

private:
    const std::string mrPrintName;
    const Properties& mrProperties;
    const std::size_t mId;
    const bool        mIncludeLower;
    const bool        mIncludeUpper;
    const double      mLowerBound = 0.0;
    const double      mUpperBound = std::numeric_limits<double>::max();
    void CheckRangeBounds(const Variable<double>& rVariable, double LowerBound, double UpperBound) const;
};

} /* namespace Kratos.*/
