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
//                   Gennady Markelov
//
#pragma once

#include "containers/flags.h"
#include "includes/ublas_interface.h"

#include <cstddef>
#include <memory>

namespace Kratos
{

class Serializer;

class InterfaceConstitutiveLawDimension
{
public:
    virtual ~InterfaceConstitutiveLawDimension() = default;

    [[nodiscard]] virtual Matrix MakeInterfaceConstitutiveMatrix(double NormalStiffness,
                                                                 double ShearStiffness,
                                                                 std::size_t TractionSize) const = 0;
    [[nodiscard]] virtual std::unique_ptr<InterfaceConstitutiveLawDimension> Clone() const = 0;
    [[nodiscard]] virtual std::size_t GetStrainSize() const                                = 0;
    [[nodiscard]] virtual std::size_t GetDimension() const                                 = 0;

private:
    friend class Serializer;
    virtual void save(Serializer& rSerializer) const = 0;
    virtual void load(Serializer& rSerializer)       = 0;
};

} // namespace Kratos
