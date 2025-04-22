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
//                   Wijtze Pieter Kikstra
//

#pragma once

#include "includes/ublas_interface.h"

namespace Kratos
{

class Serializer;

class YieldSurface
{
public:
    virtual ~YieldSurface() = default;

    [[nodiscard]] virtual double YieldFunctionValue(const Vector& rSigmaTau) const       = 0;
    [[nodiscard]] virtual Vector DerivativeOfFlowFunction(const Vector& rPrincipalStress) const = 0;

private:
    friend class Serializer;
    virtual void save(Serializer& rSerializer) const = 0;
    virtual void load(Serializer& rSerializer)       = 0;
};

} // namespace Kratos.
