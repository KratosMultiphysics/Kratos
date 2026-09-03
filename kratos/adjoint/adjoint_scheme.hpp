//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#pragma once

// --- Kratos Core Includes ---
#include "solving_strategies/schemes/scheme.h"
#include "adjoint/response_function.hpp"

// --- STL Includes ---
#include <span>


namespace Kratos {


template <class TSparse, class TDense>
class KRATOS_API(KRATOS_CORE) AdjointScheme : public Scheme<TSparse,TDense> {
public:
    KRATOS_CLASS_POINTER_DEFINITION(AdjointScheme);

    using Scheme<TSparse,TDense>::Scheme;

    AdjointScheme(
        Parameters Settings,
        ResponseFunction::Pointer pResponse)
            :   Scheme<TSparse,TDense>(Settings),
                mpResponseFunction() {
                    this->SetResponseFunction(pResponse);
    }

    void SetResponseFunction(ResponseFunction::Pointer pResponse) noexcept {
        mpResponseFunction = pResponse;
    }

    const ResponseFunction& GetResponseFunction() const {
        KRATOS_ERROR_IF_NOT(mpResponseFunction) << "no response function is set";
        return *mpResponseFunction;
    }

protected:
    template <class T>
    requires (std::is_same_v<T,Element> || std::is_same_v<T,Condition>)
    void GetEquationIds(
        std::span<Dof<double>::EquationIdType> EquationIds,
        std::span<const IAdjoint::DynamicVariable> Variables,
        const T& rElement) const {
            KRATOS_ERROR_IF_NOT(EquationIds.size() == Variables.size());
            KRATOS_TRY
                const Geometry<Node>& r_geometry = rElement.GetGeometry();
                std::transform(
                    Variables.begin(),
                    Variables.end(),
                    EquationIds.begin(),
                    [&r_geometry] (const IAdjoint::DynamicVariable& r_variable) -> Dof<double>::EquationIdType {
                        const auto i_node = r_variable.GetDynamicIndex();
                        KRATOS_DEBUG_ERROR_IF_NOT(0 <= i_node);
                        KRATOS_DEBUG_ERROR_IF_NOT(i_node < r_geometry.size());
                        return r_geometry[r_variable.GetDynamicIndex()].GetDof(r_variable).EquationId();
                    });
            KRATOS_CATCH("")
    }

private:
    ResponseFunction::Pointer mpResponseFunction;
}; // class AdjointScheme


} // namespace Kratos
