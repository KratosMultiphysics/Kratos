// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Mohamed Nabi
//                   John van Esch
//

#pragma once

#include "custom_utilities/condition_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/condition.h"
#include "includes/define.h"
#include "includes/process_info.h"
#include "includes/serializer.h"

namespace Kratos {

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoTCondition : public Condition {
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeoTCondition);

    GeoTCondition();

    GeoTCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    GeoTCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ~GeoTCondition() override;

    Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& rThisNodes,
                              PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GeoTCondition>(
            NewId, GetGeometry().Create(rThisNodes), pProperties);
    }

    using Condition::Create;

    void GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo&) const override
    {
        rConditionDofList = GetDofs();
    }

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult,
                          const ProcessInfo& rCurrentProcessInfo) const override;

protected:
    virtual void CalculateAll(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateRHS(VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo);

private:
    [[nodiscard]] DofsVectorType GetDofs() const
    {
        auto       result = DofsVectorType{};
        const auto nodes  = this->GetGeometry();
        std::transform(nodes.begin(), nodes.end(), std::back_inserter(result),
                       [](const auto& r_node) { return r_node.pGetDof(TEMPERATURE); });
        return result;
    }

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition)
    }
};

} // namespace Kratos
