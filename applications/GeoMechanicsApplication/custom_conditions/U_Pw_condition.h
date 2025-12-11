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
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#pragma once

// System includes

// Project includes
#include "includes/condition.h"
#include "includes/define.h"
#include "includes/process_info.h"
#include "includes/serializer.h"

// Application includes
#include "custom_utilities/condition_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) UPwCondition : public Condition
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(UPwCondition);

    UPwCondition();

    UPwCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    UPwCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ~UPwCondition() override = default;

    Condition::Pointer Create(IndexType               NewId,
                              NodesArrayType const&   ThisNodes,
                              PropertiesType::Pointer pProperties) const override;
    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    void GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo&) const override;

    IntegrationMethod GetIntegrationMethod() const override { return mThisIntegrationMethod; }

    void SetIntegrationMethod(IntegrationMethod method) { mThisIntegrationMethod = method; }

    void CalculateLocalSystem(Matrix&            rLeftHandSideMatrix,
                              Vector&            rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(Vector& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const override;

    std::string Info() const override;

protected:
    virtual void CalculateAll(Matrix& rLeftHandSideMatrix, Vector& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateRHS(Vector& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo);

    [[nodiscard]] DofsVectorType GetDofs() const;

private:
    GeometryData::IntegrationMethod mThisIntegrationMethod{Condition::GetIntegrationMethod()};

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

}; // class UPwCondition.

} // namespace Kratos.
