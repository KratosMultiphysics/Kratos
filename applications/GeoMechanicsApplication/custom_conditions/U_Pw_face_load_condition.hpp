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

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/U_Pw_condition.hpp"
#include "custom_utilities/condition_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) UPwFaceLoadCondition : public UPwCondition<TDim, TNumNodes>
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(UPwFaceLoadCondition);

    using IndexType      = std::size_t;
    using PropertiesType = Properties;
    using GeometryType   = Geometry<Node>;
    using NodesArrayType = GeometryType::PointsArrayType;

    // Default constructor
    UPwFaceLoadCondition() : UPwCondition<TDim, TNumNodes>() {}

    // Constructor 1
    UPwFaceLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : UPwCondition<TDim, TNumNodes>(NewId, pGeometry)
    {
    }

    // Constructor 2
    UPwFaceLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : UPwCondition<TDim, TNumNodes>(NewId, pGeometry, pProperties)
    {
    }

    Condition::Pointer Create(IndexType               NewId,
                              NodesArrayType const&   ThisNodes,
                              PropertiesType::Pointer pProperties) const override;

    std::string Info() const override;

protected:
    void CalculateRHS(Vector& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo) override;

private:
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition)
    }
}; // class UPwFaceLoadCondition.

} // namespace Kratos.