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
#include "custom_conditions/U_Pw_condition.h"
#include "custom_conditions/U_Pw_face_load_interface_condition.h"
#include "custom_utilities/element_utilities.hpp"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) UPwNormalFluxInterfaceCondition
    : public UPwFaceLoadInterfaceCondition<TDim, TNumNodes>
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(UPwNormalFluxInterfaceCondition);

    using IndexType      = std::size_t;
    using PropertiesType = Properties;
    using GeometryType   = Geometry<Node>;
    using NodesArrayType = GeometryType::PointsArrayType;

    UPwNormalFluxInterfaceCondition();

    UPwNormalFluxInterfaceCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    UPwNormalFluxInterfaceCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    Condition::Pointer Create(IndexType               NewId,
                              NodesArrayType const&   ThisNodes,
                              PropertiesType::Pointer pProperties) const override;

    std::string Info() const override;

protected:
    void CalculateRHS(Vector& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo) override;

private:
    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

}; // class UPwNormalFluxInterfaceCondition.

} // namespace Kratos.
