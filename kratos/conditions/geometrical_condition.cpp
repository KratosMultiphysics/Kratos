//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//			 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "conditions/geometrical_condition.h"

namespace Kratos
{
GeometricalCondition::GeometricalCondition(IndexType NewId)
    : BaseType(NewId)
{
}

/***********************************************************************************/
/***********************************************************************************/

GeometricalCondition::GeometricalCondition(
    IndexType NewId, 
    const NodesArrayType& rThisNodes
    ) : BaseType(NewId, rThisNodes)
{
}

/***********************************************************************************/
/***********************************************************************************/

GeometricalCondition::GeometricalCondition(
    IndexType NewId, 
    GeometryType::Pointer pGeometry
    ) : BaseType(NewId, pGeometry)
{
}

/***********************************************************************************/
/***********************************************************************************/

GeometricalCondition::GeometricalCondition(
    IndexType NewId, 
    GeometryType::Pointer pGeometry, 
    PropertiesType::Pointer pProperties
    ) : BaseType(NewId,pGeometry, pProperties)
{
}

/***********************************************************************************/
/***********************************************************************************/

GeometricalCondition::GeometricalCondition(GeometricalCondition const& rOther)
    : BaseType(rOther)
{
}

/***********************************************************************************/
/***********************************************************************************/

GeometricalCondition::~GeometricalCondition()
{
}

/***********************************************************************************/
/***********************************************************************************/

GeometricalCondition& GeometricalCondition::operator=(GeometricalCondition const& rOther)
{
    //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

    Condition::operator=(rOther);

    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer GeometricalCondition::Create(
    IndexType NewId, 
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_shared<GeometricalCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer GeometricalCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_shared<GeometricalCondition>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer GeometricalCondition::Clone (
    IndexType NewId, 
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_shared<GeometricalCondition>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void GeometricalCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
}

void GeometricalCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
}

} // Namespace Kratos


