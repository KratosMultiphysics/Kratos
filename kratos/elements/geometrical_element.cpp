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
#include "elements/geometrical_element.h"

namespace Kratos
{
GeometricalElement::GeometricalElement(IndexType NewId)
    : BaseType(NewId)
{
}

/***********************************************************************************/
/***********************************************************************************/

GeometricalElement::GeometricalElement(
    IndexType NewId, 
    const NodesArrayType& rThisNodes
    ) : BaseType(NewId, rThisNodes)
{
}

/***********************************************************************************/
/***********************************************************************************/

GeometricalElement::GeometricalElement(
    IndexType NewId, 
    GeometryType::Pointer pGeometry
    ) : BaseType(NewId, pGeometry)
{
}

/***********************************************************************************/
/***********************************************************************************/

GeometricalElement::GeometricalElement(
    IndexType NewId, 
    GeometryType::Pointer pGeometry, 
    PropertiesType::Pointer pProperties
    ) : BaseType(NewId,pGeometry, pProperties)
{
}

/***********************************************************************************/
/***********************************************************************************/

GeometricalElement::GeometricalElement(GeometricalElement const& rOther)
    : BaseType(rOther)
{
}

/***********************************************************************************/
/***********************************************************************************/

GeometricalElement::~GeometricalElement()
{
}

/***********************************************************************************/
/***********************************************************************************/

GeometricalElement& GeometricalElement::operator=(GeometricalElement const& rOther)
{
    //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

    Element::operator=(rOther);

    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer GeometricalElement::Create(
    IndexType NewId, 
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_shared<Element>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer GeometricalElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_shared<Element>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer GeometricalElement::Clone (
    IndexType NewId, 
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Element::Pointer p_new_elem = Kratos::make_shared<Element>(NewId, GetGeometry().Create(ThisNodes), pGetProperties()); 
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this)); 
    return p_new_elem; 

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void GeometricalElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
}

void GeometricalElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
}

} // Namespace Kratos


