//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Lorenzo Gracia
//

/* Project includes */
#include "custom_elements/small_displacement_thermo_mechanic_element.hpp"


namespace Kratos
{

// Default Constructor
SmallDisplacementThermoMechanicElement::SmallDisplacementThermoMechanicElement() : SmallDisplacementElement() {}

//----------------------------------------------------------------------------------------

//Constructor 1
SmallDisplacementThermoMechanicElement::SmallDisplacementThermoMechanicElement( IndexType NewId, GeometryType::Pointer pGeometry ) : SmallDisplacementElement( NewId, pGeometry ) {}

//----------------------------------------------------------------------------------------

//Constructor 2
SmallDisplacementThermoMechanicElement::SmallDisplacementThermoMechanicElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : SmallDisplacementElement( NewId, pGeometry, pProperties )
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

//----------------------------------------------------------------------------------------

//Destructor
SmallDisplacementThermoMechanicElement::~SmallDisplacementThermoMechanicElement() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Element::Pointer SmallDisplacementThermoMechanicElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new SmallDisplacementThermoMechanicElement( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

} // Namespace Kratos
