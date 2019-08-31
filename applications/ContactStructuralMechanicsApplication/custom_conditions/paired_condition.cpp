// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_conditions/paired_condition.h"
#include "contact_structural_mechanics_application_variables.h"

namespace Kratos
{
/************************************* OPERATIONS **********************************/
/***********************************************************************************/

Condition::Pointer PairedCondition::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties ) const
{
    auto p_geometry = this->GetParentGeometry().Create( rThisNodes );
    return Kratos::make_intrusive< PairedCondition >( NewId, p_geometry, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer PairedCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive< PairedCondition >( NewId, pGeometry, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer PairedCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties,
    GeometryType::Pointer pPairedGeom) const
{
    return Kratos::make_intrusive< PairedCondition >( NewId, pGeometry, pProperties, pPairedGeom);
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

PairedCondition::~PairedCondition( )
= default;

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

void PairedCondition::Initialize( )
{
    KRATOS_TRY;

    BaseType::Initialize();

    // Setting paired condition
    if (!this->Has(PAIRED_GEOMETRY)) {
        this->SetValue(PAIRED_GEOMETRY, mpPairedGeometry);
    }

    KRATOS_CATCH( "" );
}

} // Namespace Kratos
