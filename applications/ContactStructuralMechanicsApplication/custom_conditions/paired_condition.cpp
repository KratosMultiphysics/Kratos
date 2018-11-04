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
    return Kratos::make_shared< PairedCondition >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer PairedCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_shared< PairedCondition >( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer PairedCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties,
    GeometryType::Pointer pPairedGeom) const
{
    return Kratos::make_shared< PairedCondition >( NewId, pGeom, pProperties, pPairedGeom);
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
    
    if (mpPairedGeometry == nullptr) {
        if (this->Has(PAIRED_GEOMETRY)) {
            mpPairedGeometry = this->GetValue(PAIRED_GEOMETRY);
        } else {
            KRATOS_ERROR << "WARNING:: PAIRED GEOMETRY NOT DEFINED" << std::endl;
        }
    }
    
    KRATOS_CATCH( "" );
}

} // Namespace Kratos
