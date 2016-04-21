 // KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

// System includes


// External includes


// Project includes
#include "custom_conditions/mortar_contact_3D_condition.hpp"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"
#include "utilities/math_utils.h"
#include "structural_mechanics_application.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

MortarContact3DCondition::MortarContact3DCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

MortarContact3DCondition::MortarContact3DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
}

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

MortarContact3DCondition::MortarContact3DCondition( MortarContact3DCondition const& rOther )
    : Condition(rOther)
{
}

/************************************* OPERATIONS **********************************/
/***********************************************************************************/

Condition::Pointer MortarContact3DCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return boost::make_shared< MortarContact3DCondition >(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

MortarContact3DCondition::~MortarContact3DCondition()
{
}

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::Initialize()
{
    KRATOS_TRY;

    KRATOS_CATCH( "" );
}


/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;


    KRATOS_CATCH( "" );
}


/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // ADD CONTENT!!!!!

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    
    // ADD CONTENT!!!!!

    KRATOS_CATCH( "" );
}


/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{

  // ADD CONTENT!!!!!
  
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::GetDofList(DofsVectorType& rConditionalDofList,ProcessInfo& rCurrentProcessInfo)
{

  // ADD CONTENT!!!!!
  
}

//*********************************GET DOUBLE VALUE***********************************
/***********************************************************************************/

void MortarContact3DCondition::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
							  std::vector<double>& rValues,
							  const ProcessInfo& rCurrentProcessInfo )
{ 
    // ADD CONTENT!!!!!
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

     // ADD CONTENT!!!!!
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
}

void MortarContact3DCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
}

} // Namespace Kratos