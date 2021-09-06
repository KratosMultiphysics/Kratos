// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:		 BSD License
//					 license: ContactStructuralMechanicsApplication/license.txt
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

void PairedCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::Initialize(rCurrentProcessInfo);

    // The normal of the paired condition
    const auto& r_paired_geometry = GetPairedGeometry();
    GeometryType::CoordinatesArrayType aux_coords;
    r_paired_geometry.PointLocalCoordinates(aux_coords, r_paired_geometry.Center());
    mPairedNormal = r_paired_geometry.UnitNormal(aux_coords);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void PairedCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::InitializeSolutionStep(rCurrentProcessInfo);

    // The normal of the paired condition
    const auto& r_paired_geometry = GetPairedGeometry();
    GeometryType::CoordinatesArrayType aux_coords;
    r_paired_geometry.PointLocalCoordinates(aux_coords, r_paired_geometry.Center());
    mPairedNormal = r_paired_geometry.UnitNormal(aux_coords);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void PairedCondition::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::InitializeNonLinearIteration(rCurrentProcessInfo);

    // We update the normals if necessary
    const auto normal_variation = rCurrentProcessInfo.Has(CONSIDER_NORMAL_VARIATION) ? static_cast<NormalDerivativesComputation>(rCurrentProcessInfo.GetValue(CONSIDER_NORMAL_VARIATION)) : NO_DERIVATIVES_COMPUTATION;
    if (normal_variation != NO_DERIVATIVES_COMPUTATION) {
        const auto& r_paired_geometry = GetPairedGeometry();
        GeometryType::CoordinatesArrayType aux_coords;
        r_paired_geometry.PointLocalCoordinates(aux_coords, r_paired_geometry.Center());
        mPairedNormal = r_paired_geometry.UnitNormal(aux_coords);
    }

    KRATOS_CATCH( "" );
}

} // Namespace Kratos
