// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_conditions/axisym_point_load_condition.h"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************

    AxisymPointLoadCondition::AxisymPointLoadCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry
        )
            : PointLoadCondition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    AxisymPointLoadCondition::AxisymPointLoadCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
        )
            : PointLoadCondition( NewId, pGeometry, pProperties )
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************

    Condition::Pointer AxisymPointLoadCondition::Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const
    {
        return Kratos::make_shared<AxisymPointLoadCondition>( NewId, pGeom, pProperties );
    }

    //************************************************************************************
    //************************************************************************************

    Condition::Pointer AxisymPointLoadCondition::Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const
    {
        return Kratos::make_shared<AxisymPointLoadCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
    }

    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************

    AxisymPointLoadCondition::~AxisymPointLoadCondition()
    {
    }

    //************************************************************************************
    //********************************* PROTECTED ****************************************
    //************************************************************************************

    double AxisymPointLoadCondition::GetPointLoadIntegrationWeight()
    {
        // We calculate the axisymmetric coefficient
        const double Radius = StructuralMechanicsMathUtilities::CalculateRadiusPoint(GetGeometry());
        const double Thickness = (GetProperties().Has( THICKNESS ) == true) ? this->GetProperties()[THICKNESS] : 1.0;
        const double AxiSymCoefficient = 2.0 * Globals::Pi * Radius/Thickness;

        return AxiSymCoefficient;
    }


    //************************************************************************************
    //************************************************************************************

    void AxisymPointLoadCondition::save( Serializer& rSerializer ) const
    {
        rSerializer.save( "Name", "AxisymPointLoadCondition" );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, PointLoadCondition );
    }

    //************************************************************************************
    //************************************************************************************

    void AxisymPointLoadCondition::load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, PointLoadCondition );
    }

} // Namespace Kratos


