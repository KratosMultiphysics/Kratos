//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


// System includes

// External includes

// Project includes
#include "custom_conditions/mpm_axisym_point_load_condition.h"
#include "custom_utilities/particle_mechanics_math_utilities.h"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************

    MPMAxisymPointLoadCondition::MPMAxisymPointLoadCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry
        )
            : MPMPointLoadCondition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    MPMAxisymPointLoadCondition::MPMAxisymPointLoadCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
        )
            : MPMPointLoadCondition( NewId, pGeometry, pProperties )
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************

    Condition::Pointer MPMAxisymPointLoadCondition::Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const
    {
        return Kratos::make_shared<MPMAxisymPointLoadCondition>( NewId, pGeom, pProperties );
    }

    //************************************************************************************
    //************************************************************************************

    Condition::Pointer MPMAxisymPointLoadCondition::Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const
    {
        return Kratos::make_shared<MPMAxisymPointLoadCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
    }

    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************

    MPMAxisymPointLoadCondition::~MPMAxisymPointLoadCondition()
    {
    }

    //************************************************************************************
    //********************************* PROTECTED ****************************************
    //************************************************************************************

    double MPMAxisymPointLoadCondition::GetPointLoadIntegrationWeight()
    {
        // We calculate the axisymmetric coefficient
        const double Radius = StructuralMechanicsMathUtilities::CalculateRadiusPoint(GetGeometry());
        const double Thickness = (GetProperties().Has( THICKNESS ) == true) ? this->GetProperties()[THICKNESS] : 1.0;
        const double AxiSymCoefficient = 2.0 * Globals::Pi * Radius/Thickness;

        return AxiSymCoefficient;
    }


    //************************************************************************************
    //************************************************************************************

    void MPMAxisymPointLoadCondition::save( Serializer& rSerializer ) const
    {
        rSerializer.save( "Name", "MPMAxisymPointLoadCondition" );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMPointLoadCondition );
    }

    //************************************************************************************
    //************************************************************************************

    void MPMAxisymPointLoadCondition::load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMPointLoadCondition );
    }

} // Namespace Kratos


