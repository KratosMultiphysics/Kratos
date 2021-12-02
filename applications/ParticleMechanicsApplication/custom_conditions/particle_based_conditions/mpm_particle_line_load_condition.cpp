//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Veronika Singer
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/particle_based_conditions/mpm_particle_line_load_condition.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************

    MPMParticleLineLoadCondition::MPMParticleLineLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry )
        : Condition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    MPMParticleLineLoadCondition::MPMParticleLineLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : Condition( NewId, pGeometry, pProperties )
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************

    Condition::Pointer MPMParticleLineLoadCondition::Create(IndexType NewId,GeometryType::Pointer pGeometry,PropertiesType::Pointer pProperties) const
    {
        return Kratos::make_intrusive<MPMParticleLineLoadCondition>(NewId, pGeometry, pProperties);
    }

    //************************************************************************************
    //************************************************************************************

    Condition::Pointer MPMParticleLineLoadCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return Kratos::make_intrusive<MPMParticleLineLoadCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
    }

    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************

    MPMParticleLineLoadCondition::~MPMParticleLineLoadCondition()
    {
    }


    void MPMParticleLineLoadCondition::CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo)
    {
        if (rValues.size() != 1)
            rValues.resize(1);

        if (rVariable == MPC_AREA) {
            rValues[0] = m_area;
        }
        else {
            KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
        }
    }

    void MPMParticleLineLoadCondition::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
        std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo)
    {
        if (rValues.size() != 1)
            rValues.resize(1);

        if (rVariable == LINE_LOAD) {
            rValues[0] = m_line_load;
        }
        else if (rVariable == POINT_LOAD) {
            rValues[0] = m_point_load;
        }
        else if (rVariable == MPC_DISPLACEMENT) {
            rValues[0] = m_delta_xg;
        }
        else if (rVariable == MP_DISPLACEMENT) {
            rValues[0] = m_delta_xg;
        }
        else if (rVariable == MP_COORD || rVariable == MPC_COORD) {
            rValues[0] = m_xg;
        }
        else if (rVariable == MPC_ID_LIST) {
            rValues[0] = m_id_list;
        }
        else if (rVariable == MPC_AREA_LIST) {
            rValues[0] = m_area_list;
        }
        else {
            KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
        }
    }

    void MPMParticleLineLoadCondition::SetValuesOnIntegrationPoints(
        const Variable<double>& rVariable,
        const std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo) {
        KRATOS_ERROR_IF(rValues.size() > 1)
            << "Only 1 value per integration point allowed! Passed values vector size: "
            << rValues.size() << std::endl;

        if (rVariable == MPC_AREA) {
            m_area = rValues[0];
        }
        else {
            KRATOS_ERROR << "Variable " << rVariable << " is called in SetValuesOnIntegrationPoints, but is not implemented." << std::endl;
        }
    }

    
    void MPMParticleLineLoadCondition::SetValuesOnIntegrationPoints(
        const Variable<array_1d<double, 3 > >& rVariable,
        const std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR_IF(rValues.size() > 1)
            << "Only 1 value per integration point allowed! Passed values vector size: "
            << rValues.size() << std::endl;

        if (rVariable == LINE_LOAD) {
            m_line_load = rValues[0];
        }
        else if (rVariable == MPC_DISPLACEMENT) {
            m_delta_xg = rValues[0];
        }
        else if (rVariable == MP_DISPLACEMENT) {
            m_delta_xg = rValues[0];
        }
        else if (rVariable == MP_COORD || rVariable == MPC_COORD) {
            m_xg = rValues[0];
        }
        else if (rVariable == MPC_ID_LIST) {
            m_id_list = rValues[0];
        }
        else if (rVariable == MPC_AREA_LIST) {
            m_area_list = rValues[0];
        }
        else {
            KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
        }
    }

    

} // Namespace Kratos