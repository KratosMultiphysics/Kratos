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
#include "custom_conditions/particle_based_conditions/mpm_particle_point_condition.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************

    MPMParticlePointCondition::MPMParticlePointCondition( IndexType NewId, GeometryType::Pointer pGeometry )
        : Condition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    MPMParticlePointCondition::MPMParticlePointCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : Condition( NewId, pGeometry, pProperties )
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************

    Condition::Pointer MPMParticlePointCondition::Create(IndexType NewId,GeometryType::Pointer pGeometry,PropertiesType::Pointer pProperties) const
    {
        return Kratos::make_intrusive<MPMParticlePointCondition>(NewId, pGeometry, pProperties);
    }

    //************************************************************************************
    //************************************************************************************

    Condition::Pointer MPMParticlePointCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return Kratos::make_intrusive<MPMParticlePointCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
    }

    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************

    MPMParticlePointCondition::~MPMParticlePointCondition()
    {
    }


    void MPMParticlePointCondition::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        GeneralVariables Variables;

        Variables.CurrentDisp = CalculateCurrentDisp(Variables.CurrentDisp, rCurrentProcessInfo);

        array_1d<double,3> delta_xg = ZeroVector(3);
        array_1d<double, 3 > MPC_velocity = ZeroVector(3);
        
        Variables.N = row(GetGeometry().ShapeFunctionsValues(), 0);

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            if (Variables.N[i] > std::numeric_limits<double>::epsilon() )
            {
                auto r_geometry = GetGeometry();
                array_1d<double, 3 > nodal_velocity = ZeroVector(3);
                if (r_geometry[i].SolutionStepsDataHas(VELOCITY))
                    nodal_velocity = r_geometry[i].FastGetSolutionStepValue(VELOCITY);
                for ( unsigned int j = 0; j < dimension; j++ )
                {
                    delta_xg[j] += Variables.N[i] * Variables.CurrentDisp(i,j);
                    MPC_velocity[j] += Variables.N[i] * nodal_velocity[j];
                }
            }
        }

        // Update the Material Point Condition Position
        m_displacement = delta_xg;
        m_velocity = MPC_velocity;
    
    KRATOS_CATCH( "" )
}



    Matrix& MPMParticlePointCondition::CalculateCurrentDisp(Matrix & rCurrentDisp, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        GeometryType& r_geometry = GetGeometry();
        const unsigned int number_of_nodes = r_geometry.PointsNumber();
        const unsigned int dimension = r_geometry.WorkingSpaceDimension();

        rCurrentDisp = ZeroMatrix(number_of_nodes, dimension);

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            const array_1d<double, 3 > & current_displacement  = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);

            for ( unsigned int j = 0; j < dimension; j++ )
            {
                rCurrentDisp(i,j) = current_displacement[j];
            }
        }

        return rCurrentDisp;

        KRATOS_CATCH( "" )
    }

    void MPMParticlePointCondition::FinalizeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
    {
        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        GeneralVariables Variables;

        Variables.CurrentDisp = CalculateCurrentDisp(Variables.CurrentDisp, rCurrentProcessInfo);

        array_1d<double,3> delta_xg = ZeroVector(3);
        array_1d<double, 3 > MPC_velocity = ZeroVector(3);
        

        Variables.N = row(GetGeometry().ShapeFunctionsValues(), 0);

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            if (Variables.N[i] > std::numeric_limits<double>::epsilon() )
            {
                auto r_geometry = GetGeometry();
                array_1d<double, 3 > nodal_velocity = ZeroVector(3);
                if (r_geometry[i].SolutionStepsDataHas(VELOCITY))
                    nodal_velocity = r_geometry[i].FastGetSolutionStepValue(VELOCITY);
                for ( unsigned int j = 0; j < dimension; j++ )
                {
                    delta_xg[j] += Variables.N[i] * Variables.CurrentDisp(i,j);
                    MPC_velocity[j] += Variables.N[i] * nodal_velocity[j];
                }
            }
        }
        
        // Update the Material Point Condition Position
        m_xg += delta_xg ;
        m_velocity = MPC_velocity;
    }

    void MPMParticlePointCondition::CalculateOnIntegrationPoints(
        const Variable<int>& rVariable,
        std::vector<int>& rValues,
        const ProcessInfo& rCurrentProcessInfo)
    {
        if (rValues.size() != 1)
            rValues.resize(1);

        if (rVariable == MPC_CORRESPONDING_CONDITION_ID) {
            rValues[0] = m_corresponding_condition_id;
        }
        else
        {
            KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
        }
    }

    void MPMParticlePointCondition::CalculateOnIntegrationPoints(
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

    void MPMParticlePointCondition::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
        std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo)
    {
        if (rValues.size() != 1)
            rValues.resize(1);

        if (rVariable == NORMAL) {
            rValues[0] = m_normal;
        }
        else if (rVariable == MPC_CONTACT_FORCE) {
            rValues[0] = m_contact_force;
        }
        else if (rVariable == MPC_DISPLACEMENT) {
            rValues[0] = m_displacement;
        }
        else if (rVariable == MP_COORD || rVariable == MPC_COORD) {
            rValues[0] = m_xg;
        }
        else if (rVariable == MPC_IMPOSED_DISPLACEMENT) {
            rValues[0] = m_imposed_displacement;
        }
        else if (rVariable == MPC_VELOCITY) {
            rValues[0] = m_velocity;
        }
        else if (rVariable == MPC_ACCELERATION) {
            rValues[0] = m_acceleration;
        }
        else if (rVariable == MPC_IMPOSED_VELOCITY) {
            rValues[0] = m_imposed_velocity;
        }
        else if (rVariable == MPC_IMPOSED_ACCELERATION) {
            rValues[0] = m_imposed_acceleration;
        }
        else {
            KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
        }
    }

    void MPMParticlePointCondition::SetValuesOnIntegrationPoints(
        const Variable<int>& rVariable,
        const std::vector<int>& rValues,
        const ProcessInfo& rCurrentProcessInfo) {
        
        KRATOS_ERROR_IF(rValues.size() > 1)
            << "Only 1 value per integration point allowed! Passed values vector size: "
            << rValues.size() << std::endl;

        if (rVariable == MPC_CORRESPONDING_CONDITION_ID) {
            m_corresponding_condition_id = rValues[0];
        }
    }

    void MPMParticlePointCondition::SetValuesOnIntegrationPoints(
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

    
    void MPMParticlePointCondition::SetValuesOnIntegrationPoints(
        const Variable<array_1d<double, 3 > >& rVariable,
        const std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR_IF(rValues.size() > 1)
            << "Only 1 value per integration point allowed! Passed values vector size: "
            << rValues.size() << std::endl;

        if (rVariable == NORMAL) {
            m_normal = rValues[0];
        }
        else if (rVariable == MPC_DISPLACEMENT) {
            m_displacement = rValues[0];
        }
        else if (rVariable == MP_COORD || rVariable == MPC_COORD) {
            m_xg = rValues[0];
        }
        else if (rVariable == MPC_CONTACT_FORCE) {
            m_contact_force = rValues[0];
        }
        else if (rVariable == MPC_IMPOSED_DISPLACEMENT) {
            m_imposed_displacement = rValues[0];
        }
        else if (rVariable == MPC_VELOCITY) {
            m_velocity = rValues[0];
        }
        else if (rVariable == MPC_ACCELERATION) {
            m_acceleration = rValues[0];
        }
        else if (rVariable == MPC_IMPOSED_VELOCITY) {
            m_imposed_velocity = rValues[0];
        }
        else if (rVariable == MPC_IMPOSED_ACCELERATION) {
            m_imposed_acceleration = rValues[0];
        }
        else {
            KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
        }
    }

    

} // Namespace Kratos