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
#include "custom_conditions/particle_based_conditions/lagrangian_moving_particle.h"
#include "includes/kratos_flags.h"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************

    MPMLagrangianMovingParticle::MPMLagrangianMovingParticle( IndexType NewId, GeometryType::Pointer pGeometry )
        : MPMParticleBaseLoadCondition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    MPMLagrangianMovingParticle::MPMLagrangianMovingParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : MPMParticleBaseLoadCondition( NewId, pGeometry, pProperties )
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************

    Condition::Pointer MPMLagrangianMovingParticle::Create(IndexType NewId,GeometryType::Pointer pGeometry,PropertiesType::Pointer pProperties) const
    {
        return Kratos::make_intrusive<MPMLagrangianMovingParticle>(NewId, pGeometry, pProperties);
    }

    //************************************************************************************
    //************************************************************************************

    Condition::Pointer MPMLagrangianMovingParticle::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return Kratos::make_intrusive<MPMLagrangianMovingParticle>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
    }

    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************

    MPMLagrangianMovingParticle::~MPMLagrangianMovingParticle()
    {
    }
    

 
    double MPMLagrangianMovingParticle::GetPointLoadIntegrationWeight()
    {
        return 1.0;
    }
    void MPMLagrangianMovingParticle::CalculateAll(
        MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag
        )
    {
        KRATOS_TRY
            const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int block_size = this->GetBlockSize();

    // Resizing as needed the LHS
    const unsigned int matrix_size = number_of_nodes * block_size;

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != matrix_size )
        {
            rLeftHandSideMatrix.resize( matrix_size, matrix_size, false );
        }

        noalias( rLeftHandSideMatrix ) = ZeroMatrix(matrix_size,matrix_size); //resetting LHS
    }

    // Resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size( ) != matrix_size )
        {
            rRightHandSideVector.resize( matrix_size, false );
        }

        noalias( rRightHandSideVector ) = ZeroVector( matrix_size ); //resetting RHS
    }

        KRATOS_CATCH( "" )
    }

    void MPMLagrangianMovingParticle::FinalizeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
    {
        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        GeneralVariables Variables;

        Variables.CurrentDisp = CalculateCurrentDisp(Variables.CurrentDisp, rCurrentProcessInfo);

        array_1d<double,3> delta_xg = ZeroVector(3);
        array_1d<double, 3 > MPC_velocity = ZeroVector(3);
        array_1d<double, 3 > MPC_acceleration = ZeroVector(3);
        
        const double delta_time = rCurrentProcessInfo[DELTA_TIME];

        MPMShapeFunctionPointValues(Variables.N);

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            if (Variables.N[i] > std::numeric_limits<double>::epsilon() )
            {
                auto r_geometry = GetGeometry();
                array_1d<double, 3 > nodal_velocity = ZeroVector(3);
                array_1d<double, 3 > nodal_acceleration = ZeroVector(3);
                if (r_geometry[i].SolutionStepsDataHas(VELOCITY))
                    nodal_velocity = r_geometry[i].FastGetSolutionStepValue(VELOCITY);
                for ( unsigned int j = 0; j < dimension; j++ )
                {
                    delta_xg[j] += Variables.N[i] * Variables.CurrentDisp(i,j);
                    MPC_velocity[j] += Variables.N[i] * nodal_velocity[j];
                    MPC_acceleration[j] += Variables.N[i] * nodal_acceleration[j];
                }
            }
        }
        
    
        // Update the Material Point Condition Position
        KRATOS_WATCH(m_xg)
        m_xg += delta_xg ;
        // m_velocity = MPC_PreviousVelocity + 0.5 * delta_time * (MPC_acceleration + MPC_PreviousAcceleration);
        m_velocity = MPC_velocity;
        m_acceleration = MPC_acceleration;
        // total displacement of boundary particle
        m_displacement += delta_xg;
        std::cout.precision(20); 

        std::cout<<"displacement"<<m_displacement<<std::endl;
        
    }

    void MPMLagrangianMovingParticle::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
        std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo)
    {
        if (rValues.size() != 1)
            rValues.resize(1);

        if (rVariable == MP_COORD || rVariable == MPC_COORD) {
            rValues[0] = m_xg;
        }
        else if (rVariable == MPC_VELOCITY) {
            rValues[0] = m_velocity;
        }
        else if (rVariable == MPC_DISPLACEMENT) {
            rValues[0] = m_displacement;
        }
        else if (rVariable == MPC_CONTACT_FORCE) {
        rValues[0] = m_contact_force;
    }
        else {
            MPMParticleBaseLoadCondition::CalculateOnIntegrationPoints(
                rVariable, rValues, rCurrentProcessInfo);
        }
    }
    void MPMLagrangianMovingParticle::SetValuesOnIntegrationPoints(
        const Variable<array_1d<double, 3 > >& rVariable,
        const std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR_IF(rValues.size() > 1)
            << "Only 1 value per integration point allowed! Passed values vector size: "
            << rValues.size() << std::endl;

        if (rVariable == POINT_LOAD) {
            m_point_load = rValues[0];
        }
        else if (rVariable == MP_COORD || rVariable == MPC_COORD) {
            m_xg = rValues[0];
        }
        else if (rVariable == MPC_VELOCITY) {
            m_velocity = rValues[0];
        }
        else if (rVariable == MPC_DISPLACEMENT) {
            m_displacement = rValues[0];
        }
        else {
            MPMParticleBaseLoadCondition::SetValuesOnIntegrationPoints(
                rVariable, rValues, rCurrentProcessInfo);
        }
    }
} // Namespace Kratos