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
#include "custom_conditions/particle_based_conditions/mpm_particle_point_load_condition.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************

    MPMParticlePointLoadCondition::MPMParticlePointLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry )
        : MPMParticleBaseLoadCondition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    MPMParticlePointLoadCondition::MPMParticlePointLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : MPMParticleBaseLoadCondition( NewId, pGeometry, pProperties )
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************

    Condition::Pointer MPMParticlePointLoadCondition::Create(IndexType NewId,GeometryType::Pointer pGeometry,PropertiesType::Pointer pProperties) const
    {
        return Kratos::make_intrusive<MPMParticlePointLoadCondition>(NewId, pGeometry, pProperties);
    }

    //************************************************************************************
    //************************************************************************************

    Condition::Pointer MPMParticlePointLoadCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return Kratos::make_intrusive<MPMParticlePointLoadCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
    }

    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************

    MPMParticlePointLoadCondition::~MPMParticlePointLoadCondition()
    {
    }

    //*************************COMPUTE FORCE AT EACH NODE*******************************
    //************************************************************************************
    /*
    This function distributes the pointload to the nodes
    */
    Matrix& MPMParticlePointLoadCondition::CalculateNodalForce(Matrix & rNodalForce, const ProcessInfo& rCurrentProcessInfo)
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        // Prepare variables
        GeneralVariables Variables;

        // Calculating shape function
        Variables.N = this->MPMShapeFunctionPointValues(Variables.N, m_xg);

        // Here MP contribution in terms of force are added
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            for (unsigned int j = 0; j < dimension; j++)
            {
                rNodalForce(j,i) = Variables.N[i] * m_point_load[j];
            }
        }

        return rNodalForce;
    }

    void MPMParticlePointLoadCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
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

        Matrix nodal_force = ZeroMatrix(3,number_of_nodes);

        nodal_force = CalculateNodalForce(nodal_force, rCurrentProcessInfo);

        for (unsigned int ii = 0; ii < number_of_nodes; ++ii)
        {
            const unsigned int base = ii*dimension;

            for(unsigned int k = 0; k < dimension; ++k)
            {
                rRightHandSideVector[base + k] += GetPointLoadIntegrationWeight() * nodal_force(k,ii);
            }
        }
        KRATOS_CATCH( "" )
    }
    //************************************************************************************
    //************************************************************************************

    double MPMParticlePointLoadCondition::GetPointLoadIntegrationWeight()
    {
        return 1.0;
    }

    void MPMParticlePointLoadCondition::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
    {
        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        GeneralVariables Variables;

        Variables.CurrentDisp = CalculateCurrentDisp(Variables.CurrentDisp, rCurrentProcessInfo);

        array_1d<double,3> delta_xg = ZeroVector(3);

        Variables.N = this->MPMShapeFunctionPointValues(Variables.N, m_xg);

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            if (Variables.N[i] > 1e-16 )
            {
                for ( unsigned int j = 0; j < dimension; j++ )
                {
                    delta_xg[j] += Variables.N[i] * Variables.CurrentDisp(i,j);
                }
            }
        }

        // Update the Material Point Condition Position
        m_xg += delta_xg ;
    }

    void MPMParticlePointLoadCondition::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
        std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo)
    {
        if (rValues.size() != 1)
            rValues.resize(1);

        if (rVariable == POINT_LOAD) {
            rValues[0] = m_point_load;
        }
        else {
            MPMParticleBaseLoadCondition::CalculateOnIntegrationPoints(
                rVariable, rValues, rCurrentProcessInfo);
        }
    }
    void MPMParticlePointLoadCondition::SetValuesOnIntegrationPoints(
        const Variable<array_1d<double, 3 > >& rVariable,
        std::vector<array_1d<double, 3 > > rValues,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR_IF(rValues.size() > 1)
            << "Only 1 value per integration point allowed! Passed values vector size: "
            << rValues.size() << std::endl;

        if (rVariable == POINT_LOAD) {
            m_point_load = rValues[0];
        }
        else {
            MPMParticleBaseLoadCondition::SetValuesOnIntegrationPoints(
                rVariable, rValues, rCurrentProcessInfo);
        }
    }
} // Namespace Kratos