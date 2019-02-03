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
#include "includes/define.h"
#include "custom_conditions/particle_based_conditions/mpm_particle_penalty_dirichlet_condition.h"
#include "utilities/math_utils.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/particle_mechanics_math_utilities.h"

namespace Kratos
{
//******************************* CONSTRUCTOR ****************************************
//************************************************************************************

MPMParticlePenaltyDirichletCondition::MPMParticlePenaltyDirichletCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : MPMParticleBaseDirichletCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

MPMParticlePenaltyDirichletCondition::MPMParticlePenaltyDirichletCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : MPMParticleBaseDirichletCondition( NewId, pGeometry, pProperties )
{
}

//********************************* CREATE *******************************************
//************************************************************************************

Condition::Pointer MPMParticlePenaltyDirichletCondition::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
{
    return Kratos::make_shared<MPMParticlePenaltyDirichletCondition>(NewId, pGeom, pProperties);
}

//************************************************************************************
//************************************************************************************

Condition::Pointer MPMParticlePenaltyDirichletCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_shared<MPMParticlePenaltyDirichletCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************

MPMParticlePenaltyDirichletCondition::~MPMParticlePenaltyDirichletCondition()
{
}

//************************************************************************************
//************************************************************************************

void MPMParticlePenaltyDirichletCondition::CalculateAll(
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

    // Get imposed displacement and normal vector
    const array_1d<double, 3 > & xg_c = this->GetValue(MPC_COORD);
    const array_1d<double, 3 > & imposed_displacement = this->GetValue (MPC_DISPLACEMENT);

    // Prepare variables
    GeneralVariables Variables;
    const double penalty_factor = this->GetValue(PENALTY_FACTOR);

    // Calculating shape function
    Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg_c);
    Variables.CurrentDisp = this->CalculateCurrentDisp(Variables.CurrentDisp, rCurrentProcessInfo);

    // Check contact
    // NOTE: the unit_normal_vector is assumed always pointing outside the boundary
    const array_1d<double, 3 > & unit_normal_vector = this->GetValue(MPC_NORMAL);
    array_1d<double, 3 > field_displacement = ZeroVector(3);
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
        if (Variables.N[i] > std::numeric_limits<double>::epsilon() )
            for ( unsigned int j = 0; j < dimension; j++ )
                field_displacement[j] += Variables.N[i] * Variables.CurrentDisp(i,j);

    // If penetrates, apply constraint
    const double penetration = MathUtils<double>::Dot((field_displacement - imposed_displacement), unit_normal_vector);
    if (penetration < 0.0)
    {
        // Arrange shape function
        Matrix shape_function = ZeroMatrix(block_size, matrix_size);
        for (unsigned int i = 0; i < number_of_nodes; i++)
            if (Variables.N[i] > std::numeric_limits<double>::epsilon())
                for (unsigned int j = 0; j < dimension; j++)
                    shape_function(j, block_size * i + j) = Variables.N[i];

        // if (!Is(SLIP)){
        //     // Slip condition (only normal direction imposition)
        //     // Prepare Rotation Matrix
        //     Matrix local_rotation_matrix = IdentityMatrix(block_size);
        //     this->GetRotationMatrix(local_rotation_matrix, unit_normal_vector);

        //     // Calculate slip shape function
        //     for (unsigned int i = 0; i < number_of_nodes; i++){
        //         const Matrix aux_matrix = Variables.N[i] * trans(local_rotation_matrix);
        //         for (unsigned int j = 0; j < dimension; j++)
        //             shape_function(i, block_size * i + j) = aux_matrix(0, j); // The first row is the normal direction
        //     }

        //     std::cout << "Identity check:" << prod(trans(local_rotation_matrix),local_rotation_matrix) << std::endl;
        //     KRATOS_ERROR << "TEST_END" << std::endl;

        // }

        // Calculate gap_function: nodal_displacement - imposed_displacement
        Vector gap_function = ZeroVector(matrix_size);
        for (unsigned int i = 0; i < number_of_nodes; i++)
            for ( unsigned int j = 0; j < dimension; j++ )
                gap_function[block_size * i + j] = (Variables.CurrentDisp(i,j) - imposed_displacement[j]);

        // Calculate LHS Matrix and RHS Vector
        noalias(rLeftHandSideMatrix)  += prod(trans(shape_function), shape_function);
        noalias(rRightHandSideVector) -= prod(prod(trans(shape_function), shape_function), gap_function);

        rLeftHandSideMatrix  *= penalty_factor * this->GetIntegrationWeight();
        rRightHandSideVector *= penalty_factor * this->GetIntegrationWeight();
    }

    KRATOS_CATCH( "" )
}

void MPMParticlePenaltyDirichletCondition::GetRotationMatrix(
    MatrixType& rRotationMatrix, const VectorType& rNormalVector)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int block_size = this->GetBlockSize();

    if(dimension == 2){
        noalias(rRotationMatrix) = IdentityMatrix(block_size);

        double aux = rNormalVector[0]*rNormalVector[0] + rNormalVector[1]*rNormalVector[1];
        aux = sqrt(aux);
        if (std::abs(aux) < std::numeric_limits<double>::epsilon()) aux = std::numeric_limits<double>::epsilon();

        rRotationMatrix(0,0) =  rNormalVector[0]/aux;
        rRotationMatrix(0,1) =  rNormalVector[1]/aux;
        rRotationMatrix(1,0) = -rNormalVector[1]/aux;
        rRotationMatrix(1,1) =  rNormalVector[0]/aux;
    }
    else if (dimension == 3){
        noalias(rRotationMatrix) = IdentityMatrix(block_size);

        double aux = rNormalVector[0]*rNormalVector[0] + rNormalVector[1]*rNormalVector[1] + rNormalVector[2]*rNormalVector[2];
        aux = sqrt(aux);
        if (std::abs(aux) < std::numeric_limits<double>::epsilon()) aux = std::numeric_limits<double>::epsilon();

        rRotationMatrix(0,0) = rNormalVector[0]/aux;
        rRotationMatrix(0,1) = rNormalVector[1]/aux;
        rRotationMatrix(0,2) = rNormalVector[2]/aux;

        // Define the new coordinate system, where the first vector is aligned with the normal
        // To choose the remaining two vectors, we project the first component of the cartesian base to the tangent plane
        Vector rT1 = ZeroVector(3);
        rT1(0) = 1.0;
        rT1(1) = 0.0;
        rT1(2) = 0.0;
        double dot = rRotationMatrix(0,0); //this->Dot(rN,rT1);

        // It is possible that the normal is aligned with (1,0,0), resulting in norm(rT1) = 0
        // If this is the case, repeat the procedure using (0,1,0)
        if ( fabs(dot) > 0.99 )
        {
            rT1(0) = 0.0;
            rT1(1) = 1.0;
            rT1(2) = 0.0;

            dot = rRotationMatrix(0,1); //this->Dot(rN,rT1);
        }

        // Calculate projection and normalize
        rT1[0] -= dot*rRotationMatrix(0,0);
        rT1[1] -= dot*rRotationMatrix(0,1);
        rT1[2] -= dot*rRotationMatrix(0,2);
        ParticleMechanicsMathUtilities<double>::Normalize(rT1);

        rRotationMatrix(1,0) = rT1[0];
        rRotationMatrix(1,0) = rT1[1];
        rRotationMatrix(1,0) = rT1[2];

        // The third base component is choosen as N x T1, which is normalized by construction
        rRotationMatrix(2,0) = rRotationMatrix(0,1)*rT1[2] - rRotationMatrix(0,2)*rT1[1];
        rRotationMatrix(2,1) = rRotationMatrix(0,2)*rT1[0] - rRotationMatrix(0,0)*rT1[2];
        rRotationMatrix(2,2) = rRotationMatrix(0,0)*rT1[1] - rRotationMatrix(0,1)*rT1[0];
    }
    else{
        KRATOS_ERROR <<  "Dimension given is wrong: Something is wrong with the given dimension in function: GetRotationMatrix" << std::endl;
    }
}


} // Namespace Kratos


