// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/line_load_condition.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************
    
    LineLoadCondition::LineLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry )
        : BaseLoadCondition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************
    
    LineLoadCondition::LineLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : BaseLoadCondition( NewId, pGeometry, pProperties )
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************
    
    Condition::Pointer LineLoadCondition::Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const
    {
        return boost::make_shared<LineLoadCondition>(NewId, pGeom, pProperties);
    }

    //************************************************************************************
    //************************************************************************************
    
    Condition::Pointer LineLoadCondition::Create( 
        IndexType NewId, 
        NodesArrayType const& ThisNodes,  
        PropertiesType::Pointer pProperties 
        ) const
    {
        return boost::make_shared<LineLoadCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
    }

    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************
    
    LineLoadCondition::~LineLoadCondition()
    {
    }

    //************************************************************************************
    //************************************************************************************

    void LineLoadCondition::CalculateAll( 
        MatrixType& rLeftHandSideMatrix, 
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag 
        )
    {
        KRATOS_TRY;
        
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const bool bHasRotDof = this->HasRotDof();

        // Resizing as needed the LHS
        unsigned int mat_size = number_of_nodes * dimension;
        if (bHasRotDof) mat_size *= 2;

        if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
        {
            if ( rLeftHandSideMatrix.size1() != mat_size )
            {
                rLeftHandSideMatrix.resize( mat_size, mat_size, false );
            }

            noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
        }

        //resizing as needed the RHS
        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            if ( rRightHandSideVector.size( ) != mat_size )
            {
                rRightHandSideVector.resize( mat_size, false );
            }

            noalias( rRightHandSideVector ) = ZeroVector( mat_size ); //resetting RHS
        }

        IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(GetGeometry());

        // Reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(integration_method);

        const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients(integration_method);

        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(integration_method);

        // Sizing work matrices
        Vector pressure_on_nodes = ZeroVector( number_of_nodes );

        // Pressure applied to the element itself
        double pressure_on_condition = 0.0;
        if( this->Has( PRESSURE ) )
        {
            pressure_on_condition += this->GetValue( PRESSURE );
        }
        if( this->Has( NEGATIVE_FACE_PRESSURE ) )
        {
            pressure_on_condition += this->GetValue( NEGATIVE_FACE_PRESSURE );
        }
        if( this->Has( POSITIVE_FACE_PRESSURE ) )
        {
            pressure_on_condition -= this->GetValue( POSITIVE_FACE_PRESSURE );
        }

        for ( unsigned int i = 0; i < pressure_on_nodes.size(); i++ )
        {
            pressure_on_nodes[i] = pressure_on_condition;
            if( GetGeometry()[i].SolutionStepsDataHas( NEGATIVE_FACE_PRESSURE) )
            {
                pressure_on_nodes[i] += GetGeometry()[i].FastGetSolutionStepValue( NEGATIVE_FACE_PRESSURE );
            }
            if( GetGeometry()[i].SolutionStepsDataHas( POSITIVE_FACE_PRESSURE) )
            {
                pressure_on_nodes[i] -= GetGeometry()[i].FastGetSolutionStepValue( POSITIVE_FACE_PRESSURE );
            }
        }
        // Vector with a loading applied to the elemnt
        array_1d<double, 3 > line_load = ZeroVector(3);
        array_1d<double,3> gauss_load  = ZeroVector(3);
        if( this->Has( LINE_LOAD ) )
        {
            noalias(line_load) = this->GetValue( LINE_LOAD );
        }
        
        for ( unsigned int point_number = 0; point_number < integration_points.size(); point_number++ )
        {                   
            const double det_j = GetGeometry().DeterminantOfJacobian( integration_points[point_number] );

            const double integration_weight = GetIntegrationWeight(integration_points, point_number, det_j); 
            
            const array_1d<double, 3> normal = GetGeometry().Normal( integration_points[point_number] );
            
            // Calculating the pressure on the gauss point
            double gauss_pressure = 0.0;
            for ( unsigned int ii = 0; ii < number_of_nodes; ii++ )
            {
                gauss_pressure += Ncontainer( point_number, ii ) * pressure_on_nodes[ii];
            }

            if ( CalculateStiffnessMatrixFlag == true )
            {
                if ( gauss_pressure != 0.0 )
                {
                    CalculateAndSubKp( rLeftHandSideMatrix, DN_De[point_number], row( Ncontainer, point_number ), gauss_pressure, integration_weight );
                }
            }
            // Adding contributions to the residual vector
            if ( CalculateResidualVectorFlag == true )
            {
                if ( gauss_pressure != 0.0 )
                {
                    CalculateAndAddPressureForce( rRightHandSideVector, row( Ncontainer, point_number ), normal, gauss_pressure, integration_weight );
                }
            }

            gauss_load = line_load;
            for (unsigned int ii = 0; ii < number_of_nodes; ++ii)
            {
                if( GetGeometry()[ii].SolutionStepsDataHas( LINE_LOAD ) )
                {
                    noalias(gauss_load) += ( Ncontainer( point_number, ii )) * GetGeometry()[ii].FastGetSolutionStepValue( LINE_LOAD );
                }
            }

            for (unsigned int ii = 0; ii < number_of_nodes; ++ii)
            {
                unsigned int base = ii * dimension;
                if (bHasRotDof) base *= 2;
                for(unsigned int k = 0; k < dimension; ++k)
                {
                    rRightHandSideVector[base + k] += integration_weight * Ncontainer( point_number, ii ) * gauss_load[k];
                }
            }
        }

        if (bHasRotDof) this->CalculateAndAddWorkEquivalentNodalForcesLineLoad(gauss_load,rRightHandSideVector);

        KRATOS_CATCH( "" )
    }

    //***********************************************************************
    //***********************************************************************

    void LineLoadCondition::CalculateAndSubKp(
        Matrix& K,
        const Matrix& DN_De,
        const Vector& N,
        const double Pressure,
        const double IntegrationWeight
        )
    {
        KRATOS_TRY

        Matrix Kij( 2, 2 );
        Matrix Cross_gn( 2, 2 );

        //TODO: decide what to do with thickness
        //const double h0 = GetProperties()[THICKNESS];
        const double h0 = 1.00;
        Cross_gn( 0, 0 ) = 0.0;
        Cross_gn( 0, 1 ) = -h0;
        Cross_gn( 1, 0 ) = -h0;
        Cross_gn( 1, 1 ) = 0.0;

        for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
        {
            const unsigned int RowIndex = i * 2;

            for ( unsigned int j = 0; j < GetGeometry().size(); j++ )
            {
                const unsigned int ColIndex = j * 2;

                const double coeff = Pressure * N[i] * DN_De( j, 0 ) * IntegrationWeight;
                Kij = -coeff * Cross_gn;

                //TAKE CARE: the load correction matrix should be SUBTRACTED not added
                MathUtils<double>::SubtractMatrix( K, Kij, RowIndex, ColIndex );
            }
        }

        KRATOS_CATCH( "" )
    }

    //***********************************************************************
    //***********************************************************************

    void LineLoadCondition::CalculateAndAddPressureForce(
        Vector& rRightHandSideVector,
        const Vector& N,
        const array_1d<double, 3>& Normal,
        double Pressure,
        double IntegrationWeight 
        )
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = dimension * i;
            if(this->HasRotDof()) index *= 2;

            const double coeff = Pressure * N[i] * IntegrationWeight;
            
            rRightHandSideVector[index   ]  -= coeff * Normal[0];
            rRightHandSideVector[index + 1] -= coeff * Normal[1];
            if (dimension == 3)
            {
            rRightHandSideVector[index + 2] -= coeff * Normal[2];
            }
        }
    }

    //***********************************************************************
    //***********************************************************************

    void LineLoadCondition::CalculateAndAddWorkEquivalentNodalForcesLineLoad(
        const Vector ForceInput, VectorType& rRightHandSideVector)
        {
            KRATOS_TRY;
            const int dimension = this->GetGeometry().WorkingSpaceDimension();
            //calculate orthogonal load vector
            Vector GeometricOrientation = ZeroVector(dimension);
            GeometricOrientation[0] = this->GetGeometry()[1].X() 
                - this->GetGeometry()[0].X();
            GeometricOrientation[1] = this->GetGeometry()[1].Y() 
                - this->GetGeometry()[0].Y();
            if (dimension == 3)
            {
                GeometricOrientation[2] = this->GetGeometry()[1].Z() 
                    - this->GetGeometry()[0].Z();
            }
   
            const double VectorNormA = MathUtils<double>::Norm(GeometricOrientation);
            if (VectorNormA != 0.00) GeometricOrientation /= VectorNormA;
    
            Vector LineLoadDir = ZeroVector(dimension);
            for (int i = 0; i < dimension; ++i)
            {
                LineLoadDir[i] = ForceInput[i];
            }
    
            const double VectorNormB = MathUtils<double>::Norm(LineLoadDir);
            if (VectorNormB != 0.00) LineLoadDir /= VectorNormB;
    
            double cosAngle = 0.00;
            for (int i = 0; i < dimension; ++i)
            {
                cosAngle += LineLoadDir[i] * GeometricOrientation[i];
            }
    
            const double sinAngle = sqrt(1.00 - (cosAngle*cosAngle));
            const double NormForceVectorOrth = sinAngle * VectorNormB;
    
    
            Vector NodeA = ZeroVector(dimension);
            NodeA[0] = this->GetGeometry()[0].X();
            NodeA[1] = this->GetGeometry()[0].Y();
            if (dimension == 3)	NodeA[2] = this->GetGeometry()[0].Z();
    
            Vector NodeB = ZeroVector(dimension);
            NodeB = NodeA + LineLoadDir;
    
            Vector NodeC = ZeroVector(dimension);
            NodeC = NodeA + (GeometricOrientation*cosAngle);
    
            Vector LoadOrthogonalDir = ZeroVector(dimension);
            LoadOrthogonalDir = NodeB - NodeC;
            const double VectorNormC = MathUtils<double>::Norm(LoadOrthogonalDir);
            if(VectorNormC != 0.00) LoadOrthogonalDir /= VectorNormC;
    
    
    
            // now caluclate respective work equivilent nodal moments
    
            const double CustomMoment = NormForceVectorOrth *
                VectorNormA*VectorNormA / 12.00;
    
            Vector MomentNodeA = ZeroVector(dimension);
            MomentNodeA = MathUtils<double>::CrossProduct(GeometricOrientation,
                LoadOrthogonalDir);
            MomentNodeA *= CustomMoment;
    
            for (int i = 0; i < dimension; ++i)
            {
                rRightHandSideVector[(1 * dimension) + i] += MomentNodeA[i];
                rRightHandSideVector[(3 * dimension) + i] -= MomentNodeA[i];
            }
    
            KRATOS_CATCH("")            
        }


} // Namespace Kratos


