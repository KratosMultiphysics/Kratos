// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/base_solid_element.h"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
    void BaseSolidElement::Initialize()
    {
        KRATOS_TRY

        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints(  );

        //Constitutive Law initialisation

        if ( mConstitutiveLawVector.size() != IntegrationPoints.size() )
        {
            mConstitutiveLawVector.resize( IntegrationPoints.size() );
        }

        InitializeMaterial();

        KRATOS_CATCH( "" )
    }
    
    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            mConstitutiveLawVector[i]->InitializeSolutionStep( GetProperties(),
                    GetGeometry(), row( GetGeometry().ShapeFunctionsValues(  ), i ),
                    CurrentProcessInfo );
    }
    
    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::InitializeNonLinearIteration( ProcessInfo& CurrentProcessInfo )
    {
        // TODO: Add somethig if necessary
    }
    
    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::FinalizeNonLinearIteration( ProcessInfo& CurrentProcessInfo )
    {
        // TODO: Add somethig if necessary
    }

    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        //         std::cout << "in TL: calling FinalizeSolutionStep" << std::endl;
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            mConstitutiveLawVector[i]->FinalizeSolutionStep( GetProperties(),
                    GetGeometry(),
                    row( GetGeometry().ShapeFunctionsValues(  ), i ),
                    CurrentProcessInfo );
    }

    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::InitializeMaterial()
    {
        KRATOS_TRY

        if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            {
                mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
                mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(),
                        row( GetGeometry().ShapeFunctionsValues(  ), i ) );
            }
        }
        else
        {
            KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;
        }
        
        KRATOS_CATCH( "" );
    }

    //************************************************************************************
    //************************************************************************************
    
    void BaseSolidElement::ResetConstitutiveLaw()
    {
        KRATOS_TRY

        if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
                mConstitutiveLawVector[i]->ResetMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues(  ), i ) );
        }

        KRATOS_CATCH( "" )
    }
    
    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY;
        
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        
        if (rResult.size() != dimension * number_of_nodes)
        {
            rResult.resize(dimension * number_of_nodes,false);
        }

        const unsigned int pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        if(dimension == 2)
        {
            for (unsigned int i = 0; i < number_of_nodes; ++i)
            {
                const unsigned int index = i * 2;
                rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos).EquationId();
                rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos+1).EquationId();
            }
        }
        else
        {
            for (unsigned int i = 0; i < number_of_nodes; ++i)
            {
                const unsigned int index = i * 3;
                rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos).EquationId();
                rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos+1).EquationId();
                rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z,pos+2).EquationId();
            }
        }

        KRATOS_CATCH("")
    };
        
    //************************************************************************************
    //************************************************************************************
        
    void BaseSolidElement::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo 
        )
    {
        KRATOS_TRY;
        
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        rElementalDofList.resize(0);
        rElementalDofList.reserve(dimension*number_of_nodes);

        if(dimension == 2)
        {
            for (unsigned int i = 0; i < number_of_nodes; ++i)
            {
                rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
                rElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            }
        }
        else
        {
            for (unsigned int i = 0; i < number_of_nodes; ++i)
            {
                rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
                rElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
                rElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
            }
        }

        KRATOS_CATCH("")
    };
        
    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::GetValuesVector(
        Vector& rValues,
        int Step 
        )
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int mat_size = number_of_nodes * dimension;
        if (rValues.size() != mat_size)
        {
            rValues.resize(mat_size, false);
        }
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            const array_1d<double, 3 >& displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            const unsigned int index = i * dimension;
            for(unsigned int k = 0; k < dimension; ++k)
            {
                rValues[index + k] = displacement[k];
            }
        }
    }

    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::GetFirstDerivativesVector(
        Vector& rValues,
        int Step 
        )
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int mat_size = number_of_nodes * dimension;
        if (rValues.size() != mat_size)
        {
            rValues.resize(mat_size, false);
        }
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            const array_1d<double, 3 >& velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
            const unsigned int index = i * dimension;
            for(unsigned int k = 0; k < dimension; ++k)
            {
                rValues[index + k] = velocity[k];
            }
        }
    }
    
    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::GetSecondDerivativesVector(
        Vector& rValues,
        int Step 
        )
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int mat_size = number_of_nodes * dimension;
        if (rValues.size() != mat_size)
        {
            rValues.resize(mat_size, false);
        }
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            const array_1d<double, 3 >& acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
            const unsigned int index = i * dimension;
            for(unsigned int k = 0; k < dimension; ++k)
            {
                rValues[index + k] = acceleration[k];
            }
        }
    }
    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::CalculateRightHandSide( 
        VectorType& rRightHandSideVector, 
        ProcessInfo& rCurrentProcessInfo 
        )
    {
        // Calculation flags
        const bool CalculateStiffnessMatrixFlag = false;
        const bool CalculateResidualVectorFlag = true;
        MatrixType temp = Matrix();

        CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
    }

    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::CalculateLocalSystem( 
        MatrixType& rLeftHandSideMatrix, 
        VectorType& rRightHandSideVector, 
        ProcessInfo& rCurrentProcessInfo 
        )
    {
        //calculation flags
        const bool CalculateStiffnessMatrixFlag = true;
        const bool CalculateResidualVectorFlag = true;

        CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
    }
    
    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo 
        )
    {
        KRATOS_TRY;

        // Lumped
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int mat_size = dimension * number_of_nodes;

        if ( rMassMatrix.size1() != mat_size )
        {
            rMassMatrix.resize( mat_size, mat_size, false );
        }

        rMassMatrix = ZeroMatrix( mat_size, mat_size );
        
        Matrix DN_DX( number_of_nodes, dimension );
        Matrix J0(dimension,dimension), InvJ0(dimension,dimension);
        
        // Reading integration points and local gradients
        IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(GetGeometry());
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( integration_method );
        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(integration_method);
        
        const double density = GetProperties()[DENSITY];
        double thickness = 1.0;
        if ( dimension == 2 && GetProperties().Has( THICKNESS )) 
        {
            thickness = GetProperties()[THICKNESS];
        }
        
        for ( unsigned int point_number = 0; point_number < integration_points.size(); point_number++ )
        {   
            const double detJ0 = CalculateDerivativesOnReference(J0, InvJ0, DN_DX, point_number, integration_method);
            const double IntegrationWeight = GetIntegrationWeight(integration_points, point_number, detJ0) * thickness;
            const Vector& N = row(Ncontainer,point_number);
            
            for ( unsigned int i = 0; i < number_of_nodes; i++ )
            {
                const unsigned int index_i = i * dimension;
                
                for ( unsigned int j = 0; j < number_of_nodes; ++j )
                {
                    const unsigned int index_j = j * dimension;
                    const double NiNj_weight = N[i] * N[j] * IntegrationWeight * density;
                    
                    for ( unsigned int k = 0; k < dimension; k++ )
                    {
                        rMassMatrix( index_i + k, index_j + k ) += NiNj_weight;
                    }
                }
            }
        }
        
        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY;
        
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        //resizing as needed the LHS
        unsigned int mat_size = number_of_nodes * dimension;

        if ( rDampingMatrix.size1() != mat_size )
        {
            rDampingMatrix.resize( mat_size, mat_size, false );
        }

        noalias( rDampingMatrix ) = ZeroMatrix( mat_size, mat_size );

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::CalculateOnIntegrationPoints( 
        const Variable<double>& rVariable, 
        std::vector<double>& rOutput, 
        const ProcessInfo& rCurrentProcessInfo 
        )
    {
        if ( rOutput.size() != GetGeometry().IntegrationPoints(  ).size() )
        {
            rOutput.resize( GetGeometry().IntegrationPoints(  ).size() );
        }

        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
        {
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rOutput[ii] );
        }
    }

    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::CalculateOnIntegrationPoints( 
        const Variable<Vector>& rVariable, 
        std::vector<Vector>& rOutput, 
        const ProcessInfo& rCurrentProcessInfo 
        )
    {
        const unsigned int StrainSize = GetGeometry().WorkingSpaceDimension() == 2 ? 3 : 6;

        Vector StrainVector( StrainSize );

        if ( rVariable == INSITU_STRESS )
        {
            for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
            {
                if ( rOutput[ii].size() != StrainVector.size() )
                {
                    rOutput[ii].resize( StrainVector.size(), false );
                }

                rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( INSITU_STRESS, rOutput[ii] );
            }
        }
        else
        {
            if ( rOutput.size() != GetGeometry().IntegrationPoints(  ).size() )
            {
                rOutput.resize( GetGeometry().IntegrationPoints(  ).size() );
            }

            for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
            {
                rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rOutput[ii] );
            }
        }
    }

    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::CalculateOnIntegrationPoints( 
        const Variable<Matrix >& rVariable, 
        std::vector< Matrix >& rOutput, 
        const ProcessInfo& rCurrentProcessInfo 
        )
    {
        KRATOS_ERROR << "You are calle to the CalculateOnIntegrationPoints (Matrix) from the base class for solid elements" << std::endl; 
    }

    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::SetValueOnIntegrationPoints( 
        const Variable<double>& rVariable, 
        std::vector<double>& rValues, 
        const ProcessInfo& rCurrentProcessInfo 
        )
    {
        for ( unsigned int point_number = 0; point_number < GetGeometry().IntegrationPoints(  ).size(); point_number++ )
        {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,
                    rValues[point_number], rCurrentProcessInfo );
        }

    }
    
    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::SetValueOnIntegrationPoints( 
        const Variable<Vector>& rVariable, 
        std::vector<Vector>& rValues, 
        const ProcessInfo& rCurrentProcessInfo 
        )
    {

        for ( unsigned int point_number = 0; point_number < GetGeometry().IntegrationPoints(  ).size(); point_number++ )
        {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,
                    rValues[point_number], rCurrentProcessInfo );
        }
    }

    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::SetValueOnIntegrationPoints( 
        const Variable<Matrix>& rVariable, 
        std::vector<Matrix>& rValues, 
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        for ( unsigned int point_number = 0; point_number < GetGeometry().IntegrationPoints(  ).size(); point_number++ )
        {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,
                    rValues[point_number], rCurrentProcessInfo );
        }

    }

    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::GetValueOnIntegrationPoints( 
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo 
        )
    {
        if ( rValues.size() != GetGeometry().IntegrationPoints(  ).size() )
        {
            rValues.resize( GetGeometry().IntegrationPoints(  ).size(), false );
        }

        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
        {
            rValues[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rValues[ii] );
        }
    }


    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::GetValueOnIntegrationPoints( 
        const Variable<Vector>& rVariable, 
        std::vector<Vector>& rValues, 
        const ProcessInfo& rCurrentProcessInfo 
        )
    {
        const unsigned int& size = GetGeometry().IntegrationPoints(  ).size();

        if ( rValues.size() != size )
        {
            rValues.resize( size );
        }

        //TODO: decide which is the correct one
//         if ( rVariable == INSITU_STRESS || rVariable == PRESTRESS )
//         {
//             for ( unsigned int point_number = 0;
//                     point_number < GetGeometry().IntegrationPoints(  ).size();
//                     point_number++ )
//             {
//                 rValues[point_number] =
//                     mConstitutiveLawVector[point_number]->GetValue( PRESTRESS, rValues[point_number] );
//             }
//         }
//         if ( rVariable == PLASTIC_STRAIN_VECTOR )
//         {
//             for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
//             {
//                 if ( rValues[i].size() != 6 )
//                     rValues[i].resize( 6 );
//                 noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( PLASTIC_STRAIN_VECTOR, rValues[i] );
//             }
//         }
//        
        if ( rVariable == STRESSES )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            {
                if ( rValues[i].size() != 6 )
                {
                    rValues[i].resize( 6, false );
                }
                noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( STRESSES, rValues[i] );
            }
        }
        if ( rVariable == MATERIAL_PARAMETERS )
        {
            for ( unsigned int point_number = 0; point_number < GetGeometry().IntegrationPoints(  ).size(); point_number++ )
            {
                rValues[point_number] =
                    mConstitutiveLawVector[point_number]->GetValue( MATERIAL_PARAMETERS, rValues[point_number] );
            }
        }

        if ( rVariable == INTERNAL_VARIABLES )
        {
            for ( unsigned int point_number = 0; point_number < GetGeometry().IntegrationPoints(  ).size(); point_number++ )
            {
                rValues[point_number] =
                    mConstitutiveLawVector[point_number]->GetValue( INTERNAL_VARIABLES, rValues[point_number] );
            }
        }
    }

    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::GetValueOnIntegrationPoints( 
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues, 
        const ProcessInfo& rCurrentProcessInfo 
        )
    {
        if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
        {
            CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
        }

        if ( rVariable == PK2_STRESS_TENSOR )
        {
            CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
        }

//         if ( rVariable == GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR )
//         {
//             CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
//         }

    }

    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::Calculate( 
        const Variable<double>& rVariable, 
        double& Output, 
        const ProcessInfo& rCurrentProcessInfo 
        )
    {

        //HERE IT IS ESTIMATING DELTA_TIME for explicit elements
//         double lamda = 1.00; // parametro que depende del tipo de problema y del elemento pag 308 libro dinamica de Barbat
//         double c1 = 0.00; //sqrt(GetProperties()[YOUNG_MODULUS]/GetProperties()[DENSITY]); velocidad del sonido en el medio
//         double c2 = 0.00; // norma de la velocidad actual dentro del elemento
//         double c = 0.00;
//         double wmax = 0.00;
//         Vector Values( GetGeometry().IntegrationPoints(  ).size() );
//         Vector Velocities;
// 
//         GetFirstDerivativesVector( Velocities, 0 );
// 
//         if ( rVariable == DELTA_TIME )
//         {
//             for ( unsigned int point_number = 0;
//                     point_number < GetGeometry().IntegrationPoints(  ).size();
//                     point_number++ )
//             {
//                 mConstitutiveLawVector[point_number]-> GetValue( DELTA_TIME, c1 );
//                 Values[point_number] = c1;
//             }
//         }
// 
//         c1 = ( *std::max_element( Values.begin(), Values.end() ) );
// 
//         c2 = norm_2( Velocities );
// 
//         c = ( c1 > c2 ) ? c1 : c2;
// 
// 
//         double le = GetGeometry().Length();
//         //KRATOS_WATCH(le)
// 
//         /// maxima frecuencia de un elemento
//         wmax = ( lamda * c ) / le;
//         Output = 2.0 / wmax;
//         //KRATOS_WATCH(Output)

    }
    
    //************************************************************************************
    //************************************************************************************

    int  BaseSolidElement::Check( const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY;
        
        if ( DISPLACEMENT.Key() == 0 )
        {
            KRATOS_ERROR <<  "DISPLACEMENT has Key zero! (check if the application is correctly registered" << std::endl;
        }

        //verify that the dofs exist
        for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
        {
            if ( this->GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
            {
                KRATOS_ERROR << "Missing variable DISPLACEMENT on node " << this->GetGeometry()[i].Id() << std::endl;
            }

            if ( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false ||
                 this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false ||
                 this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false )
            {
                KRATOS_ERROR << "Missing one of the dofs for the variable DISPLACEMENT on node " << GetGeometry()[i].Id() << " of condition " << Id() << std::endl;
            }
        }
        
        return 0;

        KRATOS_CATCH( "" );
    }
    
    //************************************************************************************
    //************************************************************************************

    void BaseSolidElement::CalculateAll(
        MatrixType& rLeftHandSideMatrix, 
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
        )
    {
       KRATOS_ERROR << "You are calle to the CalculateAll from the base class for solid elements" << std::endl; 
    }
    
    //***********************************************************************
    //***********************************************************************

    double BaseSolidElement::GetIntegrationWeight(
        const GeometryType::IntegrationPointsArrayType& ThisIntegrationMethod,
        const unsigned int PointNumber,
        const double detJ
        )
    {
        return ThisIntegrationMethod[PointNumber].Weight() * detJ;
    }
    
    //************************************************************************************
    //************************************************************************************
    
    double BaseSolidElement::CalculateDerivativesOnReference(
        Matrix& J0, 
        Matrix& InvJ0, 
        Matrix& DN_DX, 
        const unsigned int PointNumber,
        IntegrationMethod ThisIntegrationMethod
        )
    {
        J0.clear();
        
        double detJ0;
        
        const Matrix& DN_De = GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod)[PointNumber];
        
        for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
        {
            const array_1d<double, 3>& coords = GetGeometry()[i].GetInitialPosition(); //NOTE: here we refer to the original, undeformed position!!
            for(unsigned int k = 0; k < GetGeometry().WorkingSpaceDimension(); k++)
            {
                for(unsigned int m = 0; m < GetGeometry().LocalSpaceDimension(); m++)
                {
                    J0(k,m) += coords[k]*DN_De(i,m);
                }
            }
        }
        
        MathUtils<double>::InvertMatrix( J0, InvJ0, detJ0 );
        
        noalias( DN_DX ) = prod( DN_De, InvJ0);
        
        return detJ0;
    }

} // Namespace Kratos


