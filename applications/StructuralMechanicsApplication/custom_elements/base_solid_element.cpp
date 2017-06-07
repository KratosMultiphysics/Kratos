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

    void BaseSolidElement::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY;
        
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int Dimension = GetGeometry().WorkingSpaceDimension();
        
        if (rResult.size() != Dimension * NumberOfNodes)
        {
            rResult.resize(Dimension * NumberOfNodes,false);
        }

        const unsigned int pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        if(Dimension == 2)
        {
            for (unsigned int i = 0; i < NumberOfNodes; ++i)
            {
                const unsigned int index = i * 2;
                rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos).EquationId();
                rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos+1).EquationId();
            }
        }
        else
        {
            for (unsigned int i = 0; i < NumberOfNodes; ++i)
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
        
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int Dimension = GetGeometry().WorkingSpaceDimension();
        rElementalDofList.resize(0);
        rElementalDofList.reserve(Dimension*NumberOfNodes);

        if(Dimension == 2)
        {
            for (unsigned int i = 0; i < NumberOfNodes; ++i)
            {
                rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
                rElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            }
        }
        else
        {
            for (unsigned int i = 0; i < NumberOfNodes; ++i)
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
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int Dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int MatSize = NumberOfNodes * Dimension;
        if (rValues.size() != MatSize)
        {
            rValues.resize(MatSize, false);
        }
        for (unsigned int i = 0; i < NumberOfNodes; i++)
        {
            const array_1d<double, 3 >& Displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            const unsigned int index = i * Dimension;
            for(unsigned int k = 0; k < Dimension; ++k)
            {
                rValues[index + k] = Displacement[k];
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
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int Dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int MatSize = NumberOfNodes * Dimension;
        if (rValues.size() != MatSize)
        {
            rValues.resize(MatSize, false);
        }
        for (unsigned int i = 0; i < NumberOfNodes; i++)
        {
            const array_1d<double, 3 >& Velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
            const unsigned int index = i * Dimension;
            for(unsigned int k = 0; k < Dimension; ++k)
            {
                rValues[index + k] = Velocity[k];
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
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int Dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int MatSize = NumberOfNodes * Dimension;
        if (rValues.size() != MatSize)
        {
            rValues.resize(MatSize, false);
        }
        for (unsigned int i = 0; i < NumberOfNodes; i++)
        {
            const array_1d<double, 3 >& Acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
            const unsigned int index = i * Dimension;
            for(unsigned int k = 0; k < Dimension; ++k)
            {
                rValues[index + k] = Acceleration[k];
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
        unsigned int Dimension = GetGeometry().WorkingSpaceDimension();
        unsigned int NumberOfNodes = GetGeometry().size();
        unsigned int MatSize = Dimension * NumberOfNodes;

        if ( rMassMatrix.size1() != MatSize )
        {
            rMassMatrix.resize( MatSize, MatSize, false );
        }

        rMassMatrix = ZeroMatrix( MatSize, MatSize );
        
        Matrix DN_DX( NumberOfNodes, Dimension );
        Matrix J0(Dimension,Dimension), InvJ0(Dimension,Dimension);
        
        //reading integration points and local gradients
        IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(GetGeometry());
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( integration_method );
        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(integration_method);
        
        const double density = GetProperties()[DENSITY];
        double thickness = 1.0;
        if ( Dimension == 2 && GetProperties().Has( THICKNESS )) 
        {
            thickness = GetProperties()[THICKNESS];
        }

        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {

            const Matrix& DN_De = GetGeometry().ShapeFunctionsLocalGradients(integration_method)[PointNumber];

            double detJ0;
            CalculateDerivativesOnReference(J0, InvJ0, DN_DX, detJ0, DN_De);
            const double weight = integration_points[PointNumber].Weight()*detJ0*thickness;
            const Vector& N = row(Ncontainer,PointNumber);

            for ( unsigned int i = 0; i < NumberOfNodes; i++ )
            {
                const unsigned int index_i = i * Dimension;

                for ( unsigned int j = 0; j < NumberOfNodes; ++j )
                {
                    const unsigned int index_j = j*Dimension;
                    const double NiNj_weight = N[i]*N[j]*weight*density;

                    for ( unsigned int k = 0; k < Dimension; k++ )
                    {
                        rMassMatrix( index_i+k, index_j+k ) += NiNj_weight;
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
        
        unsigned int NumberOfNodes = GetGeometry().size();
        unsigned int Dimension = GetGeometry().WorkingSpaceDimension();

        //resizing as needed the LHS
        unsigned int MatSize = NumberOfNodes * Dimension;

        if ( rDampingMatrix.size1() != MatSize )
        {
            rDampingMatrix.resize( MatSize, MatSize, false );
        }

        noalias( rDampingMatrix ) = ZeroMatrix( MatSize, MatSize );

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
        for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints(  ).size(); PointNumber++ )
        {
            mConstitutiveLawVector[PointNumber]->SetValue( rVariable,
                    rValues[PointNumber], rCurrentProcessInfo );
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

        for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints(  ).size(); PointNumber++ )
        {
            mConstitutiveLawVector[PointNumber]->SetValue( rVariable,
                    rValues[PointNumber], rCurrentProcessInfo );
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
        for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints(  ).size(); PointNumber++ )
        {
            mConstitutiveLawVector[PointNumber]->SetValue( rVariable,
                    rValues[PointNumber], rCurrentProcessInfo );
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
//             for ( unsigned int PointNumber = 0;
//                     PointNumber < GetGeometry().IntegrationPoints(  ).size();
//                     PointNumber++ )
//             {
//                 rValues[PointNumber] =
//                     mConstitutiveLawVector[PointNumber]->GetValue( PRESTRESS, rValues[PointNumber] );
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
            for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints(  ).size(); PointNumber++ )
            {
                rValues[PointNumber] =
                    mConstitutiveLawVector[PointNumber]->GetValue( MATERIAL_PARAMETERS, rValues[PointNumber] );
            }
        }

        if ( rVariable == INTERNAL_VARIABLES )
        {
            for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints(  ).size(); PointNumber++ )
            {
                rValues[PointNumber] =
                    mConstitutiveLawVector[PointNumber]->GetValue( INTERNAL_VARIABLES, rValues[PointNumber] );
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
//             for ( unsigned int PointNumber = 0;
//                     PointNumber < GetGeometry().IntegrationPoints(  ).size();
//                     PointNumber++ )
//             {
//                 mConstitutiveLawVector[PointNumber]-> GetValue( DELTA_TIME, c1 );
//                 Values[PointNumber] = c1;
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

} // Namespace Kratos


