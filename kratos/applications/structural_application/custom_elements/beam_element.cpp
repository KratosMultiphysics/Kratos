
/* *********************************************************
*
*   Last Modified by:    $Author: nelson $
*   Date:                $Date: 2009-01-21 09:56:09 $
*   Revision:            $Revision: 1.2 $
*
* ***********************************************************/

// System includes

//custom_elements
// External includes

// Project includes

#include "includes/define.h"
#include "structural_application.h"
#include "custom_elements/beam_element.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "custom_utilities/sd_math_utils.h"



namespace Kratos
{
    //*****************************************************************************
    //*****************************************************************************

    typedef GeometryData::IntegrationMethod IntegrationMethod;

    BeamElement::BeamElement( IndexType NewId, GeometryType::Pointer pGeometry )
            : Element( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
        //THIS IS THE DEFAULT CONSTRUCTOR
    }

    BeamElement::BeamElement( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
            : Element( NewId, pGeometry, pProperties )
    {
        KRATOS_TRY

        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        unsigned int number_of_nodes = GetGeometry().size();


        if ( dimension != 3 )
        {
            std::cout << "This element works only with a 2 node line and 3D dimension" << std::endl;
            return;
        }

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            GetGeometry()[i].pAddDof( DISPLACEMENT_X, REACTION_X );
            GetGeometry()[i].pAddDof( DISPLACEMENT_Y, REACTION_Y );
            GetGeometry()[i].pAddDof( DISPLACEMENT_Z, REACTION_Z );
            GetGeometry()[i].pAddDof( ROTATION_X, MOMENT_X );
            GetGeometry()[i].pAddDof( ROTATION_Y, MOMENT_Y );
            GetGeometry()[i].pAddDof( ROTATION_Z, MOMENT_Z );
        }

        KRATOS_CATCH( "" )

    }


    Element::Pointer BeamElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new BeamElement( NewId, GetGeometry().Create( ThisNodes ),
                                 pProperties ) );
    }


    BeamElement::~BeamElement()
    {
    }

    //************************************************************************************
    //THIS IS THE INITIALIZATION OF THE ELEMENT (CALLED AT THE BEGIN OF EACH CALCULATION)
    //************************************************************************************


    void BeamElement::Initialize()
    {
        KRATOS_TRY
        CalculateSectionProperties();
        KRATOS_CATCH( "" )

    }

//************************************************************************************
//************************************************************************************

    void BeamElement::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {

    }

//************************************************************************************
//************************************************************************************

    void BeamElement::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        Matrix Rotation( 12, 12 );
        Matrix LocalMatrix;
        array_1d<double, 12 > CurrentDisplacement;
        array_1d<double, 12 > LocalDisplacement;
        array_1d<double, 12 > Results;

        CurrentDisplacement( 0 )  =   GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_X );
        CurrentDisplacement( 1 )  =   GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_Y );
        CurrentDisplacement( 2 )  =   GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_Z );
        CurrentDisplacement( 3 )  =   GetGeometry()[0].GetSolutionStepValue( ROTATION_X );
        CurrentDisplacement( 4 )  =   GetGeometry()[0].GetSolutionStepValue( ROTATION_Y );
        CurrentDisplacement( 5 )  =   GetGeometry()[0].GetSolutionStepValue( ROTATION_Z );
        CurrentDisplacement( 6 )  =   GetGeometry()[1].GetSolutionStepValue( DISPLACEMENT_X );
        CurrentDisplacement( 7 )  =   GetGeometry()[1].GetSolutionStepValue( DISPLACEMENT_Y );
        CurrentDisplacement( 8 )  =   GetGeometry()[1].GetSolutionStepValue( DISPLACEMENT_Z );
        CurrentDisplacement( 9 )  =   GetGeometry()[1].GetSolutionStepValue( ROTATION_X );
        CurrentDisplacement( 10 )         =   GetGeometry()[1].GetSolutionStepValue( ROTATION_Y );
        CurrentDisplacement( 11 )         =   GetGeometry()[1].GetSolutionStepValue( ROTATION_Z );

        CalculateTransformationMatrix( Rotation );
        CalculateLocalMatrix( LocalMatrix );
        noalias( LocalDisplacement ) = prod( trans( Rotation ), CurrentDisplacement );
        noalias( Results ) = prod( LocalMatrix, LocalDisplacement );

        GetGeometry()[0].SetValue( MOMENT_X, Results( 3 ) );
        GetGeometry()[0].SetValue( MOMENT_Y, Results( 4 ) );
        GetGeometry()[0].SetValue( MOMENT_Z, Results( 5 ) );
        GetGeometry()[1].SetValue( MOMENT_X, Results( 9 ) );
        GetGeometry()[1].SetValue( MOMENT_Y, Results( 10 ) );
        GetGeometry()[1].SetValue( MOMENT_Z, Results( 11 ) );

        return;
    }


//************************************************************************************
//************************************************************************************


    void BeamElement::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo,
                                    bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag )
    {
        KRATOS_TRY

        unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        if ( dimension != 3 )
        {
            std::cout << "this element works only with a 2 node line and 3D dimension" << std::endl;
            return;
        }


        if ( CalculateStiffnessMatrixFlag == true )
        {
            CalculateLHS( rLeftHandSideMatrix );
        }


        if ( CalculateResidualVectorFlag == true )
        {
            CalculateRHS( rRightHandSideVector );
        }

        KRATOS_CATCH( "" )
    }


//************************************************************************************
//************************************************************************************

    void  BeamElement::CalculateRightHandSide( VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = false;
        bool CalculateResidualVectorFlag = true;
        MatrixType temp = Matrix();

        CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag,  CalculateResidualVectorFlag );
    }

//************************************************************************************
//************************************************************************************


    void BeamElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                            VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;
        CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                      CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

    }

//************************************************************************************
//************************************************************************************

    void BeamElement::EquationIdVector( EquationIdVectorType& rResult,
                                        ProcessInfo& CurrentProcessInfo )
    {
        if ( rResult.size() != 12 )
            rResult.resize( 12, false );



        rResult[0] = GetGeometry()[0].GetDof( DISPLACEMENT_X ).EquationId();

        rResult[1] = GetGeometry()[0].GetDof( DISPLACEMENT_Y ).EquationId();

        rResult[2] = GetGeometry()[0].GetDof( DISPLACEMENT_Z ).EquationId();

        rResult[3] = GetGeometry()[0].GetDof( ROTATION_X ).EquationId();

        rResult[4] = GetGeometry()[0].GetDof( ROTATION_Y ).EquationId();

        rResult[5] = GetGeometry()[0].GetDof( ROTATION_Z ).EquationId();

        rResult[6] = GetGeometry()[1].GetDof( DISPLACEMENT_X ).EquationId();

        rResult[7] = GetGeometry()[1].GetDof( DISPLACEMENT_Y ).EquationId();

        rResult[8] = GetGeometry()[1].GetDof( DISPLACEMENT_Z ).EquationId();

        rResult[9] = GetGeometry()[1].GetDof( ROTATION_X ).EquationId();

        rResult[10] = GetGeometry()[1].GetDof( ROTATION_Y ).EquationId();

        rResult[11] = GetGeometry()[1].GetDof( ROTATION_Z ).EquationId();

    }


    //************************************************************************************
    //************************************************************************************

    void BeamElement::GetDofList( DofsVectorType& ElementalDofList, ProcessInfo&
                                  CurrentProcessInfo )
    {
        ElementalDofList.resize( 0 );

        ElementalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_X ) );
        ElementalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_Y ) );
        ElementalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_Z ) );
        ElementalDofList.push_back( GetGeometry()[0].pGetDof( ROTATION_X ) );
        ElementalDofList.push_back( GetGeometry()[0].pGetDof( ROTATION_Y ) );
        ElementalDofList.push_back( GetGeometry()[0].pGetDof( ROTATION_Z ) );
        ElementalDofList.push_back( GetGeometry()[1].pGetDof( DISPLACEMENT_X ) );
        ElementalDofList.push_back( GetGeometry()[1].pGetDof( DISPLACEMENT_Y ) );
        ElementalDofList.push_back( GetGeometry()[1].pGetDof( DISPLACEMENT_Z ) );
        ElementalDofList.push_back( GetGeometry()[1].pGetDof( ROTATION_X ) );
        ElementalDofList.push_back( GetGeometry()[1].pGetDof( ROTATION_Y ) );
        ElementalDofList.push_back( GetGeometry()[1].pGetDof( ROTATION_Z ) );

    }

    //************************************************************************************
    //************************************************************************************

    void BeamElement::GetValuesVector( Vector& values, int Step )
    {
        if ( values.size() != 12 )
            values.resize( 12, false );

        values( 0 ) = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_X, Step );

        values( 1 ) = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_Y, Step );

        values( 2 ) = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_Z, Step );

        values( 3 ) = GetGeometry()[0].GetSolutionStepValue( ROTATION_X, Step );

        values( 4 ) = GetGeometry()[0].GetSolutionStepValue( ROTATION_Y, Step );

        values( 5 ) = GetGeometry()[0].GetSolutionStepValue( ROTATION_Z, Step );

        values( 6 ) = GetGeometry()[1].GetSolutionStepValue( DISPLACEMENT_X, Step );

        values( 7 ) = GetGeometry()[1].GetSolutionStepValue( DISPLACEMENT_Y, Step );

        values( 8 ) = GetGeometry()[1].GetSolutionStepValue( DISPLACEMENT_Z, Step );

        values( 9 ) = GetGeometry()[1].GetSolutionStepValue( ROTATION_X, Step );

        values( 10 ) = GetGeometry()[1].GetSolutionStepValue( ROTATION_Y, Step );

        values( 11 ) = GetGeometry()[1].GetSolutionStepValue( ROTATION_Z, Step );


    }

//************************************************************************************
//************************************************************************************

    void BeamElement::CalculateLHS( Matrix& rLeftHandSideMatrix )
    {

        Matrix LocalMatrix;
        Matrix Rotation;
        Matrix aux_matrix;

        LocalMatrix.resize( 12, 12, false );
        Rotation.resize( 12, 12, false );
        aux_matrix.resize( 12, 12, false );
        rLeftHandSideMatrix.resize( 12, 12, false );


        CalculateLocalMatrix( LocalMatrix );
        CalculateTransformationMatrix( Rotation );
        noalias( aux_matrix ) = prod( Rotation, LocalMatrix );
        noalias( rLeftHandSideMatrix ) = prod( aux_matrix, Matrix( trans( Rotation ) ) );
        return;
    }


//************************************************************************************
//************************************************************************************
    void BeamElement::CalculateRHS( Vector& rRightHandSideVector )

    {

        Matrix Rotation;
        Matrix GlobalMatrix;
        Vector LocalBody;

        array_1d<double, 12 > CurrentDisplacement;

        Rotation.resize( 12, 12, false );
        rRightHandSideVector = ZeroVector( 12 );
        LocalBody = ZeroVector( 12 );

        CalculateTransformationMatrix( Rotation );
        CalculateBodyForce( Rotation, LocalBody, rRightHandSideVector );


        CurrentDisplacement( 0 )  =   GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_X );
        CurrentDisplacement( 1 )  =   GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_Y );
        CurrentDisplacement( 2 )  =   GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_Z );
        CurrentDisplacement( 3 )  =   GetGeometry()[0].GetSolutionStepValue( ROTATION_X );
        CurrentDisplacement( 4 )  =   GetGeometry()[0].GetSolutionStepValue( ROTATION_Y );
        CurrentDisplacement( 5 )  =   GetGeometry()[0].GetSolutionStepValue( ROTATION_Z );
        CurrentDisplacement( 6 )  =   GetGeometry()[1].GetSolutionStepValue( DISPLACEMENT_X );
        CurrentDisplacement( 7 )  =   GetGeometry()[1].GetSolutionStepValue( DISPLACEMENT_Y );
        CurrentDisplacement( 8 )  =   GetGeometry()[1].GetSolutionStepValue( DISPLACEMENT_Z );
        CurrentDisplacement( 9 )  =   GetGeometry()[1].GetSolutionStepValue( ROTATION_X );
        CurrentDisplacement( 10 )         =   GetGeometry()[1].GetSolutionStepValue( ROTATION_Y );
        CurrentDisplacement( 11 )         =   GetGeometry()[1].GetSolutionStepValue( ROTATION_Z );

        CalculateLHS( GlobalMatrix );
        noalias( rRightHandSideVector ) -= prod( GlobalMatrix, CurrentDisplacement );
        return;
    }




//************************************************************************************
//************************************************************************************

    void BeamElement::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
    {

        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag  = false;
        Vector temp = Vector();
        CalculateAll( rLeftHandSideMatrix, temp, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

    }



    void BeamElement::CalculateSectionProperties()

    {
        KRATOS_TRY

        array_1d<double, 3> x_0;
        array_1d<double, 3> x_1; // Vector que contiene coordenadas de los nodos.
        array_1d<double, 3> length;   // Vector que contiene la direccion de la barra.

        double minimo, maximo, B;
        const double b        = GetProperties()[BASE];
        const double h        = GetProperties()[HEIGHT];


        mInertia_x     = b * h * h * h / 12.0;
        mInertia_y     = b * b * b * h / 12.0;
        minimo         = std::min( b, h );
        maximo        = std::max( b, h );
        B        = ( 1.00 - 0.63 * ( minimo / maximo ) * ( 1 - ( pow( minimo, 4 ) / ( 12 * pow( maximo, 4 ) ) ) ) ) / 3; // constante torsional. Solo para secciones rectangulares.
        mInertia_Polar = B * minimo * minimo * minimo * maximo;
        mArea        = b * h;


        x_0( 0 ) = GetGeometry()[0].X0();
        x_0( 1 ) = GetGeometry()[0].Y0();
        x_0( 2 ) = GetGeometry()[0].Z0();
        x_1( 0 ) = GetGeometry()[1].X0();
        x_1( 1 ) = GetGeometry()[1].Y0();
        x_1( 2 ) = GetGeometry()[1].Z0();

        noalias( length ) = x_1 - x_0;
        mlength = std::sqrt( inner_prod( length, length ) );

        KRATOS_CATCH( "" )

    }


//************************************************************************************
//************************************************************************************

    void BeamElement::CalculateLocalMatrix( Matrix& LocalMatrix )
    {
        KRATOS_TRY




        if ( LocalMatrix.size1() != 12 || LocalMatrix.size2() != 12 )   // Matriz local de rigidez de la Estructura.
            LocalMatrix.resize( 12, 12, false );

        noalias( LocalMatrix )   = zero_matrix<double>( 12, 12 );

        //Inicializando matriz local de la estructura

        //const double mlength = GetGeometry().Length();
        const double Poisson = GetProperties()[POISSON_RATIO];

        const double Youngs  = GetProperties()[YOUNG_MODULUS];

        const double Elasticidad_Cortante = Youngs / ( 2.0 * ( 1.0 + Poisson ) );

        const double L =  mlength;

        const double LL   =  mlength * mlength;

        const double LLL =  mlength * mlength * mlength;

        double const EA   =  mArea          * Youngs;

        double const EIx  =  mInertia_x     * Youngs;

        double const EIy  =  mInertia_y     * Youngs;

        double const JG   =  mInertia_Polar * Elasticidad_Cortante;

        LocalMatrix( 0, 0 ) = ( EA ) / ( L );

        LocalMatrix( 6, 0 ) =   -( EA ) / ( L );

        LocalMatrix( 1, 1 ) = ( 12 * EIx ) / ( LLL );

        LocalMatrix( 5, 1 ) = ( 6 * EIx ) / ( LL );

        LocalMatrix( 7, 1 ) =   -( 12 * EIx ) / ( LLL );

        LocalMatrix( 11, 1 ) = ( 6 * EIx ) / ( LL );

        LocalMatrix( 2, 2 ) = ( 12 * EIy ) / ( LLL );

        LocalMatrix( 4, 2 ) =   -( 6 * EIy ) / ( LL );

        LocalMatrix( 8, 2 ) =   -( 12 * EIy ) / ( LLL );

        LocalMatrix( 10, 2 ) =   -( 6 * EIy ) / ( LL );

        LocalMatrix( 3, 3 ) = ( JG ) / L;

        LocalMatrix( 9, 3 ) =   -( JG ) / L;

        LocalMatrix( 2, 4 ) =   -( 6 * EIy ) / ( LL );

        LocalMatrix( 4, 4 ) = ( 4 * EIy ) / L;

        LocalMatrix( 8, 4 ) = ( 6 * EIy ) / ( LL );

        LocalMatrix( 10, 4 ) = ( 2 * EIy ) / L;

        LocalMatrix( 1, 5 ) = ( 6 * EIx ) / ( LL );

        LocalMatrix( 5, 5 ) = ( 4 * EIx ) / L;

        LocalMatrix( 7, 5 ) =   -( 6 * EIx ) / ( LL );

        LocalMatrix( 11, 5 ) = ( 2 * EIx ) / L;

        LocalMatrix( 0, 6 ) =    -( EA ) / ( L );

        LocalMatrix( 6, 6 ) = ( EA ) / ( L );

        LocalMatrix( 1, 7 ) =   -( 12 * EIx ) / ( LLL );

        LocalMatrix( 5, 7 ) =   -( 6 * EIx ) / ( LL );

        LocalMatrix( 7, 7 ) = ( 12 * EIx ) / ( LLL );

        LocalMatrix( 11, 7 ) = -( 6 * EIx ) / ( LL );

        LocalMatrix( 2, 8 ) =  -( 12 * EIy ) / ( LLL );

        LocalMatrix( 4, 8 ) = ( 6 * EIy ) / ( LL );

        LocalMatrix( 8, 8 ) = ( 12 * EIy ) / ( LLL );

        LocalMatrix( 10, 8 ) = ( 6 * EIy ) / ( LL );

        LocalMatrix( 3, 9 ) =   -( JG ) / L;

        LocalMatrix( 9, 9 ) = ( JG ) / L;

        LocalMatrix( 2, 10 ) =   -( 6 * EIy ) / ( LL );

        LocalMatrix( 4, 10 ) = ( 2 * EIy ) / L;

        LocalMatrix( 8, 10 ) = ( 6 * EIy ) / ( LL );

        LocalMatrix( 10, 10 ) = ( 4 * EIy ) / L;

        LocalMatrix( 1, 11 ) = ( 6 * EIx ) / ( LL );

        LocalMatrix( 5, 11 ) = ( 2 * EIx ) / L;

        LocalMatrix( 7, 11 ) =  -( 6 * EIx ) / ( LL );

        LocalMatrix( 11, 11 ) = ( 4 * EIx ) / L;



        KRATOS_CATCH( "" )

    }


//*****************************************************************************
//*****************************************************************************

    void BeamElement::CalculateTransformationMatrix( Matrix& Rotation )

    {

        KRATOS_TRY

        Vector Normal_zero( 9 ); // vector que contiene los cosenos directores.
        Vector x_zero( 6 );
        Vector Vector_zero( 3 );
        noalias( Normal_zero ) = zero_vector<double>( 9 );
        noalias( x_zero )      = zero_vector<double>( 6 );
        noalias( Vector_zero ) = zero_vector<double>( 3 );
        noalias( Rotation )    = zero_matrix<double> ( 12, 12 );

        double nx, ny, nz, teta, phi;
        //const double& mlength = GetGeometry().Length();

        x_zero( 0 ) = GetGeometry()[0].X0();
        x_zero( 1 ) = GetGeometry()[0].Y0();
        x_zero( 2 ) = GetGeometry()[0].Z0();
        x_zero( 3 ) = GetGeometry()[1].X0();
        x_zero( 4 ) = GetGeometry()[1].Y0();
        x_zero( 5 ) = GetGeometry()[1].Z0();

        for ( unsigned int i = 0; i < 3; i++ )
        {
            Vector_zero[i] = x_zero[i+3] - x_zero[i];
        }

        noalias( Normal_zero ) = Vector_zero * ( 1.00 / mlength );


        nx = Normal_zero[0];
        ny = Normal_zero[1];
        nz = Normal_zero[2];

        if ( nx == 0.0 )
        {
            teta = PI / 2;

            if ( ny == 0.0 )
            {
                teta = 0.0;
                phi  = PI / 2;
            }
            else
            {
                phi = atan( nz / sqrt( nx * nx + ny * ny ) );
            }
        }
        else
        {
            teta = atan( ny / nx );
            phi  = atan( nz / sqrt( nx * nx + ny * ny ) );
        }

        if ( nx < 0.0 )
            teta = teta + PI;

        Normal_zero[3] = -sin( teta );

        Normal_zero[4] =  cos( teta );

        Normal_zero[5] =  0.0;

        Normal_zero[6] = -nz * cos( teta );

        Normal_zero[7] = -nz * sin( teta );

        Normal_zero[8] =  nx * cos( teta ) + ny * sin( teta );


        // Creacion de la matriz de transformacion.
        for ( unsigned int kk = 0; kk < 12; kk += 3 )
        {
            for ( unsigned int i = 0; i < 3; i++ )
            {
                for ( unsigned int j = 0; j < 3; j++ )
                {
                    Rotation( i + kk, j + kk ) = Normal_zero( 3 * j + i );
                }
            }
        }

        KRATOS_CATCH( "" )

    }


//************************************************************************************
//************************************************************************************
    void BeamElement::CalculateBodyForce( Matrix& Rotation, Vector& LocalBody, Vector& GlobalBody )

    {
        KRATOS_TRY
        //Creacion de los vectores de cargas externas.
        // Fuerzas externas Uniformente distriduida. Por lo general es un dato suministrado por el usuario.
        // Cambiaro de posicion una vez terminado el programa.

        double alpha  =  0.00;
        double signo  =  1.00;
        //const double mlength = GetGeometry().Length();
        double  sino;
        double  cose;

        array_1d<double, 3> Weight;
        Weight[0]        =  GetProperties()[BODY_FORCE]( 0 );
        Weight[1]        =  GetProperties()[BODY_FORCE]( 1 );
        Weight[2]        =  GetProperties()[BODY_FORCE]( 2 );


        array_1d<double, 12 > Cargas_X = ZeroVector( 12 );
        array_1d<double, 12 > Cargas_Y = ZeroVector( 12 );
        array_1d<double, 12 > Cargas_Z = ZeroVector( 12 );

        array_1d<double, 2 > Load;
        array_1d<double, 6 > x_zero;

        Vector Normal_Loads;
        Vector Vector_zero;

        Normal_Loads.resize( 3, false );
        Vector_zero.resize( 3, false );


        x_zero( 0 ) = GetGeometry()[0].X0();
        x_zero( 1 ) = GetGeometry()[0].Y0();
        x_zero( 2 ) = GetGeometry()[0].Z0();
        x_zero( 3 ) = GetGeometry()[1].X0();
        x_zero( 4 ) = GetGeometry()[1].Y0();
        x_zero( 5 ) = GetGeometry()[1].Z0();



        for ( unsigned int i = 0; i < 3; i++ )
        {
            Vector_zero[i] = x_zero[i+3] - x_zero[i];
        }


        //Fuerza En X
        //***********************************
        if ( Weight[0] != 0.00 )
        {
            Normal_Loads[0]   = 0.00;
            Normal_Loads[1]   = Vector_zero[1] ;
            Normal_Loads[2]   = Vector_zero[2] ;

            if ( Vector_zero[0] < 0 )
            {
                signo = -1.00;
            }

            if ( norm_2( Normal_Loads ) == 0 || norm_2( Vector_zero ) == 0 )
            {
                alpha = signo * PI / 2;
            }
            else
            {
                alpha = inner_prod( Normal_Loads, Vector_zero ) / ( norm_2( Vector_zero ) * norm_2( Normal_Loads ) );
                alpha = signo * acos( alpha );
            }

            sino = sin( alpha );

            cose = cos( alpha );

            if ( fabs( sino ) < 1E-7 ) sino = 0.00;

            if ( fabs( cose ) < 1E-7 ) cose = 0.00;

            // las fuerzas consideradas son las de peso propio.
            Load[0] = mArea * Weight[0] * sino;    // Carga Axialmente Distribuida.

            Load[1] = mArea * Weight[0] * cose;    // Carga en la Direccion gravedad

            Cargas_X[0] =   Load[0] * mlength / 2.00;     // Fuerza en X;

            Cargas_X[1] =   -( Load[1] * mlength ) / 2.00;  // Fuerza en Y; graveded

            Cargas_X[2] =   0.00;                                                          // Fuerza en Z

            Cargas_X[3] =   0.00;                                                          // Momento Tersor X;

            Cargas_X[4] =   0.00;                                                          // Momento Y

            Cargas_X[5] =  -( Load[1] ) * mlength * mlength / 12.00;; // Momento Z

            Cargas_X[6] =   Load[0] * mlength / 2.00;

            Cargas_X[7] =   -( Load[1] ) * mlength / 2.00;

            Cargas_X[8] =    0.00;

            Cargas_X[9] =    0.00;

            Cargas_X[10] =   0.00;

            Cargas_X[11] = ( Load[1] ) * mlength * mlength / 12.00;

            noalias( GlobalBody ) = prod( Rotation, Cargas_X );  // Cargas externas en coordenadas globales.

            noalias( LocalBody )  = Cargas_X;

        }




        //Fuerza En Z
        //***********************************
        if ( Weight[2] != 0.00 )
        {
            Normal_Loads[0] = Vector_zero[0] ;
            Normal_Loads[1] = Vector_zero[1] ;
            Normal_Loads[2] = 0.00;

            if ( Vector_zero[2] < 0 )
            {
                signo = -1.00;
            }

            if ( norm_2( Normal_Loads ) == 0 || norm_2( Vector_zero ) == 0 )
            {
                alpha = signo * PI / 2;
            }
            else
            {
                alpha = inner_prod( Normal_Loads, Vector_zero ) / ( norm_2( Vector_zero ) * norm_2( Normal_Loads ) );
                alpha = signo * acos( alpha );
            }

            sino = sin( alpha );

            cose = cos( alpha );

            if ( fabs( sino ) < 1E-7 ) sino = 0.00;

            if ( fabs( cose ) < 1E-7 ) cose = 0.00;

            // las fuerzas consideradas son las de peso propio.
            Load[0] = mArea * Weight[2] * sino;    // Carga Axialmente Distribuida.

            Load[1] = mArea * Weight[2] * cose;    // Carga en la Direccion gravedad

            Cargas_Z[0] =   -Load[0] * mlength / 2.00;     // Fuerza en X;

            Cargas_Z[1] =   0.00;

            Cargas_Z[2] =   -( Load[1] * mlength ) / 2.00;  // Fuerza en Z; graveded                                                           // Fuerza en Z

            Cargas_Z[3] =   0.00;                                                          // Momento Tersor X;

            Cargas_Z[4] =   -( Load[1] ) * mlength * mlength / 12.00;                                                 // Momento Y

            Cargas_Z[5] =    0.00;

            Cargas_Z[6] =   -Load[0] * mlength / 2.00;

            Cargas_Z[7] =   0.00;

            Cargas_Z[8] =    -( Load[1] ) * mlength / 2.00;

            Cargas_Z[9] =   0.00;

            Cargas_Z[10] = ( Load[1] ) * mlength * mlength / 12.00;

            Cargas_Z[11] =  0.00;

            noalias( GlobalBody ) = prod( Rotation, Cargas_Z );  // Cargas externas en coordenadas globales.

            noalias( LocalBody )  = Cargas_Z;

        }

        //Fuerza En Y
        //***********************************
        if ( Weight[1] != 0.00 )
        {
            Normal_Loads      = ZeroVector( 3 );
            Normal_Loads[0] = Vector_zero[0] ;
            Normal_Loads[1] = 0.00 ;
            Normal_Loads[2] = Vector_zero[2];

            if ( Vector_zero[1] < 0 )
            {
                signo = -1.00;
            }

            if ( norm_2( Normal_Loads ) == 0 || norm_2( Vector_zero ) == 0 )
            {
                alpha = signo * PI / 2;
            }
            else
            {
                alpha = inner_prod( Normal_Loads, Vector_zero ) / ( norm_2( Vector_zero ) * norm_2( Normal_Loads ) );
                alpha = signo * acos( alpha );
            }

            sino = sin( alpha );

            cose = cos( alpha );

            if ( fabs( sino ) < 1E-7 ) sino = 0.00;

            if ( fabs( cose ) < 1E-7 ) cose = 0.00;

            // las fuerzas consideradas son las de peso propio.
            Load[0] = mArea * Weight[1] * sino;    // Carga Axialmente Distribuida.

            Load[1] = mArea * Weight[1] * cose;    // Carga en la Direccion gravedad


            Cargas_Y[0] =   -Load[0] * mlength / 2.00;     // Fuerza en X;

            Cargas_Y[1] =   -( Load[1] * mlength ) / 2.00;  // Fuerza en Y; graveded

            Cargas_Y[2] =   0.00;                                                          // Fuerza en Z

            Cargas_Y[3] =   0.00;                                                          // Momento Tersor X;

            Cargas_Y[4] =   0.00;                                                          // Momento Y

            Cargas_Y[5] =  -( Load[1] ) * mlength * mlength / 12.00;; // Momento Z

            Cargas_Y[6] =   -Load[0] * mlength / 2.00;

            Cargas_Y[7] =   -( Load[1] ) * mlength / 2.00;

            Cargas_Y[8] =    0.00;

            Cargas_Y[9] =    0.00;

            Cargas_Y[10] =   0.00;

            Cargas_Y[11] = ( Load[1] ) * mlength * mlength / 12.00;

            noalias( GlobalBody ) = prod( Rotation, Cargas_Y );  // Cargas externas en coordenadas globales.

            noalias( LocalBody )  = Cargas_Y;
        }

        KRATOS_CATCH( "" )

    }

//************************************************************************************
//************************************************************************************

    void BeamElement::MassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
    {

        KRATOS_TRY
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        unsigned int NumberOfNodes = GetGeometry().size();
        unsigned int MatSize = dimension * NumberOfNodes;

        if ( rMassMatrix.size1() != MatSize )
            rMassMatrix.resize( MatSize, MatSize, false );

        rMassMatrix = ZeroMatrix( MatSize, MatSize );

        //const double& mlength = GetGeometry().Length();

        double TotalMass = mArea * mlength * GetProperties()[DENSITY];

        Vector LumpFact;

        LumpFact = GetGeometry().LumpingFactors( LumpFact );

        for ( unsigned int i = 0; i < NumberOfNodes; i++ )
        {
            double temp = LumpFact[i] * TotalMass;

            for ( unsigned int j = 0; j < dimension; j++ )
            {
                unsigned int index = i * dimension + j;
                rMassMatrix( index, index ) = temp;

                if ( index == 3 || index == 4 || index == 5 )
                    rMassMatrix( index, index ) = 0.00;

            }
        }

        KRATOS_CATCH( "" )
    }



//************************************************************************************
    //************************************************************************************
    void BeamElement::GetFirstDerivativesVector( Vector& values, int Step )
    {
        KRATOS_TRY
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int MatSize = 2 * number_of_nodes * dim;

        if ( values.size() != MatSize )   values.resize( MatSize, false );

        for ( unsigned int i = 0;i < number_of_nodes;i++ )
        {
            unsigned int index =  i * number_of_nodes * dim;


            values[index] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
            values[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );
            values[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
            values[index + 3] = 0.00; // GetGeometry()[i].GetSolutionStepValue(ANGULAR_VELOCITY_X,Step);
            values[index + 4] = 0.00; // GetGeometry()[i].GetSolutionStepValue(ANGULAR_VELOCITY_Y,Step);
            values[index + 5] = 0.00; // GetGeometry()[i].GetSolutionStepValue(ANGULAR_VELOCITY_Z,Step);
        }


        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************
    void BeamElement::GetSecondDerivativesVector( Vector& values, int Step )
    {
        KRATOS_TRY
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int MatSize = 2 * number_of_nodes * dim;

        if ( values.size() != MatSize ) values.resize( MatSize, false );

        for ( unsigned int i = 0;i < number_of_nodes;i++ )
        {
            unsigned int index = i * number_of_nodes * dim;
            //KRATOS_WATCH(index)
            values[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
            values[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );
            values[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
            values[index + 3] = 0.00; // GetGeometry()[i].GetSolutionStepValue(ANGULAR_ACCELERATION_X,Step);
            values[index + 4] = 0.00; // GetGeometry()[i].GetSolutionStepValue(ANGULAR_ACCELERATION_Y,Step);
            values[index + 5] = 0.00; // GetGeometry()[i].GetSolutionStepValue(ANGULAR_ACCELERATION_Z,Step);

        }

        KRATOS_CATCH( "" )

    }

    //************************************************************************************
    //************************************************************************************

    void BeamElement::CalculateOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable,
            std::vector< array_1d<double, 3> >& Output,
            const ProcessInfo& rCurrentProcessInfo )
    {


        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( GeometryData::GI_GAUSS_3 );


        if ( Output.size() != integration_points.size() )
            Output.resize( integration_points.size() );

        Vector Stress;

        Vector Load1;

        Vector Load2;

        Vector Load3;

        CalculateLocalNodalStress( Stress );

        for ( unsigned int i = 0; i < Stress.size(); i++ )
        {
            if ( std::fabs( Stress[i] ) < 1E-6 ) Stress[i] = 0.00;
        }

        if ( rVariable == MOMENT )
        {

            /// Punto Inical
            Output[0][0] = Stress[3];
            Output[0][1] = Stress[4];
            Output[0][2] = Stress[5];

            if ( Id() == 3 ) KRATOS_WATCH( Stress )
                CalculateDistrubuitedBodyForce( 1, Load1 );

            CalculateDistrubuitedBodyForce( 2, Load2 );

            CalculateDistrubuitedBodyForce( 3, Load3 );

            Output[1][0] = 0.00;

            Output[1][1] = 0.00;

            Output[1][2] = CalculateInternalForces( Stress[5], Stress[1], Load2[1], mlength / 2 );


            Output[2][0] = 0.00;

            Output[2][1] = 0.00;

            Output[2][2] = CalculateInternalForces( Stress[5], Stress[1], Load2[1], mlength );

        }

        if ( rVariable == FORCE )
        {
            Output[0][0] = Stress[0];
            Output[0][1] = Stress[1];
            Output[0][2] = Stress[2];
            Output[1][0] = 0.00;
            Output[1][1] = 0.00;
            Output[1][2] = 0.00;
            Output[2][0] = Stress[6];
            Output[2][1] = Stress[7];
            Output[2][2] = Stress[8];
        }

    }


    double BeamElement::CalculateInternalForces( const double& Mo, const double& Po, const double& Load, const double& X )
    {

        return Mo -Po*X - 0.5 * Load * X * X;
    }

    void BeamElement::GetValueOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable,
            std::vector<array_1d<double, 3> >& rValues,
            const ProcessInfo& rCurrentProcessInfo )
    {

        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

    }

    IntegrationMethod  BeamElement::GetIntegrationMethod()
    {
        return GeometryData::GI_GAUSS_3;
    }

    void BeamElement::CalculateLocalNodalStress( Vector& Stress )
    {

        Matrix Rotation;
        Matrix LocalMatrix;
        array_1d<double, 12 > CurrentDisplacement;
        array_1d<double, 12 > LocalDisplacement;
        Vector LocalBody  = ZeroVector( 12 );
        Vector GlobalBody = ZeroVector( 12 );
        Rotation.resize( 12, 12, false );
        Stress.resize( 12, false );

        CurrentDisplacement( 0 )  =   GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_X );
        CurrentDisplacement( 1 )  =   GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_Y );
        CurrentDisplacement( 2 )  =   GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_Z );
        CurrentDisplacement( 3 )  =   GetGeometry()[0].GetSolutionStepValue( ROTATION_X );
        CurrentDisplacement( 4 )  =   GetGeometry()[0].GetSolutionStepValue( ROTATION_Y );
        CurrentDisplacement( 5 )  =   GetGeometry()[0].GetSolutionStepValue( ROTATION_Z );
        CurrentDisplacement( 6 )  =   GetGeometry()[1].GetSolutionStepValue( DISPLACEMENT_X );
        CurrentDisplacement( 7 )  =   GetGeometry()[1].GetSolutionStepValue( DISPLACEMENT_Y );
        CurrentDisplacement( 8 )  =   GetGeometry()[1].GetSolutionStepValue( DISPLACEMENT_Z );
        CurrentDisplacement( 9 )  =   GetGeometry()[1].GetSolutionStepValue( ROTATION_X );
        CurrentDisplacement( 10 )         =   GetGeometry()[1].GetSolutionStepValue( ROTATION_Y );
        CurrentDisplacement( 11 )         =   GetGeometry()[1].GetSolutionStepValue( ROTATION_Z );

        CalculateTransformationMatrix( Rotation );
        CalculateLocalMatrix( LocalMatrix );

        noalias( LocalDisplacement ) = prod( Matrix( trans( Rotation ) ), CurrentDisplacement );
        CalculateBodyForce( Rotation, LocalBody, GlobalBody );

        if ( Id() == 3 )
        {
            KRATOS_WATCH( CurrentDisplacement );
            KRATOS_WATCH( LocalDisplacement );
        }

        noalias( Stress ) = -LocalBody + prod( LocalMatrix, LocalDisplacement );

        return;

    }

    void BeamElement::CalculateDistrubuitedBodyForce( const int Direction, Vector& Load )
    {

        array_1d<double, 3> Weight;
        Load.resize( 2, false );
        Weight[0]        =  GetProperties()[BODY_FORCE]( 0 );
        Weight[1]        =  GetProperties()[BODY_FORCE]( 1 );
        Weight[2]        =  GetProperties()[BODY_FORCE]( 2 );

        double alpha  =  0.00;
        double signo  =  1.00;
        double  sino;
        double  cose;

        array_1d<double, 6 > x_zero;

        Vector Normal_Loads;
        Vector Vector_zero;
        Normal_Loads.resize( 3, false );
        Vector_zero.resize( 3, false );

        x_zero( 0 ) = GetGeometry()[0].X0();
        x_zero( 1 ) = GetGeometry()[0].Y0();
        x_zero( 2 ) = GetGeometry()[0].Z0();
        x_zero( 3 ) = GetGeometry()[1].X0();
        x_zero( 4 ) = GetGeometry()[1].Y0();
        x_zero( 5 ) = GetGeometry()[1].Z0();

        for ( unsigned int i = 0; i < 3; i++ )
        {
            Vector_zero[i] = x_zero[i+3] - x_zero[i];
        }

        if ( Direction == 1 )
            Normal_Loads[0]   = 0.00;

        Normal_Loads[1]   = Vector_zero[1] ;

        Normal_Loads[2]   = Vector_zero[2] ;

        if ( Vector_zero[0] < 0 )
        {
            signo = -1.00;
        }

        if ( norm_2( Normal_Loads ) == 0 || norm_2( Vector_zero ) == 0 )
        {
            alpha = signo * PI / 2;
        }
        else
        {
            alpha = inner_prod( Normal_Loads, Vector_zero ) / ( norm_2( Vector_zero ) * norm_2( Normal_Loads ) );
            alpha = signo * acos( alpha );
        }

        sino = sin( alpha );

        cose = cos( alpha );

        if ( fabs( sino ) < 1E-7 ) sino = 0.00;

        if ( fabs( cose ) < 1E-7 ) cose = 0.00;

        // las fuerzas consideradas son las de peso propio.
        Load[0] = mArea * Weight[0] * sino;    // Carga Axialmente Distribuida.

        Load[1] = mArea * Weight[0] * cose;    // Carga en la Direccion gravedad


        if ( Direction == 2 ) // 1=x, 2=y, 3=z
        {
            Normal_Loads    = ZeroVector( 3 );
            Normal_Loads[0] = Vector_zero[0] ;
            Normal_Loads[1] = 0.00 ;
            Normal_Loads[2] = Vector_zero[2];

            if ( Vector_zero[1] < 0 )
            {
                signo = -1.00;
            }

            if ( norm_2( Normal_Loads ) == 0 || norm_2( Vector_zero ) == 0 )
            {
                alpha = signo * PI / 2;
            }
            else
            {
                alpha = inner_prod( Normal_Loads, Vector_zero ) / ( norm_2( Vector_zero ) * norm_2( Normal_Loads ) );
                alpha = signo * acos( alpha );
            }

            sino = sin( alpha );

            cose = cos( alpha );

            if ( fabs( sino ) < 1E-7 ) sino = 0.00;

            if ( fabs( cose ) < 1E-7 ) cose = 0.00;


            // las fuerzas consideradas son las de peso propio.
            Load[0] = mArea * Weight[1] * sino;    // Carga Axialmente Distribuida.

            Load[1] = mArea * Weight[1] * cose;    // Carga en la Direccion gravedad
        }

        if ( Direction == 3 ) // 1=x, 2=y, 3=z
        {
            Normal_Loads[0] = Vector_zero[0] ;
            Normal_Loads[1] = Vector_zero[1] ;
            Normal_Loads[2] = 0.00;

            if ( Vector_zero[2] < 0 )
            {
                signo = -1.00;
            }

            if ( norm_2( Normal_Loads ) == 0 || norm_2( Vector_zero ) == 0 )
            {
                alpha = signo * PI / 2;
            }
            else
            {
                alpha = inner_prod( Normal_Loads, Vector_zero ) / ( norm_2( Vector_zero ) * norm_2( Normal_Loads ) );
                alpha = signo * acos( alpha );
            }

            sino = sin( alpha );

            cose = cos( alpha );

            if ( fabs( sino ) < 1E-7 ) sino = 0.00;

            if ( fabs( cose ) < 1E-7 ) cose = 0.00;

            // las fuerzas consideradas son las de peso propio.
            Load[0] = mArea * Weight[2] * sino;    // Carga Axialmente Distribuida.

            Load[1] = mArea * Weight[2] * cose;    // Carga en la Direccion gravedad
        }



    }




} // Namespace Kratos


