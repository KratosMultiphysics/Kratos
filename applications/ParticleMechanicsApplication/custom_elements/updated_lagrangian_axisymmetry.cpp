//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//


// System includes
#include <omp.h>
#include <sstream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/updated_lagrangian_axisymmetry.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "particle_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

	UpdatedLagrangianAxisymmetry::UpdatedLagrangianAxisymmetry( )
		: UpdatedLagrangian( )
	{
	  //DO NOT CALL IT: only needed for Register and Serialization!!!
	}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
    UpdatedLagrangianAxisymmetry::UpdatedLagrangianAxisymmetry( IndexType NewId, GeometryType::Pointer pGeometry )
            : UpdatedLagrangian( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

    UpdatedLagrangianAxisymmetry::UpdatedLagrangianAxisymmetry( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
            : UpdatedLagrangian( NewId, pGeometry, pProperties )
    {
    mFinalizedStep = true;
    }

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

    UpdatedLagrangianAxisymmetry::UpdatedLagrangianAxisymmetry( UpdatedLagrangianAxisymmetry const& rOther)
        :UpdatedLagrangian(rOther)
    {
    }

//*********************************OPERATIONS*****************************************
//************************************************************************************

    Element::Pointer UpdatedLagrangianAxisymmetry::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new UpdatedLagrangianAxisymmetry( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
    }

//************************************CLONE*******************************************
//************************************************************************************

    Element::Pointer UpdatedLagrangianAxisymmetry::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
    {

        UpdatedLagrangianAxisymmetry NewElement (NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

        return Element::Pointer( new UpdatedLagrangianAxisymmetry(NewElement) );
    }

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

    UpdatedLagrangianAxisymmetry::~UpdatedLagrangianAxisymmetry()
    {
    }

//************************************************************************************
    void UpdatedLagrangianAxisymmetry::Initialize()
    {
        KRATOS_TRY
        UpdatedLagrangian::Initialize();

        const array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
        const double mp_volume = this->GetValue(MP_VOLUME);
        const double mp_mass = mp_volume * 2* GetPI() * xg[0] * GetProperties()[DENSITY];
        this->SetValue(MP_MASS, mp_mass);

        mDeterminantF0 = 1;
        mDeformationGradientF0 = IdentityMatrix(3);

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianAxisymmetry::InitializeGeneralVariables (GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int strain_size = GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
        unsigned int voigtsize  = 4;

        rVariables.detF  = 1;

        rVariables.detF0 = 1;

        rVariables.detFT = 1;

        rVariables.detJ = 1;

        rVariables.B.resize( strain_size , number_of_nodes * dimension );

        rVariables.F.resize( 3, 3, false );
		rVariables.F = IdentityMatrix(3);

        rVariables.F0.resize( 3, 3, false );
        rVariables.F0 = IdentityMatrix(3);

        rVariables.FT.resize( 3, 3, false );
		rVariables.FT = IdentityMatrix(3);

        rVariables.ConstitutiveMatrix.resize( strain_size, strain_size, false );

        rVariables.StrainVector.resize( strain_size, false );

        rVariables.StressVector.resize( strain_size, false );

        rVariables.DN_DX.resize( number_of_nodes, dimension, false );
        rVariables.DN_De.resize( number_of_nodes, dimension, false );

        array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);

        rVariables.N = this->MPMShapeFunctionPointValues(rVariables.N, xg);

        // Reading shape functions local gradients
        rVariables.DN_De = this->MPMShapeFunctionsLocalGradients( rVariables.DN_De);

        // CurrentDisp is the unknown variable. It represents the nodal delta displacement. When it is predicted is equal to zero.
        rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);

        // Calculating the current jacobian from cartesian coordinates to parent coordinates for the MP element [dx_n+1/d£]
        rVariables.j = this->MPMJacobianDelta( rVariables.j, xg, rVariables.CurrentDisp);

        // Calculating the reference jacobian from cartesian coordinates to parent coordinates for the MP element [dx_n/d£]
        rVariables.J = this->MPMJacobian( rVariables.J, xg);

    }

//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianAxisymmetry::CalculateElementalSystem( LocalSystemComponents& rLocalSystem,
                                 ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY


        //create and initialize element variables:
        GeneralVariables Variables;

        this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);


        //create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);


        //set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();

        //std::cout<<"in CalculateElementalSystem 5"<<std::endl;
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN);

        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);


        //auxiliary terms
        Vector VolumeForce;


        //compute element kinematics B, F, DN_DX ...
        this->CalculateKinematics(Variables,rCurrentProcessInfo);

        //set general variables to constitutivelaw parameters
        this->SetGeneralVariables(Variables,Values);

        mConstitutiveLawVector->CalculateMaterialResponse(Values, Variables.StressMeasure);

        //the MP density is updated
        double MP_Density = (GetProperties()[DENSITY]) / Variables.detFT;
        this->SetValue(MP_DENSITY, MP_Density);
        //if(this->Id() == 1786 || this->Id() == 1836)
            //{
                //std::cout<<"density "<<this->Id() << GetProperties()[DENSITY]<<std::endl;
            //}

        //the integration weight is evaluated

        double MP_Volume = this->GetValue(MP_MASS)/this->GetValue(MP_DENSITY);

        //const double Radius = CalculateRadius(Variables.N, GetGeometry(), "Current");

        //const double AxiSymCoefficient = 2.0 * M_PI * Radius;


        //MP_Volume = AxiSymCoefficient * MP_Volume;

        this->SetValue(MP_VOLUME, MP_Volume);

        if ( rLocalSystem.CalculationFlags.Is(UpdatedLagrangian::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
        {

        //contributions to stiffness matrix calculated on the reference config
        this->CalculateAndAddLHS ( rLocalSystem, Variables, MP_Volume );

        }

        if ( rLocalSystem.CalculationFlags.Is(UpdatedLagrangian::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
        {
        //contribution to external forces
        VolumeForce  = this->CalculateVolumeForce( VolumeForce, Variables );

        this->CalculateAndAddRHS ( rLocalSystem, Variables, VolumeForce, MP_Volume );

        }



        KRATOS_CATCH( "" )
    }

    /**
     * Calculates the radius of axisymmetry
     * @param N: The Gauss Point shape function
     * @param Geom: The geometry studied
     * @return Radius: The radius of axisymmetry
     */

    double UpdatedLagrangianAxisymmetry::CalculateRadius(
        Vector& N,
        GeometryType& Geom,
        std::string ThisConfiguration = "Current"
        )
    {
        double Radius = 0.0;

        for (unsigned int iNode = 0; iNode < Geom.size(); iNode++)
        {
            // Displacement from the reference to the current configuration
            if (ThisConfiguration == "Current")
            {
                const array_1d<double, 3 > DeltaDisplacement = Geom[iNode].FastGetSolutionStepValue(DISPLACEMENT);
                const array_1d<double, 3 > ReferencePosition = Geom[iNode].Coordinates();

                const array_1d<double, 3 > CurrentPosition = ReferencePosition + DeltaDisplacement;
                Radius += CurrentPosition[0] * N[iNode];
                //if(this->Id() == 30211)
				//{

				//std::cout<<" CurrentPosition[0]"<< CurrentPosition[0]<<std::endl;
				//std::cout<<" N[iNode]"<< N[iNode]<<std::endl;
				//std::cout<<" Radius "<< Radius<<std::endl;


				//}
            }
            else
            {

                const array_1d<double, 3 > ReferencePosition = Geom[iNode].Coordinates();
                Radius += ReferencePosition[0] * N[iNode];
            }
        }

        return Radius;
    }


//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


    void UpdatedLagrangianAxisymmetry::CalculateKinematics(GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo)

    {
        KRATOS_TRY

        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        //Define the stress measure
        rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

        //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
        Matrix InvJ;

        MathUtils<double>::InvertMatrix( rVariables.J, InvJ, rVariables.detJ);



        rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);

        //Compute cartesian derivatives [dN/dx_n]
        rVariables.DN_DX = prod( rVariables.DN_De, InvJ);
        const double current_radius = CalculateRadius(rVariables.N, GetGeometry(), "Current");
        const double initial_radius = CalculateRadius(rVariables.N, GetGeometry(), "Initial");

        CalculateDeformationGradient (rVariables.DN_DX, rVariables.F, rVariables.CurrentDisp, current_radius, initial_radius);

		if(this->Id() == 10)
    {
		//std::cout<<"DF "<<DF<<std::endl;
		std::cout<<"rVariables.DN_DX "<<rVariables.DN_DX<<std::endl;
		std::cout<<"rVariables.CurrentDisp "<<rVariables.CurrentDisp<<std::endl;


	}
		//Calculating the inverse of the jacobian and the parameters needed [d£/(dx_n+1)]
        Matrix Invj;
        MathUtils<double>::InvertMatrix( rVariables.j, Invj, rVariables.detJ ); //overwrites detJ

        //Compute cartesian derivatives [dN/dx_n+1]
        rVariables.DN_DX = prod( rVariables.DN_De, Invj); //overwrites DX now is the current position dx


        //Matrix DF = prod( rVariables.j, InvJ);

        //// Axisymmetric case

        //DF.resize(3, 3); // We keep the old values
        //for (unsigned int index = 0; index < 1; index++)
            //{
                //DF(index, 2) = 0.0;
                //DF(2, index) = 0.0;
            //}


       //const double current_radius = CalculateRadius(rVariables.N, GetGeometry(), "Current");
       //const double initial_radius = CalculateRadius(rVariables.N, GetGeometry(), "Initial");
       //DF(2, 2) = current_radius/initial_radius;
       rVariables.CurrentRadius =  current_radius;
       rVariables.ReferenceRadius = initial_radius;


        //noalias(rVariables.F) = DF;



        //noalias( rVariables.F ) = prod( rVariables.j, InvJ);



        //Determinant of the Deformation Gradient F_n

        rVariables.detF0 = mDeterminantF0;
        rVariables.F0    = mDeformationGradientF0;

        if(this->Id() == 10)
    {
		//std::cout<<"DF "<<DF<<std::endl;
		std::cout<<"rVariables.F "<<rVariables.F<<std::endl;
		std::cout<<"rVariables.detF0 "<<rVariables.detF0<<std::endl;
		std::cout<<"rVariables.F0 "<<rVariables.F0<<std::endl;
		std::cout<<"rVariables.CurrentRadius "<<rVariables.CurrentRadius <<std::endl;
		std::cout<<"rVariables.ReferenceRadius "<<rVariables.ReferenceRadius<<std::endl;
		//std::cout<<"mDeformationGradientF0 "<<mDeformationGradientF0<<std::endl;

	}

        //if(this->Id() == 365)
        //{

        //std::cout<<"rVariables.DN_DX "<<this->Id()<<rVariables.DN_DX<<std::endl;
        //std::cout<<"rVariables.DN_De "<<this->Id()<<rVariables.DN_De<<std::endl;
        //std::cout<<"rVariables.J "<<this->Id()<<rVariables.J<<std::endl;
        //std::cout<<"rVariables.j "<<this->Id()<<rVariables.j<<std::endl;
        //std::cout<<"Invj "<<this->Id()<<Invj<<std::endl;
        //}

        //Compute the deformation matrix B
        this->CalculateDeformationMatrix(rVariables.B, rVariables.F, rVariables.DN_DX, rVariables.N);

        KRATOS_CATCH( "" )
    }
//************************************************************************************

    void UpdatedLagrangianAxisymmetry::CalculateDeformationMatrix(Matrix& rB,
            Matrix& rF,
            Matrix& rDN_DX,
            Vector& rN)
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
        double Radius = 0.0f;
        Radius = CalculateRadius(rN, GetGeometry(), "Current");

        rB.clear(); //set all components to zero

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            const unsigned int index = dimension * i;


            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 0 ) = rN[i]/Radius;
            rB( 3, index + 0 ) = rDN_DX( i, 1 );
            rB( 3, index + 1 ) = rDN_DX( i, 0 );
        }
        if(this->Id() == 10)
    {
		std::cout<<"Radius "<<Radius<<std::endl;
		std::cout<<"rB "<<rB<<std::endl;

		//std::cout<<"mDeformationGradientF0 "<<mDeformationGradientF0<<std::endl;

	}

        KRATOS_CATCH( "" )
    }
//************************************************************************************
//************************************************************************************
//*************************COMPUTE DEFORMATION GRADIENT*******************************
//************************************************************************************

	void UpdatedLagrangianAxisymmetry::CalculateDeformationGradient(const Matrix& rDN_DX,
			Matrix&  rF,
			Matrix&  rDeltaPosition,
			const double & rCurrentRadius,
			const double & rReferenceRadius)
	{
		KRATOS_TRY

		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
		const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

		rF = identity_matrix<double> ( 3 );

		if( dimension == 2 )
		{

			for ( unsigned int i = 0; i < number_of_nodes; i++ )
			{
				rF ( 0 , 0 ) += rDeltaPosition(i,0)*rDN_DX ( i , 0 );
				rF ( 0 , 1 ) += rDeltaPosition(i,0)*rDN_DX ( i , 1 );
				rF ( 1 , 0 ) += rDeltaPosition(i,1)*rDN_DX ( i , 0 );
				rF ( 1 , 1 ) += rDeltaPosition(i,1)*rDN_DX ( i , 1 );
			}

			rF ( 2 , 2 ) = rCurrentRadius/rReferenceRadius;
		}
		else if( dimension == 3)
		{

			std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;
		}
		else
		{

			KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

		}

		KRATOS_CATCH( "" )
	}

    //void UpdatedLagrangianAxisymmetry::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
    //{

        ////contribution of the internal and external forces
        //if( rLocalSystem.CalculationFlags.Is( UpdatedLagrangian::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
        //{

          //std::vector<VectorType>& rRightHandSideVectors = rLocalSystem.GetRightHandSideVectors();
          //const std::vector< Variable< VectorType > >& rRightHandSideVariables = rLocalSystem.GetRightHandSideVariables();
          //for( unsigned int i=0; i<rRightHandSideVariables.size(); i++ )
        //{
          //bool calculated = false;
          //if( rRightHandSideVariables[i] == EXTERNAL_FORCES_VECTOR ){
            //// operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
            //this->CalculateAndAddExternalForces( rRightHandSideVectors[i], rVariables, rVolumeForce, rIntegrationWeight );
            //calculated = true;
          //}

          //if( rRightHandSideVariables[i] == INTERNAL_FORCES_VECTOR ){
            //// operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
            //this->CalculateAndAddInternalForces( rRightHandSideVectors[i], rVariables, rIntegrationWeight );
            //calculated = true;
          //}

          //if(calculated == false)
            //{
              //KRATOS_THROW_ERROR( std::logic_error, " ELEMENT can not supply the required local system variable: ", rRightHandSideVariables[i] )
            //}

        //}
        //}
        //else{

          //VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();

          //// operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
          //this->CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

          //// operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
          //this->CalculateAndAddInternalForces( rRightHandSideVector, rVariables, rIntegrationWeight );
          ////KRATOS_WATCH( rRightHandSideVector )

          ////if(this->Id() == 15939)
    ////{
		////KRATOS_WATCH( rRightHandSideVector )
	////}
        //}


    //}
//************************************************************************************
//*********************Calculate the contribution of external force*******************

    //void UpdatedLagrangianAxisymmetry::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
            //GeneralVariables& rVariables,
            //Vector& rVolumeForce,
            //double& rIntegrationWeight)

    //{
        //KRATOS_TRY
        //unsigned int number_of_nodes = GetGeometry().PointsNumber();

        //unsigned int dimension = GetGeometry().WorkingSpaceDimension();



        //for ( unsigned int i = 0; i < number_of_nodes; i++ )
        //{
            //int index = dimension * i;

            //for ( unsigned int j = 0; j < dimension; j++ )
            //{
            //rRightHandSideVector[index + j] += rVariables.N[i] * rVolumeForce[j];

            //}

        //}


        //KRATOS_CATCH( "" )
    //}
//************************************************************************************
//************************************************************************************

    //void UpdatedLagrangianAxisymmetry::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
            //GeneralVariables & rVariables,
            //double& rIntegrationWeight)
    //{
        //KRATOS_TRY

        //VectorType InternalForces = rIntegrationWeight * prod( trans( rVariables.B ), rVariables.StressVector );
        //noalias( rRightHandSideVector ) -= InternalForces;




        //KRATOS_CATCH( "" )
    //}
//************************************************************************************
//************************************************************************************

    //void UpdatedLagrangianAxisymmetry::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
    //{

      ////contributions of the stiffness matrix calculated on the reference configuration
    //if( rLocalSystem.CalculationFlags.Is( UpdatedLagrangian::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
        //{
          //std::vector<MatrixType>& rLeftHandSideMatrices = rLocalSystem.GetLeftHandSideMatrices();
          //const std::vector< Variable< MatrixType > >& rLeftHandSideVariables = rLocalSystem.GetLeftHandSideVariables();

          //for( unsigned int i=0; i<rLeftHandSideVariables.size(); i++ )
        //{
          //bool calculated = false;
          //if( rLeftHandSideVariables[i] == MATERIAL_STIFFNESS_MATRIX ){
            //// operation performed: add Km to the rLefsHandSideMatrix
            //this->CalculateAndAddKuum( rLeftHandSideMatrices[i], rVariables, rIntegrationWeight );
            //calculated = true;
          //}

          //if( rLeftHandSideVariables[i] == GEOMETRIC_STIFFNESS_MATRIX ){
            //// operation performed: add Kg to the rLefsHandSideMatrix
            //this->CalculateAndAddKuug( rLeftHandSideMatrices[i], rVariables, rIntegrationWeight );
            //calculated = true;
          //}

          //if(calculated == false)
            //{
              //KRATOS_THROW_ERROR(std::logic_error, " ELEMENT can not supply the required local system variable: ",rLeftHandSideVariables[i])
            //}

        //}
        //}
    //else{

        //MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

        //this->CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

        //this->CalculateAndAddKuug( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

        ////if(this->Id() == 15939)
    ////{
		////KRATOS_WATCH( rLeftHandSideMatrix )
	////}
      //}


    //}
//************************************************************************************
//************************************************************************************

    //void UpdatedLagrangianAxisymmetry::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
            //GeneralVariables& rVariables,
            //double& rIntegrationWeight
                                                      //)
    //{
        //KRATOS_TRY
        ////std::stringstream ss;



        //noalias( rLeftHandSideMatrix ) += prod( trans( rVariables.B ),  rIntegrationWeight * Matrix( prod( rVariables.ConstitutiveMatrix, rVariables.B ) ) );


        ////std::cout << ss.str();

        //KRATOS_CATCH( "" )
    //}



//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianAxisymmetry::CalculateAndAddKuug(MatrixType& rK,
            GeneralVariables& rVariables,
            double& rIntegrationWeight)

    {
        KRATOS_TRY
		MatrixType Kh=rK;
        //unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        ////std::cout<<" StressVector "<<rVariables.StressVector<<std::endl;
        //Matrix StressTensor = MathUtils<double>::StressVectorToTensor( rVariables.StressVector );
        //Matrix ReducedKg = prod( rVariables.DN_DX, rIntegrationWeight * Matrix( prod( StressTensor, trans( rVariables.DN_DX ) ) ) ); //to be optimized
        //MathUtils<double>::ExpandAndAddReducedMatrix( rLeftHandSideMatrix, ReducedKg, dimension );

		//if(this->Id() == 9794)
        //{
			//std::cout<<" Kgeo "<<rLeftHandSideMatrix-Kh<<std::endl;
		//}

		const unsigned int number_of_nodes = GetGeometry().size();

		//Matrix Kh = rK;

		// axisymmetric geometric matrix

		double alpha1 = 0;
		double alpha2 = 0;
		double alpha3 = 0;

		unsigned int indexi = 0;
		unsigned int indexj = 0;

		const double Radius = CalculateRadius(rVariables.N, GetGeometry(), "Current");

		for ( unsigned int i = 0; i < number_of_nodes; i++ )
		{
			indexj =0;
			for ( unsigned int j = 0; j < number_of_nodes; j++ )
			{
				alpha1 = rVariables.DN_DX(j,0) * ( rVariables.DN_DX(i,0) * rVariables.StressVector[0] + rVariables.DN_DX(i,1) * rVariables.StressVector[3] );
				alpha2 = rVariables.DN_DX(j,1) * ( rVariables.DN_DX(i,0) * rVariables.StressVector[3] + rVariables.DN_DX(i,1) * rVariables.StressVector[1] );
				alpha3 = rVariables.N[i] * rVariables.N[j] * rVariables.StressVector[2] * (1.0/Radius*Radius);

				rK(indexi,indexj)     += (alpha1 + alpha2 + alpha3) * rIntegrationWeight ;
				rK(indexi+1,indexj+1) += (alpha1 + alpha2) * rIntegrationWeight ;

				indexj+=2;
			}

			indexi+=2;

		}
		if(this->Id() == 10)
    {

		std::cout<<"Kgeo "<<rK-Kh<<std::endl;


	}

        KRATOS_CATCH( "" )
    }
//************************************CALCULATE VOLUME CHANGE*************************
//************************************************************************************

    double& UpdatedLagrangianAxisymmetry::CalculateVolumeChange( double& rVolumeChange, GeneralVariables& rVariables )
    {
        KRATOS_TRY

        rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0);

        return rVolumeChange;

        KRATOS_CATCH( "" )
    }
//************************************CALCULATE VOLUME ACCELERATION*******************
//************************************************************************************

    Vector& UpdatedLagrangianAxisymmetry::CalculateVolumeForce( Vector& rVolumeForce, GeneralVariables& rVariables )
    {
        KRATOS_TRY

        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        rVolumeForce = ZeroVector(dimension);


        rVolumeForce = this->GetValue(MP_VOLUME_ACCELERATION)* this->GetValue(MP_MASS);


        return rVolumeForce;

        KRATOS_CATCH( "" )
    }
//************************************************************************************
//************************************************************************************
    void UpdatedLagrangianAxisymmetry::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
    {
        //create local system components
        LocalSystemComponents LocalSystem;

        //calculation flags
        LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_RHS_VECTOR);

        MatrixType LeftHandSideMatrix = Matrix();

        //Initialize sizes for the system components:
        this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

        //Set Variables to Local system components
        LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
        LocalSystem.SetRightHandSideVector(rRightHandSideVector);
		//std::cout<<" HERE 2"<<std::endl;
        //Calculate elemental system
        CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );
    }

//************************************************************************************
//************************************************************************************


    void UpdatedLagrangianAxisymmetry::CalculateRightHandSide( std::vector< VectorType >& rRightHandSideVectors, const std::vector< Variable< VectorType > >& rRHSVariables, ProcessInfo& rCurrentProcessInfo )
    {
        //create local system components
        LocalSystemComponents LocalSystem;

        //calculation flags
        LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_RHS_VECTOR);
        LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);

        MatrixType LeftHandSideMatrix = Matrix();

        //Initialize sizes for the system components:
        if( rRHSVariables.size() != rRightHandSideVectors.size() )
          rRightHandSideVectors.resize(rRHSVariables.size());

        for( unsigned int i=0; i<rRightHandSideVectors.size(); i++ )
          {
        this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVectors[i], LocalSystem.CalculationFlags );
          }

        //Set Variables to Local system components
        LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
        LocalSystem.SetRightHandSideVectors(rRightHandSideVectors);

        LocalSystem.SetRightHandSideVariables(rRHSVariables);
		//std::cout<<" HERE 3"<<std::endl;
        //Calculate elemental system
        CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );
    }


//************************************************************************************
//************************************************************************************


    void UpdatedLagrangianAxisymmetry::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
    {
        //create local system components
        LocalSystemComponents LocalSystem;

        //calculation flags
        LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_LHS_MATRIX);

        VectorType RightHandSideVector = Vector();

        //Initialize sizes for the system components:
        this->InitializeSystemMatrices( rLeftHandSideMatrix, RightHandSideVector, LocalSystem.CalculationFlags );

        //Set Variables to Local system components
        LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
        LocalSystem.SetRightHandSideVector(RightHandSideVector);

        //Calculate elemental system
        CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

    }
//************************************************************************************
//************************************************************************************


    void UpdatedLagrangianAxisymmetry::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
    {

        //create local system components
        LocalSystemComponents LocalSystem;

        //calculation flags
        LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_LHS_MATRIX);
        LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_RHS_VECTOR);

        //Initialize sizes for the system components:
        this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

        //Set Variables to Local system components
        LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
        LocalSystem.SetRightHandSideVector(rRightHandSideVector);

        //Calculate elemental system
        //std::cout<<" in CalculateLocalSystem ends"<<std::endl;
        CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );



    }


//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianAxisymmetry::CalculateLocalSystem( std::vector< MatrixType >& rLeftHandSideMatrices,
                                 const std::vector< Variable< MatrixType > >& rLHSVariables,
                                 std::vector< VectorType >& rRightHandSideVectors,
                                 const std::vector< Variable< VectorType > >& rRHSVariables,
                                 ProcessInfo& rCurrentProcessInfo )
    {
        //create local system components
        LocalSystemComponents LocalSystem;

        //calculation flags
        LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_LHS_MATRIX_WITH_COMPONENTS);
        LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);


        //Initialize sizes for the system components:
        if( rLHSVariables.size() != rLeftHandSideMatrices.size() )
        rLeftHandSideMatrices.resize(rLHSVariables.size());

        if( rRHSVariables.size() != rRightHandSideVectors.size() )
          rRightHandSideVectors.resize(rRHSVariables.size());

        LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_LHS_MATRIX);
        for( unsigned int i=0; i<rLeftHandSideMatrices.size(); i++ )
          {
        //Note: rRightHandSideVectors.size() > 0
        this->InitializeSystemMatrices( rLeftHandSideMatrices[i], rRightHandSideVectors[0], LocalSystem.CalculationFlags );
          }

        LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_RHS_VECTOR);
        LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_LHS_MATRIX,false);

        for( unsigned int i=0; i<rRightHandSideVectors.size(); i++ )
          {
        //Note: rLeftHandSideMatrices.size() > 0
            this->InitializeSystemMatrices( rLeftHandSideMatrices[0], rRightHandSideVectors[i], LocalSystem.CalculationFlags );
          }
        LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_LHS_MATRIX,true);


        //Set Variables to Local system components
        LocalSystem.SetLeftHandSideMatrices(rLeftHandSideMatrices);
        LocalSystem.SetRightHandSideVectors(rRightHandSideVectors);

        LocalSystem.SetLeftHandSideVariables(rLHSVariables);
        LocalSystem.SetRightHandSideVariables(rRHSVariables);
		//std::cout<<" HERE 1"<<std::endl;
        //Calculate elemental system
        CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

    }

////***********************************************************************************
////***********************************************************************************
	void UpdatedLagrangianAxisymmetry::Calculate(const Variable<double>& rVariable,
                           double& Output,
                           const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

		if (rVariable == DENSITY)
		{
		unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
        GeneralVariables Variables;
		Matrix J0 = ZeroMatrix(dimension, dimension);

        J0 = this->MPMJacobian(J0, xg);

        //calculating and storing inverse and the determinant of the jacobian
        MathUtils<double>::InvertMatrix( J0, mInverseJ0, mDeterminantJ0 );

		Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);
		double MP_Mass = this->GetValue(MP_MASS);

		for (unsigned int i=0;i<number_of_nodes;i++)

            {

				GetGeometry()[i].GetSolutionStepValue(AUX_R) += Variables.N[i] * (MP_Mass);// - AUX_MP_Mass);

			}
	    }


		else if (rVariable == MP_EQUIVALENT_PLASTIC_STRAIN)
		{
		const unsigned int number_of_nodes = GetGeometry().PointsNumber();

		double SmoothMPPlasticStrain = 0.0;

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
            {
				double NodalEqPlasticStrain = GetGeometry()[i].GetSolutionStepValue(NODAL_EQ_PLASTICSTRAIN, 0);


				SmoothMPPlasticStrain += NodalEqPlasticStrain * Variables.N[i];

			}
        this->SetValue(MP_EQUIVALENT_PLASTIC_STRAIN, SmoothMPPlasticStrain);


	    }


	    KRATOS_CATCH( "" )
	}




	void UpdatedLagrangianAxisymmetry::Calculate(const Variable<array_1d<double, 3 > >& rVariable,
                           array_1d<double, 3 > & Output,
                           const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY


	    if(rVariable == VELOCITY)
	    {
		unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
        GeneralVariables Variables;
		Matrix J0 = ZeroMatrix(dimension, dimension);

        J0 = this->MPMJacobian(J0, xg);

        //calculating and storing inverse and the determinant of the jacobian
        MathUtils<double>::InvertMatrix( J0, mInverseJ0, mDeterminantJ0 );

		Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);
		array_1d<double,3>& MP_Velocity = this->GetValue(MP_VELOCITY);
		double MP_Mass = this->GetValue(MP_MASS);

		array_1d<double,3> NodalAuxRVel;
		for (unsigned int i=0;i<number_of_nodes;i++)

            {
				for (unsigned int j = 0; j < dimension; j++)
                {
				NodalAuxRVel[j] = Variables.N[i] * MP_Mass * MP_Velocity[j];
				}
				GetGeometry()[i].GetSolutionStepValue(AUX_R_VEL) += NodalAuxRVel;

				//std::cout<<" Variables.N[i] "<< Variables.N[i]<<std::endl;
				//std::cout<<" MP_Mass "<< MP_Mass<<std::endl;
				//std::cout<<" MP_Velocity "<< MP_Velocity<<std::endl;
				//std::cout<<" NODE ID "<<  GetGeometry()[i].Id()<<std::endl;
				//std::cout<<" NodalAuxRVel in the element "<< GetGeometry()[i].GetSolutionStepValue(AUX_R_VEL)<<std::endl;
			}
		}
		if(rVariable == ACCELERATION)
	    {
		unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
        GeneralVariables Variables;
		Matrix J0 = ZeroMatrix(dimension, dimension);

        J0 = this->MPMJacobian(J0, xg);

        //calculating and storing inverse and the determinant of the jacobian
        MathUtils<double>::InvertMatrix( J0, mInverseJ0, mDeterminantJ0 );

		Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);
		array_1d<double,3>& MP_Acceleration = this->GetValue(MP_ACCELERATION);
		double MP_Mass = this->GetValue(MP_MASS);

		array_1d<double,3> NodalAuxRAcc;
		for (unsigned int i=0;i<number_of_nodes;i++)

            {
				for (unsigned int j = 0; j < dimension; j++)
                {
				NodalAuxRAcc[j] = Variables.N[i] * MP_Mass * MP_Acceleration[j];
				}
				GetGeometry()[i].GetSolutionStepValue(AUX_R_ACC) += NodalAuxRAcc;
			}
		}

		KRATOS_CATCH( "" )
	}
//*******************************************************************************************


////************************************************************************************
	//void UpdatedLagrangian::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
    //{
		//unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		//array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
		//GeneralVariables Variables;
		//Matrix J0 = ZeroMatrix(dimension, dimension);

        //J0 = this->MPMJacobian(J0, xg);

        ////calculating and storing inverse and the determinant of the jacobian
        //MathUtils<double>::InvertMatrix( J0, mInverseJ0, mDeterminantJ0 );

		//Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);
		//mConstitutiveLawVector->InitializeSolutionStep( GetProperties(),
                    //GetGeometry(), Variables.N, rCurrentProcessInfo );

            //mFinalizedStep = false;
     //}




    void UpdatedLagrangianAxisymmetry::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
    {
        // In the Initialize of each time step the nodal initial conditions are evaluated
        //1. first of all I need to evaluate the MP momentum and MP_inertia



        int MP_bool = this->GetValue(MP_BOOL);

        //std::cout<<" in InitializeSolutionStep2"<<std::endl;
            unsigned int dimension = GetGeometry().WorkingSpaceDimension();
            const unsigned int number_of_nodes = GetGeometry().PointsNumber();
            array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
            GeneralVariables Variables;
            //this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);



            Matrix J0 = ZeroMatrix(dimension, dimension);

            J0 = this->MPMJacobian(J0, xg);

            //calculating and storing inverse and the determinant of the jacobian
            MathUtils<double>::InvertMatrix( J0, mInverseJ0, mDeterminantJ0 );

            Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);

            mConstitutiveLawVector->InitializeSolutionStep( GetProperties(),
                    GetGeometry(), Variables.N, rCurrentProcessInfo );

            mFinalizedStep = false;



            array_1d<double,3>& MP_Velocity = this->GetValue(MP_VELOCITY);
            array_1d<double,3>& MP_Acceleration = this->GetValue(MP_ACCELERATION);
            array_1d<double,3>& AUX_MP_Velocity = this->GetValue(AUX_MP_VELOCITY);
            array_1d<double,3>& AUX_MP_Acceleration = this->GetValue(AUX_MP_ACCELERATION);
            double MP_Mass = this->GetValue(MP_MASS);
            array_1d<double,3> MP_Momentum;
            array_1d<double,3> MP_Inertia;
            array_1d<double,3> NodalMomentum;
            array_1d<double,3> NodalInertia;

           for (unsigned int j=0;j<number_of_nodes;j++)
            {
                //these are the values of nodal velocity and nodal acceleration evaluated in the initialize solution step

              for (unsigned int l=0;l<number_of_nodes;l++)

				{
					array_1d<double, 3 > & NodalAcceleration = GetGeometry()[l].FastGetSolutionStepValue(ACCELERATION,1);
                    array_1d<double, 3 > & NodalVelocity = GetGeometry()[l].FastGetSolutionStepValue(VELOCITY,1);
                //std::cout<<"NodalVelocity "<< GetGeometry()[j].Id()<<std::endl;
					for (unsigned int k = 0; k < dimension; k++)
					{
					AUX_MP_Velocity[k] += Variables.N[j] *Variables.N[l] * NodalVelocity[k];
					AUX_MP_Acceleration[k] += Variables.N[j] *Variables.N[l] * NodalAcceleration[k];
					}
				}
            }



            // Here MP contribution in terms of momentum, inertia and mass are added

            for ( unsigned int i = 0; i < number_of_nodes; i++ )
            {
                for (unsigned int j = 0; j < dimension; j++)
                {
                    NodalMomentum[j] = Variables.N[i] * (MP_Velocity[j] - AUX_MP_Velocity[j]) * MP_Mass;
                    NodalInertia[j] = Variables.N[i] * (MP_Acceleration[j] - AUX_MP_Acceleration[j]) * MP_Mass;

                }
                GetGeometry()[i].GetSolutionStepValue(NODAL_MOMENTUM, 0) += NodalMomentum;
                GetGeometry()[i].GetSolutionStepValue(NODAL_INERTIA, 0) += NodalInertia;

                GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0) += Variables.N[i] * MP_Mass;

            }

            AUX_MP_Velocity.clear();
            AUX_MP_Acceleration.clear();






    }

////************************************************************************************
////************************************************************************************


////************************************************************************************
////************************************************************************************
    void UpdatedLagrangianAxisymmetry::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
    {

    }

////************************************************************************************
////************************************************************************************

    void UpdatedLagrangianAxisymmetry::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
    {

    }

////************************************************************************************
////************************************************************************************

    void UpdatedLagrangianAxisymmetry::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY


        //std::cout<<"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"<<std::endl;
		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
	   const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
	   unsigned int voigtsize  = 4;


        //Vector NodalStress = ZeroVector(voigtsize);
        array_1d<double,3> NodalStress;
        double EqPlasticStrain;
        //create and initialize element variables:
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

        //create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        //set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
        //ConstitutiveLawOptions.Set(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY);
        //compute element kinematics B, F, DN_DX ...
        this->CalculateKinematics(Variables, rCurrentProcessInfo);

        //set general variables to constitutivelaw parameters
        this->SetGeneralVariables(Variables,Values);

        //call the constitutive law to update material variables
        mConstitutiveLawVector->FinalizeMaterialResponse(Values, Variables.StressMeasure);

        //call the constitutive law to finalize the solution step
        mConstitutiveLawVector->FinalizeSolutionStep( GetProperties(),
                GetGeometry(),
                Variables.N,
                rCurrentProcessInfo );

        //if(this->Id() == 550 || this->Id() == 551 || this->Id() == 552 || this->Id() == 553)
        //{

        //std::cout<<" StressVector after constitutive law in the finalize"<< Variables.StressVector<<std::endl;
        //std::cout<<" total def gradient in the finalize"<< Variables.detFT<<std::endl;
		//}
        //call the element internal variables update
        this->FinalizeStepVariables(Variables, rCurrentProcessInfo);


        double EquivalentPlasticStrain = this->GetValue(MP_EQUIVALENT_PLASTIC_STRAIN );
        Vector StressVector = this->GetValue(MP_CAUCHY_STRESS_VECTOR);

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
            {
                //GetGeometry()[i].GetSolutionStepValue(STRESSES, 0) = ZeroVector(voigtsize);
                if(GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0) > 1e-4)
                {
				for ( unsigned int j = 0; j< voigtsize; j++)
					{
					NodalStress[j] = Variables.N[i] * StressVector[j] * this->GetValue(MP_MASS) / GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0);

					}
			    EqPlasticStrain = Variables.N[i] * EquivalentPlasticStrain * this->GetValue(MP_MASS) / GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0);
			    GetGeometry()[i].GetSolutionStepValue(NODAL_EQ_PLASTICSTRAIN, 0) += EqPlasticStrain;
                GetGeometry()[i].GetSolutionStepValue(NODAL_STRESSES, 0) += NodalStress;
                //std::cout<<" GetGeometry()[i].GetSolutionStepValue(NODAL_STRESSES, 0) "<<  GetGeometry()[i].GetSolutionStepValue(NODAL_STRESSES, 0)<<std::endl;
			    }
			}


        mFinalizedStep = true;

        KRATOS_CATCH( "" )
    }


////************************************************************************************************************
	void UpdatedLagrangianAxisymmetry::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY


		GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);
		if (rVariable == MP_CAUCHY_STRESS_VECTOR)
		{
		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		unsigned int voigtsize  = 3;
		array_1d<double,3> SmoothMPStress;

        if( dimension == 3 )
        {
            voigtsize  = 6;
            array_1d<double,6> SmoothMPStress;
        }
        //Vector SmoothMPStress = ZeroVector(voigtsize);

        SmoothMPStress.clear();
        //Vector NodalStress0 = GetGeometry()[0].GetSolutionStepValue(STRESSES, 0);
        //Vector NodalStress1 = GetGeometry()[1].GetSolutionStepValue(STRESSES, 0);
        //Vector NodalStress2 = GetGeometry()[2].GetSolutionStepValue(STRESSES, 0);
        //Vector NodalStress3 = GetGeometry()[3].GetSolutionStepValue(STRESSES, 0);

        //SmoothMPStress = NodalStress0 * Variables.N[0] + NodalStress1 * Variables.N[1] + NodalStress2 * Variables.N[2] + NodalStress3 * Variables.N[3];
        //if (this->Id() == 541 || this->Id() == 534 || this->Id() == 538)
        //{
			//std::cout<<" GetGeometry()[0].GetSolutionStepValue(STRESSES, 0); "<< " Id "<< GetGeometry()[0].GetSolutionStepValue(STRESSES, 0)<<std::endl;
		//}

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
            {
				array_1d<double,3> NodalStress = GetGeometry()[i].GetSolutionStepValue(NODAL_STRESSES, 0);
				if( dimension == 3 )
				{
					array_1d<double,6> NodalStress = GetGeometry()[i].GetSolutionStepValue(NODAL_STRESSES, 0);
				}
                for ( unsigned int j = 0; j< voigtsize; j++)
				{
				SmoothMPStress[j] += NodalStress[j] * Variables.N[i];
				}
			}
        this->SetValue(MP_CAUCHY_STRESS_VECTOR, SmoothMPStress);
        //std::cout<<"cccccccccccccccccccccccccccccc"<<std::endl;

        //if (this->Id() == 541 || this->Id() == 534 || this->Id() == 538)
        //{
			//std::cout<<" MP_CAUCHY_STRESS_VECTOR "<< " Id "<< this->GetValue(MP_CAUCHY_STRESS_VECTOR)<<std::endl;
		//}

	    }

		KRATOS_CATCH( "" )
	}
////************************************************************************************
////************************************************************************************

    void UpdatedLagrangianAxisymmetry::FinalizeStepVariables( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
    {
    //update internal (historical) variables
    mDeterminantF0         = rVariables.detF* rVariables.detF0;
    mDeformationGradientF0 = prod(rVariables.F, rVariables.F0);

    if(this->Id() == 10)
    {
		std::cout<<"mDeterminantF0 "<<mDeterminantF0<<std::endl;
		std::cout<<"mDeformationGradientF0 "<<mDeformationGradientF0<<std::endl;

	}

    this->SetValue(MP_CAUCHY_STRESS_VECTOR, rVariables.StressVector);




    //std::cout<<"StressVector in element "<<rVariables.StressVector<<std::endl;
    this->SetValue(MP_ALMANSI_STRAIN_VECTOR, rVariables.StrainVector);

    this->SetValue(MP_JACOBIAN, mDeterminantF0);

    //double MP_Pressure = mConstitutiveLawVector->GetValue(MP_PRESSURE, MP_Pressure );
    //this->SetValue(MP_CONSTITUTIVE_PRESSURE, MP_Pressure);
    //this->SetValue(MP_PRESSURE, MP_Pressure);

    double EquivalentPlasticStrain = mConstitutiveLawVector->GetValue(PLASTIC_STRAIN, EquivalentPlasticStrain );
    ////std::cout<<" EquivalentPlasticStrain in the element "<<EquivalentPlasticStrain<<std::endl;
    //if(EquivalentPlasticStrain < 1.0 && EquivalentPlasticStrain >= 0.0)
	this->SetValue(MP_EQUIVALENT_PLASTIC_STRAIN, EquivalentPlasticStrain);

    double Damage = mConstitutiveLawVector->GetValue(DAMAGE_VARIABLE, Damage );
    this->SetValue(MP_DAMAGE, Damage);

    double Pressure = mConstitutiveLawVector->GetValue(PRESSURE, Pressure );
    this->SetValue(MP_PRESSURE, Pressure);

    MathUtils<double>::InvertMatrix( rVariables.j, mInverseJ, rVariables.detJ );


    this->UpdateGaussPoint(rVariables, rCurrentProcessInfo);

    }

//************************************************************************************
//************************************************************************************
    /**
     * The position of the Gauss points/Material points is updated
     */

    void UpdatedLagrangianAxisymmetry::UpdateGaussPoint( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);

        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        array_1d<double,3> xg = this->GetValue(GAUSS_COORD);
        array_1d<double,3> MP_PreviousAcceleration = this->GetValue(MP_ACCELERATION);
        array_1d<double,3> MP_PreviousVelocity = this->GetValue(MP_VELOCITY);
        double MP_Mass = this->GetValue(MP_MASS);
        array_1d<double,3> delta_xg = ZeroVector(3);
        array_1d<double,3> MP_Acceleration = ZeroVector(3);
        array_1d<double,3> MP_Velocity = ZeroVector(3);
        array_1d<double,3> MP_AuxVelocity = ZeroVector(3);
        double DeltaTime = rCurrentProcessInfo[DELTA_TIME];

        array_1d<double,3> MP_DeltaVelocity = ZeroVector(3);
        rVariables.N = this->MPMShapeFunctionPointValues(rVariables.N, xg);
        int MP_number = this->GetValue(MP_NUMBER);

        //double total_nodal_mass = 0.0;
        //for ( unsigned int i = 0; i < number_of_nodes; i++ )
        //{
            //total_nodal_mass += GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0);
        //}
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            if (rVariables.N[i] > 1e-16)
            {
            array_1d<double, 3 > & NodalAcceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION);
            array_1d<double, 3 > & NodalVelocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3 > & PreviousNodalVelocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,1);
            double NodalMass = GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0);
            array_1d<double,3> NodalMomentum = NodalMass * NodalVelocity;
            array_1d<double,3> NodalInertia = NodalMass * NodalAcceleration;

            //if (this->Id() == 8752)// || this->Id() == 1513)
            //{
                //std::cout<< "Nodal ID "<< GetGeometry()[i].Id()<<std::endl;
                //std::cout<< "NodalAcceleration "<<NodalAcceleration<<std::endl;
                //std::cout<< "NodalVelocity "<<NodalVelocity<<std::endl;
                //std::cout<< "NodalMass "<<NodalMass<<std::endl;

                //std::cout<< "rVariables.N "<<rVariables.N<<std::endl;
            //}




            for ( unsigned int j = 0; j < dimension; j++ )
                {

                delta_xg[j] += rVariables.N[i] * rVariables.CurrentDisp(i,j);
                MP_Acceleration[j] += rVariables.N[i] * NodalAcceleration[j];
                MP_AuxVelocity[j] += rVariables.N[i] * NodalVelocity[j];
                //MP_DeltaVelocity[j] += rVariables.N[i] * (NodalVelocity[j]-PreviousNodalVelocity[j]);
                //MP_Acceleration[j] +=NodalInertia[j]/(rVariables.N[i] * MP_Mass * MP_number);//
                //MP_Velocity[j] += NodalMomentum[j]/(rVariables.N[i] * MP_Mass * MP_number);
                //MP_Velocity[j] += DeltaTime * rVariables.N[i] * NodalAcceleration[j];////




                }
            }

        }


        //**************************************************************************************************************************
        //Another way to update the MP velocity (see paper Guilkey and Weiss, 2003)
        MP_Velocity = MP_PreviousVelocity + 0.5 * DeltaTime * (MP_Acceleration + MP_PreviousAcceleration);
        //MP_Acceleration = 2.0/DeltaTime * (MP_Velocity - MP_PreviousVelocity) - MP_PreviousAcceleration;
        //MP_Acceleration = 4/(DeltaTime * DeltaTime) * delta_xg - 4/DeltaTime * MP_PreviousVelocity;
        //MP_Velocity = 2/DeltaTime * delta_xg - MP_PreviousVelocity;
        //MP_Velocity = MP_PreviousVelocity + MP_DeltaVelocity;
        this -> SetValue(MP_VELOCITY,MP_Velocity );
        this -> SetValue(MP_AUXVELOCITY,MP_AuxVelocity );


        const array_1d<double,3>& new_xg = xg + delta_xg ;

        //Update the MP Position
        this -> SetValue(GAUSS_COORD,new_xg);

        //Update the MP Acceleration
        this -> SetValue(MP_ACCELERATION,MP_Acceleration);

        array_1d<double,3>& MP_Displacement = this->GetValue(MP_DISPLACEMENT);

        MP_Displacement += delta_xg;

        //Update the MP Displacement
        this -> SetValue(MP_DISPLACEMENT,MP_Displacement );



        //if (this->Id() == 1518 || this->Id() == 1513)

        //{
            //std::cout<<" MP position "<<this->Id()<<this -> GetValue(GAUSS_COORD)<<std::endl;
            //std::cout<<" delta_xg "<<this->Id()<<delta_xg<<std::endl;

            //std::cout<<" MP_Velocity "<<this->Id()<<this -> GetValue(MP_VELOCITY)<<std::endl;

            //std::cout<<" MP_Acceleration "<<this->Id()<<this -> GetValue(MP_ACCELERATION)<<std::endl;

        //}

        KRATOS_CATCH( "" )
    }

 //************************************************************************************
//************************************************************************************

    void UpdatedLagrangianAxisymmetry::InitializeMaterial()
    {
        KRATOS_TRY
        array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
        GeneralVariables Variables;
        //this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);


        if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
        {

            mConstitutiveLawVector = GetProperties()[CONSTITUTIVE_LAW]->Clone();


            Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);

            mConstitutiveLawVector->InitializeMaterial( GetProperties(), GetGeometry(),
                    Variables.N );

            //}
        }
        else
            KRATOS_THROW_ERROR( std::logic_error, "a constitutive law needs to be specified for the element with ID ", this->Id() )
        //std::cout<< "in initialize material "<<std::endl;
        KRATOS_CATCH( "" )
    }


//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianAxisymmetry::ResetConstitutiveLaw()
    {
        KRATOS_TRY
        array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
        GeneralVariables Variables;
        //this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

        //create and initialize element variables:

        if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
        {

            mConstitutiveLawVector->ResetMaterial( GetProperties(), GetGeometry(), this->MPMShapeFunctionPointValues(Variables.N, xg) );
        }

        KRATOS_CATCH( "" )
    }




//*************************COMPUTE CURRENT DISPLACEMENT*******************************
//************************************************************************************


    Matrix& UpdatedLagrangianAxisymmetry::CalculateCurrentDisp(Matrix & rCurrentDisp, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        rCurrentDisp = zero_matrix<double>( number_of_nodes , dimension);

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
            {

            array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);


            for ( unsigned int j = 0; j < dimension; j++ )
                {

                rCurrentDisp(i,j) = CurrentDisplacement[j];
                }
            }

        return rCurrentDisp;

        KRATOS_CATCH( "" )
    }


//*************************COMPUTE ALMANSI STRAIN*************************************
//************************************************************************************
    void UpdatedLagrangianAxisymmetry::CalculateAlmansiStrain(const Matrix& rF,
        Vector& rStrainVector )
    {
        KRATOS_TRY

        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        //Left Cauchy-Green Calculation
        Matrix LeftCauchyGreen = prod( rF, trans( rF ) );

        //Calculating the inverse of the jacobian
        Matrix InverseLeftCauchyGreen ( dimension, dimension );
        double det_b=0;
        MathUtils<double>::InvertMatrix( LeftCauchyGreen, InverseLeftCauchyGreen, det_b);

        if( dimension == 2 )
        {

            //Almansi Strain Calculation
            rStrainVector[0] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 0, 0 ) );

            rStrainVector[1] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 1, 1 ) );

            rStrainVector[2] = - InverseLeftCauchyGreen( 0, 1 ); // xy

        }
        else if( dimension == 3 )
        {

            //Almansi Strain Calculation
            if ( rStrainVector.size() != 6 ) rStrainVector.resize( 6, false );

            rStrainVector[0] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 0, 0 ) );

            rStrainVector[1] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 1, 1 ) );

            rStrainVector[2] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 2, 2 ) );

            rStrainVector[3] = - InverseLeftCauchyGreen( 0, 1 ); // xy

            rStrainVector[4] = - InverseLeftCauchyGreen( 1, 2 ); // yz

            rStrainVector[5] = - InverseLeftCauchyGreen( 0, 2 ); // xz

        }
        else
        {

            KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

        }


        KRATOS_CATCH( "" )
    }
//*************************COMPUTE GREEN-LAGRANGE STRAIN*************************************
//************************************************************************************

    void UpdatedLagrangianAxisymmetry::CalculateGreenLagrangeStrain(const Matrix& rF,
        Vector& rStrainVector )
    {
        KRATOS_TRY

        const unsigned int dimension  = GetGeometry().WorkingSpaceDimension();

        //Right Cauchy-Green Calculation
        Matrix C ( dimension, dimension );
        noalias( C ) = prod( trans( rF ), rF );

        if( dimension == 2 )
        {

            //Green Lagrange Strain Calculation
            if ( rStrainVector.size() != 3 ) rStrainVector.resize( 3, false );

            rStrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

            rStrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

            rStrainVector[2] = C( 0, 1 ); // xy

        }
        else if( dimension == 3 )
        {

            //Green Lagrange Strain Calculation
            if ( rStrainVector.size() != 6 ) rStrainVector.resize( 6, false );

            rStrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

            rStrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

            rStrainVector[2] = 0.5 * ( C( 2, 2 ) - 1.00 );

            rStrainVector[3] = C( 0, 1 ); // xy

            rStrainVector[4] = C( 1, 2 ); // yz

            rStrainVector[5] = C( 0, 2 ); // xz

        }
        else
        {

            KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" )

        }

        KRATOS_CATCH( "" )
    }



//************************************************************************************
//************************************************************************************

    double& UpdatedLagrangianAxisymmetry::CalculateIntegrationWeight(double& rIntegrationWeight)
    {
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        if( dimension == 2 )
            rIntegrationWeight *= GetProperties()[THICKNESS];

        return rIntegrationWeight;
    }


//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianAxisymmetry::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo )
    {
        int number_of_nodes = GetGeometry().size();
        int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int dim2 = number_of_nodes * dim;

        if ( rResult.size() != dim2 )
            rResult.resize( dim2, false );

        for ( int i = 0; i < number_of_nodes; i++ )
        {
            int index = i * dim;
            rResult[index] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();

            if ( dim == 3 )
                rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
        }

    }

//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianAxisymmetry::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo )
    {
        rElementalDofList.resize( 0 );

        for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
        {
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

            if ( GetGeometry().WorkingSpaceDimension() == 3 )
            {
                rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
            }
        }
        //std::cout<< "ElementalDofList.size() "<<rElementalDofList.size()<<std::endl;
    }



//************************************************************************************
//*******************DAMPING MATRIX***************************************************

    void UpdatedLagrangianAxisymmetry::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        //0.-Initialize the DampingMatrix:
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        //resizing as needed the LHS
        unsigned int MatSize = number_of_nodes * dimension;

        if ( rDampingMatrix.size1() != MatSize )
            rDampingMatrix.resize( MatSize, MatSize, false );

        noalias( rDampingMatrix ) = ZeroMatrix( MatSize, MatSize );


        //1.-Calculate StiffnessMatrix:

        MatrixType StiffnessMatrix  = Matrix();

        this->CalculateLeftHandSide( StiffnessMatrix, rCurrentProcessInfo );

        //2.-Calculate MassMatrix:

        MatrixType MassMatrix  = Matrix();

        this->CalculateMassMatrix ( MassMatrix, rCurrentProcessInfo );


        //3.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
        double alpha = 0;
        if( GetProperties().Has(RAYLEIGH_ALPHA) ){
          alpha = GetProperties()[RAYLEIGH_ALPHA];
        }
        else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) ){
          alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
        }

        double beta  = 0;
        if( GetProperties().Has(RAYLEIGH_BETA) ){
          beta = GetProperties()[RAYLEIGH_BETA];
        }
        else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) ){
          beta = rCurrentProcessInfo[RAYLEIGH_BETA];
        }

        //4.-Compose the Damping Matrix:

        //Rayleigh Damping Matrix: alpha*M + beta*K
        rDampingMatrix  = alpha * MassMatrix;
        rDampingMatrix += beta  * StiffnessMatrix;
        //std::cout<<" rDampingMatrix "<<rDampingMatrix<<std::endl;

        KRATOS_CATCH( "" )
    }
//************************************************************************************
//****************MASS MATRIX*********************************************************

     void UpdatedLagrangianAxisymmetry::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        //I need to call the values of the shape function for the single element
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

        //lumped
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        unsigned int MatSize = dimension * number_of_nodes;

        if ( rMassMatrix.size1() != MatSize )
            rMassMatrix.resize( MatSize, MatSize, false );

        rMassMatrix = ZeroMatrix( MatSize, MatSize );

        double TotalMass = 0;

        //TOTAL MASS OF ONE MP ELEMENT

        TotalMass = this->GetValue(MP_MASS);

        //LUMPED MATRIX

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
        double temp = Variables.N[i] * TotalMass;

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            unsigned int index = i * dimension + j;
            rMassMatrix( index, index ) = temp;
        }
        }

        //CONSISTENT MATRIX
        //for ( unsigned int i = 0; i < number_of_nodes; i++ )
        //{
            //for ( unsigned int j = 0; j < number_of_nodes; j++ )
            //{

                //rMassMatrix( i*2, j*2 ) += Variables.N[i] * Variables.N[j] * TotalMass;
                //rMassMatrix( i * 2 + 1, j * 2 + 1 ) += Variables.N[i] * Variables.N[j] * TotalMass;

            //}
        //}


        //std::cout<<"rMassMatrix "<<rMassMatrix<<std::endl;

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************



    Matrix& UpdatedLagrangianAxisymmetry::MPMJacobian( Matrix& rResult, array_1d<double,3>& rPoint)
    {

        KRATOS_TRY

        //derivatives of shape functions
        Matrix shape_functions_gradients;
        shape_functions_gradients =this->MPMShapeFunctionsLocalGradients(
                                        shape_functions_gradients);
        const GeometryType& rGeom = GetGeometry();

        unsigned int number_nodes = rGeom.PointsNumber();
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        if (dimension ==2)
        {
            rResult.resize( 2, 2);
            rResult = ZeroMatrix(2,2);


            for ( unsigned int i = 0; i < number_nodes; i++ )
            {


                rResult( 0, 0 ) += ( GetGeometry().GetPoint( i ).X() *  shape_functions_gradients( i, 0 ) );
                rResult( 0, 1 ) += ( GetGeometry().GetPoint( i ).X() *  shape_functions_gradients( i, 1 ) );
                rResult( 1, 0 ) += ( GetGeometry().GetPoint( i ).Y() *  shape_functions_gradients( i, 0 ) );
                rResult( 1, 1 ) += ( GetGeometry().GetPoint( i ).Y() *  shape_functions_gradients( i, 1 ) );


            }

        }
        else if(dimension ==3)
        {


            rResult.resize( 3,3);
            rResult = ZeroMatrix(3,3);

            for ( unsigned int i = 0; i < number_nodes; i++ )
            {
                rResult( 0, 0 ) += ( GetGeometry().GetPoint( i ).X() *  shape_functions_gradients( i, 0 ) );
                rResult( 0, 1 ) += ( GetGeometry().GetPoint( i ).X() *  shape_functions_gradients( i, 1 ) );
                rResult( 0, 2 ) += ( GetGeometry().GetPoint( i ).X() *  shape_functions_gradients( i, 2 ) );
                rResult( 1, 0 ) += ( GetGeometry().GetPoint( i ).Y() *  shape_functions_gradients( i, 0 ) );
                rResult( 1, 1 ) += ( GetGeometry().GetPoint( i ).Y() *  shape_functions_gradients( i, 1 ) );
                rResult( 1, 2 ) += ( GetGeometry().GetPoint( i ).Y() *  shape_functions_gradients( i, 2 ) );
                rResult( 2, 0 ) += ( GetGeometry().GetPoint( i ).Z() *  shape_functions_gradients( i, 0 ) );
                rResult( 2, 1 ) += ( GetGeometry().GetPoint( i ).Z() *  shape_functions_gradients( i, 1 ) );
                rResult( 2, 2 ) += ( GetGeometry().GetPoint( i ).Z() *  shape_functions_gradients( i, 2 ) );

            }

        }

        return rResult;

        KRATOS_CATCH( "" )
    }
    /**
       * Jacobian in given point and given a delta position. This method calculate jacobian
       * matrix in given point and a given delta position.
       *
       * @param rPoint point which jacobians has to
    * be calculated in it.
    *
    * @return Matrix of double which is jacobian matrix \f$ J \f$ in given point and a given delta position.
    *
    * @see DeterminantOfJacobian
    * @see InverseOfJacobian
     */
    Matrix& UpdatedLagrangianAxisymmetry::MPMJacobianDelta( Matrix& rResult, array_1d<double,3>& rPoint, Matrix & rDeltaPosition )
    {
        KRATOS_TRY

        Matrix shape_functions_gradients;

        shape_functions_gradients = this->MPMShapeFunctionsLocalGradients(
                                        shape_functions_gradients );



        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        //Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)

        if (dimension ==2)
        {

            rResult.resize( 2, 2);
            rResult = ZeroMatrix(2,2);

            for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
            {
                rResult( 0, 0 ) += ( GetGeometry().GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 0 ) );
                rResult( 0, 1 ) += ( GetGeometry().GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 1 ) );
                rResult( 1, 0 ) += ( GetGeometry().GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 0 ) );
                rResult( 1, 1 ) += ( GetGeometry().GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 1 ) );
            }
        }
        else if(dimension ==3)
        {

            rResult.resize( 3,3);
            rResult = ZeroMatrix(3,3);
            for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
            {
                rResult( 0, 0 ) += ( GetGeometry().GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 0 ) );
                rResult( 0, 1 ) += ( GetGeometry().GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 1 ) );
                rResult( 0, 2 ) += ( GetGeometry().GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 2 ) );
                rResult( 1, 0 ) += ( GetGeometry().GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 0 ) );
                rResult( 1, 1 ) += ( GetGeometry().GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 1 ) );
                rResult( 1, 2 ) += ( GetGeometry().GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 2 ) );
                rResult( 2, 0 ) += ( GetGeometry().GetPoint( i ).Z() + rDeltaPosition(i,2)) * ( shape_functions_gradients( i, 0 ) );
                rResult( 2, 1 ) += ( GetGeometry().GetPoint( i ).Z() + rDeltaPosition(i,2)) * ( shape_functions_gradients( i, 1 ) );
                rResult( 2, 2 ) += ( GetGeometry().GetPoint( i ).Z() + rDeltaPosition(i,2)) * ( shape_functions_gradients( i, 2 ) );
            }
        }



        return rResult;

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    /**
       * Shape function values in given point. This method calculate the shape function
       * vector in given point.
       *
       * @param rPoint point which shape function values have to
    * be calculated in it.
    *
    * @return Vector of double which is shape function vector \f$ N \f$ in given point.
    *
     */
    Vector& UpdatedLagrangianAxisymmetry::MPMShapeFunctionPointValues( Vector& rResult, array_1d<double,3>& rPoint )
    {
        KRATOS_TRY

        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        //Vector rPointLocal = ZeroVector(dimension);
        array_1d<double,3> rPointLocal = ZeroVector(dimension);
        rPointLocal = GetGeometry().PointLocalCoordinates(rPointLocal, rPoint);

        if (dimension == 2 && GetGeometry().PointsNumber() == 3)
        {

            rResult.resize(3, false);

            //1. I evaluate the local coordinates of a point
            //rPointLocal[0] = ((GetGeometry()[2].Coordinates()[1] - GetGeometry()[0].Coordinates()[1])*(rPoint[0] - GetGeometry()[0].Coordinates()[0]) -
                    //(GetGeometry()[2].Coordinates()[0] - GetGeometry()[0].Coordinates()[0])*(rPoint[1] - GetGeometry()[0].Coordinates()[1]))/mDeterminantJ0;


            //rPointLocal[1] = (-(GetGeometry()[1].Coordinates()[1] - GetGeometry()[0].Coordinates()[1])*(rPoint[0] - GetGeometry()[0].Coordinates()[0]) +
                    //(GetGeometry()[1].Coordinates()[0] - GetGeometry()[0].Coordinates()[0])*(rPoint[1] - GetGeometry()[0].Coordinates()[1]))/mDeterminantJ0;

            //2. Shape functions
            rResult( 0 ) = 1 - rPointLocal[0] - rPointLocal[1] ;
            rResult( 1 ) = rPointLocal[0] ;
            rResult( 2 ) = rPointLocal[1];

            //if(this->Id() == 12515)
            //{
				//std::cout<<"CHECK CONNECTIVITY NODES POSITION"<<std::endl;
				//std::cout<<"position x0 "<<GetGeometry()[0].Coordinates()[0]<<std::endl;
				//std::cout<<"position y0 "<<GetGeometry()[0].Coordinates()[1]<<std::endl;
				//std::cout<<"position x1 "<<GetGeometry()[1].Coordinates()[0]<<std::endl;
				//std::cout<<"position y1 "<<GetGeometry()[1].Coordinates()[1]<<std::endl;
				//std::cout<<"position x2 "<<GetGeometry()[2].Coordinates()[0]<<std::endl;
				//std::cout<<"position y2 "<<GetGeometry()[2].Coordinates()[1]<<std::endl;

			//}
			if(rResult(0)*rResult(1)*rResult(2) < 0.0 && rResult(0)+rResult(1)+rResult(2) > 1.0)
			{

				std::cout<<"ELEMENT ID "<<this->Id()<<std::endl;
				std::cout<<"ERROR IN THE EVALUATION OF THE SHAPE FUNCTIONS"<<std::endl;
			}

        }
        else
        {
			rResult.resize(4, false);
			rResult( 0 ) = 0.25 * (1 - rPointLocal[0]) * (1 - rPointLocal[1]) ;
            rResult( 1 ) = 0.25 * (1 + rPointLocal[0]) * (1 - rPointLocal[1]) ;
            rResult( 2 ) = 0.25 * (1 + rPointLocal[0]) * (1 + rPointLocal[1]) ;
            rResult( 3 ) = 0.25 * (1 - rPointLocal[0]) * (1 + rPointLocal[1]) ;

            if(rResult(0)*rResult(1)*rResult(2)*rResult(3) < 0.0 && rResult(0)+rResult(1)+rResult(2)+rResult(3) > 1.0)
			{

				std::cout<<"ELEMENT ID "<<this->Id()<<std::endl;
				std::cout<<"ERROR IN THE EVALUATION OF THE SHAPE FUNCTIONS"<<std::endl;
			}

		}

        //else if (dimension == 3)
        //{
            //rResult.resize(4, false);
            //array_1d<double,3> rPointLocal = ZeroVector(dimension);
            //rPointLocal = GetGeometry().PointLocalCoordinates(rPointLocal, rPoint);
            ////double x10 = GetGeometry()[1].Coordinates()[0]-GetGeometry()[0].Coordinates()[0];
            ////double x20 = GetGeometry()[2].Coordinates()[0]-GetGeometry()[0].Coordinates()[0];
            ////double x30 = GetGeometry()[3].Coordinates()[0]-GetGeometry()[0].Coordinates()[0];
            ////double y10 = GetGeometry()[1].Coordinates()[1]-GetGeometry()[0].Coordinates()[1];
            ////double y20 = GetGeometry()[2].Coordinates()[1]-GetGeometry()[0].Coordinates()[1];
            ////double y30 = GetGeometry()[3].Coordinates()[1]-GetGeometry()[0].Coordinates()[1];
            ////double z10 = GetGeometry()[1].Coordinates()[2]-GetGeometry()[0].Coordinates()[2];
            ////double z20 = GetGeometry()[2].Coordinates()[2]-GetGeometry()[0].Coordinates()[2];
            ////double z30 = GetGeometry()[3].Coordinates()[2]-GetGeometry()[0].Coordinates()[2];

            ////rPointLocal[3] = ((rPoint[0] - GetGeometry()[0].Coordinates()[0])*(y10*z20 - z10*y20) - (rPoint[1] - GetGeometry()[0].Coordinates()[1])*(x10*z20-x20*z10) + (rPoint[2] - GetGeometry()[0].Coordinates()[2])*(y20*x10 - y10*x20))/mDeterminantJ0;

            ////rPointLocal[2] = ((rPoint[0] - GetGeometry()[0].Coordinates()[0])*(y30*z10-y10*z30) + (rPoint[1] - GetGeometry()[0].Coordinates()[1])*(x10*z30-x30*z10) + (rPoint[2] - GetGeometry()[0].Coordinates()[2])*(y10*x30 - y30*x10))/mDeterminantJ0;

            ////rPointLocal[1] = ((rPoint[0] - GetGeometry()[0].Coordinates()[0])*(y20*z30-y30*z20) + (rPoint[1] - GetGeometry()[0].Coordinates()[1])*(x30*z20-x20*z30) + (rPoint[2] - GetGeometry()[0].Coordinates()[2])*(y30*x20 - x30*y20))/mDeterminantJ0;



            ////rPointLocal[0] = 1 - rPointLocal[1] - rPointLocal[2] -rPointLocal[3];


            //rResult( 0 ) =  1.0-(rPointLocal[0]+rPointLocal[1]+rPointLocal[2]) ;
            //rResult( 1 ) = rPointLocal[0] ;
            //rResult( 2 ) = rPointLocal[1];
            //rResult( 3 ) = rPointLocal[2];


        //}

        return rResult;

        KRATOS_CATCH( "" )
    }



    Vector& UpdatedLagrangianAxisymmetry::MPMLocalCoordinates(Vector& rResult, array_1d<double,3>& rPoint)
    {
        KRATOS_TRY

        //Only local coordinated of a point in a tetrahedron is computed
        rResult.resize(4,false);

        double x10 = GetGeometry()[1].Coordinates()[0]-GetGeometry()[0].Coordinates()[0];
        double x20 = GetGeometry()[2].Coordinates()[0]-GetGeometry()[0].Coordinates()[0];
        double x30 = GetGeometry()[3].Coordinates()[0]-GetGeometry()[0].Coordinates()[0];
        double y10 = GetGeometry()[1].Coordinates()[1]-GetGeometry()[0].Coordinates()[1];
        double y20 = GetGeometry()[2].Coordinates()[1]-GetGeometry()[0].Coordinates()[1];
        double y30 = GetGeometry()[3].Coordinates()[1]-GetGeometry()[0].Coordinates()[1];
        double z10 = GetGeometry()[1].Coordinates()[2]-GetGeometry()[0].Coordinates()[2];
        double z20 = GetGeometry()[2].Coordinates()[2]-GetGeometry()[0].Coordinates()[2];
        double z30 = GetGeometry()[3].Coordinates()[2]-GetGeometry()[0].Coordinates()[2];

        rResult[3] = (rPoint[0] - GetGeometry()[0].Coordinates()[0])*(y10*z20 - z10*y20) - (rPoint[1] - GetGeometry()[0].Coordinates()[1])*(x10*z20-x20*z10) + (rPoint[2] - GetGeometry()[0].Coordinates()[2])*(y20*x10 - y10*x20)/mDeterminantJ0;

        rResult[2] = (rPoint[0] - GetGeometry()[0].Coordinates()[0])*(y30*z10-y10*z30) + (rPoint[1] - GetGeometry()[0].Coordinates()[1])*(x10*z30-x30*z10) + (rPoint[2] - GetGeometry()[0].Coordinates()[2])*(y10*x30 - y30*x10)/mDeterminantJ0;

        rResult[1] = (rPoint[0] - GetGeometry()[0].Coordinates()[0])*(y20*z30-y30*z20) + (rPoint[1] - GetGeometry()[0].Coordinates()[1])*(x30*z20-x20*z30) + (rPoint[2] - GetGeometry()[0].Coordinates()[2])*(y30*x20 - x30*y20)/mDeterminantJ0;

        rResult[0] = 1 - rResult[1] - rResult[2] -rResult[3];

        return rResult;

        KRATOS_CATCH( "" )
    }




    Matrix& UpdatedLagrangianAxisymmetry::MPMShapeFunctionsLocalGradients( Matrix& rResult )
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        array_1d<double,3> rPointLocal = ZeroVector(dim);

        array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
        rPointLocal = GetGeometry().PointLocalCoordinates(rPointLocal, xg);


        if (dim == 2 && GetGeometry().PointsNumber() == 3)
        {
            rResult = ZeroMatrix( 3, 2 );
            rResult( 0, 0 ) = -1.0;
            rResult( 0, 1 ) = -1.0;
            rResult( 1, 0 ) =  1.0;
            rResult( 1, 1 ) =  0.0;
            rResult( 2, 0 ) =  0.0;
            rResult( 2, 1 ) =  1.0;
        }
        else
        {

			rResult = ZeroMatrix( 4, 2 );
			rResult( 0, 0 ) = -0.25 * (1 - rPointLocal[1]);
            rResult( 0, 1 ) = -0.25 * (1 - rPointLocal[0]);
            rResult( 1, 0 ) = 0.25 * (1 - rPointLocal[1]);
            rResult( 1, 1 ) = -0.25 * (1 + rPointLocal[0]);
            rResult( 2, 0 ) = 0.25 * (1 + rPointLocal[1]);
            rResult( 2, 1 ) = 0.25 * (1 + rPointLocal[0]);
            rResult( 3, 0 ) = -0.25 * (1 + rPointLocal[1]);
            rResult( 3, 1 ) = 0.25 * (1 - rPointLocal[0]);


		}
        //else if(dim == 3)
        //{
            //rResult = ZeroMatrix( 4, 3 );
            //rResult(0,0) = -1.0;
            //rResult(0,1) = -1.0;
            //rResult(0,2) = -1.0;
            //rResult(1,0) =  1.0;
            //rResult(1,1) =  0.0;
            //rResult(1,2) =  0.0;
            //rResult(2,0) =  0.0;
            //rResult(2,1) =  1.0;
            //rResult(2,2) =  0.0;
            //rResult(3,0) =  0.0;
            //rResult(3,1) =  0.0;
            //rResult(3,2) =  1.0;
        //}

        return rResult;
    }
//************************************************************************************
//************************************************************************************

    //void UpdatedLagrangian::CalculateOnIntegrationPoints( const Variable<double>& rVariable, double& rOutput, ProcessInfo& rCurrentProcessInfo )
    //{

        //rOutput = mConstitutiveLawVector->GetValue( rVariable, rOutput );
    //}

//************************************************************************************
//************************************************************************************

    //void UpdatedLagrangian::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, Vector& rOutput, ProcessInfo& rCurrentProcessInfo )
    //{

        //KRATOS_TRY





        //if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR )
        //{
            ////create and initialize element variables:
            //GeneralVariables Variables;
            //this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

            ////create constitutive law parameters:
            //ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

            ////set constitutive law flags:
            //Flags &ConstitutiveLawOptions=Values.GetOptions();

            //ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN);
            //ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);


            ////compute element kinematics B, F, DN_DX ...
            //this->CalculateKinematics(Variables, rCurrentProcessInfo);

            ////to take in account previous step writing
            //if( mFinalizedStep ){
                //this->GetHistoricalVariables(Variables);
            //}

            ////set general variables to constitutivelaw parameters
            //this->SetGeneralVariables(Variables,Values);

            ////call the constitutive law to update material variables
            //if( rVariable == CAUCHY_STRESS_VECTOR)
                //mConstitutiveLawVector->CalculateMaterialResponseCauchy(Values);
            //else
                //mConstitutiveLawVector->CalculateMaterialResponsePK2(Values);


            //rOutput = Variables.StressVector;


        //}


        //else if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR  || rVariable == ALMANSI_STRAIN_VECTOR )
        //{
            ////create and initialize element variables:
            //GeneralVariables Variables;
            //this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);


            ////compute element kinematics B, F, DN_DX ...
            //this->CalculateKinematics(Variables, rCurrentProcessInfo);

            ////to take in account previous step writing
            //if( mFinalizedStep ){
                //this->GetHistoricalVariables(Variables);
                //Variables.FT = prod(Variables.F,Variables.F0);
            //}

            ////Compute Green-Lagrange Strain
            //if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR )
                //this->CalculateGreenLagrangeStrain( Variables.FT, Variables.StrainVector );
            //else
                //this->CalculateAlmansiStrain( Variables.FT, Variables.StrainVector );

            //if ( rOutput.size() != Variables.StrainVector.size() )
                //rOutput.resize( Variables.StrainVector.size(), false );

            //rOutput = Variables.StrainVector;



        //}
        //else
        //{

            //rOutput = mConstitutiveLawVector->GetValue( rVariable , rOutput );

        //}

        //KRATOS_CATCH( "" )
    //}

//************************************************************************************
//************************************************************************************

    //void UpdatedLagrangian::CalculateOnIntegrationPoints( const Variable<Matrix >& rVariable, Matrix& rOutput, ProcessInfo& rCurrentProcessInfo )
    //{
        //KRATOS_TRY



        //const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();



        //if ( rVariable == CAUCHY_STRESS_TENSOR || rVariable == PK2_STRESS_TENSOR )
        //{
            ////std::vector<Vector> StressVector;
            //Vector StressVector;
            //if( rVariable == CAUCHY_STRESS_TENSOR )
                //this->CalculateOnIntegrationPoints( CAUCHY_STRESS_VECTOR, StressVector, rCurrentProcessInfo );
            //else
                //this->CalculateOnIntegrationPoints( PK2_STRESS_VECTOR, StressVector, rCurrentProcessInfo );


            //if ( rOutput.size2() != dimension )
                //rOutput.resize( dimension, dimension, false );

            //rOutput = MathUtils<double>::StressVectorToTensor(StressVector);


        //}
        //else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR  || rVariable == ALMANSI_STRAIN_TENSOR)
        //{


            //Vector StrainVector;
            //if( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
                //CalculateOnIntegrationPoints( GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );
            //else
                //CalculateOnIntegrationPoints( ALMANSI_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );



            //if ( rOutput.size2() != dimension )
                //rOutput.resize( dimension, dimension, false );

            //rOutput = MathUtils<double>::StrainVectorToTensor(StrainVector);

        //}
        //else if ( rVariable == CONSTITUTIVE_MATRIX )
        //{
            ////create and initialize element variables:
            //GeneralVariables Variables;
            //this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

            ////create constitutive law parameters:
            //ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

            ////set constitutive law flags:
            //Flags &ConstitutiveLawOptions=Values.GetOptions();

            //ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
            ////ConstitutiveLawOptions.Set(ConstitutiveLaw::LAST_KNOWN_CONFIGURATION); //contact domain formulation UL

            ////compute element kinematics B, F, DN_DX ...
            //this->CalculateKinematics(Variables, rCurrentProcessInfo);

            ////set general variables to constitutivelaw parameters
            //this->SetGeneralVariables(Variables,Values);

            ////call the constitutive law to update material variables
            ////mConstitutiveLawVector[PointNumber]->CalculateMaterialResponsePK2(Values); //contact domain formulation UL
            //mConstitutiveLawVector->CalculateMaterialResponseCauchy(Values); //contact domain formulation SL

            //if( rOutput.size2() != Variables.ConstitutiveMatrix.size2() )
                //rOutput.resize( Variables.ConstitutiveMatrix.size1() , Variables.ConstitutiveMatrix.size2() , false );

            //rOutput = Variables.ConstitutiveMatrix;

            ////}


        //}
        //else if ( rVariable == DEFORMATION_GRADIENT )  // VARIABLE SET FOR TRANSFER PURPOUSES
        //{
            ////create and initialize element variables:
            //GeneralVariables Variables;
            //this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);


            ////compute element kinematics B, F, DN_DX ...
            //this->CalculateKinematics(Variables, rCurrentProcessInfo);

            //if( rOutput.size2() != Variables.F.size2() )
                //rOutput.resize( Variables.F.size1() , Variables.F.size2() , false );

            //rOutput = Variables.F;


        //}
        //else
        //{

            //rOutput = mConstitutiveLawVector->GetValue( rVariable , rOutput );

        //}

        //KRATOS_CATCH( "" )
    //}
//*********************************SET DOUBLE VALUE***********************************
//************************************************************************************

    //void UpdatedLagrangian::SetValueOnIntegrationPoints( const Variable<double>& rVariable,
        //double& rValues,
        //ProcessInfo& rCurrentProcessInfo )
    //{
        //if (rVariable == DETERMINANT_F){


            //mDeterminantF0 = rValues;
            //mConstitutiveLawVector->SetValue(rVariable, rValues, rCurrentProcessInfo);


        //}
        //else{

            //UpdatedLagrangian::SetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
        //}
    //}
//********************************SET VECTOR VALUE************************************
//************************************************************************************

    //void UpdatedLagrangian::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, Vector& rValues, ProcessInfo& rCurrentProcessInfo )
    //{

        //mConstitutiveLawVector->SetValue( rVariable, rValues, rCurrentProcessInfo );


    //}


//*******************************SET MATRIX VALUE*************************************
//************************************************************************************

    //void UpdatedLagrangian::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, Matrix& rValues, ProcessInfo& rCurrentProcessInfo )
    //{

        //mConstitutiveLawVector->SetValue( rVariable,
                    //rValues, rCurrentProcessInfo );


    //}
//********************************SET CONSTITUTIVE VALUE******************************
//************************************************************************************

    //void UpdatedLagrangian::SetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,
        //ConstitutiveLaw::Pointer& rValues,
        //ProcessInfo& rCurrentProcessInfo )
    //{
        //if(rVariable == CONSTITUTIVE_LAW)
        //{

            //mConstitutiveLawVector = rValues->Clone();

        //}

        //if(rVariable == CONSTITUTIVE_LAW_POINTER)
        //{

            //mConstitutiveLawVector = rValues;

        //}

    //}
//***************************GET DOUBLE VALUE*****************************************
//************************************************************************************

    //void UpdatedLagrangian::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
            //double& rValues,
            //ProcessInfo& rCurrentProcessInfo )
    //{
        //if (rVariable == DETERMINANT_F){


            //rValues = mDeterminantF0;

        //}
        //else{

            //rValues = mConstitutiveLawVector->GetValue( rVariable, rValues );
        //}
    //}

//**************************GET VECTOR VALUE******************************************
//************************************************************************************

    //void UpdatedLagrangian::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, Vector& rValues, ProcessInfo& rCurrentProcessInfo )
    //{

        //if ( rVariable == PK2_STRESS_VECTOR ||  rVariable == CAUCHY_STRESS_VECTOR )
        //{

            //CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

        //}
        //else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR ||  rVariable == ALMANSI_STRAIN_TENSOR )
        //{

            //CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

        //}
        //else
        //{


            //rValues = mConstitutiveLawVector->GetValue( rVariable, rValues );


        //}


    //}

//************************************************************************************
//************************************************************************************

    //void UpdatedLagrangian::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
            //Matrix& rValues, ProcessInfo& rCurrentProcessInfo )
    //{


        //if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR ||  rVariable == ALMANSI_STRAIN_TENSOR )
        //{
            //CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
        //}

        //if ( rVariable == PK2_STRESS_TENSOR ||  rVariable == CAUCHY_STRESS_TENSOR )
        //{
            //CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
        //}
        //else
        //{

            //rValues = mConstitutiveLawVector->GetValue( rVariable, rValues );

        //}

    //}
//********************************GET CONSTITUTIVE VALUE******************************
//************************************************************************************

    //void UpdatedLagrangian::GetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,
            //ConstitutiveLaw::Pointer& rValues,
            //ProcessInfo& rCurrentProcessInfo )
    //{

        //if(rVariable == CONSTITUTIVE_LAW || rVariable == CONSTITUTIVE_LAW_POINTER)
        //{

            //rValues = mConstitutiveLawVector;

        //}

    //}
//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianAxisymmetry::GetValuesVector( Vector& values, int Step )
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int MatSize = number_of_nodes * dim;

        if ( values.size() != MatSize ) values.resize( MatSize, false );

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = i * dim;
            values[index] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
            values[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

            if ( dim == 3 )
                values[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
        }
    }


//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianAxisymmetry::GetFirstDerivativesVector( Vector& values, int Step )
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int MatSize = number_of_nodes * dim;

        if ( values.size() != MatSize ) values.resize( MatSize, false );

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = i * dim;
            values[index] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
            values[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );

            if ( dim == 3 )
                values[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
        }
    }

//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianAxisymmetry::GetSecondDerivativesVector( Vector& values, int Step )
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int MatSize = number_of_nodes * dim;

        if ( values.size() != MatSize ) values.resize( MatSize, false );

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = i * dim;
            values[index] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
            values[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

            if ( dim == 3 )
                values[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
        }
    }
//************************************************************************************
//************************************************************************************
    void UpdatedLagrangianAxisymmetry::GetHistoricalVariables( GeneralVariables& rVariables )
    {
        //Deformation Gradient F ( set to identity )
        unsigned int size =  rVariables.F.size1();
        rVariables.detF  = 1;
        rVariables.F     = IdentityMatrix(size);

        rVariables.detF0 = mDeterminantF0;
        rVariables.F0    = mDeformationGradientF0;

    }

//************************************************************************************
//************************************************************************************

    //void UpdatedLagrangian::Calculate( const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo )
    //{

        //double lamda = 1.00; // parametro que depende del tipo de problema y del elemento pag 308 libro dinamica de Barbat
        //double c1 = 0.00; //sqrt(GetProperties()[YOUNG_MODULUS]/GetProperties()[DENSITY]); velocidad del sonido en el medio
        //double c2 = 0.00; // norma de la velocidad actual dentro del elemento
        //double c = 0.00;
        //double wmax = 0.00;
        //Vector Values( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() );
        //Vector Velocities;

        //GetFirstDerivativesVector( Velocities, 0 );

        //if ( rVariable == DELTA_TIME )
        //{
            //for ( unsigned int PointNumber = 0;
                    //PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size();
                    //PointNumber++ )
            //{
                //mConstitutiveLawVector[PointNumber]-> GetValue( DELTA_TIME, c1 );
                //Values[PointNumber] = c1;
            //}
        //}

        //c1 = ( *std::max_element( Values.begin(), Values.end() ) );

        //c2 = norm_2( Velocities );

        //c = ( c1 > c2 ) ? c1 : c2;


        //double le = GetGeometry().Length();
        ////KRATOS_WATCH(le)

        ///// maxima frecuencia de un elemento
        //wmax = ( lamda * c ) / le;
        //Output = 2.0 / wmax;
        ////KRATOS_WATCH(Output)

    //}

//************************************************************************************
//************************************************************************************

    //void UpdatedLagrangian::Comprobate_State_Vector( Vector& Result )
    //{
        //for ( unsigned int i = 0.00; i < Result.size(); i++ )
        //{
            //if ( fabs( Result( i ) ) < 1E-9 )
            //{
                //Result( i ) = 0.00;
            //}
        //}
    //}

//*************************DECIMAL CORRECTION OF STRAINS******************************
//************************************************************************************

    void UpdatedLagrangianAxisymmetry::DecimalCorrection(Vector& rVector)
    {
        KRATOS_TRY

        for ( unsigned int i = 0; i < rVector.size(); i++ )
        {
            if( rVector[i]*rVector[i]<1e-24 )
            {
                rVector[i]=0;
            }

        }

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************
/**
 * This function provides the place to perform checks on the completeness of the input.
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 */
    int  UpdatedLagrangianAxisymmetry::Check( const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

        //verify compatibility with the constitutive law
        ConstitutiveLaw::Features LawFeatures;

        this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);


        bool correct_strain_measure = false;
        for(unsigned int i=0; i<LawFeatures.mStrainMeasures.size(); i++)
        {
            if(LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Deformation_Gradient)
                correct_strain_measure = true;
        }

        if( correct_strain_measure == false )
            KRATOS_THROW_ERROR( std::logic_error, "constitutive law is not compatible with the element type ", " Large Displacements " )


        //verify that the variables are correctly initialized

        if ( VELOCITY.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "VELOCITY has Key zero! (check if the application is correctly registered", "" )

        if ( DISPLACEMENT.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "DISPLACEMENT has Key zero! (check if the application is correctly registered", "" )

        if ( ACCELERATION.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "ACCELERATION has Key zero! (check if the application is correctly registered", "" )

        if ( DENSITY.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "DENSITY has Key zero! (check if the application is correctly registered", "" )

        // if ( BODY_FORCE.Key() == 0 )
        //     KRATOS_THROW_ERROR( std::invalid_argument, "BODY_FORCE has Key zero! (check if the application is correctly registered", "" );

        //std::cout << " the variables have been correctly inizialized "<<std::endl;

        //verify that the dofs exist

        for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
        {
            if ( this->GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
                KRATOS_THROW_ERROR( std::invalid_argument, "missing variable DISPLACEMENT on node ", this->GetGeometry()[i].Id() )

            if ( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false )

                KRATOS_THROW_ERROR( std::invalid_argument, "missing one of the dofs for the variable DISPLACEMENT on node ", GetGeometry()[i].Id() )
        }

        //verify that the constitutive law exists
        if ( this->GetProperties().Has( CONSTITUTIVE_LAW ) == false )
        {
            KRATOS_THROW_ERROR( std::logic_error, "constitutive law not provided for property ", this->GetProperties().Id() )
        }

        //Verify that the body force is defined
        // if ( this->GetProperties().Has( BODY_FORCE ) == false )
        // {
        //     KRATOS_THROW_ERROR( std::logic_error, "BODY_FORCE not provided for property ", this->GetProperties().Id() )
        // }

        //verify that the constitutive law has the correct dimension
        if ( dimension == 2 )
        {

            if ( THICKNESS.Key() == 0 )
                KRATOS_THROW_ERROR( std::invalid_argument, "THICKNESS has Key zero! (check if the application is correctly registered", "" )

            if ( this->GetProperties().Has( THICKNESS ) == false )
                KRATOS_THROW_ERROR( std::logic_error, "THICKNESS not provided for element ", this->Id() )
        }
        else
        {
            if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() != 6 )
                KRATOS_THROW_ERROR( std::logic_error, "wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) ", this->Id() )
        }

        //check constitutive law

        if (mConstitutiveLawVector!= 0)
        {
        return mConstitutiveLawVector->Check( GetProperties(), GetGeometry(), rCurrentProcessInfo );
        }

        return 0;

        KRATOS_CATCH( "" );
    }


    double GetPI()
    {
      return std::atan(1.0)*4.0;
    }


    void UpdatedLagrangianAxisymmetry::save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
        //int IntMethod = int(mThisIntegrationMethod);
        //rSerializer.save("IntegrationMethod",IntMethod);
        rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
        rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
        rSerializer.save("DeterminantF0",mDeterminantF0);
        rSerializer.save("InverseJ0",mInverseJ0);
        rSerializer.save("DeterminantJ0",mDeterminantJ0);

    }

    void UpdatedLagrangianAxisymmetry::load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
        rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
        rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
        rSerializer.load("DeterminantF0",mDeterminantF0);
        rSerializer.load("InverseJ0",mInverseJ0);
        rSerializer.load("DeterminantJ0",mDeterminantJ0);
        //int IntMethod;
        //rSerializer.load("IntegrationMethod",IntMethod);
        //mThisIntegrationMethod = IntegrationMethod(IntMethod);


    }





} // Namespace Kratos

