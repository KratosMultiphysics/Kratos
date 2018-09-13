///    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Zhiming Guo
//                   Riccardo Rossi
//




// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/updated_lagrangian_mpm_element.h"
#include "includes/element.h"
#include "lagrangian_mpm_application_variables.h"

#include "utilities/math_utils.h"

#include "geometries/geometry.h"
//#include "custom_geometries/meshless_geometry.h"

#include <iostream>
#include <fstream>

namespace Kratos
{

//************************************************************************************
//************************************************************************************
UpdatedLagrangianMPMElement::UpdatedLagrangianMPMElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : MeshlessBaseElement(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!

}

//************************************************************************************
//************************************************************************************
UpdatedLagrangianMPMElement::UpdatedLagrangianMPMElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : MeshlessBaseElement(NewId, pGeometry, pProperties)
{
    mFinalizedStep = true;
}

Element::Pointer UpdatedLagrangianMPMElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UpdatedLagrangianMPMElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

Element::Pointer UpdatedLagrangianMPMElement::Create(Element::IndexType NewId,
        Element::GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Element::Pointer(new UpdatedLagrangianMPMElement(NewId, pGeom, pProperties));
    KRATOS_CATCH("");
}

UpdatedLagrangianMPMElement::~UpdatedLagrangianMPMElement()
{
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianMPMElement::Initialize()
{
    KRATOS_TRY


    //const unsigned int dim = GetDomainSize();
    //const unsigned int dim = GetProcessInfo()[DOMAIN_SIZE];
    const unsigned int dim = 2;
    mDeterminantF0 = 1;

    mDeformationGradientF0 = identity_matrix<double> (dim);

    InitializeMaterial();

    KRATOS_CATCH( "" )
}

////************************************************************************************
////************************************************************************************

void UpdatedLagrangianMPMElement::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    GeneralVariables Variables;

    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    mConstitutiveLawVector->InitializeSolutionStep( GetProperties(),
            GetGeometry(), Variables.N, rCurrentProcessInfo );

    mFinalizedStep = false;
    KRATOS_CATCH( "" )
}
////************************************************************************************
////************************************************************************************

void UpdatedLagrangianMPMElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    //const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    //unsigned int voigtsize  = 3;

    //if( dimension == 3 )
    //{
    //voigtsize  = 6;
    //}

    //create and initialize element variables:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    //create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    //set constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

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

    //call the element internal variables update
    this->FinalizeStepVariables(Variables, rCurrentProcessInfo);



    //****************** update the postion of the MPM particle
    array_1d<double,3> delta_displacement = ZeroVector(3);
    array_1d<double,3> previous_displacement = this->GetValue(GAUSS_DISPLACEMENT);
    for(unsigned int i=0; i<GetGeometry().size(); i++)
    {
        delta_displacement += Variables.N[i]*(GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1));

    }
    noalias(this->GetValue(GAUSS_POINT_COORDINATES)) += delta_displacement;
    array_1d<double,3> new_displacement = delta_displacement + previous_displacement;
    this->SetValue(GAUSS_DISPLACEMENT, new_displacement);
    
    
    
    //this->GetValue(GAUSS_POINT_COORDINATES) += delta_displacement;



    mFinalizedStep = true;

    KRATOS_CATCH( "" )
}

////************************************************************************************
////************************************************************************************

void UpdatedLagrangianMPMElement::FinalizeStepVariables( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    //update internal (historical) variables
    mDeterminantF0         = rVariables.detF* rVariables.detF0;
    mDeformationGradientF0 = prod(rVariables.F, rVariables.F0);

}
//**************************************************************************************
void UpdatedLagrangianMPMElement::InitializeMaterial()
{
    KRATOS_TRY;
    //array_1d<double,3>& xg = this->GetValue(GAUSS_POINT_COORDINATES);


    const unsigned int number_of_nodes = GetGeometry().size();
    //const unsigned int dimension = GetDomainSize();//GetGeometry().WorkingSpaceDimension();
    const unsigned int dimension = 2;
    //ProcessInfo& rCurrentProcessInfo = (mpModelPart)->GetProcessInfo();

    GeneralVariables Variables;

    //this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    Variables.N.resize(number_of_nodes);

    Variables.DN_DX.resize(number_of_nodes, dimension); //gradient values at the previous time step

    Variables.integration_weight = 0.0;

    this->GetGeometryData(Variables.integration_weight,Variables.N,Variables.DN_DX);



    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {

        mConstitutiveLawVector = GetProperties()[CONSTITUTIVE_LAW]->Clone();

        mConstitutiveLawVector->InitializeMaterial( GetProperties(), GetGeometry(),
                Variables.N );

    }
    else
    {
        KRATOS_THROW_ERROR( std::logic_error, "a constitutive law needs to be specified for the element with ID ", this->Id() );
        //std::cout<< "in initialize material "<<std::endl;

    }


    KRATOS_CATCH( "" );

}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianMPMElement::ResetConstitutiveLaw()
{
    KRATOS_TRY
    const unsigned int number_of_nodes = GetGeometry().size();
    //const unsigned int dimension = GetDomainSize();//GetGeometry().WorkingSpaceDimension();
    const unsigned int dimension = 2;
    
    //array_1d<double,3>& xg = this->GetValue(GAUSS_COORD_COORDINATES);

    //create and initialize element variables:
      GeneralVariables Variables;
    Variables.N.resize(number_of_nodes);
    Variables.DN_DX.resize(number_of_nodes, dimension); //gradient values at the previous time step
    Variables.integration_weight = 0.0;

    this->GetGeometryData(Variables.integration_weight,Variables.N,Variables.DN_DX);


    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {

        mConstitutiveLawVector->ResetMaterial( GetProperties(), GetGeometry(), Variables.N );
    }

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************
void UpdatedLagrangianMPMElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    //const unsigned int dim = GetDomainSize();
    const unsigned int dim = 2;



    unsigned int matrix_size = number_of_nodes * dim;
    if (rLeftHandSideMatrix.size1() != matrix_size)
        rLeftHandSideMatrix.resize(matrix_size, matrix_size, false);
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( matrix_size, matrix_size ); 
    if (rRightHandSideVector.size() != matrix_size)
        rRightHandSideVector.resize(matrix_size, false);
    rRightHandSideVector = ZeroVector( matrix_size );

    //if (this->Id() == 1)
    //{
		//std::cout<<"rRightHandSideVector in Calculate Local System"<<rRightHandSideVector<<std::endl;
	//}
    CalculateElementalSystem( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo );
    //obtain shape functions

    //obtaine delta_disp and store it in a Matrix

    //Compute deltaF = I + D delta_disp / D xn

    //compute B

    //Compute Ftot

    //Compute material response sigma_cauchy, C

    //compute RHS

    //compute LHS
    

}
void UpdatedLagrangianMPMElement::CalculateElementalSystem( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

	const unsigned int number_of_nodes = GetGeometry().size();
    //const unsigned int dim = GetDomainSize();
    const unsigned int dim = 2;
    unsigned int matrix_size = number_of_nodes * dim;

    //create and initialize element variables:
    GeneralVariables Variables;

    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    //create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);


    //set constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();

    //std::cout<<"in CalculateElementalSystem 5"<<std::endl;
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
//KRATOS_WATCH(__LINE__);

    //auxiliary terms
    Vector VolumeForce;
    
    

    //compute element kinematics B, F, DN_DX ...
    this->CalculateKinematics(Variables,rCurrentProcessInfo);
    
//KRATOS_WATCH(__LINE__);
    //set general variables to constitutivelaw parameters
    this->SetGeneralVariables(Variables,Values);
//KRATOS_WATCH(__LINE__);
    mConstitutiveLawVector->CalculateMaterialResponse(Values, Variables.StressMeasure);

    
    Vector disp0 = ZeroVector(matrix_size);
    Vector disp1 = ZeroVector(matrix_size);
    GetValuesVector(disp0,0);
    GetValuesVector(disp1,1);
    Variables.integration_weight *= Variables.detFT; 



    //update of effective radius
    //double GAUSS_AREA_NEW = Variables.integration_weight;
    //this->GetValue(EFFECTIVE_RADIUS) = sqrt(2*GAUSS_AREA_NEW)/2;
    //this->GetValue(SEARCH_RADIUS) = 4*this->GetValue(EFFECTIVE_RADIUS);

    //KRATOS_WATCH(this->GetValue(EFFECTIVE_RADIUS));

    //FILE *fp;
    std::ofstream outfile;
    outfile.open("1.txt");

    if(this->Id() == 8)
    {
         outfile<<"Element_id:"<<this->Id()<<std::endl;
        //fp = fopen("1.txt","a");
        //fprintf(fp,"%I\n",this->Id());

         //KRATOS_WATCH(Variables.integration_weight);
        for(unsigned int i=0; i<GetGeometry().size(); i++)
        {
            //KRATOS_WATCH(GetGeometry()[i].Id());
            //KRATOS_WATCH(GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT));

            //KRATOS_WATCH(GetGeometry()[i].GetDof(DISPLACEMENT_Y).IsFixed());

            //if(fp!=NULL)
            //{

                //fprintf(fp,"%I",GetGeometry()[i].Id());
                outfile<<"node_id:"<<GetGeometry()[i].Id()<<std::endl;
                outfile<<"Displacement:"<<GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT)<<std::endl;
                outfile<<"ISFixed:"<<GetGeometry()[i].GetDof(DISPLACEMENT_Y).IsFixed()<<std::endl;
                //fprintf(fp,"%d\n",GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT));
                //std::fout<<"IsFixed"<<GetGeometry()[i].GetDof(DISPLACEMENT_Y).IsFixed()<<std::endl;
            //}
        }

        //fclose(fp);
        //outfile.close();
      }


/*
    if (this->Id() == 1)
    {
		//std::cout<<"rVariables.N "<<Variables.N<<std::endl;
		//std::cout<<"rVariables.DN_DX "<<Variables.DN_DX<<std::endl;
		//std::cout<<"rVariables.integration_weight "<<Variables.integration_weight<<std::endl;
		std::cout<<"rVariables.StressVector "<<Variables.StressVector<<std::endl;
		std::cout<<"rVariables.StrainVector "<<Variables.StrainVector<<std::endl;
		std::cout<<"B*GetValuesVector "<< prod(Variables.B,(disp0-disp1))<<std::endl;
		//std::cout<<"rVariables.ConstitutiveMatrix "<<Variables.ConstitutiveMatrix<<std::endl;
    }*/
//KRATOS_WATCH(__LINE__);
    this->CalculateAndAddLHS ( rLeftHandSideMatrix, Variables);
//KRATOS_WATCH(__LINE__);
    VolumeForce  = this->CalculateVolumeForce( VolumeForce, Variables );
    //if (this->Id() == 1)
    //{
		//std::cout<<"VolumeForce "<<VolumeForce<<std::endl;
	//}
//KRATOS_WATCH(__LINE__);
    this->CalculateAndAddRHS ( rRightHandSideVector, Variables, VolumeForce);
//KRATOS_WATCH(__LINE__);
    KRATOS_CATCH( "" )
}

//**********************************CALCULATE LHS************************************************************************************
//***********************************************************************************************************************************
void UpdatedLagrangianMPMElement::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, GeneralVariables& rVariables)
{
    this->CalculateAndAddKuum( rLeftHandSideMatrix, rVariables );
    this->CalculateAndAddKuug( rLeftHandSideMatrix, rVariables );
}

//***********************************MATERIAL MATRIX**********************************
//************************************************************************************

void UpdatedLagrangianMPMElement::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables )
{
    KRATOS_TRY
    
    Matrix tmp = rVariables.integration_weight * prod( rVariables.ConstitutiveMatrix, rVariables.B );
    noalias( rLeftHandSideMatrix ) += prod( trans( rVariables.B ), tmp );

    KRATOS_CATCH( "" )
}
//**********************************GEOMETRICAL MATRIX********************************
//************************************************************************************

void UpdatedLagrangianMPMElement::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables)

{
    KRATOS_TRY

    //unsigned int dimension = GetDomainSize();//GetGeometry().WorkingSpaceDimension();
    const unsigned int dimension = 2;
    Matrix StressTensor = MathUtils<double>::StressVectorToTensor( rVariables.StressVector );
    Matrix tmp = rVariables.integration_weight *  prod( StressTensor, trans( rVariables.DN_DX ) );
    Matrix ReducedKg = prod( rVariables.DN_DX, tmp );
    MathUtils<double>::ExpandAndAddReducedMatrix( rLeftHandSideMatrix, ReducedKg, dimension );

    KRATOS_CATCH( "" )
}

//**********************************CALCULATE RHS************************************************************************************
//***********************************************************************************************************************************
void UpdatedLagrangianMPMElement::CalculateAndAddRHS(VectorType& rRightHandSideVector, GeneralVariables& rVariables, Vector& rVolumeForce)
{
    this->CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce );
    this->CalculateAndAddInternalForces( rRightHandSideVector, rVariables);
}

//*********************************************************************************************************************************
//************************************VOLUME FORCE VECTOR**************************************************************************
Vector& UpdatedLagrangianMPMElement::CalculateVolumeForce( Vector& rVolumeForce, GeneralVariables& rVariables )
{
    KRATOS_TRY

    //const unsigned int number_of_nodes = GetGeometry().size();
    //const unsigned int dimension       = GetDomainSize();//GetGeometry().WorkingSpaceDimension();
    const unsigned int dimension = 2;

    rVolumeForce = ZeroVector(dimension);



    rVolumeForce = this->GetValue(MP_VOLUME_ACCELERATION)* 1.0 / (rVariables.detFT) * GetProperties()[DENSITY];


    return rVolumeForce;

    KRATOS_CATCH( "" )


}

//************************************************************************************
//*********************Calculate the contribution of external force*******************

void UpdatedLagrangianMPMElement::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
        GeneralVariables& rVariables,
        Vector& rVolumeForce)

{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().PointsNumber();

    //unsigned int dimension = GetDomainSize();//GetGeometry().WorkingSpaceDimension();
    const unsigned int dimension = 2;



    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = dimension * i;

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rRightHandSideVector[index + j] += rVariables.N[i] * rVolumeForce[j] * rVariables.integration_weight;

        }

    }

    //KRATOS_WATCH(rVariables.N);
    //KRATOS_WATCH(rVolumeForce);
    //KRATOS_WATCH(rVariables.integration_weight);
//     if (this->Id() == 1)
//     {
// 		std::cout<<"rRightHandSideVector "<<rRightHandSideVector<<std::endl;
// 	}

    KRATOS_CATCH( "" )
}
//************************************************************************************
//*********************Calculate the contribution of internal force*******************

void UpdatedLagrangianMPMElement::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
        GeneralVariables & rVariables)
{
    KRATOS_TRY

    VectorType InternalForces = rVariables.integration_weight * prod( trans( rVariables.B ), rVariables.StressVector );

    noalias( rRightHandSideVector ) -= InternalForces;
    /*
    if (this->Id() == 1)
    {
		std::cout<<"rVariables.integration_weight "<<rVariables.integration_weight<<std::endl;
		std::cout<<"rRightHandSideVector "<<rRightHandSideVector<<std::endl;
    }*/

    KRATOS_CATCH( "" )
}
//**********************************INITIALIZATION OF ELEMENT VARIABLES***************************************************************
//************************************************************************************************************************************
void UpdatedLagrangianMPMElement::InitializeGeneralVariables (GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    //const unsigned int dimension = GetDomainSize();//GetGeometry().WorkingSpaceDimension();
    const unsigned int dimension = 2;
    unsigned int voigtsize  = 3;

    if( dimension == 3 )
    {
        voigtsize  = 6;
    }
    rVariables.detF  = 1;

    rVariables.detF0 = 1;

    rVariables.detFT = 1;

    //rVariables.detJ = 1;

    rVariables.B.resize( voigtsize , number_of_nodes * dimension );

    rVariables.F.resize( dimension, dimension );

    rVariables.F0.resize( dimension, dimension );

    rVariables.FT.resize( dimension, dimension );

    rVariables.ConstitutiveMatrix.resize( voigtsize, voigtsize );

    rVariables.StrainVector.resize( voigtsize );

    rVariables.StressVector.resize( voigtsize );



    rVariables.DN_DX.resize( number_of_nodes, dimension );


    rVariables.N.resize(number_of_nodes);

    rVariables.DN_DX.resize(number_of_nodes, dimension); //gradient values at the previous time step

    rVariables.integration_weight = 0.0;

    this->GetGeometryData(rVariables.integration_weight,rVariables.N,rVariables.DN_DX);
    
    

    //**********************************************************************************************************************

    //CurrentDisp is the variable unknown. It represents the nodal delta displacement. When it is predicted is equal to zero.

    rVariables.DeltaPosition.resize(number_of_nodes, dimension);
    //CalculateDeltaPosition(rVariables.DeltaPosition);



    //std::cout<<"The general variables are initialized"<<std::endl;
    //*************************************************************************************************************************

}


//**************************SET GENERAL CONSTITUTIVE LAW PARAMETERS****************************************************************
//*********************************************************************************************************************************
void UpdatedLagrangianMPMElement::SetGeneralVariables(GeneralVariables& rVariables,
        ConstitutiveLaw::Parameters& rValues)
{
    //Variables.detF is the determinant of the incremental total deformation gradient
    rVariables.detF  = MathUtils<double>::Det(rVariables.F);

    if(rVariables.detF<0.1)
    {

        std::cout<<" Element: "<<this->Id()<<std::endl;
        std::cout<<" Element position "<<this->GetValue(GAUSS_POINT_COORDINATES)<<std::endl;
        unsigned int number_of_nodes = GetGeometry().PointsNumber();
        Geometry< Node<3> > geom = this->GetGeometry();
        KRATOS_WATCH( number_of_nodes );

        array_1d<double, 3> DeltaDisplacement;

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            array_1d<double, 3> &CurrentPosition  = GetGeometry()[i].Coordinates();
            array_1d<double, 3> & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3> & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
            array_1d<double, 3> PreviousPosition  = CurrentPosition - (CurrentDisplacement-PreviousDisplacement);
            std::cout<<" NODE ["<<GetGeometry()[i].Id()<<"]: (Current position: "<<CurrentPosition<<") "<<std::endl;
            //std::cout<<" ---Current Disp: "<<CurrentDisplacement<<" (Previour Disp: "<<PreviousDisplacement<<")"<<std::endl;
            DeltaDisplacement  = CurrentDisplacement-PreviousDisplacement;

        }

        double MaxDelta = DeltaDisplacement[0];

        for ( unsigned int i = 1; i < number_of_nodes; i++ )
        {
            if(DeltaDisplacement[i] > MaxDelta)
            {
                MaxDelta = DeltaDisplacement[i];
            }

        }
        for ( unsigned int i = 1; i < number_of_nodes; i++ )
        {
            if(MaxDelta == DeltaDisplacement[i])
            {
                geom[i].Set(TO_ERASE, true);
            }

        }


        //KRATOS_THROW_ERROR( std::invalid_argument," MPM UPDATED LAGRANGIAN DISPLACEMENT ELEMENT INVERTED: |F|<0  detF = ", rVariables.detF )
    }



    if(rVariables.detF<0)
    {

        std::cout<<" Element: "<<this->Id()<<std::endl;
        std::cout<<" Element position "<<this->GetValue(GAUSS_POINT_COORDINATES)<<std::endl;
        unsigned int number_of_nodes = GetGeometry().PointsNumber();
        KRATOS_WATCH( number_of_nodes );

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            array_1d<double, 3> &CurrentPosition  = GetGeometry()[i].Coordinates();
            array_1d<double, 3> & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3> & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
            array_1d<double, 3> PreviousPosition  = CurrentPosition - (CurrentDisplacement-PreviousDisplacement);
            std::cout<<" NODE ["<<GetGeometry()[i].Id()<<"]: (Current position: "<<CurrentPosition<<") "<<std::endl;
            //std::cout<<" ---Current Disp: "<<CurrentDisplacement<<" (Previour Disp: "<<PreviousDisplacement<<")"<<std::endl;
        }

        //for ( unsigned int i = 0; i < number_of_nodes; i++ )
        //{
            //if( GetGeometry()[i].SolutionStepsDataHas(CONTACT_FORCE) )
            //{
                //array_1d<double, 3 > & PreContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE,1);
                //array_1d<double, 3 > & ContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE);
                //std::cout<<" ---Contact_Force: (Pre:"<<PreContactForce<<", Cur:"<<ContactForce<<") "<<std::endl;
            //}
            //else
            //{
                //std::cout<<" ---Contact_Force: NULL "<<std::endl;
            //}
        //}

        //this->Set(TO_ERASE, true);

        KRATOS_THROW_ERROR( std::invalid_argument," MPM UPDATED LAGRANGIAN DISPLACEMENT ELEMENT INVERTED: |F|<0  detF = ", rVariables.detF )
    }

    rVariables.detFT = rVariables.detF * rVariables.detF0;
    rVariables.FT    = prod( rVariables.F, rVariables.F0 );


    rValues.SetDeterminantF(rVariables.detFT);
    rValues.SetDeformationGradientF(rVariables.FT);
    rValues.SetStrainVector(rVariables.StrainVector);
    rValues.SetStressVector(rVariables.StressVector);
    rValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
    rValues.SetShapeFunctionsDerivatives(rVariables.DN_DX);
    rValues.SetShapeFunctionsValues(rVariables.N);


    //std::cout<<"The general variables are set"<<std::endl;

}



//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void UpdatedLagrangianMPMElement::CalculateKinematics(GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)

{
    KRATOS_TRY

    //const unsigned int dimension = GetDomainSize();//GetGeometry().WorkingSpaceDimension();
    const unsigned int dimension = 2;

    //Define the stress measure
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;




    //Deformation Gradient F [(dx_n+1 - dx_n)/dx_n] to be updated in constitutive law parameter as total deformation gradient
    //the increment of total deformation gradient can be evaluated in 2 ways.
    //1 way.
    //noalias( rVariables.F ) = prod( rVariables.j, InvJ);

    //2 way by means of the gradient of nodal displacement: using this second expression quadratic convergence is not guarantee

    Matrix I=identity_matrix<double>( dimension );

    Matrix GradientDisp = ZeroMatrix(dimension, dimension);
    rVariables.DeltaPosition = CalculateDeltaPosition(rVariables.DeltaPosition);
    GradientDisp = prod(trans(rVariables.DeltaPosition),this->GetValue(SHAPE_FUNCTIONS_DERIVATIVES)); //here I'm using DN_DX of the previous time step

//     KRATOS_WATCH( norm_frobenius(GradientDisp - rVariables.DN_DX) );

    noalias( rVariables.F ) = (I + GradientDisp);

    if (this->Id() == 2)
    {
        //std::cout<< "rVariables.DeltaPosition "<< rVariables.DeltaPosition<<std::endl;
        //std::cout<< "rVariables.N "<< rVariables.N<<std::endl;
        //std::cout<< "rVariables.DN_DX "<< rVariables.DN_DX<<std::endl;
        //std::cout<< "GradientDisp "<< GradientDisp <<std::endl;
        //std::cout<< "rVariables.F "<< rVariables.F <<std::endl;
    }
    Matrix InvF;

    MathUtils<double>::InvertMatrix( rVariables.F, InvF, rVariables.detF );
    //Compute cartesian derivatives [dN/dx_n+1]
    //rVariables.DN_DX = prod( rVariables.DN_DX, InvF); //overwrites DX now is the current position dx
    rVariables.DN_DX = prod( this->GetValue(SHAPE_FUNCTIONS_DERIVATIVES), InvF);
    //Determinant of the Deformation Gradient F_n

    rVariables.detF0 = mDeterminantF0;
    rVariables.F0    = mDeformationGradientF0;

    //KRATOS_WATCH(this->GetGeometry().size());

    //Compute the deformation matrix B
    this->CalculateDeformationMatrix(rVariables.B, rVariables.F, rVariables.DN_DX);


    KRATOS_CATCH( "" )
}
//***********************************************************************************************************
//                                   Deformation Matrix B
//***********************************************************************************************************
void UpdatedLagrangianMPMElement::CalculateDeformationMatrix(Matrix& rB,
        Matrix& rF,
        Matrix& rDN_DX)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    //const unsigned int dimension       = GetDomainSize();//GetGeometry().WorkingSpaceDimension();
    const unsigned int dimension = 2;

    rB.clear(); //set all components to zero

    if( dimension == 2 )
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = 2 * i;

            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 0 ) = rDN_DX( i, 1 );
            rB( 2, index + 1 ) = rDN_DX( i, 0 );

        }
        //if(this->Id() == 365)
        //{
        //std::cout<<"rB "<< this->Id()<< rB<<std::endl;
        //}

    }

    else if( dimension == 3 )
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = 3 * i;

            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 2 ) = rDN_DX( i, 2 );

            rB( 3, index + 0 ) = rDN_DX( i, 1 );
            rB( 3, index + 1 ) = rDN_DX( i, 0 );

            rB( 4, index + 1 ) = rDN_DX( i, 2 );
            rB( 4, index + 2 ) = rDN_DX( i, 1 );

            rB( 5, index + 0 ) = rDN_DX( i, 2 );
            rB( 5, index + 2 ) = rDN_DX( i, 0 );

        }
    }
    else
    {

        KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" )

    }

    KRATOS_CATCH( "" )
}

//*************************COMPUTE CURRENT DISPLACEMENT*******************************
//************************************************************************************


Matrix& UpdatedLagrangianMPMElement::CalculateDeltaPosition(Matrix & rDeltaPosition)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    //unsigned int dimension = GetDomainSize();//GetGeometry().WorkingSpaceDimension();
    const unsigned int dimension = 2;

    rDeltaPosition = zero_matrix<double>( number_of_nodes , dimension);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        const array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
        if (GetGeometry()[i].Id() ==1)
        {
            //std::cout << "PreviousDisplacement "<<PreviousDisplacement<<std::endl;
            //std::cout << "CurrentDisplacement "<<CurrentDisplacement<<std::endl;
        }
        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rDeltaPosition(i,j) = CurrentDisplacement[j]-PreviousDisplacement[j];
        }
    }

    return rDeltaPosition;

    KRATOS_CATCH( "" )
}

//*******************************************************************************************************************
void UpdatedLagrangianMPMElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    //const unsigned int dim = GetDomainSize();
    const unsigned int dim = 2;
    unsigned int matrix_size = number_of_nodes * dim;

    if (rMassMatrix.size1() != matrix_size)
        rMassMatrix.resize(matrix_size, matrix_size, false);

    //set it to zero
    rMassMatrix.clear();
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);
    double density = 1.0 / (Variables.detFT) * GetProperties()[DENSITY];



    //fill the matrix
    for(unsigned int i=0; i<number_of_nodes; i++)
    {
        for(unsigned int j=0; j<number_of_nodes; j++)
        {
            for(unsigned int k=0; k<dim; k++)
            {
                rMassMatrix(i*dim+k, j*dim+k) += density*Variables.integration_weight*Variables.N[i]*Variables.N[j]; //consistent matrix??
            }
        }
    }
//     KRATOS_WATCH(rMassMatrix);

}

//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianMPMElement::Calculate( const Variable<Matrix >& rVariable,Matrix& Output, const ProcessInfo& rCurrentProcessInfo)
    {


        KRATOS_TRY;
        //this is if you need to compute a variable
        const unsigned int number_of_nodes = GetGeometry().size();
        //const unsigned int dim = GetDomainSize();
        const unsigned int dim = 2;


        unsigned int StrainSize;

        if ( dim == 2 )
            StrainSize = 3;
        else
            StrainSize = 6;
        Geometry< Node<3> >& geom = this->GetGeometry();
        //KRATOS_WATCH("BBBBBBBBBBBBBBBBBBB");

        //create and initialize element variables:
        GeneralVariables Variables;

        this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

        //create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);


        //set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();

        //std::cout<<"in CalculateElementalSystem 5"<<std::endl;
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);



        //compute element kinematics B, F, DN_DX ...
        this->CalculateKinematics(Variables,rCurrentProcessInfo);

    //KRATOS_WATCH(__LINE__);
        //set general variables to constitutivelaw parameters
        this->SetGeneralVariables(Variables,Values);
    //KRATOS_WATCH(__LINE__);
        Matrix I=identity_matrix<double>( dim );
		Matrix InvI;
		double det;
        MathUtils<double>::InvertMatrix( I, InvI, det );
    
		Variables.FT = prod(mDeformationGradientF0, InvI);
    
		Variables.detFT = mDeterminantF0 / det;
    
        Values.SetDeterminantF(Variables.detFT);
        Values.SetDeformationGradientF(Variables.FT);
        mConstitutiveLawVector->CalculateMaterialResponse(Values, Variables.StressMeasure);

         //displacements Vector values;
         //boost::numeric::ublas::vector<double> values;

         //unsigned int MatSize = number_of_nodes * dim;


         //Matrix Output0;

        if (rVariable == NODAL_STRESS)
        {
           //for(Geometry< Node<3> >::iterator inode = rGeom.begin(); inode!=rGeom.end(); inode++)
           for(unsigned int i = 0; i < geom.size(); i++)
           {
               geom[i].FastGetSolutionStepValue(NODAL_STRESS).resize(1,StrainSize,false);



               if ( Output.size2() != StrainSize ){
                   Output.resize(1,StrainSize,false);}


                for(unsigned int k = 0; k < StrainSize; k++)
                {

                    Output(0,k) = Variables.integration_weight*Variables.N[i]*Variables.StressVector[k];
                   //row(Output,i) = A*Ng[i]*StressVector[k];

                }


                geom[i].FastGetSolutionStepValue(NODAL_STRESS) += Output;

               //std::cout << geom[i].FastGetSolutionStepValue(NODAL_STRESS) << std::endl;
            }



         }
        else if (rVariable == NODAL_STRAIN)
        {


           //Matrix Output1;
           for(unsigned int i = 0; i < geom.size(); i++)
           {
               geom[i].FastGetSolutionStepValue(NODAL_STRAIN).resize(1,StrainSize,false);
               //geom[i].FastGetSolutionStepValue(NODAL_STRESS) = ZeroMatrix(1,StrainSize);


               if ( Output.size2() != StrainSize ){
                    Output.resize(1,StrainSize,false);
                    Output = ZeroMatrix(1,StrainSize);}

                for(unsigned int k = 0; k < StrainSize; k++)
                {

                    Output(0,k) = Variables.integration_weight*Variables.N[i]*Variables.StrainVector[k];


                }


                geom[i].FastGetSolutionStepValue(NODAL_STRAIN) += Output;
                //KRATOS_WATCH(geom[i].FastGetSolutionStepValue(NODAL_STRAIN));

           }

        }


         KRATOS_CATCH("");

    }

    void UpdatedLagrangianMPMElement::Calculate( const Variable<double >& rVariable,double& Output, const ProcessInfo& rCurrentProcessInfo)
    {


        KRATOS_TRY;
        //this is if you need to compute a variable
        const unsigned int number_of_nodes = GetGeometry().size();
        //const unsigned int dim = GetDomainSize();
        const unsigned int dim = 2;
        unsigned int matrix_size = number_of_nodes * dim;

        //create and initialize element variables:
        GeneralVariables Variables;

        this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

        //create constitutive law parameters:
        //ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);


        ////set constitutive law flags:
        //Flags &ConstitutiveLawOptions=Values.GetOptions();

        ////std::cout<<"in CalculateElementalSystem 5"<<std::endl;
        //ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

        //ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

        //ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);


        ////compute element kinematics B, F, DN_DX ...
        //this->CalculateKinematics(Variables,rCurrentProcessInfo);

    ////KRATOS_WATCH(__LINE__);
        ////set general variables to constitutivelaw parameters
        //this->SetGeneralVariables(Variables,Values);
    ////KRATOS_WATCH(__LINE__);
        //mConstitutiveLawVector->CalculateMaterialResponse(Values, Variables.StressMeasure);



         Geometry< Node<3> >& geom = this->GetGeometry();

         if (rVariable == NODAL_AREA)
         {


           //for(Geometry< Node<3> >::iterator inode = rGeom.begin(); inode!=rGeom.end(); inode++)
           for(unsigned int i = 0; i < geom.size(); i++)
           {

               //geom[i].FastGetSolutionStepValue(NODAL_AREA) = 0;

           //for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
           //{
               //GetGeometry()[i].FastGetSolutionStepValue(NODAL_STRESS) += A*Ng[i]*StressVector;
               //GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += A*Ng[i];



               double aux = Variables.integration_weight*Variables.N[i];
               //KRATOS_WATCH(Variables.N[i]);

               //Output += aux;

               geom[i].FastGetSolutionStepValue(NODAL_AREA) += aux;
               //std::cout << aux << std::endl;

               //std::cout << mTotalDomainSize << std::endl;
               //std::cout << Ng[i] << std::endl;


               //KRATOS_WATCH(geom[i].FastGetSolutionStepValue(NODAL_AREA));



            }



         }

/*
         else if ( rVariable==NODE_EQUIVALENT_PLASTIC_STRAIN )
         {
             //rVariable = PLASTIC_STRAIN;
             //Output = mConstitutiveLawVector[0]->GetValue( PLASTIC_STRAIN, Output );


             for(unsigned int i = 0; i < geom.size(); i++)
             {

                 double aux = EQ_PLASTIC_STRAIN*Variables.N[i];
                 //KRATOS_WATCH(aux);

                 geom[i].FastGetSolutionStepValue(NODE_EQUIVALENT_PLASTIC_STRAIN) += aux;


              }

         }




         else if ( rVariable == DAMAGE )
         {

              Output = mConstitutiveLawVector[0]-> GetValue( DAMAGE,Output);

              for(unsigned int i = 0; i < geom.size(); i++)
              {


                  //Output = Output*Ng[i];
                  geom[i].FastGetSolutionStepValue(DAMAGE) += Output;
                  //std::cout << aux << std::endl;

                  //KRATOS_WATCH(geom[i].FastGetSolutionStepValue(DAMAGE));

               }




         }
*/

         KRATOS_CATCH("");

    }



//*********************************************************************************************************************
void UpdatedLagrangianMPMElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    //const unsigned int dim = GetDomainSize();
    const unsigned int dim = 2;
    unsigned int matrix_size = number_of_nodes * dim;

    if ( rResult.size() != matrix_size )
        rResult.resize( matrix_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
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

void UpdatedLagrangianMPMElement::GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo )
{
    ElementalDofList.resize( 0 );
    const unsigned int dim = 2;

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

        if ( dim == 3 )
        {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
        }
    }
}

void UpdatedLagrangianMPMElement::GetValuesVector( Vector& values, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    //const unsigned int dim = GetDomainSize();
    const unsigned int dim = 2;
    unsigned int matrix_size = number_of_nodes * dim;

    if ( values.size() != matrix_size ) values.resize( matrix_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dim;
        const array_1d<double,3>& disp = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT, Step );
        values[index] = disp[0];
        values[index + 1] = disp[1];

        if ( dim == 3 )
            values[index + 2] = disp[2];
    }
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianMPMElement::GetFirstDerivativesVector( Vector& values, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    //const unsigned int dim = GetDomainSize();
    const unsigned int dim = 2;
    unsigned int matrix_size = number_of_nodes * dim;

    if ( values.size() != matrix_size ) values.resize( matrix_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dim;
        const array_1d<double,3>& vel = GetGeometry()[i].GetSolutionStepValue( VELOCITY, Step );
        values[index] = vel[0];
        values[index + 1] = vel[1];

        if ( dim == 3 )
            values[index + 2] = vel[2];
    }
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianMPMElement::GetSecondDerivativesVector( Vector& values, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    //const unsigned int dim = GetDomainSize();
    const unsigned int dim = 2;
    unsigned int matrix_size = number_of_nodes * dim;

    if ( values.size() != matrix_size ) values.resize( matrix_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dim;
        const array_1d<double,3>& acc = GetGeometry()[i].GetSolutionStepValue( ACCELERATION, Step );
        values[index] = acc[0];
        values[index + 1] = acc[1];

        if ( dim == 3 )
            values[index + 2] = acc[2];
    }
}

//************************************************************************************
//************************************************************************************

//************************************************************************************
//************************************************************************************
/**
 * This function provides the place to perform checks on the completeness of the input.
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 */
int  UpdatedLagrangianMPMElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //unsigned int dimension = GetDomainSize();
    //const unsigned int dimension = GetProcessInfo()[DOMAIN_SIZE];
    const unsigned int dimension = 2;

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

void UpdatedLagrangianMPMElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    //int IntMethod = int(mThisIntegrationMethod);
    //rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.save("DeterminantF0",mDeterminantF0);
    rSerializer.save("FinalizedStep",mFinalizedStep);
    //rSerializer.save("InverseJ0",mInverseJ0);
    //rSerializer.save("DeterminantJ0",mDeterminantJ0);

}

void UpdatedLagrangianMPMElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.load("DeterminantF0",mDeterminantF0);
    rSerializer.load("FinalizedStep",mFinalizedStep);
    //rSerializer.load("InverseJ0",mInverseJ0);
    //rSerializer.load("DeterminantJ0",mDeterminantJ0);


}









} // Namespace Kratos


