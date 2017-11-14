//
//   Project Name:        KratosUmatApplication      $
//   Created by:          $Author:       LlMonforte  $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:      September 2017 $
//   Revision:            $Revision:             0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_laws/large_strain_umat_3D_law.hpp"


extern "C" void umat_wrapper_( double* STRESS, double* STATEV, double** DDSDDE, double* SSE, double* SPD, double* SCD,
			       double* RPL, double* DDSDDT, double* DRPLDE, double* DRPLDT, double* STRAN, double* DSTRAN,
			       double* TIME, double* DTIME, double* TEMP, double* DTEMP, double* PREDEF, double* DPRED,
			       char* MATERL, int* NDI, int* NSHR, int* NTENS, int* NSTATV, double* PROPS, int* NPROPS,
			       double* COORDS, double** DROT, double* PNEWDT, double* CELENT, double** DFGRD0,
			       double** DFGRD1, double* NOEL, int* NPT, double* KSLAY, double* KSPT, double* KSTEP,
			       double* KINC, int* MATERIALNUMBER );

namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  LargeStrainUmat3DLaw::LargeStrainUmat3DLaw() : Umat3DLaw()
  {
  }
  
  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  LargeStrainUmat3DLaw::LargeStrainUmat3DLaw(const LargeStrainUmat3DLaw& rOther) : Umat3DLaw(rOther)
  {
  }

  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  LargeStrainUmat3DLaw& LargeStrainUmat3DLaw::operator=(const LargeStrainUmat3DLaw& rOther)
  {
    return *this;
  } 

  
  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer LargeStrainUmat3DLaw::Clone() const
  {
    return ( LargeStrainUmat3DLaw::Pointer(new LargeStrainUmat3DLaw(*this)) );
  }



  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  LargeStrainUmat3DLaw::~LargeStrainUmat3DLaw()
  {
  }

  //*******************************OPERATIONS FROM BASE CLASS***************************
  //************************************************************************************

  // LMV: Some  variables do not need to be saved: for instance NDI; NSHR, NTENS. 
  // Additionally, looking at mises_umat it can be seen that NDI == 3,
  // since in the calculation of J2 all the three stress components are used.
  // I think it is better to use the philosophy of ConstitutiveModelsApplication:
  // with laws and models (this is a model)

  //************************************************************************************
  //************************************************************************************

  
  void LargeStrainUmat3DLaw::InitializeMaterial( const Properties& props,
						 const GeometryType& geom,
						 const Vector& ShapeFunctionsValues )
  {

    KRATOS_TRY

    mPreviousDeformationGradient = identity_matrix<double>(3);


    // 1. Define sizes of vectors
    int NDI[1]; NDI[0] = 3;
    int NSHR[1]; NSHR[0] = 3;
    int NTENS[1]; NTENS[0] = 6;

    //umat variables initialisation used by mises:
    STRESS = new double[6];


    // Finalized stress
    STRESS_FINALIZED = new double[6];
    STRAN_FINALIZED  = new double[6];

    for ( int i = 0; i < NTENS[0]; i++ )
      {
	STRESS[i] = 0.0 ;
	STRESS_FINALIZED[i] = 0.0;
	STRAN_FINALIZED[i] = 0.0;
      }

    // size of internal variables and material properties
    NSTATV = new int[1];
    NPROPS = new int[1];

    MaterialNumber = new int[1];

    // LMV: Define the number of material properties and allocate the vector
    NPROPS[0] = 2;
    PROPS = new double[NPROPS[0]];

    for (unsigned int i = 0; i < NPROPS[0]; i++)
      PROPS[i] = -1.0;

    // LMV: Define the material
    MaterialNumber[0] = ( int ) 0;

    switch ( MaterialNumber[0] )
      {

      case 0:
	//linearElastic material, mises umat, 2 matprops

	if ( NPROPS[0] != 2 )
	  KRATOS_THROW_ERROR( std::logic_error, "LinearElastic umat material number material properties failure must be 2 ", "" );

	NSTATV[0] = 13;
	STATEV = new double [ NSTATV[0] ]; //[0..5] epsilonElastic, [6..11] epsilonPlastic, [12] alpha
	STATEV_FINALIZED = new double [ NSTATV[0] ];

	for ( unsigned int i = 0; i < NSTATV[0]; i++ ) {
	  STATEV[i] = 0.0;
	  STATEV_FINALIZED[i] = 0.0;
	}

	break;


      case 1:
	//mises material, mises umat,>=4 matprops
	if ( NPROPS[0] < 4 )
	  KRATOS_THROW_ERROR( std::logic_error, "Mises umat material number material properties failure must be >=4, E,nu,Syield,EPlasticYield ", "" );

	STATEV = new double[13]; //[0..5] epsilonElastic, [6..11] epsilonPlastic, [12] alpha

	for ( unsigned int i = 0; i < 13; i++ )
	  STATEV[i] = 0.0;

	NSTATV[0] = 13;

	break;

      case 2:
	//hypoplastic material with small-strain stiffness
	STATEV = new double[14];
	for ( unsigned int i = 0; i < 14; i++ )
	  STATEV[i] = 0.0;

	NSTATV[0] = 14;

	break;


      default:
	std::cout << "No umat material with id: " << MaterialNumber[0] << " defined" << std::endl;

	KRATOS_THROW_ERROR( std::logic_error, "switch umat material error", "" );
      }


    KRATOS_CATCH("")
 }


  //************************************************************************************
  //************************************************************************************


  void LargeStrainUmat3DLaw::CalculateMaterialResponseKirchhoff(  Parameters & rValues)
  {

    KRATOS_TRY

    // 1. Define sizes of vectors
    int NDI[1]; NDI[0] = 3;
    int NSHR[1]; NSHR[0] = 3;
    int NTENS[1]; NTENS[0] = 6;


    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const ProcessInfo & rCurrentProcessInfo = rValues.GetProcessInfo();
    const Matrix& DeformationGradientF     = rValues.GetDeformationGradientF();

    const Properties& MaterialProperties   = rValues.GetMaterialProperties();
    //const Vector& rStrainVector            = rValues.GetStrainVector();
    Vector& rStressVector                  = rValues.GetStressVector();
    Matrix& rConstitutiveMatrix            = rValues.GetConstitutiveMatrix();


    if ( PROPS[0] < 0.0) {

      if ( MaterialProperties.Has(YIELD_STRESS) ) {
	if ( MaterialProperties[YIELD_STRESS]  > 0.0) {
	  NPROPS[0] = 3; // Mises
	  PROPS = new double[NPROPS[0]];
	  PROPS[2] = MaterialProperties[YIELD_STRESS];
	  //PROPS[3] = 0.1;
	}
	PROPS[0] = MaterialProperties[YOUNG_MODULUS];
	PROPS[1] = MaterialProperties[POISSON_RATIO];
      }
    }


    double DDSDDE[NTENS[0]][NTENS[0]];

    // information to remove?
    double TIM[2];
    double DTIME[1];
    TIM[0] = rCurrentProcessInfo[TIME];
    TIM[1] = rCurrentProcessInfo[DELTA_TIME];
    DTIME[0] = rCurrentProcessInfo[DELTA_TIME];
    int NPT[1]; NPT[0] = 0;

    // LOAD PREVIOUS INFORMATION
    LoadPreviousInformation();

    double DSTRAN[6];
    for ( int i = 0; i < NTENS[0]; i++ )
      {

	// DeltaEpsilon = Epsilon_n+1 - Epsilon_n^e - Epsilon_n^p
	// if you have the umat defined in this manner. xD
	if ( MaterialNumber[0] == 0 || MaterialNumber[0] == 1 )
	  {
            //DSTRAN[i] = rStrainVector[i] - STATEV_FINALIZED[i] - STATEV_FINALIZED[i+6];
            //DSTRAN[i] = rStrainVector(i) - STRAN_FINALIZED[i];
	  }

	for ( int j = 0; j < NTENS[0]; j++ )
	  {
            DDSDDE[i][j] = 0.0;
	  }
      }

    // Compute f_n^n+1
    Matrix RelativeDeformationGradient = ZeroMatrix(3);

    double det;
    MathUtils<double>::InvertMatrix( mPreviousDeformationGradient, RelativeDeformationGradient, det);
    RelativeDeformationGradient = prod( DeformationGradientF, RelativeDeformationGradient);

    // project previous stress
    Matrix PreviousStressTensor = ZeroMatrix(3);
    Vector StressVectorAUX = ZeroVector(6);
    for (unsigned int i = 0; i < 6; i++)
      StressVectorAUX(i) = STRESS_FINALIZED[i];
    PreviousStressTensor = MathUtils<double>::StressVectorToTensor( StressVectorAUX);
    PreviousStressTensor = prod( RelativeDeformationGradient, Matrix( prod( PreviousStressTensor, trans(RelativeDeformationGradient) ) ));
    StressVectorAUX = MathUtils<double>::StressTensorToVector( PreviousStressTensor, 6 );
    for (unsigned int i = 0; i < 6; i++)
      STRESS[i] = StressVectorAUX(i);

    // compute delta strain
    Matrix AuxMatrix = ZeroMatrix(3);
    MathUtils<double>::InvertMatrix( RelativeDeformationGradient, AuxMatrix, det);
    Matrix StrainMatrix = identity_matrix<double> (3);

    StrainMatrix -= prod( trans(AuxMatrix), AuxMatrix);
    Vector StrainVectorAux = ZeroVector(6);
    StrainVectorAux = MathUtils<double>::StrainTensorToVector( StrainMatrix, 6);
    for (unsigned int i = 0; i < 6;i++)
      DSTRAN[i] = StrainVectorAux(i) / 2.0;

    // NOTE: parameters that are not required by the umat implementations used so far are given as NULL pointers
    // if any new umat is implemented, please check the required parameters and add them accordingly
    // make sure that for backward compatibility the new parameters are initialized as NULL pointers for all
    // other umat materials
    double* STRAN;
    STRAN = NULL;
    double SPD[1];

    umat_wrapper_( STRESS, STATEV, ( double** ) DDSDDE, NULL, SPD, NULL, NULL, NULL, NULL, NULL, STRAN, DSTRAN,
		   TIM, DTIME, NULL, NULL, NULL, NULL, NULL, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,
		   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NPT, NULL, NULL, NULL, NULL, MaterialNumber );


    //copy variables back
    Matrix AlgorithmicTangent = ZeroMatrix(6,6);
    Vector StressVector = ZeroVector(6);
    for ( int i = 0; i < NTENS[0]; i++ )
      {
	for ( int j = 0; j < NTENS[0]; j++ )
	  {
            AlgorithmicTangent( i, j ) = DDSDDE[i][j];
	  }

	StressVector[i] = STRESS[i];
      }


    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
      {
	noalias( rStressVector ) = StressVector;
      }
    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
      {
	noalias( rConstitutiveMatrix ) = AlgorithmicTangent; 
      }

    KRATOS_CATCH("")

  }


  //***********************************FINALIZE*****************************************
  //************************************************************************************

  void LargeStrainUmat3DLaw::FinalizeMaterialResponseCauchy( Parameters & rValues)
  {
    KRATOS_TRY

    const Vector& rStrainVector  = rValues.GetStrainVector();
    for (unsigned int i = 0; i < 6; i++)
      STRAN_FINALIZED[i] = rStrainVector(i);

    const Matrix& DeformationGradientF  = rValues.GetDeformationGradientF();
    mPreviousDeformationGradient = DeformationGradientF;
    KRATOS_CATCH("")
  }

  //*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
  //************************************************************************************

  void LargeStrainUmat3DLaw::GetLawFeatures(Features& rFeatures)
  {
    //Set the type of law
    rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
    rFeatures.mOptions.Set( FINITE_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    //Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();

  }


} // Namespace Kratos
