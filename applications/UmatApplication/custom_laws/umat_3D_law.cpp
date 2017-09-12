//
//   Project Name:        KratosUmatApplication      $
//   Created by:          $Author:       JMCarbonell $
//   Last modified by:    $Co-Author:    LlMonforte  $
//   Date:                $Date:      September 2017 $
//   Revision:            $Revision:             0.0 $
//
//


// System includes

// External includes

// Project includes
#include "custom_laws/umat_3D_law.hpp"


/**
 * wrapper function for calling the UMAT fortran subroutine
 * @param STRESS ......... the vector of stresses
 * @param STATEV ......... the vector of state variables
 * @param DDSDDE ......... the material tangent
 * @param SSE ............
 * @param SPD ............
 * @param SCD ............
 * @param RPL ............
 * @param DDSDDT .........
 * @param DRPLDE .........
 * @param DRPLDT .........
 * @param STRAN .......... the vector of total strains
 * @param DSTRAN ......... the vector of incremental strains
 * @param TIME ........... current time
 * @param DTIME .......... current time increment
 * @param TEMP ........... current temperature
 * @param DTEMP .......... current increment of temperature
 * @param PREDEF .........
 * @param DPRED ..........
 * @param MATERL .........
 * @param NDI ............ number of direct strain components (3 in 3D)
 * @param NSHR ........... number if shear strain components (3 in 3D)
 * @param NTENS .......... number of stress components (6 in 3D)
 * @param NSTATV ......... number of state variables (size of STATEV)
 * @param PROPS .......... material parameters
 * @param NPROPS ......... number of material paramters (size of PROPS)
 * @param COORDS .........
 * @param DROT ...........
 * @param PNEWDT .........
 * @param CELENT .........
 * @param DFGRD0 .........
 * @param DFGRD1 .........
 * @param NOEL ...........
 * @param NPT ............ some parameter that is needed by hypoplastic material law
 * @param KSLAY ..........
 * @param KSPT ...........
 * @param KSTEP ..........
 * @param KINC ...........
 * @param MATERIALNUMBER . identifier of the UMAT subroutine to be selected
 */
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

  Umat3DLaw::Umat3DLaw() : ConstitutiveLaw()
  {
  }
  
  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  Umat3DLaw::Umat3DLaw(const Umat3DLaw& rOther) : ConstitutiveLaw(rOther)
  {
  }

  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  Umat3DLaw& Umat3DLaw::operator=(const Umat3DLaw& rOther)
  {
    return *this;
  } 

  
  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer Umat3DLaw::Clone() const
  {
    return ( Umat3DLaw::Pointer(new Umat3DLaw(*this)) );
  }



  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  Umat3DLaw::~Umat3DLaw()
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

  void Umat3DLaw::InitializeMaterial( const Properties& props,
				      const GeometryType& geom,
				      const Vector& ShapeFunctionsValues )
  {

    KRATOS_TRY

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
	STATEV = new double[13]; //[0..5] epsilonElastic, [6..11] epsilonPlastic, [12] alpha
	STATEV_FINALIZED = new double [13];

	for ( unsigned int i = 0; i < 13; i++ ) {
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

  void Umat3DLaw::InitializeSolutionStep( const Properties& props,
					  const GeometryType& geom,
					  const Vector& ShapeFunctionsValues ,
					  const ProcessInfo& CurrentProcessInfo )
  {
  }

  //************************************************************************************
  //************************************************************************************
  
  void Umat3DLaw::FinalizeSolutionStep( const Properties& props,
					const GeometryType& geom,
					const Vector& ShapeFunctionsValues ,
					const ProcessInfo& CurrentProcessInfo )
  {
    KRATOS_TRY

    //does nothing
    for (unsigned int i = 0; i < 6; i++)
      STRESS_FINALIZED[i] = STRESS[i];
    for (int i = 0; i < NSTATV[0]; i++)
      STATEV_FINALIZED[i] = STATEV[i];

    KRATOS_CATCH("")
  }

  //************************************************************************************
  //************************************************************************************

  
  void Umat3DLaw::LoadPreviousInformation()
  {
    KRATOS_TRY

    for (unsigned int i = 0; i < 6; i++)
      STRESS[i] = STRESS_FINALIZED[i];

    for (int i = 0; i < NSTATV[0]; i++)
      STATEV[i] = STATEV_FINALIZED[i];

    KRATOS_CATCH("")
  }

  
  //************************************************************************************
  //************************************************************************************

  void Umat3DLaw::CalculateMaterialResponseCauchy(  Parameters & rValues)
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

    const Properties& MaterialProperties  = rValues.GetMaterialProperties();
    const Vector& rStrainVector            = rValues.GetStrainVector();
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
            DSTRAN[i] = rStrainVector(i) - STRAN_FINALIZED[i];
	  }

	for ( int j = 0; j < NTENS[0]; j++ )
	  {
            DDSDDE[i][j] = 0.0;
	  }
      }

    Matrix AlgorithmicTangent = ZeroMatrix(6,6);
    Vector StressVector = ZeroVector(6);

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

  void Umat3DLaw::FinalizeMaterialResponseCauchy( Parameters & rValues)
  {
    KRATOS_TRY

    const Vector& rStrainVector  = rValues.GetStrainVector();
    for (unsigned int i = 0; i < 6; i++)
      STRAN_FINALIZED[i] = rStrainVector(i);

    KRATOS_CATCH("")
  }


  //*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
  //************************************************************************************

  void Umat3DLaw::GetLawFeatures(Features& rFeatures)
  {
    //Set the type of law
    rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
    rFeatures.mOptions.Set( ISOTROPIC );


    //Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    //Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
  }

  //******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW*************************
  //************************************************************************************

  bool Umat3DLaw::CheckParameters(Parameters& rValues)
  {
    KRATOS_TRY

    return rValues.CheckAllParameters();

    KRATOS_CATCH(" ")    
  }

  //************************************************************************************
  //************************************************************************************


  int Umat3DLaw::Check(const Properties& rMaterialProperties,
		       const GeometryType& rElementGeometry,
		       const ProcessInfo& rCurrentProcessInfo)
  {    
    KRATOS_TRY
      
    return 0;
    
    KRATOS_CATCH(" ")
  }


} // Namespace Kratos
