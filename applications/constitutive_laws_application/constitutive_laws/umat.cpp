#include <iostream>

#include "includes/define.h"
#include "constitutive_laws/umat.h"


//#include "includes/constitutive_law.h"

//#include "includes/variables.h"
//#include "includes/process_info.h"
#include "includes/properties.h"
//#include "geometries/geometry.h"
//#include "constitutive_laws_application.h"
//#include "includes/ublas_interface.h"



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

   //typedef:
   typedef Geometry<Node < 3 > > GeometryType;


   //constructor:
   Umat::Umat() : HyperElastic3DLaw() {}

   //destructor:
   Umat::~Umat()
   {
   }


   //Vale, el problema és que això és hypo-elastic i estic finalitzant la solució a cada iteració (amb les variables aquestes rares i tal). Per això no puc reduir l'error aquell petit que em surt. I a sobre encara tinc sort que és petit.
   //


   //init material:
   void Umat::InitializeMaterial( const Properties& props,
         const GeometryType& geom,
         const Vector& ShapeFunctionsValues )
   {

      KRATOS_TRY

      for (unsigned int i = 0; i < 30; i++)
      //umat variables initialisation used by mises:
      NTENS = new int[1];
      NTENS[0] = 6;
      NDI = new int[1];
      NDI[0] = 3;
      NSHR = new int[1];
      NSHR[0] = 3;

      STRESS = new double[6];
      STRAN = new double[6];
      DSTRAN = new double[6];

      // Finalized stress and state variables
      STRESS_END = new double[6];

      for ( int i = 0; i < NTENS[0]; i++ )
      {
         STRESS[i] = 0.0 ;
         STRESS_END[i] = 0.0;
         STRAN[i] = 0.0 ;
         DSTRAN[i] = 0.0 ;
      }

      NSTATV = new int[1];

      NPROPS = new int[1];

      MaterialNumber = new int[1];

      // material data should be readed here. I just set the number for the moment //LMV
      //Vector mdata = props[MATERIAL_PARAMETERS] ;
      //NPROPS[0] = mdata.size() - 1;
      NPROPS[0] = 2;


      PROPS = new double[NPROPS[0]];

      //choosing right Law and deleting first mdata entry, creating PROPS array:
      //MaterialNumber[0] = ( int ) mdata[0];
      MaterialNumber[0] = ( int ) 0;

      /*for ( int i = 1; i < NPROPS[0] + 1; i++ )
      {
         PROPS[i-1] = mdata[i];
      }*/


      switch ( MaterialNumber[0] )
      {

         case 0:
            //linearElastic material, mises umat, 2 matprops

            if ( NPROPS[0] != 2 )
               KRATOS_THROW_ERROR( std::logic_error, "LinearElastic umat material number material properties failure must be 2 ", "" );

            STATEV = new double[13]; //[0..5] epsilonElastic, [6..11] epsilonPlastic, [12] alpha
            STATEV_END = new double [13];

            for ( unsigned int i = 0; i < 13; i++ ) {
               STATEV[i] = 0.0;
               STATEV_END[i] = 0.0;
            }

            NSTATV[0] = 13;

            // LMV!!!!
            PROPS[0] = 1000.0; // young??
            PROPS[1] = 0.3; // nu ??
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

   void Umat::InitializeSolutionStep( const Properties& props,
         const GeometryType& geom,
         const Vector& ShapeFunctionsValues ,
         const ProcessInfo& CurrentProcessInfo )
   {
      //does nothing
   }

   /*void Umat::CalculateMaterialResponse( const Vector& StrainVector,
     const Matrix& DeformationGradient,
     Vector& StressVector,
     Matrix& AlgorithmicTangent,
     const ProcessInfo& CurrentProcessInfo,
     const Properties& props,
     const GeometryType& geom,
     const Vector& ShapeFunctionsValues,
     bool CalculateStresses,
     int CalculateTangent,
     bool SaveInternalVariables )*/

   void Umat::LoadPreviousInformation()
   {
      KRATOS_TRY

      for (unsigned int i = 0; i < 6; i++)
         STRESS[i] = STRESS_END[i];
      for (unsigned int i = 0; i < NSTATV[0]; i++)
         STATEV[i] = STATEV_END[i];

      KRATOS_CATCH("")
   }
   void Umat::CalculateMaterialResponseKirchhoff(  Parameters & rValues)
   {


      KRATOS_TRY

      //a.-Check if the constitutive parameters are passed correctly to the law calculation
      CheckParameters(rValues);

      //b.- Get Values to compute the constitutive law:
      Flags &Options=rValues.GetOptions();


      const ProcessInfo & rCurrentProcessInfo = rValues.GetProcessInfo();
      const Properties& MaterialProperties  = rValues.GetMaterialProperties();

      Vector& rStrainVector                  = rValues.GetStrainVector();
      Vector& rStressVector                  = rValues.GetStressVector();
      Matrix& rConstitutiveMatrix            = rValues.GetConstitutiveMatrix();



      Vector StrainVector(6);
      noalias( StrainVector) = rStrainVector;


      double DDSDDE[NTENS[0]][NTENS[0]];
      double DDSDDT[NTENS[0]];
      double DRPLDE[NTENS[0]];
      double TIM[2];
      double DTIME[1];
      int NPT[1];
      TIM[0] = rCurrentProcessInfo[TIME];
      TIM[1] = rCurrentProcessInfo[DELTA_TIME];
      DTIME[0] = rCurrentProcessInfo[DELTA_TIME];
      NPT[0] = 0;

      // LOAD PREVIOUS INFORMATION
      LoadPreviousInformation();

      for ( int i = 0; i < NTENS[0]; i++ )
      {
         STRAN[i] = StrainVector[i];

         //deltaEpsilon = Epsilon - STATEV  in case of mises_umat!
         DDSDDT[i] = 0.0;
         DRPLDE[i] = 0.0;

         if ( MaterialNumber[0] == 0 || MaterialNumber[0] == 1 )
         {
            DSTRAN[i] = StrainVector[i] - STATEV[i] - STATEV[i+6];
         }

         for ( int j = 0; j < NTENS[0]; j++ )
         {
            DDSDDE[i][j] = 0.0;
         }
      }

      std::cout << " STRAN " << StrainVector << std::endl;
      std::cout << " DSTRAN" << DSTRAN << std::endl;
      for (unsigned int i = 0; i < 6; i++)
         std::cout << " , " << DSTRAN[i];
      std::cout << std::endl;
      std::cout << " STATEV " << STATEV << std::endl;
      for (unsigned int i = 0; i < 12; i++)
         std::cout << "  ,  " << STATEV[i];
      std::cout << std::endl;

      Matrix AlgorithmicTangent = ZeroMatrix(6,6);
      Vector StressVector = ZeroVector(6);

      // NOTE: parameters that are not required by the umat implementations used so far are given as NULL pointers
      // if any new umat is implemented, please check the required parameters and add them accordingly
      // make sure that for backward compatibility the new parameters are initialized as NULL pointers for all
      // other umat materials
      umat_wrapper_( STRESS, STATEV, ( double** ) DDSDDE, NULL, NULL, NULL, NULL, DDSDDT, DRPLDE, NULL, STRAN, DSTRAN,
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

   void Umat::FinalizeSolutionStep( const Properties& props,
         const GeometryType& geom,
         const Vector& ShapeFunctionsValues ,
         const ProcessInfo& CurrentProcessInfo )
   {
      KRATOS_TRY

      //does nothing
      for (unsigned int i = 0; i < 6; i++)
         STRESS_END[i] = STRESS[i];
      for (unsigned int i = 0; i < NSTATV[0]; i++)
         STATEV_END[i] = STATEV[i];

      KRATOS_CATCH("")
   }

   int Umat::Check(const Properties& rMaterialProperties, const GeometryType & rElementGeometry, const ProcessInfo & rCurrentProcessInfo)
   {
      return 0;
   }

   void Umat::GetLawFeatures(Features& rFeatures)
   {
      //Set the type of law
      rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
      rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
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
