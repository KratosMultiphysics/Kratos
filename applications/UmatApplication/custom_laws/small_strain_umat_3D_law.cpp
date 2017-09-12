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
#include "custom_laws/small_strain_umat_3D_law.hpp"


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

  SmallStrainUmat3DLaw::SmallStrainUmat3DLaw() : Umat3DLaw()
  {
  }
  
  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  SmallStrainUmat3DLaw::SmallStrainUmat3DLaw(const SmallStrainUmat3DLaw& rOther) : Umat3DLaw(rOther)
  {
  }

  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  SmallStrainUmat3DLaw& SmallStrainUmat3DLaw::operator=(const SmallStrainUmat3DLaw& rOther)
  {
    return *this;
  } 

  
  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer SmallStrainUmat3DLaw::Clone() const
  {
    return ( SmallStrainUmat3DLaw::Pointer(new SmallStrainUmat3DLaw(*this)) );
  }



  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  SmallStrainUmat3DLaw::~SmallStrainUmat3DLaw()
  {
  }

  //*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
  //************************************************************************************

  void SmallStrainUmat3DLaw::GetLawFeatures(Features& rFeatures)
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
