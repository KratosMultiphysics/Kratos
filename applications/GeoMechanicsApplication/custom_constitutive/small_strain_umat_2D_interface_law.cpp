// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

// System includes

// External includes

#include "custom_constitutive/small_strain_umat_2D_interface_law.hpp"


namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallStrainUMAT2DInterfaceLaw::SmallStrainUMAT2DInterfaceLaw()
   : SmallStrainUMAT3DLaw()
   {
    KRATOS_TRY;
    //KRATOS_INFO("SmallStrainUMAT2DInterfaceLaw()") << std::endl;

    KRATOS_CATCH("")

   }

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************
SmallStrainUMAT2DInterfaceLaw::
   SmallStrainUMAT2DInterfaceLaw(const SmallStrainUMAT2DInterfaceLaw &rOther)
   : SmallStrainUMAT3DLaw(rOther)
{
   KRATOS_TRY;
   //KRATOS_INFO("SmallStrainUMAT2DInterfaceLaw(const...)") << std::endl;

   KRATOS_CATCH("");
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer SmallStrainUMAT2DInterfaceLaw::Clone() const
{
   KRATOS_TRY;
   //KRATOS_INFO("Clone()") << std::endl;

   return Kratos::make_shared<SmallStrainUMAT2DInterfaceLaw>(*this);

   KRATOS_CATCH("");
}

//********************************ASSIGNMENT******************************************
//************************************************************************************
SmallStrainUMAT2DInterfaceLaw 
  &SmallStrainUMAT2DInterfaceLaw::operator=(SmallStrainUMAT2DInterfaceLaw const &rOther)
{
   KRATOS_TRY;

   SmallStrainUMAT3DLaw::operator=(rOther);

   //KRATOS_INFO("operator=") << std::endl;

   return *this;

   KRATOS_CATCH("");
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

SmallStrainUMAT2DInterfaceLaw::~SmallStrainUMAT2DInterfaceLaw()
{
   //KRATOS_INFO("~SmallStrainUMAT3DLaw()") << std::endl;
}


//************************************************************************************
//************************************************************************************

void SmallStrainUMAT2DInterfaceLaw::UpdateInternalDeltaStrainVector(ConstitutiveLaw::Parameters &rValues)
{
   const Vector& rStrainVector = rValues.GetStrainVector();

   mDeltaStrainVector[INDEX_3D_ZZ] = rStrainVector(INDEX_2D_INTERFACE_ZZ) - mStrainVectorFinalized[INDEX_3D_ZZ];
   mDeltaStrainVector[INDEX_3D_XZ] = rStrainVector(INDEX_2D_INTERFACE_XZ) - mStrainVectorFinalized[INDEX_3D_XZ];
}

void SmallStrainUMAT2DInterfaceLaw::SetExternalStressVector(Vector& rStressVector)
{
   //KRATOS_INFO("mStressVector") << mStressVector << std::endl;

   rStressVector(INDEX_2D_INTERFACE_ZZ) = mStressVector[INDEX_3D_ZZ];
   rStressVector(INDEX_2D_INTERFACE_XZ) = mStressVector[INDEX_3D_XZ];
}


void SmallStrainUMAT2DInterfaceLaw::SetInternalStressVector(const Vector& rStressVector)
{
   // KRATOS_INFO("SetInternalStressVector:rStressVector") << rStressVector << std::endl;
   KRATOS_TRY;
   mStressVectorFinalized[INDEX_3D_ZZ] = rStressVector(INDEX_2D_INTERFACE_ZZ);
   mStressVectorFinalized[INDEX_3D_XZ] = rStressVector(INDEX_2D_INTERFACE_XZ);
   KRATOS_CATCH("")
}

void SmallStrainUMAT2DInterfaceLaw::SetInternalStrainVector(const Vector& rStrainVector)
{
   // KRATOS_INFO("SetInternalStrainVector:rStrainVector") << rStrainVector << std::endl;

   mStrainVectorFinalized[INDEX_3D_ZZ] = rStrainVector(INDEX_2D_INTERFACE_ZZ);
   mStrainVectorFinalized[INDEX_3D_XZ] = rStrainVector(INDEX_2D_INTERFACE_XZ);
}


void SmallStrainUMAT2DInterfaceLaw::CopyConstitutiveMatrix( ConstitutiveLaw::Parameters &rValues,
                                                            Matrix& rConstitutiveMatrix )
{
   if (rValues.GetMaterialProperties()[IS_FORTRAN_UDSM])
   {
      // transfer fortran style matrix to C++ style
      for (unsigned int i = 0; i < VoigtSize; i++) {
         for (unsigned int j = 0; j < VoigtSize; j++) {
            rConstitutiveMatrix(i,j) = mMatrixD[getIndex3D(static_cast<indexStress2DInterface>(j))][getIndex3D(static_cast<indexStress2DInterface>(i))];
         }
      }
   }
   else
   {
      for (unsigned int i = 0; i < VoigtSize; i++) {
         for (unsigned int j = 0; j < VoigtSize; j++) {
            rConstitutiveMatrix(i,j) = mMatrixD[getIndex3D(static_cast<indexStress2DInterface>(i))][getIndex3D(static_cast<indexStress2DInterface>(j))];
         }
      }
   }
}

indexStress3D SmallStrainUMAT2DInterfaceLaw::getIndex3D(indexStress2DInterface index2D)
{
   switch (index2D)
   {
      case INDEX_2D_INTERFACE_ZZ:
        return INDEX_3D_ZZ;
      case INDEX_2D_INTERFACE_XZ:
        return INDEX_3D_XZ;
      default:
        KRATOS_ERROR << "invalid index: " << index2D << std::endl;
   }
}


/***********************************************************************************/
/***********************************************************************************/

void SmallStrainUMAT2DInterfaceLaw::
   CalculateCauchyGreenStrain( ConstitutiveLaw::Parameters& rValues,
                               Vector& rStrainVector )
{
   KRATOS_ERROR << "CalculateCauchyGreenStrain is not implemented in SmallStrainUMAT2DInterfaceLaw" << std::endl;
}

//----------------------------------------------------------------------------------------
Vector& SmallStrainUMAT2DInterfaceLaw::
   GetValue( const Variable<Vector> &rThisVariable,
             Vector &rValue )
{
   // KRATOS_INFO("0-SmallStrainUMAT2DInterfaceLaw::GetValue()") << std::endl;

   if (rThisVariable == STATE_VARIABLES)
   {
      SmallStrainUMAT3DLaw::GetValue(rThisVariable, rValue );
   }
   else if (rThisVariable == CAUCHY_STRESS_VECTOR)
   {
      if (rValue.size() != VoigtSize)
         rValue.resize(VoigtSize);

      rValue[INDEX_2D_INTERFACE_ZZ] = mStressVectorFinalized[INDEX_3D_ZZ];
      rValue[INDEX_2D_INTERFACE_XZ] = mStressVectorFinalized[INDEX_3D_XZ];

   }

   // KRATOS_INFO("1-SmallStrainUMAT2DInterfaceLaw::GetValue()") << std::endl;

    return rValue;
}

//----------------------------------------------------------------------------------------
void SmallStrainUMAT2DInterfaceLaw::SetValue( const Variable<Vector>& rThisVariable,
                                                const Vector& rValue,
                                                const ProcessInfo& rCurrentProcessInfo )
{
   // KRATOS_INFO("02-SmallStrainUMAT2DInterfaceLaw::SetValue()") << std::endl;

   if (rThisVariable == STATE_VARIABLES)
   {
      SmallStrainUMAT3DLaw::SetValue(rThisVariable, rValue, rCurrentProcessInfo );
   }
   else if (rThisVariable == CAUCHY_STRESS_VECTOR)
   {
      if (rValue.size() == VoigtSize) 
      {
         this->SetInternalStressVector(rValue);
      }
   }

   // KRATOS_INFO("12-SmallStrainUMAT2DInterfaceLaw::SetValue()") << std::endl;
}

} // Namespace Kratos
