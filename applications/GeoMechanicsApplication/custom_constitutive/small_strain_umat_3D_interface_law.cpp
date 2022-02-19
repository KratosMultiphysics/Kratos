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

#include "custom_constitutive/small_strain_umat_3D_interface_law.hpp"


namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallStrainUMAT3DInterfaceLaw::SmallStrainUMAT3DInterfaceLaw()
   : SmallStrainUMAT3DLaw()
   {
    KRATOS_TRY;
    //KRATOS_INFO("SmallStrainUMAT3DInterfaceLaw()") << std::endl;

    KRATOS_CATCH("")

   }

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************
SmallStrainUMAT3DInterfaceLaw::
   SmallStrainUMAT3DInterfaceLaw(const SmallStrainUMAT3DInterfaceLaw &rOther)
   : SmallStrainUMAT3DLaw(rOther)
{
   KRATOS_TRY;
   //KRATOS_INFO("SmallStrainUMAT3DInterfaceLaw(const...)") << std::endl;

   KRATOS_CATCH("");
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer SmallStrainUMAT3DInterfaceLaw::Clone() const
{
   KRATOS_TRY;
   //KRATOS_INFO("Clone()") << std::endl;

   return Kratos::make_shared<SmallStrainUMAT3DInterfaceLaw>(*this);

   KRATOS_CATCH("");
}

//********************************ASSIGNMENT******************************************
//************************************************************************************
SmallStrainUMAT3DInterfaceLaw 
  &SmallStrainUMAT3DInterfaceLaw::operator=(SmallStrainUMAT3DInterfaceLaw const &rOther)
{
   KRATOS_TRY;

   SmallStrainUMAT3DLaw::operator=(rOther);

   //KRATOS_INFO("operator=") << std::endl;

   return *this;

   KRATOS_CATCH("");
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

SmallStrainUMAT3DInterfaceLaw::~SmallStrainUMAT3DInterfaceLaw()
{
   //KRATOS_INFO("~SmallStrainUMAT3DLaw()") << std::endl;
}


//************************************************************************************
//************************************************************************************

void SmallStrainUMAT3DInterfaceLaw::UpdateInternalDeltaStrainVector(ConstitutiveLaw::Parameters &rValues)
{
   const Vector& rStrainVector = rValues.GetStrainVector();

   mDeltaStrainVector[INDEX_3D_ZZ] = rStrainVector(INDEX_3D_INTERFACE_ZZ) - mStrainVectorFinalized[INDEX_3D_ZZ];
   mDeltaStrainVector[INDEX_3D_YZ] = rStrainVector(INDEX_3D_INTERFACE_YZ) - mStrainVectorFinalized[INDEX_3D_YZ];
   mDeltaStrainVector[INDEX_3D_XZ] = rStrainVector(INDEX_3D_INTERFACE_XZ) - mStrainVectorFinalized[INDEX_3D_XZ];

}

void SmallStrainUMAT3DInterfaceLaw::SetExternalStressVector(Vector& rStressVector)
{
   //KRATOS_INFO("mStressVector") << mStressVector << std::endl;

   rStressVector(INDEX_3D_INTERFACE_ZZ) = mStressVector[INDEX_3D_ZZ];
   rStressVector(INDEX_3D_INTERFACE_YZ) = mStressVector[INDEX_3D_YZ];
   rStressVector(INDEX_3D_INTERFACE_XZ) = mStressVector[INDEX_3D_XZ];
}


void SmallStrainUMAT3DInterfaceLaw::SetInternalStressVector(const Vector& rStressVector)
{
   // KRATOS_INFO("SetInternalStressVector:rStressVector") << rStressVector << std::endl;
   KRATOS_TRY;
   std::fill(mStressVectorFinalized.begin(), mStressVectorFinalized.end(), 0.0);

   mStressVectorFinalized[INDEX_3D_ZZ] = rStressVector(INDEX_3D_INTERFACE_ZZ);
   mStressVectorFinalized[INDEX_3D_YZ] = rStressVector(INDEX_3D_INTERFACE_YZ);
   mStressVectorFinalized[INDEX_3D_XZ] = rStressVector(INDEX_3D_INTERFACE_XZ);
   KRATOS_CATCH("")
}

void SmallStrainUMAT3DInterfaceLaw::SetInternalStrainVector(const Vector& rStrainVector)
{
   // KRATOS_INFO("SetInternalStrainVector:rStrainVector") << rStrainVector << std::endl;
   std::fill(mStrainVectorFinalized.begin(), mStrainVectorFinalized.end(), 0.0);

   mStrainVectorFinalized[INDEX_3D_ZZ] = rStrainVector(INDEX_3D_INTERFACE_ZZ);
   mStrainVectorFinalized[INDEX_3D_YZ] = rStrainVector(INDEX_3D_INTERFACE_YZ);
   mStrainVectorFinalized[INDEX_3D_XZ] = rStrainVector(INDEX_3D_INTERFACE_XZ);

}


void SmallStrainUMAT3DInterfaceLaw::CopyConstitutiveMatrix( ConstitutiveLaw::Parameters &rValues,
                                                            Matrix& rConstitutiveMatrix )
{
   if (rValues.GetMaterialProperties()[IS_FORTRAN_UDSM])
   {
      // transfer fortran style matrix to C++ style
      for (unsigned int i = 0; i < VoigtSize; i++) {
         for (unsigned int j = 0; j < VoigtSize; j++) {
            rConstitutiveMatrix(i,j) = mMatrixD[getIndex3D(static_cast<indexStress3DInterface>(j))][getIndex3D(static_cast<indexStress3DInterface>(i))];
         }
      }
   }
   else
   {
      for (unsigned int i = 0; i < VoigtSize; i++) {
         for (unsigned int j = 0; j < VoigtSize; j++) {
            rConstitutiveMatrix(i,j) = mMatrixD[getIndex3D(static_cast<indexStress3DInterface>(i))][getIndex3D(static_cast<indexStress3DInterface>(j))];
         }
      }
   }
}

indexStress3D SmallStrainUMAT3DInterfaceLaw::getIndex3D(indexStress3DInterface index3D)
{
   switch (index3D)
   {
      case INDEX_3D_INTERFACE_ZZ:
        return INDEX_3D_ZZ;
      case INDEX_3D_INTERFACE_YZ:
        return INDEX_3D_YZ;
      case INDEX_3D_INTERFACE_XZ:
        return INDEX_3D_XZ;
      default:
        KRATOS_ERROR << "invalid index: " << index3D << std::endl;
   }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainUMAT3DInterfaceLaw::
   CalculateCauchyGreenStrain( ConstitutiveLaw::Parameters& rValues,
                               Vector& rStrainVector )
{
   KRATOS_ERROR << "CalculateCauchyGreenStrain is not implemented in SmallStrainUMAT3DInterfaceLaw" << std::endl;
}

//----------------------------------------------------------------------------------------
Vector& SmallStrainUMAT3DInterfaceLaw::
   GetValue( const Variable<Vector> &rThisVariable,
             Vector &rValue )
{
   // KRATOS_INFO("0-SmallStrainUMAT3DInterfaceLaw::GetValue()") << std::endl;

   if (rThisVariable == STATE_VARIABLES)
   {
      SmallStrainUMAT3DLaw::GetValue(rThisVariable, rValue );
   }
   else if (rThisVariable == CAUCHY_STRESS_VECTOR)
   {
      if (rValue.size() != VoigtSize)
         rValue.resize(VoigtSize);

      rValue[INDEX_3D_INTERFACE_ZZ] = mStressVectorFinalized[INDEX_3D_ZZ];
      rValue[INDEX_3D_INTERFACE_YZ] = mStressVectorFinalized[INDEX_3D_YZ];
      rValue[INDEX_3D_INTERFACE_XZ] = mStressVectorFinalized[INDEX_3D_XZ];

   }

   // KRATOS_INFO("1-SmallStrainUMAT3DInterfaceLaw::GetValue()") << std::endl;

    return rValue;
}

//----------------------------------------------------------------------------------------
void SmallStrainUMAT3DInterfaceLaw::SetValue( const Variable<Vector>& rThisVariable,
                                                const Vector& rValue,
                                                const ProcessInfo& rCurrentProcessInfo )
{
   // KRATOS_INFO("02-SmallStrainUMAT3DInterfaceLaw::SetValue()") << std::endl;

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

   // KRATOS_INFO("12-SmallStrainUMAT3DInterfaceLaw::SetValue()") << std::endl;
}

} // Namespace Kratos
