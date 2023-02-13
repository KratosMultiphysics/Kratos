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

#include "custom_constitutive/small_strain_umat_2D_plane_strain_law.hpp"


namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallStrainUMAT2DPlaneStrainLaw::SmallStrainUMAT2DPlaneStrainLaw()
   : SmallStrainUMAT3DLaw()
   {
    KRATOS_TRY

    KRATOS_CATCH("")

   }

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************
SmallStrainUMAT2DPlaneStrainLaw::
   SmallStrainUMAT2DPlaneStrainLaw(const SmallStrainUMAT2DPlaneStrainLaw &rOther)
   : SmallStrainUMAT3DLaw(rOther)
{
   KRATOS_TRY

   KRATOS_CATCH("")
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer SmallStrainUMAT2DPlaneStrainLaw::Clone() const
{
   KRATOS_TRY

   return Kratos::make_shared<SmallStrainUMAT2DPlaneStrainLaw>(*this);

   KRATOS_CATCH("")
}

//********************************ASSIGNMENT******************************************
//************************************************************************************
SmallStrainUMAT2DPlaneStrainLaw 
  &SmallStrainUMAT2DPlaneStrainLaw::operator=(SmallStrainUMAT2DPlaneStrainLaw const &rOther)
{
   KRATOS_TRY

   SmallStrainUMAT3DLaw::operator=(rOther);

   return *this;

   KRATOS_CATCH("")
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

SmallStrainUMAT2DPlaneStrainLaw::~SmallStrainUMAT2DPlaneStrainLaw() {}


//************************************************************************************
//************************************************************************************

void SmallStrainUMAT2DPlaneStrainLaw::
   UpdateInternalDeltaStrainVector(ConstitutiveLaw::Parameters &rValues)
{
   const Vector& rStrainVector = rValues.GetStrainVector();

   for (unsigned int i = 0; i < VoigtSize; ++i) {
      mDeltaStrainVector[i] = rStrainVector(i) - mStrainVectorFinalized[i];
   }

}

void SmallStrainUMAT2DPlaneStrainLaw::
   SetExternalStressVector(Vector& rStressVector)
{
   KRATOS_TRY

   for (unsigned int i = 0; i < VoigtSize; ++i) {
      rStressVector(i) = mStressVector[i];
   }

   KRATOS_CATCH("")
}

void SmallStrainUMAT2DPlaneStrainLaw::
   SetInternalStressVector(const Vector& rStressVector)
{
   KRATOS_TRY

   for (unsigned int i = 0; i < VoigtSize; ++i) {
      mStressVectorFinalized[i] = rStressVector(i);
   }

   KRATOS_CATCH("")
}

void SmallStrainUMAT2DPlaneStrainLaw::
   SetInternalStrainVector(const Vector& rStrainVector)
{
   KRATOS_TRY

   for (unsigned int i = 0; i < VoigtSize; ++i) {
      mStrainVectorFinalized[i] = rStrainVector(i);
   }

   KRATOS_CATCH("")
}

void SmallStrainUMAT2DPlaneStrainLaw::
   CopyConstitutiveMatrix( ConstitutiveLaw::Parameters &rValues,
                           Matrix& rConstitutiveMatrix )
{
   KRATOS_TRY

   if (rValues.GetMaterialProperties()[IS_FORTRAN_UDSM])
   {
      // transfer fortran style matrix to C++ style
      for (unsigned int i = 0; i < VoigtSize; i++) {
         for (unsigned int j = 0; j < VoigtSize; j++) {
            rConstitutiveMatrix(i,j) = mMatrixD[j][i];
         }
      }
   }
   else
   {
      for (unsigned int i = 0; i < VoigtSize; i++) {
         for (unsigned int j = 0; j < VoigtSize; j++) {
            rConstitutiveMatrix(i,j) = mMatrixD[i][j];
         }
      }
   }

   KRATOS_CATCH("")
}


/***********************************************************************************/
/***********************************************************************************/

void SmallStrainUMAT2DPlaneStrainLaw::
   CalculateCauchyGreenStrain( ConstitutiveLaw::Parameters& rValues,
                               Vector& rStrainVector )
{
   //1.-Compute total deformation gradient
   const Matrix& F = rValues.GetDeformationGradientF();

   // for shells/membranes in case the DeformationGradient is of size 3x3
   BoundedMatrix<double, 2, 2> F2x2;
   for (unsigned int i = 0; i<2; ++i)
      for (unsigned int j = 0; j<2; ++j)
         F2x2(i, j) = F(i, j);

   Matrix E_tensor = prod(trans(F2x2), F2x2);

   for (unsigned int i = 0; i<2; ++i)
      E_tensor(i, i) -= 1.0;

   E_tensor *= 0.5;
   noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector(E_tensor);
}

//----------------------------------------------------------------------------------------
Vector& SmallStrainUMAT2DPlaneStrainLaw::
   GetValue( const Variable<Vector> &rThisVariable,
             Vector &rValue )
{
   // KRATOS_INFO("0-SmallStrainUMAT2DPlaneStrainLaw::GetValue()") << std::endl;

   if (rThisVariable == STATE_VARIABLES) {
      SmallStrainUMAT3DLaw::GetValue(rThisVariable, rValue );
   }
   else if (rThisVariable == CAUCHY_STRESS_VECTOR) {
      if (rValue.size() != VoigtSize)
         rValue.resize(VoigtSize);

      for (unsigned int i = 0; i < VoigtSize; ++i) {
         rValue[i] = mStressVectorFinalized[i];
      }
   }

   // KRATOS_INFO("1-SmallStrainUMAT2DPlaneStrainLaw::GetValue()") << std::endl;

    return rValue;
}

//----------------------------------------------------------------------------------------
void SmallStrainUMAT2DPlaneStrainLaw::
   SetValue( const Variable<Vector>& rThisVariable,
             const Vector& rValue,
             const ProcessInfo& rCurrentProcessInfo )
{
   // KRATOS_INFO("02-SmallStrainUMAT2DPlaneStrainLaw::SetValue()") << std::endl;

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

   // KRATOS_INFO("12-SmallStrainUMAT2DPlaneStrainLaw::SetValue()") << std::endl;
}

} // Namespace Kratos
