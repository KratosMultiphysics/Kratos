//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:       Massimo Petracca $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:             January 2014 $
//   Revision:            $Revision:                  0.0 $
//
//

#include "shell_cross_section.hpp"
#include "solid_mechanics_application.h"

namespace Kratos
{

ShellCrossSection::ShellCrossSection()
    : mThickness(0.0)
    , mOffset(0.0)
    , mStack()
    , mEditingStack(false)
    , mHasDrillingPenalty(false)
    , mDrillingPenalty(0.0)
    , mOrientation(0.0)
    , mBehavior(Thick)
    , mInitialized(false)
    , mNeedsOOPCondensation(false)
{
}

ShellCrossSection::ShellCrossSection(const ShellCrossSection & other)
{
  PrivateCopy(other);
}

ShellCrossSection::~ShellCrossSection()
{
}

ShellCrossSection & ShellCrossSection::operator = (const ShellCrossSection & other)
{
  PrivateCopy(other);
  return *this;
}

void ShellCrossSection::BeginStack()
{
  if(!mEditingStack)
  {
    mEditingStack = true;
    mThickness = 0.0;
    mStack.clear();
  }
}

void ShellCrossSection::AddPly(double thickness, double orientationAngle, int numPoints, const Properties::Pointer & pProperties)
{
  if((mEditingStack) && (pProperties != NULL) && (thickness > 0.0))
  {
    mStack.push_back( Ply( thickness, 0.0, orientationAngle, numPoints, pProperties ) );
    mThickness += thickness;
  }
}

void ShellCrossSection::EndStack()
{
  if(mEditingStack)
  {
    double currentLocation = mThickness * 0.5;

    for(PlyCollection::iterator it = mStack.begin(); it != mStack.end(); ++it)
    {
      Ply& iPly = *it;
      double iTh = iPly.GetThickness();
      iPly.SetLocation(currentLocation - iTh * 0.5 - mOffset);
      currentLocation -= iTh;
    }

    mEditingStack = false;
  }
}

std::string ShellCrossSection::GetInfo()const
{
  std::stringstream ss;
  ss << std::fixed;

  ss << std::endl;
  ss << "===============================================================" << std::endl;
  ss << "                      SellCrossSection Info:" << std::endl;
  ss << "===============================================================" << std::endl;
  ss << "Total Thickness: " << mThickness << std::endl;
  ss << "Offset from the midplane: " << mOffset << std::endl;
  ss << "Number of Plies: " << mStack.size() << std::endl;
  ss << "===============================================================" << std::endl;
  ss << "=======================       STACK      ======================" << std::endl;
  ss << "===============================================================" << std::endl;
  if(mStack.size() < 1)
  {
    ss << " EMPTY STACK" << std::endl;
    ss << "===============================================================" << std::endl;
  }
  else
  {
    for(PlyCollection::const_iterator it = mStack.begin(); it != mStack.end(); ++it)
    {
      const Ply& iPly = *it;

      ss << " - Thickness :" << iPly.GetThickness() << std::endl;
      ss << " - Location :" << iPly.GetLocation() << std::endl;
      ss << " - Orientation Angle: " << iPly.GetOrientationAngle() << " (degrees)" << std::endl;
      ss << " - Through-The-Thickness Integration Points (" << iPly.GetIntegrationPoints().size() << "):" << std::endl;
      for(unsigned int i = 0; i < iPly.GetIntegrationPoints().size(); i++) {
        const IntegrationPoint& iPoint = iPly.GetIntegrationPoints()[i];
        ss << " - - [" << i << "] "
           << "[ H: " << iPoint.GetWeight() << "; POS: " << iPoint.GetLocation() << "; C-LAW: " << iPoint.GetConstitutiveLaw() << "]"
           << std::endl;
      }
      ss << "===============================================================" << std::endl;
    }
  }
  ss << std::endl;
  return ss.str();
}

ShellCrossSection::Pointer ShellCrossSection::Clone()const
{
  ShellCrossSection::Pointer theClone( new ShellCrossSection(*this) );
  theClone->EndStack();
  return theClone;
}

bool ShellCrossSection::Has(const Variable<double>& rThisVariable)
{
  return false;
}

bool ShellCrossSection::Has(const Variable<Vector>& rThisVariable)
{
  return false;
}

bool ShellCrossSection::Has(const Variable<Matrix>& rThisVariable)
{
  return false;
}

bool ShellCrossSection::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
{
  return false;
}

bool ShellCrossSection::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
{
  return false;
}

double& ShellCrossSection::GetValue(const Variable<double>& rThisVariable, double& rValue)
{
  double meanValue = 0.0;
  double iValue = 0.0;
  double accum = 0.0;
  for(PlyCollection::iterator ply_it = mStack.begin(); ply_it != mStack.end(); ++ply_it)
  {
    Ply& iPly = *ply_it;
    //const Properties& iPlyProps = iPly.GetProperties();
    for(Ply::IntegrationPointCollection::iterator intp_it = iPly.GetIntegrationPoints().begin(); intp_it != iPly.GetIntegrationPoints().end(); ++intp_it)
    {
      IntegrationPoint& iPoint = *intp_it;
      iValue = 0.0;
      if(iPoint.GetConstitutiveLaw()->Has(rThisVariable))
      {
        iValue = iPoint.GetConstitutiveLaw()->GetValue(rThisVariable, iValue);
        meanValue += iValue * iPoint.GetWeight();
        accum += iPoint.GetWeight();
      }
    }
  }
  if(accum != 0.0)
    rValue = meanValue / accum;
  return rValue;
}

Vector& ShellCrossSection::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
  return rValue;
}

Matrix& ShellCrossSection::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
{
  return rValue;
}

array_1d<double, 3 > & ShellCrossSection::GetValue(const Variable<array_1d<double, 3 > >& rVariable,
                                                   array_1d<double, 3 > & rValue)
{
  return rValue;
}

array_1d<double, 6 > & ShellCrossSection::GetValue(const Variable<array_1d<double, 6 > >& rVariable,
                                                   array_1d<double, 6 > & rValue)
{
  return rValue;
}

void ShellCrossSection::SetValue(const Variable<double>& rVariable,
                                 const double& rValue,
                                 const ProcessInfo& rCurrentProcessInfo)
{
}

void ShellCrossSection::SetValue(const Variable<Vector >& rVariable,
                                 const Vector& rValue,
                                 const ProcessInfo& rCurrentProcessInfo)
{
}

void ShellCrossSection::SetValue(const Variable<Matrix >& rVariable,
                                 const Matrix& rValue,
                                 const ProcessInfo& rCurrentProcessInfo)
{
}

void ShellCrossSection::SetValue(const Variable<array_1d<double, 3 > >& rVariable,
                                 const array_1d<double, 3 > & rValue,
                                 const ProcessInfo& rCurrentProcessInfo)
{
}

void ShellCrossSection::SetValue(const Variable<array_1d<double, 6 > >& rVariable,
                                 const array_1d<double, 6 > & rValue,
                                 const ProcessInfo& rCurrentProcessInfo)
{
}

bool ShellCrossSection::ValidateInput(const Properties& rMaterialProperties)
{
  return true;
}

void ShellCrossSection::InitializeCrossSection(const Properties& rMaterialProperties,
                                               const GeometryType& rElementGeometry,
                                               const Vector& rShapeFunctionsValues)
{
  if(!mInitialized)
  {
    if(mEditingStack) EndStack();

    mNeedsOOPCondensation = false;

    for(PlyCollection::iterator ply_it = mStack.begin(); ply_it != mStack.end(); ++ply_it)
    {
      Ply& iPly = *ply_it;
      const Properties& iPlyProps = iPly.GetProperties();
      for(Ply::IntegrationPointCollection::iterator intp_it = iPly.GetIntegrationPoints().begin(); intp_it != iPly.GetIntegrationPoints().end(); ++intp_it)
      {
        IntegrationPoint& iPoint = *intp_it;
        iPoint.GetConstitutiveLaw()->InitializeMaterial(iPlyProps, rElementGeometry, rShapeFunctionsValues);

        if(!mNeedsOOPCondensation)
          if(iPoint.GetConstitutiveLaw()->GetStrainSize() == 6)
            mNeedsOOPCondensation = true;
      }
    }

    if(mNeedsOOPCondensation)
    {
      unsigned int condensed_strain_size = mBehavior == Thick ? 1 : 3;

      if(mOOP_CondensedStrains.size() != condensed_strain_size)
        mOOP_CondensedStrains.resize(condensed_strain_size, false);

      if(mOOP_CondensedStrains_converged.size() != condensed_strain_size)
        mOOP_CondensedStrains_converged.resize(condensed_strain_size, false);

      noalias(mOOP_CondensedStrains) = ZeroVector(condensed_strain_size);
      noalias(mOOP_CondensedStrains_converged) = ZeroVector(condensed_strain_size);
    }

    mInitialized = true;
  }
}

void ShellCrossSection::InitializeSolutionStep(const Properties& rMaterialProperties,
                                               const GeometryType& rElementGeometry,
                                               const Vector& rShapeFunctionsValues,
                                               const ProcessInfo& rCurrentProcessInfo)
{
  for(PlyCollection::iterator ply_it = mStack.begin(); ply_it != mStack.end(); ++ply_it)
  {
    Ply& iPly = *ply_it;
    const Properties& iPlyProps = iPly.GetProperties();
    for(Ply::IntegrationPointCollection::iterator intp_it = iPly.GetIntegrationPoints().begin(); intp_it != iPly.GetIntegrationPoints().end(); ++intp_it)
    {
      IntegrationPoint& iPoint = *intp_it;
      iPoint.GetConstitutiveLaw()->InitializeSolutionStep(iPlyProps, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);
    }
  }

  if(mNeedsOOPCondensation)
    mOOP_CondensedStrains = mOOP_CondensedStrains_converged;
}

void ShellCrossSection::FinalizeSolutionStep(const Properties& rMaterialProperties,
                                             const GeometryType& rElementGeometry,
                                             const Vector& rShapeFunctionsValues,
                                             const ProcessInfo& rCurrentProcessInfo)
{
  for(PlyCollection::iterator ply_it = mStack.begin(); ply_it != mStack.end(); ++ply_it)
  {
    Ply& iPly = *ply_it;
    const Properties& iPlyProps = iPly.GetProperties();
    for(Ply::IntegrationPointCollection::iterator intp_it = iPly.GetIntegrationPoints().begin(); intp_it != iPly.GetIntegrationPoints().end(); ++intp_it)
    {
      IntegrationPoint& iPoint = *intp_it;
      iPoint.GetConstitutiveLaw()->FinalizeSolutionStep(iPlyProps, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);
    }
  }

  if(mNeedsOOPCondensation)
    mOOP_CondensedStrains_converged = mOOP_CondensedStrains;
}

void ShellCrossSection::InitializeNonLinearIteration(const Properties& rMaterialProperties,
                                                     const GeometryType& rElementGeometry,
                                                     const Vector& rShapeFunctionsValues,
                                                     const ProcessInfo& rCurrentProcessInfo)
{
}

void ShellCrossSection::FinalizeNonLinearIteration(const Properties& rMaterialProperties,
                                                   const GeometryType& rElementGeometry,
                                                   const Vector& rShapeFunctionsValues,
                                                   const ProcessInfo& rCurrentProcessInfo)
{
}

void ShellCrossSection::CalculateSectionResponse(Parameters& rValues, const ConstitutiveLaw::StressMeasure& rStressMeasure)
{
  // parameters initialization
  ConstitutiveLaw::Parameters materialValues;
  ElementVariables variables;
  InitializeParameters(rValues, materialValues, variables);

  Flags& Options = rValues.GetOptions();
  bool compute_stress              = Options.Is(ConstitutiveLaw::COMPUTE_STRESS);
  bool compute_constitutive_tensor = Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
  SizeType strain_size = GetStrainSize();
  SizeType condensed_strain_size = GetCondensedStrainSize();

  if(!compute_constitutive_tensor && mNeedsOOPCondensation)
  {
    compute_constitutive_tensor = true;
    Options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    materialValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
  }

  // references
  Vector& generalizedStrainVector = rValues.GetGeneralizedStrainVector();
  Vector& generalizedStressVector = rValues.GetGeneralizedStressVector();
  Matrix& constitutiveMatrix = rValues.GetConstitutiveMatrix();

  Vector& condensedStressVector   = variables.CondensedStressVector;
  Matrix& H  = variables.H;
  Matrix& L  = variables.L;
  Matrix& LT = variables.LT;

  // working matrices to avoid re-allocations when we need to rotate
  // intput and outputs in composite layups
  Matrix R(strain_size, strain_size);
  Matrix DRT(strain_size, strain_size);
  Matrix Rc;
  Matrix HRcT;
  Matrix LRcT;
  Matrix LTRT;
  Matrix Hinv;
  if(mNeedsOOPCondensation)
  {
    Rc.resize(condensed_strain_size, condensed_strain_size, false);
    if(compute_constitutive_tensor)
    {
      HRcT.resize(condensed_strain_size, condensed_strain_size, false);
      LRcT.resize(strain_size, condensed_strain_size, false);
      LTRT.resize(condensed_strain_size, strain_size, false);
    }

    Hinv.resize(condensed_strain_size, condensed_strain_size);
    noalias(Hinv) = ZeroMatrix(condensed_strain_size, condensed_strain_size);
  }
  else if(compute_constitutive_tensor){

    Hinv.resize(condensed_strain_size, condensed_strain_size);
    noalias(Hinv) = ZeroMatrix(condensed_strain_size, condensed_strain_size);

  }


  // compute the generalized strain vector in section coordinate system
  Vector generalizedStrainVector_element;
  if(mOrientation != 0.0)
  {
    // make a copy of the one in element coordinate system. (original)
    generalizedStrainVector_element.resize(strain_size, false);
    noalias( generalizedStrainVector_element ) = generalizedStrainVector;

    // rotate the original one to the section coordinate system
    GetRotationMatrixForGeneralizedStrains( -mOrientation, R );
    noalias( generalizedStrainVector ) = prod( R, generalizedStrainVector_element );
  }

  // ************************************* NOW WE ARE IN THE CROSS SECTION COORDINATE SYSTEM *************************************

  // initialize vector and matrices to store temporary values
  Vector generalizedStrainVector_section;
  Vector generalizedStressVector_section;
  Matrix constitutiveMatrix_section;
  Matrix H_section;
  Matrix L_section;
  Matrix LT_section;
  Vector condensedStressVector_section;
  Vector condensedStrainVector_section;

  // initialize iteration data for condensation of out-of-plane strain components
  int    max_iter = 10;
  double relative_tolerance = 1.0E-8;
  double always_converged_tolerance = 1.0E-10;
  double tolerance = relative_tolerance;
  int    iter = 0;
  double oop_stress_norm = 0.0;
  bool   converged = false;

  // BEGIN LOOP: Newthon iteration to solve for condensed strains
  while(true)
  {
    noalias( generalizedStressVector ) = ZeroVector( strain_size );
    noalias( condensedStressVector ) = ZeroVector( condensed_strain_size );
    noalias( constitutiveMatrix ) = ZeroMatrix( strain_size, strain_size );
    noalias( H ) = ZeroMatrix( condensed_strain_size, condensed_strain_size );
    noalias( L ) = ZeroMatrix( strain_size, condensed_strain_size );
    noalias( LT ) = ZeroMatrix( condensed_strain_size, strain_size );

    // BEGIN LOOP: integrate the response of each ply in this cross section
    for(PlyCollection::iterator ply_it = mStack.begin(); ply_it != mStack.end(); ++ply_it)
    {
      Ply& iPly = *ply_it;
      const Properties& iPlyProps = iPly.GetProperties();
      materialValues.SetMaterialProperties( iPlyProps );
      double iPlyAngle = iPly.GetOrientationAngle();

      if(iPlyAngle == 0.0)
      {
        // BEGIN LOOP: integrate the response of each integration point in this ply
        for(Ply::IntegrationPointCollection::iterator intp_it = iPly.GetIntegrationPoints().begin(); intp_it != iPly.GetIntegrationPoints().end(); ++intp_it)
        {
          IntegrationPoint& iPoint = *intp_it;
          UpdateIntegrationPointParameters(iPoint, materialValues, variables);
          CalculateIntegrationPointResponse(iPoint, materialValues, rValues, variables, rStressMeasure);
        } // END LOOP: integrate the response of each integration point in this ply
      }
      else
      {
        // get the angle in radians of this ply w.r.t the parent section
        double alpha = Globals::Pi / 180.0 * iPlyAngle;

        // make a copy of the generalized strain vector in section coordinate system
        // and then rotate the (working) generalized strain vector in this ply coordinate system
        if(generalizedStrainVector_section.size() != strain_size)
          generalizedStrainVector_section.resize(strain_size, false);
        noalias( generalizedStrainVector_section ) = generalizedStrainVector; // make a copy
        GetRotationMatrixForGeneralizedStrains( -alpha, R );
        noalias( generalizedStrainVector ) = prod( R, generalizedStrainVector_section ); // rotate

        // make a copy of the condensed strain vector in section coordinate system
        // and then rotate the (working) condensed strain vector in this ply coordinate system.
        if(mNeedsOOPCondensation)
        {
          if(condensedStrainVector_section.size() != condensed_strain_size)
            condensedStrainVector_section.resize(condensed_strain_size, false);
          noalias( condensedStrainVector_section ) = mOOP_CondensedStrains; // make a copy
          GetRotationMatrixForCondensedStrains( -alpha, Rc );
          noalias( mOOP_CondensedStrains ) = prod( Rc, condensedStrainVector_section ); // rotate
        }

        // make a copy of the generalized stress vector in section coordinate system (which is being integrated)
        // and then set to zero the (working) generalized stress vector
        if(compute_stress)
        {
          if(generalizedStressVector_section.size() != strain_size)
            generalizedStressVector_section.resize(strain_size, false);
          noalias( generalizedStressVector_section ) = generalizedStressVector; // make a copy
          noalias( generalizedStressVector ) = ZeroVector(strain_size); // set to zero

          if(mNeedsOOPCondensation)
          {
            if(condensedStressVector_section.size() != condensed_strain_size)
              condensedStressVector_section.resize(condensed_strain_size, false);
            noalias( condensedStressVector_section ) = condensedStressVector; // make a copy
            noalias( condensedStressVector ) = ZeroVector(condensed_strain_size); // set to zero
          }
        }

        // make a copy of the section constitutive matrix in section coordinate system (which is being integrated)
        // and then set to zero the (working) section constitutive matrix
        if(compute_constitutive_tensor)
        {
          if(constitutiveMatrix_section.size1() != strain_size || constitutiveMatrix_section.size2() != strain_size)
            constitutiveMatrix_section.resize(strain_size, strain_size, false);
          noalias( constitutiveMatrix_section ) = constitutiveMatrix; // make a copy
          noalias( constitutiveMatrix ) = ZeroMatrix(strain_size, strain_size); // set to zero

          if(mNeedsOOPCondensation)
          {
            if(H_section.size1() != condensed_strain_size || H_section.size2() != condensed_strain_size)
              H_section.resize(condensed_strain_size, condensed_strain_size, false);
            noalias( H_section ) = H; // make a copy
            noalias( H ) = ZeroMatrix(condensed_strain_size, condensed_strain_size); // set to zero

            if(L_section.size1() != strain_size || L_section.size2() != condensed_strain_size)
              L_section.resize(strain_size, condensed_strain_size, false);
            noalias( L_section ) = L; // make a copy
            noalias( L ) = ZeroMatrix(strain_size, condensed_strain_size); // set to zero

            if(LT_section.size1() != condensed_strain_size || L_section.size2() != strain_size)
              LT_section.resize(condensed_strain_size, strain_size, false);
            noalias( LT_section ) = LT; // make a copy
            noalias( LT ) = ZeroMatrix(condensed_strain_size, strain_size); // set to zero
          }
        }

        // BEGIN LOOP: integrate the response of each integration point in this ply
        for(Ply::IntegrationPointCollection::iterator intp_it = iPly.GetIntegrationPoints().begin(); intp_it != iPly.GetIntegrationPoints().end(); ++intp_it)
        {
          IntegrationPoint& iPoint = *intp_it;
          UpdateIntegrationPointParameters(iPoint, materialValues, variables);
          CalculateIntegrationPointResponse(iPoint, materialValues, rValues, variables, rStressMeasure);
        } // END LOOP: integrate the response of each integration point in this ply

        // restore the (working) generalized strain vector with the one in section coordinate system
        noalias( generalizedStrainVector ) = generalizedStrainVector_section;

        // restore the (working) condensed strain vector with the one in section coordinate system
        if(mNeedsOOPCondensation)
          noalias( mOOP_CondensedStrains ) = condensedStrainVector_section;

        // transform the output stress and constitutive matrix from this ply to the parent section
        // coordinate system. then add them to the already integrated quantities.
        if(compute_stress || compute_constitutive_tensor)
        {
          GetRotationMatrixForGeneralizedStresses(alpha, R);
          if(mNeedsOOPCondensation)
            GetRotationMatrixForCondensedStresses(alpha, Rc);

          if(compute_stress)
          {
            noalias( generalizedStressVector_section ) += prod( R, generalizedStressVector );
            noalias( generalizedStressVector ) = generalizedStressVector_section;

            if(mNeedsOOPCondensation)
            {
              noalias( condensedStressVector_section ) += prod( Rc, condensedStressVector );
              noalias( condensedStressVector ) = condensedStressVector_section;
            }
          }
          if(compute_constitutive_tensor)
          {
            noalias( DRT ) = prod( constitutiveMatrix, trans( R ) );
            noalias( constitutiveMatrix_section ) += prod( R, DRT );
            constitutiveMatrix.swap(constitutiveMatrix_section);

            if(mNeedsOOPCondensation)
            {
              noalias( HRcT ) = prod( H, trans( Rc ) );
              noalias( H_section ) += prod( Rc, HRcT );
              noalias( H ) = H_section;

              noalias( LRcT ) = prod( L, trans( Rc ) );
              noalias( L_section ) += prod( R, LRcT );
              noalias( L ) = L_section;

              noalias( LTRT ) = prod( LT, trans( R ) );
              noalias( LT_section ) += prod( Rc, LTRT );
              noalias( LT ) = LT_section;
            }
          }
        }
      }
    } // END LOOP: integrate the response of each ply in this cross section

    // quick return if no static condensation is required
    if(!mNeedsOOPCondensation)
    {
      converged = true;
      break;
    }

    // compute out-of-plane stress norm
    if(mBehavior == Thick)
      oop_stress_norm = std::abs( condensedStressVector(0) );
    else
      oop_stress_norm = norm_2( condensedStressVector );

    // initialize tolerance
    if(iter == 0)
    {
      tolerance = oop_stress_norm * relative_tolerance;
      if(tolerance < always_converged_tolerance)
        tolerance = always_converged_tolerance;
    }

    // compute H^-1
    if(mBehavior == Thick)
    {
      Hinv(0, 0) = 1.0 / H(0, 0);
    }
    else
    {
      double dummy_det;
      MathUtils<double>::InvertMatrix3(H, Hinv, dummy_det);
    }

    // check convergence
    if(oop_stress_norm <= tolerance)
    {
      converged = true;
      break;
    }

    // update out-of-plane strains
    noalias( mOOP_CondensedStrains ) -= prod( Hinv, condensedStressVector );

    iter++;

    if(iter > max_iter) break;

  } // END LOOP: Newthon iteration

  if(!converged || compute_constitutive_tensor)
  {
    Matrix LHinv( prod( L, Hinv ) );
    if(!converged && compute_stress)
    {
      noalias( generalizedStressVector ) += prod( LHinv, condensedStressVector );
    }

    if(compute_constitutive_tensor)
    {
      noalias( constitutiveMatrix ) -= prod( LHinv, LT );
    }
  }

  // *********************************** NOW WE MOVE TO THE PARENT ELEMENT COORDINATE SYSTEM ************************************

  // transform the outputs back to the element coordinate system (if necessary)
  if(mOrientation != 0.0)
  {
    if(compute_stress || compute_constitutive_tensor)
    {
      GetRotationMatrixForGeneralizedStresses( mOrientation, R );
      if(compute_stress)
      {
        generalizedStressVector = prod( R, generalizedStressVector );
      }
      if(compute_constitutive_tensor)
      {
        noalias( DRT ) = prod( constitutiveMatrix, trans( R ) );
        noalias( constitutiveMatrix ) = prod( R, DRT );
      }
    }
  }

  // restore the original strain vector in element coordinate system
  if(mOrientation != 0.0)
    noalias( generalizedStrainVector ) = generalizedStrainVector_element;

  // compute the drilling stiffness parameter
  if(!mHasDrillingPenalty && compute_constitutive_tensor)
  {
    mDrillingPenalty = constitutiveMatrix(2, 2);
    mHasDrillingPenalty = true;
  }
}

void ShellCrossSection::FinalizeSectionResponse(Parameters& rValues, const ConstitutiveLaw::StressMeasure& rStressMeasure)
{
  ConstitutiveLaw::Parameters materialValues;
  ElementVariables variables;
  InitializeParameters(rValues, materialValues, variables);

  for(PlyCollection::iterator ply_it = mStack.begin(); ply_it != mStack.end(); ++ply_it)
  {
    Ply& iPly = *ply_it;
    //const Properties& iPlyProps = iPly.GetProperties();
    for(Ply::IntegrationPointCollection::iterator intp_it = iPly.GetIntegrationPoints().begin(); intp_it != iPly.GetIntegrationPoints().end(); ++intp_it)
    {
      IntegrationPoint& iPoint = *intp_it;
      UpdateIntegrationPointParameters(iPoint, materialValues, variables);
      iPoint.GetConstitutiveLaw()->FinalizeMaterialResponse(materialValues, rStressMeasure);
    }
  }
}

void ShellCrossSection::ResetCrossSection(const Properties& rMaterialProperties,
                                          const GeometryType& rElementGeometry,
                                          const Vector& rShapeFunctionsValues)
{
  mInitialized = false;

  for(PlyCollection::iterator ply_it = mStack.begin(); ply_it != mStack.end(); ++ply_it)
  {
    Ply& iPly = *ply_it;
    const Properties& iPlyProps = iPly.GetProperties();
    for(Ply::IntegrationPointCollection::iterator intp_it = iPly.GetIntegrationPoints().begin(); intp_it != iPly.GetIntegrationPoints().end(); ++intp_it)
    {
      IntegrationPoint& iPoint = *intp_it;
      iPoint.GetConstitutiveLaw()->ResetMaterial(iPlyProps, rElementGeometry, rShapeFunctionsValues);
    }
  }

  if(mNeedsOOPCondensation)
  {
    int condensed_strain_size = mBehavior == Thick ? 1 : 3;

    noalias(mOOP_CondensedStrains) = ZeroVector(condensed_strain_size);
    noalias(mOOP_CondensedStrains_converged) = ZeroVector(condensed_strain_size);
  }
}

int ShellCrossSection::Check(const Properties& rMaterialProperties,
                             const GeometryType& rElementGeometry,
                             const ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

      if(this->mEditingStack)
        KRATOS_THROW_ERROR(std::logic_error, "The Ply Stack of a ShellCrossSection is in Editing mode", "")

            if(this->mStack.size() < 1)
              KRATOS_THROW_ERROR(std::logic_error, "The Ply Stack of a ShellCrossSection cannot be empty", "")

                  if(this->mThickness <= 0.0)
                    KRATOS_THROW_ERROR(std::logic_error, "The Thickness of a ShellCrossSection should be a positive real number", this->mThickness)

                        for(PlyCollection::iterator ply_it = mStack.begin(); ply_it != mStack.end(); ++ply_it)
                        {
                          Ply& iPly = *ply_it;

                          Ply::IntegrationPointCollection::size_type numip( iPly.NumberOfIntegrationPoints() );
                          if(numip < 1)
                            KRATOS_THROW_ERROR(std::logic_error, "The number of integration points in a Ply is not set properly", numip);

                          if(iPly.GetPropertiesPointer() == NULL)
                            KRATOS_THROW_ERROR(std::logic_error, "The Properties of a Ply cannot be NULL", "");

                          const Properties & iPlyProps = iPly.GetProperties();

                          if(!iPlyProps.Has(DENSITY))
                            KRATOS_THROW_ERROR(std::logic_error, "DENSITY not provided for a Ply object","");

                          for(Ply::IntegrationPointCollection::iterator intp_it = iPly.GetIntegrationPoints().begin(); intp_it != iPly.GetIntegrationPoints().end(); ++intp_it)
                          {
                            IntegrationPoint& iPoint = *intp_it;

                            if(iPoint.GetWeight() <= 0.0)
                              KRATOS_THROW_ERROR(std::logic_error, "The Weight of a ShellCrossSection.IntegrationPoint should be a positive real number", iPoint.GetWeight());

                            const ConstitutiveLaw::Pointer& iPointLaw = iPoint.GetConstitutiveLaw();

                            if(iPointLaw == NULL)
                              KRATOS_THROW_ERROR(std::logic_error, "The Constitutive law of a ShellCrossSection.IntegrationPoint is NULL", iPointLaw);

                            ConstitutiveLaw::Features iPointLawFeatures;
                            iPointLaw->GetLawFeatures(iPointLawFeatures);

                            int correct_strain_size = iPointLawFeatures.mStrainSize;
                            if(correct_strain_size != 3 && correct_strain_size != 6)
                              KRATOS_THROW_ERROR(std::logic_error,
                                                 "The Constitutive law of a ShellCrossSection.IntegrationPoint needs a ConstitutiveLaw with 3 or 6 components, instead of ", correct_strain_size);

                            //bool correct_strain_measure = false;
                            //for(unsigned int i=0; i<iPointLawFeatures.mStrainMeasures.size(); i++)
                            //{
                            //	if(iPointLawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Infinitesimal){
                            //		correct_strain_measure = true;
                            //		break;
                            //	}
                            //}
                            //
                            //if( correct_strain_measure == false )
                            //	KRATOS_THROW_ERROR( std::logic_error,
                            //		"The Constitutive law of a ShellCrossSection.IntegrationPoint is incompatible with the strain measure required by this cross section ",
                            //		"Required strain measure: StrainMeasure_Infinitesimal" );

                            iPointLaw->Check(iPlyProps, rElementGeometry, rCurrentProcessInfo);
                          }
                        }

  return 0;

  KRATOS_CATCH("")
      }

void ShellCrossSection::InitializeParameters(Parameters& rValues, ConstitutiveLaw::Parameters& rMaterialValues, ElementVariables& rVariables)
{
  // share common data between section and materials

  rMaterialValues.SetOptions(rValues.GetOptions());
  rMaterialValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);

  rMaterialValues.SetShapeFunctionsValues( rValues.GetShapeFunctionsValues() );
  rMaterialValues.SetShapeFunctionsDerivatives( rValues.GetShapeFunctionsDerivatives() );

  rMaterialValues.SetProcessInfo( rValues.GetProcessInfo() );
  rMaterialValues.SetElementGeometry( rValues.GetElementGeometry() );

  // initialize the general variables container

  rVariables.DeterminantF = 1.0;

  rVariables.DeformationGradientF_2D = IdentityMatrix(2,2);
  rVariables.StrainVector_2D.resize(3);
  rVariables.StressVector_2D.resize(3);
  rVariables.ConstitutiveMatrix_2D.resize(3,3);
  noalias( rVariables.StrainVector_2D ) = ZeroVector(3);
  noalias( rVariables.StressVector_2D ) = ZeroVector(3);
  noalias( rVariables.ConstitutiveMatrix_2D ) = ZeroMatrix(3,3);

  if(mNeedsOOPCondensation) // avoid useless allocations
  {
    rVariables.DeformationGradientF_3D = IdentityMatrix(3,3);
    rVariables.DeformationGradientF0_3D = IdentityMatrix(3,3);
    rVariables.StrainVector_3D.resize(6);
    rVariables.StressVector_3D.resize(6);
    rVariables.ConstitutiveMatrix_3D.resize(6,6);
  }

  // by default set the 2D data for materials

  rMaterialValues.SetStrainVector(rVariables.StrainVector_2D);
  rMaterialValues.SetStressVector(rVariables.StressVector_2D);
  rMaterialValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix_2D);
  rMaterialValues.SetDeterminantF(rVariables.DeterminantF);
  rMaterialValues.SetDeformationGradientF(rVariables.DeformationGradientF_2D);

  // initialize to zero the generalized vectors / matrices because
  // they will be integrated
  int strain_size = GetStrainSize();
  int condensed_strain_size = GetCondensedStrainSize();
  noalias( rValues.GetGeneralizedStressVector() ) = ZeroVector(strain_size);
  noalias( rValues.GetConstitutiveMatrix() ) = ZeroMatrix(strain_size, strain_size);
  rVariables.CondensedStressVector = ZeroVector(condensed_strain_size);
  rVariables.H = ZeroMatrix(condensed_strain_size, condensed_strain_size);
  rVariables.L = ZeroMatrix(strain_size, condensed_strain_size);
  rVariables.LT = ZeroMatrix(condensed_strain_size, strain_size);
}

void ShellCrossSection::UpdateIntegrationPointParameters(IntegrationPoint& rPoint, ConstitutiveLaw::Parameters& rMaterialValues, ElementVariables& rVariables)
{
  if(rPoint.GetConstitutiveLaw()->GetStrainSize() == 3)
  {
    // use 2D matrices and vectors
    rMaterialValues.SetStrainVector(rVariables.StrainVector_2D);
    rMaterialValues.SetStressVector(rVariables.StressVector_2D);
    rMaterialValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix_2D);
    rMaterialValues.SetDeterminantF(rVariables.DeterminantF);
    rMaterialValues.SetDeformationGradientF(rVariables.DeformationGradientF_2D);

    if(mBehavior == Thick)
    {
      // get elastic data for the trasverse shear part (treated elastically)
      const Properties& props = rMaterialValues.GetMaterialProperties();
      if(props.Has(SHEAR_MODULUS_YZ) && props.Has(SHEAR_MODULUS_XZ))
      {
        rVariables.GYZ = props[SHEAR_MODULUS_YZ];
        rVariables.GXZ = props[SHEAR_MODULUS_XZ];
      }
      else if(props.Has(YOUNG_MODULUS) && props.Has(POISSON_RATIO))
      {
        double giso = props[YOUNG_MODULUS] / (2.0 * (1.0 + props[POISSON_RATIO]));
        rVariables.GYZ = giso;
        rVariables.GXZ = giso;
      }
      else
      {
        // This should NEVER happen!
        rVariables.GYZ = 0.0;
        rVariables.GXZ = 0.0;
      }
    }
  }
  else // 6
  {
    // use 3D matrices and vectors
    rMaterialValues.SetStrainVector(rVariables.StrainVector_3D);
    rMaterialValues.SetStressVector(rVariables.StressVector_3D);
    rMaterialValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix_3D);
    rMaterialValues.SetDeterminantF(rVariables.DeterminantF);
    rMaterialValues.SetDeformationGradientF(rVariables.DeformationGradientF_3D);
  }
}

void ShellCrossSection::CalculateIntegrationPointResponse(IntegrationPoint& rPoint,
                                                          ConstitutiveLaw::Parameters& rMaterialValues,
                                                          Parameters& rValues,
                                                          ElementVariables& rVariables,
                                                          const ConstitutiveLaw::StressMeasure& rStressMeasure)
{
  // get some data/references...

  Flags& Options = rValues.GetOptions();
  bool compute_stress              = Options.Is(ConstitutiveLaw::COMPUTE_STRESS);
  bool compute_constitutive_tensor = Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

  double h = rPoint.GetWeight();
  double z = rPoint.GetLocation();

  const Vector& generalizedStrainVector = rValues.GetGeneralizedStrainVector();
  Vector& generalizedStressVector       = rValues.GetGeneralizedStressVector();
  Matrix& sectionConstitutiveMatrix     = rValues.GetConstitutiveMatrix();

  Vector& materialStrainVector       = rMaterialValues.GetStrainVector();
  Vector& materialStressVector       = rMaterialValues.GetStressVector();
  Matrix& materialConstitutiveMatrix = rMaterialValues.GetConstitutiveMatrix();

  Vector& condensedStressVector = rVariables.CondensedStressVector;
  Matrix& H                     = rVariables.H;
  Matrix& L                     = rVariables.L;
  Matrix& LT                    = rVariables.LT;

  int material_strain_size = rPoint.GetConstitutiveLaw()->GetStrainSize();

  // shear correction factors
  // standard: but unsymmetric for anisotropic 3d materials
  //double ce = 1.0;
  //double cs = 5.0 / 6.0;
  // modified: symmetric for anisotropic 3d materials
  double ce = std::sqrt(5.0 / 6.0);
  double cs = ce;
  // modified quadratic strains: quadratic shear strains, no modification for stresses, unsymmetric for anisotropic 3d materials
  //double ce = 1.0 - 4.0 / (mThickness * mThickness) * z * z;
  //double cs = 1.0;

  // calculate the material strain vector.

  if(material_strain_size == 3) // plane-stress case
  {
    materialStrainVector(0) = generalizedStrainVector(0) + z * generalizedStrainVector(3); //  e.xx
    materialStrainVector(1) = generalizedStrainVector(1) + z * generalizedStrainVector(4); //  e.yy
    materialStrainVector(2) = generalizedStrainVector(2) + z * generalizedStrainVector(5); //  e.xy
  }
  else // full 3D case
  {
    materialStrainVector(0) = generalizedStrainVector(0) + z * generalizedStrainVector(3);	//  e.xx
    materialStrainVector(1) = generalizedStrainVector(1) + z * generalizedStrainVector(4);	//  e.yy
    materialStrainVector(2) = mOOP_CondensedStrains(0);										//  e.zz (condensed)
    materialStrainVector(3) = generalizedStrainVector(2) + z * generalizedStrainVector(5);	//  e.xy
    if(mBehavior == Thick)
    {
      materialStrainVector(4) = ce * generalizedStrainVector(6);							// 2e.yz
      materialStrainVector(5) = ce * generalizedStrainVector(7);							// 2e.xz
    }
    else // Thin
    {
      materialStrainVector(4) = ce * mOOP_CondensedStrains(1);							// 2e.yz (condensed)
      materialStrainVector(5) = ce * mOOP_CondensedStrains(2);							// 2e.xz (condensed)
    }
  }

  // calculate the deformation gradient
  // here we consider F* = R'*F = U -> approx -> I + eps

  if(material_strain_size == 3)
  {
    Matrix& F = rVariables.DeformationGradientF_2D;
    F(0, 0) = materialStrainVector(0) + 1.0;
    F(1, 1) = materialStrainVector(1) + 1.0;
    F(0, 1) = F(1, 0) = materialStrainVector(2) * 0.5;
    rVariables.DeterminantF = MathUtils<double>::Det2(F);
  }
  else // 6
  {
    Matrix& F = rVariables.DeformationGradientF_3D;
    F(0, 0) = materialStrainVector(0) + 1.0; // xx
    F(1, 1) = materialStrainVector(1) + 1.0; // yy
    F(2, 2) = materialStrainVector(2) + 1.0; // zz
    F(0, 1) = F(1, 0) = materialStrainVector(3) * 0.5; // xy - yx
    F(0, 2) = F(2, 0) = materialStrainVector(5) * 0.5; // xz - zx
    F(1, 2) = F(2, 1) = materialStrainVector(4) * 0.5; // yz - zy
    rVariables.DeterminantF = MathUtils<double>::Det3(F);
  }
  rVariables.DeterminantF0 = 1.0;

  // calculate the material response

  rPoint.GetConstitutiveLaw()->CalculateMaterialResponse(rMaterialValues, rStressMeasure);

  // compute stress resultants and stress couples

  if(compute_stress)
  {
    if(material_strain_size == 3) // plane-stress case
    {
      generalizedStressVector(0) += h * materialStressVector(0);			// N.xx
      generalizedStressVector(1) += h * materialStressVector(1);			// N.yy
      generalizedStressVector(2) += h * materialStressVector(2);			// N.xy
      generalizedStressVector(3) += h * z * materialStressVector(0);		// M.xx
      generalizedStressVector(4) += h * z * materialStressVector(1);		// M.yy
      generalizedStressVector(5) += h * z * materialStressVector(2);		// M.xy
      if(mBehavior == Thick)
      {
        // here the transverse shear is treated elastically
        generalizedStressVector(6) += cs * h * rVariables.GYZ * ce * generalizedStrainVector(6);		// V.yz
        generalizedStressVector(7) += cs * h * rVariables.GXZ * ce * generalizedStrainVector(7);		// V.xz
      }
    }
    else // full 3D case
    {
      generalizedStressVector(0) += h * materialStressVector(0);			// N.xx
      generalizedStressVector(1) += h * materialStressVector(1);			// N.yy
      condensedStressVector(0)   += h * materialStressVector(2);			// N.zz (condensed - should be 0 after integration)
      generalizedStressVector(2) += h * materialStressVector(3);			// N.xy
      generalizedStressVector(3) += h * z * materialStressVector(0);		// M.xx
      generalizedStressVector(4) += h * z * materialStressVector(1);		// M.yy
      generalizedStressVector(5) += h * z * materialStressVector(3);		// M.xy
      if(mBehavior == Thick)
      {
        generalizedStressVector(6) += cs * h * materialStressVector(4);		// V.yz
        generalizedStressVector(7) += cs * h * materialStressVector(5);		// V.xz
      }
      else // Thin
      {
        condensedStressVector(1) += cs * h * materialStressVector(4);		// V.yz (condensed - should be 0 after integration)
        condensedStressVector(2) += cs * h * materialStressVector(5);		// V.xz (condensed - should be 0 after integration)
      }
    }
  }

  // compute the section tangent matrix

  if(compute_constitutive_tensor)
  {
    Matrix & C = materialConstitutiveMatrix;
    Matrix & D = sectionConstitutiveMatrix;

    if(material_strain_size == 3) // plane-stress case
    {
      // membrane part
      D(0,0) += h*C(0,0);			D(0,1) += h*C(0,1);			D(0,2) += h*C(0,2);
      D(1,0) += h*C(1,0);			D(1,1) += h*C(1,1);			D(1,2) += h*C(1,2);
      D(2,0) += h*C(2,0);			D(2,1) += h*C(2,1);			D(2,2) += h*C(2,2);

      // bending part
      D(3,3) += h*z*z*C(0,0);		D(3,4) += h*z*z*C(0,1);		D(3,5) += h*z*z*C(0,2);
      D(4,3) += h*z*z*C(1,0);		D(4,4) += h*z*z*C(1,1);		D(4,5) += h*z*z*C(1,2);
      D(5,3) += h*z*z*C(2,0);		D(5,4) += h*z*z*C(2,1);		D(5,5) += h*z*z*C(2,2);

      // membrane-bending part
      D(0,3) += h*z*C(0,0);		D(0,4) += h*z*C(0,1);		D(0,5) += h*z*C(0,2);
      D(1,3) += h*z*C(1,0);		D(1,4) += h*z*C(1,1);		D(1,5) += h*z*C(1,2);
      D(2,3) += h*z*C(2,0);		D(2,4) += h*z*C(2,1);		D(2,5) += h*z*C(2,2);

      // bending-membrane part
      D(3,0) += h*z*C(0,0);		D(3,1) += h*z*C(0,1);		D(3,2) += h*z*C(0,2);
      D(4,0) += h*z*C(1,0);		D(4,1) += h*z*C(1,1);		D(4,2) += h*z*C(1,2);
      D(5,0) += h*z*C(2,0);		D(5,1) += h*z*C(2,1);		D(5,2) += h*z*C(2,2);

      if(mBehavior == Thick)
      {
        // here the transverse shear is treated elastically
        D(6,6) += h * cs * ce * rVariables.GYZ;
        D(7,7) += h * cs * ce * rVariables.GXZ;
      }
    }
    else // full 3D case
    {
      // membrane part
      D(0,0) += h*C(0,0);			D(0,1) += h*C(0,1);			D(0,2) += h*C(0,3);
      D(1,0) += h*C(1,0);			D(1,1) += h*C(1,1);			D(1,2) += h*C(1,3);
      D(2,0) += h*C(3,0);			D(2,1) += h*C(3,1);			D(2,2) += h*C(3,3);

      // bending part
      D(3,3) += h*z*z*C(0,0);		D(3,4) += h*z*z*C(0,1);		D(3,5) += h*z*z*C(0,3);
      D(4,3) += h*z*z*C(1,0);		D(4,4) += h*z*z*C(1,1);		D(4,5) += h*z*z*C(1,3);
      D(5,3) += h*z*z*C(3,0);		D(5,4) += h*z*z*C(3,1);		D(5,5) += h*z*z*C(3,3);

      // membrane-bending part
      D(0,3) += h*z*C(0,0);		D(0,4) += h*z*C(0,1);		D(0,5) += h*z*C(0,3);
      D(1,3) += h*z*C(1,0);		D(1,4) += h*z*C(1,1);		D(1,5) += h*z*C(1,3);
      D(2,3) += h*z*C(3,0);		D(2,4) += h*z*C(3,1);		D(2,5) += h*z*C(3,3);

      // bending-membrane part
      D(3,0) += h*z*C(0,0);		D(3,1) += h*z*C(0,1);		D(3,2) += h*z*C(0,3);
      D(4,0) += h*z*C(1,0);		D(4,1) += h*z*C(1,1);		D(4,2) += h*z*C(1,3);
      D(5,0) += h*z*C(3,0);		D(5,1) += h*z*C(3,1);		D(5,2) += h*z*C(3,3);

      if(mBehavior == Thick)
      {
        // membrane-shear part
        D(0,6) += ce*h*C(0, 4);   D(0,7) += ce*h*C(0, 5);
        D(1,6) += ce*h*C(1, 4);   D(1,7) += ce*h*C(1, 5);
        D(2,6) += ce*h*C(3, 4);   D(2,7) += ce*h*C(3, 5);

        // bending-shear part
        D(3,6) += ce*h*z*C(0, 4);   D(3,7) += ce*h*z*C(0, 5);
        D(4,6) += ce*h*z*C(1, 4);   D(4,7) += ce*h*z*C(1, 5);
        D(5,6) += ce*h*z*C(3, 4);   D(5,7) += ce*h*z*C(3, 5);

        // shear-membrane part
        D(6,0) += cs*h*C(4, 0);   D(6,1) += cs*h*C(4, 1);   D(6,2) += cs*h*C(4, 3);
        D(7,0) += cs*h*C(5, 0);   D(7,1) += cs*h*C(5, 1);   D(7,2) += cs*h*C(5, 3);

        // shear-bending part
        D(6,3) += cs*h*z*C(4, 0);   D(6,4) += cs*h*z*C(4, 1);   D(6,5) += cs*h*z*C(4, 3);
        D(7,3) += cs*h*z*C(5, 0);   D(7,4) += cs*h*z*C(5, 1);   D(7,5) += cs*h*z*C(5, 3);

        // shear part
        D(6,6) += cs*ce*h*C(4, 4);   D(6,7) += cs*ce*h*C(4, 5);
        D(7,6) += cs*ce*h*C(5, 4);   D(7,7) += cs*ce*h*C(5, 5);

        // matrices for static condensation

        H(0,0) += h*C(2,2);

        LT(0,0) += h*C(2, 0);   LT(0,1) += h*C(2, 1);   LT(0,2) += h*C(2, 3);   LT(0,3) += h*z*C(2, 0);   LT(0,4) += h*z*C(2, 1);   LT(0,5) += h*z*C(2, 3);   LT(0,6) += ce*h*C(2, 4);   LT(0,7) += ce*h*C(2, 5);

        L(0,0) +=    h*C(0, 2);
        L(1,0) +=    h*C(1, 2);
        L(2,0) +=    h*C(3, 2);
        L(3,0) +=  h*z*C(0, 2);
        L(4,0) +=  h*z*C(1, 2);
        L(5,0) +=  h*z*C(3, 2);
        L(6,0) += cs*h*C(4, 2);
        L(7,0) += cs*h*C(5, 2);
      }
      else
      {
        // matrices for static condensation

        H(0,0) +=    h*C(2, 2);   H(0,1) +=    ce*h*C(2, 4);   H(0,2) +=    ce*h*C(2, 5);
        H(1,0) += cs*h*C(4, 2);   H(1,1) += ce*cs*h*C(4, 4);   H(1,2) += ce*cs*h*C(4, 5);
        H(2,0) += cs*h*C(5, 2);   H(2,1) += ce*cs*h*C(5, 4);   H(2,2) += ce*cs*h*C(5, 5);

        LT(0,0) +=    h*C(2, 0);   LT(0,1) +=    h*C(2, 1);   LT(0,2) +=    h*C(2, 3);   LT(0,3) +=    h*z*C(2, 0);   LT(0,4) +=    h*z*C(2, 1);   LT(0,5) +=    h*z*C(2, 3);
        LT(1,0) += cs*h*C(4, 0);   LT(1,1) += cs*h*C(4, 1);   LT(1,2) += cs*h*C(4, 3);   LT(1,3) += cs*h*z*C(4, 0);   LT(1,4) += cs*h*z*C(4, 1);   LT(1,5) += cs*h*z*C(4, 3);
        LT(2,0) += cs*h*C(5, 0);   LT(2,1) += cs*h*C(5, 1);   LT(2,2) += cs*h*C(5, 3);   LT(2,3) += cs*h*z*C(5, 0);   LT(2,4) += cs*h*z*C(5, 1);   LT(2,5) += cs*h*z*C(5, 3);

        L(0,0) +=   h*C(0, 2);  L(0,1) +=   ce*h*C(0, 4);  L(0,2) +=   ce*h*C(0, 5);
        L(1,0) +=   h*C(1, 2);  L(1,1) +=   ce*h*C(1, 4);  L(1,2) +=   ce*h*C(1, 5);
        L(2,0) +=   h*C(3, 2);  L(2,1) +=   ce*h*C(3, 4);  L(2,2) +=   ce*h*C(3, 5);
        L(3,0) += h*z*C(0, 2);  L(3,1) += ce*h*z*C(0, 4);  L(3,2) += ce*h*z*C(0, 5);
        L(4,0) += h*z*C(1, 2);  L(4,1) += ce*h*z*C(1, 4);  L(4,2) += ce*h*z*C(1, 5);
        L(5,0) += h*z*C(3, 2);  L(5,1) += ce*h*z*C(3, 4);  L(5,2) += ce*h*z*C(3, 5);
      }
    }
  }
}

void ShellCrossSection::PrivateCopy(const ShellCrossSection & other)
{
  if(this != &other)
  {
    mThickness = other.mThickness;
    mOffset = other.mOffset;
    mStack = other.mStack;
    mEditingStack = other.mEditingStack;
    mHasDrillingPenalty = other.mHasDrillingPenalty;
    mDrillingPenalty = other.mDrillingPenalty;
    mOrientation = other.mOrientation;
    mBehavior = other.mBehavior;
    mInitialized = other.mInitialized;
    mNeedsOOPCondensation = other.mNeedsOOPCondensation;
    mOOP_CondensedStrains = other.mOOP_CondensedStrains;
    mOOP_CondensedStrains_converged = other.mOOP_CondensedStrains_converged;
  }
}

}
