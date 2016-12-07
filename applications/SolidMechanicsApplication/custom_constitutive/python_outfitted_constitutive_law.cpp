//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2016 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/python_outfitted_constitutive_law.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

PythonOutfittedConstitutiveLaw::PythonOutfittedConstitutiveLaw()
    : ConstitutiveLaw()
{
	
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

PythonOutfittedConstitutiveLaw::PythonOutfittedConstitutiveLaw(PyObject* pPyConstitutiveLaw)
    : ConstitutiveLaw()
{
  mpPyConstitutiveLaw = boost::make_shared<boost::python::object>(boost::python::object( boost::python::handle<>( boost::python::borrowed(pPyConstitutiveLaw) ) ));
    //boost::python::call_method<void>(mpPyConstitutiveLaw->ptr(), "CallLaw");

}

  
//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

PythonOutfittedConstitutiveLaw::PythonOutfittedConstitutiveLaw(const PythonOutfittedConstitutiveLaw& rOther)
    : ConstitutiveLaw(rOther)
    ,mInverseDeformationGradientF0(rOther.mInverseDeformationGradientF0)
    ,mDeterminantF0(rOther.mDeterminantF0)
    ,mStrainEnergy(rOther.mStrainEnergy)
    ,mpPyConstitutiveLaw(rOther.mpPyConstitutiveLaw)
{
  //deep copy
  mpPyConstitutiveLaw.reset(new boost::python::object(*(rOther.mpPyConstitutiveLaw)));
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer PythonOutfittedConstitutiveLaw::Clone() const
{    
    PythonOutfittedConstitutiveLaw::Pointer p_clone(new PythonOutfittedConstitutiveLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

PythonOutfittedConstitutiveLaw::~PythonOutfittedConstitutiveLaw()
{
}


//*******************************OPERATIONS FROM BASE CLASS***************************
//************************************************************************************

//***********************HAS : DOUBLE - VECTOR - MATRIX*******************************
//************************************************************************************

bool PythonOutfittedConstitutiveLaw::Has( const Variable<double>& rThisVariable )
{
    return false;
}

bool PythonOutfittedConstitutiveLaw::Has( const Variable<Vector>& rThisVariable )
{
    return false;
}

bool PythonOutfittedConstitutiveLaw::Has( const Variable<Matrix>& rThisVariable )
{
    return false;
}


//***********************GET VALUE: DOUBLE - VECTOR - MATRIX**************************
//************************************************************************************

double& PythonOutfittedConstitutiveLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
  if (rThisVariable == STRAIN_ENERGY)
  {
      
    rValue = mStrainEnergy;
  }
    
    return( rValue );
}

Vector& PythonOutfittedConstitutiveLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    return( rValue );
}

Matrix& PythonOutfittedConstitutiveLaw::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    return( rValue );
}


//***********************SET VALUE: DOUBLE - VECTOR - MATRIX**************************
//************************************************************************************


void PythonOutfittedConstitutiveLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                                  const ProcessInfo& rCurrentProcessInfo )
{

  if (rThisVariable == DETERMINANT_F)
    {
      mDeterminantF0 = rValue;
    }
}

void PythonOutfittedConstitutiveLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                                  const ProcessInfo& rCurrentProcessInfo )
{

}

void PythonOutfittedConstitutiveLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                                  const ProcessInfo& rCurrentProcessInfo )
{

}



//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************


void PythonOutfittedConstitutiveLaw::InitializeMaterial( const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues )
{
  mDeterminantF0                = 1;
  mInverseDeformationGradientF0 = identity_matrix<double> (3);
  mStrainEnergy                 = 0;

}

//************************************************************************************
//************************************************************************************


void PythonOutfittedConstitutiveLaw::InitializeSolutionStep( const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry, //this is just to give the array of nodes
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo)
{

}

//************************************************************************************
//************************************************************************************


void PythonOutfittedConstitutiveLaw::FinalizeSolutionStep( const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry, //this is just to give the array of nodes
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo)
{

}



//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************



//*****************************MATERIAL RESPONSES*************************************
//************************************************************************************


void  PythonOutfittedConstitutiveLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
{
//python is not thread safe
#pragma omp critical
  {
    Vector& StressVector                = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix          = rValues.GetConstitutiveMatrix();
    
    //it is important to have the correct size, any resize from python will make a local copy
    //any local copy from python will be lost out of the python local scope
    unsigned int size = this->GetStrainSize();
    if( StressVector.size() != size ) 
      StressVector.resize(size,false);
    if( ConstitutiveMatrix.size1() != size || ConstitutiveMatrix.size2() != size )
      ConstitutiveMatrix.resize(size,size,false);
	
    boost::python::call_method<void>(mpPyConstitutiveLaw->ptr(), "CalculateMaterialResponsePK2",  boost::ref(rValues));

  }
}


//************************************************************************************
//************************************************************************************


void PythonOutfittedConstitutiveLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
{
    this->CalculateMaterialResponsePK2 (rValues);

    Vector& StressVector               = rValues.GetStressVector();
    const Matrix& DeformationGradientF = rValues.GetDeformationGradientF();
    const double& DeterminantF         = rValues.GetDeterminantF();

    TransformStresses(StressVector,DeformationGradientF,DeterminantF,StressMeasure_PK2,StressMeasure_PK1);
}

//************************************************************************************
//************************************************************************************


void PythonOutfittedConstitutiveLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
{
  
//python is not thread safe
#pragma omp critical
  {
    Vector& StressVector                = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix          = rValues.GetConstitutiveMatrix();

    //it is important to have the correct size, any resize from python will make a local copy
    //any local copy from python will be lost out of the python local scope
    unsigned int size = this->GetStrainSize();
    if( StressVector.size() != size ) 
      StressVector.resize(size,false);
    if( ConstitutiveMatrix.size1() != size || ConstitutiveMatrix.size2() != size )
      ConstitutiveMatrix.resize(size,size,false);
    
    boost::python::call_method<void>(mpPyConstitutiveLaw->ptr(), "CalculateMaterialResponseKirchhoff", boost::ref(rValues));
  
  }
}


//************************************************************************************
//************************************************************************************

void PythonOutfittedConstitutiveLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{

    this->CalculateMaterialResponseKirchhoff (rValues);

    const double& DeterminantF          = rValues.GetDeterminantF();
    Vector& StressVector                = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix          = rValues.GetConstitutiveMatrix();

     //Set to cauchy Stress:
    StressVector       /= DeterminantF;
    ConstitutiveMatrix /= DeterminantF;

}


//***********************************UPDATE*******************************************
//************************************************************************************

void PythonOutfittedConstitutiveLaw::FinalizeMaterialResponsePK2 (Parameters& rValues)
{

    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponsePK2 (rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    UpdateInternalVariables( rValues );
}

//************************************************************************************
//************************************************************************************


void PythonOutfittedConstitutiveLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
{

    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponsePK1 (rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    UpdateInternalVariables( rValues );
}

//************************************************************************************
//************************************************************************************


void PythonOutfittedConstitutiveLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
{

    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponseKirchhoff (rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    UpdateInternalVariables( rValues );
}


//************************************************************************************
//************************************************************************************

void PythonOutfittedConstitutiveLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
{
    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponseCauchy (rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    UpdateInternalVariables( rValues );
}


//************************************************************************************
//************************************************************************************

void PythonOutfittedConstitutiveLaw::UpdateInternalVariables(Parameters& rValues)
{
    const Matrix& DeformationGradientF    = rValues.GetDeformationGradientF();
    const double& DeterminantF            = rValues.GetDeterminantF();

    Matrix DeformationGradientF0          = DeformationGradientF;
    DeformationGradientF0 = Transform2DTo3D(DeformationGradientF0);
    MathUtils<double>::InvertMatrix( DeformationGradientF0, this->mInverseDeformationGradientF0, mDeterminantF0);
    mDeterminantF0 = DeterminantF; //special treatment of the determinant 
}


//*************************CONSTITUTIVE LAW PARTICULAR UTILITIES**********************
//************************************************************************************

/**
 * Takes a matrix 2x2 and transforms it to a 3x3 adding a 3rd row and a 3rd column with a 1 in the diagonal
 * if the matrix passed is 3D is does nothing
 * if the matrix passed is bigger or smaller throws an error
 * @param rMatrix : usually the DeformationGradientF
 */
Matrix& PythonOutfittedConstitutiveLaw::Transform2DTo3D (Matrix& rMatrix)
{


    if (rMatrix.size1() == 2 && rMatrix.size2() == 2)
    {

        rMatrix.resize( 3, 3, true);

        rMatrix( 0 , 2 ) = 0.0;
        rMatrix( 1 , 2 ) = 0.0;

        rMatrix( 2 , 0 ) = 0.0;
        rMatrix( 2 , 1 ) = 0.0;

        rMatrix( 2 , 2 ) = 1.0;

    }
    else if(rMatrix.size1() != 3 && rMatrix.size2() != 3)
    {

        KRATOS_THROW_ERROR( std::invalid_argument,"Matrix Dimensions are not correct ", "" )

    }

    return rMatrix;
}

  
//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void PythonOutfittedConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{
  boost::python::call_method<void>(mpPyConstitutiveLaw->ptr(), "GetLawFeatures", boost::ref<Features>(rFeatures));
}


int PythonOutfittedConstitutiveLaw::Check(const Properties& rMaterialProperties,
					  const GeometryType& rElementGeometry,
					  const ProcessInfo& rCurrentProcessInfo)
{
  return boost::python::call_method<int>(mpPyConstitutiveLaw->ptr(),"Check", boost::ref<const Properties>(rMaterialProperties),boost::ref<const GeometryType>(rElementGeometry),boost::ref<const ProcessInfo>(rCurrentProcessInfo));
}

} // Namespace Kratos
