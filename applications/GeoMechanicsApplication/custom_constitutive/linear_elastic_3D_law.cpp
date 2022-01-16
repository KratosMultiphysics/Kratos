// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Hoang-Giang Bui
//

// System includes
#include <iostream>

// External includes
#include <cmath>

// Project includes

#include "custom_constitutive/linear_elastic_3D_law.h"
#include "custom_utilities/math_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{


/**
 * TO BE TESTED!!!
 */
template<int TType>
LinearElastic<TType>::LinearElastic() : ConstitutiveLaw()
{
}

template<int TType>
bool LinearElastic<TType>::Has( const Variable<int>& rThisVariable )
{
    if(rThisVariable == PARENT_ELEMENT_ID)
        return true;
    if(rThisVariable == INTEGRATION_POINT_INDEX)
        return true;
    return false;
}

template<int TType>
bool LinearElastic<TType>::Has( const Variable<double>& rThisVariable )
{
    if ( rThisVariable == PRESTRESS_FACTOR || rThisVariable == YOUNG_MODULUS || rThisVariable == POISSON_RATIO )
        return true;

    return false;
}

template<int TType>
bool LinearElastic<TType>::Has( const Variable<array_1d<double, 3> >& rThisVariable )
{
    return false;
}

template<int TType>
bool LinearElastic<TType>::Has( const Variable<Vector>& rThisVariable )
{
    if ( rThisVariable == INSITU_STRESS )
        return true;

    if ( rThisVariable == PRESTRESS )
        return true;

    if ( rThisVariable == STRESSES )
        return true;

    if ( rThisVariable == STRAIN )
        return true;

    return false;
}

template<int TType>
bool LinearElastic<TType>::Has( const Variable<Matrix>& rThisVariable )
{
    if( rThisVariable == ALGORITHMIC_TANGENT )
        return true;
    return false;
}

template<int TType>
int& LinearElastic<TType>::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    if(rThisVariable == PARENT_ELEMENT_ID)
        rValue = mElemId;
    if(rThisVariable == INTEGRATION_POINT_INDEX)
        rValue = mGaussId;
    if (rThisVariable == IS_SHAPE_FUNCTION_REQUIRED)
        rValue = 0;
    return rValue;
}

template<int TType>
double& LinearElastic<TType>::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if ( rThisVariable == PRESTRESS_FACTOR )
    {
        rValue = mPrestressFactor;
        return rValue;
    }

    if(rThisVariable == YOUNG_MODULUS )
    {
       rValue = mE;
       return rValue;
    }

    if ( rThisVariable == POISSON_RATIO )
    {
        rValue = mNU;
        return rValue;
    }

    if (rThisVariable == PRESSURE_P)
    {
        rValue = -(m_stress_n(0, 0) + m_stress_n(1, 1) + m_stress_n(2, 2)) / 3;
        return rValue;
    }

    if (rThisVariable == PRESSURE_Q || rThisVariable == VON_MISES_STRESS)
    {
        double p = (m_stress_n(0, 0) + m_stress_n(1, 1) + m_stress_n(2, 2)) / 3;

        double sxx = m_stress_n(0, 0) - p;
        double syy = m_stress_n(1, 1) - p;
        double szz = m_stress_n(2, 2) - p;
        double sxy = m_stress_n(0, 1);
        double syz = m_stress_n(1, 2);
        double sxz = m_stress_n(0, 2);

        double j2 = 0.5*( pow(sxx, 2) + pow(syy, 2) + pow(szz, 2) ) + pow(sxy, 2) + pow(syz, 2) + pow(sxz, 2);
        rValue = sqrt( 3.0 * j2 );
        return rValue;
    }

    if ( rThisVariable == EQUIVALENT_VOLUMETRIC_STRAIN )
    {
        if (TType == 0) // plane strain
        {
            rValue = (mOldStrain[0] + mOldStrain[1]);
        }
        else if (TType == 1) // plane stress
        {
            rValue = (1.0-2.0*mNU)/(1.0-mNU)*(mOldStrain[0] + mOldStrain[1]);
        }
        else if (TType == 2) // axisymmetric
        {
            rValue = (mOldStrain[0] + mOldStrain[1] + mOldStrain[3]);
        }
        else if (TType == 3) // 3D
        {
            rValue = (mOldStrain[0] + mOldStrain[1] + mOldStrain[2]);
        }
        return rValue;
    }

    if ( rThisVariable == EQUIVALENT_DEVIATORIC_STRAIN )
    {
        double sxx = mOldStrain[0];
        double syy = mOldStrain[1];
        double szz, sxy, syz, sxz;
        if (TType == 0) // plane strain
        {
            szz = sxz = syz = 0.0; sxy = mOldStrain[2];
        }
        else if (TType == 1) // plane stress
        {
            szz = -mNU/(1.0-mNU)*(mOldStrain[0] + mOldStrain[1]);
            sxy = mOldStrain[2];
            sxz = syz = 0.0;
        }
        else if (TType == 2) // axisymmetric
        {
            szz = mOldStrain[3];
            sxy = mOldStrain[2];
            sxz = syz = 0.0;
        }
        else if (TType == 3) // 3D
        {
            szz = mOldStrain[2];
            sxy = mOldStrain[3];
            syz = mOldStrain[4];
            sxz = mOldStrain[5];
        }

        // REF: https://dianafea.com/manuals/d944/Analys/node405.html
        double exx = 2*sxx/3 - syy/3 - szz/3;
        double eyy = -sxx/3 + 2*syy/3 - szz/3;
        double ezz = -sxx/3 - syy/3 + 2*szz/3;

        rValue = (2.0/3) * std::sqrt(1.5*(pow(exx, 2) + pow(eyy, 2) + pow(ezz, 2) + 3*(pow(sxy, 2) + pow(syz, 2) + pow(sxz, 2))));

        return rValue;
    }

    // REF: https://en.wikipedia.org/wiki/Elastic_modulus#:~:text=An%20elastic%20modulus%20(also%20known,stress%20is%20applied%20to%20it.
    if ( rThisVariable == BULK_MODULUS )
    {
        rValue = mE/(3.0*(1.0-2.0*mNU));
        return rValue;
    }

    if ( rThisVariable == SHEAR_MODULUS )
    {
        rValue = mE/(2.0*(1+mNU));
        return rValue;
    }

    return rValue;
}

template<int TType>
array_1d<double, 3>& LinearElastic<TType>::GetValue( const Variable<array_1d<double, 3> >& rThisVariable, array_1d<double, 3>& rValue )
{
    return rValue;
}

template<int TType>
Vector& LinearElastic<TType>::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    if ( rThisVariable == INSITU_STRESS || rThisVariable == PRESTRESS )
    {
        const unsigned int size = mPrestress.size();
        if(rValue.size() != size)
            rValue.resize(size, false );
        noalias(rValue) = mPrestress;
        return rValue;
    }

    if ( rThisVariable == STRESSES_OLD )
    {
        const unsigned int size = mPrestress.size();
        if(rValue.size() != size)
            rValue.resize(size, false );
        StructuralMechanicsMathUtilities<double>::StressTensorToVector(m_stress_n, rValue);
        return rValue;
    }

    if ( rThisVariable == STRAIN )
    {
        const unsigned int size = mOldStrain.size();
        if(rValue.size() != size)
            rValue.resize(size, false );
        noalias(rValue) = mOldStrain;
        return rValue;
    }

    if ( rThisVariable == STRESSES )
    {
        const unsigned int size = mCurrentStress.size();
        if(rValue.size() != size)
            rValue.resize(size, false );
        noalias(rValue) = mCurrentStress;
        return rValue;
    }

    if ( rThisVariable == THREED_STRESSES )
    {
        rValue = ZeroVector(6);
        StructuralMechanicsMathUtilities<double>::StressTensorToVector(m_stress_n, rValue);
        return rValue;
    }

    return rValue;
}

template<int TType>
Matrix& LinearElastic<TType>::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    if(rThisVariable == ALGORITHMIC_TANGENT || rThisVariable == ELASTIC_TANGENT)
    {
        const std::size_t strain_size = this->GetStrainSize();
        if(rValue.size1() != strain_size || rValue.size2() != strain_size)
            rValue.resize(strain_size, strain_size, false);
        LinearElastic_Helper<TType>::CalculateElasticMatrix( rValue, mE, mNU );
    }
    return rValue;
}

template<int TType>
double& LinearElastic<TType>::CalculateValue( Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue )
{
    return this->GetValue(rThisVariable, rValue);
}

template<int TType>
array_1d<double, 3>& LinearElastic<TType>::CalculateValue( Parameters& rParameterValues, const Variable<array_1d<double, 3> >& rThisVariable, array_1d<double, 3>& rValue )
{
    return this->GetValue(rThisVariable, rValue);
}

template<int TType>
Vector& LinearElastic<TType>::CalculateValue( Parameters& rParameterValues, const Variable<Vector>& rThisVariable, Vector& rValue )
{
    return this->GetValue(rThisVariable, rValue);
}

template<int TType>
Matrix& LinearElastic<TType>::CalculateValue( Parameters& rParameterValues, const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    return this->GetValue(rThisVariable, rValue);
}

template<int TType>
void LinearElastic<TType>::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if(rThisVariable == PARENT_ELEMENT_ID)
        mElemId = rValue;
    if(rThisVariable == INTEGRATION_POINT_INDEX)
        mGaussId = rValue;
}

template<int TType>
void LinearElastic<TType>::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == PRESTRESS_FACTOR )
        mPrestressFactor = rValue;
    if ( rThisVariable == YOUNG_MODULUS )
        mE = rValue;
    if ( rThisVariable == POISSON_RATIO )
        mNU = rValue;
}

template<int TType>
void LinearElastic<TType>::SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                            const array_1d<double, 3>& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}

template<int TType>
void LinearElastic<TType>::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == INSITU_STRESS || rThisVariable == PRESTRESS )
    {
        if(mPrestress.size() != rValue.size())
            mPrestress.resize(rValue.size(), false);
        noalias(mPrestress) = rValue;
    }
    else if ( rThisVariable == STRESSES )
    {
        if(mCurrentStress.size() != rValue.size())
            mCurrentStress.resize(rValue.size(), false);
        noalias(mCurrentStress) = rValue;
        StructuralMechanicsMathUtilities<double>::StressVectorToTensor(rValue, m_stress_n);
        if (TType == 0) // plane strain
            m_stress_n(2, 2) = mNU*(m_stress_n(0, 0) + m_stress_n(1, 1));
    }
    else if ( rThisVariable == INITIAL_STRESS )
    {
        if(mCurrentStress.size() != rValue.size())
            mCurrentStress.resize(rValue.size(), false);
        noalias(mCurrentStress) = rValue;
        StructuralMechanicsMathUtilities<double>::StressVectorToTensor(-rValue, m_stress_n);
        if (TType == 0) // plane strain
            m_stress_n(2, 2) = mNU*(m_stress_n(0, 0) + m_stress_n(1, 1));
    }
    else if ( rThisVariable == STRAIN )
    {
        if(mCurrentStrain.size() != rValue.size())
            mCurrentStrain.resize(rValue.size(), false);
        noalias(mCurrentStrain) = rValue;
        if(mOldStrain.size() != rValue.size())
            mOldStrain.resize(rValue.size(), false);
        noalias(mOldStrain) = rValue;
    }
}

template<int TType>
void LinearElastic<TType>::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}

template<int TType>
void LinearElastic<TType>::InitializeMaterial( const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues )
{
    std::size_t strain_size = this->GetStrainSize();
    m_stress_n = ZeroMatrix( 3, 3 );
    mCurrentStress = ZeroVector( strain_size );
    mCurrentStrain = ZeroVector( strain_size );
    mOldStrain = ZeroVector( strain_size );
    mPrestress = ZeroVector( strain_size );
    mPrestressFactor = 1.0;
    mE = props[YOUNG_MODULUS];
    mNU = props[POISSON_RATIO];
}

template<int TType>
void LinearElastic<TType>::ResetMaterial()
{
    std::size_t strain_size = this->GetStrainSize();
    m_stress_n = ZeroMatrix( 3, 3 );
    StructuralMechanicsMathUtilities<double>::StressVectorToTensor(-mPrestressFactor*mPrestress, m_stress_n);
    if (TType == 0) // plane strain
        m_stress_n(2, 2) = mNU*(m_stress_n(0, 0) + m_stress_n(1, 1));
    mCurrentStress = ZeroVector( strain_size );
    StructuralMechanicsMathUtilities<double>::StressTensorToVector(m_stress_n, mCurrentStress);
    mCurrentStrain = ZeroVector( strain_size );
    mOldStrain = ZeroVector( strain_size );
}

template<int TType>
void LinearElastic<TType>::ResetMaterial( const Properties& props,
                                 const GeometryType& geom,
                                 const Vector& ShapeFunctionsValues )
{
    this->ResetMaterial();
}

template<int TType>
void LinearElastic<TType>::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

template<int TType>
void LinearElastic<TType>::InitializeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

template<int TType>
void LinearElastic<TType>::FinalizeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

template<int TType>
void LinearElastic<TType>::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    StructuralMechanicsMathUtilities<double>::StressVectorToTensor(mCurrentStress, m_stress_n);
    if (TType == 0) // plane strain
        m_stress_n(2, 2) = mNU*(m_stress_n(0, 0) + m_stress_n(1, 1));
    noalias(mOldStrain) = mCurrentStrain;
}

template<int TType>
void LinearElastic<TType>::CalculateMaterialResponseCauchy( Parameters& parameters )
{
    this->CalculateMaterialResponse( parameters.GetStrainVector()
        , parameters.GetDeformationGradientF()
        , parameters.GetStressVector()
        , parameters.GetConstitutiveMatrix()
        , parameters.GetProcessInfo()
        , parameters.GetMaterialProperties()
        , parameters.GetElementGeometry()
        , parameters.GetShapeFunctionsValues()
        , parameters.IsSetStressVector()
        , parameters.IsSetConstitutiveMatrix()
        , true
    );
}

template<int TType>
void LinearElastic<TType>::CalculateMaterialResponsePK2( Parameters& rValues )
{
    CalculateMaterialResponseCauchy(rValues);
}

template<int TType>
void LinearElastic<TType>::CalculateMaterialResponse( const Vector& StrainVector,
        const Matrix& DeformationGradient,
        Vector& StressVector,
        Matrix& AlgorithmicTangent,
        const ProcessInfo& CurrentProcessInfo,
        const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues,
        bool CalculateStresses,
        int CalculateTangent,
        bool SaveInternalVariables )
{
    if (CalculateStresses || CalculateTangent)
        LinearElastic_Helper<TType>::CalculateElasticMatrix( AlgorithmicTangent, mE, mNU );

    if (CalculateStresses)
    {
        CalculateStress( StrainVector, AlgorithmicTangent, mCurrentStress );
        noalias(StressVector) = mCurrentStress;
        noalias(mCurrentStrain) = StrainVector;
    }
}

template<int TType>
void LinearElastic<TType>::CalculateStress( const Vector& StrainVector, Matrix& AlgorithmicTangent, Vector& StressVector )
{
    if ( StressVector.size() != this->GetStrainSize() )
    {
        StressVector.resize( this->GetStrainSize() );
    }

    StructuralMechanicsMathUtilities<double>::StressTensorToVector(m_stress_n, StressVector);
    noalias( StressVector ) += prod( AlgorithmicTangent, StrainVector - mOldStrain );
}

//**********************************************************************
template<int TType>
int LinearElastic<TType>::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY

    if ( !props.Has( YOUNG_MODULUS ) || !props.Has(POISSON_RATIO) )
    {
        KRATOS_THROW_ERROR( std::logic_error, "this constitutive law requires YOUNG_MODULUS and POISSON_RATIO given as KRATOS variables", "" );
    }

    return 0;

    KRATOS_CATCH( "" );
}

// Explicit template instantiation
template class LinearElastic<0>;
template class LinearElastic<1>;
template class LinearElastic<2>;
template class LinearElastic<3>;

} // Namespace Kratos
