//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

// Application includes
#include "custom_constitutive/history_linear_elastic_3D_law.hpp"

namespace Kratos
{

int HistoryLinearElastic3DLaw::Check(const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const ProcessInfo& rCurrentProcessInfo)
{
    int ierr = BaseType::Check(rMaterialProperties,rElementGeometry,rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    if ( INITIAL_STRESS_VECTOR.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"INITIAL_STRESS_VECTOR Key is 0. Check if all applications were correctly registered.", "" )
    if ( STEP_INITIAL_STRESS.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"STEP_INITIAL_STRESS Key is 0. Check if all applications were correctly registered.", "" )

    return ierr;
}

//----------------------------------------------------------------------------------------

void HistoryLinearElastic3DLaw::InitializeMaterial(
                        const Properties& rMaterialProperties,
                        const GeometryType& rElementGeometry,
                        const Vector& rShapeFunctionsValues
                        )
{
    BaseType::InitializeMaterial(rMaterialProperties,rElementGeometry,rShapeFunctionsValues);

    unsigned int voigt_size = GetStrainSize();
    if(mInitialStressVector.size() != voigt_size) {
        mInitialStressVector.resize(voigt_size,false);
    }
    noalias(mInitialStressVector) = ZeroVector(voigt_size);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void HistoryLinearElastic3DLaw::CalculateStress( const Vector & rStrainVector,
					  const Matrix & rConstitutiveMatrix,
					  Vector& rStressVector )
{
    BaseType::CalculateStress(rStrainVector,rConstitutiveMatrix,rStressVector);

    noalias(rStressVector) += mInitialStressVector;

    // TODO: the printed stress tensor in the STEP = STEP_INITIAL_STRESS will be wrong !
}

//----------------------------------------------------------------------------------------

void HistoryLinearElastic3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
{
    BaseType::FinalizeMaterialResponseCauchy(rValues);

    if(rValues.GetProcessInfo()[STEP] == rValues.GetProcessInfo()[STEP_INITIAL_STRESS])
    {
        mInitialStressVector = rValues.GetStressVector();
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Vector& HistoryLinearElastic3DLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    if(rThisVariable==INITIAL_STRESS_VECTOR) {
        rValue = mInitialStressVector;
    }
    // TODO: I don't understand why this is an error
    // else {
    //     return BaseType::GetValue(rThisVariable,rValue);
    // }

    return rValue;
}

//----------------------------------------------------------------------------------------

void HistoryLinearElastic3DLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                                        const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == INITIAL_STRESS_VECTOR) {
        mInitialStressVector = rValue;
    } else {
        BaseType::SetValue(rThisVariable,rValue,rCurrentProcessInfo);
    }
}

} // Namespace Kratos
